#!/usr/bin/env Rscript

# Independently rebuild a frozen-input DEC run at the selected finite candidate's estimates,
# recalculate its likelihood, and extract ancestral state posteriors.  The real
# optimizer object remains in fit_result.rds; this fixed-parameter post-fit
# object is never treated as optimizer evidence. R4 v1 uses the explicitly
# authorized Pelopidas mathias legacy-code exception and otherwise new geography.

suppressPackageStartupMessages({
  library(ape)
  library(BioGeoBEARS)
  library(methods)
  library(digest)
  library(jsonlite)
})

is_symlink <- function(path) {
  link <- Sys.readlink(path)
  !is.na(link) && nzchar(link)
}

parse_args <- function(x) {
  z <- list(); i <- 1L
  while (i <= length(x)) {
    if (!startsWith(x[[i]], "--") || i == length(x)) stop("Use --key value arguments")
    key <- sub("^--", "", x[[i]])
    if (!is.null(z[[key]])) stop("Duplicate argument key")
    z[[key]] <- x[[i + 1L]]; i <- i + 2L
  }
  z
}

script_arg <- grep("^--file=", commandArgs(FALSE), value=TRUE)
if (length(script_arg) != 1L) stop("Cannot resolve script path")
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork=TRUE)
job_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork=TRUE)
opts <- parse_args(commandArgs(trailingOnly=TRUE))
if (!setequal(names(opts),c("fit-tag", "cores", "outtag"))) stop("Required exactly: --fit-tag TAG --cores N --outtag TAG")
if (!grepl("^[A-Za-z0-9_.-]+$",opts[["fit-tag"]]) || !grepl("^[A-Za-z0-9_.-]+$",opts$outtag)) stop("Invalid fit/output tag")
fit_path <- file.path(job_root,"04_runs",opts[["fit-tag"]],"fit_result.rds")
if (!file.exists(fit_path) || is_symlink(fit_path)) stop("Expected regular fit file")
fit_fn <- normalizePath(fit_path, mustWork=TRUE)
ncores <- as.integer(opts$cores)
if (!startsWith(fit_fn, paste0(normalizePath(file.path(job_root, "04_runs")), "/"))) {
  stop("fit-rds must be inside this job's 04_runs directory")
}
if (basename(fit_fn) != "fit_result.rds") stop("fit-rds basename must be fit_result.rds")
if (!is.finite(ncores) || ncores < 1L || ncores > 6L) stop("cores must be 1..6")

outdir <- file.path(job_root,"04_runs",opts$outtag)
if (file.exists(outdir) || dir.exists(outdir) || is_symlink(outdir)) {
  stop("Refusing any existing postfit path: ", outdir)
}
run_parent <- normalizePath(dirname(outdir),mustWork=TRUE)
if (!startsWith(run_parent,paste0(job_root,"/")) || is_symlink(dirname(outdir))) stop("Output parent outside job or symlink")
if (!dir.create(outdir, recursive=FALSE, showWarnings=FALSE)) stop("Cannot create new postfit directory")
status_fn <- file.path(outdir, "STATUS.tsv")
now <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
write_status <- function(status, detail="") {
  write.table(data.frame(status=status, timestamp=now(), detail=detail), status_fn,
              sep="\t", quote=FALSE, row.names=FALSE)
}
write_status("RUNNING", "loading_real_optimizer_result")

sha256_file <- function(path) digest(file=path, algo="sha256", serialize=FALSE)
snapshot_file <- function(path) {
  canonical <- normalizePath(path, mustWork=TRUE)
  raw <- readBin(canonical, what="raw", n=file.info(canonical)$size)
  list(path=canonical, sha256=digest(raw, algo="sha256", serialize=FALSE),
       size_bytes=length(raw), raw=raw)
}
manifest_path <- file.path(job_root,"01_inputs/frozen/input_manifest.json")
if (sha256_file(manifest_path)!="f7e121a3ce4aff304767eef6df39eba6f9ef473242aad48411a3b2c204d4a306") stop("New input manifest SHA changed")
manifest <- jsonlite::read_json(manifest_path,simplifyVector=FALSE)
if (!identical(manifest$schema_version,"geography_pelopidas_legacy_exception_v1") || !identical(manifest$job_root,job_root)) stop("Manifest schema/job-root mismatch")
rec <- manifest$recipe
if (rec$n_tips!=417 || rec$n_areas!=11 || rec$max_range_size!=4 || !isTRUE(rec$include_null_range) ||
    !identical(unlist(rec$models,use.names=FALSE),c("M0","M1","M2"))) stop("Manifest recipe mismatch")
if (as.character(getRversion())!="4.4.0" || as.character(packageVersion("BioGeoBEARS"))!="1.1.3" ||
    as.character(packageVersion("optimx"))!="2025.4.9") stop("Runtime versions differ")
expected <- character()
for (item in manifest$files) {
  rel <- item$relative_path
  if (!is.character(rel) || length(rel)!=1L || startsWith(rel,"/") ||
      grepl("(^|/)\\.\\.(/|$)",rel) || grepl("\\\\",rel)) stop("Unsafe manifest path")
  p <- file.path(job_root,rel)
  if (!file.exists(p) || dir.exists(p) || is_symlink(p) ||
      !startsWith(normalizePath(p,mustWork=TRUE),paste0(job_root,"/"))) stop("Manifest file must be regular and inside job")
  if (basename(p) %in% names(expected)) stop("Duplicate manifest basename")
  if (sha256_file(p)!=item$sha256 || file.info(p)$size!=item$bytes) stop("Manifest file mismatch: ",rel)
  expected[basename(p)] <- item$sha256
}
expected["input_manifest.json"] <- sha256_file(manifest_path)
expected["04_fit_dec_model_R4_v1.R"] <- "3c196903191d9b026705e18b44a8fe2652d7987bd25c5dff227c9f5a23be2f4f"

assert_hash <- function(path) {
  key <- basename(path); observed <- sha256_file(path)
  if (!(key %in% names(expected)) || !identical(observed, unname(expected[[key]]))) {
    stop("Frozen-file hash mismatch: ", path, " observed=", observed)
  }
}

atomic_save <- function(x, path) {
  tmp <- paste0(path, ".tmp"); saveRDS(x, tmp, version=3)
  if (!file.rename(tmp, path)) stop("Atomic rename failed: ", path)
}

blank_path <- function(x) {
  if (is.null(x) || (is.character(x) && length(x) == 0L)) return(TRUE)
  if (length(x) != 1L || !is.atomic(x)) return(FALSE)
  if (is.na(x)) return(TRUE)
  is.character(x) && identical(x, "")
}

verify_provenance_snapshots <- function(snapshots, expected_paths) {
  if (!is.list(snapshots) || length(snapshots) != length(expected_paths)) {
    stop("Final fit provenance snapshot count mismatch")
  }
  normalized <- vapply(expected_paths, normalizePath, character(1), mustWork=TRUE)
  observed <- character(0)
  for (path in normalized) {
    idx <- which(vapply(snapshots, function(x) {
      is.list(x) && !is.null(x$path) &&
        identical(normalizePath(as.character(x$path), mustWork=TRUE), path)
    }, logical(1)))
    if (length(idx) != 1L) stop("Final fit provenance does not contain exactly one snapshot for ", path)
    snap <- snapshots[[idx]]
    live_sha <- sha256_file(path)
    if (!identical(as.character(snap$sha256), live_sha) || !is.raw(snap$raw) ||
        !identical(digest(snap$raw, algo="sha256", serialize=FALSE), live_sha) ||
        as.numeric(snap$size_bytes) != length(snap$raw)) {
      stop("Final fit provenance snapshot/raw/live-file mismatch for ", path)
    }
    observed <- c(observed, path)
  }
  snap_paths <- vapply(snapshots, function(x) {
    if (!is.list(x) || is.null(x$path)) return(NA_character_)
    normalizePath(as.character(x$path), mustWork=TRUE)
  }, character(1))
  if (anyNA(snap_paths) || anyDuplicated(snap_paths) || !setequal(snap_paths, observed)) {
    stop("Final fit provenance snapshot path set has duplicates, omissions, or extras")
  }
  invisible(TRUE)
}

result <- tryCatch({
  fit_sha_before <- sha256_file(fit_fn)
  fit <- readRDS(fit_fn)
  if (is.null(fit$total_loglikelihood) || !is.finite(as.numeric(fit$total_loglikelihood))) {
    stop("Accepted fit has no finite total_loglikelihood")
  }
  if (is.null(fit$outputs) || !isS4(fit$outputs) ||
      !("params_table" %in% methods::slotNames(fit$outputs))) stop("Accepted fit lacks params table")
  if (is.null(fit$inputs$provenance)) stop("Accepted fit lacks embedded provenance")
  pv <- fit$inputs$provenance
  if (!identical(pv$mode, "final") || !isTRUE(pv$final_fit) ||
      !identical(pv$speedup, FALSE) || !identical(pv$force_sparse, FALSE) ||
      !identical(fit$inputs$speedup, FALSE) || !identical(fit$inputs$force_sparse, FALSE)) {
    stop("Post-fit is allowed only for a final speedup=FALSE, force_sparse=FALSE fit")
  }
  if (!isTRUE(pv$include_null_range) || !isTRUE(fit$inputs$include_null_range)) {
    stop("Final fit does not preserve include_null_range=TRUE")
  }
  model <- pv$model
  if (!(model %in% c("M0", "M1", "M2"))) stop("Unexpected model in fit provenance")
  max_range <- as.integer(pv$max_range_size)
  if (max_range != 4L || !identical(pv$numerical_backend,"native")) stop("Expected native full-R4 final fit")

  tree_fn <- file.path(job_root, "01_inputs", "frozen", "tree_scenarioA417.tre")
  geog_fn <- file.path(job_root, "01_inputs", "frozen", "geog_scenarioA417_analysis_order.LagrangePHYLIP")
  times_fn <- file.path(job_root, "02_config", "timeperiods_5epochs.txt")
  area_order_fn <- file.path(job_root, "02_config", "area_order.tsv")
  input_manifest_fn <- file.path(job_root, "01_inputs", "frozen", "input_manifest.json")
  runner_fn <- file.path(job_root, "03_scripts", "04_fit_dec_model_R4_v1.R")
  mult_fn <- switch(model,
    M1=file.path(job_root, "02_config", "M1_conservative_dispersal_multipliers.txt"),
    M2=file.path(job_root, "02_config", "M2_permissive_dispersal_multipliers.txt"),
    M0=NA_character_)
  assert_hash(tree_fn); assert_hash(geog_fn); assert_hash(area_order_fn)
  assert_hash(input_manifest_fn); assert_hash(runner_fn)
  if (model != "M0") { assert_hash(times_fn); assert_hash(mult_fn) }

  if (!identical(normalizePath(fit$inputs$trfn, mustWork=TRUE), normalizePath(tree_fn, mustWork=TRUE)) ||
      !identical(normalizePath(fit$inputs$geogfn, mustWork=TRUE), normalizePath(geog_fn, mustWork=TRUE))) {
    stop("Final fit runtime tree/geography path mismatch")
  }
  if (model == "M0") {
    if (!blank_path(fit$inputs$timesfn) || !blank_path(fit$inputs$dispersal_multipliers_fn) ||
        length(fit$inputs$timeperiods) != 0L || length(fit$inputs$list_of_dispersal_multipliers_mats) != 0L) {
      stop("M0 final fit unexpectedly contains stratification inputs")
    }
  } else if (!identical(normalizePath(fit$inputs$timesfn, mustWork=TRUE), normalizePath(times_fn, mustWork=TRUE)) ||
             !identical(normalizePath(fit$inputs$dispersal_multipliers_fn, mustWork=TRUE), normalizePath(mult_fn, mustWork=TRUE))) {
    stop("Time-stratified final fit runtime path mismatch")
  }

  expected_snapshot_paths <- c(tree_fn, geog_fn, area_order_fn, input_manifest_fn,
                               runner_fn)
  if (model != "M0") expected_snapshot_paths <- c(expected_snapshot_paths, times_fn, mult_fn)
  verify_provenance_snapshots(pv$files, expected_snapshot_paths)

  expected_remote_sha <- "1672cc0c171b1a05effad69fa426427b3d9ef4e1"
  if (!identical(as.character(pv$BioGeoBEARS_version), "1.1.3") ||
      !identical(unname(as.character(pv$BioGeoBEARS_remote_sha)), expected_remote_sha) ||
      !identical(as.character(packageVersion("BioGeoBEARS")), "1.1.3") ||
      !identical(unname(packageDescription("BioGeoBEARS")$RemoteSha), expected_remote_sha)) {
    stop("Final fit or post-fit environment uses an unreviewed BioGeoBEARS version")
  }
  expected_package_versions <- c(
    ape=as.character(packageVersion("ape")), optimx=as.character(packageVersion("optimx")),
    minqa=as.character(packageVersion("minqa")), digest=as.character(packageVersion("digest"))
  )
  if (!identical(unname(as.character(pv$package_versions[names(expected_package_versions)])),
                 unname(expected_package_versions))) {
    stop("Final fit dependency versions differ from the post-fit environment")
  }

  pt_fit <- fit$outputs@params_table
  if (!identical(rownames(pt_fit)[pt_fit[, "type"] == "free"], c("d", "e"))) {
    stop("Final fit does not have exactly d/e free")
  }
  if (as.numeric(pt_fit["j", "est"]) != 0 || as.numeric(pt_fit["w", "est"]) != 1) {
    stop("Final fit violates j=0 or w=1")
  }
  d_hat <- as.numeric(pt_fit["d", "est"]); e_hat <- as.numeric(pt_fit["e", "est"])
  if (any(!is.finite(c(d_hat, e_hat))) ||
      any(c(d_hat,e_hat) < as.numeric(pt_fit[c("d","e"),"min"])) ||
      any(c(d_hat,e_hat) > as.numeric(pt_fit[c("d","e"),"max"]))) {
    stop("Final fit d/e estimates are not finite positive values")
  }

  run <- define_BioGeoBEARS_run()
  run$trfn <- normalizePath(tree_fn, mustWork=TRUE)
  run$geogfn <- normalizePath(geog_fn, mustWork=TRUE)
  run$max_range_size <- max_range
  run$min_branchlength <- 1e-6
  run$include_null_range <- TRUE
  run$on_NaN_error <- -1e50
  run$force_sparse <- FALSE
  run$speedup <- FALSE
  run$use_optimx <- TRUE
  run$num_cores_to_use <- ncores
  run$return_condlikes_table <- TRUE
  run$calc_TTL_loglike_from_condlikes_table <- TRUE
  run$calc_ancprobs <- TRUE
  if (model != "M0") {
    run$timesfn <- normalizePath(times_fn, mustWork=TRUE)
    run$dispersal_multipliers_fn <- normalizePath(mult_fn, mustWork=TRUE)
  }
  pt <- run$BioGeoBEARS_model_object@params_table
  pt["d", "type"] <- "free"; pt["e", "type"] <- "free"
  pt["d", c("init", "est")] <- d_hat; pt["e", c("init", "est")] <- e_hat
  pt["j", "type"] <- "fixed"; pt["j", c("init", "est")] <- 0
  pt["w", "type"] <- "fixed"; pt["w", c("init", "est")] <- 1
  run$BioGeoBEARS_model_object@params_table <- pt
  run <- readfiles_BioGeoBEARS_run(run)
  if (!is.null(run$states_list)) stop("Restricted state list not authorized")
  geog <- getranges_from_LagrangePHYLIP(geog_fn)
  g <- as.matrix(geog@df)
  if (any(!g %in% c(0,1,"0","1"))) stop("Nonbinary geography")
  storage.mode(g)<-"numeric"
  tr_check <- ape::read.tree(tree_fn)
  if (nrow(g)!=417L || length(tr_check$tip.label)!=417L || !setequal(tr_check$tip.label,rownames(g)) ||
      !identical(colnames(g),LETTERS[1:11]) || any(rowSums(g)==0) || max(rowSums(g))!=4) stop("Fresh Pelopidas-legacy exception geography QA failed")
  if (model != "M0") {
    if (!identical(as.numeric(run$timeperiods), c(3.3, 13.9, 23.03, 33.9, 45))) stop("Time periods changed")
    if (length(run$list_of_dispersal_multipliers_mats) != 5L) stop("Matrix count changed")
    run <- section_the_tree(inputs=run, make_master_table=TRUE, plot_pieces=FALSE,
                            cut_fossils=FALSE, min_dist_between_node_and_stratum_line=1e-5)
  }
  if (!isTRUE(check_BioGeoBEARS_run(run))) stop("Independent rebuilt run failed package check")
  run$postfit_provenance <- list(
    real_optimizer_fit_path=fit_fn,
    real_optimizer_fit_sha256=sha256_file(fit_fn),
    model=model, max_range_size=max_range, d=d_hat, e=e_hat,
    fixed_j=0, fixed_w=1, speedup=FALSE, force_sparse=FALSE,
    frozen_input_sha256=lapply(c(tree_fn, geog_fn, area_order_fn, input_manifest_fn,
                                 if (model == "M0") character(0) else c(times_fn, mult_fn)),
                               sha256_file),
    BioGeoBEARS_version=as.character(packageVersion("BioGeoBEARS")),
    BioGeoBEARS_remote_sha=unname(packageDescription("BioGeoBEARS")$RemoteSha),
    started=now(), script_path=script_path, script_sha256=sha256_file(script_path),
    script_snapshot=snapshot_file(script_path)
  )
  atomic_save(run, file.path(outdir, "independently_rebuilt_postfit_run.rds"))

  # Stage evidence only: this does not certify a completed likelihood or posterior.
  marker <- file.path(outdir,"POSTFIT_ENTERED.json")
  if (file.exists(marker) || dir.exists(marker) || is_symlink(marker)) stop("Postfit entry marker exists")
  jsonlite::write_json(list(schema_version="postfit_entry_v2",model=model,
    fit_tag=opts[["fit-tag"]],outtag=opts$outtag,
    postfit_runner_sha256=sha256_file(script_path),fit_sha256=sha256_file(fit_fn),
    input_manifest_sha256=sha256_file(input_manifest_fn),timestamp=now()),
    marker,auto_unbox=TRUE,pretty=TRUE,digits=NA)
  started <- Sys.time()
  postfit_warnings <- character()
  postfit <- withCallingHandlers(bears_optim_run(run, skip_optim=TRUE, skip_optim_option="return_all"),warning=function(w){
    postfit_warnings <<- c(postfit_warnings,conditionMessage(w))
  })
  postfit$postfit_warnings <- postfit_warnings
  postfit$scientific_acceptance <- "NOT_GRANTED_TECHNICAL_CONDITIONAL_RESULT"
  postfit$source_optimizer_diagnostics <- list(optim_result=fit$optim_result,bounded_optimizer=fit$bounded_optimizer,runner_warnings=fit$runner_warnings)
  elapsed <- as.numeric(difftime(Sys.time(), started, units="secs"))
  postfit$postfit_completion <- list(finished=now(), elapsed_seconds=elapsed,
                                     real_optimizer_fit_sha256=sha256_file(fit_fn))
  postfit_fn <- file.path(outdir, "postfit_recalculated_ancestral_states.rds")

  ll_fit <- as.numeric(fit$total_loglikelihood)
  ll_check <- as.numeric(postfit$total_loglikelihood)
  delta_ll <- ll_check - ll_fit
  if (length(ll_fit)!=1L || length(ll_check)!=1L || !is.finite(ll_fit) || !is.finite(ll_check) ||
      ll_fit>0 || ll_check>0 || ll_fit<=-1e49 || ll_check<=-1e49 || abs(delta_ll) > 1e-6) {
    stop(sprintf("Independent likelihood mismatch: fit=%.12f check=%.12f delta=%.12g", ll_fit, ll_check, delta_ll))
  }

  top <- as.matrix(postfit$ML_marginal_prob_each_state_at_branch_top_AT_node)
  bottom <- as.matrix(postfit$ML_marginal_prob_each_state_at_branch_bottom_below_node)
  tr <- ape::read.tree(tree_fn)
  n_total <- length(tr$tip.label) + tr$Nnode
  root_nodes <- setdiff(unique(as.integer(tr$edge[, 1])),
                        unique(as.integer(tr$edge[, 2])))
  if (length(root_nodes) != 1L ||
      !identical(root_nodes, as.integer(length(tr$tip.label) + 1L))) {
    stop("Frozen tree does not have the unique expected ape root node")
  }
  root_node <- root_nodes[[1]]
  states <- postfit$inputs$all_geog_states_list_usually_inferred_from_areas_maxareas
  expected_states <- 562L
  exact_states <- cladoRcpp::rcpp_areas_list_to_states_list(areas=LETTERS[1:11],maxareas=4L,include_null_range=TRUE)
  if (!identical(states,exact_states)) stop("Full R4 state dictionary order/content differs")
  if (!identical(dim(top), c(n_total, expected_states)) || !identical(dim(bottom), c(n_total, expected_states))) {
    stop("Unexpected ancestral probability dimensions")
  }
  check_probs <- function(x, label, undefined_root=FALSE) {
    rows <- seq_len(nrow(x))
    if (isTRUE(undefined_root)) {
      # BioGeoBEARS defines branch-bottom probabilities on branches.  The
      # root has no ancestral branch, so its entire bottom row is
      # structurally NA.  No partial root row or additional undefined row is
      # permitted.
      if (!all(is.na(x[root_node, , drop=FALSE])) ||
          any(is.nan(x[root_node, , drop=FALSE]))) {
        stop(label, " root row is not entirely structurally undefined")
      }
      rows <- setdiff(rows, root_node)
    }
    defined <- x[rows, , drop=FALSE]
    if (any(!is.finite(defined))) stop(label, " contains non-finite values outside the root")
    if (min(defined) < -1e-12 || max(defined) > 1 + 1e-12) {
      stop(label, " contains values outside [0,1]")
    }
    err <- max(abs(rowSums(defined) - 1))
    if (err > 1e-6) stop(label, " row-sum error exceeds 1e-6: ", err)
    err
  }
  top_err <- check_probs(top, "top probabilities", undefined_root=FALSE)
  bottom_err <- check_probs(bottom, "bottom probabilities", undefined_root=TRUE)
  if (length(states) != expected_states) stop("State dictionary length mismatch")

  area_map <- read.delim(file.path(job_root, "02_config", "area_order.tsv"),
                         sep="\t", stringsAsFactors=FALSE, check.names=FALSE)
  area_abbrev <- area_map$abbrev[order(area_map$analysis_bit_index)]
  internal_codes <- area_map$analysis_internal_code[order(area_map$analysis_bit_index)]
  decode <- function(s, labels) {
    if (length(s) == 0L || all(is.na(s))) return("NULL")
    idx <- as.integer(s) + 1L
    if (any(!is.finite(idx)) || any(idx < 1L) || any(idx > length(labels))) stop("Invalid 0-based state index")
    paste(labels[idx], collapse="+")
  }
  state_df <- data.frame(
    state_index=seq_along(states),
    state_zero_based=seq_along(states) - 1L,
    internal_state=vapply(states, decode, character(1), labels=internal_codes),
    semantic_state=vapply(states, decode, character(1), labels=area_abbrev),
    range_size=vapply(states, function(s) if (length(s) == 0L || all(is.na(s))) 0L else length(s), integer(1)),
    stringsAsFactors=FALSE
  )
  write.table(state_df, file.path(outdir, "state_dictionary.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

  node_depth <- ape::node.depth.edgelength(tr)
  root_age <- max(node_depth[seq_len(length(tr$tip.label))])
  node_df <- data.frame(
    ape_node=seq_len(n_total),
    node_type=ifelse(seq_len(n_total) <= length(tr$tip.label), "tip", "internal"),
    tip_label=c(tr$tip.label, rep(NA_character_, tr$Nnode)),
    age_ma=pmax(0, root_age - node_depth),
    max_state_probability=apply(top, 1, max),
    entropy_nats=apply(top, 1, function(p) -sum(ifelse(p > 0, p * log(p), 0))),
    stringsAsFactors=FALSE
  )
  write.table(node_df, file.path(outdir, "node_metadata.tsv"), sep="\t", quote=FALSE, row.names=FALSE, na="")

  top_rows <- vector("list", n_total)
  for (node in seq_len(n_total)) {
    ord <- order(top[node, ], decreasing=TRUE)[1:3]
    probs <- top[node, ord]
    top_rows[[node]] <- rbind(
      data.frame(ape_node=node, rank=1:3, state_index=ord,
                 semantic_state=state_df$semantic_state[ord], probability=probs,
                 stringsAsFactors=FALSE),
      data.frame(ape_node=node, rank=4L, state_index=NA_integer_,
                 semantic_state="OTHER", probability=max(0, 1 - sum(probs)),
                 stringsAsFactors=FALSE)
    )
  }
  top3 <- do.call(rbind, top_rows); rownames(top3) <- NULL
  write.table(top3, file.path(outdir, "node_top3_states_plus_other.tsv"), sep="\t", quote=FALSE, row.names=FALSE, na="")
  atomic_save(list(top=top, bottom=bottom, state_dictionary=state_df),
              file.path(outdir, "ancestral_probability_matrices.rds"))
  # Recheck long-running read dependencies before committing the new result.
  if (sha256_file(fit_fn)!=fit_sha_before) stop("Fit source changed during postfit")
  if (sha256_file(manifest_path)!="f7e121a3ce4aff304767eef6df39eba6f9ef473242aad48411a3b2c204d4a306") stop("Manifest changed during postfit")
  for (item in manifest$files) assert_hash(file.path(job_root,item$relative_path))
  assert_hash(runner_fn)
  # Commit only after likelihood/probability invariants pass; no optimizer or
  # scientific-acceptance conclusion is inferred from this technical check.
  atomic_save(postfit, postfit_fn)

  qa <- data.frame(
    status="TECHNICAL_PASS_CONDITIONAL_ONLY", scientific_acceptance="NOT_GRANTED",
    warning_count=length(postfit_warnings),
    model=model, max_range=max_range, n_tips=length(tr$tip.label),
    n_internal=tr$Nnode, n_states=length(states), fit_lnL=ll_fit,
    recalculated_lnL=ll_check, delta_lnL=delta_ll,
    top_max_rowsum_error=top_err, bottom_max_rowsum_error=bottom_err,
    bottom_undefined_root_node=root_node, bottom_undefined_row_count=1L,
    elapsed_seconds=elapsed, optimizer_fit_rds_sha256=sha256_file(fit_fn),
    postfit_rds_sha256=sha256_file(postfit_fn), postfit_script_sha256=sha256_file(script_path),
    stringsAsFactors=FALSE
  )
  write.table(qa, file.path(outdir, "postfit_QA.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
  write_status("TECHNICAL_PASS_CONDITIONAL_ONLY", sprintf("delta_lnL=%.12g;states=%d;top_rowsum_error=%.3g;bottom_root_na=%d",
                               delta_ll, length(states), top_err, root_node))
  cat(sprintf("POSTFIT_TECHNICAL_PASS model=%s max_range=%d fit_lnL=%.12f check_lnL=%.12f delta=%.12g states=%d bottom_root_na=%d elapsed_seconds=%.3f\n",
              model, max_range, ll_fit, ll_check, delta_ll, length(states),
              root_node, elapsed))
}, error=function(e) {
  msg <- conditionMessage(e); write_status("FAILED", msg)
  writeLines(c(now(), msg, capture.output(traceback())), file.path(outdir, "ERROR.txt"))
  message("POSTFIT_FAILED: ", msg); quit(status=1L)
})
