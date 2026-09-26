#!/usr/bin/env Rscript

# One explicit BioGeoBEARS DEC fit from frozen inputs.  This script never
# searches neighbouring directories for previous results and refuses to write
# into any existing output path. R4 v1 retains the NA-safe path guard.
# Pelopidas mathias alone retains its legacy code by explicit user instruction;
# all other evidence-supported corrections remain frozen in the new manifest.

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
  out <- list()
  i <- 1L
  while (i <= length(x)) {
    if (!startsWith(x[[i]], "--") || i == length(x)) {
      stop("Arguments must be supplied as --key value pairs")
    }
    key <- sub("^--", "", x[[i]])
    if (!is.null(out[[key]])) stop("Duplicate argument: ", key)
    out[[key]] <- x[[i + 1L]]
    i <- i + 2L
  }
  out
}

required <- function(opts, keys) {
  absent <- setdiff(keys, names(opts))
  if (length(absent) > 0L) stop("Missing arguments: ", paste(absent, collapse=", "))
}

script_arg <- grep("^--file=", commandArgs(FALSE), value=TRUE)
if (length(script_arg) != 1L) stop("Cannot resolve runner script path")
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork=TRUE)
job_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork=TRUE)

opts <- parse_args(commandArgs(trailingOnly=TRUE))
required(opts, c("model", "start", "mode", "max-range", "cores", "d", "e", "outtag", "maxit"))
if (length(setdiff(names(opts),c("model","start","mode","max-range","cores","d","e","outtag","maxit")))) stop("Unknown argument")

model <- toupper(opts$model)
start_id <- opts$start
mode <- tolower(opts$mode)
max_range <- as.integer(opts[["max-range"]])
ncores <- as.integer(opts$cores)
d0 <- as.numeric(opts$d)
e0 <- as.numeric(opts$e)
outtag <- opts$outtag
maxit <- as.integer(opts$maxit)
if (!is.finite(maxit) || maxit < 10L || maxit > 500L) stop("maxit must be 10..500")

if (!(model %in% c("M0", "M1", "M2"))) stop("model must be M0, M1, or M2")
if (!(mode %in% c("prepare", "smoke", "screen", "final"))) stop("mode must be prepare, smoke, screen, or final")
if (max_range != 4L) stop("This Pelopidas-legacy exception recipe requires max-range 4")
if (!is.finite(ncores) || ncores < 1L || ncores > 6L) stop("cores must be between 1 and 6")
if (!is.finite(d0) || !is.finite(e0) || d0 <= 0 || e0 <= 0) stop("d and e must be finite and > 0")
if (!grepl("^[A-Za-z0-9_.-]+$", outtag)) stop("outtag contains unsafe characters")

tree_fn <- file.path(job_root, "01_inputs", "frozen", "tree_scenarioA417.tre")
geog_fn <- file.path(job_root, "01_inputs", "frozen", "geog_scenarioA417_analysis_order.LagrangePHYLIP")
times_fn <- file.path(job_root, "02_config", "timeperiods_5epochs.txt")
mult_fn <- switch(model,
  M1=file.path(job_root, "02_config", "M1_conservative_dispersal_multipliers.txt"),
  M2=file.path(job_root, "02_config", "M2_permissive_dispersal_multipliers.txt"),
  M0=NA_character_
)
area_order_fn <- file.path(job_root, "02_config", "area_order.tsv")
input_manifest_fn <- file.path(job_root, "01_inputs", "frozen", "input_manifest.json")

# New immutable manifest: only Pelopidas mathias deliberately uses legacy coding.
# No old fitted result or geography-dependent cache is reused.
if (digest(file=input_manifest_fn,algo="sha256",serialize=FALSE)!="f7e121a3ce4aff304767eef6df39eba6f9ef473242aad48411a3b2c204d4a306") stop("Frozen input manifest changed")
manifest <- jsonlite::read_json(input_manifest_fn,simplifyVector=FALSE)
if (!identical(manifest$schema_version,"geography_pelopidas_legacy_exception_v1") || !identical(manifest$job_root,job_root)) stop("Manifest schema/job-root mismatch")
if (as.character(getRversion())!="4.4.0" || as.character(packageVersion("BioGeoBEARS"))!="1.1.3" ||
    packageDescription("BioGeoBEARS")$RemoteSha!="1672cc0c171b1a05effad69fa426427b3d9ef4e1" ||
    as.character(packageVersion("optimx"))!="2025.4.9") stop("Runtime identity differs from reviewed implementation")
rec <- manifest$recipe
if (is.null(rec) || rec$n_tips!=417 || rec$n_areas!=11 || rec$max_range_size!=4 ||
    !isTRUE(rec$include_null_range) ||
    !identical(unlist(rec$models,use.names=FALSE),c("M0","M1","M2")) ||
    !identical(unlist(rec$area_order,use.names=FALSE),c("AF","AUS","CAM","ENA","EPA","IND","MDG","ORI","SAM","WNA","WPA"))) stop("Manifest recipe differs")
expected_hashes <- character()
manifest_paths <- character()
if (!length(manifest$files)) stop("Manifest has no files")
for (item in manifest$files) {
  rel <- item$relative_path
  if (!is.character(rel) || length(rel)!=1L || startsWith(rel,"/") ||
      grepl("(^|/)\\.\\.(/|$)",rel) || grepl("\\\\",rel)) stop("Unsafe manifest relative path")
  path <- file.path(job_root,rel)
  if (!file.exists(path) || dir.exists(path) || is_symlink(path) ||
      !startsWith(normalizePath(path,mustWork=TRUE),paste0(job_root,"/"))) stop("Manifest file must be a regular file inside this job")
  if (rel %in% manifest_paths || basename(path) %in% names(expected_hashes)) stop("Duplicate manifest path/basename")
  actual <- digest(file=path,algo="sha256",serialize=FALSE)
  if (!is.character(item$sha256) || !grepl("^[0-9a-f]{64}$",item$sha256) ||
      actual!=item$sha256 || !is.numeric(item$bytes) || file.info(path)$size!=item$bytes) stop("Manifest binding mismatch: ",rel)
  manifest_paths <- c(manifest_paths,rel)
  expected_hashes[basename(path)] <- item$sha256
}
immutable <- c(tree_scenarioA417.tre="06725d5a0c1ef29007aeb0f7657175c2d2329c75cb4d8fb2170e39b5fdaf6962",
  timeperiods_5epochs.txt="90dadc7755385b44d9333f22e2acaeac163bd0ce4507289e1720671eccc4e51f",
  M1_conservative_dispersal_multipliers.txt="34cc3ee7f0789d978948e35085fb9ed89de09f2600b50384bb074bd11db7029c",
  M2_permissive_dispersal_multipliers.txt="aad709ab6571d541292c640efd74e0cf167959f1fa2d7cefa3a200b37953b745")
if (!identical(unname(expected_hashes[names(immutable)]),unname(immutable))) stop("Tree/time/matrix identity changed")
if (expected_hashes[["geog_scenarioA417_analysis_order.LagrangePHYLIP"]]!="3aa227e4600629eb77178ab68a30bb7cf408189bbcb9f8989ae1f184ce834ca4") stop("Not the authorized 68-bit corrected geography with Pelopidas legacy exception")

sha256_file <- function(path) digest(file=path, algo="sha256", serialize=FALSE)
assert_hash <- function(path) {
  path <- normalizePath(path, mustWork=TRUE)
  key <- basename(path)
  if (!(key %in% names(expected_hashes))) stop("No frozen expected hash for ", key)
  observed <- sha256_file(path)
  if (!identical(unname(expected_hashes[[key]]), observed)) {
    stop("SHA256 mismatch for ", path, ": ", observed)
  }
  invisible(observed)
}

assert_hash(tree_fn)
assert_hash(geog_fn)
if (model != "M0") {
  assert_hash(times_fn)
  assert_hash(mult_fn)
}

outdir <- file.path(job_root, "04_runs", outtag)
if (file.exists(outdir) || dir.exists(outdir) || is_symlink(outdir)) {
  stop("Refusing any existing run path: ", outdir)
}
run_parent <- normalizePath(dirname(outdir),mustWork=TRUE)
if (!startsWith(run_parent,paste0(job_root,"/")) || is_symlink(dirname(outdir))) stop("Run parent outside job or symlink")
if (!dir.create(outdir, recursive=FALSE, showWarnings=FALSE)) stop("Run directory creation failed")

timestamp <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
status_fn <- file.path(outdir, "STATUS.tsv")
write_status <- function(status, detail="") {
  write.table(
    data.frame(status=status, timestamp=timestamp(), model=model, start=start_id,
               mode=mode, max_range=max_range, cores=ncores, detail=detail,
               stringsAsFactors=FALSE),
    status_fn, sep="\t", quote=FALSE, row.names=FALSE
  )
}
write_status("RUNNING", "building_fresh_run_object")

snapshot_file <- function(path) {
  p <- normalizePath(path, mustWork=TRUE)
  size <- unname(file.info(p)$size)
  list(path=p, size_bytes=size, sha256=sha256_file(p),
       raw=readBin(p, what="raw", n=size))
}

build_run <- function() {
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
  run$num_cores_to_use <- as.integer(ncores)
  run$return_condlikes_table <- FALSE
  run$calc_TTL_loglike_from_condlikes_table <- TRUE
  run$calc_ancprobs <- FALSE

  if (model != "M0") {
    run$timesfn <- normalizePath(times_fn, mustWork=TRUE)
    run$dispersal_multipliers_fn <- normalizePath(mult_fn, mustWork=TRUE)
  }

  pt <- run$BioGeoBEARS_model_object@params_table
  if (d0 < as.numeric(pt["d","min"]) || d0 > as.numeric(pt["d","max"]) ||
      e0 < as.numeric(pt["e","min"]) || e0 > as.numeric(pt["e","max"])) stop("Initial d/e outside installed parameter bounds")
  pt["d", "type"] <- "free"
  pt["e", "type"] <- "free"
  pt["d", c("init", "est")] <- d0
  pt["e", c("init", "est")] <- e0
  pt["j", "type"] <- "fixed"
  pt["j", c("init", "est")] <- 0
  pt["w", "type"] <- "fixed"
  pt["w", c("init", "est")] <- 1
  run$BioGeoBEARS_model_object@params_table <- pt

  free_names <- rownames(pt)[pt[, "type"] == "free"]
  if (!identical(free_names, c("d", "e"))) {
    stop("Unexpected free parameters: ", paste(free_names, collapse=", "))
  }

  run <- readfiles_BioGeoBEARS_run(run)
  tr <- ape::read.tree(tree_fn)
  geog <- getranges_from_LagrangePHYLIP(geog_fn)
  if (length(tr$tip.label) != 417L || nrow(geog@df) != 417L) stop("Expected 417 tree/geography taxa")
  if (!setequal(tr$tip.label, rownames(geog@df))) stop("Tree/geography taxon-set mismatch")
  if (!identical(colnames(geog@df),LETTERS[1:11]) || anyDuplicated(rownames(geog@df))) stop("Area/taxon identity mismatch")
  geogmat <- as.matrix(geog@df)
  if (any(!geogmat %in% c(0,1,"0","1"))) stop("Geography is not binary")
  storage.mode(geogmat)<-"numeric"
  if (any(rowSums(geogmat)==0) || max(rowSums(geogmat))!=4) stop("Expected nonempty observed ranges and maximum four")
  states <- cladoRcpp::rcpp_areas_list_to_states_list(areas=LETTERS[1:11],maxareas=4L,include_null_range=TRUE)
  if (length(states)!=562L) stop("Unexpected R4 state count")
  if (!is.null(run$states_list)) stop("Restricted states list not authorized")

  if (model != "M0") {
    expected_times <- c(3.3, 13.9, 23.03, 33.9, 45)
    if (!isTRUE(all.equal(as.numeric(run$timeperiods), expected_times, tolerance=0))) {
      stop("Parsed timeperiods differ from frozen specification")
    }
    mats <- run$list_of_dispersal_multipliers_mats
    if (length(mats) != length(expected_times)) stop("Expected exactly five dispersal matrices")
    allowed <- if (model == "M1") c(1e-7, 0.01, 0.1, 1) else c(0.01, 0.1, 0.5, 1)
    expected_area_header <- LETTERS[1:11]
    for (i in seq_along(mats)) {
      m <- as.matrix(mats[[i]])
      if (!identical(dim(m), c(11L, 11L))) stop("Matrix ", i, " is not 11x11")
      if (!identical(colnames(m), expected_area_header)) stop("Matrix ", i, " internal header is not A:K")
      if (any(!is.finite(m)) || any(m < 0)) stop("Invalid values in matrix ", i)
      if (max(abs(m - t(m))) > 1e-15) stop("Matrix ", i, " is not symmetric")
      if (max(abs(diag(m) - 1)) > 1e-15) stop("Matrix ", i, " diagonal is not 1")
      if (any(!m %in% allowed)) stop("Unexpected multiplier value in matrix ", i)
    }
    run <- section_the_tree(
      inputs=run,
      make_master_table=TRUE,
      plot_pieces=FALSE,
      cut_fossils=FALSE,
      min_dist_between_node_and_stratum_line=1e-5
    )
  }

  if (!isTRUE(check_BioGeoBEARS_run(run))) stop("check_BioGeoBEARS_run did not return TRUE")

  provenance_paths <- c(tree_fn, geog_fn, area_order_fn,
                        input_manifest_fn, script_path)
  if (model != "M0") provenance_paths <- c(provenance_paths, times_fn, mult_fn)
  run$provenance <- list(
    schema_version="1.0",
    model=model,
    start=list(id=start_id, d=d0, e=e0),
    mode=mode,
    final_fit=identical(mode, "final"),
    max_range_size=max_range,
    include_null_range=TRUE,
    force_sparse=FALSE,
    speedup=run$speedup,
    num_cores_to_use=ncores,
    numerical_backend="native",
    state_space="ALL_RANGES_UP_TO_4_INCLUDING_NULL_562_STATES",
    optimization_maxfun=maxit,
    optimizer_scope="Single bounded bobyqa; no automatic restarts",
    input_manifest_sha256=sha256_file(input_manifest_fn),
    expected_free_parameters=c("d", "e"),
    j_fixed=0,
    w_fixed=1,
    files=lapply(provenance_paths, snapshot_file),
    BioGeoBEARS_version=as.character(packageVersion("BioGeoBEARS")),
    BioGeoBEARS_remote_sha=unname(packageDescription("BioGeoBEARS")$RemoteSha),
    package_versions=c(
      ape=as.character(packageVersion("ape")),
      optimx=as.character(packageVersion("optimx")),
      minqa=as.character(packageVersion("minqa")),
      digest=as.character(packageVersion("digest"))
    ),
    session_info=capture.output(sessionInfo()),
    started=timestamp()
  )
  run
}

atomic_save_rds <- function(x, path) {
  tmp <- paste0(path, ".tmp")
  saveRDS(x, tmp, version=3)
  if (!file.rename(tmp, path)) stop("Atomic rename failed for ", path)
}

extract_optimizer <- function(res) {
  opt <- tryCatch(as.data.frame(res$optim_result), error=function(e) data.frame())
  get1 <- function(nm) if (nm %in% names(opt) && nrow(opt) > 0L) opt[[nm]][1] else NA
  data.frame(
    optimizer_row=if (nrow(opt) > 0L) rownames(opt)[1] else NA_character_,
    value=as.numeric(get1("value")),
    convcode=as.integer(get1("convcode")),
    kkt1=as.logical(get1("kkt1")),
    kkt2=as.logical(get1("kkt2")),
    fevals=as.integer(get1("fevals")),
    optimization_budget_reached=if(is.finite(as.numeric(get1("fevals")))) as.numeric(get1("fevals")) >= maxit else NA,
    native_exit_status="UNAVAILABLE_IN_OPTIMX_RESULT",
    runtime_seconds=as.numeric(if ("xtime" %in% names(opt)) get1("xtime") else get1("xtimes")),
    stringsAsFactors=FALSE
  )
}

extract_loglikelihood <- function(res, optdf=extract_optimizer(res)) {
  explicit_names <- c("total_loglikelihood", "ML_loglik", "total_loglik",
                      "LnL", "loglik", "lnL")
  explicit <- lapply(explicit_names, function(nm) {
    if (is.null(res[[nm]])) return(NULL)
    value <- suppressWarnings(as.numeric(res[[nm]]))
    if (length(value) == 1L && is.finite(value)) value else NULL
  })
  explicit <- unlist(explicit, use.names=FALSE)

  objective <- suppressWarnings(as.numeric(optdf$value))
  if (length(objective) != 1L || !is.finite(objective)) {
    stop("Optimizer objective is not one finite number")
  }
  # This exact pinned optimx uses maximize=TRUE and returns log-likelihood.
  # Do not silently negate an invalid positive result or accept the NaN sentinel.
  if (objective > 0 || objective <= -1e49) stop("Optimizer returned invalid positive/sentinel log-likelihood")
  objective_lnL <- objective

  if (length(explicit) > 0L) {
    if (any(abs(explicit - explicit[[1]]) > 1e-6) ||
        abs(explicit[[1]] - objective_lnL) > 1e-6) {
      stop("Explicit likelihood and optimizer objective disagree")
    }
    return(as.numeric(explicit[[1]]))
  }
  as.numeric(objective_lnL)
}

# Clone only the top-level dispatcher environment. Its function body and all
# likelihood/Q/node/exponential functions remain byte-identical installed native.
make_bounded_bears <- function(budget) {
  native_bears <- BioGeoBEARS::bears_optim_run
  native_optimx <- optimx::optimx
  native_body_sha <- digest(paste(deparse(body(native_bears)),collapse="\n"),algo="sha256",serialize=FALSE)
  control_env <- new.env(parent=asNamespace("BioGeoBEARS"))
  call_record <- new.env(parent=emptyenv())
  call_record$n <- 0L
  control_env$optimx <- function(...) {
    call_record$n <- call_record$n+1L
    if (call_record$n!=1L) stop("More than one optimizer call is not authorized")
    args <- list(...)
    if (!identical(args$method,c("bobyqa"))) stop("Unexpected optimizer method")
    args$itnmax <- as.integer(budget)
    args$control$maxit <- as.integer(budget)
    args$control$all.methods <- FALSE
    args$control$maximize <- TRUE
    call_record$itnmax <- args$itnmax
    call_record$control <- args$control
    marker <- file.path(outdir,"OPTIMIZER_ENTERED.json")
    if (file.exists(marker) || dir.exists(marker) || is_symlink(marker)) stop("Optimizer entry marker already exists")
    jsonlite::write_json(list(schema_version="optimizer_entry_v2", model=model,
      start=start_id, outtag=outtag, runner_sha256=sha256_file(script_path),
      input_manifest_sha256=sha256_file(input_manifest_fn),
      timestamp=timestamp(), maxfun=budget, optimizer_calls=call_record$n),
      marker,auto_unbox=TRUE,pretty=TRUE,digits=NA)
    do.call(native_optimx,args)
  }
  bounded_bears <- native_bears
  environment(bounded_bears) <- control_env
  stopifnot(identical(body(bounded_bears),body(native_bears)),
    identical(formals(bounded_bears),formals(native_bears)),
    identical(environment(native_bears),asNamespace("BioGeoBEARS")))
  list(run=bounded_bears,record=call_record,native_body_sha256=native_body_sha)
}
bounded <- make_bounded_bears(maxit)
run_started <- Sys.time()
result <- tryCatch({
  run <- build_run()
  atomic_save_rds(run, file.path(outdir, "prefit_run.rds"))
  write.table(run$BioGeoBEARS_model_object@params_table,
              file.path(outdir, "prefit_params.tsv"), sep="\t", quote=FALSE,
              row.names=TRUE, col.names=NA)

  if (mode == "prepare") {
    if (bounded$record$n != 0L) stop("Prepare must not call the optimizer")
    write_status("VALIDATED_NOT_FITTED", "fresh_input_build_only;zero_likelihood_calls;zero_optimizer_calls")
    cat(sprintf("PREPARE_PASS model=%s max_range=%d states=562 likelihood_calls=0 optimizer_calls=0\n", model, max_range))
    quit(status=0L)
  }

  if (mode == "smoke") {
    ll <- bears_optim_run(run, skip_optim=TRUE, skip_optim_option="return_loglike")
    ll <- as.numeric(ll)
    if (length(ll) != 1L || !is.finite(ll) || ll > 0 || ll <= -1e49) stop("Initial likelihood is invalid")
    elapsed <- as.numeric(difftime(Sys.time(), run_started, units="secs"))
    smoke <- data.frame(model=model, start=start_id, d=d0, e=e0,
                        max_range=max_range, lnL=ll, elapsed_seconds=elapsed,
                        stringsAsFactors=FALSE)
    write.table(smoke, file.path(outdir, "smoke_likelihood.tsv"), sep="\t",
                quote=FALSE, row.names=FALSE)
    write_status("PASS", sprintf("finite_initial_lnL=%.12f", ll))
    cat(sprintf("SMOKE_PASS model=%s start=%s max_range=%d lnL=%.12f elapsed_seconds=%.3f\n",
                model, start_id, max_range, ll, elapsed))
    quit(status=0L)
  }

  write_status("RUNNING", "native_driver_started_no_optimizer_confirmation")
  write.table(data.frame(model=model, start=start_id, timestamp=timestamp(),
    phase="ENTERING_NATIVE_DRIVER_NOT_A_COMPLETED_OPTIMIZATION"),
    file.path(outdir,"optimization_entered.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
  fit_warnings <- character()
  fit <- withCallingHandlers(bounded$run(run), warning=function(w) {
    fit_warnings <<- c(fit_warnings,conditionMessage(w))
  })
  fit$runner_warnings <- fit_warnings
  if (bounded$record$n!=1L || bounded$record$itnmax!=maxit) stop("Bounded optimizer control did not execute")
  fit$bounded_optimizer <- list(method="bobyqa",maxfun=maxit,optimizer_calls=bounded$record$n,
    native_exit_status="UNAVAILABLE_IN_OPTIMX_RESULT",scientific_acceptance="NOT_GRANTED",
    native_bears_body_sha256=bounded$native_body_sha256,control=bounded$record$control,
    interpretation="maxfun bounds bobyqa evaluations; initial test, diagnostics and final likelihood add finite evaluations")
  fit$run_completion <- list(finished=timestamp(), elapsed_seconds=as.numeric(difftime(Sys.time(), run_started, units="secs")))

  if (is.null(fit$outputs) || !isS4(fit$outputs) ||
      !("params_table" %in% methods::slotNames(fit$outputs))) stop("Result lacks outputs@params_table")
  pt <- fit$outputs@params_table
  fitted_de <- as.numeric(pt[c("d","e"),"est"])
  if (any(!is.finite(fitted_de)) || any(fitted_de < as.numeric(pt[c("d","e"),"min"])) ||
      any(fitted_de > as.numeric(pt[c("d","e"),"max"]))) stop("Fitted d/e outside installed parameter bounds")
  optdf <- extract_optimizer(fit)
  lnL <- extract_loglikelihood(fit, optdf)
  fit$total_loglikelihood <- lnL
  if (sha256_file(input_manifest_fn)!="f7e121a3ce4aff304767eef6df39eba6f9ef473242aad48411a3b2c204d4a306") stop("Manifest changed during fit")
  for (rel in manifest_paths) assert_hash(file.path(job_root,rel))
  atomic_save_rds(fit, file.path(outdir, "fit_result.rds"))
  write.table(pt, file.path(outdir, "fit_params.tsv"), sep="\t", quote=FALSE,
              row.names=TRUE, col.names=NA)
  d_hat <- as.numeric(pt["d", "est"])
  e_hat <- as.numeric(pt["e", "est"])
  k <- sum(pt[, "type"] == "free")
  n <- 417L
  aic <- 2 * k - 2 * lnL
  aicc <- aic + 2 * k * (k + 1) / (n - k - 1)
  write.table(optdf, file.path(outdir, "optimizer_diagnostics.tsv"), sep="\t",
              quote=FALSE, row.names=FALSE)
  summary <- data.frame(
    model=model, start=start_id, mode=mode, max_range=max_range,
    lnL=lnL, k=k, AIC=aic, AICc=aicc, n_tips=n,
    d=d_hat, e=e_hat, j=as.numeric(pt["j", "est"]), w=as.numeric(pt["w", "est"]),
    convcode=optdf$convcode, kkt1=optdf$kkt1, kkt2=optdf$kkt2,
    fevals=optdf$fevals, optimization_budget=maxit,
    optimization_budget_reached=optdf$optimization_budget_reached,
    warning_count=length(fit_warnings),
    native_exit_status=optdf$native_exit_status,
    elapsed_seconds=fit$run_completion$elapsed_seconds,
    stringsAsFactors=FALSE
  )
  write.table(summary, file.path(outdir, "fit_summary.tsv"), sep="\t",
              quote=FALSE, row.names=FALSE)
  preliminary_pass <- isTRUE(optdf$convcode == 0L) && is.finite(optdf$fevals) && length(fit_warnings)==0L &&
    !isTRUE(optdf$optimization_budget_reached) && all(is.finite(c(d_hat, e_hat))) &&
    identical(k, 2L) && isTRUE(as.numeric(pt["j", "est"]) == 0) &&
    isTRUE(as.numeric(pt["w", "est"]) == 1)
  # Pinned optimx does not retain minqa's native ierr/msg; reported convcode=0
  # is not proof of a successful native exit, even in the absence of warnings.
  write_status(if (preliminary_pass) "FINITE_CANDIDATE" else "FINITE_CANDIDATE_WITH_DIAGNOSTIC_FLAGS",
               sprintf("lnL=%.12f;convcode=%s;kkt1=%s;kkt2=%s",
                       lnL, optdf$convcode, optdf$kkt1, optdf$kkt2))
  cat(sprintf("FIT_DONE model=%s start=%s mode=%s max_range=%d lnL=%.12f d=%.12g e=%.12g convcode=%s kkt1=%s kkt2=%s elapsed_seconds=%.3f\n",
              model, start_id, mode, max_range, lnL, d_hat, e_hat,
              optdf$convcode, optdf$kkt1, optdf$kkt2, fit$run_completion$elapsed_seconds))
  invisible(summary)
}, error=function(e) {
  msg <- conditionMessage(e)
  write_status("FAILED", msg)
  writeLines(c(timestamp(), msg, capture.output(traceback())), file.path(outdir, "ERROR.txt"))
  message("RUN_FAILED: ", msg)
  quit(status=1L)
})
