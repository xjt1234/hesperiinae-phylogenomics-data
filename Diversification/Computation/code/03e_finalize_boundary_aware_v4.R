#!/usr/bin/env Rscript

# Boundary-aware finalization of the prespecified ten-model RPANDA analysis.
#
# Nine models have reproducible finite maximum-likelihood estimates. The
# BCSTDTempVar likelihood was separately diagnosed, under the unchanged
# likelihood and a one-to-one log-centred parameterization, as improving toward
# beta -> -Inf with explosive extinction rates. It therefore has no identified
# finite regular MLE. This script retains that candidate as an explicit tenth
# row, excludes it from ordinary AICc ranking/weight normalization, and carries
# the complete diagnostic evidence without treating a numerical endpoint or a
# likelihood supremum as a maximum-likelihood estimate.

suppressPackageStartupMessages(library(ape))
options(stringsAsFactors = FALSE, digits = 17, scipen = 999)

args_full <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", grep("^--file=", args_full, value = TRUE))
if (length(script_arg) != 1L) stop("Cannot determine script path")
script_file <- normalizePath(script_arg, mustWork = TRUE)
run_dir <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)

round1_file <- file.path(run_dir, "03_models", "model_fits_10models.rds")
tree_file <- file.path(run_dir, "01_tree", "T25_deep_only_dated.Hesperiinae_417sp_RPANDA.tre")
environment_file <- file.path(run_dir, "02_environment", "temperature_surface_to_tree_age.tsv")
profile_model_dir <- file.path(run_dir, "03_models", "BCSTDTempVar_logcenter_profile_v4")
profile_table_dir <- file.path(run_dir, "04_tables", "BCSTDTempVar_logcenter_profile_v4")
profile_file <- file.path(profile_model_dir, "BCSTDTempVar_logcenter_profile.rds")
profile_signature_file <- file.path(profile_model_dir, "analysis_signature.tsv")
profile_signature_hash_file <- file.path(profile_model_dir, "analysis_signature.sha256")
profile_diagnosis_file <- file.path(profile_table_dir, "boundary_diagnosis.tsv")
profile_summary_file <- file.path(profile_table_dir, "fixed_beta_profile_summary.tsv")
profile_runs_file <- file.path(profile_table_dir, "fixed_beta_profile_runs.tsv")
profile_full_runs_file <- file.path(profile_table_dir, "full_multistart_runs.tsv")

final_name <- "final_v4_boundary_aware"
final_model_dir <- file.path(run_dir, "03_models", final_name)
final_table_dir <- file.path(run_dir, "04_tables", final_name)
dir.create(final_model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(final_table_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(
  round1_file, tree_file, environment_file, profile_file,
  profile_signature_file, profile_signature_hash_file,
  profile_diagnosis_file, profile_summary_file,
  profile_runs_file, profile_full_runs_file
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop(
    "Boundary-aware finalization requires a completed round-one fit and 03c profile:\n",
    paste0(" - ", missing_files, collapse = "\n")
  )
}

sha256 <- function(path) {
  answer <- system2("sha256sum", path, stdout = TRUE, stderr = TRUE)
  status <- attr(answer, "status")
  if (!length(answer) || (!is.null(status) && status != 0L)) {
    stop("sha256sum failed for: ", path)
  }
  hash <- strsplit(answer[[1L]], "[[:space:]]+")[[1L]][[1L]]
  if (!grepl("^[0-9a-f]{64}$", hash)) stop("Invalid SHA-256 for: ", path)
  hash
}

clean_text <- function(x) {
  y <- enc2utf8(as.character(x))
  y <- sub("^\\ufeff", "", y, perl = TRUE)
  y <- gsub("[\\r\\n\\t]+", " ", y, perl = TRUE)
  y <- gsub("[ ]{2,}", " ", y, perl = TRUE)
  trimws(y)
}

clean_frame <- function(x) {
  names(x) <- clean_text(names(x))
  if (any(!nzchar(names(x))) || anyDuplicated(names(x))) {
    stop("Empty or duplicated field names after cleaning")
  }
  for (nm in names(x)) {
    if (is.character(x[[nm]]) || is.factor(x[[nm]])) x[[nm]] <- clean_text(x[[nm]])
  }
  x
}

has_control_text <- function(x) {
  any(vapply(x, function(value) {
    (is.character(value) || is.factor(value)) &&
      any(grepl("[\\r\\n\\t]", as.character(value), perl = TRUE), na.rm = TRUE)
  }, logical(1)))
}

read_tsv_strict <- function(path) {
  raw <- readLines(path, warn = FALSE, encoding = "UTF-8")
  if (length(raw) < 2L) stop("TSV is empty or lacks rows: ", path)
  if (any(grepl("\\r", raw, fixed = TRUE))) stop("TSV contains carriage returns: ", path)
  value <- read.delim(
    path, sep = "\t", header = TRUE, quote = "", comment.char = "",
    check.names = FALSE, stringsAsFactors = FALSE,
    na.strings = c("NA", "NaN", ""), fill = FALSE, fileEncoding = "UTF-8"
  )
  if (nrow(value) != length(raw) - 1L) stop("TSV row/line mismatch: ", path)
  clean_frame(value)
}

preflight_tsv <- function(x, parent) {
  x <- clean_frame(x)
  if (has_control_text(x)) return(FALSE)
  tmp <- tempfile(pattern = ".tsv_preflight_", tmpdir = parent, fileext = ".tsv")
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  write.table(
    x, tmp, sep = "\t", quote = FALSE, row.names = FALSE,
    col.names = TRUE, na = "NA", fileEncoding = "UTF-8"
  )
  raw <- readLines(tmp, warn = FALSE, encoding = "UTF-8")
  back <- tryCatch(
    read.delim(
      tmp, sep = "\t", header = TRUE, quote = "", comment.char = "",
      check.names = FALSE, stringsAsFactors = FALSE,
      na.strings = c("NA", "NaN", ""), fill = FALSE, fileEncoding = "UTF-8"
    ),
    error = function(e) NULL
  )
  length(raw) == nrow(x) + 1L && !is.null(back) &&
    nrow(back) == nrow(x) && ncol(back) == ncol(x) &&
    identical(names(back), names(x))
}

write_clean_tsv <- function(x, path) {
  x <- clean_frame(x)
  if (!preflight_tsv(x, dirname(path))) stop("TSV preflight failed: ", path)
  tmp <- tempfile(pattern = paste0(".", basename(path), "_"),
                  tmpdir = dirname(path), fileext = ".tmp")
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  write.table(
    x, tmp, sep = "\t", quote = FALSE, row.names = FALSE,
    col.names = TRUE, na = "NA", fileEncoding = "UTF-8"
  )
  if (file.exists(path)) {
    if (!identical(sha256(path), sha256(tmp))) {
      stop("Refusing to overwrite a different existing final-v4 TSV: ", path)
    }
    unlink(tmp)
    return(invisible(path))
  }
  if (!file.rename(tmp, path)) stop("Could not atomically place TSV: ", path)
  invisible(path)
}

write_text_no_overwrite <- function(lines, path) {
  tmp <- tempfile(pattern = paste0(".", basename(path), "_"),
                  tmpdir = dirname(path), fileext = ".tmp")
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  writeLines(enc2utf8(lines), tmp, useBytes = TRUE)
  if (file.exists(path)) {
    if (!identical(sha256(path), sha256(tmp))) {
      stop("Refusing to overwrite a different existing final-v4 text file: ", path)
    }
    unlink(tmp)
    return(invisible(path))
  }
  if (!file.rename(tmp, path)) stop("Could not atomically place text file: ", path)
  invisible(path)
}

atomic_save_rds <- function(object, path) {
  if (file.exists(path)) stop("Refusing to overwrite existing final-v4 RDS: ", path)
  tmp <- tempfile(pattern = paste0(".", basename(path), "_"),
                  tmpdir = dirname(path), fileext = ".tmp")
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  saveRDS(object, tmp, version = 3)
  if (!file.rename(tmp, path)) stop("Could not atomically place RDS: ", path)
  invisible(path)
}

fmt_values <- function(x, fallback_names = NULL) {
  if (is.null(x) || !length(x)) return(NA_character_)
  x <- as.numeric(x)
  nm <- names(x)
  if ((is.null(nm) || any(!nzchar(nm))) && !is.null(fallback_names) &&
      length(fallback_names) == length(x)) nm <- fallback_names
  values <- sprintf("%.17g", x)
  if (!is.null(nm) && length(nm) == length(values)) values <- paste0(nm, "=", values)
  paste(values, collapse = ";")
}

collapse_messages <- function(x) {
  if (is.null(x) || !length(x) || all(is.na(x))) return(NA_character_)
  clean_text(paste(as.character(x), collapse = " | "))
}

count_value <- function(fit, which) {
  if (is.null(fit) || is.null(fit$counts) || !which %in% names(fit$counts)) return(NA_integer_)
  value <- unname(fit$counts[[which]])
  if (!length(value) || is.na(value)) NA_integer_ else as.integer(value)
}

same_number <- function(x, target, tolerance = 1e-12) {
  length(x) > 0L && all(is.finite(x)) && all(abs(x - target) <= tolerance)
}

round1 <- readRDS(round1_file)
profile <- readRDS(profile_file)

required_round1 <- c(
  "analysis_signature", "signature_components", "selected_fits", "all_results",
  "model_table", "tree", "env", "total_time", "sampling_fraction", "df",
  "cond", "input_sha256", "code_sha256"
)
if (!all(required_round1 %in% names(round1))) {
  stop("Round-one RDS lacks fields: ", paste(setdiff(required_round1, names(round1)), collapse = ", "))
}
required_profile <- c(
  "analysis_signature", "signature_components", "model", "parameterization",
  "diagnosis", "diagnosis_table", "full_runs", "full_best_i",
  "profile_payloads", "profile_summary", "input_sha256", "code_sha256"
)
if (!all(required_profile %in% names(profile))) {
  stop("Profile RDS lacks fields: ", paste(setdiff(required_profile, names(profile)), collapse = ", "))
}

model_names <- c(
  "BCST", "BCSTDCST", "BTimeVar", "BTimeVarDCST", "BCSTDTimeVar",
  "BTimeVarDTimeVar", "BTempVar", "BTempVarDCST", "BCSTDTempVar",
  "BTempVarDTempVar"
)
boundary_model <- "BCSTDTempVar"
regular_models <- setdiff(model_names, boundary_model)
if (!setequal(names(round1$all_results), model_names) ||
    !setequal(round1$model_table$model, model_names) || nrow(round1$model_table) != 10L) {
  stop("Round-one RDS does not contain exactly the ten prespecified models")
}
if (!identical(profile$model, boundary_model) || !identical(profile$diagnosis, "boundary_nonidentifiable")) {
  stop("Completed 03c profile did not classify BCSTDTempVar as boundary_nonidentifiable")
}
if (!identical(profile$parameterization$likelihood_changed, FALSE)) {
  stop("Profile diagnostic is not certified as using the unchanged likelihood")
}

profile_signature_sha <- sha256(profile_signature_file)
profile_signature_record <- readLines(profile_signature_hash_file, warn = FALSE, encoding = "UTF-8")
if (length(profile_signature_record) != 1L ||
    strsplit(trimws(profile_signature_record), "[[:space:]]+")[[1L]][[1L]] != profile_signature_sha ||
    !identical(profile$analysis_signature, profile_signature_sha)) {
  stop("Profile signature chain failed")
}

round1_sha <- sha256(round1_file)
profile_sha <- sha256(profile_file)
tree_sha <- sha256(tree_file)
environment_sha <- sha256(environment_file)
if (!identical(unname(round1$input_sha256[c("tree", "environment")]),
               unname(c(tree_sha, environment_sha)))) {
  stop("Round-one tree/environment SHA-256 mismatch")
}
if (!identical(unname(profile$input_sha256[c("tree", "environment", "round1_rds")]),
               unname(c(tree_sha, environment_sha, round1_sha)))) {
  stop("Profile input SHA-256 chain does not match current formal inputs")
}

diagnosis <- profile$diagnosis_table
if (!is.data.frame(diagnosis) || nrow(diagnosis) != 1L) stop("Invalid profile diagnosis table")
required_diagnosis <- c(
  "model", "diagnosis", "finite_MLE_accepted", "boundary_nonidentifiable",
  "full_best_logLik", "full_best_beta", "full_converged_reproduction_count",
  "most_negative_profile_beta", "most_negative_profile_logLik",
  "all_profile_points_converged", "profile_has_clear_interior_peak",
  "best_at_negative_grid_edge", "negative_edge_within_1e_4_of_profile_best",
  "last_four_profile_values_non_decreasing", "analysis_signature"
)
if (!all(required_diagnosis %in% names(diagnosis))) {
  stop("Profile diagnosis lacks fields: ", paste(setdiff(required_diagnosis, names(diagnosis)), collapse = ", "))
}
boundary_evidence_complete <-
  identical(as.character(diagnosis$model), boundary_model) &&
  identical(as.character(diagnosis$diagnosis), "boundary_nonidentifiable") &&
  identical(as.logical(diagnosis$finite_MLE_accepted), FALSE) &&
  identical(as.logical(diagnosis$boundary_nonidentifiable), TRUE) &&
  identical(as.logical(diagnosis$all_profile_points_converged), TRUE) &&
  identical(as.logical(diagnosis$profile_has_clear_interior_peak), FALSE) &&
  identical(as.logical(diagnosis$best_at_negative_grid_edge), TRUE) &&
  identical(as.logical(diagnosis$negative_edge_within_1e_4_of_profile_best), TRUE) &&
  identical(as.logical(diagnosis$last_four_profile_values_non_decreasing), TRUE) &&
  identical(as.character(diagnosis$analysis_signature), profile$analysis_signature)
if (!boundary_evidence_complete) stop("Profile boundary evidence is incomplete or internally inconsistent")

full_runs <- profile$full_runs
if (!is.list(full_runs) || length(full_runs) != 10L) stop("Expected exactly ten free log-centred profile runs")
full_ll <- vapply(full_runs, function(x) {
  if (is.null(x$fit) || !is.finite(x$fit$LH)) NA_real_ else x$fit$LH
}, numeric(1))
if (!any(is.finite(full_ll))) stop("No finite diagnostic free-run likelihood")
diagnostic_supremum <- max(full_ll, na.rm = TRUE)
diagnostic_best_i <- which.max(full_ll)
diagnostic_best <- full_runs[[diagnostic_best_i]]
canonical_native_start <- function(x) paste(sprintf("%.10g", as.numeric(x)), collapse = ";")
diagnostic_reproducing <- vapply(full_runs, function(x) {
  !is.null(x$fit) && is.finite(x$fit$LH) && x$fit$convergence == 0L &&
    abs(x$fit$LH - diagnostic_supremum) <= 1e-4
}, logical(1))
diagnostic_replication <- length(unique(vapply(
  full_runs[diagnostic_reproducing], function(x) canonical_native_start(x$native_start), character(1)
)))
if (diagnostic_replication < 2L ||
    abs(diagnostic_supremum - as.numeric(diagnosis$full_best_logLik)) > 1e-8 ||
    diagnostic_replication != as.integer(diagnosis$full_converged_reproduction_count)) {
  stop("Free log-centred diagnostic likelihood/replication does not match diagnosis")
}

profile_summary <- profile$profile_summary
if (!is.data.frame(profile_summary) || !all(c("beta", "logLik") %in% names(profile_summary))) {
  stop("Invalid fixed-beta profile summary")
}
edge_i <- which.min(profile_summary$beta)
profile_edge_beta <- as.numeric(profile_summary$beta[edge_i])
profile_edge_logLik <- as.numeric(profile_summary$logLik[edge_i])
if (abs(profile_edge_beta - as.numeric(diagnosis$most_negative_profile_beta)) > 1e-12 ||
    abs(profile_edge_logLik - as.numeric(diagnosis$most_negative_profile_logLik)) > 1e-8) {
  stop("Fixed-beta profile edge does not match diagnosis table")
}

tree <- round1$tree
if (!inherits(tree, "phylo") || Ntip(tree) != 417L || tree$Nnode != 416L ||
    !is.rooted(tree) || !is.binary(tree) ||
    any(!is.finite(tree$edge.length)) || any(tree$edge.length <= 0)) {
  stop("Round-one formal tree failed structural checks")
}
n_obs <- Ntip(tree)
sampling_fraction <- 417 / 2100
conditioning <- "crown"
environment_df <- 80L
total_time <- max(node.depth.edgelength(tree)[seq_len(Ntip(tree))])
if (!same_number(round1$sampling_fraction, sampling_fraction) ||
    !identical(round1$cond, conditioning) || as.integer(round1$df) != environment_df ||
    abs(round1$total_time - total_time) > 1e-12) {
  stop("Round-one fixed settings do not match the formal specification")
}

old_table <- clean_frame(round1$model_table)
old_table <- old_table[match(model_names, old_table$model), , drop = FALSE]
regular_best_count <- setNames(integer(length(regular_models)), regular_models)
regular_selected_runs <- vector("list", length(regular_models))
names(regular_selected_runs) <- regular_models
for (model in regular_models) {
  result <- round1$all_results[[model]]
  if (!is.list(result) || is.na(result$best_i) || result$best_i < 1L ||
      result$best_i > length(result$runs)) stop("Invalid selected run for ", model)
  selected <- result$runs[[result$best_i]]
  fit <- selected$fit
  rate <- selected$rate_qa
  if (is.null(fit) || !is.finite(fit$LH) || fit$convergence != 0L ||
      is.null(rate) || !isTRUE(rate$finite) || !isTRUE(rate$nonnegative) ||
      isTRUE(rate$explosive)) stop("Regular selected fit failed diagnostics: ", model)
  candidates <- vapply(result$runs, function(x) {
    !is.null(x$fit) && is.finite(x$fit$LH) && x$fit$convergence == 0L &&
      abs(x$fit$LH - fit$LH) <= 1e-4
  }, logical(1))
  regular_best_count[[model]] <- length(unique(vapply(
    result$runs[candidates], function(x) canonical_native_start(x$start), character(1)
  )))
  if (regular_best_count[[model]] < 2L) stop("Selected regular likelihood was not reproduced twice: ", model)
  regular_selected_runs[[model]] <- selected
}

# Recalculate inferential quantities only for the nine regular finite-MLE fits.
model_table <- old_table
model_table$model_status <- "finite_MLE"
model_table$eligible_for_AICc <- TRUE
model_table$diagnostic_supremum_logLik <- NA_real_
model_table$diagnostic_supremum_source <- NA_character_
model_table$diagnostic_naive_AICc <- NA_real_
model_table$diagnostic_replication_count <- NA_integer_
model_table$profile_edge_beta <- NA_real_
model_table$profile_edge_logLik <- NA_real_
model_table$diagnostic_full_beta <- NA_real_
model_table$diagnostic_max_mu_per_Ma <- NA_real_

for (model in regular_models) {
  i <- match(model, model_table$model)
  selected <- regular_selected_runs[[model]]
  model_table$logLik[i] <- selected$fit$LH
  model_table$convergence[i] <- selected$fit$convergence
  model_table$best_start_count[i] <- regular_best_count[[model]]
  model_table$AIC[i] <- -2 * model_table$logLik[i] + 2 * model_table$k[i]
  model_table$AICc[i] <- model_table$AIC[i] +
    2 * model_table$k[i] * (model_table$k[i] + 1) /
    (n_obs - model_table$k[i] - 1)
  model_table$notes[i] <- paste0(
    model_table$notes[i],
    "; eligible finite MLE in boundary-aware nine-model AICc normalization"
  )
}

boundary_i <- match(boundary_model, model_table$model)
inferential_numeric <- c(
  "logLik", "AIC", "AICc", "deltaAICc", "Akaike_weight",
  "lambda0", "alpha", "mu0", "beta", "convergence", "best_start_count"
)
for (field in inferential_numeric) model_table[[field]][boundary_i] <- NA
model_table$model_status[boundary_i] <- "boundary_nonidentifiable"
model_table$eligible_for_AICc[boundary_i] <- FALSE
model_table$boundary_flag[boundary_i] <- "boundary_nonidentifiable_beta_to_minus_infinity"
model_table$diagnostic_supremum_logLik[boundary_i] <- diagnostic_supremum
model_table$diagnostic_supremum_source[boundary_i] <- "03c_free_logcenter_boundary_sequence"
diagnostic_naive_aic <- -2 * diagnostic_supremum + 2 * model_table$k[boundary_i]
diagnostic_naive_aicc <- diagnostic_naive_aic +
  2 * model_table$k[boundary_i] * (model_table$k[boundary_i] + 1) /
  (n_obs - model_table$k[boundary_i] - 1)
model_table$diagnostic_naive_AICc[boundary_i] <- diagnostic_naive_aicc
model_table$diagnostic_replication_count[boundary_i] <- diagnostic_replication
model_table$profile_edge_beta[boundary_i] <- profile_edge_beta
model_table$profile_edge_logLik[boundary_i] <- profile_edge_logLik
model_table$diagnostic_full_beta[boundary_i] <- as.numeric(diagnostic_best$native_final["beta"])
model_table$diagnostic_max_mu_per_Ma[boundary_i] <- diagnostic_best$rate_qa$max_mu
model_table$notes[boundary_i] <- paste(
  "Attempted as prespecified; same-likelihood log-centred profile is boundary_nonidentifiable.",
  "No finite MLE, ordinary AICc, delta AICc, Akaike weight, or finite parameter estimate is reported.",
  "diagnostic_naive_AICc is shown only to expose why treating the extreme numerical endpoint as regular would be invalid; it is excluded from ranking and weights."
)

eligible <- which(model_table$eligible_for_AICc)
model_table$deltaAICc <- NA_real_
model_table$Akaike_weight <- NA_real_
model_table$deltaAICc[eligible] <- model_table$AICc[eligible] - min(model_table$AICc[eligible])
relative <- exp(-0.5 * model_table$deltaAICc[eligible])
model_table$Akaike_weight[eligible] <- relative / sum(relative)
model_table$rank <- NA_integer_
eligible_order <- eligible[order(model_table$AICc[eligible], model_table$model[eligible])]
model_table$rank[eligible_order] <- seq_along(eligible_order)
model_table <- model_table[c(eligible_order, boundary_i), , drop = FALSE]
rownames(model_table) <- NULL

front <- c(
  "rank", "model", "model_status", "eligible_for_AICc", "family",
  "lambda_formula", "mu_formula", "k", "logLik", "AIC", "AICc",
  "deltaAICc", "Akaike_weight", "lambda0", "alpha", "mu0", "beta",
  "convergence", "boundary_flag", "best_start_count",
  "diagnostic_supremum_logLik", "diagnostic_supremum_source",
  "diagnostic_naive_AICc", "diagnostic_replication_count",
  "profile_edge_beta", "profile_edge_logLik", "diagnostic_full_beta",
  "diagnostic_max_mu_per_Ma", "notes"
)
model_table <- clean_frame(model_table[, front, drop = FALSE])

formatted <- model_table
format_fixed <- function(x, digits) ifelse(is.finite(x), formatC(x, format = "f", digits = digits), NA_character_)
format_sig <- function(x, digits = 7L) ifelse(is.finite(x), formatC(x, format = "g", digits = digits), NA_character_)
for (field in c("logLik", "AIC", "AICc", "deltaAICc", "diagnostic_supremum_logLik",
                "diagnostic_naive_AICc", "profile_edge_logLik")) {
  formatted[[field]] <- format_fixed(as.numeric(formatted[[field]]), 3L)
}
formatted$Akaike_weight <- format_fixed(as.numeric(formatted$Akaike_weight), 4L)
for (field in c("lambda0", "alpha", "mu0", "beta", "profile_edge_beta",
                "diagnostic_full_beta", "diagnostic_max_mu_per_Ma")) {
  formatted[[field]] <- format_sig(as.numeric(formatted[[field]]))
}

# Build one clean optimization ledger: all original starts for the nine regular
# models plus the ten free log-centred BCSTDTempVar diagnostic runs.
regular_run_row <- function(x, spec, selected) {
  fit <- x$fit
  rate <- x$rate_qa
  parameter_names <- spec$par_names
  data.frame(
    model = x$model, model_status = "finite_MLE", eligible_for_AICc = TRUE,
    run_source = "first_round_v2", start_id = x$start_id,
    start_label = sprintf("first_round_start%02d", x$start_id),
    parameterization = "raw_original", seed = x$seed, method = x$method,
    maxit = 500L, reltol = sqrt(.Machine$double.eps),
    parscale = paste(rep(1, length(x$start)), collapse = ";"),
    starting_values_native = fmt_values(x$start, parameter_names),
    starting_values_raw = fmt_values(x$start, parameter_names),
    final_values_native = if (is.null(fit)) NA_character_ else fmt_values(fit$optim_par, parameter_names),
    final_values_raw = if (is.null(fit)) NA_character_ else fmt_values(fit$optim_par, parameter_names),
    logLik = if (is.null(fit)) NA_real_ else fit$LH,
    AICc_from_vendor = if (is.null(fit)) NA_real_ else fit$aicc,
    diagnostic_naive_AICc = NA_real_,
    convergence = if (is.null(fit)) NA_integer_ else fit$convergence,
    optimizer_message = if (is.null(fit)) NA_character_ else fit$message,
    n_function = count_value(fit, "function"), n_gradient = count_value(fit, "gradient"),
    elapsed_sec = x$elapsed_sec, warnings = collapse_messages(x$warnings), error = collapse_messages(x$error),
    rates_finite = !is.null(rate) && isTRUE(rate$finite),
    rates_nonnegative = !is.null(rate) && isTRUE(rate$nonnegative),
    rate_explosion = if (is.null(rate)) NA else isTRUE(rate$explosive),
    min_lambda = if (is.null(rate)) NA_real_ else rate$min_lambda,
    max_lambda = if (is.null(rate)) NA_real_ else rate$max_lambda,
    min_mu = if (is.null(rate)) NA_real_ else rate$min_mu,
    max_mu = if (is.null(rate)) NA_real_ else rate$max_mu,
    boundary_flag = if (is.null(rate)) "fit_failed" else rate$boundary_flag,
    selected_for_formal_MLE = selected, diagnostic_only = FALSE,
    diagnostic_sequence_best = FALSE,
    total_time = x$total_time, sampling_fraction = x$sampling_fraction,
    conditioning = x$conditioning, environment_df = x$environment_df,
    source_analysis_signature = x$analysis_signature, run_signature = x$run_signature,
    tree_sha256 = x$tree_sha256, environment_sha256 = x$env_sha256,
    source_driver_sha256 = x$driver_sha256, stringsAsFactors = FALSE
  )
}

logcenter_to_original <- function(native, temperature_reference) {
  native <- as.numeric(native)
  names(native) <- c("log_lambda", "log_mu_ref", "beta")
  log_mu0 <- native["log_mu_ref"] - native["beta"] * temperature_reference
  mu0 <- if (is.finite(log_mu0) && log_mu0 <= log(.Machine$double.xmax)) exp(log_mu0) else Inf
  c(lambda0 = exp(native["log_lambda"]), mu0 = mu0,
    beta = native["beta"], log_mu_ref = native["log_mu_ref"], log_mu0 = log_mu0)
}

profile_run_row <- function(x) {
  fit <- x$fit
  rate <- x$rate_qa
  raw_start <- logcenter_to_original(x$native_start, profile$parameterization$temperature_reference)
  raw_final <- x$original_final
  data.frame(
    model = boundary_model, model_status = "boundary_nonidentifiable",
    eligible_for_AICc = FALSE, run_source = "03c_free_logcenter_boundary_sequence",
    start_id = x$start_id, start_label = x$label,
    parameterization = "log_lambda_log_mu_ref_beta", seed = x$seed,
    method = if (is.null(fit)) "Nelder-Mead" else fit$method,
    maxit = x$control$maxit, reltol = x$control$reltol,
    parscale = fmt_values(x$control$parscale),
    starting_values_native = fmt_values(x$native_start),
    starting_values_raw = fmt_values(raw_start),
    final_values_native = if (is.null(fit)) NA_character_ else fmt_values(x$native_final),
    final_values_raw = if (is.null(fit)) NA_character_ else fmt_values(raw_final),
    logLik = if (is.null(fit)) NA_real_ else fit$LH,
    AICc_from_vendor = NA_real_,
    diagnostic_naive_AICc = if (is.null(fit)) NA_real_ else fit$aicc,
    convergence = if (is.null(fit)) NA_integer_ else fit$convergence,
    optimizer_message = if (is.null(fit)) NA_character_ else fit$message,
    n_function = count_value(fit, "function"), n_gradient = count_value(fit, "gradient"),
    elapsed_sec = x$elapsed_sec, warnings = collapse_messages(x$warnings), error = collapse_messages(x$error),
    rates_finite = !is.null(rate) && isTRUE(rate$finite),
    rates_nonnegative = !is.null(rate) && isTRUE(rate$nonnegative),
    rate_explosion = if (is.null(rate)) NA else isTRUE(rate$explosive),
    min_lambda = if (is.null(rate)) NA_real_ else rate$min_lambda,
    max_lambda = if (is.null(rate)) NA_real_ else rate$max_lambda,
    min_mu = if (is.null(rate)) NA_real_ else rate$min_mu,
    max_mu = if (is.null(rate)) NA_real_ else rate$max_mu,
    boundary_flag = if (is.null(rate)) "fit_failed" else rate$boundary_flag,
    selected_for_formal_MLE = FALSE, diagnostic_only = TRUE,
    diagnostic_sequence_best = !is.null(fit) && is.finite(fit$LH) &&
      abs(fit$LH - diagnostic_supremum) <= 1e-4,
    total_time = x$total_time, sampling_fraction = x$sampling_fraction,
    conditioning = x$conditioning, environment_df = x$environment_df,
    source_analysis_signature = x$analysis_signature, run_signature = x$run_signature,
    tree_sha256 = x$tree_sha256, environment_sha256 = x$environment_sha256,
    source_driver_sha256 = unname(profile$code_sha256["driver"]), stringsAsFactors = FALSE
  )
}

run_rows <- list()
for (model in regular_models) {
  result <- round1$all_results[[model]]
  for (i in seq_along(result$runs)) {
    run_rows[[length(run_rows) + 1L]] <- regular_run_row(
      result$runs[[i]], result$spec, i == result$best_i
    )
  }
}
for (run in full_runs) run_rows[[length(run_rows) + 1L]] <- profile_run_row(run)
optimization_runs <- clean_frame(do.call(rbind, run_rows))

all_round_runs <- unlist(lapply(round1$all_results, `[[`, "runs"), recursive = FALSE)
same_tree_all <- all(vapply(all_round_runs, function(x) identical(x$tree_sha256, tree_sha), logical(1))) &&
  all(vapply(full_runs, function(x) identical(x$tree_sha256, tree_sha), logical(1)))
same_environment_all <- all(vapply(all_round_runs, function(x) identical(x$env_sha256, environment_sha), logical(1))) &&
  all(vapply(full_runs, function(x) identical(x$environment_sha256, environment_sha), logical(1)))
same_f_all <- same_number(vapply(all_round_runs, `[[`, numeric(1), "sampling_fraction"), sampling_fraction) &&
  same_number(vapply(full_runs, `[[`, numeric(1), "sampling_fraction"), sampling_fraction)
same_time_all <- same_number(vapply(all_round_runs, `[[`, numeric(1), "total_time"), total_time) &&
  same_number(vapply(full_runs, `[[`, numeric(1), "total_time"), total_time)
same_df_all <- all(vapply(all_round_runs, `[[`, integer(1), "environment_df") == environment_df) &&
  all(vapply(full_runs, `[[`, integer(1), "environment_df") == environment_df)
same_cond_all <- all(vapply(all_round_runs, `[[`, character(1), "conditioning") == conditioning) &&
  all(vapply(full_runs, `[[`, character(1), "conditioning") == conditioning)

driver_sha <- sha256(script_file)
source_hashes <- vapply(required_files, sha256, character(1))
signature_components <- data.frame(
  component = c(
    "signature_schema", "overall_status", "driver_sha256", "round1_rds_sha256",
    "round1_analysis_signature", "profile_rds_sha256", "profile_analysis_signature",
    "tree_sha256", "environment_sha256", "sampling_fraction", "conditioning",
    "environment_df", "total_time", "n_observed", "prespecified_models",
    paste0("source_sha256:", basename(required_files))
  ),
  value = c(
    "03e_finalize_boundary_aware_v4_signature_1",
    "complete_with_one_prespecified_nonregular_boundary", driver_sha, round1_sha,
    round1$analysis_signature, profile_sha, profile$analysis_signature,
    tree_sha, environment_sha, sprintf("%.17g", sampling_fraction), conditioning,
    as.character(environment_df), sprintf("%.17g", total_time), as.character(n_obs),
    paste(model_names, collapse = ";"), unname(source_hashes)
  ),
  stringsAsFactors = FALSE
)
signature_components <- clean_frame(signature_components)
signature_file <- file.path(final_model_dir, "analysis_signature_v4_boundary_aware_final.tsv")
write_clean_tsv(signature_components, signature_file)
analysis_signature <- sha256(signature_file)
signature_hash_file <- file.path(final_model_dir, "analysis_signature_v4_boundary_aware_final.sha256")
write_text_no_overwrite(
  paste(analysis_signature, basename(signature_file)), signature_hash_file
)
optimization_runs$analysis_signature <- analysis_signature

# Preserve the completed profile evidence inside the final table directory.
profile_copy_map <- data.frame(
  source = c(profile_diagnosis_file, profile_summary_file, profile_runs_file, profile_full_runs_file),
  destination = file.path(final_table_dir, c(
    "BCSTDTempVar_boundary_diagnosis.tsv",
    "BCSTDTempVar_fixed_beta_profile_summary.tsv",
    "BCSTDTempVar_fixed_beta_profile_runs.tsv",
    "BCSTDTempVar_full_logcenter_multistart_runs.tsv"
  )),
  stringsAsFactors = FALSE
)
for (i in seq_len(nrow(profile_copy_map))) {
  source <- profile_copy_map$source[i]
  destination <- profile_copy_map$destination[i]
  if (file.exists(destination)) {
    if (!identical(sha256(source), sha256(destination))) {
      stop("Existing copied profile table differs from source: ", destination)
    }
  } else if (!file.copy(source, destination, overwrite = FALSE, copy.mode = TRUE, copy.date = TRUE)) {
    stop("Could not copy profile evidence: ", source)
  }
}
profile_copy_map$source_relative <- substring(profile_copy_map$source, nchar(run_dir) + 2L)
profile_copy_map$destination_relative <- substring(profile_copy_map$destination, nchar(run_dir) + 2L)
profile_copy_map$source_sha256 <- vapply(profile_copy_map$source, sha256, character(1))
profile_copy_map$destination_sha256 <- vapply(profile_copy_map$destination, sha256, character(1))
profile_copy_map$byte_identical <- profile_copy_map$source_sha256 == profile_copy_map$destination_sha256
profile_copy_map$analysis_signature <- analysis_signature
profile_copy_table <- profile_copy_map[, c(
  "source_relative", "destination_relative", "source_sha256",
  "destination_sha256", "byte_identical", "analysis_signature"
)]

regular_rows <- model_table$model_status == "finite_MLE"
boundary_rows <- model_table$model_status == "boundary_nonidentifiable"
selected_ledger <- optimization_runs$selected_for_formal_MLE
regular_selected_ledger <- selected_ledger & optimization_runs$model_status == "finite_MLE"
regular_aicc_check <- -2 * model_table$logLik[regular_rows] + 2 * model_table$k[regular_rows] +
  2 * model_table$k[regular_rows] * (model_table$k[regular_rows] + 1) /
  (n_obs - model_table$k[regular_rows] - 1)
boundary_inferential_na <- all(vapply(
  model_table[boundary_rows, inferential_numeric, drop = FALSE],
  function(x) all(is.na(x)), logical(1)
))

qa <- data.frame(
  check = c(
    "prespecified_models_attempted", "final_table_rows", "finite_regular_MLE_models",
    "boundary_nonidentifiable_models", "ten_models_have_finite_MLE",
    "boundary_model_is_BCSTDTempVar", "boundary_not_eligible_for_AICc",
    "boundary_inferential_fields_are_NA", "boundary_diagnostic_supremum_not_ranked",
    "boundary_profile_evidence_complete", "boundary_free_logcenter_runs",
    "boundary_diagnostic_sequence_reproduced_twice_1e4",
    "boundary_profile_edge_is_most_negative", "boundary_extinction_rate_explosion_detected",
    "eligible_Akaike_weights_sum", "boundary_Akaike_weight_is_NA",
    "regular_AICc_matches_full_precision_recalculation",
    "regular_selected_convergence_zero", "regular_selected_best_reproduced_twice_1e4",
    "regular_selected_rates_finite", "regular_selected_rates_nonnegative",
    "regular_selected_no_rate_explosion_gt_100_per_Ma",
    "same_tree_sha256_all_attempts", "same_environment_sha256_all_attempts",
    "same_sampling_fraction_all_attempts", "same_conditioning_all_attempts",
    "same_temperature_df_all_attempts", "same_total_time_all_attempts",
    "profile_tables_copied_byte_identical", "model_table_TSV_roundtrip",
    "optimization_ledger_TSV_roundtrip", "no_control_characters_in_final_tables",
    "overall_status"
  ),
  observed = c(
    length(unique(c(names(round1$all_results), profile$model))), nrow(model_table),
    sum(regular_rows), sum(boundary_rows), FALSE,
    identical(model_table$model[boundary_rows], boundary_model),
    all(!model_table$eligible_for_AICc[boundary_rows]), boundary_inferential_na,
    all(is.na(model_table$rank[boundary_rows])) && all(is.na(model_table$Akaike_weight[boundary_rows])),
    boundary_evidence_complete, length(full_runs), diagnostic_replication,
    profile_edge_beta == min(profile_summary$beta),
    isTRUE(diagnostic_best$rate_qa$explosive) && diagnostic_best$rate_qa$max_mu > 100,
    format(sum(model_table$Akaike_weight[regular_rows]), digits = 17),
    all(is.na(model_table$Akaike_weight[boundary_rows])),
    max(abs(model_table$AICc[regular_rows] - regular_aicc_check)),
    all(optimization_runs$convergence[regular_selected_ledger] == 0L),
    min(model_table$best_start_count[regular_rows]),
    all(optimization_runs$rates_finite[regular_selected_ledger]),
    all(optimization_runs$rates_nonnegative[regular_selected_ledger]),
    all(!optimization_runs$rate_explosion[regular_selected_ledger]),
    same_tree_all, same_environment_all, same_f_all, same_cond_all,
    same_df_all, same_time_all, all(profile_copy_map$byte_identical),
    preflight_tsv(model_table, final_table_dir),
    preflight_tsv(optimization_runs, final_model_dir),
    !has_control_text(model_table) && !has_control_text(optimization_runs) && !has_control_text(profile_copy_table),
    "complete_with_one_prespecified_nonregular_boundary"
  ),
  expected = c(
    "10", "10", "9", "1", "FALSE", "TRUE", "TRUE", "TRUE", "TRUE",
    "TRUE", "10", ">=2", "TRUE", "TRUE", "1 within floating-point error",
    "TRUE", "<=1e-12", "TRUE", ">=2", "TRUE", "TRUE", "TRUE",
    "TRUE", "TRUE", "TRUE", "TRUE", "TRUE", "TRUE", "TRUE", "TRUE",
    "TRUE", "TRUE", "complete_with_one_prespecified_nonregular_boundary"
  ),
  pass = c(
    length(unique(c(names(round1$all_results), profile$model))) == 10L,
    nrow(model_table) == 10L, sum(regular_rows) == 9L, sum(boundary_rows) == 1L,
    TRUE, identical(model_table$model[boundary_rows], boundary_model),
    all(!model_table$eligible_for_AICc[boundary_rows]), boundary_inferential_na,
    all(is.na(model_table$rank[boundary_rows])) && all(is.na(model_table$Akaike_weight[boundary_rows])),
    boundary_evidence_complete, length(full_runs) == 10L, diagnostic_replication >= 2L,
    profile_edge_beta == min(profile_summary$beta),
    isTRUE(diagnostic_best$rate_qa$explosive) && diagnostic_best$rate_qa$max_mu > 100,
    abs(sum(model_table$Akaike_weight[regular_rows]) - 1) < 1e-12,
    all(is.na(model_table$Akaike_weight[boundary_rows])),
    max(abs(model_table$AICc[regular_rows] - regular_aicc_check)) <= 1e-12,
    all(optimization_runs$convergence[regular_selected_ledger] == 0L),
    min(model_table$best_start_count[regular_rows]) >= 2L,
    all(optimization_runs$rates_finite[regular_selected_ledger]),
    all(optimization_runs$rates_nonnegative[regular_selected_ledger]),
    all(!optimization_runs$rate_explosion[regular_selected_ledger]),
    same_tree_all, same_environment_all, same_f_all, same_cond_all,
    same_df_all, same_time_all, all(profile_copy_map$byte_identical),
    preflight_tsv(model_table, final_table_dir),
    preflight_tsv(optimization_runs, final_model_dir),
    !has_control_text(model_table) && !has_control_text(optimization_runs) && !has_control_text(profile_copy_table),
    TRUE
  ),
  stringsAsFactors = FALSE
)
qa <- clean_frame(qa)
if (!all(qa$pass)) stop("Boundary-aware final QA failed: ", paste(qa$check[!qa$pass], collapse = ", "))
if (!preflight_tsv(qa, final_table_dir)) stop("QA table failed TSV preflight")

write_clean_tsv(optimization_runs, file.path(final_model_dir, "optimization_runs.tsv"))
write_clean_tsv(model_table, file.path(final_table_dir, "model_table_10models_full_precision.tsv"))
write_clean_tsv(formatted, file.path(final_table_dir, "model_table_10models_formatted.tsv"))
write_clean_tsv(profile_copy_table, file.path(final_table_dir, "BCSTDTempVar_profile_evidence_provenance.tsv"))
write_clean_tsv(qa, file.path(final_table_dir, "model_selection_QA.tsv"))

best <- model_table[which.min(model_table$AICc), , drop = FALSE]
runner <- model_table[order(model_table$AICc)[2L], , drop = FALSE]
report <- c(
  "# Boundary-aware finalization report",
  "",
  "## Formal status",
  "",
  "`complete_with_one_prespecified_nonregular_boundary`",
  "",
  paste0(
    "All ten prespecified models were attempted. Nine models yielded reproducible finite regular MLEs and are eligible for ordinary AICc comparison. `BCSTDTempVar` was retained as the tenth row but classified `boundary_nonidentifiable`; it is not silently dropped."
  ),
  "",
  "## Why BCSTDTempVar is not assigned ordinary AICc or weight",
  "",
  paste0(
    "Under the unchanged likelihood, ten independent free log-centred runs converged to an indistinguishable extreme sequence (largest evaluated log likelihood ",
    format(diagnostic_supremum, digits = 15), "; beta ",
    format(as.numeric(diagnostic_best$native_final["beta"]), digits = 10),
    "; maximum evaluated extinction rate ",
    format(diagnostic_best$rate_qa$max_mu, digits = 10),
    " per Ma). The fixed-beta profile continued improving to its most-negative evaluated edge (beta = ",
    format(profile_edge_beta, digits = 8), "; log likelihood ",
    format(profile_edge_logLik, digits = 15),
    ") and had no finite interior peak. Thus beta is not finitely identified and the ordinary regular-MLE assumptions needed for its AICc and Akaike weight do not hold."
  ),
  "",
  paste0(
    "For transparency only, applying the ordinary k=3 formula to the largest extreme diagnostic endpoint would give `diagnostic_naive_AICc = ",
    format(diagnostic_naive_aicc, digits = 15),
    "`. This value is invalid for formal ranking, delta AICc, or weight and is never included in their calculation."
  ),
  "",
  "## Eligible comparison",
  "",
  paste0(
    "Among the nine candidates with finite regular MLEs, ", best$model,
    " ranked first (AICc = ", format(best$AICc, digits = 15),
    ", delta AICc = 0, weight = ", format(best$Akaike_weight, digits = 10),
    "). The runner-up was ", runner$model, " (delta AICc = ",
    format(runner$deltaAICc, digits = 10), ", weight = ",
    format(runner$Akaike_weight, digits = 10), ")."
  ),
  "",
  "These weights sum to one only across the nine eligible regular-MLE candidates. They measure relative support inside that eligible subset and are not causal probabilities.",
  "",
  "## Evidence files",
  "",
  "- `model_table_10models_full_precision.tsv`: ten rows, explicit eligibility/status, and diagnostic-only boundary fields.",
  "- `model_selection_QA.tsv`: strict expected-failure-aware QA; `ten_models_have_finite_MLE = FALSE` is an expected, passing fact.",
  "- `BCSTDTempVar_fixed_beta_profile_summary.tsv`: byte-identical copy of the completed profile summary.",
  "- `BCSTDTempVar_boundary_diagnosis.tsv`: byte-identical copy of the signed 03c diagnosis.",
  "- `../03_models/final_v4_boundary_aware/optimization_runs.tsv`: nine-model original multi-start runs plus ten free log-centred boundary-diagnostic runs.",
  "",
  paste0("Final analysis signature: `", analysis_signature, "`.")
)
write_text_no_overwrite(report, file.path(final_table_dir, "boundary_exception_report.md"))

selected_fits <- round1$selected_fits[regular_models]
selected_fits[boundary_model] <- list(list(
  model = boundary_model, status = "boundary_nonidentifiable", fit = NULL,
  finite_MLE = FALSE, eligible_for_AICc = FALSE,
  diagnostic_supremum_logLik = diagnostic_supremum,
  diagnostic_naive_AICc = diagnostic_naive_aicc,
  diagnostic_replication_count = diagnostic_replication,
  profile_edge_beta = profile_edge_beta, profile_edge_logLik = profile_edge_logLik
))
selected_fits <- selected_fits[model_names]
all_results <- round1$all_results[regular_models]
all_results[boundary_model] <- list(list(
  model = boundary_model, status = "boundary_nonidentifiable",
  eligible_for_AICc = FALSE, profile_diagnostic = profile
))
all_results <- all_results[model_names]

final_payload <- list(
  overall_status = "complete_with_one_prespecified_nonregular_boundary",
  analysis_signature = analysis_signature,
  signature_components = signature_components,
  prespecified_models = model_names, regular_finite_MLE_models = regular_models,
  boundary_nonidentifiable_model = boundary_model,
  selected_fits = selected_fits, all_results = all_results,
  model_table = model_table, optimization_runs = optimization_runs,
  model_selection_QA = qa,
  boundary_exception = list(
    diagnosis = diagnosis, diagnostic_supremum_logLik = diagnostic_supremum,
    diagnostic_naive_AICc = diagnostic_naive_aicc,
    diagnostic_replication_count = diagnostic_replication,
    profile_edge_beta = profile_edge_beta, profile_edge_logLik = profile_edge_logLik,
    instruction = "Do not use the BCSTDTempVar diagnostic supremum/naive AICc in formal ranking or Akaike-weight normalization."
  ),
  source_round1 = list(path = round1_file, sha256 = round1_sha,
                       analysis_signature = round1$analysis_signature),
  source_profile = list(path = profile_file, sha256 = profile_sha,
                        analysis_signature = profile$analysis_signature),
  tree = tree, env = round1$env, total_time = total_time,
  sampling_fraction = sampling_fraction, df = environment_df, cond = conditioning,
  input_sha256 = c(tree = tree_sha, environment = environment_sha,
                   round1_rds = round1_sha, profile_rds = profile_sha),
  code_sha256 = c(driver = driver_sha)
)
atomic_save_rds(final_payload, file.path(final_model_dir, "model_fits_10models.rds"))
atomic_save_rds(final_payload$boundary_exception,
                file.path(final_model_dir, "BCSTDTempVar_boundary_exception.rds"))

cat(
  "Boundary-aware v4 finalization complete: 10 attempted; 9 finite-MLE eligible; ",
  "BCSTDTempVar boundary_nonidentifiable and excluded from ordinary AICc weights; ",
  "best eligible model = ", best$model, "; signature = ", analysis_signature, "\n",
  sep = ""
)
