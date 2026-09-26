#!/usr/bin/env Rscript

# Finalize the BCSTDTempVar dt-sensitivity diagnostic from immutable v5 and v6
# RDS/checkpoints.  No likelihood is evaluated here.  v7 exists because the v6
# fits completed successfully but its TSV guard used a non-PCRE control-character
# test.  This script sanitizes actual CR/LF/TAB characters and writes new outputs.

suppressPackageStartupMessages({
  library(ape)
  library(pspline)
})
options(stringsAsFactors = FALSE, digits = 17, scipen = 999)

args_full <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", grep("^--file=", args_full, value = TRUE))
if (length(script_arg) != 1L) stop("Cannot determine script path")
script_file <- normalizePath(script_arg, mustWork = TRUE)
run_dir <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)

tree_file <- file.path(run_dir, "01_tree", "T25_deep_only_dated.Hesperiinae_417sp_RPANDA.tre")
env_file <- file.path(run_dir, "02_environment", "temperature_surface_to_tree_age.tsv")
v5_dir <- file.path(run_dir, "03_models", "BCSTDTempVar_dt_sensitivity_v5")
v5_file <- file.path(v5_dir, "BCSTDTempVar_dt_sensitivity.rds")
v5_signature_file <- file.path(v5_dir, "analysis_signature.tsv")
v6_dir <- file.path(run_dir, "03_models", "BCSTDTempVar_dt_sensitivity_rescue_v6")
v6_signature_file <- file.path(v6_dir, "analysis_signature.tsv")
v6_signature_hash_file <- file.path(v6_dir, "analysis_signature.sha256")
v6_payload_files <- list.files(v6_dir, pattern = "present_center_rescue[.]rds$", full.names = TRUE)
out_model_dir <- file.path(run_dir, "03_models", "BCSTDTempVar_dt_sensitivity_final_v7")
out_table_dir <- file.path(run_dir, "04_tables", "BCSTDTempVar_dt_sensitivity_final_v7")
dir.create(out_model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)

required <- c(tree_file, env_file, v5_file, v5_signature_file,
              v6_signature_file, v6_signature_hash_file, v6_payload_files)
if (!length(v6_payload_files) || any(!file.exists(required))) {
  stop("Missing v5/v6 sensitivity evidence:\n", paste(required[!file.exists(required)], collapse = "\n"))
}

sha256 <- function(path) {
  z <- system2("sha256sum", path, stdout = TRUE, stderr = TRUE)
  status <- attr(z, "status")
  if (!length(z) || (!is.null(status) && status != 0L)) stop("sha256sum failed: ", path)
  strsplit(z[1], "[[:space:]]+")[[1]][1]
}
fmt_num <- function(x) paste(sprintf("%.17g", as.numeric(x)), collapse = ";")
fmt_vec <- function(x) {
  if (is.null(x) || !length(x)) return(NA_character_)
  paste(sprintf("%.17g", as.numeric(x)), collapse = ";")
}
clean_text <- function(x) {
  y <- enc2utf8(as.character(x))
  y <- gsub("[\\r\\n\\t]+", " ", y, perl = TRUE)
  trimws(gsub("[ ]{2,}", " ", y, perl = TRUE))
}
clean_frame <- function(x) {
  names(x) <- clean_text(names(x))
  for (nm in names(x)) {
    if (is.character(x[[nm]]) || is.factor(x[[nm]])) x[[nm]] <- clean_text(x[[nm]])
  }
  x
}
has_actual_control <- function(x) {
  any(vapply(x, function(v) {
    (is.character(v) || is.factor(v)) &&
      any(grepl("[\\r\\n\\t]", as.character(v), perl = TRUE), na.rm = TRUE)
  }, logical(1)))
}
write_clean_tsv <- function(x, path) {
  x <- clean_frame(x)
  if (has_actual_control(x)) stop("Actual CR/LF/TAB remains: ", path)
  tmp <- tempfile(pattern = paste0(".", basename(path)), tmpdir = dirname(path), fileext = ".tmp")
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  write.table(x, tmp, sep = "\t", quote = FALSE, row.names = FALSE,
              col.names = TRUE, na = "NA", fileEncoding = "UTF-8")
  raw <- readLines(tmp, warn = FALSE, encoding = "UTF-8")
  back <- read.delim(tmp, sep = "\t", header = TRUE, quote = "", comment.char = "",
                     check.names = FALSE, stringsAsFactors = FALSE,
                     na.strings = c("NA", "NaN", ""), fill = FALSE, fileEncoding = "UTF-8")
  if (length(raw) != nrow(x) + 1L || nrow(back) != nrow(x) || ncol(back) != ncol(x) ||
      !identical(names(back), names(x))) stop("TSV roundtrip failed: ", path)
  if (file.exists(path)) {
    if (!identical(sha256(path), sha256(tmp))) stop("Refusing to overwrite different v7 TSV: ", path)
    return(invisible(path))
  }
  if (!file.rename(tmp, path)) stop("Could not atomically place TSV: ", path)
  invisible(path)
}
atomic_save_rds <- function(object, path) {
  if (file.exists(path)) stop("Refusing to overwrite v7 RDS: ", path)
  tmp <- tempfile(pattern = paste0(".", basename(path)), tmpdir = dirname(path), fileext = ".tmp")
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  saveRDS(object, tmp, version = 3)
  if (!file.rename(tmp, path)) stop("Could not atomically place RDS: ", path)
  invisible(path)
}

tree <- read.tree(tree_file)
if (Ntip(tree) != 417L || tree$Nnode != 416L || !is.rooted(tree) || !is.binary(tree) ||
    any(!is.finite(tree$edge.length)) || any(tree$edge.length <= 0)) stop("Formal tree QA failed")
total_time <- max(node.depth.edgelength(tree)[seq_len(Ntip(tree))])
env_all <- read.delim(env_file, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
env_data <- data.frame(time = as.numeric(env_all[[1]]), temp = as.numeric(env_all[[2]]))
if (anyNA(env_data) || any(!is.finite(as.matrix(env_data)))) stop("Environment QA failed")
env_spline <- sm.spline(env_data$time, env_data$temp, df = 80L)
time_grid <- sort(unique(c(seq(0, total_time, length.out = 20001L),
                           env_data$time[env_data$time >= 0 & env_data$time <= total_time])))
temperature_grid <- as.numeric(predict(env_spline, time_grid))
temperature_ref <- median(temperature_grid)
temperature_min <- min(temperature_grid)

tree_sha <- sha256(tree_file)
env_sha <- sha256(env_file)
v5_sha <- sha256(v5_file)
v5_sig_sha <- sha256(v5_signature_file)
v6_sig_sha <- sha256(v6_signature_file)
v6_payload_sha <- setNames(vapply(v6_payload_files, sha256, character(1)), basename(v6_payload_files))
driver_sha <- sha256(script_file)

v5 <- readRDS(v5_file)
if (!identical(v5$analysis_signature, v5_sig_sha) ||
    !identical(unname(v5$input_sha256[c("tree", "environment")]), unname(c(tree_sha, env_sha))) ||
    !identical(v5$formal_inference_dt, 0.005) || !identical(v5$formal_likelihood_changed, FALSE)) {
  stop("v5 signature/input/formal-status chain failed")
}
v6_sig_record <- strsplit(trimws(readLines(v6_signature_hash_file, warn = FALSE)), "[[:space:]]+")[[1]][1]
if (!identical(v6_sig_record, v6_sig_sha)) stop("v6 signature hash chain failed")
v6_payloads <- lapply(v6_payload_files, readRDS)
if (any(vapply(v6_payloads, function(x) !identical(x$analysis_signature, v6_sig_sha), logical(1)))) {
  stop("v6 payload signature mismatch")
}
v6_runs <- unlist(lapply(v6_payloads, `[[`, "runs"), recursive = FALSE)
if (length(v6_runs) < 4L || any(vapply(v6_runs, function(x) {
  !identical(x$tree_sha256, tree_sha) || !identical(x$environment_sha256, env_sha) ||
    !identical(x$sampling_fraction, 417 / 2100) || !identical(x$conditioning, "crown") ||
    !identical(as.integer(x$environment_df), 80L) || abs(x$total_time - total_time) > 1e-12
}, logical(1)))) stop("v6 run settings/input mismatch")

v6_conversion_error <- max(vapply(v6_runs, function(x) {
  expected <- x$present_native_final[2] + x$present_native_final[3] *
    (temperature_ref - temperature_min)
  abs(unname(x$native_final[2]) - unname(expected))
}, numeric(1)))
if (!is.finite(v6_conversion_error) || v6_conversion_error > 1e-10) {
  stop("v6 present-centred to reference-centred conversion failed")
}

signature_components <- data.frame(
  component = c("signature_schema", "driver_sha256", "tree_sha256", "environment_sha256",
                "v5_rds_sha256", "v5_analysis_signature", "v6_analysis_signature",
                "v6_signature_record_sha256", "sampling_fraction", "conditioning",
                "environment_df", "total_time", "formal_inference_dt", "diagnostic_dt_values",
                "v6_conversion_max_abs_log_scale", "replication_tolerance",
                "scientific_rule", "R_version", paste0("v6_payload_sha256:", names(v6_payload_sha))),
  value = c("03d3_BCSTDTempVar_dt_sensitivity_final_signature_1", driver_sha, tree_sha, env_sha,
            v5_sha, v5$analysis_signature, v6_sig_sha, sha256(v6_signature_hash_file),
            sprintf("%.17g", 417 / 2100), "crown", "80", sprintf("%.17g", total_time),
            sprintf("%.17g", 0.005), fmt_num(sort(v5$diagnostic_dt_values)),
            sprintf("%.17g", v6_conversion_error), sprintf("%.17g", 1e-4),
            paste("inverse-grid scaling requires slope [-1.05,-0.95] and pulse-mass ratio <=1.05;",
                  "likelihood invariance is reported separately at range <=0.001; formal inference remains official dt=0.005"),
            R.version.string, unname(v6_payload_sha)), stringsAsFactors = FALSE
)
signature_path <- file.path(out_model_dir, "analysis_signature.tsv")
write_clean_tsv(signature_components, signature_path)
analysis_signature <- sha256(signature_path)
signature_hash_path <- file.path(out_model_dir, "analysis_signature.sha256")
signature_hash_line <- paste(analysis_signature, basename(signature_path))
if (file.exists(signature_hash_path)) {
  if (!identical(readLines(signature_hash_path, warn = FALSE), signature_hash_line)) stop("Existing v7 signature hash differs")
} else writeLines(signature_hash_line, signature_hash_path, useBytes = TRUE)

run_id <- function(source, x) paste(source, sprintf("%.17g", x$dt), x$start_id, sep = ":")
run_row <- function(x, source) {
  is_v6 <- identical(source, "v6_present_center_rescue")
  native_start <- if (is_v6) x$native_start else x$native_start
  native_final <- if (is_v6) x$present_native_final else x$native_final
  fit <- x$fit
  rate <- x$rate_qa
  cell_width <- total_time / (1L + as.integer(total_time / x$dt))
  data.frame(
    run_id = run_id(source, x), model = "BCSTDTempVar", dt = x$dt,
    source = source, source_parameterization = if (is_v6) "log_lambda|log_mu_present|beta" else "log_lambda|log_mu_ref|beta",
    start_id = x$start_id, label = x$label,
    starting_values_native = fmt_vec(native_start), final_values_native = fmt_vec(native_final),
    final_log_lambda = x$native_final[1], final_log_mu_ref = x$native_final[2],
    final_beta = x$native_final[3], mu_present = x$mu_present,
    logLik = if (is.null(fit)) NA_real_ else fit$LH,
    convergence = if (is.null(fit)) NA_integer_ else fit$convergence,
    optimizer_message = if (is.null(fit)) NA_character_ else fit$message,
    n_function = if (is.null(fit) || is.null(fit$counts)) NA_integer_ else unname(fit$counts["function"]),
    elapsed_sec = x$elapsed_sec,
    warnings = if (length(x$warnings) && !all(is.na(x$warnings))) paste(x$warnings, collapse = " | ") else NA_character_,
    error = if (length(x$error) && !all(is.na(x$error))) paste(x$error, collapse = " | ") else NA_character_,
    rates_finite = !is.null(rate) && isTRUE(rate$finite),
    rates_nonnegative = !is.null(rate) && isTRUE(rate$nonnegative),
    rate_gt_100 = if (is.null(rate)) NA else if (!is.null(rate$rate_gt_100)) rate$rate_gt_100 else rate$max_mu > 100,
    rate_gt_1000 = if (is.null(rate)) NA else if (!is.null(rate$rate_gt_1000)) rate$rate_gt_1000 else rate$max_mu > 1000,
    min_lambda = if (is.null(rate)) NA_real_ else rate$min_lambda,
    max_lambda = if (is.null(rate)) NA_real_ else rate$max_lambda,
    min_mu = if (is.null(rate)) NA_real_ else rate$min_mu,
    max_mu = if (is.null(rate)) NA_real_ else rate$max_mu,
    first_cell_width_total_time = cell_width,
    mu_present_times_first_cell_width_total = x$mu_present * cell_width,
    sampling_fraction = x$sampling_fraction, conditioning = x$conditioning,
    environment_df = x$environment_df, total_time = x$total_time,
    tree_sha256 = x$tree_sha256, environment_sha256 = x$environment_sha256,
    source_analysis_signature = x$analysis_signature,
    run_signature = x$run_signature, analysis_signature = analysis_signature,
    stringsAsFactors = FALSE
  )
}

v5_runs <- unlist(lapply(v5$dt_payloads, `[[`, "runs"), recursive = FALSE)
all_runs <- c(v5_runs, v6_runs)
all_sources <- c(rep("v5_reference_center", length(v5_runs)),
                 rep("v6_present_center_rescue", length(v6_runs)))
runs_table <- clean_frame(do.call(rbind, Map(run_row, all_runs, all_sources)))

eligible_run <- function(x) !is.null(x$fit) && is.finite(x$fit$LH) && x$fit$convergence == 0L
dt_values <- sort(unique(vapply(all_runs, `[[`, numeric(1), "dt")), decreasing = TRUE)
summary_rows <- lapply(dt_values, function(dt) {
  at <- which(vapply(all_runs, function(x) identical(as.numeric(x$dt), dt), logical(1)))
  runs <- all_runs[at]
  sources <- all_sources[at]
  good <- which(vapply(runs, eligible_run, logical(1)))
  best <- if (length(good)) good[which.max(vapply(runs[good], function(x) x$fit$LH, numeric(1)))] else NA_integer_
  chosen <- if (is.na(best)) NULL else runs[[best]]
  best_ll <- if (is.null(chosen)) NA_real_ else chosen$fit$LH
  reproduced <- if (is.null(chosen)) logical(length(runs)) else vapply(runs, function(x) {
    eligible_run(x) && abs(x$fit$LH - best_ll) <= 1e-4
  }, logical(1))
  nonconv_ll <- vapply(runs, function(x) {
    if (!is.null(x$fit) && is.finite(x$fit$LH) && x$fit$convergence != 0L) x$fit$LH else NA_real_
  }, numeric(1))
  max_nonconv <- if (any(is.finite(nonconv_ll))) max(nonconv_ll, na.rm = TRUE) else NA_real_
  cell_width <- total_time / (1L + as.integer(total_time / dt))
  data.frame(
    dt = dt, n_runs = length(runs), n_converged = length(good),
    best_start_count_1e_4 = sum(reproduced), selected_run_id = if (is.null(chosen)) NA_character_ else run_id(sources[best], chosen),
    selected_source = if (is.null(chosen)) NA_character_ else sources[best],
    logLik = best_ll, maximum_finite_nonconverged_logLik = max_nonconv,
    no_higher_finite_nonconverged = is.finite(best_ll) && (!is.finite(max_nonconv) || best_ll + 1e-8 >= max_nonconv),
    beta = if (is.null(chosen)) NA_real_ else chosen$native_final[3],
    lambda = if (is.null(chosen)) NA_real_ else exp(chosen$native_final[1]),
    mu_present = if (is.null(chosen)) NA_real_ else chosen$mu_present,
    min_mu_observed_interval = if (is.null(chosen) || is.null(chosen$rate_qa)) NA_real_ else chosen$rate_qa$min_mu,
    max_mu_observed_interval = if (is.null(chosen) || is.null(chosen$rate_qa)) NA_real_ else chosen$rate_qa$max_mu,
    selected_rates_finite = !is.null(chosen) && !is.null(chosen$rate_qa) && chosen$rate_qa$finite,
    first_cell_width_total_time = cell_width,
    mu_present_times_first_cell_width_total = if (is.null(chosen)) NA_real_ else chosen$mu_present * cell_width,
    analysis_signature = analysis_signature, stringsAsFactors = FALSE
  )
})
summary_table <- clean_frame(do.call(rbind, summary_rows))

strict <- summary_table$n_converged >= 2L & summary_table$best_start_count_1e_4 >= 2L &
  summary_table$no_higher_finite_nonconverged & summary_table$selected_rates_finite
logLik_range <- if (all(is.finite(summary_table$logLik))) diff(range(summary_table$logLik)) else Inf
slope <- if (all(is.finite(summary_table$mu_present)) && all(summary_table$mu_present > 0)) {
  unname(coef(lm(log(mu_present) ~ log(dt), data = summary_table))[2])
} else NA_real_
pulse_ratio <- if (all(is.finite(summary_table$mu_present_times_first_cell_width_total)) &&
                   all(summary_table$mu_present_times_first_cell_width_total > 0)) {
  max(summary_table$mu_present_times_first_cell_width_total) /
    min(summary_table$mu_present_times_first_cell_width_total)
} else Inf
inverse_grid_scaling <- all(strict) && is.finite(slope) && slope >= -1.05 && slope <= -0.95 && pulse_ratio <= 1.05
likelihood_invariant_1e3 <- logLik_range <= 1e-3
diagnosis <- if (!all(strict)) {
  "inconclusive_dt_optimization"
} else if (inverse_grid_scaling && likelihood_invariant_1e3) {
  "endpoint_grid_ridge_with_invariant_likelihood"
} else if (inverse_grid_scaling) {
  "endpoint_grid_scaling_with_residual_logLik_discretization"
} else {
  "dt_sensitivity_without_inverse_grid_scaling"
}
diagnosis_table <- data.frame(
  model = "BCSTDTempVar", diagnosis = diagnosis,
  all_dt_strictly_optimized_and_replicated = all(strict),
  inverse_grid_scaling_supported = inverse_grid_scaling,
  likelihood_invariant_within_1e_3 = likelihood_invariant_1e3,
  logLik_range_across_dt = logLik_range,
  expected_mu_present_loglog_slope = -1,
  observed_mu_present_loglog_slope = slope,
  present_pulse_mass_max_min_ratio = pulse_ratio,
  implication = paste(
    "The present-endpoint rate follows the integration-grid scale; residual log-likelihood",
    "variation across dt is reported rather than hidden. This supports numerical endpoint",
    "grid dependence and does not define a finite biological beta estimate."
  ),
  formal_inference_dt = 0.005,
  formal_inference_statement = paste(
    "Diagnostic only. The signed 03c profile under the archived Toussaint/RPANDA dt=0.005",
    "likelihood remains the formal basis for boundary_nonidentifiable status."
  ),
  analysis_signature = analysis_signature, stringsAsFactors = FALSE
)
diagnosis_table <- clean_frame(diagnosis_table)

qa <- data.frame(
  check = c(
    "v5_signature_chain", "v6_signature_chain", "tree_environment_hash_chain",
    "exact_four_dt_values", "v5_sixteen_runs_retained", "v6_four_rescue_runs_retained",
    "at_least_four_runs_per_dt", "all_dt_have_converged_solution",
    "all_dt_best_reproduced_at_least_twice_1e_4", "no_higher_finite_nonconverged_each_dt",
    "all_selected_rates_finite", "present_to_reference_center_conversion",
    "sampling_fraction_conditioning_df_total_time", "formal_inference_remains_official_dt_0.005",
    "actual_control_characters_removed", "runs_TSV_roundtrip", "summary_TSV_roundtrip",
    "diagnosis_TSV_roundtrip", "signature_embedded_all_rows"
  ),
  observed = c(
    v5$analysis_signature, v6_sig_sha, paste(tree_sha, env_sha, sep = ";"),
    paste(sort(dt_values), collapse = ";"), length(v5_runs), length(v6_runs),
    min(summary_table$n_runs), all(summary_table$n_converged > 0L),
    min(summary_table$best_start_count_1e_4), all(summary_table$no_higher_finite_nonconverged),
    all(summary_table$selected_rates_finite), sprintf("%.17g", v6_conversion_error),
    sprintf("f=%.17g;cond=crown;df=80;time=%.17g", 417 / 2100, total_time),
    TRUE, TRUE, TRUE, TRUE, TRUE,
    all(runs_table$analysis_signature == analysis_signature) && all(summary_table$analysis_signature == analysis_signature)
  ),
  expected = c(v5_sig_sha, v6_sig_record, "current formal tree/environment SHA", "0.001;0.0025;0.005;0.01",
               "16", ">=4", ">=4", "TRUE", ">=2", "TRUE", "TRUE", "<=1e-10",
               "417/2100;crown;80;formal total_time", "TRUE", "TRUE", "TRUE", "TRUE", "TRUE", "TRUE"),
  pass = c(
    identical(v5$analysis_signature, v5_sig_sha), identical(v6_sig_sha, v6_sig_record),
    identical(unname(v5$input_sha256[c("tree", "environment")]), unname(c(tree_sha, env_sha))),
    identical(sort(dt_values), c(0.001, 0.0025, 0.005, 0.01)), length(v5_runs) == 16L,
    length(v6_runs) >= 4L, min(summary_table$n_runs) >= 4L,
    all(summary_table$n_converged > 0L), min(summary_table$best_start_count_1e_4) >= 2L,
    all(summary_table$no_higher_finite_nonconverged), all(summary_table$selected_rates_finite),
    v6_conversion_error <= 1e-10,
    all(runs_table$sampling_fraction == 417 / 2100) && all(runs_table$conditioning == "crown") &&
      all(runs_table$environment_df == 80L) && all(abs(runs_table$total_time - total_time) <= 1e-12),
    TRUE, !has_actual_control(runs_table) && !has_actual_control(summary_table) && !has_actual_control(diagnosis_table),
    TRUE, TRUE, TRUE,
    all(runs_table$analysis_signature == analysis_signature) && all(summary_table$analysis_signature == analysis_signature)
  ), stringsAsFactors = FALSE
)
qa$analysis_signature <- analysis_signature
qa <- clean_frame(qa)
if (!all(qa$pass)) stop("v7 QA failed: ", paste(qa$check[!qa$pass], collapse = ", "))

runs_path <- file.path(out_table_dir, "dt_sensitivity_runs_all_v7.tsv")
summary_path <- file.path(out_table_dir, "dt_sensitivity_summary_v7.tsv")
diagnosis_path <- file.path(out_table_dir, "dt_sensitivity_diagnosis_v7.tsv")
qa_path <- file.path(out_table_dir, "dt_sensitivity_QA_v7.tsv")
write_clean_tsv(runs_table, runs_path)
write_clean_tsv(summary_table, summary_path)
write_clean_tsv(diagnosis_table, diagnosis_path)
write_clean_tsv(qa, qa_path)

final_payload <- list(
  analysis_signature = analysis_signature, signature_components = signature_components,
  model = "BCSTDTempVar", diagnosis = diagnosis, diagnosis_table = diagnosis_table,
  formal_inference_dt = 0.005, formal_likelihood_changed = FALSE,
  v5_analysis_signature = v5$analysis_signature, v5_rds_sha256 = v5_sha,
  v6_analysis_signature = v6_sig_sha, v6_payload_sha256 = v6_payload_sha,
  all_runs = all_runs, all_run_sources = all_sources,
  summary = summary_table, qa = qa,
  input_sha256 = c(tree = tree_sha, environment = env_sha, v5_rds = v5_sha,
                   v6_signature = v6_sig_sha, v6_payload_sha),
  code_sha256 = c(driver = driver_sha)
)
final_path <- file.path(out_model_dir, "BCSTDTempVar_dt_sensitivity_final_v7.rds")
if (!file.exists(final_path)) atomic_save_rds(final_payload, final_path)
cat(sprintf(
  "v7 dt sensitivity complete: diagnosis=%s strict=%s slope=%.9g pulse_ratio=%.9g LL_range=%.9g signature=%s\n",
  diagnosis, all(strict), slope, pulse_ratio, logLik_range, analysis_signature
))
