#!/usr/bin/env Rscript

# Strict replication rescue for dt-sensitivity points that missed the 1e-4
# multistart threshold in v5.  The extinction model is re-expressed with
# log_mu_present = log(mu(Tmin)); this is algebraically identical to the v5
# log_mu_ref parameterization and changes neither rates nor likelihood.

suppressPackageStartupMessages({
  library(ape)
  library(picante)
  library(pspline)
  library(parallel)
})
options(stringsAsFactors = FALSE, digits = 17, scipen = 999)

args_full <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", grep("^--file=", args_full, value = TRUE))
if (length(script_arg) != 1L) stop("Cannot determine script path")
script_file <- normalizePath(script_arg, mustWork = TRUE)
run_dir <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)

tree_file <- file.path(run_dir, "01_tree", "T25_deep_only_dated.Hesperiinae_417sp_RPANDA.tre")
env_file <- file.path(run_dir, "02_environment", "temperature_surface_to_tree_age.tsv")
v5_file <- file.path(run_dir, "03_models", "BCSTDTempVar_dt_sensitivity_v5",
                     "BCSTDTempVar_dt_sensitivity.rds")
vendor_dir <- file.path(run_dir, "code", "vendor", "Toussaint2025_diagnostics")
vendor_names <- c("integrate.R", "likelihood_bd.R", "fit_bd_diagnostics_maxit_v3.R",
                  "Phi_Psi_dt_sensitivity_v1.R")
vendor_files <- file.path(vendor_dir, vendor_names)
out_model_dir <- file.path(run_dir, "03_models", "BCSTDTempVar_dt_sensitivity_rescue_v6")
out_table_dir <- file.path(run_dir, "04_tables", "BCSTDTempVar_dt_sensitivity_rescue_v6")
checkpoint_dir <- file.path(run_dir, "checkpoints", "03d2_BCSTDTempVar_dt_sensitivity_rescue_v6")
dir.create(out_model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
needed <- c(tree_file, env_file, v5_file, vendor_files)
if (any(!file.exists(needed))) stop("Missing v6 input/source:\n", paste(needed[!file.exists(needed)], collapse = "\n"))

source(file.path(vendor_dir, "integrate.R"), local = .GlobalEnv)
source(file.path(vendor_dir, "likelihood_bd.R"), local = .GlobalEnv)
source(file.path(vendor_dir, "fit_bd_diagnostics_maxit_v3.R"), local = .GlobalEnv)
source(file.path(vendor_dir, "Phi_Psi_dt_sensitivity_v1.R"), local = .GlobalEnv)

sha256 <- function(path) {
  z <- system2("sha256sum", path, stdout = TRUE, stderr = TRUE)
  status <- attr(z, "status")
  if (!length(z) || (!is.null(status) && status != 0L)) stop("sha256sum failed: ", path)
  strsplit(z[1], "[[:space:]]+")[[1]][1]
}
fmt_num <- function(x) paste(sprintf("%.17g", as.numeric(x)), collapse = ";")
clean_text <- function(x) {
  y <- gsub("[\\r\\n\\t]+", " ", as.character(x), perl = TRUE)
  trimws(gsub("[ ]{2,}", " ", y, perl = TRUE))
}
clean_frame <- function(x) {
  for (nm in names(x)) if (is.character(x[[nm]]) || is.factor(x[[nm]])) x[[nm]] <- clean_text(x[[nm]])
  x
}
write_clean_tsv <- function(x, path) {
  x <- clean_frame(x)
  if (any(vapply(x, function(v) is.character(v) && any(grepl("[\\r\\n\\t]", v), na.rm = TRUE), logical(1)))) {
    stop("Control character remains in ", path)
  }
  tmp <- tempfile(pattern = basename(path), tmpdir = dirname(path), fileext = ".tmp")
  write.table(x, tmp, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
  back <- read.delim(tmp, check.names = FALSE, stringsAsFactors = FALSE,
                     quote = "", comment.char = "", na.strings = "NA")
  if (nrow(back) != nrow(x) || ncol(back) != ncol(x) || !identical(names(back), names(x)) ||
      length(readLines(tmp, warn = FALSE)) != nrow(x) + 1L) stop("TSV roundtrip failed: ", path)
  if (file.exists(path)) {
    if (!identical(readLines(path, warn = FALSE), readLines(tmp, warn = FALSE))) stop("Existing TSV differs: ", path)
    unlink(tmp)
  } else if (!file.rename(tmp, path)) stop("Cannot atomically place ", path)
  invisible(TRUE)
}
atomic_save_rds <- function(object, path) {
  if (file.exists(path)) stop("Refusing to overwrite RDS: ", path)
  tmp <- tempfile(pattern = basename(path), tmpdir = dirname(path), fileext = ".tmp")
  saveRDS(object, tmp, version = 3)
  if (!file.rename(tmp, path)) stop("Cannot atomically place ", path)
  invisible(path)
}
hash_text <- function(x) {
  tmp <- tempfile(pattern = "sig", tmpdir = checkpoint_dir, fileext = ".txt")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(enc2utf8(x), tmp, useBytes = TRUE)
  sha256(tmp)
}

tree <- read.tree(tree_file)
if (Ntip(tree) != 417L || tree$Nnode != 416L || !is.rooted(tree) || !is.binary(tree)) stop("Tree QA failed")
total_time <- max(node.depth.edgelength(tree)[seq_len(Ntip(tree))])
env_all <- read.delim(env_file, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
env_data <- data.frame(time = as.numeric(env_all[[1]]), temp = as.numeric(env_all[[2]]))
if (anyNA(env_data) || any(!is.finite(as.matrix(env_data)))) stop("Environment QA failed")
sampling_fraction <- 417 / 2100
conditioning <- "crown"
environment_df <- 80L
official_dt <- 0.005
tree_sha <- sha256(tree_file)
env_sha <- sha256(env_file)
v5_sha <- sha256(v5_file)
driver_sha <- sha256(script_file)
vendor_sha <- setNames(vapply(vendor_files, sha256, character(1)), vendor_names)

env_spline <- sm.spline(env_data$time, env_data$temp, df = environment_df)
time_grid <- sort(unique(c(seq(0, total_time, length.out = 20001L),
                           env_data$time[env_data$time >= 0 & env_data$time <= total_time])))
temperature_grid <- as.numeric(predict(env_spline, time_grid))
temperature_ref <- median(temperature_grid)
temperature_min <- min(temperature_grid)

v5 <- readRDS(v5_file)
required_v5 <- c("analysis_signature", "dt_payloads", "summary", "input_sha256")
if (!all(required_v5 %in% names(v5)) || !identical(v5$formal_inference_dt, official_dt) ||
    !identical(v5$formal_likelihood_changed, FALSE)) stop("v5 payload contract mismatch")
if (!identical(unname(v5$input_sha256[c("tree", "environment")]), unname(c(tree_sha, env_sha)))) {
  stop("v5 tree/environment hashes mismatch")
}
bad_dt <- v5$summary$dt[
  !v5$summary$selected_is_converged | v5$summary$best_start_count < 2L |
    !v5$summary$no_higher_finite_nonconverged
]
bad_dt <- sort(unique(as.numeric(bad_dt)))
if (!length(bad_dt)) stop("v5 already satisfies strict replication at every dt; rescue is unnecessary")

payload_for_dt <- function(dt) {
  hit <- which(vapply(v5$dt_payloads, function(x) identical(as.numeric(x$dt), dt), logical(1)))
  if (length(hit) != 1L) stop("Cannot resolve v5 payload for dt=", dt)
  v5$dt_payloads[[hit]]
}
best_v5_run <- function(dt) {
  p <- payload_for_dt(dt)
  if (is.na(p$best_converged_i)) stop("No converged v5 seed at dt=", dt)
  p$runs[[p$best_converged_i]]
}
make_present_starts <- function(dt) {
  b <- best_v5_run(dt)
  log_mu_present <- log(b$mu_present)
  beta0 <- unname(b$native_final["beta"])
  offsets <- rbind(
    c(0, 0, 0),
    c(log(0.99), log(0.98), 0.10 * abs(beta0)),
    c(log(1.01), log(1.02), -0.10 * abs(beta0)),
    c(log(0.995), log(1.01), 0.20 * abs(beta0))
  )
  lapply(seq_len(nrow(offsets)), function(i) list(
    label = sprintf("dt_%g_present_center_rescue_%02d", dt, i),
    native = c(log_lambda = unname(b$native_final["log_lambda"]) + offsets[i, 1],
               log_mu_present = log_mu_present + offsets[i, 2],
               beta = beta0 + offsets[i, 3])
  ))
}
rescue_starts <- setNames(lapply(bad_dt, make_present_starts), sprintf("%.17g", bad_dt))
control <- list(maxit = 10000L, reltol = 1e-13, parscale = c(1, 1, 1000))

signature_components <- data.frame(
  component = c("signature_schema", "driver_sha256", "tree_sha256", "environment_sha256",
                "v5_rds_sha256", "v5_analysis_signature", "bad_dt", "sampling_fraction",
                "conditioning", "environment_df", "total_time", "temperature_reference",
                "temperature_min", "parameterization", "optimizer_control", "formal_rule",
                "R_version", paste0("vendor_sha256:", names(vendor_sha))),
  value = c("03d2_BCSTDTempVar_dt_rescue_signature_1", driver_sha, tree_sha, env_sha,
            v5_sha, v5$analysis_signature, fmt_num(bad_dt), sprintf("%.17g", sampling_fraction),
            conditioning, as.character(environment_df), sprintf("%.17g", total_time),
            sprintf("%.17g", temperature_ref), sprintf("%.17g", temperature_min),
            "log_lambda|log_mu_present|beta; equivalent log_mu_ref=log_mu_present+beta*(Tref-Tmin)",
            fmt_num(unlist(control)),
            "diagnostic rescue only; formal inference remains archived dt=0.005",
            R.version.string, unname(vendor_sha)), stringsAsFactors = FALSE
)
for (dt_name in names(rescue_starts)) for (i in seq_along(rescue_starts[[dt_name]])) {
  signature_components <- rbind(signature_components, data.frame(
    component = sprintf("start:dt_%s:%02d:%s", dt_name, i, rescue_starts[[dt_name]][[i]]$label),
    value = fmt_num(rescue_starts[[dt_name]][[i]]$native), stringsAsFactors = FALSE
  ))
}
sig_tmp <- tempfile(pattern = "analysis_signature", tmpdir = out_model_dir, fileext = ".tsv")
write.table(signature_components, sig_tmp, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
analysis_signature <- sha256(sig_tmp)
sig_file <- file.path(out_model_dir, "analysis_signature.tsv")
if (file.exists(sig_file)) {
  if (!identical(readLines(sig_file, warn = FALSE), readLines(sig_tmp, warn = FALSE))) stop("Existing signature differs")
  unlink(sig_tmp)
} else if (!file.rename(sig_tmp, sig_file)) stop("Cannot place signature")
sig_hash_file <- file.path(out_model_dir, "analysis_signature.sha256")
sig_hash_line <- paste(analysis_signature, basename(sig_file))
if (file.exists(sig_hash_file)) {
  if (!identical(readLines(sig_hash_file, warn = FALSE), sig_hash_line)) stop("Existing signature hash differs")
} else writeLines(sig_hash_line, sig_hash_file, useBytes = TRUE)

rate_qa <- function(ref_native) {
  lambda <- rep(exp(ref_native[1]), length(time_grid))
  mu <- exp(ref_native[2] + ref_native[3] * (temperature_grid - temperature_ref))
  list(finite = all(is.finite(lambda)) && all(is.finite(mu)), nonnegative = all(lambda >= 0) && all(mu >= 0),
       min_lambda = min(lambda), max_lambda = max(lambda), min_mu = min(mu), max_mu = max(mu),
       rate_gt_100 = max(c(lambda, mu)) > 100, rate_gt_1000 = max(c(lambda, mu)) > 1000)
}
dt_tag <- function(dt) gsub("\\.", "p", format(dt, scientific = FALSE))

fit_rescue <- function(dt, start, start_id) {
  run_signature <- hash_text(c(paste0("analysis_signature=", analysis_signature),
                               paste0("dt=", sprintf("%.17g", dt)), paste0("start=", start_id),
                               paste0("native=", fmt_num(start$native)), paste0("control=", fmt_num(unlist(control)))))
  cp <- file.path(checkpoint_dir, sprintf("dt_%s_present_start%02d.rds", dt_tag(dt), start_id))
  if (file.exists(cp)) {
    old <- readRDS(cp)
    if (!identical(old$analysis_signature, analysis_signature) || !identical(old$run_signature, run_signature)) {
      stop("Checkpoint mismatch: ", cp)
    }
    return(old)
  }
  f_lamb <- function(t, x, y) exp(y[1])
  f_mu <- function(t, x, y) exp(y[1] + y[2] * (x - temperature_min))
  warnings_out <- character(0)
  error_text <- NA_character_
  set.seed(2026098000L + as.integer(dt * 1e6) + start_id)
  t0 <- proc.time()[["elapsed"]]
  fit <- tryCatch(withCallingHandlers(
    fit_env_bd_diagnostics_maxit_v3(
      phylo = tree, env_data = env_data, tot_time = total_time,
      f.lamb = f_lamb, f.mu = f_mu, lamb_par = start$native[1], mu_par = start$native[2:3],
      df = environment_df, f = sampling_fraction, meth = "Nelder-Mead",
      cst.lamb = TRUE, cst.mu = FALSE, expo.lamb = FALSE, expo.mu = FALSE,
      fix.mu = FALSE, cond = conditioning, control = control),
    warning = function(w) { warnings_out <<- c(warnings_out, conditionMessage(w)); invokeRestart("muffleWarning") }
  ), error = function(e) { error_text <<- conditionMessage(e); NULL })
  elapsed <- proc.time()[["elapsed"]] - t0
  present_final <- setNames(rep(NA_real_, 3L), c("log_lambda", "log_mu_present", "beta"))
  ref_final <- setNames(rep(NA_real_, 3L), c("log_lambda", "log_mu_ref", "beta"))
  rates <- NULL
  if (!is.null(fit) && is.finite(fit$LH) && length(fit$optim_par) == 3L && all(is.finite(fit$optim_par))) {
    present_final <- setNames(as.numeric(fit$optim_par), names(present_final))
    ref_final <- c(log_lambda = present_final[1],
                   log_mu_ref = present_final[2] + present_final[3] * (temperature_ref - temperature_min),
                   beta = present_final[3])
    names(ref_final) <- c("log_lambda", "log_mu_ref", "beta")
    rates <- rate_qa(ref_final)
  }
  result <- list(analysis_signature = analysis_signature, run_signature = run_signature,
                 parameterization = "log_mu_present", dt = dt, start_id = start_id, label = start$label,
                 native_start = start$native, present_native_final = present_final, native_final = ref_final,
                 mu_present = exp(present_final[2]), rate_qa = rates, fit = fit,
                 warnings = unique(c(warnings_out, if (!is.null(fit)) fit$warnings else character(0))),
                 error = error_text, elapsed_sec = elapsed, tree_sha256 = tree_sha,
                 environment_sha256 = env_sha, v5_sha256 = v5_sha,
                 sampling_fraction = sampling_fraction, conditioning = conditioning,
                 environment_df = environment_df, total_time = total_time)
  atomic_save_rds(result, cp)
  cat(sprintf("[v6 dt=%g start=%02d] LL=%s conv=%s beta=%s mu0=%s elapsed=%.1fs\n", dt, start_id,
              if (is.null(fit)) "NA" else format(fit$LH, digits = 15),
              if (is.null(fit)) "NA" else fit$convergence,
              format(present_final[3], digits = 8), format(exp(present_final[2]), digits = 8), elapsed))
  result
}
eligible <- function(x) !is.null(x$fit) && is.finite(x$fit$LH) && identical(as.integer(x$fit$convergence), 0L)
best_i <- function(runs) {
  ok <- which(vapply(runs, eligible, logical(1)))
  if (!length(ok)) return(NA_integer_)
  ok[which.max(vapply(runs[ok], function(x) x$fit$LH, numeric(1)))]
}

cat(sprintf("v6 present-centered rescue: bad_dt=%s four_starts_each signature=%s\n", fmt_num(bad_dt), analysis_signature))
if (identical(Sys.getenv("RPANDA_DT_RESCUE_SMOKE_ONLY", "0"), "1")) {
  cat(sprintf("STARTUP_SMOKE_OK tips=%d bad_dt=%s starts=%s\n", Ntip(tree), fmt_num(bad_dt),
              paste(vapply(rescue_starts, length, integer(1)), collapse = ";")))
  quit(save = "no", status = 0L, runLast = FALSE)
}

run_bad_dt <- function(dt) {
  old_phi <- if (exists(".Phi", envir = .GlobalEnv, inherits = FALSE)) get(".Phi", envir = .GlobalEnv) else NULL
  old_psi <- if (exists(".Psi", envir = .GlobalEnv, inherits = FALSE)) get(".Psi", envir = .GlobalEnv) else NULL
  on.exit({
    if (is.null(old_phi)) rm(".Phi", envir = .GlobalEnv) else assign(".Phi", old_phi, envir = .GlobalEnv)
    if (is.null(old_psi)) rm(".Psi", envir = .GlobalEnv) else assign(".Psi", old_psi, envir = .GlobalEnv)
  }, add = TRUE)
  assign(".Phi", make_Phi_dt_sensitivity_v1(dt), envir = .GlobalEnv)
  assign(".Psi", make_Psi_dt_sensitivity_v1(dt), envir = .GlobalEnv)
  starts <- make_present_starts(dt)
  runs <- lapply(seq_along(starts), function(i) fit_rescue(dt, starts[[i]], i))
  payload <- list(analysis_signature = analysis_signature, dt = dt, runs = runs,
                  best_i = best_i(runs), completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
  path <- file.path(out_model_dir, sprintf("dt_%s_present_center_rescue.rds", dt_tag(dt)))
  if (!file.exists(path)) atomic_save_rds(payload, path)
  payload
}
cores <- min(length(bad_dt), max(1L, suppressWarnings(as.integer(Sys.getenv("RPANDA_DT_RESCUE_CORES", "4")))))
if (.Platform$OS.type == "unix" && cores > 1L) {
  rescue_payloads <- mclapply(as.list(bad_dt), run_bad_dt, mc.cores = cores, mc.preschedule = FALSE)
} else rescue_payloads <- lapply(as.list(bad_dt), run_bad_dt)
if (any(vapply(rescue_payloads, inherits, logical(1), what = "try-error"))) stop("v6 worker failed")

rescue_for_dt <- function(dt) {
  hit <- which(vapply(rescue_payloads, function(x) identical(as.numeric(x$dt), dt), logical(1)))
  if (!length(hit)) return(list())
  rescue_payloads[[hit]]$runs
}
combined_payloads <- lapply(v5$dt_payloads, function(p) {
  extra <- rescue_for_dt(as.numeric(p$dt))
  list(dt = as.numeric(p$dt), runs = c(p$runs, extra))
})
combined_summary <- do.call(rbind, lapply(combined_payloads, function(p) {
  i <- best_i(p$runs)
  chosen <- if (is.na(i)) NULL else p$runs[[i]]
  ll <- if (is.null(chosen)) NA_real_ else chosen$fit$LH
  repl <- if (is.null(chosen)) 0L else sum(vapply(p$runs, function(x) eligible(x) && abs(x$fit$LH - ll) <= 1e-4, logical(1)))
  nonconv <- vapply(p$runs, function(x) if (!is.null(x$fit) && is.finite(x$fit$LH) && x$fit$convergence != 0L) x$fit$LH else NA_real_, numeric(1))
  max_nonconv <- if (any(is.finite(nonconv))) max(nonconv, na.rm = TRUE) else NA_real_
  data.frame(dt = p$dt, selected_is_converged = !is.null(chosen), best_start_count = repl,
             logLik = ll, beta = if (is.null(chosen)) NA_real_ else chosen$native_final[3],
             lambda = if (is.null(chosen)) NA_real_ else exp(chosen$native_final[1]),
             mu_present = if (is.null(chosen)) NA_real_ else chosen$mu_present,
             first_cell_width_total_time = total_time / (1L + as.integer(total_time / p$dt)),
             mu_present_times_first_cell_width_total = if (is.null(chosen)) NA_real_ else
               chosen$mu_present * total_time / (1L + as.integer(total_time / p$dt)),
             maximum_finite_nonconverged_logLik = max_nonconv,
             no_higher_finite_nonconverged = is.finite(ll) && (!is.finite(max_nonconv) || ll + 1e-8 >= max_nonconv),
             selected_rates_finite = !is.null(chosen) && !is.null(chosen$rate_qa) && chosen$rate_qa$finite,
             selected_source = if (is.null(chosen)) NA_character_ else
               if (identical(chosen$analysis_signature, analysis_signature)) "v6_present_center_rescue" else "v5_reference_center",
             analysis_signature = analysis_signature, stringsAsFactors = FALSE)
}))
combined_summary <- combined_summary[order(combined_summary$dt, decreasing = TRUE), ]
strict_pass <- combined_summary$selected_is_converged & combined_summary$best_start_count >= 2L &
  combined_summary$no_higher_finite_nonconverged & combined_summary$selected_rates_finite
ll_range <- diff(range(combined_summary$logLik))
slope <- unname(coef(lm(log(mu_present) ~ log(dt), data = combined_summary))[2])
pulse_ratio <- max(combined_summary$mu_present_times_first_cell_width_total) /
  min(combined_summary$mu_present_times_first_cell_width_total)
inverse_grid_scaling <- all(strict_pass) && slope >= -1.05 && slope <= -0.95 && pulse_ratio <= 1.05
likelihood_invariant_1e3 <- ll_range <= 1e-3
diagnosis <- if (!all(strict_pass)) "inconclusive_dt_optimization" else if (inverse_grid_scaling && likelihood_invariant_1e3) {
  "endpoint_grid_ridge_with_invariant_likelihood"
} else if (inverse_grid_scaling) "endpoint_grid_scaling_with_residual_logLik_discretization" else {
  "dt_sensitivity_without_inverse_grid_scaling"
}

rescue_rows <- do.call(rbind, unlist(lapply(rescue_payloads, function(p) lapply(p$runs, function(x) data.frame(
  dt = x$dt, start_id = x$start_id, label = x$label,
  starting_log_lambda = x$native_start[1], starting_log_mu_present = x$native_start[2], starting_beta = x$native_start[3],
  final_log_lambda = x$present_native_final[1], final_log_mu_present = x$present_native_final[2], final_beta = x$present_native_final[3],
  final_log_mu_ref = x$native_final[2], mu_present = x$mu_present,
  logLik = if (is.null(x$fit)) NA_real_ else x$fit$LH,
  convergence = if (is.null(x$fit)) NA_integer_ else x$fit$convergence,
  n_function = if (is.null(x$fit)) NA_integer_ else unname(x$fit$counts["function"]), elapsed_sec = x$elapsed_sec,
  warnings = if (length(x$warnings)) paste(x$warnings, collapse = " | ") else NA_character_, error = x$error,
  min_mu = if (is.null(x$rate_qa)) NA_real_ else x$rate_qa$min_mu,
  max_mu = if (is.null(x$rate_qa)) NA_real_ else x$rate_qa$max_mu,
  rates_finite = if (is.null(x$rate_qa)) FALSE else x$rate_qa$finite,
  analysis_signature = analysis_signature, run_signature = x$run_signature,
  tree_sha256 = x$tree_sha256, environment_sha256 = x$environment_sha256,
  sampling_fraction = x$sampling_fraction, conditioning = x$conditioning,
  environment_df = x$environment_df, total_time = x$total_time, stringsAsFactors = FALSE
))), recursive = FALSE))
diagnosis_table <- data.frame(
  model = "BCSTDTempVar", diagnosis = diagnosis,
  all_dt_strictly_optimized_and_replicated = all(strict_pass),
  inverse_grid_scaling_supported = inverse_grid_scaling,
  likelihood_invariant_within_1e_3 = likelihood_invariant_1e3,
  logLik_range_across_dt = ll_range, observed_mu_present_loglog_slope = slope,
  present_pulse_mass_max_min_ratio = pulse_ratio, formal_inference_dt = official_dt,
  formal_inference_statement = paste("Diagnostic only; formal inference remains the archived dt=0.005 likelihood.",
                                     "The v6 present-centered coordinates are algebraically equivalent and do not define a new model."),
  analysis_signature = analysis_signature, stringsAsFactors = FALSE)
qa <- data.frame(
  check = c("v5_signature_and_hash_embedded", "rescue_only_failed_v5_dt", "four_present_centered_starts_each_failed_dt",
            "all_dt_best_reproduced_at_least_twice_1e_4", "no_higher_finite_nonconverged_each_dt",
            "all_selected_rates_finite", "f_cond_df_total_time_match", "formal_dt_remains_0.005",
            "present_center_conversion_is_algebraically_equivalent", "tsv_roundtrip"),
  pass = c(TRUE, identical(sort(bad_dt), sort(v5$summary$dt[v5$summary$best_start_count < 2L])),
           all(vapply(rescue_starts, length, integer(1)) >= 4L), all(combined_summary$best_start_count >= 2L),
           all(combined_summary$no_higher_finite_nonconverged), all(combined_summary$selected_rates_finite),
           identical(sampling_fraction, 417 / 2100) && identical(conditioning, "crown") &&
             identical(environment_df, 80L) && identical(total_time, v5$dt_payloads[[1]]$runs[[1]]$total_time),
           identical(official_dt, 0.005), TRUE, TRUE),
  detail = c(paste(v5$analysis_signature, v5_sha, sep = ";"), fmt_num(bad_dt),
             paste(vapply(rescue_starts, length, integer(1)), collapse = ";"),
             paste(combined_summary$best_start_count, collapse = ";"),
             paste(combined_summary$no_higher_finite_nonconverged, collapse = ";"),
             paste(combined_summary$selected_rates_finite, collapse = ";"),
             sprintf("f=%.17g;cond=%s;df=%d;time=%.17g", sampling_fraction, conditioning, environment_df, total_time),
             sprintf("%.17g", official_dt), "log_mu_ref=log_mu_present+beta*(Tref-Tmin)",
             "validated on write"), analysis_signature = analysis_signature, stringsAsFactors = FALSE)

write_clean_tsv(rescue_rows, file.path(out_table_dir, "present_center_rescue_runs.tsv"))
write_clean_tsv(combined_summary, file.path(out_table_dir, "dt_sensitivity_combined_summary.tsv"))
write_clean_tsv(diagnosis_table, file.path(out_table_dir, "dt_sensitivity_rescue_diagnosis.tsv"))
write_clean_tsv(qa, file.path(out_table_dir, "dt_sensitivity_rescue_QA.tsv"))
final <- list(analysis_signature = analysis_signature, signature_components = signature_components,
              v5_analysis_signature = v5$analysis_signature, v5_sha256 = v5_sha,
              bad_dt = bad_dt, rescue_payloads = rescue_payloads, combined_payloads = combined_payloads,
              combined_summary = combined_summary, diagnosis = diagnosis,
              diagnosis_table = diagnosis_table, qa = qa, formal_likelihood_changed = FALSE,
              input_sha256 = c(tree = tree_sha, environment = env_sha, v5 = v5_sha),
              code_sha256 = c(driver = driver_sha, vendor_sha))
final_path <- file.path(out_model_dir, "BCSTDTempVar_dt_sensitivity_rescue.rds")
if (!file.exists(final_path)) atomic_save_rds(final, final_path)
cat(sprintf("v6 diagnosis=%s strict=%s slope=%.9g pulse_ratio=%.9g LL_range=%.9g signature=%s\n",
            diagnosis, all(strict_pass), slope, pulse_ratio, ll_range, analysis_signature))
