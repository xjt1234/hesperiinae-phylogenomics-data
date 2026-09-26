#!/usr/bin/env Rscript

# Diagnose the BCSTDTempVar likelihood surface without changing its likelihood.
# The free fit uses [log(lambda0), log(mu_ref), beta], where
# mu_ref = mu(Tref). Fixed-beta profiles optimize the other two parameters.
# A finite MLE is accepted only with replicated convergence and an interior
# profile peak; a monotonically improving or asymptotically flat negative-beta
# tail is explicitly classified as boundary_nonidentifiable.

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
round1_file <- file.path(run_dir, "03_models", "model_fits_10models.rds")
vendor_dir <- file.path(run_dir, "code", "vendor", "Toussaint2025_diagnostics")
out_model_dir <- file.path(run_dir, "03_models", "BCSTDTempVar_logcenter_profile_v4")
out_table_dir <- file.path(run_dir, "04_tables", "BCSTDTempVar_logcenter_profile_v4")
checkpoint_dir <- file.path(run_dir, "checkpoints", "03c_BCSTDTempVar_logcenter_profile_v4")
dir.create(out_model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

vendor_names <- c(
  "integrate.R", "Phi.R", "Psi.R", "likelihood_bd.R",
  "fit_bd_diagnostics_maxit_v3.R"
)
vendor_files <- file.path(vendor_dir, vendor_names)
needed <- c(tree_file, env_file, round1_file, vendor_files)
if (any(!file.exists(needed))) {
  stop("Missing required profile input/source:\n", paste(needed[!file.exists(needed)], collapse = "\n"))
}
source(file.path(vendor_dir, "integrate.R"), local = .GlobalEnv)
source(file.path(vendor_dir, "Phi.R"), local = .GlobalEnv)
source(file.path(vendor_dir, "Psi.R"), local = .GlobalEnv)
source(file.path(vendor_dir, "likelihood_bd.R"), local = .GlobalEnv)
source(file.path(vendor_dir, "fit_bd_diagnostics_maxit_v3.R"), local = .GlobalEnv)

sha256 <- function(path) {
  out <- system2("sha256sum", path, stdout = TRUE, stderr = TRUE)
  status <- attr(out, "status")
  if (!length(out) || (!is.null(status) && status != 0L)) stop("sha256sum failed: ", path)
  strsplit(out[1], "[[:space:]]+")[[1]][1]
}
hash_text <- function(x) {
  tmp <- tempfile(pattern = "signature_", tmpdir = checkpoint_dir, fileext = ".txt")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(enc2utf8(x), tmp, useBytes = TRUE)
  sha256(tmp)
}
fmt_num <- function(x) paste(sprintf("%.17g", as.numeric(x)), collapse = ";")
canonical_start <- function(x) paste(sprintf("%.10g", as.numeric(x)), collapse = ";")
clean_text <- function(x) {
  y <- as.character(x)
  y <- gsub("[\\r\\n\\t]+", " ", y, perl = TRUE)
  y <- gsub("[ ]{2,}", " ", y, perl = TRUE)
  trimws(y)
}
clean_frame <- function(x) {
  for (nm in names(x)) {
    if (is.character(x[[nm]]) || is.factor(x[[nm]])) x[[nm]] <- clean_text(x[[nm]])
  }
  x
}
write_clean_tsv <- function(x, path) {
  x <- clean_frame(x)
  bad <- any(vapply(x, function(v) {
    is.character(v) && any(grepl("[\\r\\n\\t]", v, perl = TRUE), na.rm = TRUE)
  }, logical(1)))
  if (bad) stop("Unclean control character remains: ", path)
  if (file.exists(path)) stop("Refusing to overwrite profile output: ", path)
  tmp <- tempfile(pattern = paste0(basename(path), "_"), tmpdir = dirname(path), fileext = ".tmp")
  write.table(x, tmp, sep = "\t", quote = FALSE, row.names = FALSE,
              na = "NA", fileEncoding = "UTF-8")
  back <- read.delim(tmp, check.names = FALSE, stringsAsFactors = FALSE,
                     quote = "", comment.char = "", na.strings = "NA")
  if (nrow(back) != nrow(x) || ncol(back) != ncol(x) ||
      !identical(names(back), names(x)) ||
      length(readLines(tmp, warn = FALSE)) != nrow(x) + 1L) {
    stop("TSV roundtrip failed: ", path)
  }
  if (!file.rename(tmp, path)) stop("Could not atomically place TSV: ", path)
  invisible(TRUE)
}
atomic_save_rds <- function(object, path) {
  if (file.exists(path)) stop("Refusing to overwrite profile RDS: ", path)
  tmp <- tempfile(pattern = paste0(basename(path), "_"), tmpdir = dirname(path), fileext = ".tmp")
  saveRDS(object, tmp, version = 3)
  if (!file.rename(tmp, path)) stop("Could not atomically place RDS: ", path)
  invisible(path)
}

tree <- read.tree(tree_file)
if (Ntip(tree) != 417L || tree$Nnode != 416L || !is.rooted(tree) || !is.binary(tree)) {
  stop("Formal tree failed 417-tip rooted-binary validation")
}
if (any(!is.finite(tree$edge.length)) || any(tree$edge.length <= 0)) {
  stop("Formal tree contains non-positive/non-finite branches")
}
tip_depths <- node.depth.edgelength(tree)[seq_len(Ntip(tree))]
if (diff(range(tip_depths)) >= 1e-5) stop("Root-to-tip range exceeds 1e-5 Ma")
total_time <- max(tip_depths)

env_all <- read.delim(env_file, check.names = FALSE, stringsAsFactors = FALSE,
                      quote = "", comment.char = "")
env_data <- data.frame(time = as.numeric(env_all[[1]]), temp = as.numeric(env_all[[2]]))
if (anyNA(env_data) || any(!is.finite(as.matrix(env_data))) ||
    anyDuplicated(env_data$time) || is.unsorted(env_data$time, strictly = TRUE)) {
  stop("Environment table is not finite, unique, and strictly increasing")
}
if (min(env_data$time) > 0 || max(env_data$time) < total_time) {
  stop("Environment does not cover formal tree age")
}

sampling_fraction <- 417 / 2100
conditioning <- "crown"
environment_df <- 80L
n_obs <- Ntip(tree)
tree_sha <- sha256(tree_file)
env_sha <- sha256(env_file)
round1_sha <- sha256(round1_file)
driver_sha <- sha256(script_file)
vendor_sha <- setNames(vapply(vendor_files, sha256, character(1)), vendor_names)

env_spline <- sm.spline(env_data$time, env_data$temp, df = environment_df)
time_grid <- sort(unique(c(
  seq(0, total_time, length.out = 20001L),
  env_data$time[env_data$time >= 0 & env_data$time <= total_time]
)))
temperature_grid <- as.numeric(predict(env_spline, time_grid))
if (any(!is.finite(temperature_grid))) stop("Non-finite df=80 spline prediction")
temperature_ref <- median(temperature_grid)
temperature_min <- min(temperature_grid)
temperature_max <- max(temperature_grid)

round1 <- readRDS(round1_file)
if (!all(c("analysis_signature", "all_results", "total_time", "sampling_fraction",
           "df", "cond", "input_sha256") %in% names(round1))) {
  stop("First-round RDS lacks required fields")
}
if (!identical(round1$sampling_fraction, sampling_fraction) ||
    !identical(round1$cond, conditioning) ||
    !identical(as.integer(round1$df), environment_df) ||
    !identical(round1$total_time, total_time) ||
    !identical(unname(round1$input_sha256[c("tree", "environment")]),
               unname(c(tree_sha, env_sha)))) {
  stop("First-round RDS configuration/input mismatch")
}
old_result <- round1$all_results$BCSTDTempVar
if (is.null(old_result) || length(old_result$runs) < 1L) {
  stop("First-round BCSTDTempVar runs unavailable")
}
old_nonconv <- Filter(function(x) {
  !is.null(x$fit) && is.finite(x$fit$LH) && x$fit$convergence != 0L
}, old_result$runs)
if (!length(old_nonconv)) stop("No finite nonconverged endpoint available for log-centered starts")

raw_to_logcenter <- function(raw) {
  lambda0 <- abs(unname(raw[1]))
  mu0 <- abs(unname(raw[2]))
  beta <- unname(raw[3])
  mu_ref <- mu0 * exp(beta * temperature_ref)
  if (!is.finite(lambda0) || !is.finite(mu_ref) || lambda0 <= 0 || mu_ref <= 0) {
    stop("Cannot convert endpoint to finite log-centered parameters")
  }
  c(log_lambda = log(lambda0), log_mu_ref = log(mu_ref), beta = beta)
}
logcenter_to_original <- function(native) {
  log_lambda <- unname(native[1])
  log_mu_ref <- unname(native[2])
  beta <- unname(native[3])
  log_mu0 <- log_mu_ref - beta * temperature_ref
  c(lambda0 = exp(log_lambda), mu0 = exp(log_mu0), beta = beta,
    log_mu_ref = log_mu_ref, log_mu0 = log_mu0)
}
profile_to_original <- function(native, beta) {
  logcenter_to_original(c(log_lambda = native[1], log_mu_ref = native[2], beta = beta))
}
equivalence_qa <- function(native) {
  converted <- logcenter_to_original(native)
  log_mu_center <- native[2] + native[3] * (temperature_grid - temperature_ref)
  mu_center <- exp(log_mu_center)
  mu_original <- converted["mu0"] * exp(converted["beta"] * temperature_grid)
  max_abs <- if (all(is.finite(mu_original))) max(abs(mu_center - mu_original)) else Inf
  exponent_identity <- max(abs(
    log_mu_center - (converted["log_mu0"] + converted["beta"] * temperature_grid)
  ))
  list(max_abs_rate_difference = max_abs,
       max_abs_log_exponent_difference = exponent_identity)
}
rate_qa <- function(native) {
  lambda <- rep(exp(native[1]), length(time_grid))
  mu <- exp(native[2] + native[3] * (temperature_grid - temperature_ref))
  finite <- all(is.finite(lambda)) && all(is.finite(mu))
  nonnegative <- all(lambda >= 0) && all(mu >= 0)
  max_rate <- if (finite) max(c(lambda, mu)) else Inf
  flags <- character(0)
  if (finite && max(mu) < 1e-7) flags <- c(flags, "mu_rate_near_zero_over_observed_T")
  if (!finite) flags <- c(flags, "nonfinite_rate")
  if (!nonnegative) flags <- c(flags, "negative_rate")
  if (max_rate > 100) flags <- c(flags, "rate_gt_100_per_Ma")
  list(
    finite = finite, nonnegative = nonnegative, explosive = max_rate > 100,
    min_lambda = suppressWarnings(min(lambda)), max_lambda = suppressWarnings(max(lambda)),
    min_mu = suppressWarnings(min(mu)), max_mu = suppressWarnings(max(mu)),
    boundary_flag = if (length(flags)) paste(flags, collapse = ";") else "none"
  )
}

make_full_start <- function(label, native) {
  native <- as.numeric(native)
  names(native) <- c("log_lambda", "log_mu_ref", "beta")
  list(label = label, native = native)
}
full_starts <- lapply(seq_along(old_nonconv), function(i) {
  p <- as.numeric(old_nonconv[[i]]$fit$optim_par)
  make_full_start(paste0("round1_nonconverged_endpoint_", sprintf("%02d", i)),
                  raw_to_logcenter(p))
})
boundary_betas <- c(-2, -4, -6, -8, -10, -12, -16)
cold_mu_targets <- c(0.05, 0.20, 0.75, 2.0, 5.0, 1.0, 3.0)
lambda_targets <- c(0.15, 0.16, 0.17, 0.16, 0.165, 0.155, 0.17)
for (i in seq_along(boundary_betas)) {
  beta <- boundary_betas[i]
  log_mu_ref <- log(cold_mu_targets[i]) - beta * (temperature_min - temperature_ref)
  full_starts[[length(full_starts) + 1L]] <- make_full_start(
    sprintf("boundary_calibrated_beta_%g", beta),
    c(log(lambda_targets[i]), log_mu_ref, beta)
  )
}
if (length(full_starts) < 6L ||
    length(unique(vapply(full_starts, function(x) canonical_start(x$native), character(1)))) < 6L) {
  stop("Fewer than six distinct full-model log-centered starts")
}

base_profile_betas <- c(0, -0.5, -1, -2, -4, -6, -8, -10, -12, -16, -20)
extension_betas_1 <- c(-24, -32)
extension_betas_2 <- c(-48)
profile_lambda_targets <- c(0.14, 0.16, 0.18, 0.21)
profile_cold_mu_targets <- c(0.01, 0.10, 1.0, 5.0)
make_profile_starts <- function(beta) {
  lapply(seq_along(profile_lambda_targets), function(i) {
    log_mu_ref <- log(profile_cold_mu_targets[i]) - beta * (temperature_min - temperature_ref)
    p <- c(log_lambda = log(profile_lambda_targets[i]), log_mu_ref = log_mu_ref)
    list(label = sprintf("beta_%g_profile_start%02d", beta, i), native = p)
  })
}

full_control <- list(maxit = 15000L, reltol = 1e-11, parscale = c(1, 10, 1))
profile_control <- list(maxit = 8000L, reltol = 1e-11, parscale = c(1, 10))
signature_components <- data.frame(
  component = c(
    "signature_schema", "driver_sha256", "tree_sha256", "environment_sha256",
    "round1_rds_sha256", "round1_analysis_signature", "sampling_fraction",
    "conditioning", "environment_df", "total_time", "temperature_reference",
    "temperature_min", "temperature_max", "parameterization",
    "optimizer", "full_control", "profile_control", "base_profile_betas",
    "extension_betas_1", "extension_betas_2", "classification_rule", "R_version",
    paste0("package_version:", c("ape", "picante", "pspline")),
    paste0("vendor_sha256:", names(vendor_sha))
  ),
  value = c(
    "03c_BCSTDTempVar_logcenter_profile_signature_1", driver_sha, tree_sha, env_sha,
    round1_sha, round1$analysis_signature, sprintf("%.17g", sampling_fraction),
    conditioning, as.character(environment_df), sprintf("%.17g", total_time),
    sprintf("%.17g", temperature_ref), sprintf("%.17g", temperature_min),
    sprintf("%.17g", temperature_max), "log_lambda|log_mu_ref|beta",
    "Nelder-Mead; archived likelihood unchanged", fmt_num(unlist(full_control)),
    fmt_num(unlist(profile_control)), fmt_num(base_profile_betas),
    fmt_num(extension_betas_1), fmt_num(extension_betas_2),
    "finite only if replicated full convergence and clear interior profile peak; improving/flat tail is boundary_nonidentifiable",
    R.version.string,
    vapply(c("ape", "picante", "pspline"), function(x) as.character(packageVersion(x)), character(1)),
    unname(vendor_sha)
  ),
  stringsAsFactors = FALSE
)
for (i in seq_along(full_starts)) {
  signature_components <- rbind(signature_components, data.frame(
    component = sprintf("full_start:%02d:%s", i, full_starts[[i]]$label),
    value = fmt_num(full_starts[[i]]$native), stringsAsFactors = FALSE
  ))
}
for (beta in c(base_profile_betas, extension_betas_1, extension_betas_2)) {
  starts <- make_profile_starts(beta)
  for (i in seq_along(starts)) {
    signature_components <- rbind(signature_components, data.frame(
      component = sprintf("profile_start:beta_%g:%02d", beta, i),
      value = fmt_num(starts[[i]]$native), stringsAsFactors = FALSE
    ))
  }
}

signature_file <- file.path(out_model_dir, "analysis_signature.tsv")
signature_tmp <- tempfile(pattern = "analysis_signature_", tmpdir = out_model_dir, fileext = ".tsv")
write.table(signature_components, signature_tmp, sep = "\t", quote = FALSE,
            row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
analysis_signature <- sha256(signature_tmp)
if (file.exists(signature_file)) {
  if (!identical(readLines(signature_file, warn = FALSE), readLines(signature_tmp, warn = FALSE))) {
    stop("Existing profile signature components differ")
  }
  unlink(signature_tmp)
} else if (!file.rename(signature_tmp, signature_file)) {
  stop("Could not atomically place profile signature components")
}
signature_hash_file <- file.path(out_model_dir, "analysis_signature.sha256")
signature_hash_text <- paste(analysis_signature, basename(signature_file))
if (file.exists(signature_hash_file)) {
  if (!identical(readLines(signature_hash_file, warn = FALSE), signature_hash_text)) {
    stop("Existing profile signature hash differs")
  }
} else {
  writeLines(signature_hash_text, signature_hash_file, useBytes = TRUE)
}

fit_logcenter_full <- function(start, start_id) {
  run_signature <- hash_text(c(
    paste0("analysis_signature=", analysis_signature), "type=full",
    paste0("start_id=", start_id), paste0("label=", start$label),
    paste0("native=", fmt_num(start$native)),
    paste0("control=", fmt_num(unlist(full_control)))
  ))
  cp <- file.path(checkpoint_dir, sprintf("full_start%02d.rds", start_id))
  if (file.exists(cp)) {
    old <- readRDS(cp)
    if (!identical(old$analysis_signature, analysis_signature) ||
        !identical(old$run_signature, run_signature) ||
        !isTRUE(all.equal(old$native_start, start$native, tolerance = 0))) {
      stop("Full-fit checkpoint mismatch: ", cp)
    }
    cat(sprintf("[full %02d] resumed\n", start_id))
    return(old)
  }

  f_lamb <- function(t, x, y) exp(y[1])
  f_mu <- function(t, x, y) exp(y[1] + y[2] * (x - temperature_ref))
  seed <- 2026095000L + start_id
  set.seed(seed)
  warnings_out <- character(0)
  error_text <- NA_character_
  t0 <- proc.time()[["elapsed"]]
  fit <- tryCatch(
    withCallingHandlers(
      fit_env_bd_diagnostics_maxit_v3(
        phylo = tree, env_data = env_data, tot_time = total_time,
        f.lamb = f_lamb, f.mu = f_mu,
        lamb_par = start$native[1], mu_par = start$native[2:3],
        df = environment_df, f = sampling_fraction, meth = "Nelder-Mead",
        cst.lamb = TRUE, cst.mu = FALSE, expo.lamb = FALSE, expo.mu = FALSE,
        fix.mu = FALSE, cond = conditioning, control = full_control
      ),
      warning = function(w) {
        warnings_out <<- c(warnings_out, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      error_text <<- conditionMessage(e)
      NULL
    }
  )
  elapsed <- proc.time()[["elapsed"]] - t0
  native_final <- rep(NA_real_, 3L)
  names(native_final) <- c("log_lambda", "log_mu_ref", "beta")
  original_final <- rep(NA_real_, 5L)
  names(original_final) <- c("lambda0", "mu0", "beta", "log_mu_ref", "log_mu0")
  eq <- list(max_abs_rate_difference = NA_real_, max_abs_log_exponent_difference = NA_real_)
  rates <- NULL
  if (!is.null(fit) && is.finite(fit$LH)) {
    native_final <- setNames(as.numeric(fit$optim_par), names(native_final))
    original_final <- logcenter_to_original(native_final)
    eq <- equivalence_qa(native_final)
    rates <- rate_qa(native_final)
  }
  result <- list(
    analysis_signature = analysis_signature, run_signature = run_signature,
    type = "full", start_id = start_id, label = start$label,
    native_start = start$native, native_final = native_final,
    original_final = original_final, equivalence_qa = eq, rate_qa = rates,
    fit = fit, seed = seed, control = full_control,
    warnings = unique(c(warnings_out, if (!is.null(fit)) fit$warnings else character(0))),
    error = error_text, elapsed_sec = elapsed,
    tree_sha256 = tree_sha, environment_sha256 = env_sha,
    sampling_fraction = sampling_fraction, conditioning = conditioning,
    environment_df = environment_df, total_time = total_time
  )
  atomic_save_rds(result, cp)
  cat(sprintf("[full %02d] logLik=%s conv=%s beta=%s elapsed=%.1fs\n",
              start_id, if (is.null(fit)) "NA" else format(fit$LH, digits = 14),
              if (is.null(fit)) "NA" else fit$convergence,
              if (is.null(fit)) "NA" else format(native_final["beta"], digits = 8), elapsed))
  result
}

beta_tag <- function(beta) gsub("-", "m", gsub("\\.", "p", format(beta, scientific = FALSE)))
fit_profile_once <- function(beta, start, start_id) {
  run_signature <- hash_text(c(
    paste0("analysis_signature=", analysis_signature), "type=profile",
    paste0("beta=", sprintf("%.17g", beta)), paste0("start_id=", start_id),
    paste0("native=", fmt_num(start$native)),
    paste0("control=", fmt_num(unlist(profile_control)))
  ))
  cp <- file.path(checkpoint_dir, sprintf("profile_beta_%s_start%02d.rds", beta_tag(beta), start_id))
  if (file.exists(cp)) {
    old <- readRDS(cp)
    if (!identical(old$analysis_signature, analysis_signature) ||
        !identical(old$run_signature, run_signature) || !identical(old$fixed_beta, beta) ||
        !isTRUE(all.equal(old$native_start, start$native, tolerance = 0))) {
      stop("Profile checkpoint mismatch: ", cp)
    }
    cat(sprintf("[profile beta=%g start=%02d] resumed\n", beta, start_id))
    return(old)
  }

  f_lamb <- function(t, x, y) exp(y[1])
  f_mu <- function(t, x, y) exp(y[1] + beta * (x - temperature_ref))
  seed <- 2026096000L + as.integer(abs(beta) * 100) + start_id
  set.seed(seed)
  warnings_out <- character(0)
  error_text <- NA_character_
  t0 <- proc.time()[["elapsed"]]
  fit <- tryCatch(
    withCallingHandlers(
      fit_env_bd_diagnostics_maxit_v3(
        phylo = tree, env_data = env_data, tot_time = total_time,
        f.lamb = f_lamb, f.mu = f_mu,
        lamb_par = start$native[1], mu_par = start$native[2],
        df = environment_df, f = sampling_fraction, meth = "Nelder-Mead",
        cst.lamb = TRUE, cst.mu = FALSE, expo.lamb = FALSE, expo.mu = FALSE,
        fix.mu = FALSE, cond = conditioning, control = profile_control
      ),
      warning = function(w) {
        warnings_out <<- c(warnings_out, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      error_text <<- conditionMessage(e)
      NULL
    }
  )
  elapsed <- proc.time()[["elapsed"]] - t0
  native_final <- rep(NA_real_, 2L)
  names(native_final) <- c("log_lambda", "log_mu_ref")
  full_native <- c(log_lambda = NA_real_, log_mu_ref = NA_real_, beta = beta)
  original_final <- rep(NA_real_, 5L)
  names(original_final) <- c("lambda0", "mu0", "beta", "log_mu_ref", "log_mu0")
  eq <- list(max_abs_rate_difference = NA_real_, max_abs_log_exponent_difference = NA_real_)
  rates <- NULL
  if (!is.null(fit) && is.finite(fit$LH)) {
    native_final <- setNames(as.numeric(fit$optim_par), names(native_final))
    full_native <- c(native_final, beta = beta)
    original_final <- logcenter_to_original(full_native)
    eq <- equivalence_qa(full_native)
    rates <- rate_qa(full_native)
  }
  result <- list(
    analysis_signature = analysis_signature, run_signature = run_signature,
    type = "profile", fixed_beta = beta, start_id = start_id, label = start$label,
    native_start = start$native, native_final = native_final,
    full_native_final = full_native, original_final = original_final,
    equivalence_qa = eq, rate_qa = rates, fit = fit,
    seed = seed, control = profile_control,
    warnings = unique(c(warnings_out, if (!is.null(fit)) fit$warnings else character(0))),
    error = error_text, elapsed_sec = elapsed,
    tree_sha256 = tree_sha, environment_sha256 = env_sha,
    sampling_fraction = sampling_fraction, conditioning = conditioning,
    environment_df = environment_df, total_time = total_time
  )
  atomic_save_rds(result, cp)
  cat(sprintf("[profile beta=%g start=%02d] logLik=%s conv=%s elapsed=%.1fs\n",
              beta, start_id, if (is.null(fit)) "NA" else format(fit$LH, digits = 14),
              if (is.null(fit)) "NA" else fit$convergence, elapsed))
  result
}

eligible <- function(x) {
  !is.null(x$fit) && is.finite(x$fit$LH) && identical(as.integer(x$fit$convergence), 0L)
}
best_converged <- function(runs) {
  ok <- which(vapply(runs, eligible, logical(1)))
  if (!length(ok)) return(NA_integer_)
  ok[which.max(vapply(runs[ok], function(x) x$fit$LH, numeric(1)))]
}
best_finite <- function(runs) {
  ok <- which(vapply(runs, function(x) !is.null(x$fit) && is.finite(x$fit$LH), logical(1)))
  if (!length(ok)) return(NA_integer_)
  ok[which.max(vapply(runs[ok], function(x) x$fit$LH, numeric(1)))]
}
replication_count <- function(runs, best_i, tolerance = 1e-4) {
  if (is.na(best_i)) return(0L)
  best_ll <- runs[[best_i]]$fit$LH
  ok <- vapply(runs, function(x) eligible(x) && abs(x$fit$LH - best_ll) <= tolerance, logical(1))
  length(unique(vapply(runs[ok], function(x) canonical_start(x$native_start), character(1))))
}

requested_cores <- suppressWarnings(as.integer(Sys.getenv("RPANDA_PROFILE_CORES", "4")))
if (!is.finite(requested_cores) || requested_cores < 1L) requested_cores <- 1L
detected_cores <- detectCores(logical = FALSE)
if (!is.finite(detected_cores) || detected_cores < 1L) detected_cores <- 1L
cores <- min(requested_cores, detected_cores)
cat(sprintf(
  "BCSTDTempVar log-centered diagnostic: full_starts=%d beta_points=%d f=%.17g cond=%s df=%d Tref=%.17g signature=%s\n",
  length(full_starts), length(base_profile_betas), sampling_fraction,
  conditioning, environment_df, temperature_ref, analysis_signature
))

full_stage_cp <- file.path(checkpoint_dir, "full_multistart_complete.rds")
full_stage_rds <- file.path(out_model_dir, "full_multistart.rds")
if (file.exists(full_stage_cp)) {
  full_payload <- readRDS(full_stage_cp)
  if (!identical(full_payload$analysis_signature, analysis_signature)) stop("Full-stage checkpoint mismatch")
  full_runs <- full_payload$runs
  if (!file.exists(full_stage_rds)) atomic_save_rds(full_payload, full_stage_rds)
} else if (file.exists(full_stage_rds)) {
  full_payload <- readRDS(full_stage_rds)
  if (!identical(full_payload$analysis_signature, analysis_signature)) stop("Full-stage RDS mismatch")
  full_runs <- full_payload$runs
  atomic_save_rds(full_payload, full_stage_cp)
} else {
  if (.Platform$OS.type == "unix" && cores > 1L) {
    full_runs <- mclapply(seq_along(full_starts), function(i) {
      fit_logcenter_full(full_starts[[i]], i)
    }, mc.cores = min(cores, length(full_starts)), mc.preschedule = FALSE)
  } else {
    full_runs <- lapply(seq_along(full_starts), function(i) fit_logcenter_full(full_starts[[i]], i))
  }
  if (any(vapply(full_runs, inherits, logical(1), what = "try-error"))) stop("Full multistart worker failed")
  full_payload <- list(
    analysis_signature = analysis_signature, runs = full_runs,
    best_converged_i = best_converged(full_runs),
    completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  )
  atomic_save_rds(full_payload, full_stage_rds)
  atomic_save_rds(full_payload, full_stage_cp)
}

run_profile_beta <- function(beta) {
  beta_cp <- file.path(checkpoint_dir, sprintf("profile_beta_%s_complete.rds", beta_tag(beta)))
  beta_rds <- file.path(out_model_dir, sprintf("profile_beta_%s.rds", beta_tag(beta)))
  if (file.exists(beta_cp)) {
    payload <- readRDS(beta_cp)
    if (!identical(payload$analysis_signature, analysis_signature) || !identical(payload$beta, beta)) {
      stop("Profile-beta completion checkpoint mismatch: ", beta)
    }
    if (!file.exists(beta_rds)) atomic_save_rds(payload, beta_rds)
    return(payload)
  }
  if (file.exists(beta_rds)) {
    payload <- readRDS(beta_rds)
    if (!identical(payload$analysis_signature, analysis_signature) || !identical(payload$beta, beta)) {
      stop("Profile-beta model RDS mismatch: ", beta)
    }
    atomic_save_rds(payload, beta_cp)
    return(payload)
  }
  starts <- make_profile_starts(beta)
  runs <- lapply(seq_along(starts), function(i) fit_profile_once(beta, starts[[i]], i))
  payload <- list(
    analysis_signature = analysis_signature, beta = beta, runs = runs,
    best_converged_i = best_converged(runs), best_finite_i = best_finite(runs),
    completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  )
  atomic_save_rds(payload, beta_rds)
  atomic_save_rds(payload, beta_cp)
  payload
}

run_beta_set <- function(betas) {
  if (.Platform$OS.type == "unix" && cores > 1L) {
    out <- mclapply(as.list(betas), run_profile_beta,
                    mc.cores = min(cores, length(betas)), mc.preschedule = FALSE)
  } else {
    out <- lapply(as.list(betas), run_profile_beta)
  }
  if (any(vapply(out, inherits, logical(1), what = "try-error"))) stop("Profile-beta worker failed")
  out
}

profile_payloads <- run_beta_set(base_profile_betas)
profile_summary_from_payloads <- function(payloads) {
  rows <- lapply(payloads, function(p) {
    i <- p$best_converged_i
    j <- p$best_finite_i
    chosen <- if (!is.na(i)) p$runs[[i]] else if (!is.na(j)) p$runs[[j]] else NULL
    data.frame(
      beta = p$beta, any_converged = !is.na(i), selected_start_id = if (is.null(chosen)) NA_integer_ else chosen$start_id,
      selected_is_converged = if (is.null(chosen)) FALSE else eligible(chosen),
      logLik = if (is.null(chosen)) NA_real_ else chosen$fit$LH,
      log_lambda = if (is.null(chosen)) NA_real_ else chosen$full_native_final["log_lambda"],
      log_mu_ref = if (is.null(chosen)) NA_real_ else chosen$full_native_final["log_mu_ref"],
      lambda0 = if (is.null(chosen)) NA_real_ else chosen$original_final["lambda0"],
      mu0 = if (is.null(chosen)) NA_real_ else chosen$original_final["mu0"],
      log_mu0 = if (is.null(chosen)) NA_real_ else chosen$original_final["log_mu0"],
      min_mu = if (is.null(chosen) || is.null(chosen$rate_qa)) NA_real_ else chosen$rate_qa$min_mu,
      max_mu = if (is.null(chosen) || is.null(chosen$rate_qa)) NA_real_ else chosen$rate_qa$max_mu,
      equivalence_max_abs = if (is.null(chosen)) NA_real_ else chosen$equivalence_qa$max_abs_rate_difference,
      best_start_count = if (is.na(i)) 0L else replication_count(p$runs, i),
      stringsAsFactors = FALSE
    )
  })
  z <- do.call(rbind, rows)
  z[order(z$beta, decreasing = TRUE), ]
}

profile_summary <- profile_summary_from_payloads(profile_payloads)
tail_requires_extension <- function(summary, tolerance = 1e-4) {
  z <- summary[summary$selected_is_converged & is.finite(summary$logLik), ]
  if (nrow(z) < 3L) return(TRUE)
  z <- z[order(z$beta, decreasing = TRUE), ]
  tail <- tail(z, min(4L, nrow(z)))
  best_at_negative_edge <- which.max(z$logLik) == nrow(z)
  tail_non_decreasing <- all(diff(tail$logLik) >= -tolerance)
  best_at_negative_edge || tail_non_decreasing
}
if (tail_requires_extension(profile_summary)) {
  profile_payloads <- c(profile_payloads, run_beta_set(extension_betas_1))
  profile_summary <- profile_summary_from_payloads(profile_payloads)
}
if (tail_requires_extension(profile_summary)) {
  profile_payloads <- c(profile_payloads, run_beta_set(extension_betas_2))
  profile_summary <- profile_summary_from_payloads(profile_payloads)
}

full_best_i <- best_converged(full_runs)
full_finite_i <- best_finite(full_runs)
full_best <- if (!is.na(full_best_i)) full_runs[[full_best_i]] else NULL
full_best_ll <- if (is.null(full_best)) NA_real_ else full_best$fit$LH
full_replication <- replication_count(full_runs, full_best_i)
full_nonconv_ll <- vapply(full_runs, function(x) {
  if (!is.null(x$fit) && is.finite(x$fit$LH) && x$fit$convergence != 0L) x$fit$LH else NA_real_
}, numeric(1))
max_full_nonconv_ll <- if (any(is.finite(full_nonconv_ll))) max(full_nonconv_ll, na.rm = TRUE) else NA_real_
no_higher_full_nonconv <- is.finite(full_best_ll) &&
  (!is.finite(max_full_nonconv_ll) || full_best_ll + 1e-8 >= max_full_nonconv_ll)

valid_profile <- profile_summary$selected_is_converged & is.finite(profile_summary$logLik)
profile_valid <- profile_summary[valid_profile, ]
profile_valid <- profile_valid[order(profile_valid$beta, decreasing = TRUE), ]
all_profile_points_converged <- nrow(profile_valid) == nrow(profile_summary)
profile_best_index <- if (nrow(profile_valid)) which.max(profile_valid$logLik) else NA_integer_
profile_best_beta <- if (is.na(profile_best_index)) NA_real_ else profile_valid$beta[profile_best_index]
profile_best_ll <- if (is.na(profile_best_index)) NA_real_ else profile_valid$logLik[profile_best_index]
negative_edge_beta <- if (nrow(profile_valid)) min(profile_valid$beta) else NA_real_
negative_edge_ll <- if (nrow(profile_valid)) profile_valid$logLik[which.min(profile_valid$beta)] else NA_real_
tail <- if (nrow(profile_valid)) tail(profile_valid, min(4L, nrow(profile_valid))) else profile_valid
tail_deltas <- if (nrow(tail) >= 2L) diff(tail$logLik) else numeric(0)
tail_non_decreasing <- length(tail_deltas) > 0L && all(tail_deltas >= -1e-4)
tail_total_change <- if (nrow(tail) >= 2L) tail$logLik[nrow(tail)] - tail$logLik[1] else NA_real_
tail_flat <- length(tail_deltas) > 0L && max(abs(tail_deltas)) <= 1e-4
best_at_negative_edge <- !is.na(profile_best_index) && profile_best_index == nrow(profile_valid)
edge_indistinguishable_from_best <- is.finite(profile_best_ll) && is.finite(negative_edge_ll) &&
  profile_best_ll - negative_edge_ll <= 1e-4
profile_has_clear_interior_peak <- !is.na(profile_best_index) &&
  profile_best_index > 1L && profile_best_index < nrow(profile_valid) &&
  profile_best_ll - negative_edge_ll > 1e-3

full_equivalence <- if (is.null(full_best)) Inf else full_best$equivalence_qa$max_abs_rate_difference
full_rates_valid <- !is.null(full_best) && !is.null(full_best$rate_qa) &&
  isTRUE(full_best$rate_qa$finite) && isTRUE(full_best$rate_qa$nonnegative) &&
  !isTRUE(full_best$rate_qa$explosive)
boundary_evidence <- best_at_negative_edge || edge_indistinguishable_from_best ||
  tail_non_decreasing || tail_flat
finite_peak_supported <- all_profile_points_converged && profile_has_clear_interior_peak &&
  full_replication >= 2L && no_higher_full_nonconv && full_equivalence < 1e-12 &&
  full_rates_valid && is.finite(full_best_ll) && full_best_ll + 1e-4 >= profile_best_ll
diagnosis <- if (!all_profile_points_converged) {
  "inconclusive_profile_optimization"
} else if (boundary_evidence) {
  "boundary_nonidentifiable"
} else if (finite_peak_supported) {
  "finite_peak_supported"
} else {
  "inconclusive_no_boundary_or_replicated_finite_peak"
}

diagnosis_table <- data.frame(
  model = "BCSTDTempVar", diagnosis = diagnosis,
  finite_MLE_accepted = identical(diagnosis, "finite_peak_supported"),
  boundary_nonidentifiable = identical(diagnosis, "boundary_nonidentifiable"),
  full_best_logLik = full_best_ll,
  full_best_beta = if (is.null(full_best)) NA_real_ else full_best$native_final["beta"],
  full_best_log_mu_ref = if (is.null(full_best)) NA_real_ else full_best$native_final["log_mu_ref"],
  full_best_log_mu0 = if (is.null(full_best)) NA_real_ else full_best$original_final["log_mu0"],
  full_converged_reproduction_count = full_replication,
  maximum_finite_nonconverged_full_logLik = max_full_nonconv_ll,
  no_higher_finite_nonconverged_full_result = no_higher_full_nonconv,
  profile_best_beta = profile_best_beta, profile_best_logLik = profile_best_ll,
  most_negative_profile_beta = negative_edge_beta,
  most_negative_profile_logLik = negative_edge_ll,
  all_profile_points_converged = all_profile_points_converged,
  profile_has_clear_interior_peak = profile_has_clear_interior_peak,
  best_at_negative_grid_edge = best_at_negative_edge,
  negative_edge_within_1e_4_of_profile_best = edge_indistinguishable_from_best,
  last_four_profile_values_non_decreasing = tail_non_decreasing,
  last_four_profile_values_flat_within_1e_4 = tail_flat,
  last_four_profile_total_logLik_change = tail_total_change,
  logcenter_to_original_max_abs_rate_difference = full_equivalence,
  likelihood_statement = paste(
    "Same archived likelihood; log-centered coordinates are a one-to-one",
    "reparameterization for positive lambda0 and mu_ref"
  ),
  interpretation = if (identical(diagnosis, "boundary_nonidentifiable")) {
    paste(
      "The likelihood improves or remains indistinguishable toward beta -> -Inf;",
      "a finite beta estimate is not identified and must not be reported as a converged MLE."
    )
  } else if (identical(diagnosis, "finite_peak_supported")) {
    "The profile supports a replicated finite interior maximum under the tested model."
  } else {
    "The diagnostic is insufficient for a finite-MLE claim; retain the explicit inconclusive status."
  },
  analysis_signature = analysis_signature,
  stringsAsFactors = FALSE
)

full_rows <- lapply(full_runs, function(x) {
  data.frame(
    start_id = x$start_id, label = x$label,
    starting_log_lambda = x$native_start[1], starting_log_mu_ref = x$native_start[2],
    starting_beta = x$native_start[3],
    final_log_lambda = x$native_final[1], final_log_mu_ref = x$native_final[2],
    final_beta = x$native_final[3],
    final_lambda0 = x$original_final["lambda0"], final_mu0 = x$original_final["mu0"],
    final_log_mu0 = x$original_final["log_mu0"],
    logLik = if (is.null(x$fit)) NA_real_ else x$fit$LH,
    AICc = if (is.null(x$fit)) NA_real_ else x$fit$aicc,
    convergence = if (is.null(x$fit)) NA_integer_ else x$fit$convergence,
    optimizer_message = if (is.null(x$fit)) NA_character_ else x$fit$message,
    n_function = if (is.null(x$fit)) NA_integer_ else unname(x$fit$counts["function"]),
    elapsed_sec = x$elapsed_sec,
    warnings = if (length(x$warnings)) paste(x$warnings, collapse = " | ") else NA_character_,
    error = x$error,
    min_mu = if (is.null(x$rate_qa)) NA_real_ else x$rate_qa$min_mu,
    max_mu = if (is.null(x$rate_qa)) NA_real_ else x$rate_qa$max_mu,
    rate_explosion = if (is.null(x$rate_qa)) NA else x$rate_qa$explosive,
    equivalence_max_abs = x$equivalence_qa$max_abs_rate_difference,
    selected_finite_candidate = !is.na(full_best_i) && x$start_id == full_runs[[full_best_i]]$start_id,
    analysis_signature = analysis_signature, run_signature = x$run_signature,
    sampling_fraction = x$sampling_fraction, conditioning = x$conditioning,
    environment_df = x$environment_df, total_time = x$total_time,
    stringsAsFactors = FALSE
  )
})
full_table <- do.call(rbind, full_rows)

profile_run_rows <- unlist(lapply(profile_payloads, function(p) {
  lapply(p$runs, function(x) {
    data.frame(
      fixed_beta = x$fixed_beta, start_id = x$start_id, label = x$label,
      starting_log_lambda = x$native_start[1], starting_log_mu_ref = x$native_start[2],
      final_log_lambda = x$native_final[1], final_log_mu_ref = x$native_final[2],
      final_lambda0 = x$original_final["lambda0"], final_mu0 = x$original_final["mu0"],
      final_log_mu0 = x$original_final["log_mu0"],
      logLik = if (is.null(x$fit)) NA_real_ else x$fit$LH,
      convergence = if (is.null(x$fit)) NA_integer_ else x$fit$convergence,
      optimizer_message = if (is.null(x$fit)) NA_character_ else x$fit$message,
      n_function = if (is.null(x$fit)) NA_integer_ else unname(x$fit$counts["function"]),
      elapsed_sec = x$elapsed_sec,
      warnings = if (length(x$warnings)) paste(x$warnings, collapse = " | ") else NA_character_,
      error = x$error,
      min_mu = if (is.null(x$rate_qa)) NA_real_ else x$rate_qa$min_mu,
      max_mu = if (is.null(x$rate_qa)) NA_real_ else x$rate_qa$max_mu,
      rate_explosion = if (is.null(x$rate_qa)) NA else x$rate_qa$explosive,
      equivalence_max_abs = x$equivalence_qa$max_abs_rate_difference,
      selected_for_beta = !is.na(p$best_converged_i) && x$start_id == p$runs[[p$best_converged_i]]$start_id,
      analysis_signature = analysis_signature, run_signature = x$run_signature,
      stringsAsFactors = FALSE
    )
  })
}), recursive = FALSE)
profile_runs_table <- do.call(rbind, profile_run_rows)

write_clean_tsv(full_table, file.path(out_table_dir, "full_multistart_runs.tsv"))
write_clean_tsv(profile_runs_table, file.path(out_table_dir, "fixed_beta_profile_runs.tsv"))
write_clean_tsv(profile_summary, file.path(out_table_dir, "fixed_beta_profile_summary.tsv"))
write_clean_tsv(diagnosis_table, file.path(out_table_dir, "boundary_diagnosis.tsv"))

diagnostic_payload <- list(
  analysis_signature = analysis_signature,
  signature_components = signature_components,
  model = "BCSTDTempVar",
  parameterization = list(
    native = c("log_lambda", "log_mu_ref", "beta"),
    temperature_reference = temperature_ref,
    original_conversion = "lambda0=exp(log_lambda); mu0=exp(log_mu_ref-beta*Tref)",
    likelihood_changed = FALSE
  ),
  diagnosis = diagnosis, diagnosis_table = diagnosis_table,
  full_runs = full_runs, full_best_i = full_best_i,
  profile_payloads = profile_payloads, profile_summary = profile_summary,
  input_sha256 = c(tree = tree_sha, environment = env_sha, round1_rds = round1_sha),
  code_sha256 = c(driver = driver_sha, vendor_sha)
)
atomic_save_rds(diagnostic_payload, file.path(out_model_dir, "BCSTDTempVar_logcenter_profile.rds"))
if (identical(diagnosis, "finite_peak_supported")) {
  atomic_save_rds(list(
    status = diagnosis, fit = full_best$fit, native_parameters = full_best$native_final,
    original_parameters = full_best$original_final,
    logLik = full_best$fit$LH, reproduction_count = full_replication,
    analysis_signature = analysis_signature
  ), file.path(out_model_dir, "finite_mle_candidate_for_final_integration.rds"))
} else if (identical(diagnosis, "boundary_nonidentifiable")) {
  atomic_save_rds(list(
    status = diagnosis, finite_MLE = FALSE,
    instruction = "Do not report a finite beta MLE or silently treat a local optimizer convergence code as the global solution.",
    diagnosis_table = diagnosis_table, profile_summary = profile_summary,
    analysis_signature = analysis_signature
  ), file.path(out_model_dir, "boundary_nonidentifiable_diagnostic.rds"))
}
cat("BCSTDTempVar log-centered profile diagnosis: ", diagnosis,
    "; profile beta range ", min(profile_summary$beta), " to ", max(profile_summary$beta),
    "; signature ", analysis_signature, "\n", sep = "")
