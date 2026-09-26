#!/usr/bin/env Rscript

# Formal 10-model RPANDA analysis for the Scenario A 417-tip tree.
# Version 2 strengthens checkpoint provenance and writes a model-level RDS
# immediately after every model. The likelihood implementation is unchanged.

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
vendor <- file.path(run_dir, "code", "vendor", "Toussaint2025_diagnostics")
checkpoint_dir <- file.path(run_dir, "checkpoints", "03_fit_models_v2")
model_dir <- file.path(run_dir, "03_models")
table_dir <- file.path(run_dir, "04_tables")
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

vendor_names <- c("integrate.R", "Phi.R", "Psi.R", "likelihood_bd.R",
                  "fit_bd_diagnostics.R", "fit_env_bd.R")
vendor_files <- file.path(vendor, vendor_names)
needed <- c(tree_file, env_file, vendor_files)
if (any(!file.exists(needed))) {
  stop("Missing required input/source:\n", paste(needed[!file.exists(needed)], collapse = "\n"))
}

source(file.path(vendor, "integrate.R"), local = .GlobalEnv)
source(file.path(vendor, "Phi.R"), local = .GlobalEnv)
source(file.path(vendor, "Psi.R"), local = .GlobalEnv)
source(file.path(vendor, "likelihood_bd.R"), local = .GlobalEnv)
source(file.path(vendor, "fit_bd_diagnostics.R"), local = .GlobalEnv)
source(file.path(vendor, "fit_env_bd.R"), local = .GlobalEnv)

tree <- read.tree(tree_file)
if (Ntip(tree) != 417L || tree$Nnode != 416L) stop("Formal tree is not 417-tip/416-node")
if (!is.rooted(tree) || !is.binary(tree)) stop("Formal tree must be rooted and fully binary")
if (any(!is.finite(tree$edge.length)) || any(tree$edge.length <= 0)) {
  stop("Tree has non-positive/non-finite branches")
}
tip_depths <- node.depth.edgelength(tree)[seq_len(Ntip(tree))]
if (diff(range(tip_depths)) >= 1e-5) stop("Root-to-tip range exceeds 1e-5 Ma")
total_time <- max(tip_depths)

env_all <- read.delim(env_file, check.names = FALSE, stringsAsFactors = FALSE)
if (ncol(env_all) < 2L) stop("Environment table must contain at least two columns")
env_data <- data.frame(time = as.numeric(env_all[[1]]), temp = as.numeric(env_all[[2]]))
if (anyNA(env_data) || any(!is.finite(as.matrix(env_data)))) stop("Invalid environment table")
if (anyDuplicated(env_data$time) || is.unsorted(env_data$time, strictly = TRUE)) {
  stop("Environment time must be strictly increasing and unique")
}
if (min(env_data$time) > 0 || max(env_data$time) < total_time) {
  stop("Environment does not cover tree age")
}

sampling_fraction <- 417 / 2100
conditioning <- "crown"
environment_df <- 80L
n_obs <- Ntip(tree)

# Used only to calibrate starts and to audit fitted rates. Each environmental
# fit still calls the archived fit_env_bd.R and independently fits df=80.
env_spline_qa <- sm.spline(env_data$time, env_data$temp, df = environment_df)
time_grid_qa <- sort(unique(c(
  seq(0, total_time, length.out = 10001L),
  env_data$time[env_data$time >= 0 & env_data$time <= total_time]
)))
temperature_grid_qa <- as.numeric(predict(env_spline_qa, time_grid_qa))
if (any(!is.finite(temperature_grid_qa))) stop("Non-finite df=80 temperature spline prediction")
temperature_ref <- stats::median(temperature_grid_qa)

rlam <- c(0.08, 0.12, 0.15, 0.20, 0.28, 0.10)
rmu <- c(0.001, 0.005, 0.010, 0.030, 0.080, 0.015)
time_slope <- c(-0.020, -0.010, 0.005, 0.010, 0.020, 0.001)
temp_slope <- c(-0.060, -0.030, 0.000, 0.030, 0.060, 0.010)

make_starts <- function(model) {
  switch(model,
    BCST = matrix(c(0.05, 0.08, 0.12, 0.16, 0.22, 0.30), ncol = 1L),
    BCSTDCST = cbind(rlam, rmu),
    BTimeVar = cbind(rlam, time_slope),
    BTimeVarDCST = cbind(rlam, time_slope, rmu),
    BCSTDTimeVar = cbind(rlam, rmu, rev(time_slope)),
    BTimeVarDTimeVar = cbind(rlam, time_slope, rmu, -time_slope),
    BTempVar = cbind(rlam * exp(-temp_slope * temperature_ref), temp_slope),
    BTempVarDCST = cbind(rlam * exp(-temp_slope * temperature_ref), temp_slope, rmu),
    BCSTDTempVar = cbind(rlam, rmu * exp(temp_slope * temperature_ref), -temp_slope),
    BTempVarDTempVar = cbind(
      rlam * exp(-temp_slope * temperature_ref), temp_slope,
      rmu * exp(temp_slope * temperature_ref), -temp_slope
    ),
    stop("Unknown model: ", model)
  )
}

specs <- list(
  BCST = list(family = "constant", lambda_type = "constant", mu_type = "zero", k = 1L,
              n_lamb = 1L, n_mu = 0L, cst_lamb = TRUE, cst_mu = TRUE,
              expo_lamb = FALSE, expo_mu = FALSE, fix_mu = TRUE,
              par_names = c("lambda0"),
              lambda_formula = "lambda(t) = lambda0", mu_formula = "mu(t) = 0"),
  BCSTDCST = list(family = "constant", lambda_type = "constant", mu_type = "constant", k = 2L,
                  n_lamb = 1L, n_mu = 1L, cst_lamb = TRUE, cst_mu = TRUE,
                  expo_lamb = FALSE, expo_mu = FALSE, fix_mu = FALSE,
                  par_names = c("lambda0", "mu0"),
                  lambda_formula = "lambda(t) = lambda0", mu_formula = "mu(t) = mu0"),
  BTimeVar = list(family = "time", lambda_type = "time", mu_type = "zero", k = 2L,
                  n_lamb = 2L, n_mu = 0L, cst_lamb = FALSE, cst_mu = TRUE,
                  expo_lamb = TRUE, expo_mu = FALSE, fix_mu = TRUE,
                  par_names = c("lambda0", "alpha"),
                  lambda_formula = "lambda(t) = lambda0 * exp(alpha * t)", mu_formula = "mu(t) = 0"),
  BTimeVarDCST = list(family = "time", lambda_type = "time", mu_type = "constant", k = 3L,
                      n_lamb = 2L, n_mu = 1L, cst_lamb = FALSE, cst_mu = TRUE,
                      expo_lamb = TRUE, expo_mu = FALSE, fix_mu = FALSE,
                      par_names = c("lambda0", "alpha", "mu0"),
                      lambda_formula = "lambda(t) = lambda0 * exp(alpha * t)", mu_formula = "mu(t) = mu0"),
  BCSTDTimeVar = list(family = "time", lambda_type = "constant", mu_type = "time", k = 3L,
                      n_lamb = 1L, n_mu = 2L, cst_lamb = TRUE, cst_mu = FALSE,
                      expo_lamb = FALSE, expo_mu = TRUE, fix_mu = FALSE,
                      par_names = c("lambda0", "mu0", "beta"),
                      lambda_formula = "lambda(t) = lambda0", mu_formula = "mu(t) = mu0 * exp(beta * t)"),
  BTimeVarDTimeVar = list(family = "time", lambda_type = "time", mu_type = "time", k = 4L,
                          n_lamb = 2L, n_mu = 2L, cst_lamb = FALSE, cst_mu = FALSE,
                          expo_lamb = TRUE, expo_mu = TRUE, fix_mu = FALSE,
                          par_names = c("lambda0", "alpha", "mu0", "beta"),
                          lambda_formula = "lambda(t) = lambda0 * exp(alpha * t)", mu_formula = "mu(t) = mu0 * exp(beta * t)"),
  BTempVar = list(family = "temperature", lambda_type = "temperature", mu_type = "zero", k = 2L,
                  n_lamb = 2L, n_mu = 0L, cst_lamb = FALSE, cst_mu = TRUE,
                  expo_lamb = FALSE, expo_mu = FALSE, fix_mu = TRUE,
                  par_names = c("lambda0", "alpha"),
                  lambda_formula = "lambda(t) = lambda0 * exp(alpha * T(t))", mu_formula = "mu(t) = 0"),
  BTempVarDCST = list(family = "temperature", lambda_type = "temperature", mu_type = "constant", k = 3L,
                      n_lamb = 2L, n_mu = 1L, cst_lamb = FALSE, cst_mu = TRUE,
                      expo_lamb = FALSE, expo_mu = FALSE, fix_mu = FALSE,
                      par_names = c("lambda0", "alpha", "mu0"),
                      lambda_formula = "lambda(t) = lambda0 * exp(alpha * T(t))", mu_formula = "mu(t) = mu0"),
  BCSTDTempVar = list(family = "temperature", lambda_type = "constant", mu_type = "temperature", k = 3L,
                      n_lamb = 1L, n_mu = 2L, cst_lamb = TRUE, cst_mu = FALSE,
                      expo_lamb = FALSE, expo_mu = FALSE, fix_mu = FALSE,
                      par_names = c("lambda0", "mu0", "beta"),
                      lambda_formula = "lambda(t) = lambda0", mu_formula = "mu(t) = mu0 * exp(beta * T(t))"),
  BTempVarDTempVar = list(family = "temperature", lambda_type = "temperature", mu_type = "temperature", k = 4L,
                          n_lamb = 2L, n_mu = 2L, cst_lamb = FALSE, cst_mu = FALSE,
                          expo_lamb = FALSE, expo_mu = FALSE, fix_mu = FALSE,
                          par_names = c("lambda0", "alpha", "mu0", "beta"),
                          lambda_formula = "lambda(t) = lambda0 * exp(alpha * T(t))", mu_formula = "mu(t) = mu0 * exp(beta * T(t))")
)

make_functions <- function(spec) {
  if (spec$family == "temperature") {
    f_lamb <- switch(spec$lambda_type,
      constant = function(t, x, y) y[1],
      temperature = function(t, x, y) y[1] * exp(y[2] * x))
    f_mu <- switch(spec$mu_type,
      zero = function(t, x, y) 0,
      constant = function(t, x, y) y[1],
      temperature = function(t, x, y) y[1] * exp(y[2] * x))
  } else {
    f_lamb <- switch(spec$lambda_type,
      constant = function(t, y) y[1],
      time = function(t, y) y[1] * exp(y[2] * t))
    f_mu <- switch(spec$mu_type,
      zero = function(t, y) 0,
      constant = function(t, y) y[1],
      time = function(t, y) y[1] * exp(y[2] * t))
  }
  list(lamb = f_lamb, mu = f_mu)
}

rate_qa <- function(spec, par) {
  names(par) <- spec$par_names
  lambda0 <- abs(unname(par["lambda0"]))
  alpha <- if ("alpha" %in% names(par)) unname(par["alpha"]) else NA_real_
  mu0 <- if ("mu0" %in% names(par)) abs(unname(par["mu0"])) else 0
  beta <- if ("beta" %in% names(par)) unname(par["beta"]) else NA_real_
  lambda <- switch(spec$lambda_type,
    constant = rep(lambda0, length(time_grid_qa)),
    time = lambda0 * exp(alpha * time_grid_qa),
    temperature = lambda0 * exp(alpha * temperature_grid_qa))
  mu <- switch(spec$mu_type,
    zero = rep(0, length(time_grid_qa)),
    constant = rep(mu0, length(time_grid_qa)),
    time = mu0 * exp(beta * time_grid_qa),
    temperature = mu0 * exp(beta * temperature_grid_qa))
  finite <- all(is.finite(lambda)) && all(is.finite(mu))
  nonnegative <- all(lambda >= 0) && all(mu >= 0)
  max_rate <- if (finite) max(c(lambda, mu)) else Inf
  flags <- character(0)
  if (lambda0 < 1e-7) flags <- c(flags, "lambda0_near_zero")
  if (spec$mu_type != "zero" && mu0 < 1e-7) flags <- c(flags, "mu0_near_zero")
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
driver_sha <- sha256(script_file)
tree_sha <- sha256(tree_file)
env_sha <- sha256(env_file)
vendor_sha <- setNames(vapply(vendor_files, sha256, character(1)), vendor_names)
package_sha <- vapply(c("ape", "picante", "pspline"), function(x) as.character(packageVersion(x)), character(1))

signature_components <- data.frame(
  component = c(
    "signature_schema", "driver_sha256", "tree_sha256", "environment_sha256",
    "sampling_fraction", "conditioning", "environment_df", "total_time",
    "temperature_reference", "R_version",
    paste0("package_version:", names(package_sha)),
    paste0("vendor_sha256:", names(vendor_sha))
  ),
  value = c(
    "03_fit_models_v2_signature_1", driver_sha, tree_sha, env_sha,
    sprintf("%.17g", sampling_fraction), conditioning, as.character(environment_df),
    sprintf("%.17g", total_time), sprintf("%.17g", temperature_ref), R.version.string,
    unname(package_sha), unname(vendor_sha)
  ),
  stringsAsFactors = FALSE
)

for (model in names(specs)) {
  spec <- specs[[model]]
  starts <- make_starts(model)
  colnames(starts) <- spec$par_names
  signature_components <- rbind(
    signature_components,
    data.frame(
      component = paste0("model_spec:", model),
      value = paste(
        spec$family, spec$lambda_type, spec$mu_type, spec$k, spec$n_lamb, spec$n_mu,
        spec$cst_lamb, spec$cst_mu, spec$expo_lamb, spec$expo_mu, spec$fix_mu,
        paste(spec$par_names, collapse = ","), spec$lambda_formula, spec$mu_formula,
        sep = "|"
      ),
      stringsAsFactors = FALSE
    )
  )
  for (i in seq_len(nrow(starts))) {
    signature_components <- rbind(
      signature_components,
      data.frame(component = sprintf("start:%s:%02d", model, i),
                 value = fmt_num(starts[i, ]), stringsAsFactors = FALSE)
    )
  }
}

signature_file <- file.path(model_dir, "analysis_signature_v2.tsv")
signature_tmp <- tempfile(pattern = "analysis_signature_", tmpdir = model_dir, fileext = ".tsv")
write.table(signature_components, signature_tmp, sep = "\t", quote = FALSE,
            row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
analysis_signature <- sha256(signature_tmp)
if (file.exists(signature_file)) {
  if (!identical(readLines(signature_file, warn = FALSE), readLines(signature_tmp, warn = FALSE))) {
    stop("Existing analysis signature components differ: ", signature_file)
  }
  unlink(signature_tmp)
} else if (!file.rename(signature_tmp, signature_file)) {
  stop("Could not atomically place analysis signature components")
}
signature_hash_file <- file.path(model_dir, "analysis_signature_v2.sha256")
signature_hash_text <- paste(analysis_signature, basename(signature_file))
if (file.exists(signature_hash_file)) {
  if (!identical(readLines(signature_hash_file, warn = FALSE), signature_hash_text)) {
    stop("Existing analysis signature hash differs: ", signature_hash_file)
  }
} else {
  writeLines(signature_hash_text, signature_hash_file, useBytes = TRUE)
}

run_signature <- function(model, start_id, start) {
  hash_text(c(
    paste0("analysis_signature=", analysis_signature),
    paste0("model=", model),
    paste0("start_id=", start_id),
    paste0("start=", fmt_num(start))
  ))
}

checkpoint_path <- function(model, start_id) {
  file.path(checkpoint_dir, sprintf("%s_start%02d.rds", model, start_id))
}
model_checkpoint_path <- function(model) file.path(checkpoint_dir, paste0(model, "_complete.rds"))
model_rds_path <- function(model) file.path(model_dir, paste0("model_fit_", model, ".rds"))

atomic_save_rds <- function(object, path) {
  if (file.exists(path)) stop("Refusing to overwrite existing RDS: ", path)
  tmp <- tempfile(pattern = paste0(basename(path), "_"), tmpdir = dirname(path), fileext = ".tmp")
  saveRDS(object, tmp, version = 3)
  if (!file.rename(tmp, path)) stop("Could not atomically place RDS: ", path)
  invisible(path)
}

validate_model_payload <- function(payload, model, path) {
  if (!is.list(payload) || !identical(payload$analysis_signature, analysis_signature) ||
      !identical(payload$model, model) || is.null(payload$result)) {
    stop("Model-level checkpoint signature/content mismatch: ", path)
  }
  invisible(TRUE)
}

fit_once <- function(model, spec, start, start_id) {
  cp <- checkpoint_path(model, start_id)
  expected_run_signature <- run_signature(model, start_id, start)
  if (file.exists(cp)) {
    old <- readRDS(cp)
    valid <- identical(old$analysis_signature, analysis_signature) &&
      identical(old$run_signature, expected_run_signature) &&
      identical(old$tree_sha256, tree_sha) && identical(old$env_sha256, env_sha) &&
      identical(old$driver_sha256, driver_sha) && identical(old$vendor_sha256, vendor_sha) &&
      identical(old$sampling_fraction, sampling_fraction) &&
      identical(old$conditioning, conditioning) && identical(old$environment_df, environment_df) &&
      identical(old$total_time, total_time) &&
      isTRUE(all.equal(as.numeric(old$start), as.numeric(start), tolerance = 0))
    if (!valid) stop("Checkpoint signature mismatch: ", cp)
    cat(sprintf("[%s] start %02d resumed; signature=%s\n", model, start_id, expected_run_signature))
    return(old)
  }

  seed <- 2026090400L + match(model, names(specs)) * 100L + start_id
  set.seed(seed)
  fun <- make_functions(spec)
  lamb_start <- start[seq_len(spec$n_lamb)]
  mu_start <- if (spec$n_mu) start[spec$n_lamb + seq_len(spec$n_mu)] else numeric(0)
  outside_warnings <- character(0)
  error_text <- NA_character_
  t0 <- proc.time()[["elapsed"]]
  fit <- tryCatch(
    withCallingHandlers({
      if (spec$family == "temperature") {
        fit_env_bd(
          phylo = tree, env_data = env_data, tot_time = total_time,
          f.lamb = fun$lamb, f.mu = fun$mu,
          lamb_par = lamb_start, mu_par = mu_start, df = environment_df,
          f = sampling_fraction, meth = "Nelder-Mead",
          cst.lamb = spec$cst_lamb, cst.mu = spec$cst_mu,
          expo.lamb = spec$expo_lamb, expo.mu = spec$expo_mu,
          fix.mu = spec$fix_mu, cond = conditioning
        )
      } else {
        fit_bd(
          phylo = tree, tot_time = total_time,
          f.lamb = fun$lamb, f.mu = fun$mu,
          lamb_par = lamb_start, mu_par = mu_start,
          f = sampling_fraction, meth = "Nelder-Mead",
          cst.lamb = spec$cst_lamb, cst.mu = spec$cst_mu,
          expo.lamb = spec$expo_lamb, expo.mu = spec$expo_mu,
          fix.mu = spec$fix_mu, cond = conditioning
        )
      }
    }, warning = function(w) {
      outside_warnings <<- c(outside_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }),
    error = function(e) {
      error_text <<- conditionMessage(e)
      NULL
    }
  )
  elapsed <- proc.time()[["elapsed"]] - t0

  qa <- NULL
  if (!is.null(fit) && is.finite(fit$LH)) qa <- rate_qa(spec, fit$optim_par)
  result <- list(
    model = model, start_id = start_id, seed = seed, method = "Nelder-Mead",
    start = as.numeric(start), fit = fit, rate_qa = qa,
    warnings = unique(c(outside_warnings, if (!is.null(fit)) fit$warnings else character(0))),
    error = error_text, elapsed_sec = elapsed,
    analysis_signature = analysis_signature, run_signature = expected_run_signature,
    tree_sha256 = tree_sha, env_sha256 = env_sha,
    driver_sha256 = driver_sha, vendor_sha256 = vendor_sha,
    total_time = total_time, sampling_fraction = sampling_fraction,
    conditioning = conditioning, environment_df = environment_df
  )
  atomic_save_rds(result, cp)
  cat(sprintf("[%s] start %02d complete: logLik=%s conv=%s elapsed=%.1fs signature=%s\n",
              model, start_id,
              if (is.null(fit)) "NA" else format(fit$LH, digits = 12),
              if (is.null(fit)) "NA" else fit$convergence, elapsed, expected_run_signature))
  result
}

eligible <- function(run) {
  !is.null(run$fit) && is.finite(run$fit$LH) && identical(as.integer(run$fit$convergence), 0L)
}

select_best <- function(runs) {
  ok <- which(vapply(runs, eligible, logical(1)))
  if (!length(ok)) {
    ok <- which(vapply(runs, function(x) !is.null(x$fit) && is.finite(x$fit$LH), logical(1)))
  }
  if (!length(ok)) return(NA_integer_)
  ok[which.max(vapply(runs[ok], function(x) x$fit$LH, numeric(1)))]
}

fit_model <- function(model) {
  complete_cp <- model_checkpoint_path(model)
  model_rds <- model_rds_path(model)
  if (file.exists(complete_cp)) {
    payload <- readRDS(complete_cp)
    validate_model_payload(payload, model, complete_cp)
    if (file.exists(model_rds)) {
      validate_model_payload(readRDS(model_rds), model, model_rds)
    } else {
      atomic_save_rds(payload, model_rds)
    }
    cat(sprintf("[%s] complete model resumed; analysis_signature=%s\n", model, analysis_signature))
    return(payload$result)
  }
  if (file.exists(model_rds)) {
    payload <- readRDS(model_rds)
    validate_model_payload(payload, model, model_rds)
    atomic_save_rds(payload, complete_cp)
    cat(sprintf("[%s] complete model recovered from model RDS; analysis_signature=%s\n", model, analysis_signature))
    return(payload$result)
  }

  spec <- specs[[model]]
  starts <- make_starts(model)
  colnames(starts) <- spec$par_names
  cat(sprintf("[%s] begin: f=%.17g cond=%s df=%d total_time=%.17g analysis_signature=%s\n",
              model, sampling_fraction, conditioning, environment_df, total_time, analysis_signature))
  runs <- lapply(seq_len(nrow(starts)), function(i) fit_once(model, spec, starts[i, ], i))

  best_i <- select_best(runs)
  if (!is.na(best_i)) {
    best_ll <- runs[[best_i]]$fit$LH
    replicated <- sum(vapply(runs, function(x) {
      eligible(x) && abs(x$fit$LH - best_ll) <= 1e-4
    }, logical(1)))
    if (replicated < 2L) {
      best_par <- as.numeric(runs[[best_i]]$fit$optim_par)
      rate_idx <- which(spec$par_names %in% c("lambda0", "mu0"))
      slope_idx <- which(spec$par_names %in% c("alpha", "beta"))
      mult <- c(1.02, 0.98, 1.05, 0.95)
      slope_shift <- c(0.0007, -0.0007, 0.0015, -0.0015)
      for (j in seq_along(mult)) {
        extra <- best_par
        extra[rate_idx] <- extra[rate_idx] * mult[j]
        extra[slope_idx] <- extra[slope_idx] + slope_shift[j]
        runs[[length(runs) + 1L]] <- fit_once(model, spec, extra, nrow(starts) + j)
        best_i <- select_best(runs)
        best_ll <- runs[[best_i]]$fit$LH
        replicated <- sum(vapply(runs, function(x) {
          eligible(x) && abs(x$fit$LH - best_ll) <= 1e-4
        }, logical(1)))
        if (replicated >= 2L) break
      }
    }
  }
  best_i <- select_best(runs)
  result <- list(model = model, spec = spec, runs = runs, best_i = best_i,
                 analysis_signature = analysis_signature)
  payload <- list(
    analysis_signature = analysis_signature, model = model, result = result,
    completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    configuration = list(
      tree_sha256 = tree_sha, environment_sha256 = env_sha,
      driver_sha256 = driver_sha, vendor_sha256 = vendor_sha,
      sampling_fraction = sampling_fraction, conditioning = conditioning,
      environment_df = environment_df, total_time = total_time
    )
  )
  atomic_save_rds(payload, model_rds)
  atomic_save_rds(payload, complete_cp)
  cat(sprintf("[%s] model-level RDS and completion checkpoint saved\n", model))
  result
}

model_names <- names(specs)
requested_cores <- suppressWarnings(as.integer(Sys.getenv("RPANDA_CORES", "4")))
if (!is.finite(requested_cores) || requested_cores < 1L) requested_cores <- 1L
detected_cores <- detectCores(logical = FALSE)
if (!is.finite(detected_cores) || detected_cores < 1L) detected_cores <- 1L
cores <- min(requested_cores, length(model_names), detected_cores)
cat(sprintf("Formal RPANDA v2: n=%d total_time=%.17g f=%.17g cond=%s df=%d cores=%d signature=%s\n",
            n_obs, total_time, sampling_fraction, conditioning, environment_df, cores,
            analysis_signature))

if (.Platform$OS.type == "unix" && cores > 1L) {
  results <- mclapply(model_names, fit_model, mc.cores = cores, mc.preschedule = FALSE)
} else {
  results <- lapply(model_names, fit_model)
}
names(results) <- model_names
worker_errors <- vapply(results, inherits, logical(1), what = "try-error")
if (any(worker_errors)) {
  stop("One or more parallel model workers failed: ",
       paste(model_names[worker_errors], collapse = ", "))
}

run_row <- function(x, selected = FALSE) {
  fit <- x$fit
  qa <- x$rate_qa
  data.frame(
    model = x$model, start_id = x$start_id, seed = x$seed, method = x$method,
    starting_values_raw = paste(format(x$start, digits = 17), collapse = ";"),
    final_values_raw = if (is.null(fit)) NA_character_ else paste(format(fit$optim_par, digits = 17), collapse = ";"),
    logLik = if (is.null(fit)) NA_real_ else fit$LH,
    AICc_from_vendor = if (is.null(fit)) NA_real_ else fit$aicc,
    convergence = if (is.null(fit)) NA_integer_ else fit$convergence,
    optimizer_message = if (is.null(fit)) NA_character_ else fit$message,
    n_function = if (is.null(fit) || is.null(fit$counts)) NA_integer_ else unname(fit$counts["function"]),
    n_gradient = if (is.null(fit) || is.null(fit$counts) || is.na(fit$counts["gradient"])) NA_integer_ else unname(fit$counts["gradient"]),
    elapsed_sec = x$elapsed_sec,
    warnings = if (length(x$warnings)) paste(x$warnings, collapse = " | ") else NA_character_,
    error = x$error,
    rates_finite = if (is.null(qa)) FALSE else qa$finite,
    rates_nonnegative = if (is.null(qa)) FALSE else qa$nonnegative,
    rate_explosion = if (is.null(qa)) NA else qa$explosive,
    min_lambda = if (is.null(qa)) NA_real_ else qa$min_lambda,
    max_lambda = if (is.null(qa)) NA_real_ else qa$max_lambda,
    min_mu = if (is.null(qa)) NA_real_ else qa$min_mu,
    max_mu = if (is.null(qa)) NA_real_ else qa$max_mu,
    boundary_flag = if (is.null(qa)) "fit_failed" else qa$boundary_flag,
    selected = selected,
    total_time = x$total_time, sampling_fraction = x$sampling_fraction,
    conditioning = x$conditioning, environment_df = x$environment_df,
    analysis_signature = x$analysis_signature, run_signature = x$run_signature,
    tree_sha256 = x$tree_sha256, environment_sha256 = x$env_sha256,
    driver_sha256 = x$driver_sha256,
    vendor_sha256 = paste(paste(names(x$vendor_sha256), x$vendor_sha256, sep = "="), collapse = ";"),
    stringsAsFactors = FALSE
  )
}

all_run_rows <- list()
selected_fits <- list()
table_rows <- list()
for (model in model_names) {
  z <- results[[model]]
  if (is.na(z$best_i)) {
    all_run_rows <- c(all_run_rows, lapply(z$runs, run_row, selected = FALSE))
    next
  }
  for (i in seq_along(z$runs)) {
    all_run_rows[[length(all_run_rows) + 1L]] <- run_row(z$runs[[i]], selected = i == z$best_i)
  }
  best <- z$runs[[z$best_i]]
  fit <- best$fit
  selected_fits[[model]] <- fit
  par <- rep(NA_real_, 4L)
  names(par) <- c("lambda0", "alpha", "mu0", "beta")
  raw <- fit$optim_par
  names(raw) <- z$spec$par_names
  par[names(raw)] <- raw
  par[c("lambda0", "mu0")] <- abs(par[c("lambda0", "mu0")])
  best_start_count <- sum(vapply(z$runs, function(x) {
    eligible(x) && abs(x$fit$LH - fit$LH) <= 1e-4
  }, logical(1)))
  notes <- c(if (fit$convergence == 0L) {
    "maximum finite converged logLik across deterministic starts"
  } else {
    "no converged finite fit; maximum finite nonconverged logLik retained for diagnosis"
  })
  if (best$rate_qa$boundary_flag != "none") notes <- c(notes, best$rate_qa$boundary_flag)
  table_rows[[model]] <- data.frame(
    model = model, family = z$spec$family,
    lambda_formula = z$spec$lambda_formula, mu_formula = z$spec$mu_formula,
    k = z$spec$k, logLik = fit$LH,
    AIC = -2 * fit$LH + 2 * z$spec$k,
    AICc = -2 * fit$LH + 2 * z$spec$k +
      (2 * z$spec$k * (z$spec$k + 1)) / (n_obs - z$spec$k - 1),
    lambda0 = par["lambda0"], alpha = par["alpha"],
    mu0 = par["mu0"], beta = par["beta"],
    convergence = fit$convergence, boundary_flag = best$rate_qa$boundary_flag,
    best_start_count = best_start_count, notes = paste(notes, collapse = "; "),
    stringsAsFactors = FALSE
  )
}

optimization_runs <- do.call(rbind, all_run_rows)
write.table(optimization_runs, file.path(model_dir, "optimization_runs.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")

if (length(table_rows) != length(model_names)) {
  incomplete <- list(
    analysis_signature = analysis_signature, results = results,
    selected_fits = selected_fits, tree = tree, env = env_data,
    total_time = total_time, sampling_fraction = sampling_fraction,
    df = environment_df, cond = conditioning
  )
  incomplete_path <- file.path(model_dir, "model_fits_10models.INCOMPLETE.v2.rds")
  if (!file.exists(incomplete_path)) atomic_save_rds(incomplete, incomplete_path)
  stop("At least one prespecified model has no finite fit; see optimization_runs.tsv")
}

model_table <- do.call(rbind, table_rows)
rownames(model_table) <- NULL
model_table$deltaAICc <- model_table$AICc - min(model_table$AICc)
rel <- exp(-0.5 * model_table$deltaAICc)
model_table$Akaike_weight <- rel / sum(rel)
model_table <- model_table[order(model_table$AICc, model_table$model), ]
rownames(model_table) <- NULL
front <- c("model", "lambda_formula", "mu_formula", "k", "logLik", "AIC", "AICc",
           "deltaAICc", "Akaike_weight", "lambda0", "alpha", "mu0", "beta",
           "convergence", "boundary_flag", "best_start_count", "notes", "family")
model_table <- model_table[front]
write.table(model_table, file.path(table_dir, "model_table_10models_full_precision.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")

formatted <- model_table
for (nm in c("logLik", "AIC", "AICc", "deltaAICc")) formatted[[nm]] <- sprintf("%.3f", formatted[[nm]])
formatted$Akaike_weight <- sprintf("%.4f", formatted$Akaike_weight)
for (nm in c("lambda0", "alpha", "mu0", "beta")) {
  formatted[[nm]] <- ifelse(is.na(formatted[[nm]]), "—", sprintf("%.6g", formatted[[nm]]))
}
write.table(formatted, file.path(table_dir, "model_table_10models_formatted.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")

same_tree_all_runs <- all(optimization_runs$tree_sha256 == tree_sha)
same_env_all_runs <- all(optimization_runs$environment_sha256 == env_sha)
same_f_all_runs <- all(optimization_runs$sampling_fraction == sampling_fraction)
same_cond_all_runs <- all(optimization_runs$conditioning == conditioning)
same_df_all_runs <- all(optimization_runs$environment_df == environment_df)
same_time_all_runs <- all(optimization_runs$total_time == total_time)
same_signature_all_runs <- all(optimization_runs$analysis_signature == analysis_signature)
selected_rows <- optimization_runs$selected
multi_models <- names(specs)[vapply(specs, function(x) x$k > 1L, logical(1))]
unique_starts <- vapply(multi_models, function(m) {
  length(unique(optimization_runs$starting_values_raw[optimization_runs$model == m &
                                                        optimization_runs$start_id <= 6L]))
}, integer(1))
run_k <- vapply(optimization_runs$model, function(x) specs[[x]]$k, integer(1))
finite_fit_rows <- is.finite(optimization_runs$logLik) & is.finite(optimization_runs$AICc_from_vendor)
recomputed_run_aicc <- -2 * optimization_runs$logLik + 2 * run_k +
  (2 * run_k * (run_k + 1)) / (n_obs - run_k - 1)
vendor_aicc_agrees <- all(abs(optimization_runs$AICc_from_vendor[finite_fit_rows] -
                               recomputed_run_aicc[finite_fit_rows]) < 1e-10)

qa <- data.frame(
  check = c(
    "prespecified_models_fitted", "all_selected_convergence_zero",
    "all_selected_rates_finite", "all_selected_rates_nonnegative",
    "no_selected_rate_explosion_gt_100_per_Ma",
    "all_models_best_reproduced_twice_1e-4", "at_least_5_unique_starts_each_multiparameter_model",
    "akaike_weights_sum", "same_tree_all_runs", "same_environment_all_runs",
    "same_sampling_fraction_all_runs", "same_conditioning_all_runs",
    "same_temperature_df_all_runs", "same_total_time_all_runs",
    "same_full_analysis_signature_all_runs", "vendor_AICc_matches_full_precision_recalculation",
    "unrounded_logLik_used_for_AICc"
  ),
  observed = c(
    nrow(model_table), all(model_table$convergence == 0L),
    all(optimization_runs$rates_finite[selected_rows]),
    all(optimization_runs$rates_nonnegative[selected_rows]),
    all(!optimization_runs$rate_explosion[selected_rows]),
    all(model_table$best_start_count >= 2L), paste(unique_starts, collapse = ";"),
    format(sum(model_table$Akaike_weight), digits = 17), same_tree_all_runs,
    same_env_all_runs, same_f_all_runs, same_cond_all_runs, same_df_all_runs,
    same_time_all_runs, same_signature_all_runs, vendor_aicc_agrees, TRUE
  ),
  expected = c(
    "10", "TRUE", "TRUE", "TRUE", "TRUE", "TRUE", ">=5 for every k>1 model",
    "1 within floating-point error", "TRUE", "TRUE", "TRUE", "TRUE", "TRUE",
    "TRUE", "TRUE", "TRUE", "TRUE"
  ),
  pass = c(
    nrow(model_table) == 10L, all(model_table$convergence == 0L),
    all(optimization_runs$rates_finite[selected_rows]),
    all(optimization_runs$rates_nonnegative[selected_rows]),
    all(!optimization_runs$rate_explosion[selected_rows]),
    all(model_table$best_start_count >= 2L), all(unique_starts >= 5L),
    abs(sum(model_table$Akaike_weight) - 1) < 1e-12,
    same_tree_all_runs, same_env_all_runs, same_f_all_runs, same_cond_all_runs,
    same_df_all_runs, same_time_all_runs, same_signature_all_runs,
    vendor_aicc_agrees, TRUE
  ),
  stringsAsFactors = FALSE
)
write.table(qa, file.path(table_dir, "model_selection_QA.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")

final_payload <- list(
  analysis_signature = analysis_signature,
  signature_components = signature_components,
  selected_fits = selected_fits, all_results = results, model_table = model_table,
  tree = tree, env = env_data, total_time = total_time,
  sampling_fraction = sampling_fraction, df = environment_df, cond = conditioning,
  temperature_reference_for_starts = temperature_ref,
  input_sha256 = c(tree = tree_sha, environment = env_sha),
  code_sha256 = c(driver = driver_sha, vendor_sha)
)
final_rds <- file.path(model_dir, "model_fits_10models.rds")
if (file.exists(final_rds)) {
  old_final <- readRDS(final_rds)
  if (!identical(old_final$analysis_signature, analysis_signature)) {
    stop("Existing final RDS has a different analysis signature: ", final_rds)
  }
} else {
  atomic_save_rds(final_payload, final_rds)
}

if (!all(qa$pass)) stop("One or more formal model-selection QA checks failed")
cat("All 10 models complete and QA-passing. Best model: ", model_table$model[1],
    "; deltaAICc(second): ", format(model_table$deltaAICc[2], digits = 7),
    "; analysis_signature: ", analysis_signature, "\n", sep = "")

