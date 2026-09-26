#!/usr/bin/env Rscript

# Numerical-grid sensitivity diagnostic for the BCSTDTempVar boundary ridge.
# This script does not alter or replace the formal dt=0.005 fit.  It uses
# versioned copies of Phi/Psi whose only change is a configurable integration
# grid, and writes all results to independent v5 directories.

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
profile_seed_file <- file.path(
  run_dir, "03_models", "BCSTDTempVar_logcenter_profile_v4", "full_multistart.rds"
)
vendor_dir <- file.path(run_dir, "code", "vendor", "Toussaint2025_diagnostics")
out_model_dir <- file.path(run_dir, "03_models", "BCSTDTempVar_dt_sensitivity_v5")
out_table_dir <- file.path(run_dir, "04_tables", "BCSTDTempVar_dt_sensitivity_v5")
checkpoint_dir <- file.path(run_dir, "checkpoints", "03d_BCSTDTempVar_dt_sensitivity_v5")
dir.create(out_model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

vendor_names <- c(
  "integrate.R", "Phi.R", "Psi.R", "likelihood_bd.R",
  "fit_bd_diagnostics_maxit_v3.R", "Phi_Psi_dt_sensitivity_v1.R"
)
vendor_files <- file.path(vendor_dir, vendor_names)
needed <- c(tree_file, env_file, round1_file, profile_seed_file, vendor_files)
if (any(!file.exists(needed))) {
  stop("Missing dt-sensitivity input/source:\n", paste(needed[!file.exists(needed)], collapse = "\n"))
}

source(file.path(vendor_dir, "integrate.R"), local = .GlobalEnv)
source(file.path(vendor_dir, "Phi.R"), local = .GlobalEnv)
source(file.path(vendor_dir, "Psi.R"), local = .GlobalEnv)
archived_Phi <- .Phi
archived_Psi <- .Psi
source(file.path(vendor_dir, "likelihood_bd.R"), local = .GlobalEnv)
source(file.path(vendor_dir, "fit_bd_diagnostics_maxit_v3.R"), local = .GlobalEnv)
source(file.path(vendor_dir, "Phi_Psi_dt_sensitivity_v1.R"), local = .GlobalEnv)

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
  y <- gsub("[\\r\\n\\t]+", " ", as.character(x), perl = TRUE)
  trimws(gsub("[ ]{2,}", " ", y, perl = TRUE))
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
  tmp <- tempfile(pattern = paste0(basename(path), "_"), tmpdir = dirname(path), fileext = ".tmp")
  write.table(x, tmp, sep = "\t", quote = FALSE, row.names = FALSE,
              na = "NA", fileEncoding = "UTF-8")
  back <- read.delim(tmp, check.names = FALSE, stringsAsFactors = FALSE,
                     quote = "", comment.char = "", na.strings = "NA")
  ok <- nrow(back) == nrow(x) && ncol(back) == ncol(x) &&
    identical(names(back), names(x)) &&
    length(readLines(tmp, warn = FALSE)) == nrow(x) + 1L
  if (!ok) stop("TSV roundtrip failed: ", path)
  if (file.exists(path)) {
    if (!identical(readLines(path, warn = FALSE), readLines(tmp, warn = FALSE))) {
      stop("Refusing to overwrite non-identical dt-sensitivity TSV: ", path)
    }
    unlink(tmp)
  } else if (!file.rename(tmp, path)) {
    stop("Could not atomically place TSV: ", path)
  }
  invisible(TRUE)
}
atomic_save_rds <- function(object, path) {
  if (file.exists(path)) stop("Refusing to overwrite dt-sensitivity RDS: ", path)
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
official_dt <- 0.005
dt_values <- c(0.01, 0.005, 0.0025, 0.001)
optimizer_control <- list(maxit = 15000L, reltol = 1e-11, parscale = c(1, 10, 1))

tree_sha <- sha256(tree_file)
env_sha <- sha256(env_file)
round1_sha <- sha256(round1_file)
profile_seed_sha <- sha256(profile_seed_file)
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

# Reproduce the exact likelihood query times so the actual (not nominal) first
# cell widths can be reported for every Psi/Phi integration.
nbtips <- Ntip(tree)
from_past <- cbind(tree$edge, node.age(tree)$ages)
ages <- rbind(from_past[, 2:3], c(nbtips + 1, 0))
ages <- ages[order(ages[, 1]), ]
tree_age_for_queries <- max(ages[, 2])
psi_times <- vapply(seq_len(nbtips - 1L), function(j) {
  node <- nbtips + j
  edges <- tree$edge[tree$edge[, 1] == node, , drop = FALSE]
  tree_age_for_queries - ages[edges[1, 1], 2]
}, numeric(1))
likelihood_query_times <- c(psi_times, total_time)
if (any(!is.finite(likelihood_query_times)) || any(likelihood_query_times <= 0)) {
  stop("Could not reconstruct positive likelihood query times")
}
first_cell_width <- function(t, dt) t / (1L + as.integer(t / dt))

round1 <- readRDS(round1_file)
if (!all(c("analysis_signature", "total_time", "sampling_fraction", "df", "cond",
           "input_sha256") %in% names(round1))) {
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

profile_seed <- readRDS(profile_seed_file)
if (!all(c("analysis_signature", "runs") %in% names(profile_seed)) ||
    length(profile_seed$runs) < 4L) {
  stop("03c full-multistart RDS is incomplete")
}
seed_ok <- vapply(profile_seed$runs, function(x) {
  !is.null(x$fit) && is.finite(x$fit$LH) && identical(as.integer(x$fit$convergence), 0L) &&
    length(x$native_final) == 3L && all(is.finite(x$native_final)) &&
    identical(x$tree_sha256, tree_sha) && identical(x$environment_sha256, env_sha) &&
    identical(x$sampling_fraction, sampling_fraction) &&
    identical(x$conditioning, conditioning) &&
    identical(as.integer(x$environment_df), environment_df) &&
    identical(x$total_time, total_time)
}, logical(1))
seed_runs <- profile_seed$runs[seed_ok]
if (length(seed_runs) < 4L) stop("Fewer than four valid converged 03c seeds")
seed_order <- order(vapply(seed_runs, function(x) x$native_final[3], numeric(1)))
seed_runs <- seed_runs[seed_order]
seed_pick <- unique(as.integer(round(seq(1, length(seed_runs), length.out = 4L))))
if (length(seed_pick) != 4L) stop("Could not select four distinct deterministic seeds")
seed_runs <- seed_runs[seed_pick]

make_starts <- function(dt) {
  h_ref <- first_cell_width(total_time, official_dt)
  h_now <- first_cell_width(total_time, dt)
  lapply(seq_along(seed_runs), function(i) {
    p <- setNames(as.numeric(seed_runs[[i]]$native_final),
                  c("log_lambda", "log_mu_ref", "beta"))
    old_log_mu_present <- p["log_mu_ref"] + p["beta"] * (temperature_min - temperature_ref)
    target_log_mu_present <- old_log_mu_present + log(h_ref / h_now)
    p["log_mu_ref"] <- target_log_mu_present - p["beta"] * (temperature_min - temperature_ref)
    list(
      label = sprintf("dt_%g_scaled_03c_seed_%02d", dt, i),
      native = p,
      source_run_label = seed_runs[[i]]$label,
      source_logLik = seed_runs[[i]]$fit$LH,
      target_present_pulse_mass = exp(target_log_mu_present) * h_now
    )
  })
}
all_starts <- setNames(lapply(dt_values, make_starts), sprintf("%.17g", dt_values))
if (any(vapply(all_starts, function(z) {
  length(z) < 4L || length(unique(vapply(z, function(x) canonical_start(x$native), character(1)))) < 4L
}, logical(1)))) stop("Each dt must have at least four distinct starts")

# Prove that the copied general-environment formulas reproduce the archived
# implementation at the formal dt before any expensive optimization begins.
test_lamb <- function(t) rep(0.17, length(t))
test_mu <- function(t) 0.03 * exp(-0.02 * t)
copy_phi <- make_Phi_dt_sensitivity_v1(official_dt)
copy_psi <- make_Psi_dt_sensitivity_v1(official_dt)
copy_test_times <- c(1.25, 7.75, total_time)
copy_differences <- c(
  vapply(copy_test_times, function(t) abs(
    copy_phi(t, test_lamb, test_mu, sampling_fraction) -
      archived_Phi(t, test_lamb, test_mu, sampling_fraction)
  ), numeric(1)),
  vapply(copy_test_times, function(t) abs(
    copy_psi(0, t, test_lamb, test_mu, sampling_fraction) -
      archived_Psi(0, t, test_lamb, test_mu, sampling_fraction)
  ), numeric(1))
)
copy_equivalence_max_abs <- max(copy_differences)
if (!is.finite(copy_equivalence_max_abs) || copy_equivalence_max_abs >= 1e-12) {
  stop("Parameterized Phi/Psi copy does not reproduce archived dt=0.005 formulas")
}

signature_components <- data.frame(
  component = c(
    "signature_schema", "driver_sha256", "tree_sha256", "environment_sha256",
    "round1_rds_sha256", "round1_analysis_signature", "profile_seed_rds_sha256",
    "profile_seed_analysis_signature", "sampling_fraction", "conditioning",
    "environment_df", "total_time", "temperature_reference", "temperature_min",
    "temperature_max", "formal_official_dt", "diagnostic_dt_values",
    "parameterization", "optimizer", "optimizer_control", "parallelization",
    "copy_equivalence_max_abs", "formal_inference_rule", "R_version",
    paste0("package_version:", c("ape", "picante", "pspline")),
    paste0("vendor_sha256:", names(vendor_sha))
  ),
  value = c(
    "03d_BCSTDTempVar_dt_sensitivity_signature_1", driver_sha, tree_sha, env_sha,
    round1_sha, round1$analysis_signature, profile_seed_sha,
    profile_seed$analysis_signature, sprintf("%.17g", sampling_fraction), conditioning,
    as.character(environment_df), sprintf("%.17g", total_time),
    sprintf("%.17g", temperature_ref), sprintf("%.17g", temperature_min),
    sprintf("%.17g", temperature_max), sprintf("%.17g", official_dt),
    fmt_num(dt_values), "log_lambda|log_mu_ref|beta",
    "Nelder-Mead; archived likelihood with diagnostic-only parameterized Phi/Psi dt",
    fmt_num(unlist(optimizer_control)),
    "outer parallel across four dt values; starts sequential within dt",
    sprintf("%.17g", copy_equivalence_max_abs),
    "formal inference remains the archived dt=0.005 analysis; other dt values are numerical diagnostics only",
    R.version.string,
    vapply(c("ape", "picante", "pspline"), function(x) as.character(packageVersion(x)), character(1)),
    unname(vendor_sha)
  ), stringsAsFactors = FALSE
)
for (dt_name in names(all_starts)) {
  starts <- all_starts[[dt_name]]
  for (i in seq_along(starts)) {
    signature_components <- rbind(signature_components, data.frame(
      component = sprintf("start:dt_%s:%02d:%s", dt_name, i, starts[[i]]$label),
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
    stop("Existing dt-sensitivity signature components differ")
  }
  unlink(signature_tmp)
} else if (!file.rename(signature_tmp, signature_file)) {
  stop("Could not atomically place signature components")
}
signature_hash_file <- file.path(out_model_dir, "analysis_signature.sha256")
signature_hash_text <- paste(analysis_signature, basename(signature_file))
if (file.exists(signature_hash_file)) {
  if (!identical(readLines(signature_hash_file, warn = FALSE), signature_hash_text)) {
    stop("Existing dt-sensitivity signature hash differs")
  }
} else {
  writeLines(signature_hash_text, signature_hash_file, useBytes = TRUE)
}

dt_tag <- function(dt) gsub("-", "m", gsub("\\.", "p", format(dt, scientific = FALSE)))
rate_qa <- function(native) {
  log_lambda <- rep(native[1], length(time_grid))
  log_mu <- native[2] + native[3] * (temperature_grid - temperature_ref)
  lambda <- exp(log_lambda)
  mu <- exp(log_mu)
  list(
    finite = all(is.finite(lambda)) && all(is.finite(mu)),
    nonnegative = all(lambda >= 0) && all(mu >= 0),
    min_lambda = suppressWarnings(min(lambda)),
    max_lambda = suppressWarnings(max(lambda)),
    min_mu = suppressWarnings(min(mu)),
    max_mu = suppressWarnings(max(mu)),
    min_log_mu = suppressWarnings(min(log_mu)),
    max_log_mu = suppressWarnings(max(log_mu)),
    rate_gt_100 = suppressWarnings(max(c(lambda, mu))) > 100,
    rate_gt_1000 = suppressWarnings(max(c(lambda, mu))) > 1000
  )
}

fit_one <- function(dt, start, start_id) {
  run_signature <- hash_text(c(
    paste0("analysis_signature=", analysis_signature),
    paste0("dt=", sprintf("%.17g", dt)), paste0("start_id=", start_id),
    paste0("label=", start$label), paste0("native=", fmt_num(start$native)),
    paste0("control=", fmt_num(unlist(optimizer_control)))
  ))
  cp <- file.path(checkpoint_dir, sprintf("dt_%s_start%02d.rds", dt_tag(dt), start_id))
  if (file.exists(cp)) {
    old <- readRDS(cp)
    if (!identical(old$analysis_signature, analysis_signature) ||
        !identical(old$run_signature, run_signature) || !identical(old$dt, dt) ||
        !isTRUE(all.equal(old$native_start, start$native, tolerance = 0))) {
      stop("Run checkpoint mismatch: ", cp)
    }
    cat(sprintf("[dt=%g start=%02d] resumed\n", dt, start_id))
    return(old)
  }

  f_lamb <- function(t, x, y) exp(y[1])
  f_mu <- function(t, x, y) exp(y[1] + y[2] * (x - temperature_ref))
  seed <- 2026097000L + as.integer(dt * 1e6) + start_id
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
        fix.mu = FALSE, cond = conditioning, control = optimizer_control
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
  native_final <- setNames(rep(NA_real_, 3L), c("log_lambda", "log_mu_ref", "beta"))
  rates <- NULL
  log_mu0_original <- NA_real_
  mu0_original <- NA_real_
  mu_present <- NA_real_
  cell_widths <- first_cell_width(likelihood_query_times, dt)
  pulse_products <- rep(NA_real_, length(cell_widths))
  if (!is.null(fit) && is.finite(fit$LH) && length(fit$optim_par) == 3L &&
      all(is.finite(fit$optim_par))) {
    native_final <- setNames(as.numeric(fit$optim_par), names(native_final))
    rates <- rate_qa(native_final)
    log_mu0_original <- native_final["log_mu_ref"] - native_final["beta"] * temperature_ref
    mu0_original <- exp(log_mu0_original)
    mu_present <- exp(native_final["log_mu_ref"] +
                        native_final["beta"] * (temperature_min - temperature_ref))
    pulse_products <- mu_present * cell_widths
  }
  result <- list(
    analysis_signature = analysis_signature, run_signature = run_signature,
    dt = dt, start_id = start_id, label = start$label,
    source_run_label = start$source_run_label, source_logLik = start$source_logLik,
    native_start = start$native, native_final = native_final,
    log_mu0_original = log_mu0_original, mu0_original = mu0_original,
    mu_present = mu_present, rate_qa = rates,
    first_cell_width_total_time = first_cell_width(total_time, dt),
    first_cell_width_query_min = min(cell_widths),
    first_cell_width_query_median = median(cell_widths),
    first_cell_width_query_max = max(cell_widths),
    mu_present_times_first_cell_width_total = mu_present * first_cell_width(total_time, dt),
    pulse_product_query_min = suppressWarnings(min(pulse_products)),
    pulse_product_query_median = suppressWarnings(median(pulse_products)),
    pulse_product_query_max = suppressWarnings(max(pulse_products)),
    fit = fit, seed = seed, control = optimizer_control,
    warnings = unique(c(warnings_out, if (!is.null(fit)) fit$warnings else character(0))),
    error = error_text, elapsed_sec = elapsed,
    tree_sha256 = tree_sha, environment_sha256 = env_sha,
    profile_seed_sha256 = profile_seed_sha,
    sampling_fraction = sampling_fraction, conditioning = conditioning,
    environment_df = environment_df, total_time = total_time
  )
  atomic_save_rds(result, cp)
  cat(sprintf("[dt=%g start=%02d] logLik=%s conv=%s beta=%s mu_present=%s pulse=%s elapsed=%.1fs\n",
              dt, start_id, if (is.null(fit)) "NA" else format(fit$LH, digits = 14),
              if (is.null(fit)) "NA" else fit$convergence,
              if (is.null(fit)) "NA" else format(native_final["beta"], digits = 8),
              format(mu_present, digits = 8),
              format(mu_present * first_cell_width(total_time, dt), digits = 8), elapsed))
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

run_dt <- function(dt) {
  old_phi <- get(".Phi", envir = .GlobalEnv)
  old_psi <- get(".Psi", envir = .GlobalEnv)
  on.exit({
    assign(".Phi", old_phi, envir = .GlobalEnv)
    assign(".Psi", old_psi, envir = .GlobalEnv)
  }, add = TRUE)
  assign(".Phi", make_Phi_dt_sensitivity_v1(dt), envir = .GlobalEnv)
  assign(".Psi", make_Psi_dt_sensitivity_v1(dt), envir = .GlobalEnv)

  dt_cp <- file.path(checkpoint_dir, sprintf("dt_%s_complete.rds", dt_tag(dt)))
  dt_rds <- file.path(out_model_dir, sprintf("dt_%s_multistart.rds", dt_tag(dt)))
  if (file.exists(dt_cp)) {
    payload <- readRDS(dt_cp)
    if (!identical(payload$analysis_signature, analysis_signature) || !identical(payload$dt, dt)) {
      stop("dt completion checkpoint mismatch: ", dt)
    }
    if (!file.exists(dt_rds)) atomic_save_rds(payload, dt_rds)
    return(payload)
  }
  if (file.exists(dt_rds)) {
    payload <- readRDS(dt_rds)
    if (!identical(payload$analysis_signature, analysis_signature) || !identical(payload$dt, dt)) {
      stop("dt model RDS mismatch: ", dt)
    }
    atomic_save_rds(payload, dt_cp)
    return(payload)
  }

  starts <- make_starts(dt)
  runs <- lapply(seq_along(starts), function(i) fit_one(dt, starts[[i]], i))
  payload <- list(
    analysis_signature = analysis_signature, dt = dt, runs = runs,
    best_converged_i = best_converged(runs), best_finite_i = best_finite(runs),
    completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  )
  atomic_save_rds(payload, dt_rds)
  atomic_save_rds(payload, dt_cp)
  payload
}

requested_cores <- suppressWarnings(as.integer(Sys.getenv("RPANDA_DT_CORES", "4")))
if (!is.finite(requested_cores) || requested_cores < 1L) requested_cores <- 1L
detected_cores <- detectCores(logical = FALSE)
if (!is.finite(detected_cores) || detected_cores < 1L) detected_cores <- 1L
cores <- min(requested_cores, detected_cores, length(dt_values))
cat(sprintf(
  "BCSTDTempVar dt sensitivity: dt=%s starts_per_dt=4 cores=%d f=%.17g cond=%s df=%d signature=%s\n",
  fmt_num(dt_values), cores, sampling_fraction, conditioning, environment_df,
  analysis_signature
))

if (.Platform$OS.type == "unix" && cores > 1L) {
  dt_payloads <- mclapply(as.list(dt_values), run_dt, mc.cores = cores, mc.preschedule = FALSE)
} else {
  dt_payloads <- lapply(as.list(dt_values), run_dt)
}
if (any(vapply(dt_payloads, inherits, logical(1), what = "try-error"))) {
  stop("At least one dt worker failed")
}

run_rows <- unlist(lapply(dt_payloads, function(p) {
  lapply(p$runs, function(x) data.frame(
    dt = x$dt, start_id = x$start_id, label = x$label,
    source_run_label = x$source_run_label, source_logLik = x$source_logLik,
    starting_log_lambda = x$native_start[1], starting_log_mu_ref = x$native_start[2],
    starting_beta = x$native_start[3], final_log_lambda = x$native_final[1],
    final_log_mu_ref = x$native_final[2], final_beta = x$native_final[3],
    final_lambda = exp(x$native_final[1]), final_log_mu0_original = x$log_mu0_original,
    final_mu0_original = x$mu0_original, mu_present = x$mu_present,
    logLik = if (is.null(x$fit)) NA_real_ else x$fit$LH,
    AICc_diagnostic_only = if (is.null(x$fit)) NA_real_ else x$fit$aicc,
    convergence = if (is.null(x$fit)) NA_integer_ else x$fit$convergence,
    optimizer_message = if (is.null(x$fit)) NA_character_ else x$fit$message,
    n_function = if (is.null(x$fit)) NA_integer_ else unname(x$fit$counts["function"]),
    elapsed_sec = x$elapsed_sec,
    warnings = if (length(x$warnings)) paste(x$warnings, collapse = " | ") else NA_character_,
    error = x$error,
    min_lambda = if (is.null(x$rate_qa)) NA_real_ else x$rate_qa$min_lambda,
    max_lambda = if (is.null(x$rate_qa)) NA_real_ else x$rate_qa$max_lambda,
    min_mu = if (is.null(x$rate_qa)) NA_real_ else x$rate_qa$min_mu,
    max_mu = if (is.null(x$rate_qa)) NA_real_ else x$rate_qa$max_mu,
    min_log_mu = if (is.null(x$rate_qa)) NA_real_ else x$rate_qa$min_log_mu,
    max_log_mu = if (is.null(x$rate_qa)) NA_real_ else x$rate_qa$max_log_mu,
    all_rates_finite = if (is.null(x$rate_qa)) FALSE else x$rate_qa$finite,
    rate_gt_100 = if (is.null(x$rate_qa)) NA else x$rate_qa$rate_gt_100,
    rate_gt_1000 = if (is.null(x$rate_qa)) NA else x$rate_qa$rate_gt_1000,
    first_cell_width_total_time = x$first_cell_width_total_time,
    first_cell_width_query_min = x$first_cell_width_query_min,
    first_cell_width_query_median = x$first_cell_width_query_median,
    first_cell_width_query_max = x$first_cell_width_query_max,
    mu_present_times_first_cell_width_total = x$mu_present_times_first_cell_width_total,
    pulse_product_query_min = x$pulse_product_query_min,
    pulse_product_query_median = x$pulse_product_query_median,
    pulse_product_query_max = x$pulse_product_query_max,
    selected_for_dt = !is.na(p$best_converged_i) && x$start_id == p$runs[[p$best_converged_i]]$start_id,
    analysis_signature = analysis_signature, run_signature = x$run_signature,
    tree_sha256 = x$tree_sha256, environment_sha256 = x$environment_sha256,
    profile_seed_sha256 = x$profile_seed_sha256,
    sampling_fraction = x$sampling_fraction, conditioning = x$conditioning,
    environment_df = x$environment_df, total_time = x$total_time,
    stringsAsFactors = FALSE
  ))
}), recursive = FALSE)
runs_table <- do.call(rbind, run_rows)

summary_rows <- lapply(dt_payloads, function(p) {
  i <- p$best_converged_i
  j <- p$best_finite_i
  chosen <- if (!is.na(i)) p$runs[[i]] else if (!is.na(j)) p$runs[[j]] else NULL
  nonconv_ll <- vapply(p$runs, function(x) {
    if (!is.null(x$fit) && is.finite(x$fit$LH) && x$fit$convergence != 0L) x$fit$LH else NA_real_
  }, numeric(1))
  max_nonconv <- if (any(is.finite(nonconv_ll))) max(nonconv_ll, na.rm = TRUE) else NA_real_
  chosen_ll <- if (is.null(chosen) || is.null(chosen$fit)) NA_real_ else chosen$fit$LH
  data.frame(
    dt = p$dt, any_converged = !is.na(i),
    selected_start_id = if (is.null(chosen)) NA_integer_ else chosen$start_id,
    selected_is_converged = if (is.null(chosen)) FALSE else eligible(chosen),
    best_start_count = if (is.na(i)) 0L else replication_count(p$runs, i),
    logLik = chosen_ll,
    maximum_finite_nonconverged_logLik = max_nonconv,
    no_higher_finite_nonconverged = is.finite(chosen_ll) &&
      (!is.finite(max_nonconv) || chosen_ll + 1e-8 >= max_nonconv),
    log_lambda = if (is.null(chosen)) NA_real_ else chosen$native_final[1],
    log_mu_ref = if (is.null(chosen)) NA_real_ else chosen$native_final[2],
    beta = if (is.null(chosen)) NA_real_ else chosen$native_final[3],
    lambda = if (is.null(chosen)) NA_real_ else exp(chosen$native_final[1]),
    log_mu0_original = if (is.null(chosen)) NA_real_ else chosen$log_mu0_original,
    mu0_original = if (is.null(chosen)) NA_real_ else chosen$mu0_original,
    mu_present = if (is.null(chosen)) NA_real_ else chosen$mu_present,
    min_mu_observed_interval = if (is.null(chosen) || is.null(chosen$rate_qa)) NA_real_ else chosen$rate_qa$min_mu,
    max_mu_observed_interval = if (is.null(chosen) || is.null(chosen$rate_qa)) NA_real_ else chosen$rate_qa$max_mu,
    all_rates_finite = if (is.null(chosen) || is.null(chosen$rate_qa)) FALSE else chosen$rate_qa$finite,
    rate_gt_100 = if (is.null(chosen) || is.null(chosen$rate_qa)) NA else chosen$rate_qa$rate_gt_100,
    rate_gt_1000 = if (is.null(chosen) || is.null(chosen$rate_qa)) NA else chosen$rate_qa$rate_gt_1000,
    first_cell_width_total_time = if (is.null(chosen)) NA_real_ else chosen$first_cell_width_total_time,
    mu_present_times_first_cell_width_total = if (is.null(chosen)) NA_real_ else chosen$mu_present_times_first_cell_width_total,
    pulse_product_query_min = if (is.null(chosen)) NA_real_ else chosen$pulse_product_query_min,
    pulse_product_query_median = if (is.null(chosen)) NA_real_ else chosen$pulse_product_query_median,
    pulse_product_query_max = if (is.null(chosen)) NA_real_ else chosen$pulse_product_query_max,
    analysis_signature = analysis_signature,
    stringsAsFactors = FALSE
  )
})
summary_table <- do.call(rbind, summary_rows)
summary_table <- summary_table[order(summary_table$dt, decreasing = TRUE), ]

strict_optimization <- summary_table$selected_is_converged &
  summary_table$best_start_count >= 2L & summary_table$no_higher_finite_nonconverged
all_rates_finite <- all(summary_table$all_rates_finite)
logLik_range <- if (all(is.finite(summary_table$logLik))) diff(range(summary_table$logLik)) else Inf
mu_dt_slope <- if (all(is.finite(summary_table$mu_present)) && all(summary_table$mu_present > 0)) {
  unname(coef(lm(log(mu_present) ~ log(dt), data = summary_table))[2])
} else NA_real_
pulse_mass_ratio <- if (all(is.finite(summary_table$mu_present_times_first_cell_width_total)) &&
                        all(summary_table$mu_present_times_first_cell_width_total > 0)) {
  max(summary_table$mu_present_times_first_cell_width_total) /
    min(summary_table$mu_present_times_first_cell_width_total)
} else Inf
grid_ridge_supported <- all(strict_optimization) && all_rates_finite &&
  logLik_range <= 1e-3 && is.finite(mu_dt_slope) &&
  mu_dt_slope >= -1.25 && mu_dt_slope <= -0.75 && pulse_mass_ratio <= 1.20
diagnosis <- if (!all(strict_optimization)) {
  "inconclusive_dt_optimization"
} else if (grid_ridge_supported) {
  "endpoint_grid_ridge_supported"
} else {
  "dt_sensitive_without_complete_inverse_grid_scaling"
}

diagnosis_table <- data.frame(
  model = "BCSTDTempVar", diagnosis = diagnosis,
  endpoint_grid_ridge_supported = grid_ridge_supported,
  formal_inference_dt = official_dt,
  diagnostic_dt_values = fmt_num(dt_values),
  all_dt_strictly_optimized_and_replicated = all(strict_optimization),
  logLik_range_across_dt = logLik_range,
  expected_mu_present_loglog_slope = -1,
  observed_mu_present_loglog_slope = mu_dt_slope,
  present_pulse_mass_max_min_ratio = pulse_mass_ratio,
  all_selected_rates_finite = all_rates_finite,
  interpretation = if (grid_ridge_supported) {
    paste(
      "The fitted present-day extinction maximum scales approximately as 1/dt",
      "while its first-cell mass and likelihood remain stable. This supports an",
      "endpoint discretization ridge, not a biologically interpretable finite beta MLE."
    )
  } else if (identical(diagnosis, "inconclusive_dt_optimization")) {
    "At least one dt did not satisfy strict replicated-optimization checks; no grid-scaling claim is accepted."
  } else {
    "All dt fits passed optimization checks, but the prespecified inverse-grid scaling criteria were not all met."
  },
  formal_inference_statement = paste(
    "This is a numerical sensitivity analysis only. Formal inference and any",
    "model-status decision retain the archived Toussaint/RPANDA dt=0.005 likelihood."
  ),
  analysis_signature = analysis_signature,
  stringsAsFactors = FALSE
)

runs_path <- file.path(out_table_dir, "dt_sensitivity_runs.tsv")
summary_path <- file.path(out_table_dir, "dt_sensitivity_summary.tsv")
diagnosis_path <- file.path(out_table_dir, "dt_sensitivity_diagnosis.tsv")
write_clean_tsv(runs_table, runs_path)
write_clean_tsv(summary_table, summary_path)
write_clean_tsv(diagnosis_table, diagnosis_path)

qa_table <- data.frame(
  check = c(
    "exact_dt_grid_0.01_0.005_0.0025_0.001",
    "at_least_four_distinct_starts_per_dt",
    "all_dt_have_converged_solution",
    "all_dt_best_reproduced_at_least_twice_1e_4",
    "no_higher_finite_nonconverged_result_each_dt",
    "selected_rates_finite_and_nonnegative",
    "parameterized_copy_matches_archived_at_dt_0.005",
    "sampling_fraction_is_417_over_2100",
    "conditioning_is_crown", "environment_df_is_80",
    "tree_is_417_tip_rooted_binary", "tree_environment_seed_hashes_embedded",
    "tsv_roundtrip_runs", "tsv_roundtrip_summary", "tsv_roundtrip_diagnosis",
    "formal_inference_remains_official_dt_0.005",
    "dt_sensitivity_is_diagnostic_only"
  ),
  pass = c(
    identical(dt_values, c(0.01, 0.005, 0.0025, 0.001)),
    all(vapply(all_starts, function(z) length(unique(vapply(z, function(x) canonical_start(x$native), character(1)))) >= 4L, logical(1))),
    all(summary_table$selected_is_converged),
    all(summary_table$best_start_count >= 2L),
    all(summary_table$no_higher_finite_nonconverged),
    all(summary_table$all_rates_finite) && all(summary_table$min_mu_observed_interval >= 0),
    copy_equivalence_max_abs < 1e-12,
    identical(sampling_fraction, 417 / 2100), identical(conditioning, "crown"),
    identical(environment_df, 80L),
    Ntip(tree) == 417L && tree$Nnode == 416L && is.rooted(tree) && is.binary(tree),
    all(runs_table$tree_sha256 == tree_sha) && all(runs_table$environment_sha256 == env_sha) &&
      all(runs_table$profile_seed_sha256 == profile_seed_sha),
    TRUE, TRUE, TRUE, identical(official_dt, 0.005), TRUE
  ),
  detail = c(
    fmt_num(dt_values),
    paste(vapply(all_starts, length, integer(1)), collapse = ";"),
    paste(summary_table$selected_is_converged, collapse = ";"),
    paste(summary_table$best_start_count, collapse = ";"),
    paste(summary_table$no_higher_finite_nonconverged, collapse = ";"),
    paste(summary_table$all_rates_finite, collapse = ";"),
    sprintf("max_abs=%.17g", copy_equivalence_max_abs),
    sprintf("%.17g", sampling_fraction), conditioning, as.character(environment_df),
    sprintf("tips=%d nodes=%d", Ntip(tree), tree$Nnode),
    paste(tree_sha, env_sha, profile_seed_sha, sep = ";"),
    basename(runs_path), basename(summary_path), basename(diagnosis_path),
    sprintf("official_dt=%.17g", official_dt),
    "non-official dt values cannot replace the formal likelihood or enter AICc weights"
  ),
  analysis_signature = analysis_signature,
  stringsAsFactors = FALSE
)
qa_path <- file.path(out_table_dir, "dt_sensitivity_QA.tsv")
write_clean_tsv(qa_table, qa_path)

diagnostic_payload <- list(
  analysis_signature = analysis_signature,
  signature_components = signature_components,
  model = "BCSTDTempVar", diagnosis = diagnosis,
  formal_inference_dt = official_dt, diagnostic_dt_values = dt_values,
  formal_likelihood_changed = FALSE,
  parameterized_helper_copy_equivalence_max_abs = copy_equivalence_max_abs,
  dt_payloads = dt_payloads, summary = summary_table,
  diagnosis_table = diagnosis_table, qa = qa_table,
  input_sha256 = c(tree = tree_sha, environment = env_sha, round1_rds = round1_sha,
                   profile_seed_rds = profile_seed_sha),
  code_sha256 = c(driver = driver_sha, vendor_sha)
)
final_rds <- file.path(out_model_dir, "BCSTDTempVar_dt_sensitivity.rds")
if (!file.exists(final_rds)) atomic_save_rds(diagnostic_payload, final_rds)
if (grid_ridge_supported) {
  ridge_rds <- file.path(out_model_dir, "endpoint_grid_ridge_diagnostic.rds")
  if (!file.exists(ridge_rds)) atomic_save_rds(list(
    status = diagnosis, finite_MLE = FALSE,
    instruction = paste(
      "Treat the extreme-beta solution as a numerical endpoint-grid ridge;",
      "do not report it as a finite biological MLE or use non-official dt fits in AICc."
    ),
    diagnosis_table = diagnosis_table, summary = summary_table,
    analysis_signature = analysis_signature
  ), ridge_rds)
}
cat("BCSTDTempVar dt sensitivity diagnosis: ", diagnosis,
    "; LL range=", format(logLik_range, digits = 8),
    "; slope=", format(mu_dt_slope, digits = 8),
    "; pulse ratio=", format(pulse_mass_ratio, digits = 8),
    "; signature=", analysis_signature, "\n", sep = "")
