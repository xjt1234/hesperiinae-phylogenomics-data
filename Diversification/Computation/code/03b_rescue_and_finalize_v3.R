#!/usr/bin/env Rscript

# Rescue BCSTDTempVar with long Nelder-Mead optimization, then combine it with
# the nine accepted first-round fits. All v3 outputs are written to new
# final_v3 directories; first-round artifacts are read-only.

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
final_model_dir <- file.path(run_dir, "03_models", "final_v3")
final_table_dir <- file.path(run_dir, "04_tables", "final_v3")
checkpoint_dir <- file.path(run_dir, "checkpoints", "03b_rescue_v3")
dir.create(final_model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(final_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

vendor_names <- c(
  "integrate.R", "Phi.R", "Psi.R", "likelihood_bd.R",
  "fit_bd_diagnostics_maxit_v3.R"
)
vendor_files <- file.path(vendor_dir, vendor_names)
needed <- c(tree_file, env_file, round1_file, vendor_files)
if (any(!file.exists(needed))) {
  stop("Missing required v3 input/source:\n", paste(needed[!file.exists(needed)], collapse = "\n"))
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
has_control_text <- function(x) {
  any(vapply(x, function(v) {
    if (!(is.character(v) || is.factor(v))) return(FALSE)
    any(grepl("[\\r\\n\\t]", as.character(v), perl = TRUE), na.rm = TRUE)
  }, logical(1)))
}

preflight_tsv <- function(x) {
  x <- clean_frame(x)
  if (has_control_text(x)) return(FALSE)
  tmp <- tempfile(pattern = "tsv_preflight_", tmpdir = checkpoint_dir, fileext = ".tsv")
  on.exit(unlink(tmp), add = TRUE)
  write.table(x, tmp, sep = "\t", quote = FALSE, row.names = FALSE,
              na = "NA", fileEncoding = "UTF-8")
  lines_ok <- length(readLines(tmp, warn = FALSE)) == nrow(x) + 1L
  back <- tryCatch(
    read.delim(tmp, check.names = FALSE, stringsAsFactors = FALSE,
               quote = "", comment.char = "", na.strings = "NA"),
    error = function(e) NULL
  )
  lines_ok && !is.null(back) && nrow(back) == nrow(x) &&
    ncol(back) == ncol(x) && identical(names(back), names(x))
}

write_clean_tsv <- function(x, path) {
  x <- clean_frame(x)
  if (!preflight_tsv(x)) stop("TSV preflight/roundtrip failed: ", path)
  if (file.exists(path)) stop("Refusing to overwrite existing v3 TSV: ", path)
  tmp <- tempfile(pattern = paste0(basename(path), "_"), tmpdir = dirname(path), fileext = ".tmp")
  write.table(x, tmp, sep = "\t", quote = FALSE, row.names = FALSE,
              na = "NA", fileEncoding = "UTF-8")
  if (!file.rename(tmp, path)) stop("Could not atomically place TSV: ", path)
  back <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE,
                     quote = "", comment.char = "", na.strings = "NA")
  if (nrow(back) != nrow(x) || ncol(back) != ncol(x) ||
      length(readLines(path, warn = FALSE)) != nrow(x) + 1L) {
    stop("Written TSV failed roundtrip: ", path)
  }
  invisible(TRUE)
}

atomic_save_rds <- function(object, path) {
  if (file.exists(path)) stop("Refusing to overwrite existing v3 RDS: ", path)
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
if (ncol(env_all) < 2L) stop("Environment table must contain at least two columns")
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
  seq(0, total_time, length.out = 10001L),
  env_data$time[env_data$time >= 0 & env_data$time <= total_time]
)))
temperature_grid <- as.numeric(predict(env_spline, time_grid))
if (any(!is.finite(temperature_grid))) stop("Non-finite df=80 spline prediction")
temperature_ref <- median(temperature_grid)

model_names <- c(
  "BCST", "BCSTDCST", "BTimeVar", "BTimeVarDCST", "BCSTDTimeVar",
  "BTimeVarDTimeVar", "BTempVar", "BTempVarDCST", "BCSTDTempVar",
  "BTempVarDTempVar"
)
rescue_model <- "BCSTDTempVar"
carry_models <- setdiff(model_names, rescue_model)

round1 <- readRDS(round1_file)
required_round1 <- c("analysis_signature", "selected_fits", "all_results", "model_table",
                     "total_time", "sampling_fraction", "df", "cond", "input_sha256")
if (!all(required_round1 %in% names(round1))) stop("First-round RDS lacks required fields")
if (!setequal(names(round1$all_results), model_names) || nrow(round1$model_table) != 10L) {
  stop("First-round RDS does not contain the prescribed 10 models")
}
if (!identical(round1$sampling_fraction, sampling_fraction) ||
    !identical(round1$cond, conditioning) || !identical(as.integer(round1$df), environment_df) ||
    !identical(round1$total_time, total_time)) {
  stop("First-round f/cond/df/total_time does not match v3 configuration")
}
if (!identical(unname(round1$input_sha256[c("tree", "environment")]),
               unname(c(tree_sha, env_sha)))) {
  stop("First-round tree/environment SHA-256 mismatch")
}

for (m in model_names) {
  payload_file <- file.path(run_dir, "03_models", paste0("model_fit_", m, ".rds"))
  if (!file.exists(payload_file)) stop("Missing first-round model RDS: ", payload_file)
  payload <- readRDS(payload_file)
  if (!identical(payload$analysis_signature, round1$analysis_signature) ||
      !identical(payload$model, m)) stop("First-round model payload mismatch: ", m)
}
carry_table <- round1$model_table[match(carry_models, round1$model_table$model), , drop = FALSE]
if (anyNA(carry_table$model) || any(carry_table$convergence != 0L) ||
    any(carry_table$best_start_count < 2L)) {
  stop("One of the nine carried first-round models is not convergence/reproduction qualified")
}

old_rescue <- round1$all_results[[rescue_model]]
old_finite_nonconv <- which(vapply(old_rescue$runs, function(x) {
  !is.null(x$fit) && is.finite(x$fit$LH) && x$fit$convergence != 0L
}, logical(1)))
if (length(old_finite_nonconv) < 1L) stop("No finite nonconverged BCSTDTempVar endpoint found")
endpoint_order <- old_finite_nonconv[order(vapply(old_rescue$runs[old_finite_nonconv],
                                                  function(x) x$start_id, integer(1)))]
endpoint_raw <- lapply(old_rescue$runs[endpoint_order], function(x) {
  p <- as.numeric(x$fit$optim_par)
  names(p) <- c("lambda0", "mu0", "beta")
  c(lambda0 = abs(p[1]), mu0 = abs(p[2]), beta = p[3])
})
endpoint_labels <- paste0("round1_nonconverged_endpoint_start", sprintf("%02d", endpoint_order))
names(endpoint_raw) <- endpoint_labels
endpoint_ll <- vapply(old_rescue$runs[endpoint_order], function(x) x$fit$LH, numeric(1))
top_raw <- endpoint_raw[[which.max(endpoint_ll)]]

raw_to_centered <- function(raw) {
  c(lambda0 = abs(raw[1]), mu_ref = abs(raw[2]) * exp(raw[3] * temperature_ref), beta = raw[3])
}
centered_to_raw <- function(centered) {
  c(lambda0 = abs(centered[1]),
    mu0 = abs(centered[2]) * exp(-centered[3] * temperature_ref),
    beta = centered[3])
}
rate_equivalence <- function(centered, raw = centered_to_raw(centered)) {
  mu_center <- abs(centered[2] * exp(centered[3] * (temperature_grid - temperature_ref)))
  mu_raw <- abs(raw[2] * exp(raw[3] * temperature_grid))
  max(abs(mu_center - mu_raw))
}

make_start <- function(label, parameterization, native_start) {
  native_start <- as.numeric(native_start)
  names(native_start) <- if (parameterization == "centered") {
    c("lambda0", "mu_ref", "beta")
  } else {
    c("lambda0", "mu0", "beta")
  }
  raw_start <- if (parameterization == "centered") centered_to_raw(native_start) else {
    c(lambda0 = abs(native_start[1]), mu0 = abs(native_start[2]), beta = native_start[3])
  }
  list(label = label, parameterization = parameterization,
       native_start = native_start, raw_start = raw_start)
}

rescue_starts <- list()
for (i in seq_along(endpoint_raw)) {
  rescue_starts[[length(rescue_starts) + 1L]] <- make_start(
    paste0(names(endpoint_raw)[i], "_centered"), "centered", raw_to_centered(endpoint_raw[[i]]))
}
base_center <- raw_to_centered(top_raw)
center_perturb <- rbind(
  c(0.98, 0.50, -0.050), c(1.02, 2.00, 0.050),
  c(1.00, 0.25, -0.100), c(1.00, 4.00, 0.100),
  c(0.95, 1.00, -0.020), c(1.05, 1.00, 0.020),
  c(1.00, 1.50, -0.075)
)
colnames(center_perturb) <- c("lambda_mult", "mu_mult", "beta_delta")
for (i in seq_len(nrow(center_perturb))) {
  p <- base_center
  p[1] <- p[1] * center_perturb[i, "lambda_mult"]
  p[2] <- p[2] * center_perturb[i, "mu_mult"]
  p[3] <- p[3] + center_perturb[i, "beta_delta"]
  rescue_starts[[length(rescue_starts) + 1L]] <- make_start(
    sprintf("top_endpoint_centered_perturb%02d", i), "centered", p)
}
for (i in seq_along(endpoint_raw)) {
  rescue_starts[[length(rescue_starts) + 1L]] <- make_start(
    paste0(names(endpoint_raw)[i], "_raw"), "raw", endpoint_raw[[i]])
}
raw_perturb <- rbind(c(0.99, 0.50, -0.030), c(1.01, 2.00, 0.030))
colnames(raw_perturb) <- c("lambda_mult", "mu_mult", "beta_delta")
for (i in seq_len(nrow(raw_perturb))) {
  p <- top_raw
  p[1] <- p[1] * raw_perturb[i, "lambda_mult"]
  p[2] <- p[2] * raw_perturb[i, "mu_mult"]
  p[3] <- p[3] + raw_perturb[i, "beta_delta"]
  rescue_starts[[length(rescue_starts) + 1L]] <- make_start(
    sprintf("top_endpoint_raw_perturb%02d", i), "raw", p)
}
if (length(rescue_starts) < 5L ||
    length(unique(vapply(rescue_starts, function(x) fmt_num(x$raw_start), character(1)))) < 5L) {
  stop("Fewer than five distinct deterministic rescue starts")
}

control_for <- function(item) {
  floor_scale <- if (item$parameterization == "centered") c(0.1, 0.001, 0.1) else c(0.1, 1, 0.1)
  list(maxit = 10000L, reltol = 1e-10, parscale = pmax(abs(item$native_start), floor_scale))
}

signature_components <- data.frame(
  component = c(
    "signature_schema", "driver_sha256", "tree_sha256", "environment_sha256",
    "round1_rds_sha256", "round1_analysis_signature", "sampling_fraction",
    "conditioning", "environment_df", "total_time", "temperature_reference",
    "optimizer", "maxit", "reltol", "adaptive_recipe", "R_version",
    paste0("package_version:", c("ape", "picante", "pspline")),
    paste0("vendor_sha256:", names(vendor_sha))
  ),
  value = c(
    "03b_rescue_and_finalize_v3_signature_1", driver_sha, tree_sha, env_sha,
    round1_sha, round1$analysis_signature, sprintf("%.17g", sampling_fraction),
    conditioning, as.character(environment_df), sprintf("%.17g", total_time),
    sprintf("%.17g", temperature_ref), "Nelder-Mead", "10000", "1e-10",
    "if best not duplicated: six centered perturbations of best raw-equivalent endpoint",
    R.version.string,
    vapply(c("ape", "picante", "pspline"), function(x) as.character(packageVersion(x)), character(1)),
    unname(vendor_sha)
  ),
  stringsAsFactors = FALSE
)
for (i in seq_along(rescue_starts)) {
  item <- rescue_starts[[i]]
  ctl <- control_for(item)
  signature_components <- rbind(signature_components, data.frame(
    component = sprintf("rescue_start:%02d:%s", i, item$label),
    value = paste(item$parameterization, fmt_num(item$native_start), fmt_num(item$raw_start),
                  fmt_num(ctl$parscale), sep = "|"), stringsAsFactors = FALSE
  ))
}

signature_file <- file.path(final_model_dir, "analysis_signature_v3.tsv")
signature_tmp <- tempfile(pattern = "analysis_signature_v3_", tmpdir = final_model_dir, fileext = ".tsv")
write.table(signature_components, signature_tmp, sep = "\t", quote = FALSE,
            row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
analysis_signature <- sha256(signature_tmp)
if (file.exists(signature_file)) {
  if (!identical(readLines(signature_file, warn = FALSE), readLines(signature_tmp, warn = FALSE))) {
    stop("Existing v3 signature components differ")
  }
  unlink(signature_tmp)
} else if (!file.rename(signature_tmp, signature_file)) {
  stop("Could not atomically place v3 signature components")
}
signature_hash_file <- file.path(final_model_dir, "analysis_signature_v3.sha256")
signature_hash_text <- paste(analysis_signature, basename(signature_file))
if (file.exists(signature_hash_file)) {
  if (!identical(readLines(signature_hash_file, warn = FALSE), signature_hash_text)) {
    stop("Existing v3 signature hash differs")
  }
} else {
  writeLines(signature_hash_text, signature_hash_file, useBytes = TRUE)
}

rescue_run_signature <- function(item, start_id, control) {
  hash_text(c(
    paste0("analysis_signature=", analysis_signature),
    paste0("model=", rescue_model), paste0("start_id=", start_id),
    paste0("label=", item$label), paste0("parameterization=", item$parameterization),
    paste0("native_start=", fmt_num(item$native_start)),
    paste0("raw_start=", fmt_num(item$raw_start)),
    paste0("maxit=", control$maxit), paste0("reltol=", fmt_num(control$reltol)),
    paste0("parscale=", fmt_num(control$parscale))
  ))
}

rate_qa <- function(raw_par) {
  lambda <- rep(abs(raw_par[1]), length(temperature_grid))
  mu <- abs(raw_par[2] * exp(raw_par[3] * temperature_grid))
  finite <- all(is.finite(lambda)) && all(is.finite(mu))
  nonnegative <- all(lambda >= 0) && all(mu >= 0)
  max_rate <- if (finite) max(c(lambda, mu)) else Inf
  flags <- character(0)
  if (abs(raw_par[1]) < 1e-7) flags <- c(flags, "lambda0_near_zero")
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

rescue_checkpoint <- function(id, parameterization) {
  file.path(checkpoint_dir, sprintf("BCSTDTempVar_%s_start%02d.rds", parameterization, id))
}

fit_rescue_once <- function(item, start_id) {
  control <- control_for(item)
  run_sig <- rescue_run_signature(item, start_id, control)
  cp <- rescue_checkpoint(start_id, item$parameterization)
  if (file.exists(cp)) {
    old <- readRDS(cp)
    valid <- identical(old$analysis_signature, analysis_signature) &&
      identical(old$run_signature, run_sig) && identical(old$model, rescue_model) &&
      identical(old$parameterization, item$parameterization) &&
      identical(old$tree_sha256, tree_sha) && identical(old$environment_sha256, env_sha) &&
      identical(old$sampling_fraction, sampling_fraction) &&
      identical(old$conditioning, conditioning) && identical(old$environment_df, environment_df) &&
      identical(old$total_time, total_time) &&
      isTRUE(all.equal(old$native_start, item$native_start, tolerance = 0)) &&
      isTRUE(all.equal(old$raw_start, item$raw_start, tolerance = 0)) &&
      identical(old$control, control)
    if (!valid) stop("Rescue checkpoint signature mismatch: ", cp)
    cat(sprintf("[rescue %02d %s] resumed signature=%s\n", start_id, item$parameterization, run_sig))
    return(old)
  }

  f_lamb <- function(t, x, y) y[1]
  f_mu <- if (item$parameterization == "centered") {
    function(t, x, y) y[1] * exp(y[2] * (x - temperature_ref))
  } else {
    function(t, x, y) y[1] * exp(y[2] * x)
  }
  seed <- 2026094300L + start_id
  set.seed(seed)
  outside_warnings <- character(0)
  error_text <- NA_character_
  t0 <- proc.time()[["elapsed"]]
  fit <- tryCatch(
    withCallingHandlers(
      fit_env_bd_diagnostics_maxit_v3(
        phylo = tree, env_data = env_data, tot_time = total_time,
        f.lamb = f_lamb, f.mu = f_mu,
        lamb_par = item$native_start[1], mu_par = item$native_start[2:3],
        df = environment_df, f = sampling_fraction, meth = "Nelder-Mead",
        cst.lamb = TRUE, cst.mu = FALSE, expo.lamb = FALSE, expo.mu = FALSE,
        fix.mu = FALSE, cond = conditioning, control = control
      ),
      warning = function(w) {
        outside_warnings <<- c(outside_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      error_text <<- conditionMessage(e)
      NULL
    }
  )
  elapsed <- proc.time()[["elapsed"]] - t0

  native_final <- raw_final <- rep(NA_real_, 3L)
  names(native_final) <- names(item$native_start)
  names(raw_final) <- c("lambda0", "mu0", "beta")
  equivalence_max_abs <- NA_real_
  qa <- NULL
  if (!is.null(fit) && is.finite(fit$LH)) {
    native_final <- as.numeric(fit$optim_par)
    names(native_final) <- names(item$native_start)
    raw_final <- if (item$parameterization == "centered") {
      centered_to_raw(native_final)
    } else {
      c(lambda0 = abs(native_final[1]), mu0 = abs(native_final[2]), beta = native_final[3])
    }
    equivalence_max_abs <- if (item$parameterization == "centered") {
      rate_equivalence(native_final, raw_final)
    } else 0
    fit$optim_par_native <- fit$optim_par
    fit$lamb_par_native <- fit$lamb_par
    fit$mu_par_native <- fit$mu_par
    fit$optim_par <- unname(raw_final)
    fit$lamb_par <- unname(raw_final[1])
    fit$mu_par <- unname(raw_final[2:3])
    fit$parameterization <- item$parameterization
    fit$temperature_reference <- temperature_ref
    fit$rate_parameterization_equivalence_max_abs <- equivalence_max_abs
    qa <- rate_qa(raw_final)
  }
  result <- list(
    model = rescue_model, start_id = start_id, label = item$label,
    parameterization = item$parameterization, seed = seed, method = "Nelder-Mead",
    native_start = item$native_start, raw_start = item$raw_start,
    native_final = native_final, raw_final = raw_final,
    parameterization_equivalence_max_abs = equivalence_max_abs,
    control = control, fit = fit, rate_qa = qa,
    warnings = unique(c(outside_warnings, if (!is.null(fit)) fit$warnings else character(0))),
    error = error_text, elapsed_sec = elapsed,
    analysis_signature = analysis_signature, run_signature = run_sig,
    tree_sha256 = tree_sha, environment_sha256 = env_sha,
    driver_sha256 = driver_sha, vendor_sha256 = vendor_sha,
    total_time = total_time, sampling_fraction = sampling_fraction,
    conditioning = conditioning, environment_df = environment_df
  )
  atomic_save_rds(result, cp)
  cat(sprintf("[rescue %02d %s] logLik=%s conv=%s eval=%s elapsed=%.1fs signature=%s\n",
              start_id, item$parameterization,
              if (is.null(fit)) "NA" else format(fit$LH, digits = 14),
              if (is.null(fit)) "NA" else fit$convergence,
              if (is.null(fit)) "NA" else unname(fit$counts["function"]),
              elapsed, run_sig))
  result
}

eligible <- function(x) {
  !is.null(x$fit) && is.finite(x$fit$LH) && identical(as.integer(x$fit$convergence), 0L)
}
select_best <- function(runs) {
  ok <- which(vapply(runs, eligible, logical(1)))
  if (!length(ok)) return(NA_integer_)
  ok[which.max(vapply(runs[ok], function(x) x$fit$LH, numeric(1)))]
}
replication_count <- function(runs, best_i) {
  if (is.na(best_i)) return(0L)
  best_ll <- runs[[best_i]]$fit$LH
  starts <- vapply(runs, function(x) fmt_num(x$raw_start), character(1))
  ok <- vapply(runs, function(x) eligible(x) && abs(x$fit$LH - best_ll) <= 1e-4, logical(1))
  length(unique(starts[ok]))
}

model_cp <- file.path(checkpoint_dir, "BCSTDTempVar_rescue_complete.rds")
model_rds <- file.path(final_model_dir, "model_fit_BCSTDTempVar_rescue_v3.rds")
validate_rescue_payload <- function(x, path) {
  if (!is.list(x) || !identical(x$analysis_signature, analysis_signature) ||
      !identical(x$model, rescue_model) || is.null(x$result)) {
    stop("Rescue model payload mismatch: ", path)
  }
  invisible(TRUE)
}

if (file.exists(model_cp)) {
  rescue_payload <- readRDS(model_cp)
  validate_rescue_payload(rescue_payload, model_cp)
  if (file.exists(model_rds)) {
    validate_rescue_payload(readRDS(model_rds), model_rds)
  } else {
    atomic_save_rds(rescue_payload, model_rds)
  }
  rescue_result <- rescue_payload$result
  cat("[BCSTDTempVar rescue] resumed completed model-level checkpoint\n")
} else if (file.exists(model_rds)) {
  rescue_payload <- readRDS(model_rds)
  validate_rescue_payload(rescue_payload, model_rds)
  atomic_save_rds(rescue_payload, model_cp)
  rescue_result <- rescue_payload$result
  cat("[BCSTDTempVar rescue] recovered completion checkpoint from model RDS\n")
} else {
  requested_cores <- suppressWarnings(as.integer(Sys.getenv("RPANDA_RESCUE_CORES", "4")))
  if (!is.finite(requested_cores) || requested_cores < 1L) requested_cores <- 1L
  detected_cores <- detectCores(logical = FALSE)
  if (!is.finite(detected_cores) || detected_cores < 1L) detected_cores <- 1L
  cores <- min(requested_cores, length(rescue_starts), detected_cores)
  cat(sprintf(
    "BCSTDTempVar v3 rescue: initial_starts=%d maxit=10000 f=%.17g cond=%s df=%d Tref=%.17g cores=%d signature=%s\n",
    length(rescue_starts), sampling_fraction, conditioning, environment_df,
    temperature_ref, cores, analysis_signature
  ))
  if (.Platform$OS.type == "unix" && cores > 1L) {
    rescue_runs <- mclapply(seq_along(rescue_starts), function(i) {
      fit_rescue_once(rescue_starts[[i]], i)
    }, mc.cores = cores, mc.preschedule = FALSE)
  } else {
    rescue_runs <- lapply(seq_along(rescue_starts), function(i) fit_rescue_once(rescue_starts[[i]], i))
  }
  worker_errors <- vapply(rescue_runs, inherits, logical(1), what = "try-error")
  if (any(worker_errors)) stop("Rescue worker failure at starts: ", paste(which(worker_errors), collapse = ","))

  best_i <- select_best(rescue_runs)
  if (!is.na(best_i) && replication_count(rescue_runs, best_i) < 2L) {
    best_center <- raw_to_centered(rescue_runs[[best_i]]$raw_final)
    adaptive <- rbind(
      c(0.99, 0.80, -0.020), c(1.01, 1.20, 0.020),
      c(0.98, 1.00, -0.010), c(1.02, 1.00, 0.010),
      c(1.00, 0.60, -0.040), c(1.00, 1.60, 0.040)
    )
    colnames(adaptive) <- c("lambda_mult", "mu_mult", "beta_delta")
    adaptive_starts <- lapply(seq_len(nrow(adaptive)), function(j) {
      p <- best_center
      p[1] <- p[1] * adaptive[j, "lambda_mult"]
      p[2] <- p[2] * adaptive[j, "mu_mult"]
      p[3] <- p[3] + adaptive[j, "beta_delta"]
      make_start(sprintf("adaptive_centered_%02d", j), "centered", p)
    })
    first_id <- length(rescue_runs) + 1L
    ids <- first_id + seq_along(adaptive_starts) - 1L
    if (.Platform$OS.type == "unix" && cores > 1L) {
      extra_runs <- mclapply(seq_along(adaptive_starts), function(j) {
        fit_rescue_once(adaptive_starts[[j]], ids[j])
      }, mc.cores = min(cores, length(adaptive_starts)), mc.preschedule = FALSE)
    } else {
      extra_runs <- lapply(seq_along(adaptive_starts), function(j) {
        fit_rescue_once(adaptive_starts[[j]], ids[j])
      })
    }
    worker_errors <- vapply(extra_runs, inherits, logical(1), what = "try-error")
    if (any(worker_errors)) stop("Adaptive rescue worker failure")
    rescue_runs <- c(rescue_runs, extra_runs)
  }
  best_i <- select_best(rescue_runs)
  rescue_result <- list(
    model = rescue_model, runs = rescue_runs, best_i = best_i,
    analysis_signature = analysis_signature,
    best_start_count = replication_count(rescue_runs, best_i),
    temperature_reference = temperature_ref,
    centered_parameterization_statement = paste(
      "mu_ref=mu(Tref); mu(T)=mu_ref*exp(beta*(T-Tref));",
      "mu0=mu_ref*exp(-beta*Tref); likelihood and rate function unchanged"
    )
  )
  rescue_payload <- list(
    analysis_signature = analysis_signature, model = rescue_model,
    result = rescue_result, completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    configuration = list(
      tree_sha256 = tree_sha, environment_sha256 = env_sha,
      round1_sha256 = round1_sha, driver_sha256 = driver_sha,
      vendor_sha256 = vendor_sha, sampling_fraction = sampling_fraction,
      conditioning = conditioning, environment_df = environment_df,
      total_time = total_time, temperature_reference = temperature_ref,
      maxit = 10000L, reltol = 1e-10
    )
  )
  atomic_save_rds(rescue_payload, model_rds)
  atomic_save_rds(rescue_payload, model_cp)
  cat("[BCSTDTempVar rescue] model-level RDS and completion checkpoint saved\n")
}

rescue_best_i <- rescue_result$best_i
if (is.na(rescue_best_i)) stop("BCSTDTempVar rescue produced no finite converged fit")
rescue_best <- rescue_result$runs[[rescue_best_i]]

run_row_round1 <- function(x, selected) {
  fit <- x$fit
  qa <- x$rate_qa
  data.frame(
    model = x$model, start_id = x$start_id, run_phase = "first_round_v2",
    start_label = paste0("first_round_start", sprintf("%02d", x$start_id)),
    parameterization = "raw_original", seed = x$seed, method = x$method,
    maxit = 500L, reltol = sqrt(.Machine$double.eps),
    parscale = paste(rep(1, length(x$start)), collapse = ";"),
    starting_values_native = fmt_num(x$start), starting_values_raw = fmt_num(x$start),
    final_values_native = if (is.null(fit)) NA_character_ else fmt_num(fit$optim_par),
    final_values_raw = if (is.null(fit)) NA_character_ else fmt_num(fit$optim_par),
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
    rate_explosion = if (is.null(qa) || is.null(qa$explosive)) NA else qa$explosive,
    min_lambda = if (is.null(qa)) NA_real_ else qa$min_lambda,
    max_lambda = if (is.null(qa)) NA_real_ else qa$max_lambda,
    min_mu = if (is.null(qa)) NA_real_ else qa$min_mu,
    max_mu = if (is.null(qa)) NA_real_ else qa$max_mu,
    boundary_flag = if (is.null(qa)) "fit_failed" else qa$boundary_flag,
    parameterization_equivalence_max_abs = 0, selected = selected,
    total_time = total_time, sampling_fraction = sampling_fraction,
    conditioning = conditioning, environment_df = environment_df,
    analysis_signature = analysis_signature,
    source_analysis_signature = x$analysis_signature,
    run_signature = x$run_signature, tree_sha256 = x$tree_sha256,
    environment_sha256 = x$env_sha256, source_driver_sha256 = x$driver_sha256,
    stringsAsFactors = FALSE
  )
}

run_row_rescue <- function(x, selected) {
  fit <- x$fit
  qa <- x$rate_qa
  data.frame(
    model = x$model, start_id = x$start_id, run_phase = "rescue_v3",
    start_label = x$label, parameterization = x$parameterization,
    seed = x$seed, method = x$method, maxit = x$control$maxit,
    reltol = x$control$reltol, parscale = fmt_num(x$control$parscale),
    starting_values_native = fmt_num(x$native_start), starting_values_raw = fmt_num(x$raw_start),
    final_values_native = if (is.null(fit)) NA_character_ else fmt_num(x$native_final),
    final_values_raw = if (is.null(fit)) NA_character_ else fmt_num(x$raw_final),
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
    parameterization_equivalence_max_abs = x$parameterization_equivalence_max_abs,
    selected = selected, total_time = x$total_time,
    sampling_fraction = x$sampling_fraction, conditioning = x$conditioning,
    environment_df = x$environment_df, analysis_signature = analysis_signature,
    source_analysis_signature = x$analysis_signature,
    run_signature = x$run_signature, tree_sha256 = x$tree_sha256,
    environment_sha256 = x$environment_sha256, source_driver_sha256 = x$driver_sha256,
    stringsAsFactors = FALSE
  )
}

run_rows <- list()
for (m in carry_models) {
  z <- round1$all_results[[m]]
  for (i in seq_along(z$runs)) {
    run_rows[[length(run_rows) + 1L]] <- run_row_round1(z$runs[[i]], i == z$best_i)
  }
}
for (i in seq_along(rescue_result$runs)) {
  run_rows[[length(run_rows) + 1L]] <- run_row_rescue(rescue_result$runs[[i]], i == rescue_best_i)
}
optimization_runs <- clean_frame(do.call(rbind, run_rows))

old_seed_rows <- lapply(endpoint_order, function(i) {
  x <- old_rescue$runs[[i]]
  data.frame(
    source_start_id = x$start_id, starting_values = fmt_num(x$start),
    endpoint_values = fmt_num(x$fit$optim_par), logLik = x$fit$LH,
    convergence = x$fit$convergence,
    warning = if (length(x$warnings)) paste(x$warnings, collapse = " | ") else NA_character_,
    analysis_signature = round1$analysis_signature, stringsAsFactors = FALSE
  )
})
old_seed_table <- clean_frame(do.call(rbind, old_seed_rows))

rescue_replication <- replication_count(rescue_result$runs, rescue_best_i)
rescue_raw <- rescue_best$raw_final
rescue_qa <- rescue_best$rate_qa
rescue_row <- data.frame(
  model = rescue_model, family = "temperature",
  lambda_formula = "lambda(t) = lambda0",
  mu_formula = "mu(t) = mu0 * exp(beta * T(t))",
  k = 3L, logLik = rescue_best$fit$LH,
  AIC = -2 * rescue_best$fit$LH + 6,
  AICc = -2 * rescue_best$fit$LH + 6 + (2 * 3 * 4) / (n_obs - 3 - 1),
  lambda0 = abs(rescue_raw[1]), alpha = NA_real_,
  mu0 = abs(rescue_raw[2]), beta = rescue_raw[3],
  convergence = rescue_best$fit$convergence,
  boundary_flag = rescue_qa$boundary_flag,
  best_start_count = rescue_replication,
  notes = paste(
    "maximum finite converged logLik across long-run rescue starts;",
    "centered parameterization is likelihood-equivalent and estimates were converted to original mu0/beta form"
  ),
  stringsAsFactors = FALSE
)

model_table <- rbind(carry_table, rescue_row)
model_table$AIC <- -2 * model_table$logLik + 2 * model_table$k
model_table$AICc <- model_table$AIC +
  (2 * model_table$k * (model_table$k + 1)) / (n_obs - model_table$k - 1)
model_table$deltaAICc <- model_table$AICc - min(model_table$AICc)
rel <- exp(-0.5 * model_table$deltaAICc)
model_table$Akaike_weight <- rel / sum(rel)
model_table <- model_table[order(model_table$AICc, model_table$model), ]
rownames(model_table) <- NULL
front <- c(
  "model", "lambda_formula", "mu_formula", "k", "logLik", "AIC", "AICc",
  "deltaAICc", "Akaike_weight", "lambda0", "alpha", "mu0", "beta",
  "convergence", "boundary_flag", "best_start_count", "notes", "family"
)
model_table <- clean_frame(model_table[front])
formatted <- model_table
for (nm in c("logLik", "AIC", "AICc", "deltaAICc")) formatted[[nm]] <- sprintf("%.3f", formatted[[nm]])
formatted$Akaike_weight <- sprintf("%.4f", formatted$Akaike_weight)
for (nm in c("lambda0", "alpha", "mu0", "beta")) {
  formatted[[nm]] <- ifelse(is.na(formatted[[nm]]), "—", sprintf("%.6g", formatted[[nm]]))
}
formatted <- clean_frame(formatted)

dominance_rows <- lapply(model_names, function(m) {
  selected_ll <- model_table$logLik[match(m, model_table$model)]
  formal_nonconv <- optimization_runs$logLik[
    optimization_runs$model == m & is.finite(optimization_runs$logLik) &
      !is.na(optimization_runs$convergence) & optimization_runs$convergence != 0L
  ]
  source_nonconv <- if (m == rescue_model) endpoint_ll else numeric(0)
  vals <- c(formal_nonconv, source_nonconv)
  max_nonconv <- if (length(vals)) max(vals) else NA_real_
  data.frame(
    model = m, selected_converged_logLik = selected_ll,
    maximum_finite_nonconverged_logLik = max_nonconv,
    n_finite_nonconverged_considered = length(vals),
    selected_not_lower_than_nonconverged = if (length(vals)) selected_ll + 1e-8 >= max_nonconv else TRUE,
    stringsAsFactors = FALSE
  )
})
dominance_qa <- do.call(rbind, dominance_rows)

selected_rows <- optimization_runs$selected
run_k <- vapply(optimization_runs$model, function(m) model_table$k[match(m, model_table$model)], integer(1))
finite_aicc <- is.finite(optimization_runs$logLik) & is.finite(optimization_runs$AICc_from_vendor)
recomputed_run_aicc <- -2 * optimization_runs$logLik + 2 * run_k +
  (2 * run_k * (run_k + 1)) / (n_obs - run_k - 1)
vendor_aicc_agrees <- all(abs(optimization_runs$AICc_from_vendor[finite_aicc] -
                               recomputed_run_aicc[finite_aicc]) < 1e-10)
centered_equiv <- optimization_runs$parameterization_equivalence_max_abs[
  optimization_runs$run_phase == "rescue_v3" & optimization_runs$parameterization == "centered" &
    is.finite(optimization_runs$parameterization_equivalence_max_abs)
]
max_centered_equiv <- if (length(centered_equiv)) max(centered_equiv) else Inf

same_f <- all(optimization_runs$sampling_fraction == sampling_fraction)
same_cond <- all(optimization_runs$conditioning == conditioning)
same_df <- all(optimization_runs$environment_df == environment_df)
same_time <- all(optimization_runs$total_time == total_time)
same_sig <- all(optimization_runs$analysis_signature == analysis_signature)
same_tree <- all(optimization_runs$tree_sha256 == tree_sha)
same_env <- all(optimization_runs$environment_sha256 == env_sha)
no_higher_nonconv <- all(dominance_qa$selected_not_lower_than_nonconverged)
no_control_chars <- !has_control_text(optimization_runs) && !has_control_text(model_table) &&
  !has_control_text(formatted) && !has_control_text(old_seed_table) && !has_control_text(dominance_qa)
preflight_outputs <- all(
  preflight_tsv(optimization_runs), preflight_tsv(model_table), preflight_tsv(formatted),
  preflight_tsv(old_seed_table), preflight_tsv(dominance_qa)
)

qa <- data.frame(
  check = c(
    "prespecified_models_fitted", "all_selected_convergence_zero",
    "all_selected_rates_finite", "all_selected_rates_nonnegative",
    "no_selected_rate_explosion_gt_100_per_Ma",
    "BCSTDTempVar_best_reproduced_twice_1e-4_distinct_raw_starts",
    "no_higher_finite_nonconverged_result", "centered_to_raw_rate_equivalence_lt_1e-12",
    "rescue_Nelder_Mead_maxit_at_least_5000", "centered_and_raw_rescue_crosscheck_present",
    "vendor_AICc_matches_full_precision_recalculation", "akaike_weights_sum",
    "same_sampling_fraction_all_runs", "same_conditioning_all_runs",
    "same_temperature_df_all_runs", "same_total_time_all_runs",
    "same_final_analysis_signature_all_runs", "same_tree_sha256_all_runs",
    "same_environment_sha256_all_runs", "TSV_strings_have_no_CR_LF_TAB",
    "TSV_roundtrip_preflight", "unrounded_logLik_used_for_AICc"
  ),
  observed = c(
    nrow(model_table), all(model_table$convergence == 0L),
    all(optimization_runs$rates_finite[selected_rows]),
    all(optimization_runs$rates_nonnegative[selected_rows]),
    all(!optimization_runs$rate_explosion[selected_rows]), rescue_replication,
    no_higher_nonconv, format(max_centered_equiv, digits = 17),
    min(optimization_runs$maxit[optimization_runs$run_phase == "rescue_v3"]),
    all(c("centered", "raw") %in% optimization_runs$parameterization[optimization_runs$run_phase == "rescue_v3"]),
    vendor_aicc_agrees, format(sum(model_table$Akaike_weight), digits = 17),
    same_f, same_cond, same_df, same_time, same_sig, same_tree, same_env,
    no_control_chars, preflight_outputs, TRUE
  ),
  expected = c(
    "10", "TRUE", "TRUE", "TRUE", "TRUE", ">=2", "TRUE", "<1e-12",
    ">=5000", "TRUE", "TRUE", "1 within floating-point error",
    "TRUE", "TRUE", "TRUE", "TRUE", "TRUE", "TRUE", "TRUE", "TRUE", "TRUE", "TRUE"
  ),
  pass = c(
    nrow(model_table) == 10L, all(model_table$convergence == 0L),
    all(optimization_runs$rates_finite[selected_rows]),
    all(optimization_runs$rates_nonnegative[selected_rows]),
    all(!optimization_runs$rate_explosion[selected_rows]), rescue_replication >= 2L,
    no_higher_nonconv, max_centered_equiv < 1e-12,
    all(optimization_runs$maxit[optimization_runs$run_phase == "rescue_v3"] >= 5000L),
    all(c("centered", "raw") %in% optimization_runs$parameterization[optimization_runs$run_phase == "rescue_v3"]),
    vendor_aicc_agrees, abs(sum(model_table$Akaike_weight) - 1) < 1e-12,
    same_f, same_cond, same_df, same_time, same_sig, same_tree, same_env,
    no_control_chars, preflight_outputs, TRUE
  ),
  stringsAsFactors = FALSE
)
qa <- clean_frame(qa)
if (!preflight_tsv(qa)) stop("QA table itself failed TSV roundtrip preflight")

write_clean_tsv(optimization_runs, file.path(final_model_dir, "optimization_runs.tsv"))
write_clean_tsv(old_seed_table, file.path(final_model_dir, "BCSTDTempVar_round1_nonconverged_seed_endpoints.tsv"))
write_clean_tsv(model_table, file.path(final_table_dir, "model_table_10models_full_precision.tsv"))
write_clean_tsv(formatted, file.path(final_table_dir, "model_table_10models_formatted.tsv"))
write_clean_tsv(dominance_qa, file.path(final_table_dir, "nonconverged_dominance_QA.tsv"))
write_clean_tsv(qa, file.path(final_table_dir, "model_selection_QA.tsv"))

selected_fits <- round1$selected_fits[carry_models]
selected_fits[[rescue_model]] <- rescue_best$fit
selected_fits <- selected_fits[model_names]
all_results <- round1$all_results[carry_models]
all_results[[rescue_model]] <- rescue_result
all_results <- all_results[model_names]
final_payload <- list(
  analysis_signature = analysis_signature,
  signature_components = signature_components,
  source_round1 = list(
    path = round1_file, sha256 = round1_sha,
    analysis_signature = round1$analysis_signature,
    carried_models = carry_models,
    excluded_model = rescue_model,
    exclusion_reason = "first-round higher finite nonconverged endpoints and best_start_count=1"
  ),
  centered_parameterization = list(
    temperature_reference = temperature_ref,
    statement = rescue_result$centered_parameterization_statement,
    maximum_rate_difference = max_centered_equiv
  ),
  selected_fits = selected_fits, all_results = all_results,
  model_table = model_table, optimization_runs = optimization_runs,
  model_selection_QA = qa, nonconverged_dominance_QA = dominance_qa,
  tree = tree, env = env_data, total_time = total_time,
  sampling_fraction = sampling_fraction, df = environment_df, cond = conditioning,
  input_sha256 = c(tree = tree_sha, environment = env_sha, round1_rds = round1_sha),
  code_sha256 = c(driver = driver_sha, vendor_sha)
)
atomic_save_rds(final_payload, file.path(final_model_dir, "model_fits_10models.rds"))

if (!all(qa$pass)) {
  failed <- qa$check[!qa$pass]
  stop("v3 final QA failed: ", paste(failed, collapse = ", "))
}
cat("v3 rescue/finalization complete and QA-passing. Best model: ", model_table$model[1],
    "; rescued BCSTDTempVar logLik: ", format(rescue_best$fit$LH, digits = 14),
    "; reproduced starts: ", rescue_replication,
    "; signature: ", analysis_signature, "\n", sep = "")
