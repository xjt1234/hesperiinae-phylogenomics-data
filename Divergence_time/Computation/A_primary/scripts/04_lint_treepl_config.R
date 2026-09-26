#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("Usage: 04_lint_treepl_config.R <config> <prime|cv|cv_extension|cv_full|cv_finalgrid|final_baseline|final> <run_root> <output.tsv>")
}

cfg_path <- normalizePath(args[[1]], mustWork = TRUE)
phase <- args[[2]]
run_root <- normalizePath(args[[3]], mustWork = TRUE)
out_path <- args[[4]]
if (!phase %in% c("prime", "cv", "cv_extension", "cv_full", "cv_finalgrid", "final_baseline", "final")) stop("Unknown phase: ", phase)

raw <- readLines(cfg_path, warn = FALSE)
active <- trimws(sub("#.*$", "", raw))
active <- active[nzchar(active)]
key <- trimws(sub("=.*$", "", active))
key[!grepl("=", active, fixed = TRUE)] <- active[!grepl("=", active, fixed = TRUE)]
value <- ifelse(grepl("=", active, fixed = TRUE), trimws(sub("^[^=]*=", "", active)), "")

values_for <- function(k) value[key == k]
has_key <- function(k) any(key == k)
exact_line <- function(x) sum(active == x) == 1L
within_run <- function(p) {
  if (!length(p) || !nzchar(p[[1]])) return(FALSE)
  resolved <- normalizePath(p[[1]], mustWork = FALSE)
  startsWith(paste0(resolved, "/"), paste0(run_root, "/"))
}

checks <- list()
add <- function(name, pass, observed, expected) {
  checks[[length(checks) + 1L]] <<- data.frame(
    check = name, pass = isTRUE(pass), observed = as.character(observed),
    expected = as.character(expected), stringsAsFactors = FALSE
  )
}

add("treefile_exact", identical(values_for("treefile"), paste0(run_root, "/input/verified_495_tip_input.treefile")),
    paste(values_for("treefile"), collapse = ";"), paste0(run_root, "/input/verified_495_tip_input.treefile"))
add("numsites_182682", identical(values_for("numsites"), "182682"), paste(values_for("numsites"), collapse = ";"), "182682")
add("exactly_two_mrca", sum(key == "mrca") == 2L, sum(key == "mrca"), 2)
add("exactly_two_min", sum(key == "min") == 2L, sum(key == "min"), 2)
add("exactly_two_max", sum(key == "max") == 2L, sum(key == "max"), 2)
add("root_mrca_exact", exact_line("mrca = CROWN_PAPILIONOIDEA Parnassius_apollo_ncbi Aeromachus_catocyanea_mydata"),
    sum(active == "mrca = CROWN_PAPILIONOIDEA Parnassius_apollo_ncbi Aeromachus_catocyanea_mydata"), 1)
add("root_bounds_exact",
    exact_line("min = CROWN_PAPILIONOIDEA 91.5046") && exact_line("max = CROWN_PAPILIONOIDEA 100.8925"),
    paste(values_for("min")[grepl("CROWN_PAPILIONOIDEA", values_for("min"))],
          values_for("max")[grepl("CROWN_PAPILIONOIDEA", values_for("max"))], collapse = ";"),
    "91.5046-100.8925")
add("papilionidae_mrca_exact",
    exact_line("mrca = CROWN_PAPILIONIDAE_NONBARONIA Parnassius_apollo_ncbi Graphium_cloanthus_mydata"),
    sum(active == "mrca = CROWN_PAPILIONIDAE_NONBARONIA Parnassius_apollo_ncbi Graphium_cloanthus_mydata"), 1)
add("papilionidae_bounds_exact",
    exact_line("min = CROWN_PAPILIONIDAE_NONBARONIA 44.1968") && exact_line("max = CROWN_PAPILIONIDAE_NONBARONIA 52.9473"),
    paste(values_for("min")[grepl("PAPILIONIDAE", values_for("min"))],
          values_for("max")[grepl("PAPILIONIDAE", values_for("max"))], collapse = ";"),
    "44.1968-52.9473")
calibration_lines <- active[key %in% c("mrca", "min", "max")]
add("no_hesperiinae_calibration", !any(grepl("HESPERIINAE", calibration_lines, ignore.case = TRUE)),
    paste(calibration_lines[grepl("HESPERIINAE", calibration_lines, ignore.case = TRUE)], collapse = ";"), "none")
for (bad in c("101.4", "36.2", "40.7", "36.210861", "40.662537")) {
  add(paste0("forbidden_value_absent_", gsub("\\.", "_", bad)),
      !any(grepl(bad, active, fixed = TRUE)), sum(grepl(bad, active, fixed = TRUE)), 0)
}
add("log_pen_absent", !has_key("log_pen"), has_key("log_pen"), FALSE)
seed_value <- suppressWarnings(as.numeric(values_for("seed")))
add("explicit_positive_seed", length(seed_value) == 1L && is.finite(seed_value) && seed_value > 0,
    paste(values_for("seed"), collapse = ";"), "one positive integer")

if (phase == "prime") {
  add("prime_enabled", has_key("prime"), has_key("prime"), TRUE)
  add("cv_disabled", !has_key("cv") && !has_key("randomcv"), paste(key[key %in% c("cv", "randomcv")], collapse = ";"), "none")
  add("thorough_disabled_for_prime", !has_key("thorough"), has_key("thorough"), FALSE)
  add("no_outfile_for_prime", !has_key("outfile") && !has_key("cvoutfile"), paste(key[key %in% c("outfile", "cvoutfile")], collapse = ";"), "none")
} else if (phase %in% c("cv", "cv_extension", "cv_full", "cv_finalgrid")) {
  add("randomcv_enabled", has_key("randomcv"), has_key("randomcv"), TRUE)
  add("prime_disabled", !has_key("prime"), has_key("prime"), FALSE)
  add("thorough_disabled_for_cv", !has_key("thorough"), has_key("thorough"), FALSE)
  expected_start <- if (phase == "cv_extension") "0.00001" else "100000"
  expected_stop <- if (phase == "cv") "0.00001" else if (phase %in% c("cv_extension", "cv_full")) "0.000000001" else "0.0000000000001"
  add("cv_grid_start", identical(values_for("cvstart"), expected_start), paste(values_for("cvstart"), collapse = ";"), expected_start)
  add("cv_grid_stop", identical(values_for("cvstop"), expected_stop), paste(values_for("cvstop"), collapse = ";"), expected_stop)
  add("cv_grid_multiplier", identical(values_for("cvmultstep"), "0.1"), paste(values_for("cvmultstep"), collapse = ";"), "0.1")
  add("cviter_3", identical(values_for("cviter"), "3"), paste(values_for("cviter"), collapse = ";"), "3")
  add("nthreads_10", identical(values_for("nthreads"), "10"), paste(values_for("nthreads"), collapse = ";"), "10")
  add("cvoutfile_within_run", length(values_for("cvoutfile")) == 1L && within_run(values_for("cvoutfile")),
      paste(values_for("cvoutfile"), collapse = ";"), paste0(run_root, "/..."))
  add("outfile_within_run", length(values_for("outfile")) == 1L && within_run(values_for("outfile")),
      paste(values_for("outfile"), collapse = ";"), paste0(run_root, "/..."))
} else {
  cv_keys <- c("prime", "cv", "randomcv", "cvstart", "cvstop", "cvmultstep", "cvoutfile")
  if (phase == "final") {
    add("thorough_enabled", has_key("thorough"), has_key("thorough"), TRUE)
  } else {
    add("thorough_disabled_for_baseline", !has_key("thorough"), has_key("thorough"), FALSE)
  }
  add("prime_and_cv_disabled", !any(key %in% cv_keys), paste(key[key %in% cv_keys], collapse = ";"), "none")
  smooth_value <- suppressWarnings(as.numeric(values_for("smooth")))
  add("one_positive_smoothing", length(smooth_value) == 1L && is.finite(smooth_value) && smooth_value > 0,
      paste(values_for("smooth"), collapse = ";"), "one positive value")
  add("outfile_within_run", length(values_for("outfile")) == 1L && within_run(values_for("outfile")),
      paste(values_for("outfile"), collapse = ";"), paste0(run_root, "/..."))
}

result <- do.call(rbind, checks)
write.table(result, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Config lint %s: %d/%d PASS\n", basename(cfg_path), sum(result$pass), nrow(result)))
if (!all(result$pass)) stop("Config lint failed; see ", out_path)
