#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: 10_compile_config_lint.R <run_root>", call. = FALSE)
}
run_root <- normalizePath(args[[1]], mustWork = TRUE)
config_dir <- file.path(run_root, "configs")
qa_dir <- file.path(run_root, "qa")
out_path <- file.path(qa_dir, "config_lint.txt")
if (file.exists(out_path)) {
  stop("Refusing to overwrite existing output: ", out_path, call. = FALSE)
}

config_names <- c(
  "prime.cfg",
  sprintf("cv_run_%02d.cfg", 1:3),
  sprintf("cv_ext_run_%02d.cfg", 1:3),
  sprintf("cv_full_run_%02d.cfg", 1:3),
  sprintf("cv_finalgrid_run_%02d.cfg", 1:3),
  "final.cfg",
  "final_baseline_no_thorough.cfg",
  "final_repeat_same_seed.cfg",
  "final_sensitivity_1e-12.cfg",
  "final_sensitivity_1e-13.cfg"
)

rows <- lapply(config_names, function(config_name) {
  config_path <- file.path(config_dir, config_name)
  if (!file.exists(config_path)) stop("Missing config: ", config_path)
  stem <- sub("\\.cfg$", "", config_name)
  lint_path <- file.path(qa_dir, paste0("config_lint_", stem, ".tsv"))
  if (!file.exists(lint_path)) stop("Missing lint table: ", lint_path)

  lint <- read.delim(lint_path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("check", "pass") %in% names(lint))) {
    stop("Malformed lint table: ", lint_path)
  }
  lines <- readLines(config_path, warn = FALSE)
  active <- trimws(sub("#.*$", "", lines))
  active <- active[nzchar(active)]
  calibration_lines <- active[
    grepl("^(mrca|min|max)[[:space:]]*=", active)
  ]
  numsites_lines <- active[grepl("^numsites[[:space:]]*=", active)]
  treefile_lines <- active[grepl("^treefile[[:space:]]*=", active)]
  outfile_lines <- active[grepl("^outfile[[:space:]]*=", active)]
  outfile_ok <- if (identical(config_name, "prime.cfg")) {
    length(outfile_lines) == 0L
  } else {
    length(outfile_lines) == 1L &&
      grepl(run_root, outfile_lines, fixed = TRUE)
  }
  forbidden_numeric <- c("101.4", "36.2", "40.7", "36.210861", "40.662537")

  checks <- c(
    lint_table_all_pass = all(lint$pass),
    exactly_two_mrca = sum(grepl("^mrca[[:space:]]*=", active)) == 2L,
    exactly_two_min = sum(grepl("^min[[:space:]]*=", active)) == 2L,
    exactly_two_max = sum(grepl("^max[[:space:]]*=", active)) == 2L,
    no_hesperiinae_calibration =
      !any(grepl("HESPERIINAE", calibration_lines, ignore.case = TRUE)),
    forbidden_old_ages_absent =
      !any(vapply(forbidden_numeric, function(x) {
        any(grepl(x, calibration_lines, fixed = TRUE))
      }, logical(1))),
    numsites_exact =
      length(numsites_lines) == 1L &&
      identical(gsub("[[:space:]]", "", numsites_lines), "numsites=182682"),
    treefile_exact =
      length(treefile_lines) == 1L &&
      identical(
        trimws(sub("^[^=]+=", "", treefile_lines)),
        file.path(run_root, "input", "verified_495_tip_input.treefile")
      ),
    outfile_policy = outfile_ok,
    no_log_pen = !any(grepl("^log_pen([[:space:]]*=|$)", active)),
    no_rescale_or_reanchor =
      !any(grepl("rescal|re-?anchor", active, ignore.case = TRUE))
  )

  data.frame(
    config = config_name,
    lint_table = basename(lint_path),
    lint_pass = sum(lint$pass),
    lint_total = nrow(lint),
    status = if (all(checks)) "PASS" else "FAIL",
    failed_consolidated_checks =
      if (all(checks)) "" else paste(names(checks)[!checks], collapse = ","),
    stringsAsFactors = FALSE
  )
})

tab <- do.call(rbind, rows)
lines <- c(
  "treePL configuration lint — consolidated final report",
  paste0("run_root\t", run_root),
  paste0("generated_at\t", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  paste0("configs_checked\t", nrow(tab)),
  paste0("overall_status\t", if (all(tab$status == "PASS")) "PASS" else "FAIL"),
  "global_requirements\texactly 2 calibrations; no Hesperiinae calibration; no old ages 101.4/36.2/40.7; numsites=182682; run-local outputs; no rescale/re-anchor",
  "",
  paste(names(tab), collapse = "\t"),
  apply(tab, 1, paste, collapse = "\t")
)
writeLines(lines, out_path, useBytes = TRUE)
cat("Consolidated config lint:", if (all(tab$status == "PASS")) "PASS" else "FAIL",
    "-", nrow(tab), "configs\n")
if (!all(tab$status == "PASS")) quit(status = 2L)
