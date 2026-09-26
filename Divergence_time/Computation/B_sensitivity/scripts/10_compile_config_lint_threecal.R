#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: 10_compile_config_lint_threecal.R <run_root>", call. = FALSE)
}

fail <- function(...) stop(paste0(...), call. = FALSE)

run_root <- normalizePath(args[[1]], mustWork = TRUE)
config_dir <- file.path(run_root, "configs")
log_dir <- file.path(run_root, "logs")
qa_dir <- file.path(run_root, "qa")
required_dirs <- c(config_dir, log_dir, qa_dir)
missing_dirs <- required_dirs[!dir.exists(required_dirs)]
if (length(missing_dirs)) {
  fail("Missing required directory/directories: ", paste(missing_dirs, collapse = ", "))
}

# This manifest gate is intentionally before any report construction. Final
# identities are data-dependent, so a pre-final invocation must never emit PASS.
manifest_path <- file.path(qa_dir, "final_variant_manifest.tsv")
if (!file.exists(manifest_path) || isTRUE(file.info(manifest_path)$isdir)) {
  fail(
    "Missing qa/final_variant_manifest.tsv. Run this consolidation only after ",
    "CV selection, every final run, and final-manifest creation are complete."
  )
}

out_path <- file.path(qa_dir, "config_lint_consolidated.txt")
if (file.exists(out_path)) {
  fail("Refusing to overwrite existing output: ", out_path)
}

inside_run <- function(path) {
  startsWith(paste0(path, "/"), paste0(run_root, "/"))
}

resolve_run_file <- function(path, label, must_exist = TRUE) {
  if (length(path) != 1L || is.na(path) || !nzchar(trimws(path))) {
    fail(label, " is empty")
  }
  candidate <- if (startsWith(path, "/")) path else file.path(run_root, path)
  resolved <- normalizePath(candidate, mustWork = must_exist)
  if (!inside_run(resolved)) fail(label, " resolves outside run root: ", resolved)
  if (must_exist && isTRUE(file.info(resolved)$isdir)) {
    fail(label, " must resolve to a file, not a directory: ", resolved)
  }
  resolved
}

relative_run_path <- function(path) {
  resolved <- normalizePath(path, mustWork = TRUE)
  if (!inside_run(resolved)) fail("Path is outside run root: ", resolved)
  substring(resolved, nchar(run_root) + 2L)
}

expected_treepl_binary <- resolve_run_file(
  file.path(run_root, "software", "treePL_randomcv_fix_src", "treePL"),
  "required run-local patched treePL binary"
)

sha256_file <- function(path) {
  resolved <- normalizePath(path, mustWork = TRUE)
  result <- system2(
    "/usr/bin/sha256sum",
    args = c("--", shQuote(resolved)),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(result, "status")
  if (is.null(status)) status <- 0L
  if (status != 0L || length(result) != 1L) {
    fail("sha256sum failed for ", resolved, ": ", paste(result, collapse = " | "))
  }
  hash <- sub("[[:space:]].*$", "", result[[1]])
  if (!grepl("^[0-9a-f]{64}$", hash)) fail("Malformed SHA-256 for ", resolved)
  hash
}

read_one_row_decision <- function(path, label) {
  if (!file.exists(path) || isTRUE(file.info(path)$isdir) ||
      file.info(path)$size <= 0) {
    fail("Missing or empty ", label, ": ", path)
  }
  tab <- read.delim(
    path,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    colClasses = "character",
    check.names = FALSE,
    na.strings = character(),
    strip.white = FALSE,
    stringsAsFactors = FALSE
  )
  if (nrow(tab) != 1L || anyDuplicated(names(tab))) {
    fail(label, " must contain exactly one row and unique column names: ", path)
  }
  tab
}

require_decision_columns <- function(tab, columns, label) {
  missing <- setdiff(columns, names(tab))
  if (length(missing)) {
    fail(label, " is missing required column(s): ", paste(missing, collapse = ", "))
  }
}

decision_value <- function(tab, key, label) {
  value <- tab[[key]][[1L]]
  if (is.na(value) || !identical(value, trimws(value))) {
    fail(label, " has NA or surrounding whitespace in field ", key)
  }
  value
}

decision_bool <- function(tab, key, label) {
  value <- decision_value(tab, key, label)
  if (!value %in% c("TRUE", "FALSE")) {
    fail(label, " field ", key, " must be exact TRUE or FALSE")
  }
  identical(value, "TRUE")
}

decision_uint <- function(tab, key, label) {
  value <- decision_value(tab, key, label)
  if (!grepl("^(0|[1-9][0-9]*)$", value)) {
    fail(label, " field ", key, " must be an unsigned integer")
  }
  numeric <- suppressWarnings(as.numeric(value))
  if (!is.finite(numeric) || numeric > 2^53) {
    fail(label, " field ", key, " is outside the exact numeric range")
  }
  numeric
}

source_decision_path <- file.path(qa_dir, "cv_decision.tsv")
source_decision <- read_one_row_decision(
  source_decision_path, "source 19-point CV decision"
)
source_required <- c(
  "status", "decision_basis", "selected_smoothing",
  "aggregate_winner_at_grid_boundary", "formal_runs",
  "formal_grid_points_per_run", "formal_rows", "formal_seed_values"
)
require_decision_columns(
  source_decision, source_required, "source 19-point CV decision"
)
if (!identical(
      decision_value(source_decision, "status", "source 19-point CV decision"),
      "REQUIRES_GRID_EXTENSION"
    ) ||
    !identical(
      decision_value(source_decision, "decision_basis", "source 19-point CV decision"),
      "lowest_median_raw_chisq_across_three_19_point_randomCV_runs"
    ) ||
    !decision_bool(
      source_decision, "aggregate_winner_at_grid_boundary",
      "source 19-point CV decision"
    ) ||
    !decision_value(
      source_decision, "selected_smoothing", "source 19-point CV decision"
    ) %in% c("", "NA") ||
    !identical(
      decision_value(source_decision, "formal_runs", "source 19-point CV decision"),
      "cv_finalgrid_run_01,cv_finalgrid_run_02,cv_finalgrid_run_03"
    ) ||
    !identical(
      decision_uint(
        source_decision, "formal_grid_points_per_run", "source 19-point CV decision"
      ),
      19
    ) ||
    !identical(
      decision_uint(source_decision, "formal_rows", "source 19-point CV decision"),
      57
    ) ||
    !identical(
      decision_value(source_decision, "formal_seed_values", "source 19-point CV decision"),
      "2026090441,2026090442,2026090443"
    )) {
  fail("Source qa/cv_decision.tsv does not document the required 19-point boundary extension")
}

extended_decision_path <- file.path(qa_dir, "cv_23point_decision.tsv")
extended_decision <- read_one_row_decision(
  extended_decision_path, "formal 23-point CV decision"
)
extended_required <- c(
  "status", "decision_basis", "selected_smoothing", "selected_median_score",
  "aggregate_winner_at_grid_boundary", "boundary_requires_extension",
  "aggregate_winner_count", "finite_precision_plateau", "stability",
  "sensitivity_required", "sensitivity_smoothing_values",
  "formal_runs", "formal_grid_points_per_run", "formal_rows",
  "formal_seed_values", "formal_exit_codes_all_zero",
  "formal_stderr_all_empty", "formal_dated_trees_all_nonempty",
  "formal_dated_trees_valid_branch_length_newick",
  "formal_dated_tree_tip_sets_match_input",
  "formal_dated_tree_topologies_match_input",
  "formal_dated_trees_ultrametric_within_tolerance",
  "formal_config_contracts_validated", "formal_config_lints_all_pass",
  "formal_metadata_validated", "formal_stdout_score_sets_validated",
  "formal_resource_logs_validated", "source_19point_decision_validated",
  "diagnostic_19point_full_stage_validation_passed",
  "diagnostic_19point_runs_used_for_selection",
  "diagnostic_19point_all_prefix_points_present",
  "diagnostic_19point_lower_boundary_in_aggregate_winner_set",
  "diagnostic_19point_aggregate_winner_values"
)
require_decision_columns(
  extended_decision, extended_required, "formal 23-point CV decision"
)
diagnostic_19point_winner_text <- decision_value(
  extended_decision, "diagnostic_19point_aggregate_winner_values",
  "formal 23-point CV decision"
)
diagnostic_19point_winner_tokens <- strsplit(
  diagnostic_19point_winner_text, ",", fixed = TRUE
)[[1L]]
if (!decision_bool(
      extended_decision,
      "diagnostic_19point_lower_boundary_in_aggregate_winner_set",
      "formal 23-point CV decision"
    ) ||
    !"0.0000000000001" %in% diagnostic_19point_winner_tokens) {
  fail(
    "Formal 23-point decision does not prove the required downward extension: ",
    "the 19-point lower boundary must be an aggregate winner and ",
    "diagnostic_19point_aggregate_winner_values must contain 0.0000000000001"
  )
}
extended_status <- decision_value(
  extended_decision, "status", "formal 23-point CV decision"
)
if (!extended_status %in% c("ACCEPTED", "ACCEPTED_WITH_SENSITIVITIES")) {
  fail("Formal 23-point CV decision is not qualified for final runs: status=", extended_status)
}
extended_boundary <- decision_bool(
  extended_decision, "aggregate_winner_at_grid_boundary",
  "formal 23-point CV decision"
)
extended_boundary_requires_extension <- decision_bool(
  extended_decision, "boundary_requires_extension",
  "formal 23-point CV decision"
)
extended_winner_count <- decision_uint(
  extended_decision, "aggregate_winner_count", "formal 23-point CV decision"
)
extended_finite_precision_plateau <- decision_bool(
  extended_decision, "finite_precision_plateau",
  "formal 23-point CV decision"
)
if (!identical(
      decision_value(
        extended_decision, "decision_basis", "formal 23-point CV decision"
      ),
      "lowest_median_raw_chisq_across_three_23_point_randomCV_runs"
    ) ||
    extended_boundary_requires_extension ||
    !identical(
      decision_uint(
        extended_decision, "formal_grid_points_per_run", "formal 23-point CV decision"
      ),
      23
    ) ||
    !identical(
      decision_uint(extended_decision, "formal_rows", "formal 23-point CV decision"),
      69
    ) ||
    !identical(
      decision_value(extended_decision, "formal_runs", "formal 23-point CV decision"),
      "cv_extendedgrid_run_01,cv_extendedgrid_run_02,cv_extendedgrid_run_03"
    ) ||
    !identical(
      decision_value(
        extended_decision, "formal_seed_values", "formal 23-point CV decision"
      ),
      "2026090441,2026090442,2026090443"
    )) {
  fail("Formal qa/cv_23point_decision.tsv does not describe the required 23-point design")
}
if (extended_finite_precision_plateau) {
  if (!identical(extended_status, "ACCEPTED_WITH_SENSITIVITIES") ||
      !extended_boundary || extended_winner_count < 3) {
    fail("Malformed qualified finite-precision plateau in formal 23-point decision")
  }
} else if (extended_boundary || !identical(extended_winner_count, 1)) {
  fail("Formal 23-point decision lacks a unique non-boundary winner")
}
extended_true_fields <- c(
  "formal_exit_codes_all_zero", "formal_stderr_all_empty",
  "formal_dated_trees_all_nonempty",
  "formal_dated_trees_valid_branch_length_newick",
  "formal_dated_tree_tip_sets_match_input",
  "formal_dated_tree_topologies_match_input",
  "formal_dated_trees_ultrametric_within_tolerance",
  "formal_config_contracts_validated", "formal_config_lints_all_pass",
  "formal_metadata_validated", "formal_stdout_score_sets_validated",
  "formal_resource_logs_validated", "source_19point_decision_validated",
  "diagnostic_19point_full_stage_validation_passed",
  "diagnostic_19point_all_prefix_points_present"
)
failed_extended_qualifications <- extended_true_fields[!vapply(
  extended_true_fields,
  function(key) decision_bool(extended_decision, key, "formal 23-point CV decision"),
  logical(1)
)]
if (length(failed_extended_qualifications) || decision_bool(
    extended_decision, "diagnostic_19point_runs_used_for_selection",
    "formal 23-point CV decision"
  )) {
  fail(
    "Formal 23-point CV decision failed qualification field(s): ",
    if (length(failed_extended_qualifications)) {
      paste(failed_extended_qualifications, collapse = ", ")
    } else {
      "diagnostic_19point_runs_used_for_selection"
    }
  )
}
extended_selected_smoothing <- decision_value(
  extended_decision, "selected_smoothing", "formal 23-point CV decision"
)
extended_selected_numeric <- suppressWarnings(as.numeric(extended_selected_smoothing))
extended_median_numeric <- suppressWarnings(as.numeric(decision_value(
  extended_decision, "selected_median_score", "formal 23-point CV decision"
)))
if (!nzchar(extended_selected_smoothing) ||
    !is.finite(extended_selected_numeric) || extended_selected_numeric <= 0 ||
    !is.finite(extended_median_numeric) || extended_median_numeric < 0) {
  fail("Formal 23-point decision has invalid selected smoothing or median score")
}
extended_sensitivity_required <- decision_bool(
  extended_decision, "sensitivity_required", "formal 23-point CV decision"
)
expected_extended_status <- if (extended_sensitivity_required) {
  "ACCEPTED_WITH_SENSITIVITIES"
} else {
  "ACCEPTED"
}
if (!identical(extended_status, expected_extended_status)) {
  fail("Formal 23-point status and sensitivity_required are inconsistent")
}
extended_sensitivity_text <- decision_value(
  extended_decision, "sensitivity_smoothing_values",
  "formal 23-point CV decision"
)
if (extended_sensitivity_required) {
  extended_sensitivity_tokens <- strsplit(
    extended_sensitivity_text, ",", fixed = TRUE
  )[[1L]]
  sensitivity_numeric <- suppressWarnings(as.numeric(extended_sensitivity_tokens))
  if (!nzchar(extended_sensitivity_text) ||
      any(!nzchar(extended_sensitivity_tokens)) ||
      !identical(extended_sensitivity_tokens, trimws(extended_sensitivity_tokens)) ||
      anyDuplicated(extended_sensitivity_tokens) ||
      any(!is.finite(sensitivity_numeric) | sensitivity_numeric <= 0)) {
    fail(
      "Formal 23-point sensitivity_smoothing_values must be unique, ",
      "comma-delimited positive numeric tokens without whitespace"
    )
  }
} else {
  if (nzchar(extended_sensitivity_text)) {
    fail(
      "Formal 23-point sensitivity_smoothing_values must be empty when ",
      "sensitivity_required=FALSE"
    )
  }
  extended_sensitivity_tokens <- character()
}

manifest <- read.delim(
  manifest_path,
  sep = "\t",
  header = TRUE,
  quote = "",
  comment.char = "",
  colClasses = "character",
  check.names = FALSE,
  na.strings = character(),
  strip.white = FALSE,
  stringsAsFactors = FALSE
)
manifest_required <- c(
  "variant_id", "variant_role", "smooth", "thorough", "seed",
  "config_path", "tree_path", "stage"
)
missing_manifest_columns <- setdiff(manifest_required, names(manifest))
if (length(missing_manifest_columns)) {
  fail(
    "Final manifest is missing required column(s): ",
    paste(missing_manifest_columns, collapse = ", ")
  )
}
if (!nrow(manifest)) fail("Final manifest has no variant rows")
manifest_required_data <- manifest[, manifest_required, drop = FALSE]
if (anyNA(manifest_required_data) ||
    any(vapply(manifest_required_data, function(x) any(!nzchar(x)), logical(1)))) {
  fail("Final manifest required fields must be non-empty and non-NA")
}
if (anyDuplicated(manifest$variant_id)) fail("Manifest variant_id values must be unique")
if (anyDuplicated(manifest$stage)) fail("Manifest stage values must be unique")
safe_id <- grepl("^[A-Za-z0-9_.-]+$", manifest$variant_id)
safe_stage <- grepl("^[A-Za-z0-9_.-]+$", manifest$stage)
if (!all(safe_id)) {
  fail("Unsafe manifest variant_id value(s): ", paste(manifest$variant_id[!safe_id], collapse = ", "))
}
if (!all(safe_stage)) {
  fail("Unsafe manifest stage value(s): ", paste(manifest$stage[!safe_stage], collapse = ", "))
}
allowed_variant_roles <- c(
  "primary", "optimization_diagnostic", "same_seed_repeat",
  "smoothing_sensitivity"
)
unknown_variant_roles <- setdiff(unique(manifest$variant_role), allowed_variant_roles)
if (length(unknown_variant_roles)) {
  fail(
    "Manifest contains disallowed variant_role value(s): ",
    paste(unknown_variant_roles, collapse = ", ")
  )
}
if (sum(manifest$variant_role == "primary") != 1L) {
  fail("Manifest must contain exactly one variant_role=primary row")
}
if (sum(manifest$variant_role == "optimization_diagnostic") != 1L) {
  fail("Manifest must contain exactly one variant_role=optimization_diagnostic row")
}
if (sum(manifest$variant_role == "same_seed_repeat") < 1L) {
  fail("Manifest must contain at least one variant_role=same_seed_repeat row")
}
if (!all(manifest$thorough %in% c("TRUE", "FALSE"))) {
  fail("Manifest thorough must use exact uppercase TRUE or FALSE values")
}
manifest$smooth_value <- suppressWarnings(as.numeric(manifest$smooth))
if (any(!is.finite(manifest$smooth_value) | manifest$smooth_value <= 0)) {
  fail("Manifest smooth values must be finite positive numbers")
}
if (any(!grepl("^[0-9]+$", manifest$seed))) {
  fail("Manifest seed values must be unsigned base-10 integers")
}
manifest$seed_value <- suppressWarnings(as.numeric(manifest$seed))
if (any(!is.finite(manifest$seed_value) | manifest$seed_value <= 0)) {
  fail("Manifest seed values must be positive and within the supported numeric range")
}

manifest$config_resolved <- vapply(
  seq_len(nrow(manifest)),
  function(i) resolve_run_file(
    manifest$config_path[[i]],
    paste0("config_path for manifest variant ", manifest$variant_id[[i]])
  ),
  character(1)
)
manifest$tree_resolved <- vapply(
  seq_len(nrow(manifest)),
  function(i) resolve_run_file(
    manifest$tree_path[[i]],
    paste0("tree_path for manifest variant ", manifest$variant_id[[i]])
  ),
  character(1)
)
if (any(dirname(manifest$config_resolved) != config_dir) ||
    any(!grepl("[.]cfg$", basename(manifest$config_resolved)))) {
  fail("Every manifest config_path must be a .cfg file directly under configs/")
}
if (anyDuplicated(manifest$config_resolved)) {
  fail("Manifest config_path values must resolve to unique files")
}
if (anyDuplicated(manifest$tree_resolved)) {
  fail("Manifest tree_path values must resolve to unique files")
}

primary <- manifest[manifest$variant_role == "primary", , drop = FALSE]
repeat_rows <- manifest[manifest$variant_role == "same_seed_repeat", , drop = FALSE]
sensitivity_rows <- manifest[
  manifest$variant_role == "smoothing_sensitivity", , drop = FALSE
]
if (extended_sensitivity_required) {
  if (nrow(sensitivity_rows) != length(extended_sensitivity_tokens) ||
      anyDuplicated(sensitivity_rows$smooth) ||
      !setequal(sensitivity_rows$smooth, extended_sensitivity_tokens)) {
    fail(
      "Manifest smoothing_sensitivity smooth set does not exactly match ",
      "cv_23point_decision.tsv sensitivity_smoothing_values"
    )
  }
} else if (nrow(sensitivity_rows) != 0L) {
  fail(
    "Manifest must contain zero smoothing_sensitivity rows when ",
    "sensitivity_required=FALSE"
  )
}
if (!identical(primary$smooth[[1]], extended_selected_smoothing)) {
  fail(
    "Manifest primary smooth does not exactly match the qualified 23-point selection: ",
    primary$smooth[[1]], " versus ", extended_selected_smoothing
  )
}
if (any(repeat_rows$smooth != primary$smooth[[1]]) ||
    any(repeat_rows$thorough != primary$thorough[[1]]) ||
    any(repeat_rows$seed != primary$seed[[1]])) {
  fail("Every same_seed_repeat row must exactly match primary smooth, thorough, and seed")
}

base_specs <- data.frame(
  stage = c(
    "prime",
    sprintf("cv_full_run_%02d", 1:3),
    sprintf("cv_finalgrid_run_%02d", 1:3),
    sprintf("cv_extendedgrid_run_%02d", 1:3)
  ),
  stage_group = c(
    "prime", rep("cv_full", 3L), rep("cv_finalgrid", 3L),
    rep("cv_extendedgrid", 3L)
  ),
  phase = c(
    "prime", rep("cv_full", 3L), rep("cv_finalgrid", 3L),
    rep("cv_extendedgrid", 3L)
  ),
  config_path = file.path(
    config_dir,
    c(
      "prime.cfg",
      sprintf("cv_full_run_%02d.cfg", 1:3),
      sprintf("cv_finalgrid_run_%02d.cfg", 1:3),
      sprintf("cv_extendedgrid_run_%02d.cfg", 1:3)
    )
  ),
  variant_id = c(
    "prime", sprintf("cv_full_run_%02d", 1:3),
    sprintf("cv_finalgrid_run_%02d", 1:3),
    sprintf("cv_extendedgrid_run_%02d", 1:3)
  ),
  manifest_smooth = NA_character_,
  manifest_thorough = NA_character_,
  manifest_seed = c(
    "2026090431",
    as.character(2026090441:2026090443),
    as.character(2026090441:2026090443),
    as.character(2026090441:2026090443)
  ),
  manifest_tree = NA_character_,
  stringsAsFactors = FALSE
)
final_specs <- data.frame(
  stage = manifest$stage,
  stage_group = "final",
  phase = ifelse(manifest$thorough == "TRUE", "final", "final_baseline"),
  config_path = manifest$config_resolved,
  variant_id = manifest$variant_id,
  manifest_smooth = manifest$smooth,
  manifest_thorough = manifest$thorough,
  manifest_seed = manifest$seed,
  manifest_tree = manifest$tree_resolved,
  stringsAsFactors = FALSE
)
specs <- rbind(base_specs, final_specs)
if (anyDuplicated(specs$stage)) fail("Stage names overlap across predefined and manifest stages")
if (anyDuplicated(specs$config_path)) fail("Config paths overlap across predefined and manifest stages")

missing_configs <- specs$config_path[!file.exists(specs$config_path)]
if (length(missing_configs)) {
  fail("Missing required config(s): ", paste(missing_configs, collapse = ", "))
}
specs$config_path <- vapply(specs$config_path, normalizePath, character(1), mustWork = TRUE)

actual_configs <- sort(list.files(
  config_dir,
  pattern = "[.]cfg$",
  full.names = TRUE,
  recursive = FALSE
))
actual_configs <- vapply(actual_configs, normalizePath, character(1), mustWork = TRUE)
expected_configs <- sort(unique(specs$config_path))
if (!setequal(actual_configs, expected_configs)) {
  unexpected <- setdiff(actual_configs, expected_configs)
  missing <- setdiff(expected_configs, actual_configs)
  fail(
    "The configs/ set is not exactly the predefined ten plus manifest final configs; ",
    "unexpected: ", if (length(unexpected)) paste(vapply(unexpected, relative_run_path, character(1)), collapse = ",") else "none",
    "; missing: ", if (length(missing)) paste(missing, collapse = ",") else "none"
  )
}

expected_treefile <- file.path(run_root, "input", "verified_495_tip_input.treefile")
expected_calibration_lines <- c(
  "mrca = CROWN_PAPILIONOIDEA Parnassius_apollo_ncbi Aeromachus_catocyanea_mydata",
  "min = CROWN_PAPILIONOIDEA 91.5046",
  "max = CROWN_PAPILIONOIDEA 100.8925",
  "mrca = CROWN_PAPILIONIDAE_NONBARONIA Parnassius_apollo_ncbi Graphium_cloanthus_mydata",
  "min = CROWN_PAPILIONIDAE_NONBARONIA 44.1968",
  "max = CROWN_PAPILIONIDAE_NONBARONIA 52.9473",
  "mrca = CROWN_HESPERIINAE Aeromachus_catocyanea_mydata Acada_biseriata_kawahara2023",
  "min = CROWN_HESPERIINAE 36.210861",
  "max = CROWN_HESPERIINAE 40.662537"
)
forbidden_exact_values <- c("101.4", "36.2", "40.7")

common_lint_checks <- c(
  "treefile_exact", "numsites_182682",
  "exactly_three_mrca", "exactly_three_min", "exactly_three_max",
  "root_mrca_exact", "root_bounds_exact",
  "papilionidae_mrca_exact", "papilionidae_bounds_exact",
  "hesperiinae_mrca_exact", "hesperiinae_bounds_exact",
  "forbidden_value_absent_101_4", "forbidden_value_absent_36_2",
  "forbidden_value_absent_40_7", "log_pen_absent"
)
phase_lint_checks <- function(phase) {
  extras <- switch(
    phase,
    prime = c(
      "prime_enabled", "cv_disabled", "thorough_disabled_for_prime",
      "no_outfile_for_prime"
    ),
    cv_full = c(
      "randomcv_enabled", "prime_disabled", "thorough_disabled_for_cv",
      "cv_grid_start", "cv_grid_stop", "cv_grid_multiplier", "cviter_3",
      "nthreads_10", "cvoutfile_within_run", "outfile_within_run"
    ),
    cv_finalgrid = c(
      "randomcv_enabled", "prime_disabled", "thorough_disabled_for_cv",
      "cv_grid_start", "cv_grid_stop", "cv_grid_multiplier", "cviter_3",
      "nthreads_10", "cvoutfile_within_run", "outfile_within_run"
    ),
    cv_extendedgrid = c(
      "randomcv_enabled", "prime_disabled", "thorough_disabled_for_cv",
      "cv_grid_start", "cv_grid_stop", "cv_grid_multiplier", "cviter_3",
      "nthreads_10", "cvoutfile_within_run", "outfile_within_run"
    ),
    final = c(
      "thorough_enabled", "prime_and_cv_disabled", "one_positive_smoothing",
      "outfile_within_run"
    ),
    final_baseline = c(
      "thorough_disabled_for_baseline", "prime_and_cv_disabled",
      "one_positive_smoothing", "outfile_within_run"
    ),
    fail("Unsupported phase: ", phase)
  )
  c(common_lint_checks, extras)
}

parse_config <- function(path) {
  raw <- readLines(path, warn = FALSE)
  active <- trimws(sub("#.*$", "", raw))
  active <- active[nzchar(active)]
  has_equals <- grepl("=", active, fixed = TRUE)
  key <- active
  key[has_equals] <- trimws(sub("=.*$", "", active[has_equals]))
  value <- rep("", length(active))
  value[has_equals] <- trimws(sub("^[^=]*=", "", active[has_equals]))
  list(
    active = active,
    key = key,
    value = value,
    values_for = function(k) value[key == k]
  )
}

read_strict_lint <- function(path, phase) {
  if (!file.exists(path) || isTRUE(file.info(path)$isdir)) {
    fail("Missing strict lint TSV: ", path)
  }
  lint <- read.delim(
    path,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  expected_columns <- c("check", "pass", "observed", "expected")
  if (!identical(names(lint), expected_columns)) {
    fail(
      "Malformed strict lint TSV columns in ", path, "; observed: ",
      paste(names(lint), collapse = ",")
    )
  }
  if (!nrow(lint) || anyNA(lint$check) || any(!nzchar(lint$check)) ||
      anyDuplicated(lint$check)) {
    fail("Strict lint TSV has empty or duplicate check identifiers: ", path)
  }
  expected_checks <- phase_lint_checks(phase)
  seed_check_present <- any(c(
    "explicit_positive_seed", "explicit_positive_integer_seed"
  ) %in% lint$check)
  required_checks_present <- all(expected_checks %in% lint$check) &&
    seed_check_present
  pass_text <- toupper(trimws(as.character(lint$pass)))
  all_pass <- length(pass_text) == nrow(lint) && !anyNA(pass_text) &&
    all(pass_text == "TRUE")
  list(
    table = lint,
    required_checks_present = required_checks_present,
    all_pass = all_pass,
    sha256 = sha256_file(path)
  )
}

read_metadata <- function(path) {
  if (!file.exists(path) || isTRUE(file.info(path)$isdir)) {
    fail("Missing stage metadata: ", path)
  }
  meta <- read.delim(
    path,
    sep = "\t",
    header = FALSE,
    quote = "",
    comment.char = "",
    col.names = c("key", "value"),
    colClasses = "character",
    fill = FALSE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  required <- c(
    "stage", "pid", "start_time", "binary", "config",
    "end_time", "exit_code", "elapsed_s"
  )
  if (ncol(meta) != 2L || !nrow(meta) ||
      any(vapply(required, function(k) sum(meta$key == k) != 1L, logical(1)))) {
    fail("Malformed or incomplete stage metadata: ", path)
  }
  value <- function(k) meta$value[meta$key == k][[1]]
  list(
    stage = value("stage"),
    pid = value("pid"),
    start_time = value("start_time"),
    binary = value("binary"),
    config = value("config"),
    end_time = value("end_time"),
    exit_code = value("exit_code"),
    elapsed_s = value("elapsed_s")
  )
}

check_rows <- vector("list", nrow(specs))
outfile_rows <- list()
add_output <- function(stage, config, directive, path) {
  outfile_rows[[length(outfile_rows) + 1L]] <<- data.frame(
    stage = stage,
    config = config,
    directive = directive,
    path = path,
    stringsAsFactors = FALSE
  )
}

for (i in seq_len(nrow(specs))) {
  spec <- specs[i, , drop = FALSE]
  cfg_path <- spec$config_path[[1]]
  cfg_name <- basename(cfg_path)
  cfg <- parse_config(cfg_path)
  values_for <- cfg$values_for
  active <- cfg$active
  key <- cfg$key
  value <- cfg$value

  exact_calibrations <- vapply(
    expected_calibration_lines,
    function(expected) sum(active == expected) == 1L,
    logical(1)
  )
  bound_values <- value[key %in% c("min", "max")]
  calibration_numeric_tokens <- sub("^.*[[:space:]]+", "", bound_values)

  checks <- c(
    exactly_three_mrca = sum(key == "mrca") == 3L,
    exactly_three_min = sum(key == "min") == 3L,
    exactly_three_max = sum(key == "max") == 3L,
    all_nine_calibration_lines_exact = all(exact_calibrations),
    numsites_exact = identical(values_for("numsites"), "182682"),
    treefile_exact = identical(values_for("treefile"), expected_treefile),
    log_pen_absent = !any(key == "log_pen"),
    forbidden_old_exact_values_absent =
      !any(calibration_numeric_tokens %in% forbidden_exact_values),
    no_rescale_or_reanchor = !any(grepl(
      "rescal|re[-_[:space:]]*anchor", active, ignore.case = TRUE
    )),
    one_exact_expected_seed = identical(values_for("seed"), spec$manifest_seed[[1]])
  )

  if (spec$phase[[1]] == "prime") {
    checks <- c(
      checks,
      phase_prime_exact = sum(key == "prime") == 1L &&
        !any(key %in% c("cv", "randomcv", "smooth", "thorough")),
      prime_has_no_outfile = !any(key %in% c("outfile", "cvoutfile"))
    )
  } else if (spec$phase[[1]] %in% c(
      "cv_full", "cv_finalgrid", "cv_extendedgrid"
    )) {
    expected_stop <- switch(
      spec$phase[[1]],
      cv_full = "0.000000001",
      cv_finalgrid = "0.0000000000001",
      cv_extendedgrid = "0.00000000000000001"
    )
    checks <- c(
      checks,
      phase_randomcv_exact = sum(key == "randomcv") == 1L &&
        !any(key %in% c("prime", "cv", "smooth", "thorough")),
      cv_grid_start_exact = identical(values_for("cvstart"), "100000"),
      cv_grid_stop_exact = identical(values_for("cvstop"), expected_stop),
      cv_grid_multiplier_exact = identical(values_for("cvmultstep"), "0.1"),
      cviter_exact = identical(values_for("cviter"), "3"),
      nthreads_exact = identical(values_for("nthreads"), "10"),
      cv_has_one_outfile = length(values_for("outfile")) == 1L,
      cv_has_one_cvoutfile = length(values_for("cvoutfile")) == 1L
    )
  } else {
    expected_thorough <- identical(spec$manifest_thorough[[1]], "TRUE")
    smooth_text <- values_for("smooth")
    smooth_numeric <- suppressWarnings(as.numeric(smooth_text))
    checks <- c(
      checks,
      phase_final_controls_absent = !any(key %in% c(
        "prime", "cv", "randomcv", "cvstart", "cvstop", "cvmultstep", "cvoutfile"
      )),
      phase_thorough_matches_manifest =
        sum(active == "thorough") == as.integer(expected_thorough),
      final_has_one_positive_smooth = length(smooth_numeric) == 1L &&
        is.finite(smooth_numeric) && smooth_numeric > 0,
      final_smooth_text_matches_manifest = identical(smooth_text, spec$manifest_smooth[[1]]),
      final_has_one_outfile = length(values_for("outfile")) == 1L
    )
  }

  config_stem <- sub("[.]cfg$", "", cfg_name)
  lint_path <- if (spec$stage_group[[1]] == "final") {
    file.path(qa_dir, paste0("config_lint_", config_stem, ".tsv"))
  } else {
    file.path(qa_dir, paste0("lint_", config_stem, ".strict.tsv"))
  }
  lint <- read_strict_lint(lint_path, spec$phase[[1]])
  checks <- c(
    checks,
    strict_lint_required_checks_present = lint$required_checks_present,
    strict_lint_all_pass = lint$all_pass
  )

  metadata_path <- file.path(log_dir, paste0(spec$stage[[1]], ".metadata.tsv"))
  metadata <- read_metadata(metadata_path)
  metadata_config <- tryCatch(
    resolve_run_file(metadata$config, paste0(spec$stage[[1]], " metadata config")),
    error = function(e) NA_character_
  )
  metadata_binary <- tryCatch(
    resolve_run_file(metadata$binary, paste0(spec$stage[[1]], " metadata binary")),
    error = function(e) NA_character_
  )
  checks <- c(
    checks,
    metadata_stage_exact = identical(metadata$stage, spec$stage[[1]]),
    metadata_config_exact = identical(metadata_config, cfg_path),
    metadata_binary_exact_patched_treepl =
      identical(metadata_binary, expected_treepl_binary),
    metadata_exit_zero = identical(metadata$exit_code, "0"),
    metadata_pid_positive_integer = grepl("^[0-9]+$", metadata$pid) &&
      suppressWarnings(as.numeric(metadata$pid)) > 0,
    metadata_elapsed_nonnegative_integer = grepl("^[0-9]+$", metadata$elapsed_s),
    metadata_has_start_and_end = nzchar(metadata$start_time) && nzchar(metadata$end_time)
  )

  stdout_path <- file.path(log_dir, paste0(spec$stage[[1]], ".stdout"))
  stderr_path <- file.path(log_dir, paste0(spec$stage[[1]], ".stderr"))
  resource_path <- file.path(log_dir, paste0(spec$stage[[1]], ".resources.txt"))
  stdout_nonempty <- file.exists(stdout_path) &&
    !isTRUE(file.info(stdout_path)$isdir) && file.info(stdout_path)$size > 0
  stderr_empty <- file.exists(stderr_path) &&
    !isTRUE(file.info(stderr_path)$isdir) && file.info(stderr_path)$size == 0
  resource_nonempty <- file.exists(resource_path) &&
    !isTRUE(file.info(resource_path)$isdir) && file.info(resource_path)$size > 0
  resource_lines <- if (resource_nonempty) {
    readLines(resource_path, warn = FALSE)
  } else {
    character()
  }
  resource_exit_zero <- sum(trimws(resource_lines) == "Exit status: 0") == 1L
  checks <- c(
    checks,
    wrapper_stdout_nonempty = stdout_nonempty,
    wrapper_stderr_empty = stderr_empty,
    wrapper_resource_nonempty = resource_nonempty,
    wrapper_resource_exactly_one_exit_status_zero = resource_exit_zero
  )

  outfile_values <- values_for("outfile")
  cvoutfile_values <- values_for("cvoutfile")
  outfile_resolved <- character()
  cvoutfile_resolved <- character()
  if (length(outfile_values) == 1L) {
    outfile_resolved <- tryCatch(
      resolve_run_file(outfile_values[[1]], paste0(spec$stage[[1]], " outfile")),
      error = function(e) NA_character_
    )
    checks <- c(
      checks,
      configured_outfile_exists_nonempty = length(outfile_resolved) == 1L &&
        !is.na(outfile_resolved) && file.info(outfile_resolved)$size > 0
    )
    if (!is.na(outfile_resolved)) {
      add_output(spec$stage[[1]], cfg_name, "outfile", outfile_resolved)
    }
  }
  if (length(cvoutfile_values) == 1L) {
    cvoutfile_resolved <- tryCatch(
      resolve_run_file(cvoutfile_values[[1]], paste0(spec$stage[[1]], " cvoutfile")),
      error = function(e) NA_character_
    )
    checks <- c(
      checks,
      configured_cvoutfile_exists_nonempty = length(cvoutfile_resolved) == 1L &&
        !is.na(cvoutfile_resolved) && file.info(cvoutfile_resolved)$size > 0
    )
    if (!is.na(cvoutfile_resolved)) {
      add_output(spec$stage[[1]], cfg_name, "cvoutfile", cvoutfile_resolved)
    }
  }
  if (spec$stage_group[[1]] == "final") {
    checks <- c(
      checks,
      final_outfile_matches_manifest_tree = length(outfile_resolved) == 1L &&
        identical(outfile_resolved, spec$manifest_tree[[1]])
    )
  }

  failed <- names(checks)[is.na(checks) | !checks]
  check_rows[[i]] <- data.frame(
    config = cfg_name,
    phase = spec$phase[[1]],
    stage = spec$stage[[1]],
    stage_group = spec$stage_group[[1]],
    variant_id = spec$variant_id[[1]],
    config_sha256 = sha256_file(cfg_path),
    config_hash_timing = "post_run_collection_snapshot_not_launch_time_attestation",
    lint_table = basename(lint_path),
    lint_sha256 = lint$sha256,
    lint_pass = sum(toupper(trimws(as.character(lint$table$pass))) == "TRUE", na.rm = TRUE),
    lint_total = nrow(lint$table),
    outfile = if (length(outfile_resolved) == 1L && !is.na(outfile_resolved)) {
      relative_run_path(outfile_resolved)
    } else {
      ""
    },
    cvoutfile = if (length(cvoutfile_resolved) == 1L && !is.na(cvoutfile_resolved)) {
      relative_run_path(cvoutfile_resolved)
    } else {
      ""
    },
    status = if (!length(failed)) "PASS" else "FAIL",
    failed_consolidated_checks = paste(failed, collapse = ","),
    stringsAsFactors = FALSE
  )
}

config_table <- do.call(rbind, check_rows)
output_table <- if (length(outfile_rows)) {
  do.call(rbind, outfile_rows)
} else {
  data.frame(stage = character(), config = character(), directive = character(), path = character())
}

expected_outfile_count <- 9L + nrow(manifest)
expected_cvoutfile_count <- 9L
global_checks <- c(
  manifest_present = file.exists(manifest_path),
  manifest_roles_allowed = all(manifest$variant_role %in% allowed_variant_roles),
  manifest_exactly_one_primary = sum(manifest$variant_role == "primary") == 1L,
  manifest_exactly_one_optimization_diagnostic =
    sum(manifest$variant_role == "optimization_diagnostic") == 1L,
  manifest_at_least_one_same_seed_repeat =
    sum(manifest$variant_role == "same_seed_repeat") >= 1L,
  manifest_smoothing_sensitivity_set_exact = if (extended_sensitivity_required) {
    nrow(sensitivity_rows) == length(extended_sensitivity_tokens) &&
      !anyDuplicated(sensitivity_rows$smooth) &&
      setequal(sensitivity_rows$smooth, extended_sensitivity_tokens)
  } else {
    nrow(sensitivity_rows) == 0L
  },
  source_19point_status_requires_extension = identical(
    decision_value(source_decision, "status", "source 19-point CV decision"),
    "REQUIRES_GRID_EXTENSION"
  ),
  formal_23point_decision_qualified = extended_status %in%
    c("ACCEPTED", "ACCEPTED_WITH_SENSITIVITIES"),
  predefined_stage_count_exact = sum(specs$stage_group != "final") == 10L,
  extendedgrid_stage_count_exact = sum(specs$stage_group == "cv_extendedgrid") == 3L,
  manifest_final_stage_count_exact = sum(specs$stage_group == "final") == nrow(manifest),
  all_stage_names_unique = !anyDuplicated(specs$stage),
  expected_config_set_exact = setequal(actual_configs, expected_configs),
  all_config_rows_pass = all(config_table$status == "PASS"),
  outfile_count_exact = sum(output_table$directive == "outfile") == expected_outfile_count,
  cvoutfile_count_exact = sum(output_table$directive == "cvoutfile") == expected_cvoutfile_count,
  outfile_paths_unique = !anyDuplicated(output_table$path[output_table$directive == "outfile"]),
  cvoutfile_paths_unique = !anyDuplicated(output_table$path[output_table$directive == "cvoutfile"]),
  outfile_and_cvoutfile_paths_globally_unique = !anyDuplicated(output_table$path)
)
overall_pass <- all(global_checks) && all(config_table$status == "PASS")

report <- c(
  "treePL three-calibration configuration lint — final consolidated report",
  paste0("run_root\t", run_root),
  paste0("generated_at\t", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  paste0("manifest\t", relative_run_path(manifest_path)),
  paste0(
    "required_patched_treepl_binary\t",
    relative_run_path(expected_treepl_binary),
    "\tsha256=", sha256_file(expected_treepl_binary)
  ),
  paste0(
    "source_19point_decision\t", relative_run_path(source_decision_path),
    "\tstatus=", decision_value(
      source_decision, "status", "source 19-point CV decision"
    ),
    "\tsha256=", sha256_file(source_decision_path)
  ),
  paste0(
    "formal_23point_decision\t", relative_run_path(extended_decision_path),
    "\tstatus=", extended_status,
    "\tselected_smoothing=", extended_selected_smoothing,
    "\tsha256=", sha256_file(extended_decision_path)
  ),
  paste0("manifest_final_variants\t", nrow(manifest)),
  paste0("configs_checked\t", nrow(config_table)),
  paste0("stages_checked\t", nrow(specs)),
  paste0("overall_status\t", if (overall_pass) "PASS" else "FAIL"),
  paste0(
    "config_hash_interpretation\t",
    "SHA-256 values were computed after completed stages; they are reproducibility snapshots, ",
    "not launch-time attestations because wrapper metadata does not record hashes"
  ),
  "",
  "GLOBAL_CHECK\tSTATUS",
  paste0(names(global_checks), "\t", ifelse(global_checks, "PASS", "FAIL")),
  "",
  paste(names(config_table), collapse = "\t"),
  vapply(
    seq_len(nrow(config_table)),
    function(i) paste(config_table[i, , drop = TRUE], collapse = "\t"),
    character(1)
  )
)

# Exclusive creation protects the no-overwrite contract even if another process
# creates the report between the initial existence check and this write.
connection <- file(out_path, open = "wx", encoding = "UTF-8")
writeLines(report, connection, useBytes = TRUE)
close(connection)

cat(
  "Consolidated three-calibration config/stage lint:",
  if (overall_pass) "PASS" else "FAIL",
  "-", nrow(config_table), "configs and", nrow(specs), "stages\n"
)
if (!overall_pass) quit(status = 2L)
