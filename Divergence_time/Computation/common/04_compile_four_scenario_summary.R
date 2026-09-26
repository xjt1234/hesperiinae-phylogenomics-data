#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, digits = 15)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop(
    "Usage: 04_compile_four_scenario_summary.R <validation_root>",
    call. = FALSE
  )
}

validation_root <- normalizePath(args[[1L]], mustWork = TRUE)
expected_root_name <- "treePL_R2_3_minimal_validation_20260904_135356"
if (!identical(basename(validation_root), expected_root_name)) {
  stop(
    "Validation-root basename mismatch: observed=", basename(validation_root),
    " expected=", expected_root_name, call. = FALSE
  )
}

qa_root <- file.path(validation_root, "qa")
if (!dir.exists(qa_root)) stop("Missing root QA directory", call. = FALSE)

output_four <- file.path(
  validation_root, "calibration_sensitivity_summary_with_leave_one_out.tsv"
)
output_diagnostic <- file.path(qa_root, "leave_one_out_diagnostic_summary.tsv")
output_checks <- file.path(qa_root, "four_scenario_validation_checks.tsv")
output_report <- file.path(qa_root, "four_scenario_report.md")
output_paths <- c(output_four, output_diagnostic, output_checks, output_report)
existing_outputs <- output_paths[file.exists(output_paths)]
if (length(existing_outputs)) {
  stop(
    "Refusing to overwrite existing output(s): ",
    paste(existing_outputs, collapse = ", "), call. = FALSE
  )
}

EXPECTED_INPUT_SHA256 <- paste0(
  "4d35d553ed5bbf08022daef2361ffa0c",
  "342f8694a53aaf08d5350692094439e5"
)
EXPECTED_TREEPL_SHA256 <- paste0(
  "5d6bfac8d18d835f24ebb79f0373b070a",
  "948b58a179f437f0aacb8fed74ccc5d"
)
EXPECTED_TIPS <- 495L
EXPECTED_NUMSITES <- "182682"
scenario_order <- c(
  "A_two_node_main",
  "B_three_node_sensitivity",
  "C_root_only",
  "D_papilionidae_only"
)
loo_scenarios <- c("C_root_only", "D_papilionidae_only")
node_order <- c(
  "CROWN_PAPILIONOIDEA",
  "CROWN_PAPILIONIDAE_NONBARONIA",
  "CROWN_HESPERIINAE"
)
calibrated_node_by_loo <- c(
  C_root_only = "CROWN_PAPILIONOIDEA",
  D_papilionidae_only = "CROWN_PAPILIONIDAE_NONBARONIA"
)

read_tsv_character <- function(path, required_columns = character()) {
  if (!file.exists(path) || file.info(path)$size <= 0) {
    stop("Missing/empty input table: ", path, call. = FALSE)
  }
  table <- read.delim(
    path, sep = "\t", quote = "", comment.char = "",
    check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character",
    na.strings = c("NA", "")
  )
  missing <- setdiff(required_columns, names(table))
  if (length(missing)) {
    stop(
      "Missing column(s) in ", path, ": ", paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  table
}

parse_bool <- function(x) {
  ifelse(
    is.na(x), NA,
    ifelse(toupper(trimws(x)) %in% c("TRUE", "T", "1"), TRUE,
           ifelse(toupper(trimws(x)) %in% c("FALSE", "F", "0"), FALSE, NA))
  )
}

as_finite_numeric <- function(x, label) {
  value <- suppressWarnings(as.numeric(x))
  if (any(!is.finite(value))) {
    stop("Non-finite numeric value in ", label, call. = FALSE)
  }
  value
}

format_number <- function(x) {
  if (length(x) != 1L || is.na(x) || !is.finite(x)) return("NA")
  format(x, digits = 15, scientific = FALSE, trim = TRUE)
}

parse_numeric_csv <- function(x) {
  if (length(x) != 1L || is.na(x) || !nzchar(trimws(x)) ||
      identical(tolower(trimws(x)), "none")) {
    return(numeric())
  }
  values <- suppressWarnings(as.numeric(trimws(strsplit(x, ",", fixed = TRUE)[[1L]])))
  if (any(!is.finite(values))) stop("Malformed numeric CSV: ", x, call. = FALSE)
  values
}

numeric_set_equal <- function(a, b, tolerance = 1e-12) {
  if (length(a) != length(b)) return(FALSE)
  if (!length(a)) return(TRUE)
  a <- sort(a)
  b <- sort(b)
  all(abs(a - b) <= tolerance * pmax(abs(a), abs(b), .Machine$double.xmin))
}

checks_list <- list()
add_check <- function(check_id, scope, passed, observed, expected, note = "") {
  passed <- isTRUE(passed)
  checks_list[[length(checks_list) + 1L]] <<- data.frame(
    check_id = check_id,
    scope = scope,
    critical = TRUE,
    status = if (passed) "PASS" else "FAIL",
    observed = as.character(observed),
    expected = as.character(expected),
    note = as.character(note),
    stringsAsFactors = FALSE
  )
}

ab_path <- file.path(validation_root, "calibration_sensitivity_summary.tsv")
summary_columns <- c(
  "scenario", "input_tree_path", "input_tree_sha256", "final_tree_path",
  "final_tree_sha256", "tip_count", "numsites", "smoothing", "root_age_ma",
  "papilionidae_age_ma", "hesperiinae_age_ma",
  "hesperiinae_shift_from_main_ma", "root_boundary_status",
  "papilionidae_boundary_status", "hesperiinae_boundary_status",
  "ultrametric_range_ma", "negative_branch_count", "zero_branch_count",
  "nonfinite_branch_count", "rooted_rf", "qa_status"
)
ab <- read_tsv_character(ab_path, summary_columns)
if (!identical(names(ab), summary_columns)) {
  stop(
    "A/B summary columns/order changed; refusing to reinterpret source table",
    call. = FALSE
  )
}
add_check(
  "ab_exact_scenarios_and_order", "A_B_source",
  nrow(ab) == 2L && identical(ab$scenario, scenario_order[1:2]),
  paste(ab$scenario, collapse = ","), paste(scenario_order[1:2], collapse = ",")
)
add_check(
  "ab_input_hash_fixed", "A_B_source",
  all(ab$input_tree_sha256 == EXPECTED_INPUT_SHA256),
  paste(unique(ab$input_tree_sha256), collapse = ","), EXPECTED_INPUT_SHA256
)
add_check(
  "ab_tip_count_495_numsites_182682", "A_B_source",
  all(ab$tip_count == as.character(EXPECTED_TIPS)) &&
    all(ab$numsites == EXPECTED_NUMSITES),
  paste0("tips=", paste(unique(ab$tip_count), collapse = ","),
         ";numsites=", paste(unique(ab$numsites), collapse = ",")),
  "tips=495;numsites=182682"
)
add_check(
  "ab_topology_branch_and_qa_pass", "A_B_source",
  all(as_finite_numeric(ab$rooted_rf, "A/B rooted_rf") == 0) &&
    all(as_finite_numeric(ab$negative_branch_count, "A/B negative branches") == 0) &&
    all(as_finite_numeric(ab$zero_branch_count, "A/B zero branches") == 0) &&
    all(as_finite_numeric(ab$nonfinite_branch_count, "A/B nonfinite branches") == 0) &&
    all(ab$qa_status == "PASS"),
  paste(ab$qa_status, collapse = ","), "all topology/branch gates zero and QA PASS"
)

a_hesperiinae_age <- as_finite_numeric(
  ab$hesperiinae_age_ma[ab$scenario == "A_two_node_main"],
  "A Hesperiinae age"
)
if (length(a_hesperiinae_age) != 1L) {
  stop("Could not identify exactly one A Hesperiinae age", call. = FALSE)
}

loo_data <- list()
for (scenario in loo_scenarios) {
  scenario_qa <- file.path(validation_root, scenario, "qa")
  paths <- list(
    node_ages = file.path(scenario_qa, "node_ages.tsv"),
    variant_checks = file.path(scenario_qa, "final_variant_checks.tsv"),
    decision = file.path(scenario_qa, "cv_decision.tsv"),
    variant_ages = file.path(scenario_qa, "final_variant_node_ages.tsv")
  )
  node_ages <- read_tsv_character(paths$node_ages, c(
    "scenario", "variant_role", "stage", "smoothing", "node_id",
    "descendant_tip_count", "expected_descendant_tip_count",
    "descendant_count_exact", "calibrated_in_scenario", "node_age_ma",
    "calibrated_node_within_bounds_1e_4", "boundary_status_1e_4"
  ))
  variant_checks <- read_tsv_character(paths$variant_checks, c(
    "scenario", "variant_role", "stage", "smoothing", "output_path",
    "output_sha256", "tip_count", "unique_tip_count", "rooted",
    "fully_binary", "tip_set_exact", "rooted_rf", "rooted_clades_exact",
    "finite_branch_lengths", "negative_branch_count", "zero_branch_count",
    "root_to_tip_range_ma", "root_to_tip_range_below_1e_4",
    "focal_anchor_descendant_counts_exact", "exactly_one_scenario_calibration",
    "calibrated_node_within_bounds_1e_4",
    "metadata_stage_config_binary_exit0", "stderr_empty",
    "stdout_nonempty_no_failure_token", "resources_exit_zero",
    "all_hard_gates_pass"
  ))
  decision <- read_tsv_character(paths$decision, c(
    "scenario", "status", "selected_smoothing", "sensitivity_required",
    "sensitivity_smoothing_values", "input_tree_sha256", "treepl_binary_sha256"
  ))
  variant_ages <- read_tsv_character(paths$variant_ages, c(
    "scenario", "variant_role", "stage", "smoothing", "node_id",
    "descendant_tip_count", "expected_descendant_tip_count",
    "descendant_count_exact", "calibrated_in_scenario", "node_age_ma",
    "calibrated_node_within_bounds_1e_4", "boundary_status_1e_4"
  ))

  if (nrow(decision) != 1L) {
    stop("Expected one CV decision row for ", scenario, call. = FALSE)
  }
  primary_checks <- variant_checks[
    variant_checks$variant_role == "primary" &
      variant_checks$stage == "final_primary_thorough", , drop = FALSE
  ]
  repeat_checks <- variant_checks[
    variant_checks$variant_role == "same_seed_repeat" &
      variant_checks$stage == "final_repeat_same_seed", , drop = FALSE
  ]
  primary_nodes <- node_ages[
    node_ages$variant_role == "primary" &
      node_ages$stage == "final_primary_thorough", , drop = FALSE
  ]

  add_check(
    "loo_source_scenario_labels_exact", scenario,
    all(node_ages$scenario == scenario) && all(variant_checks$scenario == scenario) &&
      identical(decision$scenario[[1L]], scenario) &&
      all(variant_ages$scenario == scenario),
    scenario, scenario
  )
  add_check(
    "loo_primary_repeat_unique", scenario,
    nrow(primary_checks) == 1L && nrow(repeat_checks) == 1L,
    paste0("primary=", nrow(primary_checks), ";repeat=", nrow(repeat_checks)),
    "primary=1;repeat=1"
  )
  add_check(
    "loo_primary_three_nodes_exact", scenario,
    nrow(primary_nodes) == 3L && setequal(primary_nodes$node_id, node_order),
    paste(primary_nodes$node_id, collapse = ","), paste(node_order, collapse = ",")
  )
  add_check(
    "loo_fixed_hashes_exact", scenario,
    identical(decision$input_tree_sha256[[1L]], EXPECTED_INPUT_SHA256) &&
      identical(decision$treepl_binary_sha256[[1L]], EXPECTED_TREEPL_SHA256),
    paste(decision$input_tree_sha256[[1L]], decision$treepl_binary_sha256[[1L]], sep = ";"),
    paste(EXPECTED_INPUT_SHA256, EXPECTED_TREEPL_SHA256, sep = ";")
  )
  add_check(
    "loo_decision_accepted", scenario,
    decision$status[[1L]] %in% c("ACCEPTED", "ACCEPTED_WITH_SENSITIVITIES"),
    decision$status[[1L]], "ACCEPTED or ACCEPTED_WITH_SENSITIVITIES"
  )
  add_check(
    "loo_all_variant_hard_gates_pass", scenario,
    nrow(variant_checks) >= 2L &&
      all(parse_bool(variant_checks$all_hard_gates_pass) %in% TRUE),
    paste(variant_checks$all_hard_gates_pass, collapse = ","), "all TRUE"
  )
  add_check(
    "loo_all_variants_tip_topology_exact", scenario,
    all(variant_checks$tip_count == as.character(EXPECTED_TIPS)) &&
      all(variant_checks$unique_tip_count == as.character(EXPECTED_TIPS)) &&
      all(parse_bool(variant_checks$rooted) %in% TRUE) &&
      all(parse_bool(variant_checks$fully_binary) %in% TRUE) &&
      all(parse_bool(variant_checks$tip_set_exact) %in% TRUE) &&
      all(as_finite_numeric(variant_checks$rooted_rf, paste(scenario, "RF")) == 0) &&
      all(parse_bool(variant_checks$rooted_clades_exact) %in% TRUE),
    paste0("variants=", nrow(variant_checks)),
    "all 495-tip rooted binary; exact tips/clades; RF=0"
  )
  add_check(
    "loo_all_variants_branch_and_ultrametric_gates", scenario,
    all(parse_bool(variant_checks$finite_branch_lengths) %in% TRUE) &&
      all(as_finite_numeric(variant_checks$negative_branch_count, paste(scenario, "negative")) == 0) &&
      all(as_finite_numeric(variant_checks$zero_branch_count, paste(scenario, "zero")) == 0) &&
      all(as_finite_numeric(variant_checks$root_to_tip_range_ma, paste(scenario, "root-to-tip")) < 1e-4) &&
      all(parse_bool(variant_checks$root_to_tip_range_below_1e_4) %in% TRUE),
    paste0("variants=", nrow(variant_checks)),
    "finite; negative=0; zero=0; root-to-tip range <1e-4 Ma"
  )
  add_check(
    "loo_all_variants_execution_provenance_pass", scenario,
    all(parse_bool(variant_checks$metadata_stage_config_binary_exit0) %in% TRUE) &&
      all(parse_bool(variant_checks$stderr_empty) %in% TRUE) &&
      all(parse_bool(variant_checks$stdout_nonempty_no_failure_token) %in% TRUE) &&
      all(parse_bool(variant_checks$resources_exit_zero) %in% TRUE),
    paste0("variants=", nrow(variant_checks)),
    "metadata/config/binary/exit0, empty stderr, valid stdout/resources"
  )
  add_check(
    "loo_same_seed_repeat_byte_identical", scenario,
    nrow(primary_checks) == 1L && nrow(repeat_checks) == 1L &&
      identical(primary_checks$output_sha256[[1L]], repeat_checks$output_sha256[[1L]]),
    if (nrow(primary_checks) == 1L && nrow(repeat_checks) == 1L) {
      paste(primary_checks$output_sha256[[1L]], repeat_checks$output_sha256[[1L]], sep = " vs ")
    } else "missing primary/repeat",
    "identical SHA-256"
  )

  selected <- as_finite_numeric(decision$selected_smoothing[[1L]], paste(scenario, "selected smooth"))
  sensitivity_declared <- parse_numeric_csv(decision$sensitivity_smoothing_values[[1L]])
  sensitivity_observed <- as_finite_numeric(
    variant_checks$smoothing[variant_checks$variant_role == "sensitivity"],
    paste(scenario, "observed sensitivity smoothings")
  )
  add_check(
    "loo_primary_smoothing_matches_cv", scenario,
    nrow(primary_checks) == 1L &&
      numeric_set_equal(as.numeric(primary_checks$smoothing[[1L]]), selected),
    if (nrow(primary_checks) == 1L) primary_checks$smoothing[[1L]] else "missing",
    decision$selected_smoothing[[1L]]
  )
  add_check(
    "loo_sensitivity_set_matches_cv", scenario,
    numeric_set_equal(sensitivity_observed, sensitivity_declared),
    paste(sort(sensitivity_observed), collapse = ","),
    paste(sort(sensitivity_declared), collapse = ",")
  )

  calibrated_node <- calibrated_node_by_loo[[scenario]]
  expected_calibrated <- primary_nodes$node_id == calibrated_node
  calibrated_flags <- parse_bool(primary_nodes$calibrated_in_scenario)
  expected_status <- ifelse(expected_calibrated, NA_character_, "not_calibrated")
  add_check(
    "loo_primary_calibration_roles_exact", scenario,
    identical(calibrated_flags, expected_calibrated) &&
      all(primary_nodes$boundary_status_1e_4[!expected_calibrated] == "not_calibrated") &&
      all(parse_bool(primary_nodes$calibrated_node_within_bounds_1e_4) %in% TRUE),
    paste(paste(primary_nodes$node_id, primary_nodes$calibrated_in_scenario,
                primary_nodes$boundary_status_1e_4, sep = ":"), collapse = ";"),
    paste0(calibrated_node, " only; all others not_calibrated")
  )

  variant_keys <- paste(variant_checks$scenario, variant_checks$stage, sep = "\034")
  age_keys <- unique(paste(variant_ages$scenario, variant_ages$stage, sep = "\034"))
  add_check(
    "loo_variant_age_stage_set_exact", scenario,
    setequal(variant_keys, age_keys) &&
      nrow(variant_ages) == 3L * nrow(variant_checks),
    paste0("check_stages=", length(unique(variant_keys)),
           ";age_rows=", nrow(variant_ages)),
    paste0("same stages; age_rows=", 3L * nrow(variant_checks))
  )

  loo_data[[scenario]] <- list(
    node_ages = node_ages,
    variant_checks = variant_checks,
    decision = decision,
    variant_ages = variant_ages,
    primary_checks = primary_checks,
    primary_nodes = primary_nodes,
    paths = paths
  )
}

loo_summary_rows <- lapply(loo_scenarios, function(scenario) {
  dat <- loo_data[[scenario]]
  pc <- dat$primary_checks
  pn <- dat$primary_nodes
  node_value <- function(node_id, column) {
    value <- pn[[column]][pn$node_id == node_id]
    if (length(value) != 1L) {
      stop("Missing/non-unique primary node value: ", scenario, " ", node_id)
    }
    value[[1L]]
  }
  finite_flag <- parse_bool(pc$finite_branch_lengths[[1L]])
  row <- as.list(setNames(rep(NA_character_, length(summary_columns)), summary_columns))
  row$scenario <- scenario
  row$input_tree_path <- ab$input_tree_path[[1L]]
  row$input_tree_sha256 <- dat$decision$input_tree_sha256[[1L]]
  row$final_tree_path <- pc$output_path[[1L]]
  row$final_tree_sha256 <- pc$output_sha256[[1L]]
  row$tip_count <- pc$tip_count[[1L]]
  row$numsites <- EXPECTED_NUMSITES
  row$smoothing <- pc$smoothing[[1L]]
  row$root_age_ma <- node_value("CROWN_PAPILIONOIDEA", "node_age_ma")
  row$papilionidae_age_ma <- node_value(
    "CROWN_PAPILIONIDAE_NONBARONIA", "node_age_ma"
  )
  row$hesperiinae_age_ma <- node_value("CROWN_HESPERIINAE", "node_age_ma")
  row$hesperiinae_shift_from_main_ma <- format_number(
    as.numeric(row$hesperiinae_age_ma) - a_hesperiinae_age
  )
  row$root_boundary_status <- node_value(
    "CROWN_PAPILIONOIDEA", "boundary_status_1e_4"
  )
  row$papilionidae_boundary_status <- node_value(
    "CROWN_PAPILIONIDAE_NONBARONIA", "boundary_status_1e_4"
  )
  row$hesperiinae_boundary_status <- node_value(
    "CROWN_HESPERIINAE", "boundary_status_1e_4"
  )
  row$ultrametric_range_ma <- pc$root_to_tip_range_ma[[1L]]
  row$negative_branch_count <- pc$negative_branch_count[[1L]]
  row$zero_branch_count <- pc$zero_branch_count[[1L]]
  row$nonfinite_branch_count <- if (isTRUE(finite_flag)) "0" else "NA"
  row$rooted_rf <- pc$rooted_rf[[1L]]
  row$qa_status <- if (isTRUE(parse_bool(pc$all_hard_gates_pass[[1L]])) &&
    all(parse_bool(dat$variant_checks$all_hard_gates_pass) %in% TRUE)) {
    "PASS"
  } else "FAIL"
  as.data.frame(row, stringsAsFactors = FALSE, check.names = FALSE)
})
loo_summary <- do.call(rbind, loo_summary_rows)
four_summary <- rbind(ab, loo_summary)
rownames(four_summary) <- NULL

add_check(
  "four_scenario_order_exact", "four_scenario_summary",
  identical(four_summary$scenario, scenario_order),
  paste(four_summary$scenario, collapse = ","), paste(scenario_order, collapse = ",")
)
add_check(
  "ab_rows_preserved_character_exact", "four_scenario_summary",
  identical(four_summary[1:2, , drop = FALSE], ab),
  "first two rows compared cell-for-cell", "identical to source A/B table"
)
add_check(
  "four_scenario_fixed_input_hash", "four_scenario_summary",
  all(four_summary$input_tree_sha256 == EXPECTED_INPUT_SHA256),
  paste(unique(four_summary$input_tree_sha256), collapse = ","), EXPECTED_INPUT_SHA256
)
add_check(
  "four_scenario_primary_tip_topology_branch_qa", "four_scenario_summary",
  all(four_summary$tip_count == as.character(EXPECTED_TIPS)) &&
    all(as_finite_numeric(four_summary$rooted_rf, "four-scenario RF") == 0) &&
    all(as_finite_numeric(four_summary$negative_branch_count, "four-scenario negative") == 0) &&
    all(as_finite_numeric(four_summary$zero_branch_count, "four-scenario zero") == 0) &&
    all(as_finite_numeric(four_summary$nonfinite_branch_count, "four-scenario nonfinite") == 0) &&
    all(four_summary$qa_status == "PASS"),
  paste(four_summary$qa_status, collapse = ","),
  "all 495 tips, RF=0, branch counts=0, QA PASS"
)
add_check(
  "loo_unconstrained_boundaries_marked_not_calibrated", "four_scenario_summary",
  four_summary$papilionidae_boundary_status[four_summary$scenario == "C_root_only"] ==
      "not_calibrated" &&
    four_summary$hesperiinae_boundary_status[four_summary$scenario == "C_root_only"] ==
      "not_calibrated" &&
    four_summary$root_boundary_status[four_summary$scenario == "D_papilionidae_only"] ==
      "not_calibrated" &&
    four_summary$hesperiinae_boundary_status[four_summary$scenario == "D_papilionidae_only"] ==
      "not_calibrated",
  "C: pap+hesp; D: root+hesp", "not_calibrated"
)

diagnostic_rows <- list()
for (scenario in loo_scenarios) {
  dat <- loo_data[[scenario]]
  vc <- dat$variant_checks
  va <- dat$variant_ages
  ages_by_node <- lapply(node_order, function(node_id) {
    as_finite_numeric(va$node_age_ma[va$node_id == node_id], paste(scenario, node_id))
  })
  names(ages_by_node) <- node_order
  ranges <- lapply(ages_by_node, function(x) c(min = min(x), max = max(x), range = diff(range(x))))

  for (i in seq_len(nrow(vc))) {
    stage <- vc$stage[[i]]
    rows <- va[va$stage == stage, , drop = FALSE]
    age_one <- function(node_id) {
      x <- as_finite_numeric(rows$node_age_ma[rows$node_id == node_id],
                             paste(scenario, stage, node_id))
      if (length(x) != 1L) stop("Non-unique diagnostic node age", call. = FALSE)
      x
    }
    diagnostic_rows[[length(diagnostic_rows) + 1L]] <- data.frame(
      scenario = scenario,
      variant_role = vc$variant_role[[i]],
      stage = stage,
      smoothing = vc$smoothing[[i]],
      output_path = vc$output_path[[i]],
      output_sha256 = vc$output_sha256[[i]],
      root_age_ma = age_one("CROWN_PAPILIONOIDEA"),
      papilionidae_age_ma = age_one("CROWN_PAPILIONIDAE_NONBARONIA"),
      hesperiinae_age_ma = age_one("CROWN_HESPERIINAE"),
      root_age_min_across_final_variants_ma = ranges[["CROWN_PAPILIONOIDEA"]][["min"]],
      root_age_max_across_final_variants_ma = ranges[["CROWN_PAPILIONOIDEA"]][["max"]],
      root_age_range_across_final_variants_ma = ranges[["CROWN_PAPILIONOIDEA"]][["range"]],
      papilionidae_age_min_across_final_variants_ma = ranges[["CROWN_PAPILIONIDAE_NONBARONIA"]][["min"]],
      papilionidae_age_max_across_final_variants_ma = ranges[["CROWN_PAPILIONIDAE_NONBARONIA"]][["max"]],
      papilionidae_age_range_across_final_variants_ma = ranges[["CROWN_PAPILIONIDAE_NONBARONIA"]][["range"]],
      hesperiinae_age_min_across_final_variants_ma = ranges[["CROWN_HESPERIINAE"]][["min"]],
      hesperiinae_age_max_across_final_variants_ma = ranges[["CROWN_HESPERIINAE"]][["max"]],
      hesperiinae_age_range_across_final_variants_ma = ranges[["CROWN_HESPERIINAE"]][["range"]],
      ultrametric_range_ma = as.numeric(vc$root_to_tip_range_ma[[i]]),
      rooted_rf = as.numeric(vc$rooted_rf[[i]]),
      all_hard_gates_pass = parse_bool(vc$all_hard_gates_pass[[i]]),
      diagnostic_extreme_flag = scenario == "D_papilionidae_only" &&
        vc$stage[[i]] == "final_sensitivity_1e-6",
      interpretation_role = "leave-one-out calibration diagnostic only",
      stringsAsFactors = FALSE
    )
  }
}
diagnostic <- do.call(rbind, diagnostic_rows)
rownames(diagnostic) <- NULL

expected_diagnostic_rows <- sum(vapply(
  loo_data, function(x) nrow(x$variant_checks), integer(1)
))
add_check(
  "diagnostic_retains_every_loo_final_variant", "leave_one_out_diagnostic",
  nrow(diagnostic) == expected_diagnostic_rows,
  nrow(diagnostic), expected_diagnostic_rows,
  "No smoothing sensitivity is filtered by its inferred ages."
)
d_extreme <- diagnostic[
  diagnostic$scenario == "D_papilionidae_only" &
    diagnostic$stage == "final_sensitivity_1e-6", , drop = FALSE
]
add_check(
  "d_smooth_1e_6_extreme_retained", "leave_one_out_diagnostic",
  nrow(d_extreme) == 1L &&
    abs(as.numeric(d_extreme$smoothing[[1L]]) - 1e-6) <= 1e-18 &&
    d_extreme$root_age_ma[[1L]] > 500 &&
    d_extreme$hesperiinae_age_ma[[1L]] > 300 &&
    isTRUE(d_extreme$all_hard_gates_pass[[1L]]),
  if (nrow(d_extreme) == 1L) {
    paste0("smooth=", d_extreme$smoothing[[1L]],
           ";root=", format_number(d_extreme$root_age_ma[[1L]]),
           ";hesp=", format_number(d_extreme$hesperiinae_age_ma[[1L]]),
           ";hard_gates=", d_extreme$all_hard_gates_pass[[1L]])
  } else "missing",
  "one unfiltered smooth=1e-6 row with root>500 Ma, Hesperiinae>300 Ma, hard gates PASS"
)

checks <- do.call(rbind, checks_list)
rownames(checks) <- NULL
critical_failures <- checks$status == "FAIL"
overall_status <- if (any(critical_failures)) "FAIL" else "PASS"

fmt_md <- function(x) {
  if (length(x) != 1L || is.na(x)) return("NA")
  value <- suppressWarnings(as.numeric(x))
  if (is.finite(value)) sprintf("%.12g", value) else as.character(x)
}

four_table <- c(
  "| Scenario | Role | Smooth | Root (Ma) | Papilionidae* (Ma) | Hesperiinae (Ma) | Hesperiinae shift from A (Ma) | QA |",
  "|---|---|---:|---:|---:|---:|---:|---:|",
  vapply(seq_len(nrow(four_summary)), function(i) {
    role <- c(
      A_two_node_main = "predeclared main analysis",
      B_three_node_sensitivity = "three-node sensitivity",
      C_root_only = "leave-one-out diagnostic",
      D_papilionidae_only = "leave-one-out diagnostic"
    )[[four_summary$scenario[[i]]]]
    paste0(
      "| ", four_summary$scenario[[i]], " | ", role,
      " | ", fmt_md(four_summary$smoothing[[i]]),
      " | ", fmt_md(four_summary$root_age_ma[[i]]),
      " | ", fmt_md(four_summary$papilionidae_age_ma[[i]]),
      " | ", fmt_md(four_summary$hesperiinae_age_ma[[i]]),
      " | ", fmt_md(four_summary$hesperiinae_shift_from_main_ma[[i]]),
      " | ", four_summary$qa_status[[i]], " |"
    )
  }, character(1))
)

d_extreme_line <- if (nrow(d_extreme) == 1L) {
  paste0(
    "The D sensitivity at `smooth=1e-6` is retained without filtering. It gives ",
    "root `", fmt_md(d_extreme$root_age_ma[[1L]]), " Ma`, Papilionidae* `",
    fmt_md(d_extreme$papilionidae_age_ma[[1L]]), " Ma`, and Hesperiinae `",
    fmt_md(d_extreme$hesperiinae_age_ma[[1L]]),
    " Ma`. Its technical tree/provenance gates pass, but the very large ",
    "unconstrained ages demonstrate smoothing instability and are diagnostic, ",
    "not a biologically preferred estimate."
  )
} else {
  "The expected D `smooth=1e-6` diagnostic row is missing."
}

report <- c(
  "# Four-scenario treePL calibration validation",
  "",
  paste0("**Overall validation status:** ", overall_status),
  "",
  "The analysis roles remain fixed by the study design: `A_two_node_main` is the manuscript's primary chronogram; `B_three_node_sensitivity` is the direct Hesperiinae-calibration sensitivity; C and D are leave-one-out diagnostics only. C or D must not replace A based on which numerical ages appear preferable.",
  "",
  "## Primary results",
  "",
  four_table,
  "",
  "Boundary labels are assigned only to nodes actually calibrated in a scenario. Unconstrained nodes are reported as `not_calibrated`, even when their inferred age lies outside the corresponding reference interval.",
  "",
  "## Leave-one-out smoothing diagnostics",
  "",
  paste0(
    "C final variants span root `",
    fmt_md(min(diagnostic$root_age_ma[diagnostic$scenario == "C_root_only"])),
    "-", fmt_md(max(diagnostic$root_age_ma[diagnostic$scenario == "C_root_only"])),
    " Ma` and Hesperiinae `",
    fmt_md(min(diagnostic$hesperiinae_age_ma[diagnostic$scenario == "C_root_only"])),
    "-", fmt_md(max(diagnostic$hesperiinae_age_ma[diagnostic$scenario == "C_root_only"])),
    " Ma`."
  ),
  paste0(
    "D final variants span root `",
    fmt_md(min(diagnostic$root_age_ma[diagnostic$scenario == "D_papilionidae_only"])),
    "-", fmt_md(max(diagnostic$root_age_ma[diagnostic$scenario == "D_papilionidae_only"])),
    " Ma` and Hesperiinae `",
    fmt_md(min(diagnostic$hesperiinae_age_ma[diagnostic$scenario == "D_papilionidae_only"])),
    "-", fmt_md(max(diagnostic$hesperiinae_age_ma[diagnostic$scenario == "D_papilionidae_only"])),
    " Ma`."
  ),
  "",
  d_extreme_line,
  "",
  "The full per-variant values and within-scenario three-node ranges are retained in `qa/leave_one_out_diagnostic_summary.tsv`.",
  "",
  "## Validation scope",
  "",
  paste0("- Fixed input SHA-256: `", EXPECTED_INPUT_SHA256, "`"),
  paste0("- Fixed patched treePL SHA-256 for C/D decisions: `", EXPECTED_TREEPL_SHA256, "`"),
  "- Every C/D primary, same-seed repeat, and declared smoothing sensitivity passed the stored 495-tip, tip-set, rooted topology, RF=0, finite/positive branch, ultrametric, calibration, execution-log, and provenance hard gates.",
  "- Same-seed primary/repeat output SHA-256 values are identical within C and within D.",
  "- The first two rows of the expanded summary are preserved cell-for-cell from the existing A/B summary; no A/B file was overwritten.",
  "",
  "## Critical failures",
  "",
  if (any(critical_failures)) {
    paste0("- `", checks$check_id[critical_failures], "`")
  } else "None."
)

if (any(critical_failures)) {
  stop(
    "Four-scenario validation hard gate failed before output creation: ",
    paste(checks$check_id[critical_failures], collapse = ", "),
    call. = FALSE
  )
}

write.table(
  four_summary, output_four, sep = "\t", quote = FALSE,
  row.names = FALSE, na = "NA"
)
write.table(
  diagnostic, output_diagnostic, sep = "\t", quote = FALSE,
  row.names = FALSE, na = "NA"
)
write.table(
  checks, output_checks, sep = "\t", quote = FALSE,
  row.names = FALSE, na = "NA"
)
writeLines(report, output_report, useBytes = TRUE)

cat(
  "Four-scenario summary PASS: ", nrow(checks), " hard checks; ",
  nrow(diagnostic), " leave-one-out final variants retained\n", sep = ""
)
cat("Expanded summary: ", output_four, "\n", sep = "")
cat("Diagnostic summary: ", output_diagnostic, "\n", sep = "")
cat("Report: ", output_report, "\n", sep = "")

