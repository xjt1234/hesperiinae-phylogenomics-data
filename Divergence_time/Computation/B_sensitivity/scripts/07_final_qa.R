#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ape))

options(stringsAsFactors = FALSE, digits = 15)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) {
  stop(
    paste(
      "Usage: 07_final_qa.R <input_ml_tree> <final_dated_tree> <run_root>",
      "<primary_final_config> <primary_stage_metadata>",
      "<cv_23point_decision_tsv> <patched_treepl_binary>"
    ),
    call. = FALSE
  )
}

input_path <- normalizePath(args[[1]], mustWork = TRUE)
final_path <- normalizePath(args[[2]], mustWork = TRUE)
run_root <- normalizePath(args[[3]], mustWork = TRUE)
config_path <- normalizePath(args[[4]], mustWork = TRUE)
metadata_path <- normalizePath(args[[5]], mustWork = TRUE)
cv_decision_path <- normalizePath(args[[6]], mustWork = TRUE)
binary_path <- normalizePath(args[[7]], mustWork = TRUE)

canonical_input_path <- file.path(
  run_root, "input", "verified_495_tip_input.treefile"
)
canonical_final_path <- file.path(
  run_root, "output", "T25_three_calibrations_dated.tre"
)
canonical_config_path <- file.path(run_root, "configs", "final.cfg")
canonical_metadata_path <- file.path(
  run_root, "logs", "final_primary_thorough.metadata.tsv"
)
canonical_cv_decision_path <- file.path(
  run_root, "qa", "cv_23point_decision.tsv"
)
source_19point_decision_path <- file.path(run_root, "qa", "cv_decision.tsv")
plateau_rule_path <- file.path(
  run_root,
  "qa",
  "cv_tie_plateau_rule_declared_before_low_extension_scores.md"
)
canonical_binary_path <- file.path(
  run_root, "software", "treePL_randomcv_fix_src", "treePL"
)
prime_recommendations_path <- file.path(
  run_root, "qa", "prime_recommendations.cfg"
)
EXPECTED_FINAL_SEED <- "2026090488"
EXPECTED_BINARY_SHA256 <- paste0(
  "5d6bfac8d18d835f24ebb79f0373b070a",
  "948b58a179f437f0aacb8fed74ccc5d"
)
EXPECTED_PATCH_SHA256 <- paste0(
  "d12cf27a8a6a7ce98583990e1f64f51c",
  "e11b26375f2121f6f005b90574a2e485"
)
EXPECTED_INPUT_SHA256 <- paste0(
  "4d35d553ed5bbf08022daef2361ffa0c",
  "342f8694a53aaf08d5350692094439e5"
)

qa_dir <- file.path(run_root, "qa")
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)

node_ages_path <- file.path(qa_dir, "node_ages.tsv")
checks_path <- file.path(qa_dir, "final_tree_checks.tsv")
report_path <- file.path(qa_dir, "qa_report_core.md")
sha256_path <- file.path(qa_dir, "final_provenance_sha256.tsv")
requested_outputs <- c(node_ages_path, checks_path, report_path, sha256_path)
existing_outputs <- requested_outputs[file.exists(requested_outputs)]
if (length(existing_outputs)) {
  stop(
    "Refusing to overwrite existing QA output(s): ",
    paste(existing_outputs, collapse = ", "),
    call. = FALSE
  )
}

STRICT_ULTRAMETRIC_TARGET_MA <- 1e-6
SERIALIZATION_ULTRAMETRIC_TOL_MA <- 1e-5
STRICT_BOUND_DIAGNOSTIC_TOL_MA <- 1e-6
SERIALIZATION_BOUND_TOL_MA <- 1e-5

fmt_num <- function(x) {
  if (length(x) != 1L || is.na(x) || !is.finite(x)) return("NA")
  sprintf("%.12g", x)
}


read_active_config <- function(path) {
  raw <- readLines(path, warn = FALSE)
  active <- trimws(sub("#.*$", "", raw))
  active <- active[nzchar(active)]
  has_equals <- grepl("=", active, fixed = TRUE)
  key <- active
  key[has_equals] <- trimws(sub("=.*$", "", active[has_equals]))
  value <- rep("", length(active))
  value[has_equals] <- trimws(sub("^[^=]*=", "", active[has_equals]))
  list(active = active, key = key, value = value)
}

values_for <- function(parsed, key) parsed$value[parsed$key == key]

read_metadata <- function(path) {
  raw <- readLines(path, warn = FALSE)
  pieces <- strsplit(raw, "\t", fixed = TRUE)
  valid <- lengths(pieces) >= 2L
  if (!length(raw) || !all(valid)) {
    return(data.frame(key = character(), value = character()))
  }
  data.frame(
    key = vapply(pieces, `[[`, character(1), 1L),
    value = vapply(
      pieces,
      function(x) paste(x[-1L], collapse = "\t"),
      character(1)
    ),
    stringsAsFactors = FALSE
  )
}

metadata_values <- function(metadata, key) metadata$value[metadata$key == key]

parse_bool <- function(x) {
  if (length(x) != 1L || is.na(x)) return(NA)
  normalized <- toupper(trimws(as.character(x)))
  if (normalized %in% c("TRUE", "T", "1")) return(TRUE)
  if (normalized %in% c("FALSE", "F", "0")) return(FALSE)
  NA
}

parse_numeric_csv <- function(x) {
  if (length(x) != 1L || is.na(x) || !nzchar(trimws(x))) return(numeric())
  tokens <- trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
  values <- suppressWarnings(as.numeric(tokens))
  if (any(!is.finite(values))) return(numeric())
  values
}

sha256_file <- function(path) {
  if (!file.exists(path) || dir.exists(path)) return(NA_character_)
  result <- suppressWarnings(system2(
    "/usr/bin/sha256sum",
    args = c("--", shQuote(path)),
    stdout = TRUE,
    stderr = TRUE
  ))
  status <- attr(result, "status")
  if ((!is.null(status) && status != 0L) || length(result) != 1L) {
    return(NA_character_)
  }
  token <- strsplit(result[[1]], "[[:space:]]+")[[1]][[1]]
  if (!grepl("^[0-9a-f]{64}$", token)) return(NA_character_)
  token
}

make_context <- function(tree, tree_name) {
  if (is.null(tree$edge.length)) {
    stop(tree_name, " has no branch lengths", call. = FALSE)
  }

  n_tip <- Ntip(tree)
  root <- setdiff(unique(tree$edge[, 1]), unique(tree$edge[, 2]))
  if (length(root) != 1L) {
    stop("Could not identify exactly one root in ", tree_name, call. = FALSE)
  }
  root <- as.integer(root[[1]])

  descendant_ids <- local({
    edge <- tree$edge
    nt <- n_tip
    walk <- function(node) {
      if (length(node) != 1L || is.na(node)) return(integer())
      node <- as.integer(node)
      if (node <= nt) return(node)
      children <- edge[edge[, 1] == node, 2]
      if (!length(children)) return(integer())
      as.integer(unlist(lapply(children, walk), use.names = FALSE))
    }
    walk
  })

  descendant_labels <- function(node) {
    ids <- descendant_ids(node)
    if (!length(ids)) return(character())
    sort(unique(tree$tip.label[ids]))
  }

  depth <- node.depth.edgelength(tree)
  root_children <- as.integer(tree$edge[tree$edge[, 1] == root, 2])

  list(
    tree = tree,
    n_tip = n_tip,
    root = root,
    depth = depth,
    tip_depth = depth[seq_len(n_tip)],
    root_children = root_children,
    descendant_ids = descendant_ids,
    descendant_labels = descendant_labels
  )
}

node_age_stats <- function(context, node) {
  if (length(node) != 1L || is.na(node)) {
    return(c(mean = NA_real_, min = NA_real_, max = NA_real_, range = NA_real_))
  }
  ids <- context$descendant_ids(node)
  if (!length(ids)) {
    return(c(mean = NA_real_, min = NA_real_, max = NA_real_, range = NA_real_))
  }
  ages <- context$depth[ids] - context$depth[node]
  c(
    mean = mean(ages),
    min = min(ages),
    max = max(ages),
    range = diff(range(ages))
  )
}

label_counts <- function(tree, labels) {
  vapply(labels, function(label) sum(tree$tip.label == label), integer(1))
}

safe_mrca <- function(tree, pair) {
  if (!all(label_counts(tree, pair) == 1L)) return(NA_integer_)
  as.integer(getMRCA(tree, pair))
}

safe_descendants <- function(context, node) {
  if (length(node) != 1L || is.na(node)) return(character())
  context$descendant_labels(node)
}

rooted_cluster_keys <- function(context) {
  internal_nodes <- sort(unique(context$tree$edge[, 1]))
  sort(vapply(
    internal_nodes,
    function(node) paste(context$descendant_labels(node), collapse = "\034"),
    character(1)
  ))
}

input_tree <- read.tree(input_path)
final_tree <- read.tree(final_path)
input <- make_context(input_tree, "input ML tree")
final <- make_context(final_tree, "final dated tree")

papilionidae_expected <- sort(c(
  "Parnassius_apollo_ncbi",
  "Parnassius_glacialis_ncbi",
  "Graphium_cloanthus_mydata",
  "Graphium_sarpedon_mydata",
  "Papilio_helenus_ncbi",
  "Papilio_machaon_ncbi",
  "Papilio_xuthus_ncbi"
))

node_specs <- data.frame(
  node_id = c(
    "CROWN_PAPILIONOIDEA",
    "CROWN_PAPILIONIDAE_NONBARONIA",
    "CROWN_HESPERIINAE"
  ),
  biological_definition = c(
    "crown Papilionoidea (total-tree root)",
    "crown Papilionidae excluding Baronia",
    "crown Hesperiinae"
  ),
  anchor_tip_1 = c(
    "Parnassius_apollo_ncbi",
    "Parnassius_apollo_ncbi",
    "Aeromachus_catocyanea_mydata"
  ),
  anchor_tip_2 = c(
    "Aeromachus_catocyanea_mydata",
    "Graphium_cloanthus_mydata",
    "Acada_biseriata_kawahara2023"
  ),
  expected_descendant_tip_count = c(495L, 7L, 488L),
  calibrated = c("yes", "yes", "yes"),
  interval_role = c(
    "treePL hard secondary calibration",
    "treePL hard secondary calibration",
    "treePL hard secondary calibration"
  ),
  lower_bound_ma = c(91.5046, 44.1968, 36.210861),
  upper_bound_ma = c(100.8925, 52.9473, 40.662537),
  stringsAsFactors = FALSE
)

required_labels <- unique(c(node_specs$anchor_tip_1, node_specs$anchor_tip_2))
input_nodes <- vapply(
  seq_len(nrow(node_specs)),
  function(i) safe_mrca(
    input_tree,
    c(node_specs$anchor_tip_1[[i]], node_specs$anchor_tip_2[[i]])
  ),
  integer(1)
)
final_nodes <- vapply(
  seq_len(nrow(node_specs)),
  function(i) safe_mrca(
    final_tree,
    c(node_specs$anchor_tip_1[[i]], node_specs$anchor_tip_2[[i]])
  ),
  integer(1)
)

input_desc <- lapply(input_nodes, function(node) safe_descendants(input, node))
final_desc <- lapply(final_nodes, function(node) safe_descendants(final, node))
input_desc_n <- lengths(input_desc)
final_desc_n <- lengths(final_desc)

final_age_matrix <- t(vapply(
  final_nodes,
  function(node) node_age_stats(final, node),
  numeric(4)
))
colnames(final_age_matrix) <- c("mean", "min", "max", "range")

node_age <- final_age_matrix[, "mean"]
distance_above_lower <- node_age - node_specs$lower_bound_ma
distance_below_upper <- node_specs$upper_bound_ma - node_age
nearest_boundary_distance <- pmin(
  abs(distance_above_lower),
  abs(distance_below_upper)
)
within_interval_no_tolerance <- (
  distance_above_lower >= 0 & distance_below_upper >= 0
)
within_interval_strict <- (
  distance_above_lower >= -STRICT_BOUND_DIAGNOSTIC_TOL_MA &
    distance_below_upper >= -STRICT_BOUND_DIAGNOSTIC_TOL_MA
)
within_interval_serialization <- (
  distance_above_lower >= -SERIALIZATION_BOUND_TOL_MA &
    distance_below_upper >= -SERIALIZATION_BOUND_TOL_MA
)
lower_boundary_hit_strict <- (
  abs(distance_above_lower) <= STRICT_BOUND_DIAGNOSTIC_TOL_MA
)
upper_boundary_hit_strict <- (
  abs(distance_below_upper) <= STRICT_BOUND_DIAGNOSTIC_TOL_MA
)
lower_boundary_near_serialization <- (
  abs(distance_above_lower) <= SERIALIZATION_BOUND_TOL_MA
)
upper_boundary_near_serialization <- (
  abs(distance_below_upper) <= SERIALIZATION_BOUND_TOL_MA
)

node_ages <- data.frame(
  node_id = node_specs$node_id,
  biological_definition = node_specs$biological_definition,
  anchor_tip_1 = node_specs$anchor_tip_1,
  anchor_tip_2 = node_specs$anchor_tip_2,
  input_mrca_node = input_nodes,
  final_mrca_node = final_nodes,
  final_is_root = final_nodes == final$root,
  input_descendant_tip_count = input_desc_n,
  final_descendant_tip_count = final_desc_n,
  calibrated = node_specs$calibrated,
  interval_role = node_specs$interval_role,
  final_node_age_ma = node_age,
  final_node_age_descendant_tip_min_ma = final_age_matrix[, "min"],
  final_node_age_descendant_tip_max_ma = final_age_matrix[, "max"],
  final_node_age_descendant_tip_range_ma = final_age_matrix[, "range"],
  lower_bound_ma = node_specs$lower_bound_ma,
  upper_bound_ma = node_specs$upper_bound_ma,
  within_interval = within_interval_serialization,
  within_interval_no_tolerance = within_interval_no_tolerance,
  within_interval_strict_1e_6_diagnostic = within_interval_strict,
  within_interval_serialization_tolerance_1e_5 =
    within_interval_serialization,
  signed_distance_above_lower_ma = distance_above_lower,
  signed_distance_below_upper_ma = distance_below_upper,
  nearest_boundary_distance_ma = nearest_boundary_distance,
  lower_boundary_hit_strict_1e_6 = lower_boundary_hit_strict,
  upper_boundary_hit_strict_1e_6 = upper_boundary_hit_strict,
  any_boundary_hit_strict_1e_6 =
    lower_boundary_hit_strict | upper_boundary_hit_strict,
  lower_boundary_near_serialization_1e_5 =
    lower_boundary_near_serialization,
  upper_boundary_near_serialization_1e_5 =
    upper_boundary_near_serialization,
  strict_boundary_diagnostic_tolerance_ma = STRICT_BOUND_DIAGNOSTIC_TOL_MA,
  age_definition = paste0(
    "mean node-to-descendant-tip path length; descendant-tip range retained ",
    "as serialization/ultrametric diagnostic"
  ),
  stringsAsFactors = FALSE
)

checks_list <- list()
add_check <- function(
    check_id,
    scope,
    passed,
    observed,
    expected,
    critical = TRUE,
    note = "") {
  passed <- isTRUE(passed)
  checks_list[[length(checks_list) + 1L]] <<- data.frame(
    check_id = check_id,
    scope = scope,
    critical = critical,
    status = if (passed) "PASS" else if (critical) "FAIL" else "WARN",
    observed = as.character(observed),
    expected = as.character(expected),
    note = as.character(note),
    stringsAsFactors = FALSE
  )
}

config <- read_active_config(config_path)
metadata <- read_metadata(metadata_path)
cv_decision <- tryCatch(
  read.delim(cv_decision_path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("NA", "")),
  error = function(e) NULL
)
source_19point_decision <- tryCatch(
  read.delim(source_19point_decision_path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("NA", "")),
  error = function(e) NULL
)
cv_required_columns <- c(
  "status", "decision_basis", "selected_smoothing_rule", "selected_smoothing",
  "aggregate_winner_at_grid_boundary", "boundary_requires_extension",
  "aggregate_winner_count", "aggregate_winner_values",
  "lower_boundary_in_aggregate_winner_set",
  "lower_boundary_strictly_below_all_interior",
  "lowest_three_grid_values_all_aggregate_winners",
  "aggregate_winner_set_contains_interior", "finite_precision_plateau",
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
  "diagnostic_19point_lower_boundary_in_aggregate_winner_set",
  "diagnostic_19point_aggregate_winner_values",
  "tie_relative_tolerance"
)
source_19point_required_columns <- c(
  "status", "decision_basis", "selected_smoothing",
  "aggregate_winner_at_grid_boundary", "formal_runs",
  "formal_grid_points_per_run", "formal_rows",
  "formal_seed_values", "formal_exit_codes_all_zero",
  "formal_stderr_all_empty", "formal_dated_trees_all_nonempty"
)
cv_decision_readable <- (
  !is.null(cv_decision) &&
    nrow(cv_decision) == 1L &&
    all(cv_required_columns %in% names(cv_decision))
)
source_19point_decision_readable <- (
  !is.null(source_19point_decision) &&
    nrow(source_19point_decision) == 1L &&
    all(source_19point_required_columns %in% names(source_19point_decision))
)
decision_value <- function(table, readable, name) {
  if (!readable || !name %in% names(table)) return(character())
  as.character(table[[name]][[1L]])
}
cv_value <- function(name) decision_value(cv_decision, cv_decision_readable, name)
source_19point_value <- function(name) decision_value(source_19point_decision, source_19point_decision_readable, name)
selected_smoothing <- suppressWarnings(as.numeric(cv_value("selected_smoothing")))
config_smoothing <- suppressWarnings(as.numeric(values_for(config, "smooth")))
aggregate_winner_values <- parse_numeric_csv(cv_value("aggregate_winner_values"))
sensitivity_smoothing_values <- parse_numeric_csv(cv_value("sensitivity_smoothing_values"))
diagnostic_19point_aggregate_winner_values <- parse_numeric_csv(
  cv_value("diagnostic_19point_aggregate_winner_values")
)
aggregate_winner_count <- suppressWarnings(as.integer(cv_value("aggregate_winner_count")))
cv_status <- cv_value("status")
cv_selected_rule <- cv_value("selected_smoothing_rule")
cv_aggregate_boundary <- parse_bool(cv_value("aggregate_winner_at_grid_boundary"))
cv_boundary_requires_extension <- parse_bool(cv_value("boundary_requires_extension"))
cv_finite_precision_plateau <- parse_bool(cv_value("finite_precision_plateau"))
cv_sensitivity_required <- parse_bool(cv_value("sensitivity_required"))
cv_diagnostic_19point_lower_boundary_in_aggregate_winner_set <- parse_bool(
  cv_value("diagnostic_19point_lower_boundary_in_aggregate_winner_set")
)
cv_tie_tolerance <- suppressWarnings(as.numeric(cv_value("tie_relative_tolerance")))

unique_nonboundary_selection_contract <- (
  identical(cv_finite_precision_plateau, FALSE) &&
    identical(cv_boundary_requires_extension, FALSE) &&
    identical(cv_aggregate_boundary, FALSE) &&
    length(aggregate_winner_count) == 1L &&
    aggregate_winner_count == 1L &&
    length(aggregate_winner_values) == 1L &&
    length(selected_smoothing) == 1L &&
    is.finite(selected_smoothing) &&
    identical(selected_smoothing, aggregate_winner_values) &&
    identical(cv_selected_rule, "unique_nonboundary_aggregate_minimum")
)

finite_precision_plateau_selection_contract <- (
  identical(cv_finite_precision_plateau, TRUE) &&
    identical(cv_boundary_requires_extension, FALSE) &&
    identical(cv_aggregate_boundary, TRUE) &&
    length(aggregate_winner_count) == 1L &&
    aggregate_winner_count >= 3L &&
    length(aggregate_winner_values) == aggregate_winner_count &&
    anyDuplicated(aggregate_winner_values) == 0L &&
    length(selected_smoothing) == 1L &&
    is.finite(selected_smoothing) &&
    identical(selected_smoothing, max(aggregate_winner_values)) &&
    identical(cv_selected_rule, "largest_smoothing_in_aggregate_tied_minimum_plateau") &&
    identical(parse_bool(cv_value("lower_boundary_in_aggregate_winner_set")), TRUE) &&
    identical(parse_bool(cv_value("lower_boundary_strictly_below_all_interior")), FALSE) &&
    identical(parse_bool(cv_value("lowest_three_grid_values_all_aggregate_winners")), TRUE) &&
    identical(parse_bool(cv_value("aggregate_winner_set_contains_interior")), TRUE) &&
    all(c(1e-15, 1e-16, 1e-17) %in% aggregate_winner_values)
)
cv_selection_contract_resolved <- (
  unique_nonboundary_selection_contract ||
    finite_precision_plateau_selection_contract
)
cv_sensitivity_status_contract <- (
  (identical(cv_status, "ACCEPTED") &&
     identical(cv_sensitivity_required, FALSE)) ||
    (identical(cv_status, "ACCEPTED_WITH_SENSITIVITIES") &&
       identical(cv_sensitivity_required, TRUE))
)
plateau_sensitivity_values_contract <- if (!isTRUE(cv_finite_precision_plateau)) {
  TRUE
} else {
  length(selected_smoothing) == 1L &&
    is.finite(selected_smoothing) &&
    length(aggregate_winner_values) >= 3L &&
    length(sensitivity_smoothing_values) > 0L &&
    !selected_smoothing %in% sensitivity_smoothing_values &&
    all(setdiff(aggregate_winner_values, selected_smoothing) %in% sensitivity_smoothing_values)
}

prime_recommendations <- if (file.exists(prime_recommendations_path)) {
  read_active_config(prime_recommendations_path)
} else {
  list(active = character(), key = character(), value = character())
}
expected_prime_lines <- c(
  "opt = 4",
  "moredetail",
  "optad = 4",
  "moredetailad",
  "optcvad = 1",
  "moredetailcvad"
)
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

metadata_stage <- metadata_values(metadata, "stage")
stage_stdout_path <- file.path(run_root, "logs", "final_primary_thorough.stdout")
stage_stderr_path <- file.path(run_root, "logs", "final_primary_thorough.stderr")
stage_resource_path <- file.path(run_root, "logs", "final_primary_thorough.resources.txt")
mapping_path <- file.path(run_root, "mapping", "calibration_node_mapping.tsv")
patch_path <- file.path(run_root, "software", "randomcv_parent_index_fix.patch")

provenance <- data.frame(
  role = c(
    "input_ml_tree",
    "final_dated_tree",
    "primary_final_config",
    "primary_stage_metadata",
    "cv_23point_decision",
    "source_19point_cv_decision",
    "predeclared_cv_tie_plateau_rule",
    "prime_recommendations",
    "calibration_mapping",
    "patched_treepl_binary",
    "randomcv_parent_index_patch",
    "primary_stage_stdout",
    "primary_stage_stderr",
    "primary_stage_resources"
  ),
  path = c(
    input_path,
    final_path,
    config_path,
    metadata_path,
    cv_decision_path,
    source_19point_decision_path,
    plateau_rule_path,
    prime_recommendations_path,
    mapping_path,
    binary_path,
    patch_path,
    stage_stdout_path,
    stage_stderr_path,
    stage_resource_path
  ),
  stringsAsFactors = FALSE
)
provenance$file_exists <- file.exists(provenance$path)
provenance$file_size_bytes <- vapply(
  provenance$path,
  function(path) {
    if (!file.exists(path) || dir.exists(path)) return(NA_real_)
    as.numeric(file.info(path)$size)
  },
  numeric(1)
)
provenance$sha256 <- vapply(provenance$path, sha256_file, character(1))
provenance$hash_timing <- paste0(
  "computed independently by 07_final_qa.R after the final stage; ",
  "not an execution-time hash embedded in wrapper metadata"
)
provenance_sha256 <- setNames(provenance$sha256, provenance$role)

input_root_sizes <- sort(vapply(
  input$root_children,
  function(node) length(input$descendant_ids(node)),
  integer(1)
))
final_root_sizes <- sort(vapply(
  final$root_children,
  function(node) length(final$descendant_ids(node)),
  integer(1)
))

same_tip_set <- (
  length(input_tree$tip.label) == length(final_tree$tip.label) &&
    setequal(input_tree$tip.label, final_tree$tip.label)
)
rf_ph85 <- if (same_tip_set) {
  suppressWarnings(as.numeric(dist.topo(input_tree, final_tree, method = "PH85")))
} else {
  NA_real_
}
rooted_clusters_equal <- if (same_tip_set) {
  identical(rooted_cluster_keys(input), rooted_cluster_keys(final))
} else {
  FALSE
}
all_equal_rooted_topology <- if (same_tip_set) {
  isTRUE(all.equal.phylo(
    input_tree,
    final_tree,
    use.edge.length = FALSE,
    use.tip.label = TRUE
  ))
} else {
  FALSE
}

input_tip_range <- diff(range(input$tip_depth))
final_tip_range <- diff(range(final$tip_depth))
ape_final_ultrametric <- is.ultrametric(final_tree)

add_check("canonical_input_path", "provenance", identical(input_path, canonical_input_path), input_path, canonical_input_path)
add_check("canonical_primary_output_path", "provenance", identical(final_path, canonical_final_path), final_path, canonical_final_path)
add_check("canonical_primary_config_path", "provenance", identical(config_path, canonical_config_path), config_path, canonical_config_path)
add_check("canonical_primary_metadata_path", "provenance", identical(metadata_path, canonical_metadata_path), metadata_path, canonical_metadata_path)
add_check("canonical_cv_decision_path", "provenance", identical(cv_decision_path, canonical_cv_decision_path), cv_decision_path, canonical_cv_decision_path)
add_check("source_19point_decision_exists_nonzero", "provenance", file.exists(source_19point_decision_path) && file.info(source_19point_decision_path)$size > 0, if (file.exists(source_19point_decision_path)) file.info(source_19point_decision_path)$size else NA, ">0 bytes")
add_check("predeclared_plateau_rule_exists_nonzero", "provenance", file.exists(plateau_rule_path) && file.info(plateau_rule_path)$size > 0, if (file.exists(plateau_rule_path)) file.info(plateau_rule_path)$size else NA, ">0 bytes")
add_check("canonical_patched_binary_path", "provenance", identical(binary_path, canonical_binary_path), binary_path, canonical_binary_path)
add_check("prime_recommendations_file_exists", "provenance", file.exists(prime_recommendations_path), file.exists(prime_recommendations_path), TRUE)
add_check("calibration_mapping_file_exists", "provenance", file.exists(mapping_path), file.exists(mapping_path), TRUE)
add_check("primary_stage_stdout_exists", "provenance", file.exists(stage_stdout_path), file.exists(stage_stdout_path), TRUE)
add_check("primary_stage_stderr_exists", "provenance", file.exists(stage_stderr_path), file.exists(stage_stderr_path), TRUE)
add_check("primary_stage_resources_exists", "provenance", file.exists(stage_resource_path), file.exists(stage_resource_path), TRUE)
add_check("all_provenance_sha256_computed", "provenance", all(!is.na(provenance$sha256)), sum(is.na(provenance$sha256)), 0L)
add_check("input_sha256_matches_predeclared", "provenance", identical(unname(provenance_sha256[["input_ml_tree"]]), EXPECTED_INPUT_SHA256), provenance_sha256[["input_ml_tree"]], EXPECTED_INPUT_SHA256)
add_check("patched_binary_sha256_matches_predeclared", "provenance", identical(unname(provenance_sha256[["patched_treepl_binary"]]), EXPECTED_BINARY_SHA256), provenance_sha256[["patched_treepl_binary"]], EXPECTED_BINARY_SHA256)
add_check("randomcv_patch_sha256_matches_predeclared", "provenance", identical(unname(provenance_sha256[["randomcv_parent_index_patch"]]), EXPECTED_PATCH_SHA256), provenance_sha256[["randomcv_parent_index_patch"]], EXPECTED_PATCH_SHA256)

metadata_required_keys <- c("stage", "pid", "start_time", "binary", "config", "end_time", "exit_code", "elapsed_s")
metadata_required_counts <- vapply(metadata_required_keys, function(key) sum(metadata$key == key), integer(1))
add_check("metadata_required_keys_exactly_once", "stage_metadata", all(metadata_required_counts == 1L), paste(paste(metadata_required_keys, metadata_required_counts, sep = "="), collapse = ";"), paste(paste(metadata_required_keys, 1L, sep = "="), collapse = ";"))
add_check("metadata_stage_is_primary_thorough", "stage_metadata", identical(metadata_stage, "final_primary_thorough"), paste(metadata_stage, collapse = ";"), "final_primary_thorough")
add_check("metadata_config_matches_primary_config", "stage_metadata", identical(metadata_values(metadata, "config"), config_path), paste(metadata_values(metadata, "config"), collapse = ";"), config_path)
add_check("metadata_binary_matches_patched_binary", "stage_metadata", identical(metadata_values(metadata, "binary"), binary_path), paste(metadata_values(metadata, "binary"), collapse = ";"), binary_path)
add_check("metadata_exit_code_zero", "stage_metadata", identical(metadata_values(metadata, "exit_code"), "0"), paste(metadata_values(metadata, "exit_code"), collapse = ";"), "0")
add_check("metadata_has_end_time", "stage_metadata", length(metadata_values(metadata, "end_time")) == 1L && nzchar(metadata_values(metadata, "end_time")), paste(metadata_values(metadata, "end_time"), collapse = ";"), "one nonempty end_time")
add_check("primary_stage_stderr_empty", "stage_metadata", file.exists(stage_stderr_path) && file.info(stage_stderr_path)$size == 0, if (file.exists(stage_stderr_path)) file.info(stage_stderr_path)$size else NA, "0 bytes")

config_calibration_lines <- config$active[config$key %in% c("mrca", "min", "max")]
config_optimization_lines <- config$active[config$key %in% c("opt", "moredetail", "optad", "moredetailad", "optcvad", "moredetailcvad")]
prohibited_final_keys <- c("prime", "cv", "randomcv", "cvstart", "cvstop", "cvmultstep", "cviter", "cvoutfile", "log_pen")
metadata_elapsed <- suppressWarnings(as.numeric(metadata_values(metadata, "elapsed_s")))
add_check("metadata_elapsed_nonnegative", "stage_metadata", length(metadata_elapsed) == 1L && is.finite(metadata_elapsed) && metadata_elapsed >= 0, paste(metadata_values(metadata, "elapsed_s"), collapse = ";"), "one finite value >=0 seconds")
add_check("config_treefile_matches_input", "primary_config", identical(values_for(config, "treefile"), input_path), paste(values_for(config, "treefile"), collapse = ";"), input_path)
add_check("config_outfile_matches_primary_tree", "primary_config", identical(values_for(config, "outfile"), final_path), paste(values_for(config, "outfile"), collapse = ";"), final_path)
add_check("config_numsites_182682", "primary_config", identical(values_for(config, "numsites"), "182682"), paste(values_for(config, "numsites"), collapse = ";"), "182682")
add_check("config_exactly_three_mrca_min_max", "primary_config", sum(config$key == "mrca") == 3L && sum(config$key == "min") == 3L && sum(config$key == "max") == 3L, paste0("mrca=", sum(config$key == "mrca"), ";min=", sum(config$key == "min"), ";max=", sum(config$key == "max")), "mrca=3;min=3;max=3")
add_check("config_three_calibrations_exact", "primary_config", identical(config_calibration_lines, expected_calibration_lines), paste(config_calibration_lines, collapse = " | "), paste(expected_calibration_lines, collapse = " | "))
add_check("config_final_seed_exact", "primary_config", identical(values_for(config, "seed"), EXPECTED_FINAL_SEED), paste(values_for(config, "seed"), collapse = ";"), EXPECTED_FINAL_SEED)
add_check("config_thorough_exactly_once", "primary_config", sum(config$active == "thorough") == 1L, sum(config$active == "thorough"), 1L)
add_check("config_prime_cv_and_logpen_absent", "primary_config", !any(config$key %in% prohibited_final_keys), paste(config$key[config$key %in% prohibited_final_keys], collapse = ";"), "none")
add_check("config_one_positive_smoothing", "primary_config", length(config_smoothing) == 1L && is.finite(config_smoothing) && config_smoothing > 0, paste(values_for(config, "smooth"), collapse = ";"), "one positive value")
source_19point_provenance_flags <- c("formal_exit_codes_all_zero", "formal_stderr_all_empty", "formal_dated_trees_all_nonempty")
add_check("source_19point_decision_readable_single_row", "source_19point_cv_decision", source_19point_decision_readable, if (is.null(source_19point_decision)) "unreadable" else paste0("rows=", nrow(source_19point_decision), ";columns=", ncol(source_19point_decision)), paste0("one row with columns: ", paste(source_19point_required_columns, collapse = ",")))
add_check("source_19point_requires_grid_extension", "source_19point_cv_decision", identical(source_19point_value("status"), "REQUIRES_GRID_EXTENSION"), paste(source_19point_value("status"), collapse = ";"), "REQUIRES_GRID_EXTENSION")
add_check("source_19point_decision_basis_exact", "source_19point_cv_decision", identical(source_19point_value("decision_basis"), "lowest_median_raw_chisq_across_three_19_point_randomCV_runs"), paste(source_19point_value("decision_basis"), collapse = ";"), "lowest_median_raw_chisq_across_three_19_point_randomCV_runs")
add_check("source_19point_selected_smoothing_missing", "source_19point_cv_decision", length(source_19point_value("selected_smoothing")) == 1L && is.na(source_19point_value("selected_smoothing")), paste(source_19point_value("selected_smoothing"), collapse = ";"), "NA")
add_check("source_19point_aggregate_winner_at_boundary", "source_19point_cv_decision", identical(parse_bool(source_19point_value("aggregate_winner_at_grid_boundary")), TRUE), paste(source_19point_value("aggregate_winner_at_grid_boundary"), collapse = ";"), TRUE)
add_check("source_19point_formal_runs_exact", "source_19point_cv_decision", identical(source_19point_value("formal_runs"), "cv_finalgrid_run_01,cv_finalgrid_run_02,cv_finalgrid_run_03"), paste(source_19point_value("formal_runs"), collapse = ";"), "cv_finalgrid_run_01,cv_finalgrid_run_02,cv_finalgrid_run_03")
add_check("source_19point_grid_19_points_57_rows", "source_19point_cv_decision", identical(source_19point_value("formal_grid_points_per_run"), "19") && identical(source_19point_value("formal_rows"), "57"), paste0("points=", paste(source_19point_value("formal_grid_points_per_run"), collapse = ";"), ";rows=", paste(source_19point_value("formal_rows"), collapse = ";")), "points=19;rows=57")
add_check("source_19point_formal_seeds_exact", "source_19point_cv_decision", identical(source_19point_value("formal_seed_values"), "2026090441,2026090442,2026090443"), paste(source_19point_value("formal_seed_values"), collapse = ";"), "2026090441,2026090442,2026090443")
add_check("source_19point_formal_provenance_flags_true", "source_19point_cv_decision", all(vapply(source_19point_provenance_flags, function(key) isTRUE(parse_bool(source_19point_value(key))), logical(1))), paste(paste(source_19point_provenance_flags, vapply(source_19point_provenance_flags, function(key) paste(source_19point_value(key), collapse = ";"), character(1)), sep = "="), collapse = ";"), paste(paste(source_19point_provenance_flags, TRUE, sep = "="), collapse = ";"))

cv_formal_validation_flags <- c(
  "formal_exit_codes_all_zero", "formal_stderr_all_empty",
  "formal_dated_trees_all_nonempty", "formal_dated_trees_valid_branch_length_newick",
  "formal_dated_tree_tip_sets_match_input", "formal_dated_tree_topologies_match_input",
  "formal_dated_trees_ultrametric_within_tolerance",
  "formal_config_contracts_validated", "formal_config_lints_all_pass",
  "formal_metadata_validated", "formal_stdout_score_sets_validated",
  "formal_resource_logs_validated"
)
add_check("cv_23point_decision_readable_single_row", "cv_23point_decision", cv_decision_readable, if (is.null(cv_decision)) "unreadable" else paste0("rows=", nrow(cv_decision), ";columns=", ncol(cv_decision)), paste0("one row with columns: ", paste(cv_required_columns, collapse = ",")))
add_check("cv_23point_status_accepted", "cv_23point_decision", cv_status %in% c("ACCEPTED", "ACCEPTED_WITH_SENSITIVITIES"), paste(cv_status, collapse = ";"), "ACCEPTED or ACCEPTED_WITH_SENSITIVITIES")
add_check("cv_23point_decision_basis_exact", "cv_23point_decision", identical(cv_value("decision_basis"), "lowest_median_raw_chisq_across_three_23_point_randomCV_runs"), paste(cv_value("decision_basis"), collapse = ";"), "lowest_median_raw_chisq_across_three_23_point_randomCV_runs")
add_check("cv_diagnostic_19point_lower_boundary_in_aggregate_winner_set", "cv_23point_decision", identical(cv_diagnostic_19point_lower_boundary_in_aggregate_winner_set, TRUE), paste(cv_value("diagnostic_19point_lower_boundary_in_aggregate_winner_set"), collapse = ";"), TRUE)
add_check("cv_diagnostic_19point_aggregate_winners_include_1e_13", "cv_23point_decision", 1e-13 %in% diagnostic_19point_aggregate_winner_values, paste(cv_value("diagnostic_19point_aggregate_winner_values"), collapse = ";"), "contains 0.0000000000001")
add_check("cv_23point_selected_smoothing_positive", "cv_23point_decision", length(selected_smoothing) == 1L && is.finite(selected_smoothing) && selected_smoothing > 0, paste(cv_value("selected_smoothing"), collapse = ";"), "one positive selected smoothing")
add_check("config_smoothing_equals_cv_23point_selected", "config_vs_cv", length(config_smoothing) == 1L && length(selected_smoothing) == 1L && is.finite(config_smoothing) && is.finite(selected_smoothing) && identical(config_smoothing, selected_smoothing), paste(values_for(config, "smooth"), cv_value("selected_smoothing"), sep = " vs ", collapse = ";"), "numeric equality")
add_check("cv_23point_boundary_fully_resolved", "cv_23point_decision", identical(cv_boundary_requires_extension, FALSE), paste(cv_value("boundary_requires_extension"), collapse = ";"), FALSE)
add_check("cv_23point_selection_contract_resolved", "cv_23point_decision", cv_selection_contract_resolved, paste0("unique_nonboundary=", unique_nonboundary_selection_contract, ";finite_precision_plateau=", finite_precision_plateau_selection_contract, ";rule=", paste(cv_selected_rule, collapse = ";")), "unique nonboundary winner OR valid predeclared finite-precision plateau")
add_check("cv_plateau_selected_is_largest_tied_smoothing", "cv_23point_decision", !isTRUE(cv_finite_precision_plateau) || (length(selected_smoothing) == 1L && length(aggregate_winner_values) >= 3L && identical(selected_smoothing, max(aggregate_winner_values))), paste0("selected=", paste(cv_value("selected_smoothing"), collapse = ";"), ";aggregate_winners=", paste(cv_value("aggregate_winner_values"), collapse = ";")), "for plateau, selected=max(aggregate_winner_values)")
add_check("cv_aggregate_winner_count_matches_values", "cv_23point_decision", length(aggregate_winner_count) == 1L && aggregate_winner_count > 0L && length(aggregate_winner_values) == aggregate_winner_count && anyDuplicated(aggregate_winner_values) == 0L, paste0("declared_count=", paste(cv_value("aggregate_winner_count"), collapse = ";"), ";parsed_count=", length(aggregate_winner_values)), "positive equal counts with no duplicates")
add_check("cv_sensitivity_status_consistent", "cv_23point_decision", cv_sensitivity_status_contract, paste0("status=", paste(cv_status, collapse = ";"), ";sensitivity_required=", paste(cv_value("sensitivity_required"), collapse = ";")), "ACCEPTED iff FALSE; ACCEPTED_WITH_SENSITIVITIES iff TRUE")
add_check("cv_plateau_sensitivity_values_declared", "cv_23point_decision", plateau_sensitivity_values_contract, paste(cv_value("sensitivity_smoothing_values"), collapse = ";"), "for plateau, all nonprimary aggregate tied winners included and primary excluded")
add_check("cv_tie_relative_tolerance_exact", "cv_23point_decision", length(cv_tie_tolerance) == 1L && is.finite(cv_tie_tolerance) && identical(cv_tie_tolerance, 1e-10), paste(cv_value("tie_relative_tolerance"), collapse = ";"), "1e-10")
add_check("cv_23point_formal_runs_exact", "cv_23point_decision", identical(cv_value("formal_runs"), "cv_extendedgrid_run_01,cv_extendedgrid_run_02,cv_extendedgrid_run_03"), paste(cv_value("formal_runs"), collapse = ";"), "cv_extendedgrid_run_01,cv_extendedgrid_run_02,cv_extendedgrid_run_03")
add_check("cv_23point_grid_23_points_69_rows", "cv_23point_decision", identical(cv_value("formal_grid_points_per_run"), "23") && identical(cv_value("formal_rows"), "69"), paste0("points=", paste(cv_value("formal_grid_points_per_run"), collapse = ";"), ";rows=", paste(cv_value("formal_rows"), collapse = ";")), "points=23;rows=69")
add_check("cv_23point_formal_seeds_exact", "cv_23point_decision", identical(cv_value("formal_seed_values"), "2026090441,2026090442,2026090443"), paste(cv_value("formal_seed_values"), collapse = ";"), "2026090441,2026090442,2026090443")
add_check("cv_23point_formal_validation_flags_true", "cv_23point_decision", all(vapply(cv_formal_validation_flags, function(key) isTRUE(parse_bool(cv_value(key))), logical(1))), paste(paste(cv_formal_validation_flags, vapply(cv_formal_validation_flags, function(key) paste(cv_value(key), collapse = ";"), character(1)), sep = "="), collapse = ";"), paste(paste(cv_formal_validation_flags, TRUE, sep = "="), collapse = ";"))
add_check("prime_recommendations_exact", "prime", identical(prime_recommendations$active, expected_prime_lines), paste(prime_recommendations$active, collapse = " | "), paste(expected_prime_lines, collapse = " | "))
add_check("config_uses_exact_prime_recommendations", "primary_config", identical(config_optimization_lines, expected_prime_lines), paste(config_optimization_lines, collapse = " | "), paste(expected_prime_lines, collapse = " | "))

add_check("input_tip_count_495", "input", input$n_tip == 495L, input$n_tip, 495L)
add_check("final_file_exists_nonzero", "final", file.info(final_path)$size > 0, file.info(final_path)$size, ">0 bytes")
add_check("final_parent_not_younger_than_child", "final", all(final_tree$edge.length >= 0), sum(final_tree$edge.length < 0), 0L, TRUE, "Nonnegative temporal edge duration is equivalent to parent age not being younger than child age.")
add_check("final_tip_count_495", "final", final$n_tip == 495L, final$n_tip, 495L)
add_check(
  "same_tip_set",
  "input_vs_final",
  same_tip_set,
  sprintf(
    "input_only=%d;final_only=%d",
    length(setdiff(input_tree$tip.label, final_tree$tip.label)),
    length(setdiff(final_tree$tip.label, input_tree$tip.label))
  ),
  "input_only=0;final_only=0"
)
add_check(
  "input_unique_tip_labels",
  "input",
  anyDuplicated(input_tree$tip.label) == 0L,
  anyDuplicated(input_tree$tip.label),
  0L
)
add_check(
  "final_unique_tip_labels",
  "final",
  anyDuplicated(final_tree$tip.label) == 0L,
  anyDuplicated(final_tree$tip.label),
  0L
)
add_check("input_rooted", "input", is.rooted(input_tree), is.rooted(input_tree), TRUE)
add_check("final_rooted", "final", is.rooted(final_tree), is.rooted(final_tree), TRUE)
add_check(
  "input_binary",
  "input",
  is.binary.phylo(input_tree),
  is.binary.phylo(input_tree),
  TRUE
)
add_check(
  "final_binary",
  "final",
  is.binary.phylo(final_tree),
  is.binary.phylo(final_tree),
  TRUE
)
add_check(
  "input_branch_lengths_finite",
  "input",
  all(is.finite(input_tree$edge.length)),
  sum(!is.finite(input_tree$edge.length)),
  0L
)
add_check(
  "final_branch_lengths_finite",
  "final",
  all(is.finite(final_tree$edge.length)),
  sum(!is.finite(final_tree$edge.length)),
  0L
)
add_check(
  "input_branch_lengths_nonnegative",
  "input",
  all(input_tree$edge.length >= 0),
  sum(input_tree$edge.length < 0),
  0L
)
add_check(
  "final_branch_lengths_nonnegative",
  "final",
  all(final_tree$edge.length >= 0),
  sum(final_tree$edge.length < 0),
  0L
)
add_check(
  "input_expected_non_ultrametric",
  "input",
  !is.ultrametric(input_tree),
  paste0("root_to_tip_range=", fmt_num(input_tip_range)),
  "non-ultrametric ML phylogram"
)
add_check(
  "final_ultrametric_ape_default",
  "final",
  ape_final_ultrametric,
  ape_final_ultrametric,
  TRUE,
  FALSE,
  "treePL serializes time branches to six decimal places, so ape default can warn despite an acceptable absolute range."
)
add_check(
  "final_ultrametric_absolute_tol_1e-6_ma",
  "final",
  final_tip_range <= STRICT_ULTRAMETRIC_TARGET_MA,
  fmt_num(final_tip_range),
  paste0("<=", fmt_num(STRICT_ULTRAMETRIC_TARGET_MA), " Ma"),
  FALSE,
  "Strict target reported as a warning because the output Newick uses six-decimal time branches."
)
add_check(
  "final_ultrametric_serialization_tol_1e-5_ma",
  "final",
  final_tip_range <= SERIALIZATION_ULTRAMETRIC_TOL_MA,
  fmt_num(final_tip_range),
  paste0("<=", fmt_num(SERIALIZATION_ULTRAMETRIC_TOL_MA), " Ma"),
  TRUE,
  "Critical absolute tolerance chosen before final QA to accommodate six-decimal serialization rounding."
)
add_check(
  "input_root_split_7_488",
  "input",
  identical(input_root_sizes, c(7L, 488L)),
  paste(input_root_sizes, collapse = ","),
  "7,488"
)
add_check(
  "final_root_split_7_488",
  "final",
  identical(final_root_sizes, c(7L, 488L)),
  paste(final_root_sizes, collapse = ","),
  "7,488"
)
add_check("ape_PH85_RF_zero", "input_vs_final", identical(rf_ph85, 0), rf_ph85, 0)
add_check(
  "rooted_cluster_sets_identical",
  "input_vs_final",
  rooted_clusters_equal,
  rooted_clusters_equal,
  TRUE,
  TRUE,
  "Root-aware clade-set comparison complements ape PH85 RF."
)
add_check(
  "ape_all_equal_topology_no_lengths",
  "input_vs_final",
  all_equal_rooted_topology,
  all_equal_rooted_topology,
  TRUE,
  FALSE,
  "Supplementary ape serialization check; RF=0 and exact rooted-cluster equality remain critical."
)
add_check(
  "required_anchor_labels_exactly_once_input",
  "input",
  all(label_counts(input_tree, required_labels) == 1L),
  paste(label_counts(input_tree, required_labels), collapse = ","),
  paste(rep(1L, length(required_labels)), collapse = ",")
)
add_check(
  "required_anchor_labels_exactly_once_final",
  "final",
  all(label_counts(final_tree, required_labels) == 1L),
  paste(label_counts(final_tree, required_labels), collapse = ","),
  paste(rep(1L, length(required_labels)), collapse = ",")
)
add_check(
  "obsolete_Graphium_ncbi_absent_final",
  "final",
  sum(final_tree$tip.label == "Graphium_cloanthus_ncbi") == 0L,
  sum(final_tree$tip.label == "Graphium_cloanthus_ncbi"),
  0L
)
add_check(
  "Baronia_absent_final",
  "final",
  !any(grepl("^Baronia_", final_tree$tip.label)),
  sum(grepl("^Baronia_", final_tree$tip.label)),
  0L
)

add_check(
  "Papilionoidea_anchor_is_input_root",
  "input",
  !is.na(input_nodes[[1]]) && input_nodes[[1]] == input$root,
  input_nodes[[1]],
  input$root
)
add_check(
  "Papilionoidea_anchor_is_final_root",
  "final",
  !is.na(final_nodes[[1]]) && final_nodes[[1]] == final$root,
  final_nodes[[1]],
  final$root
)
add_check(
  "Papilionoidea_descendants_495_input",
  "input",
  input_desc_n[[1]] == 495L,
  input_desc_n[[1]],
  495L
)
add_check(
  "Papilionoidea_descendants_495_final",
  "final",
  final_desc_n[[1]] == 495L,
  final_desc_n[[1]],
  495L
)
add_check(
  "Papilionidae_descendants_exact_7_input",
  "input",
  identical(input_desc[[2]], papilionidae_expected),
  paste(input_desc[[2]], collapse = ";"),
  paste(papilionidae_expected, collapse = ";")
)
add_check(
  "Papilionidae_descendants_exact_7_final",
  "final",
  identical(final_desc[[2]], papilionidae_expected),
  paste(final_desc[[2]], collapse = ";"),
  paste(papilionidae_expected, collapse = ";")
)
add_check(
  "Papilionidae_is_input_root_child",
  "input",
  !is.na(input_nodes[[2]]) && input_nodes[[2]] %in% input$root_children,
  input_nodes[[2]],
  paste(input$root_children, collapse = ",")
)
add_check(
  "Papilionidae_is_final_root_child",
  "final",
  !is.na(final_nodes[[2]]) && final_nodes[[2]] %in% final$root_children,
  final_nodes[[2]],
  paste(final$root_children, collapse = ",")
)

input_hesperiinae_expected <- sort(setdiff(input_tree$tip.label, papilionidae_expected))
final_hesperiinae_expected <- sort(setdiff(final_tree$tip.label, papilionidae_expected))
add_check(
  "Hesperiinae_descendants_exact_488_input",
  "input",
  identical(input_desc[[3]], input_hesperiinae_expected),
  length(input_desc[[3]]),
  488L
)
add_check(
  "Hesperiinae_descendants_exact_488_final",
  "final",
  identical(final_desc[[3]], final_hesperiinae_expected),
  length(final_desc[[3]]),
  488L
)
add_check(
  "Hesperiinae_is_input_root_child",
  "input",
  !is.na(input_nodes[[3]]) && input_nodes[[3]] %in% input$root_children,
  input_nodes[[3]],
  paste(input$root_children, collapse = ",")
)
add_check(
  "Hesperiinae_is_final_root_child",
  "final",
  !is.na(final_nodes[[3]]) && final_nodes[[3]] %in% final$root_children,
  final_nodes[[3]],
  paste(final$root_children, collapse = ",")
)

for (i in which(node_specs$calibrated == "yes")) {
  add_check(
    paste0(node_specs$node_id[[i]], "_age_within_calibration_bounds"),
    "final_age",
    within_interval_serialization[[i]],
    fmt_num(node_age[[i]]),
    paste0(
      fmt_num(node_specs$lower_bound_ma[[i]]),
      "-",
      fmt_num(node_specs$upper_bound_ma[[i]]),
      " Ma"
    ),
    TRUE,
    "Hard gate allows only the prespecified 1e-5 Ma six-decimal Newick serialization tolerance."
  )
  add_check(paste0(node_specs$node_id[[i]], "_age_within_bounds_strict_1e_6_diagnostic"), "final_age", within_interval_strict[[i]], paste0("signed_above_lower=", fmt_num(distance_above_lower[[i]]), ";signed_below_upper=", fmt_num(distance_below_upper[[i]])), "both signed distances >= -1e-6 Ma", FALSE, "Strict diagnostic only; a warning can arise solely from six-decimal Newick serialization.")
  add_check(paste0(node_specs$node_id[[i]], "_age_within_bounds_without_tolerance_diagnostic"), "final_age", within_interval_no_tolerance[[i]], paste0("signed_above_lower=", fmt_num(distance_above_lower[[i]]), ";signed_below_upper=", fmt_num(distance_below_upper[[i]])), "both signed distances >= 0 Ma", FALSE, "No-tolerance diagnostic preserves the sign of any serialized boundary excursion; it is not the hard gate.")
}
mapping <- if (file.exists(mapping_path)) {
  tryCatch(
    read.delim(mapping_path, sep = "\t", check.names = FALSE),
    error = function(e) NULL
  )
} else {
  NULL
}
mapping_columns <- c(
  "node_id", "anchor_tip_1", "anchor_tip_2", "calibrated",
  "lower_bound_ma", "upper_bound_ma"
)
mapping_readable <- !is.null(mapping) && all(mapping_columns %in% names(mapping))
add_check(
  "calibration_mapping_readable",
  "run_root",
  mapping_readable,
  if (mapping_readable) mapping_path else "missing/unreadable/invalid columns",
  mapping_path
)
if (mapping_readable) {
  mapping_node_counts <- vapply(node_specs$node_id, function(z) sum(mapping$node_id == z), integer(1))
  mapping_exactly_three_unique <- nrow(mapping) == 3L && all(mapping_node_counts == 1L) && setequal(mapping$node_id, node_specs$node_id)
  add_check(
    "calibration_mapping_exactly_one_row_per_focal_node",
    "run_root",
    mapping_exactly_three_unique,
    paste(paste(node_specs$node_id, mapping_node_counts, sep = "="), collapse = ";"),
    paste(paste(node_specs$node_id, 1L, sep = "="), collapse = ";")
  )
  mapping_rows <- match(node_specs$node_id, mapping$node_id)
  mapping_ids_complete <- all(!is.na(mapping_rows))
  add_check(
    "calibration_mapping_has_all_focal_nodes",
    "run_root",
    mapping_ids_complete,
    paste(mapping$node_id, collapse = ","),
    paste(node_specs$node_id, collapse = ",")
  )
  if (mapping_ids_complete) {
    mm <- mapping[mapping_rows, , drop = FALSE]
    add_check(
      "calibration_flags_exact_yes_yes_yes",
      "run_root",
      identical(tolower(as.character(mm$calibrated)), c("yes", "yes", "yes")),
      paste(tolower(as.character(mm$calibrated)), collapse = ","),
      "yes,yes,yes"
    )
    add_check(
      "calibration_mapping_anchors_exact",
      "run_root",
      identical(as.character(mm$anchor_tip_1), node_specs$anchor_tip_1) &&
        identical(as.character(mm$anchor_tip_2), node_specs$anchor_tip_2),
      paste(paste(mm$anchor_tip_1, mm$anchor_tip_2, sep = "+"), collapse = ";"),
      paste(
        paste(node_specs$anchor_tip_1, node_specs$anchor_tip_2, sep = "+"),
        collapse = ";"
      )
    )
    add_check(
      "Hesperiinae_mapping_explicitly_calibrated",
      "run_root",
      tolower(as.character(mm$calibrated[[3]])) == "yes",
      mm$calibrated[[3]],
      "yes"
    )
    add_check(
      "three_calibration_bounds_exact",
      "run_root",
      isTRUE(all.equal(
        as.numeric(mm$lower_bound_ma),
        node_specs$lower_bound_ma,
        tolerance = 0
      )) && isTRUE(all.equal(
        as.numeric(mm$upper_bound_ma),
        node_specs$upper_bound_ma,
        tolerance = 0
      )),
      paste(
        paste(mm$lower_bound_ma, mm$upper_bound_ma, sep = "-"),
        collapse = ";"
      ),
      "91.5046-100.8925;44.1968-52.9473;36.210861-40.662537"
    )
  }
}

checks <- do.call(rbind, checks_list)
critical_failures <- checks$critical & checks$status == "FAIL"
overall_status <- if (any(critical_failures)) "FAIL" else "PASS"

write.table(
  node_ages,
  node_ages_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)
write.table(
  checks,
  checks_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)
write.table(
  provenance,
  sha256_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

md_escape <- function(x) gsub("\\|", "\\\\|", as.character(x))
check_table <- c(
  "| Check | Scope | Critical | Status | Observed | Expected |",
  "|---|---|---:|---|---|---|",
  vapply(seq_len(nrow(checks)), function(i) {
    paste0(
      "| ", md_escape(checks$check_id[[i]]),
      " | ", md_escape(checks$scope[[i]]),
      " | ", checks$critical[[i]],
      " | ", checks$status[[i]],
      " | ", md_escape(checks$observed[[i]]),
      " | ", md_escape(checks$expected[[i]]), " |"
    )
  }, character(1))
)
node_table <- c(
  "| Node | Calibrated | Final age (Ma) | Interval (Ma) | 1e-5 hard gate | No-tolerance interval | Signed above lower (Ma) | Signed below upper (Ma) | Strict boundary hit (1e-6) | Descendants |",
  "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
  vapply(seq_len(nrow(node_ages)), function(i) {
    paste0(
      "| ", node_ages$node_id[[i]],
      " | ", node_ages$calibrated[[i]],
      " | ", fmt_num(node_ages$final_node_age_ma[[i]]),
      " | ", fmt_num(node_ages$lower_bound_ma[[i]]), "-",
      fmt_num(node_ages$upper_bound_ma[[i]]),
      " | ", node_ages$within_interval[[i]],
      " | ", node_ages$within_interval_no_tolerance[[i]],
      " | ", fmt_num(node_ages$signed_distance_above_lower_ma[[i]]),
      " | ", fmt_num(node_ages$signed_distance_below_upper_ma[[i]]),
      " | ", node_ages$any_boundary_hit_strict_1e_6[[i]],
      " | ", node_ages$final_descendant_tip_count[[i]], " |"
    )
  }, character(1))
)

failed_ids <- checks$check_id[critical_failures]
strict_boundary_hits <- node_ages$node_id[which(node_ages$any_boundary_hit_strict_1e_6 %in% TRUE)]
serialization_boundary_near <- node_ages$node_id[which(node_ages$lower_boundary_near_serialization_1e_5 | node_ages$upper_boundary_near_serialization_1e_5)]
report_lines <- c(
  "# Final chronogram core QA",
  "",
  paste0("**Overall status:** ", overall_status),
  "",
  paste0("- Input ML tree: `", input_path, "`"),
  paste0("- Final dated tree: `", final_path, "`"),
  paste0("- Primary final config: `", config_path, "`"),
  paste0("- Primary stage metadata: `", metadata_path, "`"),
  paste0("- Formal 23-point CV decision: `", cv_decision_path, "`"),
  paste0("- Source 19-point boundary decision: `", source_19point_decision_path, "`"),
  paste0("- Predeclared finite-precision plateau rule: `", plateau_rule_path, "`"),
  paste0("- Patched treePL binary: `", binary_path, "`"),
  paste0("- CV-selected smoothing: `", paste(cv_value("selected_smoothing"), collapse = ";"), "`"),
  paste0("- CV selection rule: `", paste(cv_selected_rule, collapse = ";"), "`"),
  paste0("- CV decision status: `", paste(cv_status, collapse = ";"), "`"),
  paste0("- Finite-precision plateau: `", paste(cv_value("finite_precision_plateau"), collapse = ";"), "`"),
  paste0("- Boundary still requires extension: `", paste(cv_value("boundary_requires_extension"), collapse = ";"), "`"),
  paste0("- Source 19-point lower boundary in aggregate winner set: `", paste(cv_value("diagnostic_19point_lower_boundary_in_aggregate_winner_set"), collapse = ";"), "`"),
  paste0("- Source 19-point aggregate winner values: `", paste(cv_value("diagnostic_19point_aggregate_winner_values"), collapse = ";"), "`"),
  paste0("- Final seed: `", paste(values_for(config, "seed"), collapse = ";"), "`"),
  paste0("- Independently computed SHA-256 table: `", sha256_path, "`"),
  paste0("- Strict ultrametric target (warning): `", fmt_num(STRICT_ULTRAMETRIC_TARGET_MA), " Ma`"),
  paste0("- Critical ultrametric serialization tolerance: `", fmt_num(SERIALIZATION_ULTRAMETRIC_TOL_MA), " Ma`"),
  paste0("- Strict calibration-bound diagnostic tolerance (warning): `", fmt_num(STRICT_BOUND_DIAGNOSTIC_TOL_MA), " Ma`"),
  paste0("- Critical calibration-bound serialization tolerance: `", fmt_num(SERIALIZATION_BOUND_TOL_MA), " Ma`"),
  paste0(
    "- Input root-to-tip min/max/range: `",
    paste(vapply(c(min(input$tip_depth), max(input$tip_depth), input_tip_range), fmt_num, character(1)), collapse = " / "),
    "`"
  ),
  paste0(
    "- Final root-to-tip min/max/range (Ma): `",
    paste(vapply(c(min(final$tip_depth), max(final$tip_depth), final_tip_range), fmt_num, character(1)), collapse = " / "),
    "`"
  ),
  paste0("- ape PH85 RF distance: `", fmt_num(rf_ph85), "`"),
  "",
  "Node ages are the mean node-to-descendant-tip path length. `node_ages.tsv` retains descendant-tip min/max/range, the two unmodified signed distances to the interval bounds, a no-tolerance result, a strict 1e-6 Ma diagnostic, and the prespecified 1e-5 Ma serialization hard gate.",
  "",
  "## Provenance hash caveat",
  "",
  "The stage wrapper metadata records the executed binary path, config path, timing, and exit code, but it does not embed execution-time hashes. No hash was fabricated in that metadata. `final_provenance_sha256.tsv` contains hashes computed independently by this QA script after stage completion; the predeclared input-tree, randomCV patch, and patched-binary hashes are also enforced as critical checks.",
  "",
  "## CV interpretation",
  "",
  if (isTRUE(cv_finite_precision_plateau)) {
    paste0("The 23-point result is a predeclared finite-precision tied-minimum plateau, not a unique optimum. The primary smooth is the largest (most strongly regularized) tied aggregate winner: `", paste(cv_value("selected_smoothing"), collapse = ";"), "`; aggregate tied winners: `", paste(cv_value("aggregate_winner_values"), collapse = ";"), "`.")
  } else {
    "The accepted 23-point primary selection must satisfy the unique, nonboundary aggregate-minimum contract."
  },
  paste0("Sensitivity-required flag: `", paste(cv_value("sensitivity_required"), collapse = ";"), "`. Required sensitivity executions and comparisons are delegated to `qa/final_variant_manifest.tsv` and `scripts/09_compare_final_variants_threecal.R`; this core QA validates the decision contract and primary tree only."),
  "",
  "## Focal node ages",
  "",
  node_table,
  "",
  "All three intervals are hard secondary treePL calibrations from the same source and are therefore correlated constraints, not three independent fossil observations.",
  "`CROWN_HESPERIINAE` is calibrated with the hard secondary interval 36.210861-40.662537 Ma. Its estimate is constraint-conditioned and is not an independent test of that Toussaint et al. interval.",
  "",
  paste0(
    "Strict calibration-boundary hits within ", fmt_num(STRICT_BOUND_DIAGNOSTIC_TOL_MA), " Ma: ",
    if (length(strict_boundary_hits)) paste(strict_boundary_hits, collapse = ", ") else "none"
  ),
  paste0(
    "Nodes within the ", fmt_num(SERIALIZATION_BOUND_TOL_MA), " Ma serialization neighborhood of either bound: ",
    if (length(serialization_boundary_near)) paste(serialization_boundary_near, collapse = ", ") else "none"
  ),
  "",
  "## Checks",
  "",
  check_table,
  "",
  "## Critical failures",
  "",
  if (length(failed_ids)) paste0("- `", failed_ids, "`") else "None."
)
writeLines(report_lines, report_path, useBytes = TRUE)

if (any(critical_failures)) {
  stop(
    "Final QA hard gate failed: ",
    paste(failed_ids, collapse = ", "),
    ". See ", checks_path,
    call. = FALSE
  )
}
cat("Provenance SHA-256: ", sha256_path, "\n", sep = "")

cat(sprintf("Final QA PASS: %d checks; %d warnings\n", nrow(checks), sum(checks$status == "WARN")))
cat("Node ages: ", node_ages_path, "\n", sep = "")
cat("Checks: ", checks_path, "\n", sep = "")
cat("Core report: ", report_path, "\n", sep = "")
