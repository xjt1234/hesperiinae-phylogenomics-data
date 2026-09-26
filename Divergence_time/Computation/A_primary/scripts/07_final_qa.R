#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ape))

options(stringsAsFactors = FALSE, digits = 15)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "Usage: 07_final_qa.R <input_ml_tree> <final_dated_tree> <run_root>",
    call. = FALSE
  )
}

input_path <- normalizePath(args[[1]], mustWork = TRUE)
final_path <- normalizePath(args[[2]], mustWork = TRUE)
run_root <- normalizePath(args[[3]], mustWork = TRUE)

qa_dir <- file.path(run_root, "qa")
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)

node_ages_path <- file.path(qa_dir, "node_ages.tsv")
checks_path <- file.path(qa_dir, "final_tree_checks.tsv")
report_path <- file.path(qa_dir, "qa_report_core.md")
requested_outputs <- c(node_ages_path, checks_path, report_path)
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
BOUND_TOL_MA <- 1e-6

fmt_num <- function(x) {
  if (length(x) != 1L || is.na(x) || !is.finite(x)) return("NA")
  sprintf("%.12g", x)
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
  calibrated = c("yes", "yes", "no"),
  interval_role = c(
    "treePL calibration",
    "treePL calibration",
    "external comparison only; not calibrated"
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
within_interval <- (
  node_age >= node_specs$lower_bound_ma - BOUND_TOL_MA &
    node_age <= node_specs$upper_bound_ma + BOUND_TOL_MA
)
lower_boundary_hit <- abs(distance_above_lower) <= BOUND_TOL_MA
upper_boundary_hit <- abs(distance_below_upper) <= BOUND_TOL_MA

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
  within_interval = within_interval,
  distance_above_lower_ma = distance_above_lower,
  distance_below_upper_ma = distance_below_upper,
  nearest_boundary_distance_ma = nearest_boundary_distance,
  lower_boundary_hit = lower_boundary_hit,
  upper_boundary_hit = upper_boundary_hit,
  any_boundary_hit = lower_boundary_hit | upper_boundary_hit,
  boundary_hit_tolerance_ma = BOUND_TOL_MA,
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
    within_interval[[i]],
    fmt_num(node_age[[i]]),
    paste0(
      fmt_num(node_specs$lower_bound_ma[[i]]),
      "-",
      fmt_num(node_specs$upper_bound_ma[[i]]),
      " Ma"
    )
  )
}
add_check(
  "CROWN_HESPERIINAE_external_range_comparison",
  "final_age",
  within_interval[[3]],
  fmt_num(node_age[[3]]),
  "36.210861-40.662537 Ma (comparison only)",
  FALSE,
  "This interval is not a treePL calibration and must not be treated as a hard age constraint."
)

mapping_path <- file.path(run_root, "mapping", "calibration_node_mapping.tsv")
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
      "calibration_flags_exact_yes_yes_no",
      "run_root",
      identical(tolower(as.character(mm$calibrated)), c("yes", "yes", "no")),
      paste(tolower(as.character(mm$calibrated)), collapse = ","),
      "yes,yes,no"
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
      "Hesperiinae_mapping_explicitly_uncalibrated",
      "run_root",
      tolower(as.character(mm$calibrated[[3]])) == "no",
      mm$calibrated[[3]],
      "no"
    )
    add_check(
      "deep_calibration_bounds_exact",
      "run_root",
      isTRUE(all.equal(
        as.numeric(mm$lower_bound_ma[1:2]),
        node_specs$lower_bound_ma[1:2],
        tolerance = 0
      )) && isTRUE(all.equal(
        as.numeric(mm$upper_bound_ma[1:2]),
        node_specs$upper_bound_ma[1:2],
        tolerance = 0
      )),
      paste(
        paste(mm$lower_bound_ma[1:2], mm$upper_bound_ma[1:2], sep = "-"),
        collapse = ";"
      ),
      "91.5046-100.8925;44.1968-52.9473"
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
  "| Node | Calibrated | Final age (Ma) | Interval (Ma) | Within interval | Boundary hit | Descendants |",
  "|---|---:|---:|---:|---:|---:|---:|",
  vapply(seq_len(nrow(node_ages)), function(i) {
    paste0(
      "| ", node_ages$node_id[[i]],
      " | ", node_ages$calibrated[[i]],
      " | ", fmt_num(node_ages$final_node_age_ma[[i]]),
      " | ", fmt_num(node_ages$lower_bound_ma[[i]]), "-",
      fmt_num(node_ages$upper_bound_ma[[i]]),
      " | ", node_ages$within_interval[[i]],
      " | ", node_ages$any_boundary_hit[[i]],
      " | ", node_ages$final_descendant_tip_count[[i]], " |"
    )
  }, character(1))
)

failed_ids <- checks$check_id[critical_failures]
boundary_hits <- node_ages$node_id[which(node_ages$any_boundary_hit %in% TRUE)]
report_lines <- c(
  "# Final chronogram core QA",
  "",
  paste0("**Overall status:** ", overall_status),
  "",
  paste0("- Input ML tree: `", input_path, "`"),
  paste0("- Final dated tree: `", final_path, "`"),
  paste0("- Strict ultrametric target (warning): `", fmt_num(STRICT_ULTRAMETRIC_TARGET_MA), " Ma`"),
  paste0("- Critical serialization tolerance: `", fmt_num(SERIALIZATION_ULTRAMETRIC_TOL_MA), " Ma`"),
  paste0("- Calibration-bound tolerance: `", fmt_num(BOUND_TOL_MA), " Ma`"),
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
  "Node ages are the mean node-to-descendant-tip path length. The corresponding min, max, and range are retained in `node_ages.tsv` so Newick rounding cannot be mistaken for biological age variation.",
  "",
  "## Focal node ages",
  "",
  node_table,
  "",
  "`CROWN_HESPERIINAE` is not calibrated. Its 36.210861-40.662537 Ma interval is an external comparison range only; an out-of-range value is a warning, not a hard calibration failure.",
  "",
  paste0(
    "Calibration-boundary hits within ", fmt_num(BOUND_TOL_MA), " Ma: ",
    if (length(boundary_hits)) paste(boundary_hits, collapse = ", ") else "none"
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

cat(sprintf("Final QA PASS: %d checks; %d warnings\n", nrow(checks), sum(checks$status == "WARN")))
cat("Node ages: ", node_ages_path, "\n", sep = "")
cat("Checks: ", checks_path, "\n", sep = "")
cat("Core report: ", report_path, "\n", sep = "")
