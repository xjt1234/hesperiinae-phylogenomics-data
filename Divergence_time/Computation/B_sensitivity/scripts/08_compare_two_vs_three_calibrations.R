#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ape))

options(stringsAsFactors = FALSE, digits = 15)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 0L) {
  stop(
    "Usage: 08_compare_two_vs_three_calibrations.R",
    call. = FALSE
  )
}

RUN_ROOT <- paste0(
  "/home/data/t200301/xjt/Hesperiinae_review/",
  "treePL_R2_3_T25_three_calibrations_20260904_075723"
)
TWO_CALIBRATION_TREE <- paste0(
  "/home/data/t200301/xjt/Hesperiinae_review/",
  "treePL_R2_3_T25_deep_only_20260904_000707/",
  "output/T25_deep_only_dated.tre"
)
THREE_CALIBRATION_TREE <- file.path(
  RUN_ROOT,
  "output",
  "T25_three_calibrations_dated.tre"
)

run_root <- normalizePath(RUN_ROOT, mustWork = TRUE)
two_path <- normalizePath(TWO_CALIBRATION_TREE, mustWork = TRUE)
three_path <- normalizePath(THREE_CALIBRATION_TREE, mustWork = TRUE)

qa_dir <- file.path(run_root, "qa")
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)

node_ages_path <- file.path(
  qa_dir,
  "two_vs_three_calibration_node_ages.tsv"
)
checks_path <- file.path(
  qa_dir,
  "two_vs_three_calibration_checks.tsv"
)
summary_path <- file.path(
  qa_dir,
  "two_vs_three_calibration_summary.md"
)
requested_outputs <- c(node_ages_path, checks_path, summary_path)
existing_outputs <- requested_outputs[file.exists(requested_outputs)]
if (length(existing_outputs)) {
  stop(
    "Refusing to overwrite existing two-versus-three calibration output(s): ",
    paste(existing_outputs, collapse = ", "),
    call. = FALSE
  )
}

EXPECTED_TIP_COUNT <- 495L
ULTRAMETRIC_SERIALIZATION_TOL_MA <- 1e-5
CALIBRATION_BOUND_TOL_MA <- 1e-6
CLUSTER_SEPARATOR <- "\034"

fmt_num <- function(x) {
  if (length(x) != 1L || is.na(x) || !is.finite(x)) return("NA")
  sprintf("%.12g", x)
}

fmt_md_num <- function(x) {
  if (length(x) != 1L || is.na(x) || !is.finite(x)) return("NA")
  sprintf("%.9f", x)
}

fmt_bool <- function(x) {
  if (length(x) != 1L || is.na(x)) return("NA")
  if (isTRUE(x)) "TRUE" else "FALSE"
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
    tree_name = tree_name,
    n_tip = n_tip,
    root = root,
    root_children = root_children,
    depth = depth,
    tip_depth = depth[seq_len(n_tip)],
    descendant_ids = descendant_ids,
    descendant_labels = descendant_labels
  )
}

label_counts <- function(tree, labels) {
  vapply(labels, function(label) sum(tree$tip.label == label), integer(1))
}

safe_mrca <- function(tree, labels) {
  if (!all(label_counts(tree, labels) == 1L)) return(NA_integer_)
  as.integer(getMRCA(tree, labels))
}

safe_descendants <- function(context, node) {
  if (length(node) != 1L || is.na(node)) return(character())
  context$descendant_labels(node)
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

internal_node_table <- function(context) {
  nodes <- sort(unique(as.integer(context$tree$edge[, 1])))
  pieces <- lapply(nodes, function(node) {
    labels <- context$descendant_labels(node)
    age <- node_age_stats(context, node)
    data.frame(
      node = node,
      descendant_tip_count = length(labels),
      cluster_key = paste(labels, collapse = CLUSTER_SEPARATOR),
      age_ma = unname(age[["mean"]]),
      age_tip_min_ma = unname(age[["min"]]),
      age_tip_max_ma = unname(age[["max"]]),
      age_tip_range_ma = unname(age[["range"]]),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, pieces)
}

boundary_label <- function(age, lower, upper) {
  if (!is.finite(age)) return("NA")
  lower_hit <- abs(age - lower) <= CALIBRATION_BOUND_TOL_MA
  upper_hit <- abs(upper - age) <= CALIBRATION_BOUND_TOL_MA
  if (lower_hit && upper_hit) return("both")
  if (lower_hit) return("lower")
  if (upper_hit) return("upper")
  "none"
}

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

add_info <- function(check_id, scope, observed, note = "") {
  checks_list[[length(checks_list) + 1L]] <<- data.frame(
    check_id = check_id,
    scope = scope,
    critical = FALSE,
    status = "INFO",
    observed = as.character(observed),
    expected = "reported; no pass/fail threshold",
    note = as.character(note),
    stringsAsFactors = FALSE
  )
}

two_tree <- read.tree(two_path)
three_tree <- read.tree(three_path)
two <- make_context(two_tree, "existing two-calibration primary tree")
three <- make_context(three_tree, "new three-calibration primary tree")

tree_specs <- list(
  old_two_calibrations = list(tree = two_tree, context = two, path = two_path),
  new_three_calibrations = list(
    tree = three_tree,
    context = three,
    path = three_path
  )
)

for (scenario_id in names(tree_specs)) {
  spec <- tree_specs[[scenario_id]]
  tree <- spec$tree
  context <- spec$context
  root_split <- sort(vapply(
    context$root_children,
    function(node) length(context$descendant_ids(node)),
    integer(1)
  ))
  root_to_tip_range <- diff(range(context$tip_depth))

  add_check(
    paste0(scenario_id, "_file_exists_nonzero"), scenario_id,
    file.info(spec$path)$size > 0,
    file.info(spec$path)$size, ">0 bytes"
  )
  add_check(
    paste0(scenario_id, "_tip_count_495"), scenario_id,
    context$n_tip == EXPECTED_TIP_COUNT,
    context$n_tip, EXPECTED_TIP_COUNT
  )
  add_check(
    paste0(scenario_id, "_unique_tip_labels"), scenario_id,
    anyDuplicated(tree$tip.label) == 0L,
    anyDuplicated(tree$tip.label), 0L
  )
  add_check(
    paste0(scenario_id, "_rooted"), scenario_id,
    is.rooted(tree), is.rooted(tree), TRUE
  )
  add_check(
    paste0(scenario_id, "_binary"), scenario_id,
    is.binary.phylo(tree), is.binary.phylo(tree), TRUE
  )
  add_check(
    paste0(scenario_id, "_branch_lengths_finite"), scenario_id,
    all(is.finite(tree$edge.length)),
    sum(!is.finite(tree$edge.length)), 0L
  )
  add_check(
    paste0(scenario_id, "_branch_lengths_nonnegative"), scenario_id,
    all(tree$edge.length >= 0),
    sum(tree$edge.length < 0), 0L
  )
  add_check(
    paste0(scenario_id, "_root_split_7_488"), scenario_id,
    identical(root_split, c(7L, 488L)),
    paste(root_split, collapse = ","), "7,488"
  )
  add_check(
    paste0(scenario_id, "_root_to_tip_range_le_1e-5_ma"), scenario_id,
    is.finite(root_to_tip_range) &&
      root_to_tip_range <= ULTRAMETRIC_SERIALIZATION_TOL_MA,
    paste0(
      "min=", fmt_num(min(context$tip_depth)),
      ";max=", fmt_num(max(context$tip_depth)),
      ";range=", fmt_num(root_to_tip_range)
    ),
    paste0(
      "range<=", fmt_num(ULTRAMETRIC_SERIALIZATION_TOL_MA), " Ma"
    ),
    note = "Absolute tolerance accommodates six-decimal Newick serialization."
  )
}

same_tip_set <- (
  length(two_tree$tip.label) == length(three_tree$tip.label) &&
    setequal(two_tree$tip.label, three_tree$tip.label)
)
add_check(
  "same_495_tip_label_set", "cross_tree",
  same_tip_set && length(two_tree$tip.label) == EXPECTED_TIP_COUNT,
  paste0(
    "old_only=", length(setdiff(two_tree$tip.label, three_tree$tip.label)),
    ";new_only=", length(setdiff(three_tree$tip.label, two_tree$tip.label)),
    ";old_n=", length(two_tree$tip.label),
    ";new_n=", length(three_tree$tip.label)
  ),
  "old_only=0;new_only=0;old_n=495;new_n=495"
)

rf_ph85 <- if (
    same_tip_set &&
      !anyDuplicated(two_tree$tip.label) &&
      !anyDuplicated(three_tree$tip.label)) {
  suppressWarnings(as.numeric(dist.topo(two_tree, three_tree, method = "PH85")))
} else {
  NA_real_
}
add_check(
  "ape_PH85_RF_zero", "cross_tree",
  is.finite(rf_ph85) && rf_ph85 == 0,
  fmt_num(rf_ph85), 0,
  note = "RF=0 is required but is supplemented by the root-aware clade check."
)

two_internal <- internal_node_table(two)
three_internal <- internal_node_table(three)
two_cluster_keys <- sort(two_internal$cluster_key)
three_cluster_keys <- sort(three_internal$cluster_key)
rooted_clades_equal <- (
  same_tip_set &&
    identical(two_cluster_keys, three_cluster_keys)
)
add_check(
  "rooted_internal_clade_sets_identical", "cross_tree",
  rooted_clades_equal,
  paste0(
    "old_internal=", nrow(two_internal),
    ";new_internal=", nrow(three_internal),
    ";old_only=", length(setdiff(two_cluster_keys, three_cluster_keys)),
    ";new_only=", length(setdiff(three_cluster_keys, two_cluster_keys))
  ),
  "old_internal=494;new_internal=494;old_only=0;new_only=0",
  note = paste0(
    "Each rooted clade is keyed by its exact sorted descendant-tip label set; ",
    "this distinguishes root placement even when unrooted RF is zero."
  )
)

all_equal_rooted_topology <- if (same_tip_set) {
  isTRUE(all.equal.phylo(
    two_tree,
    three_tree,
    use.edge.length = FALSE,
    use.tip.label = TRUE
  ))
} else {
  FALSE
}
add_check(
  "ape_all_equal_rooted_topology", "cross_tree",
  all_equal_rooted_topology,
  all_equal_rooted_topology, TRUE,
  note = "Independent topology check with edge lengths ignored."
)

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
  lower_bound_ma = c(91.5046, 44.1968, 36.210861),
  upper_bound_ma = c(100.8925, 52.9473, 40.662537),
  old_two_calibrated = c(TRUE, TRUE, FALSE),
  new_three_calibrated = c(TRUE, TRUE, TRUE),
  stringsAsFactors = FALSE
)

required_anchors <- unique(c(node_specs$anchor_tip_1, node_specs$anchor_tip_2))
for (scenario_id in names(tree_specs)) {
  tree <- tree_specs[[scenario_id]]$tree
  counts <- label_counts(tree, required_anchors)
  add_check(
    paste0(scenario_id, "_required_focal_anchors_exactly_once"), scenario_id,
    all(counts == 1L),
    paste(paste(names(counts), counts, sep = "="), collapse = ";"),
    paste(paste(names(counts), rep(1L, length(counts)), sep = "="), collapse = ";")
  )
}

expected_descendants <- list(
  sort(two_tree$tip.label),
  papilionidae_expected,
  sort(setdiff(two_tree$tip.label, papilionidae_expected))
)

scenario_nodes <- list()
scenario_descendants <- list()
scenario_ages <- list()

for (scenario_id in names(tree_specs)) {
  tree <- tree_specs[[scenario_id]]$tree
  context <- tree_specs[[scenario_id]]$context
  nodes <- vapply(seq_len(nrow(node_specs)), function(i) {
    safe_mrca(
      tree,
      c(node_specs$anchor_tip_1[[i]], node_specs$anchor_tip_2[[i]])
    )
  }, integer(1))
  descendants <- lapply(nodes, function(node) safe_descendants(context, node))
  ages <- t(vapply(
    nodes,
    function(node) node_age_stats(context, node),
    numeric(4)
  ))
  colnames(ages) <- c("mean", "min", "max", "range")

  scenario_nodes[[scenario_id]] <- nodes
  scenario_descendants[[scenario_id]] <- descendants
  scenario_ages[[scenario_id]] <- ages

  for (i in seq_len(nrow(node_specs))) {
    node_id <- node_specs$node_id[[i]]
    descendant_set_exact <- identical(descendants[[i]], expected_descendants[[i]])
    expected_position <- if (i == 1L) context$root else context$root_children
    position_exact <- if (i == 1L) {
      !is.na(nodes[[i]]) && nodes[[i]] == context$root
    } else {
      !is.na(nodes[[i]]) && nodes[[i]] %in% context$root_children
    }

    add_check(
      paste0(scenario_id, "_", node_id, "_descendant_set_exact"),
      scenario_id,
      descendant_set_exact &&
        length(descendants[[i]]) == node_specs$expected_descendant_tip_count[[i]],
      paste0("descendants=", length(descendants[[i]])),
      paste0(
        "exact_set;descendants=",
        node_specs$expected_descendant_tip_count[[i]]
      )
    )
    add_check(
      paste0(scenario_id, "_", node_id, "_position_exact"),
      scenario_id,
      position_exact,
      if (is.na(nodes[[i]])) "NA" else as.character(nodes[[i]]),
      paste(expected_position, collapse = ",")
    )
    add_check(
      paste0(scenario_id, "_", node_id, "_age_finite"),
      scenario_id,
      is.finite(ages[i, "mean"]),
      fmt_num(ages[i, "mean"]),
      "finite Ma"
    )

    is_calibrated <- if (scenario_id == "old_two_calibrations") {
      node_specs$old_two_calibrated[[i]]
    } else {
      node_specs$new_three_calibrated[[i]]
    }
    within_bounds <- (
      is.finite(ages[i, "mean"]) &&
        ages[i, "mean"] >=
          node_specs$lower_bound_ma[[i]] - CALIBRATION_BOUND_TOL_MA &&
        ages[i, "mean"] <=
          node_specs$upper_bound_ma[[i]] + CALIBRATION_BOUND_TOL_MA
    )
    if (is_calibrated) {
      add_check(
        paste0(scenario_id, "_", node_id, "_within_calibration_bounds"),
        scenario_id,
        within_bounds,
        fmt_num(ages[i, "mean"]),
        paste0(
          fmt_num(node_specs$lower_bound_ma[[i]]), "-",
          fmt_num(node_specs$upper_bound_ma[[i]]), " Ma"
        ),
        note = paste0(
          "Inclusive hard bounds with ",
          fmt_num(CALIBRATION_BOUND_TOL_MA), " Ma numeric tolerance."
        )
      )
    } else {
      add_info(
        paste0(scenario_id, "_", node_id, "_within_reference_interval"),
        scenario_id,
        paste0(
          "within=", fmt_bool(within_bounds),
          ";age_ma=", fmt_num(ages[i, "mean"]),
          ";interval=", fmt_num(node_specs$lower_bound_ma[[i]]), "-",
          fmt_num(node_specs$upper_bound_ma[[i]])
        ),
        note = paste0(
          "This node was not calibrated in the two-calibration run; the ",
          "Toussaint interval is reported only as an external comparison."
        )
      )
    }
    add_info(
      paste0(scenario_id, "_", node_id, "_boundary_hit"),
      scenario_id,
      boundary_label(
        ages[i, "mean"],
        node_specs$lower_bound_ma[[i]],
        node_specs$upper_bound_ma[[i]]
      ),
      note = paste0(
        "Boundary hit uses mean node age and absolute tolerance ",
        fmt_num(CALIBRATION_BOUND_TOL_MA), " Ma; calibration status=",
        is_calibrated, "."
      )
    )
  }
}

old_nodes <- scenario_nodes$old_two_calibrations
new_nodes <- scenario_nodes$new_three_calibrations
old_descendants <- scenario_descendants$old_two_calibrations
new_descendants <- scenario_descendants$new_three_calibrations
old_ages <- scenario_ages$old_two_calibrations
new_ages <- scenario_ages$new_three_calibrations

node_rows <- lapply(seq_len(nrow(node_specs)), function(i) {
  old_age <- old_ages[i, "mean"]
  new_age <- new_ages[i, "mean"]
  lower <- node_specs$lower_bound_ma[[i]]
  upper <- node_specs$upper_bound_ma[[i]]
  old_within <- is.finite(old_age) &&
    old_age >= lower - CALIBRATION_BOUND_TOL_MA &&
    old_age <= upper + CALIBRATION_BOUND_TOL_MA
  new_within <- is.finite(new_age) &&
    new_age >= lower - CALIBRATION_BOUND_TOL_MA &&
    new_age <= upper + CALIBRATION_BOUND_TOL_MA

  data.frame(
    node_id = node_specs$node_id[[i]],
    biological_definition = node_specs$biological_definition[[i]],
    anchor_tip_1 = node_specs$anchor_tip_1[[i]],
    anchor_tip_2 = node_specs$anchor_tip_2[[i]],
    expected_descendant_tip_count =
      node_specs$expected_descendant_tip_count[[i]],
    old_two_calibrated = node_specs$old_two_calibrated[[i]],
    new_three_calibrated = node_specs$new_three_calibrated[[i]],
    lower_bound_ma = lower,
    upper_bound_ma = upper,
    old_two_calibration_mrca_node = old_nodes[[i]],
    new_three_calibration_mrca_node = new_nodes[[i]],
    old_two_calibration_descendant_tip_count = length(old_descendants[[i]]),
    new_three_calibration_descendant_tip_count = length(new_descendants[[i]]),
    old_two_calibration_node_age_ma = old_age,
    new_three_calibration_node_age_ma = new_age,
    new_minus_old_ma = new_age - old_age,
    absolute_new_minus_old_ma = abs(new_age - old_age),
    old_two_calibration_node_age_tip_min_ma = old_ages[i, "min"],
    old_two_calibration_node_age_tip_max_ma = old_ages[i, "max"],
    old_two_calibration_node_age_tip_range_ma = old_ages[i, "range"],
    new_three_calibration_node_age_tip_min_ma = new_ages[i, "min"],
    new_three_calibration_node_age_tip_max_ma = new_ages[i, "max"],
    new_three_calibration_node_age_tip_range_ma = new_ages[i, "range"],
    old_two_calibration_within_interval = old_within,
    new_three_calibration_within_interval = new_within,
    old_two_calibration_distance_above_lower_ma = old_age - lower,
    old_two_calibration_distance_below_upper_ma = upper - old_age,
    new_three_calibration_distance_above_lower_ma = new_age - lower,
    new_three_calibration_distance_below_upper_ma = upper - new_age,
    old_two_calibration_lower_boundary_hit =
      is.finite(old_age) && abs(old_age - lower) <= CALIBRATION_BOUND_TOL_MA,
    old_two_calibration_upper_boundary_hit =
      is.finite(old_age) && abs(upper - old_age) <= CALIBRATION_BOUND_TOL_MA,
    old_two_calibration_boundary_hit = boundary_label(old_age, lower, upper),
    new_three_calibration_lower_boundary_hit =
      is.finite(new_age) && abs(new_age - lower) <= CALIBRATION_BOUND_TOL_MA,
    new_three_calibration_upper_boundary_hit =
      is.finite(new_age) && abs(upper - new_age) <= CALIBRATION_BOUND_TOL_MA,
    new_three_calibration_boundary_hit = boundary_label(new_age, lower, upper),
    calibration_boundary_hit_tolerance_ma = CALIBRATION_BOUND_TOL_MA,
    age_definition = paste0(
      "mean node-to-descendant-tip path length; min/max/range retained for ",
      "serialization diagnostics"
    ),
    stringsAsFactors = FALSE
  )
})
node_ages <- do.call(rbind, node_rows)

new_index_for_old_clade <- match(two_internal$cluster_key, three_internal$cluster_key)
all_internal_clades_matched <- (
  rooted_clades_equal &&
    length(new_index_for_old_clade) == nrow(two_internal) &&
    all(!is.na(new_index_for_old_clade))
)
add_check(
  "all_internal_clades_matched_for_age_comparison", "cross_tree",
  all_internal_clades_matched,
  paste0(
    "matched=", sum(!is.na(new_index_for_old_clade)),
    ";old_internal=", nrow(two_internal),
    ";new_internal=", nrow(three_internal)
  ),
  "matched=494;old_internal=494;new_internal=494"
)

if (all_internal_clades_matched) {
  internal_delta <-
    three_internal$age_ma[new_index_for_old_clade] - two_internal$age_ma
  max_index <- which.max(abs(internal_delta))
  max_abs_internal_delta <- abs(internal_delta[[max_index]])
  signed_delta_at_max <- internal_delta[[max_index]]
  mean_abs_internal_delta <- mean(abs(internal_delta))
  mean_signed_internal_delta <- mean(internal_delta)
  median_abs_internal_delta <- median(abs(internal_delta))
  old_node_at_max <- two_internal$node[[max_index]]
  new_node_at_max <- three_internal$node[new_index_for_old_clade[[max_index]]]
  descendants_at_max <- two_internal$descendant_tip_count[[max_index]]
} else {
  internal_delta <- numeric()
  max_abs_internal_delta <- NA_real_
  signed_delta_at_max <- NA_real_
  mean_abs_internal_delta <- NA_real_
  mean_signed_internal_delta <- NA_real_
  median_abs_internal_delta <- NA_real_
  old_node_at_max <- NA_integer_
  new_node_at_max <- NA_integer_
  descendants_at_max <- NA_integer_
}

add_info(
  "all_internal_node_age_difference_summary", "cross_tree",
  paste0(
    "n=", length(internal_delta),
    ";max_abs_ma=", fmt_num(max_abs_internal_delta),
    ";mean_abs_ma=", fmt_num(mean_abs_internal_delta),
    ";mean_signed_new_minus_old_ma=", fmt_num(mean_signed_internal_delta),
    ";median_abs_ma=", fmt_num(median_abs_internal_delta)
  ),
  note = paste0(
    "Internal nodes were paired only by identical rooted descendant-tip sets; ",
    "new-minus-old is three-calibration minus two-calibration age."
  )
)
add_info(
  "max_abs_all_internal_node_age_difference_ma", "cross_tree",
  fmt_num(max_abs_internal_delta),
  note = "Maximum absolute new-minus-old age change across exact matched rooted clades."
)
add_info(
  "mean_abs_all_internal_node_age_difference_ma", "cross_tree",
  fmt_num(mean_abs_internal_delta),
  note = "Mean absolute new-minus-old age change across exact matched rooted clades."
)
add_info(
  "mean_signed_all_internal_node_age_difference_ma", "cross_tree",
  fmt_num(mean_signed_internal_delta),
  note = "Mean signed age change; positive values mean the three-calibration tree is older."
)
add_info(
  "max_absolute_internal_node_age_difference_location", "cross_tree",
  paste0(
    "max_abs_ma=", fmt_num(max_abs_internal_delta),
    ";signed_new_minus_old_ma=", fmt_num(signed_delta_at_max),
    ";old_node=", old_node_at_max,
    ";new_node=", new_node_at_max,
    ";descendant_tips=", descendants_at_max
  ),
  note = "Location of the maximum absolute age change among matched internal clades."
)

checks <- do.call(rbind, checks_list)

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

critical_failures <- checks$critical & checks$status == "FAIL"
overall_status <- if (any(critical_failures)) "FAIL" else "PASS"
failed_ids <- checks$check_id[critical_failures]

md_table_rows <- vapply(seq_len(nrow(node_ages)), function(i) {
  row <- node_ages[i, , drop = FALSE]
  paste0(
    "| ", row$node_id, " | ",
    fmt_md_num(row$old_two_calibration_node_age_ma), " | ",
    fmt_md_num(row$new_three_calibration_node_age_ma), " | ",
    fmt_md_num(row$new_minus_old_ma), " | ",
    fmt_bool(row$old_two_calibration_within_interval), " | ",
    row$old_two_calibration_boundary_hit, " | ",
    fmt_bool(row$new_three_calibration_within_interval), " | ",
    row$new_three_calibration_boundary_hit, " |"
  )
}, character(1))

summary_lines <- c(
  "# Two- versus three-calibration treePL comparison",
  "",
  paste0("Overall structural/constraint QA: **", overall_status, "**"),
  "",
  "## Compared trees",
  "",
  paste0("- Existing two-calibration primary tree: `", two_path, "`"),
  paste0("- New three-calibration primary tree: `", three_path, "`"),
  "- Both trees are expected to contain the same 495 uniquely labelled tips.",
  paste0(
    "- Node age is the mean node-to-descendant-tip path length; descendant-tip ",
    "minimum, maximum, and range are retained in the TSV as serialization diagnostics."
  ),
  "",
  "## Focal-node ages and calibration diagnostics",
  "",
  paste0(
    "Boundary hits use an absolute tolerance of ",
    fmt_num(CALIBRATION_BOUND_TOL_MA), " Ma. In the existing two-calibration ",
    "tree, the Hesperiinae interval is an external reference only; in the new ",
    "tree it is a hard treePL calibration."
  ),
  "",
  paste0(
    "| Node | Two-calibration age (Ma) | Three-calibration age (Ma) | ",
    "New-old (Ma) | Old within interval | Old boundary | ",
    "New within interval | New boundary |"
  ),
  "|---|---:|---:|---:|:---:|:---:|:---:|:---:|",
  md_table_rows,
  "",
  "Intervals assessed:",
  "",
  "- crown Papilionoidea: 91.5046-100.8925 Ma",
  "- crown Papilionidae excluding *Baronia*: 44.1968-52.9473 Ma",
  "- crown Hesperiinae: 36.210861-40.662537 Ma",
  "",
  "## Topology and all-internal-node age comparison",
  "",
  paste0("- ape PH85 RF distance: ", fmt_num(rf_ph85)),
  paste0("- Exact rooted internal-clade sets identical: ", rooted_clades_equal),
  paste0("- Matched internal clades: ", length(internal_delta)),
  paste0(
    "- Maximum absolute internal-node age change: ",
    fmt_md_num(max_abs_internal_delta), " Ma"
  ),
  paste0(
    "- Signed change at that maximum (new minus old): ",
    fmt_md_num(signed_delta_at_max), " Ma"
  ),
  paste0(
    "- Mean absolute internal-node age change: ",
    fmt_md_num(mean_abs_internal_delta), " Ma"
  ),
  paste0(
    "- Mean signed internal-node age change (new minus old): ",
    fmt_md_num(mean_signed_internal_delta), " Ma"
  ),
  paste0(
    "- Median absolute internal-node age change: ",
    fmt_md_num(median_abs_internal_delta), " Ma"
  ),
  paste0(
    "- Maximum-change clade: old node ", old_node_at_max,
    ", new node ", new_node_at_max,
    ", descendant tips ", descendants_at_max
  ),
  "",
  "Internal nodes are paired by exact sorted descendant-tip label sets. No node-number equivalence is assumed.",
  "",
  "## Output files",
  "",
  paste0("- Focal-node table: `", node_ages_path, "`"),
  paste0("- QA checks and internal-node summaries: `", checks_path, "`"),
  paste0("- This summary: `", summary_path, "`")
)

if (any(critical_failures)) {
  summary_lines <- c(
    summary_lines,
    "",
    "## Critical failures",
    "",
    paste0("- ", failed_ids)
  )
}

writeLines(summary_lines, summary_path, useBytes = TRUE)

if (any(critical_failures)) {
  stop(
    "Two-versus-three calibration comparison hard gate failed: ",
    paste(failed_ids, collapse = ", "),
    ". Diagnostic outputs were retained in ", qa_dir,
    call. = FALSE
  )
}

cat("Two-versus-three calibration comparison PASS\n")
cat("Focal-node ages: ", node_ages_path, "\n", sep = "")
cat("Checks: ", checks_path, "\n", sep = "")
cat("Summary: ", summary_path, "\n", sep = "")
