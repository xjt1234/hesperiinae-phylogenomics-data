#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ape))

options(stringsAsFactors = FALSE, digits = 15)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop(
    "Usage: 09_compare_final_variants.R <run_root>",
    call. = FALSE
  )
}

run_root <- normalizePath(args[[1]], mustWork = TRUE)
qa_dir <- file.path(run_root, "qa")
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)

node_ages_path <- file.path(qa_dir, "final_variant_node_ages.tsv")
checks_path <- file.path(qa_dir, "final_variant_checks.tsv")
requested_outputs <- c(node_ages_path, checks_path)
existing_outputs <- requested_outputs[file.exists(requested_outputs)]
if (length(existing_outputs)) {
  stop(
    "Refusing to overwrite existing final-variant QA output(s): ",
    paste(existing_outputs, collapse = ", "),
    call. = FALSE
  )
}

ULTRAMETRIC_SERIALIZATION_TOL_MA <- 1e-5
CALIBRATION_BOUND_TOL_MA <- 1e-6

variant_specs <- data.frame(
  variant_id = c(
    "primary_thorough_smooth_1e-10",
    "baseline_no_thorough_smooth_1e-10",
    "thorough_smooth_1e-12",
    "thorough_smooth_1e-13",
    "repeat_same_seed_primary_thorough_smooth_1e-10"
  ),
  variant_role = c("primary", "optimization_diagnostic", "smoothing_sensitivity", "smoothing_sensitivity", "reproducibility"),
  is_primary = c(TRUE, FALSE, FALSE, FALSE, FALSE),
  thorough = c(TRUE, FALSE, TRUE, TRUE, TRUE),
  smooth = c("1e-10", "1e-10", "1e-12", "1e-13", "1e-10"),
  relative_tree_path = c(
    "output/T25_deep_only_dated.tre",
    "output/T25_deep_only_dated_baseline_no_thorough.tre",
    "output/T25_deep_only_dated_sensitivity_1e-12.tre",
    "output/T25_deep_only_dated_sensitivity_1e-13.tre",
    "output/T25_deep_only_dated_repeat_same_seed.tre"
  ),
  stringsAsFactors = FALSE
)

input_path <- normalizePath(
  file.path(run_root, "input", "verified_495_tip_input.treefile"),
  mustWork = TRUE
)
variant_specs$tree_path <- vapply(
  file.path(run_root, variant_specs$relative_tree_path),
  normalizePath,
  character(1),
  mustWork = TRUE
)

files_byte_identical <- function(path_a, path_b) {
  size_a <- file.info(path_a)$size
  size_b <- file.info(path_b)$size
  if (is.na(size_a) || is.na(size_b) || size_a != size_b) return(FALSE)
  identical(
    readBin(path_a, what = "raw", n = size_a),
    readBin(path_b, what = "raw", n = size_b)
  )
}
primary_tree_path <- variant_specs$tree_path[variant_specs$is_primary][[1]]
variant_specs$byte_identical_to_primary <- vapply(
  variant_specs$tree_path, files_byte_identical, logical(1),
  path_b = primary_tree_path
)

make_context <- function(tree, tree_name) {
  if (is.null(tree$edge.length)) {
    stop(tree_name, " has no branch lengths", call. = FALSE)
  }
  root <- setdiff(unique(tree$edge[, 1]), unique(tree$edge[, 2]))
  if (length(root) != 1L) {
    stop("Could not identify exactly one root in ", tree_name, call. = FALSE)
  }
  root <- as.integer(root[[1]])
  n_tip <- Ntip(tree)

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
  list(
    tree = tree,
    n_tip = n_tip,
    root = root,
    depth = depth,
    tip_depth = depth[seq_len(n_tip)],
    root_children = as.integer(tree$edge[tree$edge[, 1] == root, 2]),
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
      cluster_key = paste(labels, collapse = "\034"),
      age_ma = unname(age[["mean"]]),
      age_tip_min_ma = unname(age[["min"]]),
      age_tip_max_ma = unname(age[["max"]]),
      age_tip_range_ma = unname(age[["range"]]),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, pieces)
}

fmt_num <- function(x) {
  if (length(x) != 1L || is.na(x) || !is.finite(x)) return("NA")
  sprintf("%.12g", x)
}

input_tree <- read.tree(input_path)
input_context <- make_context(input_tree, "verified input ML tree")
if (Ntip(input_tree) != 495L || anyDuplicated(input_tree$tip.label)) {
  stop("Verified input tree does not contain 495 unique tips", call. = FALSE)
}

trees <- setNames(
  lapply(variant_specs$tree_path, read.tree),
  variant_specs$variant_id
)
contexts <- setNames(
  Map(
    function(tree, variant_id) make_context(tree, variant_id),
    trees,
    variant_specs$variant_id
  ),
  variant_specs$variant_id
)

input_internal <- internal_node_table(input_context)
variant_internal <- lapply(contexts, internal_node_table)
input_cluster_keys <- sort(input_internal$cluster_key)

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
  calibrated = c(TRUE, TRUE, FALSE),
  lower_bound_ma = c(91.5046, 44.1968, NA_real_),
  upper_bound_ma = c(100.8925, 52.9473, NA_real_),
  interval_role = c(
    "treePL calibration",
    "treePL calibration",
    "uncalibrated focal comparison"
  ),
  stringsAsFactors = FALSE
)
required_anchors <- unique(c(node_specs$anchor_tip_1, node_specs$anchor_tip_2))

checks_list <- list()
add_check <- function(
    variant_id,
    check_id,
    passed,
    observed,
    expected,
    critical = TRUE,
    note = "") {
  passed <- isTRUE(passed)
  checks_list[[length(checks_list) + 1L]] <<- data.frame(
    variant_id = variant_id,
    check_id = check_id,
    critical = critical,
    status = if (passed) "PASS" else if (critical) "FAIL" else "WARN",
    observed = as.character(observed),
    expected = as.character(expected),
    note = as.character(note),
    stringsAsFactors = FALSE
  )
}

add_info <- function(variant_id, check_id, observed, note = "") {
  checks_list[[length(checks_list) + 1L]] <<- data.frame(
    variant_id = variant_id,
    check_id = check_id,
    critical = FALSE,
    status = "INFO",
    observed = as.character(observed),
    expected = "reported; no pass/fail threshold",
    note = as.character(note),
    stringsAsFactors = FALSE
  )
}

focal_rows <- list()
focal_by_variant <- list()

for (i in seq_len(nrow(variant_specs))) {
  spec <- variant_specs[i, , drop = FALSE]
  variant_id <- spec$variant_id[[1]]
  tree <- trees[[variant_id]]
  context <- contexts[[variant_id]]
  internal <- variant_internal[[variant_id]]

  same_tip_set <- length(tree$tip.label) == length(input_tree$tip.label) &&
    setequal(tree$tip.label, input_tree$tip.label)
  rf_ph85 <- if (same_tip_set && !anyDuplicated(tree$tip.label)) {
    suppressWarnings(as.numeric(dist.topo(input_tree, tree, method = "PH85")))
  } else {
    NA_real_
  }
  rooted_clusters_equal <- if (same_tip_set) {
    identical(sort(internal$cluster_key), input_cluster_keys)
  } else {
    FALSE
  }
  root_split <- sort(vapply(
    context$root_children,
    function(node) length(context$descendant_ids(node)),
    integer(1)
  ))
  root_tip_min <- min(context$tip_depth)
  root_tip_max <- max(context$tip_depth)
  root_tip_range <- root_tip_max - root_tip_min

  add_check(variant_id, "tip_count_495", context$n_tip == 495L, context$n_tip, 495L)
  add_check(
    variant_id,
    "unique_tip_labels",
    anyDuplicated(tree$tip.label) == 0L,
    anyDuplicated(tree$tip.label),
    0L
  )
  add_check(
    variant_id,
    "same_tip_set_as_input",
    same_tip_set,
    paste0(
      "input_only=", length(setdiff(input_tree$tip.label, tree$tip.label)),
      ";variant_only=", length(setdiff(tree$tip.label, input_tree$tip.label))
    ),
    "input_only=0;variant_only=0"
  )
  add_check(variant_id, "rooted", is.rooted(tree), is.rooted(tree), TRUE)
  add_check(variant_id, "binary", is.binary.phylo(tree), is.binary.phylo(tree), TRUE)
  add_check(
    variant_id,
    "branch_lengths_finite",
    all(is.finite(tree$edge.length)),
    sum(!is.finite(tree$edge.length)),
    0L
  )
  add_check(
    variant_id,
    "branch_lengths_nonnegative",
    all(tree$edge.length >= 0),
    sum(tree$edge.length < 0),
    0L
  )
  add_check(
    variant_id,
    "ape_PH85_RF_zero_vs_input",
    is.finite(rf_ph85) && rf_ph85 == 0,
    fmt_num(rf_ph85),
    0
  )
  add_check(
    variant_id,
    "rooted_cluster_sets_identical_to_input",
    rooted_clusters_equal,
    rooted_clusters_equal,
    TRUE,
    note = "Root-aware topology gate complements ape PH85 RF."
  )
  add_check(
    variant_id,
    "root_split_7_488",
    identical(root_split, c(7L, 488L)),
    paste(root_split, collapse = ","),
    "7,488"
  )
  add_check(
    variant_id,
    "root_to_tip_range_le_1e-5_ma",
    is.finite(root_tip_range) && root_tip_range <= ULTRAMETRIC_SERIALIZATION_TOL_MA,
    paste0(
      "min=", fmt_num(root_tip_min),
      ";max=", fmt_num(root_tip_max),
      ";range=", fmt_num(root_tip_range)
    ),
    paste0("range<=", fmt_num(ULTRAMETRIC_SERIALIZATION_TOL_MA), " Ma"),
    note = "Absolute tolerance accommodates six-decimal Newick serialization."
  )
  add_check(
    variant_id,
    "required_focal_anchors_exactly_once",
    all(label_counts(tree, required_anchors) == 1L),
    paste(label_counts(tree, required_anchors), collapse = ","),
    paste(rep(1L, length(required_anchors)), collapse = ",")
  )
  if (spec$variant_role[[1]] == "reproducibility") {
    add_check(
      variant_id, "byte_identical_to_primary",
      spec$byte_identical_to_primary[[1]],
      paste0("byte_identical_to_primary=", spec$byte_identical_to_primary[[1]]),
      "byte_identical_to_primary=TRUE",
      note = "Same-seed repeat is a reproducibility check, not a sensitivity variant."
    )
  }

  nodes <- vapply(seq_len(nrow(node_specs)), function(j) {
    safe_mrca(tree, c(node_specs$anchor_tip_1[[j]], node_specs$anchor_tip_2[[j]]))
  }, integer(1))
  descendants <- lapply(nodes, function(node) safe_descendants(context, node))
  ages <- t(vapply(nodes, function(node) node_age_stats(context, node), numeric(4)))
  colnames(ages) <- c("mean", "min", "max", "range")

  pap_complement <- sort(setdiff(tree$tip.label, papilionidae_expected))
  descendant_sets_ok <- c(
    identical(descendants[[1]], sort(tree$tip.label)),
    identical(descendants[[2]], papilionidae_expected),
    identical(descendants[[3]], pap_complement)
  )
  node_positions_ok <- c(
    !is.na(nodes[[1]]) && nodes[[1]] == context$root,
    !is.na(nodes[[2]]) && nodes[[2]] %in% context$root_children,
    !is.na(nodes[[3]]) && nodes[[3]] %in% context$root_children
  )

  for (j in seq_len(nrow(node_specs))) {
    node_id <- node_specs$node_id[[j]]
    desc_n <- length(descendants[[j]])
    add_check(
      variant_id,
      paste0(node_id, "_descendant_set_exact"),
      descendant_sets_ok[[j]] && desc_n == node_specs$expected_descendant_tip_count[[j]],
      paste0("descendants=", desc_n),
      paste0("exact_set;descendants=", node_specs$expected_descendant_tip_count[[j]])
    )
    add_check(
      variant_id,
      paste0(node_id, "_position_exact"),
      node_positions_ok[[j]],
      if (is.na(nodes[[j]])) "NA" else as.character(nodes[[j]]),
      if (j == 1L) as.character(context$root) else paste(context$root_children, collapse = ",")
    )
    add_check(
      variant_id,
      paste0(node_id, "_age_finite"),
      is.finite(ages[j, "mean"]),
      fmt_num(ages[j, "mean"]),
      "finite Ma"
    )

    within_bounds <- if (node_specs$calibrated[[j]]) {
      is.finite(ages[j, "mean"]) &&
        ages[j, "mean"] >= node_specs$lower_bound_ma[[j]] - CALIBRATION_BOUND_TOL_MA &&
        ages[j, "mean"] <= node_specs$upper_bound_ma[[j]] + CALIBRATION_BOUND_TOL_MA
    } else {
      NA
    }
    if (node_specs$calibrated[[j]]) {
      add_check(
        variant_id,
        paste0(node_id, "_within_calibration_bounds"),
        within_bounds,
        fmt_num(ages[j, "mean"]),
        paste0(
          fmt_num(node_specs$lower_bound_ma[[j]]), "-",
          fmt_num(node_specs$upper_bound_ma[[j]]), " Ma"
        ),
        note = paste0("Inclusive bounds with ", CALIBRATION_BOUND_TOL_MA, " Ma numeric tolerance.")
      )
    }

    focal_rows[[length(focal_rows) + 1L]] <- data.frame(
      variant_id = variant_id,
      variant_role = spec$variant_role[[1]],
      is_primary = spec$is_primary[[1]],
      byte_identical_to_primary = spec$byte_identical_to_primary[[1]],
      thorough = spec$thorough[[1]],
      smooth = spec$smooth[[1]],
      tree_path = spec$tree_path[[1]],
      node_id = node_id,
      anchor_tip_1 = node_specs$anchor_tip_1[[j]],
      anchor_tip_2 = node_specs$anchor_tip_2[[j]],
      mrca_node = nodes[[j]],
      node_position_exact = node_positions_ok[[j]],
      descendant_set_exact = descendant_sets_ok[[j]],
      descendant_tip_count = desc_n,
      node_age_ma = ages[j, "mean"],
      node_age_descendant_tip_min_ma = ages[j, "min"],
      node_age_descendant_tip_max_ma = ages[j, "max"],
      node_age_descendant_tip_range_ma = ages[j, "range"],
      calibrated = node_specs$calibrated[[j]],
      interval_role = node_specs$interval_role[[j]],
      lower_bound_ma = node_specs$lower_bound_ma[[j]],
      upper_bound_ma = node_specs$upper_bound_ma[[j]],
      within_calibration_bounds = within_bounds,
      distance_above_lower_ma = if (node_specs$calibrated[[j]]) {
        ages[j, "mean"] - node_specs$lower_bound_ma[[j]]
      } else {
        NA_real_
      },
      distance_below_upper_ma = if (node_specs$calibrated[[j]]) {
        node_specs$upper_bound_ma[[j]] - ages[j, "mean"]
      } else {
        NA_real_
      },
      lower_boundary_hit = if (node_specs$calibrated[[j]]) {
        abs(ages[j, "mean"] - node_specs$lower_bound_ma[[j]]) <= CALIBRATION_BOUND_TOL_MA
      } else {
        NA
      },
      upper_boundary_hit = if (node_specs$calibrated[[j]]) {
        abs(node_specs$upper_bound_ma[[j]] - ages[j, "mean"]) <= CALIBRATION_BOUND_TOL_MA
      } else {
        NA
      },
      root_to_tip_min_ma = root_tip_min,
      root_to_tip_max_ma = root_tip_max,
      root_to_tip_range_ma = root_tip_range,
      age_definition = paste0(
        "mean node-to-descendant-tip path length; min/max/range retained for ",
        "serialization diagnostics"
      ),
      stringsAsFactors = FALSE
    )
  }

  focal_by_variant[[variant_id]] <- ages[, "mean"]
}

node_ages <- do.call(rbind, focal_rows)
primary_id <- variant_specs$variant_id[variant_specs$is_primary][[1]]
primary_rows <- node_ages[node_ages$variant_id == primary_id, , drop = FALSE]
primary_age_by_node <- setNames(primary_rows$node_age_ma, primary_rows$node_id)
node_ages$primary_node_age_ma <- unname(primary_age_by_node[node_ages$node_id])
node_ages$age_difference_from_primary_ma <-
  node_ages$node_age_ma - node_ages$primary_node_age_ma
node_ages$absolute_age_difference_from_primary_ma <-
  abs(node_ages$age_difference_from_primary_ma)

primary_internal <- variant_internal[[primary_id]]
max_internal_summary <- data.frame(
  variant_id = variant_specs$variant_id,
  max_abs_all_internal_node_age_difference_from_primary_ma = NA_real_,
  signed_difference_at_max_ma = NA_real_,
  primary_node_at_max = NA_integer_,
  variant_node_at_max = NA_integer_,
  descendant_tip_count_at_max = NA_integer_,
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(variant_specs))) {
  variant_id <- variant_specs$variant_id[[i]]
  internal <- variant_internal[[variant_id]]
  idx <- match(primary_internal$cluster_key, internal$cluster_key)
  if (all(!is.na(idx)) && length(idx) == nrow(primary_internal)) {
    delta <- internal$age_ma[idx] - primary_internal$age_ma
    if (length(delta) && all(is.finite(delta))) {
      at_max <- which.max(abs(delta))
      max_internal_summary$max_abs_all_internal_node_age_difference_from_primary_ma[[i]] <-
        abs(delta[[at_max]])
      max_internal_summary$signed_difference_at_max_ma[[i]] <- delta[[at_max]]
      max_internal_summary$primary_node_at_max[[i]] <- primary_internal$node[[at_max]]
      max_internal_summary$variant_node_at_max[[i]] <- internal$node[idx[[at_max]]]
      max_internal_summary$descendant_tip_count_at_max[[i]] <-
        primary_internal$descendant_tip_count[[at_max]]
    }
  }
}

max_focal_by_variant <- vapply(
  variant_specs$variant_id,
  function(variant_id) {
    values <- node_ages$absolute_age_difference_from_primary_ma[
      node_ages$variant_id == variant_id
    ]
    if (!length(values) || any(!is.finite(values))) return(NA_real_)
    max(values)
  },
  numeric(1)
)
node_ages$max_abs_focal_node_age_difference_from_primary_ma <-
  unname(max_focal_by_variant[node_ages$variant_id])
node_ages$max_abs_all_internal_node_age_difference_from_primary_ma <-
  max_internal_summary$max_abs_all_internal_node_age_difference_from_primary_ma[
    match(node_ages$variant_id, max_internal_summary$variant_id)
  ]

for (i in seq_len(nrow(variant_specs))) {
  variant_id <- variant_specs$variant_id[[i]]
  summary_row <- max_internal_summary[
    max_internal_summary$variant_id == variant_id,
    ,
    drop = FALSE
  ]
  add_info(
    variant_id,
    "max_abs_focal_node_age_difference_from_primary_ma",
    fmt_num(max_focal_by_variant[[variant_id]]),
    "Maximum absolute difference across root, seven-tip Papilionidae, and 488-tip Hesperiinae."
  )
  add_info(
    variant_id,
    "max_abs_all_internal_node_age_difference_from_primary_ma",
    fmt_num(summary_row$max_abs_all_internal_node_age_difference_from_primary_ma[[1]]),
    paste0(
      "Compared by exact rooted descendant-tip set across all internal nodes; signed difference=",
      fmt_num(summary_row$signed_difference_at_max_ma[[1]]),
      "; primary_node=", summary_row$primary_node_at_max[[1]],
      "; variant_node=", summary_row$variant_node_at_max[[1]],
      "; descendants=", summary_row$descendant_tip_count_at_max[[1]], "."
    )
  )
}

sensitivity_variant <- variant_specs$variant_role == "smoothing_sensitivity"
global_internal_values <-
  max_internal_summary$max_abs_all_internal_node_age_difference_from_primary_ma[sensitivity_variant]
global_focal_values <- max_focal_by_variant[variant_specs$variant_id[sensitivity_variant]]
if (length(global_internal_values) && all(is.finite(global_internal_values))) {
  winner <- which.max(global_internal_values)
  winner_id <- variant_specs$variant_id[sensitivity_variant][[winner]]
  add_info(
    "ALL_SENSITIVITY_VARIANTS",
    "largest_all_internal_node_age_difference_from_primary_ma",
    fmt_num(global_internal_values[[winner]]),
    paste0("Largest value occurs in ", winner_id, "; reported without a sensitivity threshold.")
  )
}
if (length(global_focal_values) && all(is.finite(global_focal_values))) {
  winner <- which.max(global_focal_values)
  winner_id <- names(global_focal_values)[[winner]]
  add_info(
    "ALL_SENSITIVITY_VARIANTS",
    "largest_focal_node_age_difference_from_primary_ma",
    fmt_num(global_focal_values[[winner]]),
    paste0("Largest value occurs in ", winner_id, "; reported without a sensitivity threshold.")
  )
}

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
if (any(critical_failures)) {
  failed <- paste(
    paste(checks$variant_id[critical_failures], checks$check_id[critical_failures], sep = ":"),
    collapse = ", "
  )
  stop(
    "Final-variant QA hard gate failed: ", failed,
    ". Diagnostic outputs were retained at ", checks_path,
    call. = FALSE
  )
}

cat("Final-variant QA PASS\n")
cat("Focal node ages: ", node_ages_path, "\n", sep = "")
cat("Checks and sensitivity summaries: ", checks_path, "\n", sep = "")
