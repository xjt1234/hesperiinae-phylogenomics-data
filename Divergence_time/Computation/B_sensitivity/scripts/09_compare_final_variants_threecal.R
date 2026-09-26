#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ape))

options(stringsAsFactors = FALSE, digits = 15)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop(
    "Usage: 09_compare_final_variants_threecal.R <run_root>",
    call. = FALSE
  )
}

run_root <- normalizePath(args[[1]], mustWork = TRUE)
qa_dir <- file.path(run_root, "qa")
manifest_path <- file.path(qa_dir, "final_variant_manifest.tsv")
if (!file.exists(manifest_path)) {
  stop("Missing final-variant manifest: ", manifest_path, call. = FALSE)
}
decision_path <- file.path(qa_dir, "cv_23point_decision.tsv")
if (!file.exists(decision_path)) {
  stop("Missing qualified 23-point CV decision: ", decision_path, call. = FALSE)
}

CALIBRATION_BOUND_TOL_MA <- 1e-5
ULTRAMETRIC_SERIALIZATION_TOL_MA <- 1e-5

output_paths <- c(
  checks = file.path(qa_dir, "final_variant_threecal_checks.tsv"),
  calibration_ages = file.path(qa_dir, "final_variant_threecal_calibration_ages.tsv"),
  internal_differences = file.path(qa_dir, "final_variant_threecal_internal_node_age_differences.tsv"),
  difference_summary = file.path(qa_dir, "final_variant_threecal_age_difference_summary.tsv"),
  report = file.path(qa_dir, "final_variant_threecal_report.md")
)
existing_outputs <- output_paths[file.exists(output_paths)]
if (length(existing_outputs)) {
  stop(
    "Refusing to overwrite existing final-variant QA output(s): ",
    paste(existing_outputs, collapse = ", "),
    call. = FALSE
  )
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
  strip.white = FALSE
)

required_columns <- c(
  "variant_id", "variant_role", "smooth", "thorough", "seed",
  "config_path", "tree_path", "stage"
)
missing_columns <- setdiff(required_columns, names(manifest))
if (length(missing_columns)) {
  stop(
    "Manifest is missing required column(s): ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}
manifest <- manifest[, required_columns, drop = FALSE]
if (!nrow(manifest)) stop("Manifest has no variant rows", call. = FALSE)
if (anyNA(manifest) || any(vapply(manifest, function(x) any(!nzchar(x)), logical(1)))) {
  stop("Manifest required fields must be non-empty and non-NA", call. = FALSE)
}
if (anyDuplicated(manifest$variant_id)) {
  stop("Manifest variant_id values must be unique", call. = FALSE)
}
if (anyDuplicated(manifest$stage)) {
  stop("Manifest stage values must be unique", call. = FALSE)
}
safe_id <- grepl("^[A-Za-z0-9_.-]+$", manifest$variant_id)
safe_stage <- grepl("^[A-Za-z0-9_.-]+$", manifest$stage)
if (!all(safe_id)) {
  stop("Unsafe variant_id value(s): ", paste(manifest$variant_id[!safe_id], collapse = ", "), call. = FALSE)
}
if (!all(safe_stage)) {
  stop("Unsafe stage value(s): ", paste(manifest$stage[!safe_stage], collapse = ", "), call. = FALSE)
}
allowed_variant_roles <- c(
  "primary",
  "optimization_diagnostic",
  "same_seed_repeat",
  "smoothing_sensitivity"
)
unknown_variant_roles <- setdiff(unique(manifest$variant_role), allowed_variant_roles)
if (length(unknown_variant_roles)) {
  stop(
    "Manifest contains unsupported variant_role value(s): ",
    paste(unknown_variant_roles, collapse = ", "),
    call. = FALSE
  )
}
if (sum(manifest$variant_role == "primary") != 1L) {
  stop("Manifest must contain exactly one variant_role=primary row", call. = FALSE)
}
if (sum(manifest$variant_role == "optimization_diagnostic") != 1L) {
  stop(
    "Manifest must contain exactly one variant_role=optimization_diagnostic row",
    call. = FALSE
  )
}
if (sum(manifest$variant_role == "same_seed_repeat") < 1L) {
  stop("Manifest must contain at least one variant_role=same_seed_repeat row", call. = FALSE)
}

parse_exact_bool <- function(x, field) {
  if (!all(x %in% c("TRUE", "FALSE"))) {
    stop(field, " must use exact uppercase TRUE or FALSE values", call. = FALSE)
  }
  x == "TRUE"
}

manifest$thorough_value <- parse_exact_bool(manifest$thorough, "Manifest thorough")
manifest$smooth_value <- suppressWarnings(as.numeric(manifest$smooth))
if (any(!is.finite(manifest$smooth_value) | manifest$smooth_value <= 0)) {
  stop("Manifest smooth values must be finite positive numbers", call. = FALSE)
}
if (any(!grepl("^[0-9]+$", manifest$seed))) {
  stop("Manifest seed values must be unsigned base-10 integers", call. = FALSE)
}
manifest$seed_value <- suppressWarnings(as.numeric(manifest$seed))
if (any(!is.finite(manifest$seed_value) | manifest$seed_value < 0)) {
  stop("Manifest seed values are outside the supported numeric range", call. = FALSE)
}

decision <- read.delim(
  decision_path,
  sep = "\t",
  header = TRUE,
  quote = "",
  comment.char = "",
  colClasses = "character",
  check.names = FALSE,
  na.strings = character(),
  strip.white = FALSE
)
if (nrow(decision) != 1L) {
  stop("Expected exactly one row in ", decision_path, call. = FALSE)
}
if (anyDuplicated(names(decision))) {
  stop("Duplicate column name(s) in ", decision_path, call. = FALSE)
}
required_decision_columns <- c(
  "sensitivity_required",
  "sensitivity_smoothing_values"
)
missing_decision_columns <- setdiff(required_decision_columns, names(decision))
if (length(missing_decision_columns)) {
  stop(
    "23-point CV decision is missing required column(s): ",
    paste(missing_decision_columns, collapse = ", "),
    call. = FALSE
  )
}

decision_value <- function(key) {
  value <- decision[[key]][[1L]]
  if (is.na(value)) value <- ""
  if (!identical(value, trimws(value))) {
    stop("23-point CV decision field has surrounding whitespace: ", key, call. = FALSE)
  }
  value
}

sensitivity_required_text <- decision_value("sensitivity_required")
if (!sensitivity_required_text %in% c("TRUE", "FALSE")) {
  stop(
    "23-point CV decision sensitivity_required must be exact TRUE or FALSE",
    call. = FALSE
  )
}
sensitivity_required <- identical(sensitivity_required_text, "TRUE")
sensitivity_field <- decision_value("sensitivity_smoothing_values")
decision_sensitivity_tokens <- character()

if (sensitivity_required) {
  if (!nzchar(sensitivity_field)) {
    stop(
      "23-point CV decision requires sensitivity runs but supplies no smoothing values",
      call. = FALSE
    )
  }
  raw_tokens <- strsplit(sensitivity_field, ",", fixed = TRUE)[[1L]]
  decision_sensitivity_tokens <- trimws(raw_tokens)
  if (any(!nzchar(decision_sensitivity_tokens)) ||
      !identical(decision_sensitivity_tokens, raw_tokens)) {
    stop(
      "23-point sensitivity_smoothing_values contains whitespace or an empty token",
      call. = FALSE
    )
  }
  decision_sensitivity_numeric <- suppressWarnings(as.numeric(decision_sensitivity_tokens))
  if (any(!is.finite(decision_sensitivity_numeric) | decision_sensitivity_numeric <= 0)) {
    stop(
      "23-point sensitivity_smoothing_values must contain only finite positive numbers",
      call. = FALSE
    )
  }
  if (anyDuplicated(decision_sensitivity_tokens) ||
      anyDuplicated(decision_sensitivity_numeric)) {
    stop(
      "23-point sensitivity_smoothing_values must be a deduplicated set",
      call. = FALSE
    )
  }
} else if (nzchar(sensitivity_field)) {
  stop(
    "23-point CV decision has sensitivity_required=FALSE but a non-empty ",
    "sensitivity_smoothing_values field",
    call. = FALSE
  )
}

manifest_sensitivity_tokens <- manifest$smooth[
  manifest$variant_role == "smoothing_sensitivity"
]
if (anyDuplicated(manifest_sensitivity_tokens)) {
  stop("Manifest contains duplicate smoothing_sensitivity smooth values", call. = FALSE)
}
if (!sensitivity_required && length(manifest_sensitivity_tokens) != 0L) {
  stop(
    "Manifest must contain zero smoothing_sensitivity rows when ",
    "sensitivity_required=FALSE",
    call. = FALSE
  )
}
if (sensitivity_required &&
    (!identical(
      sort(manifest_sensitivity_tokens),
      sort(decision_sensitivity_tokens)
    ) || length(manifest_sensitivity_tokens) != length(decision_sensitivity_tokens))) {
  stop(
    "Manifest smoothing_sensitivity smooth set does not exactly match the ",
    "23-point decision: manifest={",
    paste(sort(manifest_sensitivity_tokens), collapse = ","),
    "}; decision={",
    paste(sort(decision_sensitivity_tokens), collapse = ","),
    "}",
    call. = FALSE
  )
}

is_absolute_path <- function(path) startsWith(path, .Platform$file.sep)

resolve_run_file <- function(path, field, variant_id) {
  candidate <- if (is_absolute_path(path)) path else file.path(run_root, path)
  resolved <- normalizePath(candidate, mustWork = TRUE)
  run_prefix <- paste0(run_root, .Platform$file.sep)
  if (!startsWith(resolved, run_prefix) || isTRUE(file.info(resolved)$isdir)) {
    stop(
      field, " for ", variant_id,
      " must resolve to a file inside the current run root: ", resolved,
      call. = FALSE
    )
  }
  resolved
}

manifest$config_path_resolved <- Map(
  resolve_run_file,
  manifest$config_path,
  MoreArgs = list(field = "config_path"),
  variant_id = manifest$variant_id
) |>
  unlist(use.names = FALSE)
manifest$tree_path_resolved <- Map(
  resolve_run_file,
  manifest$tree_path,
  MoreArgs = list(field = "tree_path"),
  variant_id = manifest$variant_id
) |>
  unlist(use.names = FALSE)
if (anyDuplicated(manifest$config_path_resolved)) {
  stop("Manifest config_path values must resolve to unique files", call. = FALSE)
}
if (anyDuplicated(manifest$tree_path_resolved)) {
  stop("Manifest tree_path values must resolve to unique files", call. = FALSE)
}

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
    variant_id = as.character(variant_id),
    check_id = as.character(check_id),
    critical = isTRUE(critical),
    status = if (passed) "PASS" else if (critical) "FAIL" else "WARN",
    observed = as.character(observed),
    expected = as.character(expected),
    note = as.character(note),
    stringsAsFactors = FALSE
  )
}

add_info <- function(variant_id, check_id, observed, note = "") {
  checks_list[[length(checks_list) + 1L]] <<- data.frame(
    variant_id = as.character(variant_id),
    check_id = as.character(check_id),
    critical = FALSE,
    status = "INFO",
    observed = as.character(observed),
    expected = "reported; no pass/fail threshold",
    note = as.character(note),
    stringsAsFactors = FALSE
  )
}

fmt_num <- function(x) {
  if (length(x) != 1L || is.na(x) || !is.finite(x)) return("NA")
  sprintf("%.12g", x)
}

read_active_config <- function(path) {
  lines <- readLines(path, warn = FALSE)
  active <- trimws(sub("#.*$", "", lines))
  active[nzchar(active)]
}

config_values <- function(active, key) {
  pattern <- paste0("^", key, "[[:space:]]*=")
  hits <- active[grepl(pattern, active)]
  trimws(sub(pattern, "", hits))
}

resolve_config_output <- function(path) {
  if (length(path) != 1L || !nzchar(path)) return(NA_character_)
  candidate <- if (is_absolute_path(path)) path else file.path(run_root, path)
  normalizePath(candidate, mustWork = FALSE)
}

required_calibration_lines <- c(
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

read_stage_metadata <- function(path) {
  if (!file.exists(path)) return(setNames(character(), character()))
  values <- read.delim(
    path,
    sep = "\t",
    header = FALSE,
    quote = "",
    comment.char = "",
    fill = TRUE,
    stringsAsFactors = FALSE,
    col.names = c("key", "value")
  )
  if (!nrow(values)) return(setNames(character(), character()))
  setNames(values$value, values$key)
}

metadata_value <- function(metadata, key) {
  values <- unname(metadata[names(metadata) == key])
  if (length(values) != 1L) return(NA_character_)
  values[[1]]
}

for (i in seq_len(nrow(manifest))) {
  variant_id <- manifest$variant_id[[i]]
  active <- read_active_config(manifest$config_path_resolved[[i]])
  smooth_cfg <- config_values(active, "smooth")
  seed_cfg <- config_values(active, "seed")
  outfile_cfg <- config_values(active, "outfile")
  cfg_smooth_value <- suppressWarnings(as.numeric(smooth_cfg))
  cfg_seed_value <- suppressWarnings(as.numeric(seed_cfg))
  thorough_cfg <- sum(active == "thorough") == 1L

  add_check(
    variant_id,
    "config_exact_three_calibration_blocks",
    all(vapply(required_calibration_lines, function(x) sum(active == x) == 1L, logical(1))) &&
      sum(grepl("^mrca[[:space:]]*=", active)) == 3L &&
      sum(grepl("^min[[:space:]]*=", active)) == 3L &&
      sum(grepl("^max[[:space:]]*=", active)) == 3L,
    paste0(
      "required_exact=", sum(vapply(required_calibration_lines, function(x) sum(active == x) == 1L, logical(1))),
      "/9;mrca=", sum(grepl("^mrca[[:space:]]*=", active)),
      ";min=", sum(grepl("^min[[:space:]]*=", active)),
      ";max=", sum(grepl("^max[[:space:]]*=", active))
    ),
    "required_exact=9/9;mrca=3;min=3;max=3"
  )
  add_check(
    variant_id,
    "config_smooth_matches_manifest",
    length(cfg_smooth_value) == 1L && is.finite(cfg_smooth_value) &&
      identical(cfg_smooth_value, manifest$smooth_value[[i]]),
    if (length(smooth_cfg) == 1L) smooth_cfg else paste(smooth_cfg, collapse = ","),
    manifest$smooth[[i]]
  )
  add_check(
    variant_id,
    "config_seed_matches_manifest",
    length(cfg_seed_value) == 1L && is.finite(cfg_seed_value) &&
      identical(cfg_seed_value, manifest$seed_value[[i]]),
    if (length(seed_cfg) == 1L) seed_cfg else paste(seed_cfg, collapse = ","),
    manifest$seed[[i]]
  )
  add_check(
    variant_id,
    "config_thorough_matches_manifest",
    identical(thorough_cfg, manifest$thorough_value[[i]]) && sum(active == "thorough") <= 1L,
    paste0("thorough_directive_count=", sum(active == "thorough")),
    paste0("thorough=", manifest$thorough[[i]])
  )
  outfile_resolved <- if (length(outfile_cfg) == 1L) resolve_config_output(outfile_cfg) else NA_character_
  add_check(
    variant_id,
    "config_outfile_matches_manifest_tree",
    length(outfile_cfg) == 1L && identical(outfile_resolved, manifest$tree_path_resolved[[i]]),
    if (length(outfile_cfg) == 1L) outfile_resolved else paste(outfile_cfg, collapse = ","),
    manifest$tree_path_resolved[[i]]
  )

  metadata_path <- file.path(run_root, "logs", paste0(manifest$stage[[i]], ".metadata.tsv"))
  metadata <- read_stage_metadata(metadata_path)
  metadata_config <- metadata_value(metadata, "config")
  metadata_config_resolved <- if (!is.na(metadata_config) && file.exists(metadata_config)) {
    normalizePath(metadata_config, mustWork = TRUE)
  } else {
    metadata_config
  }
  add_check(
    variant_id,
    "stage_metadata_exists",
    file.exists(metadata_path),
    file.exists(metadata_path),
    TRUE
  )
  add_check(
    variant_id,
    "stage_metadata_stage_matches_manifest",
    identical(metadata_value(metadata, "stage"), manifest$stage[[i]]),
    metadata_value(metadata, "stage"),
    manifest$stage[[i]]
  )
  add_check(
    variant_id,
    "stage_metadata_config_matches_manifest",
    identical(metadata_config_resolved, manifest$config_path_resolved[[i]]),
    metadata_config_resolved,
    manifest$config_path_resolved[[i]]
  )
  add_check(
    variant_id,
    "stage_metadata_exit_code_zero",
    identical(metadata_value(metadata, "exit_code"), "0"),
    metadata_value(metadata, "exit_code"),
    0
  )
}

files_byte_identical <- function(path_a, path_b) {
  info_a <- file.info(path_a)
  info_b <- file.info(path_b)
  if (is.na(info_a$size) || is.na(info_b$size) || info_a$size != info_b$size) return(FALSE)
  identical(
    readBin(path_a, what = "raw", n = info_a$size),
    readBin(path_b, what = "raw", n = info_b$size)
  )
}

make_context <- function(tree, tree_name) {
  root <- setdiff(unique(tree$edge[, 1]), unique(tree$edge[, 2]))
  if (length(root) != 1L) {
    stop("Could not identify exactly one root in ", tree_name, call. = FALSE)
  }
  root <- as.integer(root[[1]])
  n_tip <- Ntip(tree)
  cache <- new.env(parent = emptyenv())
  descendant_ids <- function(node) {
    key <- as.character(node)
    if (exists(key, envir = cache, inherits = FALSE)) return(get(key, envir = cache))
    result <- if (node <= n_tip) {
      as.integer(node)
    } else {
      children <- tree$edge[tree$edge[, 1] == node, 2]
      as.integer(unlist(lapply(children, descendant_ids), use.names = FALSE))
    }
    assign(key, result, envir = cache)
    result
  }
  descendant_labels <- function(node) sort(unique(tree$tip.label[descendant_ids(node)]))
  depth <- if (is.null(tree$edge.length)) {
    rep(NA_real_, n_tip + Nnode(tree))
  } else {
    node.depth.edgelength(tree)
  }
  list(
    tree = tree,
    name = tree_name,
    n_tip = n_tip,
    root = root,
    root_children = as.integer(tree$edge[tree$edge[, 1] == root, 2]),
    depth = depth,
    descendant_ids = descendant_ids,
    descendant_labels = descendant_labels
  )
}

node_age_stats <- function(context, node) {
  if (length(node) != 1L || is.na(node)) {
    return(c(mean = NA_real_, min = NA_real_, max = NA_real_, range = NA_real_))
  }
  tips <- context$descendant_ids(node)
  ages <- context$depth[tips] - context$depth[node]
  if (!length(ages) || any(!is.finite(ages))) {
    return(c(mean = NA_real_, min = NA_real_, max = NA_real_, range = NA_real_))
  }
  c(mean = mean(ages), min = min(ages), max = max(ages), range = diff(range(ages)))
}

internal_node_table <- function(context) {
  nodes <- sort(unique(as.integer(context$tree$edge[, 1])))
  do.call(rbind, lapply(nodes, function(node) {
    labels <- context$descendant_labels(node)
    ages <- node_age_stats(context, node)
    data.frame(
      node = node,
      descendant_tip_count = length(labels),
      clade_tip_key = paste(labels, collapse = "|"),
      age_ma = unname(ages[["mean"]]),
      age_tip_min_ma = unname(ages[["min"]]),
      age_tip_max_ma = unname(ages[["max"]]),
      age_tip_range_ma = unname(ages[["range"]]),
      stringsAsFactors = FALSE
    )
  }))
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

trees <- setNames(
  lapply(manifest$tree_path_resolved, function(path) {
    tree <- read.tree(path)
    if (!inherits(tree, "phylo")) stop("Expected exactly one phylo tree in ", path, call. = FALSE)
    tree
  }),
  manifest$variant_id
)
contexts <- setNames(
  Map(function(tree, id) make_context(tree, id), trees, manifest$variant_id),
  manifest$variant_id
)
internal_tables <- lapply(contexts, internal_node_table)

primary_index <- which(manifest$variant_role == "primary")
primary_id <- manifest$variant_id[[primary_index]]
primary_tree <- trees[[primary_id]]
primary_context <- contexts[[primary_id]]
primary_internal <- internal_tables[[primary_id]]
primary_tree_path <- manifest$tree_path_resolved[[primary_index]]
primary_tip_set <- sort(primary_tree$tip.label)
primary_clade_set <- sort(primary_internal$clade_tip_key)

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
  lower_bound_ma = c(91.5046, 44.1968, 36.210861),
  upper_bound_ma = c(100.8925, 52.9473, 40.662537),
  stringsAsFactors = FALSE
)
required_anchors <- unique(c(node_specs$anchor_tip_1, node_specs$anchor_tip_2))

calibration_rows <- list()
internal_difference_rows <- list()
summary_rows <- list()

for (i in seq_len(nrow(manifest))) {
  spec <- manifest[i, , drop = FALSE]
  variant_id <- spec$variant_id[[1]]
  tree <- trees[[variant_id]]
  context <- contexts[[variant_id]]
  internal <- internal_tables[[variant_id]]

  unique_tips <- anyDuplicated(tree$tip.label) == 0L
  same_tip_set <- unique_tips && length(tree$tip.label) == length(primary_tree$tip.label) &&
    identical(sort(tree$tip.label), primary_tip_set)
  rf_ph85 <- if (same_tip_set) {
    tryCatch(
      suppressWarnings(as.numeric(dist.topo(primary_tree, tree, method = "PH85"))),
      error = function(e) NA_real_
    )
  } else {
    NA_real_
  }
  rooted_clades_equal <- same_tip_set && identical(sort(internal$clade_tip_key), primary_clade_set)
  root_tip_range <- unname(node_age_stats(context, context$root)[["range"]])

  add_check(variant_id, "tip_count_495", Ntip(tree) == 495L, Ntip(tree), 495L)
  add_check(variant_id, "unique_tip_labels", unique_tips, anyDuplicated(tree$tip.label), 0L)
  add_check(
    variant_id,
    "same_tip_labels_as_primary",
    same_tip_set,
    paste0(
      "primary_only=", length(setdiff(primary_tree$tip.label, tree$tip.label)),
      ";variant_only=", length(setdiff(tree$tip.label, primary_tree$tip.label))
    ),
    "primary_only=0;variant_only=0"
  )
  add_check(variant_id, "rooted", is.rooted(tree), is.rooted(tree), TRUE)
  add_check(variant_id, "binary", is.binary.phylo(tree), is.binary.phylo(tree), TRUE)
  add_check(
    variant_id,
    "branch_lengths_present_finite",
    !is.null(tree$edge.length) && all(is.finite(tree$edge.length)),
    if (is.null(tree$edge.length)) "missing" else sum(!is.finite(tree$edge.length)),
    "present;nonfinite=0"
  )
  add_check(
    variant_id,
    "branch_lengths_nonnegative",
    !is.null(tree$edge.length) && all(tree$edge.length >= 0),
    if (is.null(tree$edge.length)) "missing" else sum(tree$edge.length < 0),
    0L
  )
  add_check(
    variant_id,
    "ape_PH85_RF_zero_vs_primary",
    is.finite(rf_ph85) && rf_ph85 == 0,
    fmt_num(rf_ph85),
    0,
    note = "Unrooted RF gate is paired with an exact rooted-clade-set gate."
  )
  add_check(
    variant_id,
    "rooted_clade_sets_identical_to_primary",
    rooted_clades_equal,
    rooted_clades_equal,
    TRUE,
    note = "Every internal node is keyed by its exact sorted descendant-tip set."
  )
  add_check(
    variant_id,
    "root_to_tip_range_le_1e-5_ma",
    is.finite(root_tip_range) && root_tip_range <= ULTRAMETRIC_SERIALIZATION_TOL_MA,
    fmt_num(root_tip_range),
    paste0("<=", fmt_num(ULTRAMETRIC_SERIALIZATION_TOL_MA), " Ma"),
    note = "Tolerance accommodates six-decimal Newick serialization."
  )
  add_check(
    variant_id,
    "calibration_anchors_exactly_once",
    all(label_counts(tree, required_anchors) == 1L),
    paste(label_counts(tree, required_anchors), collapse = ","),
    paste(rep(1L, length(required_anchors)), collapse = ",")
  )

  byte_identical <- files_byte_identical(spec$tree_path_resolved[[1]], primary_tree_path)
  if (spec$variant_role[[1]] == "same_seed_repeat") {
    add_check(
      variant_id,
      "same_seed_repeat_seed_matches_primary",
      identical(spec$seed_value[[1]], manifest$seed_value[[primary_index]]),
      spec$seed[[1]],
      manifest$seed[[primary_index]]
    )
    add_check(
      variant_id,
      "same_seed_repeat_smooth_matches_primary",
      identical(spec$smooth_value[[1]], manifest$smooth_value[[primary_index]]),
      spec$smooth[[1]],
      manifest$smooth[[primary_index]]
    )
    add_check(
      variant_id,
      "same_seed_repeat_thorough_matches_primary",
      identical(spec$thorough_value[[1]], manifest$thorough_value[[primary_index]]),
      spec$thorough[[1]],
      manifest$thorough[[primary_index]]
    )
    add_check(
      variant_id,
      "same_seed_repeat_tree_byte_identical_to_primary",
      byte_identical,
      byte_identical,
      TRUE,
      note = "Exact raw-byte comparison of the complete Newick files."
    )
  } else {
    add_info(
      variant_id,
      "tree_byte_identical_to_primary",
      byte_identical,
      "Informational except for variant_role=same_seed_repeat."
    )
  }

  nodes <- vapply(seq_len(nrow(node_specs)), function(j) {
    safe_mrca(tree, c(node_specs$anchor_tip_1[[j]], node_specs$anchor_tip_2[[j]]))
  }, integer(1))
  descendants <- lapply(nodes, function(node) safe_descendants(context, node))
  age_stats <- t(vapply(nodes, function(node) node_age_stats(context, node), numeric(4)))
  colnames(age_stats) <- c("mean", "min", "max", "range")

  expected_descendants <- list(
    sort(tree$tip.label),
    papilionidae_expected,
    sort(setdiff(tree$tip.label, papilionidae_expected))
  )
  position_exact <- c(
    !is.na(nodes[[1]]) && nodes[[1]] == context$root,
    !is.na(nodes[[2]]) && nodes[[2]] %in% context$root_children,
    !is.na(nodes[[3]]) && nodes[[3]] %in% context$root_children
  )

  for (j in seq_len(nrow(node_specs))) {
    node_id <- node_specs$node_id[[j]]
    descendants_exact <- identical(descendants[[j]], expected_descendants[[j]]) &&
      length(descendants[[j]]) == node_specs$expected_descendant_tip_count[[j]]
    age <- unname(age_stats[j, "mean"])
    within_bounds <- is.finite(age) &&
      age >= node_specs$lower_bound_ma[[j]] - CALIBRATION_BOUND_TOL_MA &&
      age <= node_specs$upper_bound_ma[[j]] + CALIBRATION_BOUND_TOL_MA

    add_check(
      variant_id,
      paste0(node_id, "_descendant_set_exact"),
      descendants_exact,
      paste0("descendants=", length(descendants[[j]])),
      paste0("exact_set;descendants=", node_specs$expected_descendant_tip_count[[j]])
    )
    add_check(
      variant_id,
      paste0(node_id, "_position_exact"),
      position_exact[[j]],
      if (is.na(nodes[[j]])) "NA" else nodes[[j]],
      if (j == 1L) context$root else paste(context$root_children, collapse = ",")
    )
    add_check(
      variant_id,
      paste0(node_id, "_within_exact_calibration_interval"),
      within_bounds,
      fmt_num(age),
      paste0(
        fmt_num(node_specs$lower_bound_ma[[j]]), "-",
        fmt_num(node_specs$upper_bound_ma[[j]]), " Ma"
      ),
      note = paste0("Inclusive bounds with ", fmt_num(CALIBRATION_BOUND_TOL_MA), " Ma tolerance.")
    )

    calibration_rows[[length(calibration_rows) + 1L]] <- data.frame(
      variant_id = variant_id,
      variant_role = spec$variant_role[[1]],
      smooth = spec$smooth[[1]],
      thorough = spec$thorough_value[[1]],
      seed = spec$seed[[1]],
      stage = spec$stage[[1]],
      node_id = node_id,
      anchor_tip_1 = node_specs$anchor_tip_1[[j]],
      anchor_tip_2 = node_specs$anchor_tip_2[[j]],
      mrca_node = nodes[[j]],
      descendant_tip_count = length(descendants[[j]]),
      descendant_set_exact = descendants_exact,
      node_position_exact = position_exact[[j]],
      node_age_ma = age,
      node_age_descendant_tip_min_ma = unname(age_stats[j, "min"]),
      node_age_descendant_tip_max_ma = unname(age_stats[j, "max"]),
      node_age_descendant_tip_range_ma = unname(age_stats[j, "range"]),
      lower_bound_ma = node_specs$lower_bound_ma[[j]],
      upper_bound_ma = node_specs$upper_bound_ma[[j]],
      bound_tolerance_ma = CALIBRATION_BOUND_TOL_MA,
      within_calibration_bounds = within_bounds,
      distance_above_lower_ma = age - node_specs$lower_bound_ma[[j]],
      distance_below_upper_ma = node_specs$upper_bound_ma[[j]] - age,
      tree_path = spec$tree_path_resolved[[1]],
      stringsAsFactors = FALSE
    )
  }

  matched_index <- match(primary_internal$clade_tip_key, internal$clade_tip_key)
  matched <- which(!is.na(matched_index))
  differences <- if (length(matched)) {
    variant_rows <- matched_index[matched]
    data.frame(
      variant_id = variant_id,
      variant_role = spec$variant_role[[1]],
      primary_node = primary_internal$node[matched],
      variant_node = internal$node[variant_rows],
      descendant_tip_count = primary_internal$descendant_tip_count[matched],
      clade_tip_key = primary_internal$clade_tip_key[matched],
      primary_age_ma = primary_internal$age_ma[matched],
      variant_age_ma = internal$age_ma[variant_rows],
      signed_age_difference_from_primary_ma = internal$age_ma[variant_rows] - primary_internal$age_ma[matched],
      absolute_age_difference_from_primary_ma = abs(internal$age_ma[variant_rows] - primary_internal$age_ma[matched]),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      variant_id = character(), variant_role = character(), primary_node = integer(),
      variant_node = integer(), descendant_tip_count = integer(), clade_tip_key = character(),
      primary_age_ma = numeric(), variant_age_ma = numeric(),
      signed_age_difference_from_primary_ma = numeric(),
      absolute_age_difference_from_primary_ma = numeric(),
      stringsAsFactors = FALSE
    )
  }
  internal_difference_rows[[length(internal_difference_rows) + 1L]] <- differences

  finite_delta <- length(differences$absolute_age_difference_from_primary_ma) &&
    all(is.finite(differences$absolute_age_difference_from_primary_ma))
  at_max <- if (finite_delta) which.max(differences$absolute_age_difference_from_primary_ma) else NA_integer_
  variant_calibration_ages <- age_stats[, "mean"]
  primary_nodes <- vapply(seq_len(nrow(node_specs)), function(j) {
    safe_mrca(primary_tree, c(node_specs$anchor_tip_1[[j]], node_specs$anchor_tip_2[[j]]))
  }, integer(1))
  primary_calibration_ages <- vapply(
    primary_nodes,
    function(node) unname(node_age_stats(primary_context, node)[["mean"]]),
    numeric(1)
  )
  calibration_delta <- variant_calibration_ages - primary_calibration_ages

  summary_rows[[length(summary_rows) + 1L]] <- data.frame(
    variant_id = variant_id,
    variant_role = spec$variant_role[[1]],
    smooth = spec$smooth[[1]],
    thorough = spec$thorough_value[[1]],
    seed = spec$seed[[1]],
    stage = spec$stage[[1]],
    byte_identical_to_primary = byte_identical,
    primary_internal_node_count = nrow(primary_internal),
    variant_internal_node_count = nrow(internal),
    matched_rooted_internal_node_count = nrow(differences),
    max_abs_all_internal_node_age_difference_from_primary_ma = if (finite_delta) {
      differences$absolute_age_difference_from_primary_ma[[at_max]]
    } else NA_real_,
    mean_abs_all_internal_node_age_difference_from_primary_ma = if (finite_delta) {
      mean(differences$absolute_age_difference_from_primary_ma)
    } else NA_real_,
    median_abs_all_internal_node_age_difference_from_primary_ma = if (finite_delta) {
      median(differences$absolute_age_difference_from_primary_ma)
    } else NA_real_,
    signed_difference_at_max_ma = if (finite_delta) {
      differences$signed_age_difference_from_primary_ma[[at_max]]
    } else NA_real_,
    primary_node_at_max = if (finite_delta) differences$primary_node[[at_max]] else NA_integer_,
    variant_node_at_max = if (finite_delta) differences$variant_node[[at_max]] else NA_integer_,
    descendant_tip_count_at_max = if (finite_delta) {
      differences$descendant_tip_count[[at_max]]
    } else NA_integer_,
    max_abs_calibration_node_age_difference_from_primary_ma = if (all(is.finite(calibration_delta))) {
      max(abs(calibration_delta))
    } else NA_real_,
    all_three_calibration_nodes_within_bounds = all(vapply(seq_len(nrow(node_specs)), function(j) {
      is.finite(variant_calibration_ages[[j]]) &&
        variant_calibration_ages[[j]] >= node_specs$lower_bound_ma[[j]] - CALIBRATION_BOUND_TOL_MA &&
        variant_calibration_ages[[j]] <= node_specs$upper_bound_ma[[j]] + CALIBRATION_BOUND_TOL_MA
    }, logical(1))),
    config_path = spec$config_path_resolved[[1]],
    tree_path = spec$tree_path_resolved[[1]],
    stringsAsFactors = FALSE
  )
}

calibration_ages <- do.call(rbind, calibration_rows)
primary_age_lookup <- setNames(
  calibration_ages$node_age_ma[calibration_ages$variant_id == primary_id],
  calibration_ages$node_id[calibration_ages$variant_id == primary_id]
)
calibration_ages$primary_node_age_ma <- unname(primary_age_lookup[calibration_ages$node_id])
calibration_ages$signed_age_difference_from_primary_ma <-
  calibration_ages$node_age_ma - calibration_ages$primary_node_age_ma
calibration_ages$absolute_age_difference_from_primary_ma <-
  abs(calibration_ages$signed_age_difference_from_primary_ma)

internal_differences <- do.call(rbind, internal_difference_rows)
difference_summary <- do.call(rbind, summary_rows)

for (i in seq_len(nrow(difference_summary))) {
  add_info(
    difference_summary$variant_id[[i]],
    "max_abs_all_internal_node_age_difference_from_primary_ma",
    fmt_num(difference_summary$max_abs_all_internal_node_age_difference_from_primary_ma[[i]]),
    paste0(
      "Exact rooted-clade matching; matched=",
      difference_summary$matched_rooted_internal_node_count[[i]], "/",
      difference_summary$primary_internal_node_count[[i]], "."
    )
  )
  add_info(
    difference_summary$variant_id[[i]],
    "max_abs_calibration_node_age_difference_from_primary_ma",
    fmt_num(difference_summary$max_abs_calibration_node_age_difference_from_primary_ma[[i]]),
    "Maximum across the three exact calibration nodes."
  )
}

checks <- do.call(rbind, checks_list)
write.table(checks, output_paths[["checks"]], sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(calibration_ages, output_paths[["calibration_ages"]], sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(internal_differences, output_paths[["internal_differences"]], sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(difference_summary, output_paths[["difference_summary"]], sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

critical_failures <- checks$critical & checks$status == "FAIL"
report_lines <- c(
  "# Three-calibration final-variant QA",
  "",
  paste0("- Manifest: `", manifest_path, "`"),
  paste0("- Qualified 23-point CV decision: `", decision_path, "`"),
  paste0("- Primary variant: `", primary_id, "`"),
  paste0("- Variants evaluated: ", nrow(manifest)),
  paste0("- Sensitivity required: ", sensitivity_required_text),
  paste0(
    "- Required sensitivity smoothing set: ",
    if (length(decision_sensitivity_tokens)) {
      paste(decision_sensitivity_tokens, collapse = ",")
    } else {
      "none"
    }
  ),
  paste0("- Calibration-bound tolerance: ", fmt_num(CALIBRATION_BOUND_TOL_MA), " Ma"),
  paste0("- Critical failures: ", sum(critical_failures)),
  "",
  "All topologies are evaluated against the manifest-designated primary tree. Node ages are",
  "mean node-to-descendant-tip path lengths; descendant-tip minima, maxima, and ranges are",
  "retained to expose Newick serialization effects. Internal nodes are paired only by exact",
  "rooted descendant-tip sets.",
  "",
  "The three calibrated ages are constraint-conditioned and must not be interpreted as",
  "independent validation of their source intervals."
)
writeLines(report_lines, output_paths[["report"]], useBytes = TRUE)

if (any(critical_failures)) {
  failed <- paste(
    paste(checks$variant_id[critical_failures], checks$check_id[critical_failures], sep = ":"),
    collapse = ", "
  )
  stop(
    "Three-calibration final-variant QA hard gate failed: ", failed,
    ". Diagnostic outputs were retained at ", output_paths[["checks"]],
    call. = FALSE
  )
}

cat("Three-calibration final-variant QA PASS\n")
cat("Checks: ", output_paths[["checks"]], "\n", sep = "")
cat("Calibration ages: ", output_paths[["calibration_ages"]], "\n", sep = "")
cat("Internal-node differences: ", output_paths[["internal_differences"]], "\n", sep = "")
cat("Difference summary: ", output_paths[["difference_summary"]], "\n", sep = "")

