#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ape))

options(stringsAsFactors = FALSE, digits = 15)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop(
    paste(
      "Usage: 03_final_qa_leave_one_out.R <scenario_root>",
      "<C_root_only|D_papilionidae_only>"
    ),
    call. = FALSE
  )
}

scenario_root <- normalizePath(args[[1L]], mustWork = TRUE)
scenario <- args[[2L]]
allowed_scenarios <- c("C_root_only", "D_papilionidae_only")
if (!scenario %in% allowed_scenarios) {
  stop("Unknown scenario: ", scenario, call. = FALSE)
}
if (!identical(basename(scenario_root), scenario)) {
  stop(
    "Scenario argument/path mismatch: argument=", scenario,
    " path_basename=", basename(scenario_root), call. = FALSE
  )
}

config_dir <- file.path(scenario_root, "configs")
log_dir <- file.path(scenario_root, "logs")
output_dir <- file.path(scenario_root, "output")
qa_dir <- file.path(scenario_root, "qa")
if (!all(dir.exists(c(config_dir, log_dir, output_dir, qa_dir)))) {
  stop("Scenario is missing one or more required directories", call. = FALSE)
}

output_paths <- file.path(
  qa_dir,
  c(
    "node_ages.tsv",
    "final_tree_checks.tsv",
    "qa_report.md",
    "final_variant_node_ages.tsv",
    "final_variant_checks.tsv",
    "provenance_sha256.tsv"
  )
)
existing_outputs <- output_paths[file.exists(output_paths)]
if (length(existing_outputs)) {
  stop(
    "Refusing to overwrite existing output(s): ",
    paste(existing_outputs, collapse = ", "), call. = FALSE
  )
}

INPUT_TREE <- paste0(
  "/home/data/t200301/xjt/Hesperiinae_review/",
  "treePL_R2_3_T25_deep_only_20260904_000707/",
  "input/verified_495_tip_input.treefile"
)
INPUT_TREE_SHA256 <- paste0(
  "4d35d553ed5bbf08022daef2361ffa0c",
  "342f8694a53aaf08d5350692094439e5"
)
TREEPL_BINARY <- paste0(
  "/home/data/t200301/xjt/Hesperiinae_review/",
  "treePL_R2_3_T25_deep_only_20260904_000707/",
  "software/treePL_randomcv_fix_src/treePL"
)
TREEPL_BINARY_SHA256 <- paste0(
  "5d6bfac8d18d835f24ebb79f0373b070a",
  "948b58a179f437f0aacb8fed74ccc5d"
)
EXPECTED_NUMSITES <- "182682"
EXPECTED_FINAL_SEED <- "2026090599"
ULTRAMETRIC_TOL_MA <- 1e-4
CALIBRATION_TOL_MA <- 1e-4

calibrations <- list(
  C_root_only = list(
    node_id = "CROWN_PAPILIONOIDEA",
    lines = c(
      paste(
        "mrca = CROWN_PAPILIONOIDEA",
        "Parnassius_apollo_ncbi Aeromachus_catocyanea_mydata"
      ),
      "min = CROWN_PAPILIONOIDEA 91.5046",
      "max = CROWN_PAPILIONOIDEA 100.8925"
    )
  ),
  D_papilionidae_only = list(
    node_id = "CROWN_PAPILIONIDAE_NONBARONIA",
    lines = c(
      paste(
        "mrca = CROWN_PAPILIONIDAE_NONBARONIA",
        "Parnassius_apollo_ncbi Graphium_cloanthus_mydata"
      ),
      "min = CROWN_PAPILIONIDAE_NONBARONIA 44.1968",
      "max = CROWN_PAPILIONIDAE_NONBARONIA 52.9473"
    )
  )
)
scenario_calibration <- calibrations[[scenario]]

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
  stringsAsFactors = FALSE
)
node_specs$calibrated_in_scenario <-
  node_specs$node_id == scenario_calibration$node_id

canonical <- function(path, must_work = FALSE) {
  normalizePath(path, mustWork = must_work)
}

sha256_file <- function(path) {
  if (!file.exists(path) || dir.exists(path)) {
    stop("Missing file for SHA-256: ", path, call. = FALSE)
  }
  executable <- Sys.which("sha256sum")
  if (!nzchar(executable)) stop("sha256sum executable not found", call. = FALSE)
  result <- system2(executable, c("--", path), stdout = TRUE, stderr = TRUE)
  status <- attr(result, "status")
  if ((!is.null(status) && status != 0L) || length(result) != 1L) {
    stop("sha256sum failed for ", path, call. = FALSE)
  }
  token <- strsplit(result[[1L]], "[[:space:]]+", perl = TRUE)[[1L]][[1L]]
  if (!grepl("^[0-9a-f]{64}$", token)) {
    stop("Malformed SHA-256 result for ", path, call. = FALSE)
  }
  token
}

read_active_config <- function(path) {
  if (!file.exists(path) || file.info(path)$size <= 0) {
    stop("Missing/empty config: ", path, call. = FALSE)
  }
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

one_config_value <- function(parsed, key, path) {
  values <- parsed$value[parsed$key == key]
  if (length(values) != 1L) {
    stop("Expected exactly one ", key, " in ", path, call. = FALSE)
  }
  values[[1L]]
}

read_metadata <- function(path) {
  if (!file.exists(path) || file.info(path)$size <= 0) {
    stop("Missing/empty metadata: ", path, call. = FALSE)
  }
  tab <- read.delim(
    path, header = FALSE, sep = "\t", quote = "", comment.char = "",
    col.names = c("key", "value"), stringsAsFactors = FALSE,
    colClasses = "character"
  )
  one <- function(key) {
    value <- tab$value[tab$key == key]
    if (length(value) != 1L) {
      stop("Metadata key missing or non-unique: ", key, " in ", path)
    }
    value[[1L]]
  }
  list(table = tab, one = one)
}

parse_bool <- function(x) {
  if (length(x) != 1L || is.na(x)) return(NA)
  x <- toupper(trimws(as.character(x)))
  if (x %in% c("TRUE", "T", "1")) return(TRUE)
  if (x %in% c("FALSE", "F", "0")) return(FALSE)
  NA
}

parse_numeric_csv <- function(x) {
  if (length(x) != 1L || is.na(x) || !nzchar(trimws(x)) ||
      identical(tolower(trimws(x)), "none")) {
    return(numeric())
  }
  values <- suppressWarnings(as.numeric(trimws(strsplit(x, ",", fixed = TRUE)[[1L]])))
  if (any(!is.finite(values)) || any(values <= 0)) {
    stop("Malformed positive numeric CSV value: ", x, call. = FALSE)
  }
  values
}

numeric_set_equal <- function(a, b, tolerance = 1e-12) {
  if (length(a) != length(b)) return(FALSE)
  if (!length(a)) return(TRUE)
  a <- sort(a)
  b <- sort(b)
  all(abs(a - b) <= tolerance * pmax(abs(a), abs(b), .Machine$double.xmin))
}

read_nonempty_lines <- function(path, label) {
  if (!file.exists(path) || file.info(path)$size <= 0) {
    stop("Missing/empty ", label, ": ", path, call. = FALSE)
  }
  readLines(path, warn = FALSE)
}

validate_resource_log <- function(path) {
  lines <- read_nonempty_lines(path, "resource log")
  hits <- grep("^[[:space:]]*Exit status:[[:space:]]*", lines, value = TRUE)
  length(hits) == 1L &&
    grepl("^[[:space:]]*Exit status:[[:space:]]*0[[:space:]]*$", hits)
}

descendant_context <- function(tree) {
  n_tip <- Ntip(tree)
  roots <- setdiff(unique(tree$edge[, 1L]), unique(tree$edge[, 2L]))
  if (length(roots) != 1L) stop("Tree does not have exactly one root", call. = FALSE)
  root <- as.integer(roots[[1L]])
  children <- split(tree$edge[, 2L], tree$edge[, 1L])
  memo <- new.env(parent = emptyenv())
  descendant_ids <- function(node) {
    key <- as.character(node)
    if (exists(key, envir = memo, inherits = FALSE)) {
      return(get(key, envir = memo, inherits = FALSE))
    }
    result <- if (node <= n_tip) {
      as.integer(node)
    } else {
      kids <- children[[key]]
      if (is.null(kids) || length(kids) != 2L) {
        stop("Non-binary or malformed internal node: ", node, call. = FALSE)
      }
      sort(unique(unlist(lapply(kids, descendant_ids), use.names = FALSE)))
    }
    assign(key, result, envir = memo)
    result
  }
  descendant_labels <- function(node) {
    sort(tree$tip.label[descendant_ids(node)])
  }
  depth <- node.depth.edgelength(tree)
  internal_nodes <- seq.int(n_tip + 1L, n_tip + tree$Nnode)
  signatures <- sort(vapply(
    internal_nodes,
    function(node) paste(descendant_labels(node), collapse = "\034"),
    character(1)
  ))
  list(
    root = root,
    depth = depth,
    descendant_ids = descendant_ids,
    descendant_labels = descendant_labels,
    signatures = signatures
  )
}

read_tree_strict <- function(path, label) {
  if (!file.exists(path) || file.info(path)$size <= 0) {
    stop("Missing/empty ", label, ": ", path, call. = FALSE)
  }
  tree <- tryCatch(
    read.tree(path),
    error = function(e) stop("Cannot parse ", label, ": ", conditionMessage(e))
  )
  if (!inherits(tree, "phylo") || length(tree$tip.label) != 495L ||
      anyNA(tree$tip.label) || any(!nzchar(tree$tip.label)) ||
      anyDuplicated(tree$tip.label)) {
    stop(label, " must contain exactly 495 unique nonempty tips", call. = FALSE)
  }
  if (!is.rooted(tree) || !is.binary.tree(tree) || tree$Nnode != 494L) {
    stop(label, " must be a rooted fully binary 495-tip tree", call. = FALSE)
  }
  if (is.null(tree$edge.length)) {
    stop(label, " has no branch lengths", call. = FALSE)
  }
  tree
}

safe_mrca <- function(tree, tip_1, tip_2) {
  if (sum(tree$tip.label == tip_1) != 1L || sum(tree$tip.label == tip_2) != 1L) {
    return(NA_integer_)
  }
  as.integer(getMRCA(tree, c(tip_1, tip_2)))
}

node_age_record <- function(tree, context, spec, role, stage, smoothing) {
  node <- safe_mrca(tree, spec$anchor_tip_1, spec$anchor_tip_2)
  if (length(node) != 1L || is.na(node)) {
    stop("Could not map focal MRCA: ", spec$node_id, call. = FALSE)
  }
  ids <- context$descendant_ids(node)
  path_ages <- context$depth[ids] - context$depth[node]
  age <- mean(path_ages)
  distance_lower <- age - spec$lower_bound_ma
  distance_upper <- spec$upper_bound_ma - age
  calibrated <- isTRUE(spec$calibrated_in_scenario)
  boundary_status <- if (!calibrated) {
    "not_calibrated"
  } else if (abs(distance_lower) <= CALIBRATION_TOL_MA) {
    "lower_boundary_hit"
  } else if (abs(distance_upper) <= CALIBRATION_TOL_MA) {
    "upper_boundary_hit"
  } else if (distance_lower >= -CALIBRATION_TOL_MA &&
             distance_upper >= -CALIBRATION_TOL_MA) {
    "within_bounds_not_boundary"
  } else {
    "outside_calibration_bounds"
  }
  data.frame(
    scenario = scenario,
    variant_role = role,
    stage = stage,
    smoothing = smoothing,
    node_id = spec$node_id,
    biological_definition = spec$biological_definition,
    anchor_tip_1 = spec$anchor_tip_1,
    anchor_tip_2 = spec$anchor_tip_2,
    mrca_node = node,
    is_root = node == context$root,
    descendant_tip_count = length(ids),
    expected_descendant_tip_count = spec$expected_descendant_tip_count,
    descendant_count_exact = length(ids) == spec$expected_descendant_tip_count,
    calibrated_in_scenario = calibrated,
    node_age_ma = age,
    descendant_tip_age_min_ma = min(path_ages),
    descendant_tip_age_max_ma = max(path_ages),
    descendant_tip_age_range_ma = diff(range(path_ages)),
    reference_lower_bound_ma = spec$lower_bound_ma,
    reference_upper_bound_ma = spec$upper_bound_ma,
    signed_distance_above_lower_ma = distance_lower,
    signed_distance_below_upper_ma = distance_upper,
    calibrated_node_within_bounds_1e_4 =
      !calibrated || (distance_lower >= -CALIBRATION_TOL_MA &&
                        distance_upper >= -CALIBRATION_TOL_MA),
    boundary_status_1e_4 = boundary_status,
    age_definition = paste0(
      "mean node-to-descendant-tip path length; min/max/range retained ",
      "as serialization and ultrametric diagnostics"
    ),
    stringsAsFactors = FALSE
  )
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

input_tree_path <- canonical(INPUT_TREE, must_work = TRUE)
binary_path <- canonical(TREEPL_BINARY, must_work = TRUE)
input_hash <- sha256_file(input_tree_path)
binary_hash <- sha256_file(binary_path)
add_check(
  "fixed_input_tree_sha256", "provenance",
  identical(input_hash, INPUT_TREE_SHA256), input_hash, INPUT_TREE_SHA256
)
add_check(
  "fixed_treepl_binary_sha256", "provenance",
  identical(binary_hash, TREEPL_BINARY_SHA256), binary_hash, TREEPL_BINARY_SHA256
)

decision_path <- file.path(qa_dir, "cv_decision.tsv")
decision <- read.delim(
  decision_path, sep = "\t", quote = "", comment.char = "",
  check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character",
  na.strings = c("NA", "")
)
required_decision_columns <- c(
  "scenario", "status", "selected_smoothing", "grid_extension_direction",
  "aggregate_winner_at_boundary", "sensitivity_required",
  "sensitivity_smoothing_values", "config_validation_passed",
  "prime_six_line_validation_passed", "metadata_validation_passed",
  "stderr_validation_passed", "stdout_and_score_validation_passed",
  "resource_exit_zero_validation_passed",
  "dated_tree_495_tip_rooted_topology_validation_passed",
  "input_tree_sha256", "treepl_binary_sha256"
)
if (nrow(decision) != 1L ||
    !all(required_decision_columns %in% names(decision))) {
  stop(
    "cv_decision.tsv must be one row and contain all required columns",
    call. = FALSE
  )
}
dv <- function(key) decision[[key]][[1L]]
selected_smoothing <- suppressWarnings(as.numeric(dv("selected_smoothing")))
if (length(selected_smoothing) != 1L || !is.finite(selected_smoothing) ||
    selected_smoothing <= 0) {
  stop("CV decision has no valid positive selected_smoothing", call. = FALSE)
}
sensitivity_required <- parse_bool(dv("sensitivity_required"))
sensitivity_values <- parse_numeric_csv(dv("sensitivity_smoothing_values"))
decision_validation_flags <- c(
  "config_validation_passed", "prime_six_line_validation_passed",
  "metadata_validation_passed", "stderr_validation_passed",
  "stdout_and_score_validation_passed",
  "resource_exit_zero_validation_passed",
  "dated_tree_495_tip_rooted_topology_validation_passed"
)
add_check("decision_scenario_exact", "cv_decision", identical(dv("scenario"), scenario), dv("scenario"), scenario)
add_check(
  "decision_status_accepted", "cv_decision",
  dv("status") %in% c("ACCEPTED", "ACCEPTED_WITH_SENSITIVITIES"),
  dv("status"), "ACCEPTED or ACCEPTED_WITH_SENSITIVITIES"
)
add_check("decision_grid_extension_none", "cv_decision", identical(dv("grid_extension_direction"), "none"), dv("grid_extension_direction"), "none")
add_check("decision_winner_not_at_boundary", "cv_decision", identical(parse_bool(dv("aggregate_winner_at_boundary")), FALSE), dv("aggregate_winner_at_boundary"), FALSE)
add_check("decision_input_hash_exact", "cv_decision", identical(dv("input_tree_sha256"), INPUT_TREE_SHA256), dv("input_tree_sha256"), INPUT_TREE_SHA256)
add_check("decision_binary_hash_exact", "cv_decision", identical(dv("treepl_binary_sha256"), TREEPL_BINARY_SHA256), dv("treepl_binary_sha256"), TREEPL_BINARY_SHA256)
add_check(
  "decision_prior_validation_flags_true", "cv_decision",
  all(vapply(decision_validation_flags, function(x) isTRUE(parse_bool(dv(x))), logical(1))),
  paste(paste(decision_validation_flags, vapply(decision_validation_flags, dv, character(1)), sep = "="), collapse = ";"),
  "all TRUE"
)
add_check(
  "decision_sensitivity_status_consistent", "cv_decision",
  (identical(dv("status"), "ACCEPTED") && identical(sensitivity_required, FALSE)) ||
    (identical(dv("status"), "ACCEPTED_WITH_SENSITIVITIES") &&
       identical(sensitivity_required, TRUE)),
  paste0("status=", dv("status"), ";required=", dv("sensitivity_required")),
  "ACCEPTED iff FALSE; ACCEPTED_WITH_SENSITIVITIES iff TRUE"
)

prime_path <- file.path(qa_dir, "prime_recommendations.cfg")
prime <- read_active_config(prime_path)
expected_prime_keys <- c(
  "opt", "moredetail", "optad", "moredetailad", "optcvad", "moredetailcvad"
)
add_check(
  "scenario_prime_has_exact_six_directives", "prime",
  length(prime$active) == 6L && identical(prime$key, expected_prime_keys),
  paste(prime$active, collapse = " | "),
  paste(expected_prime_keys, collapse = " | ")
)
prime_integer_text <- prime$value[prime$key %in% c("opt", "optad", "optcvad")]
prime_integer_valid <- length(prime_integer_text) == 3L &&
  all(grepl("^[1-9][0-9]*$", prime_integer_text))
add_check(
  "scenario_prime_integer_values_positive", "prime", prime_integer_valid,
  paste(prime_integer_text, collapse = ","), "three positive integers"
)

primary_config <- file.path(config_dir, "final.cfg")
repeat_config <- file.path(config_dir, "final_repeat_same_seed.cfg")
sensitivity_configs <- sort(list.files(
  config_dir,
  pattern = "^final_sensitivity_.+\\.cfg$",
  full.names = TRUE
))
if (!file.exists(primary_config) || !file.exists(repeat_config)) {
  stop("Missing final.cfg or final_repeat_same_seed.cfg", call. = FALSE)
}

stage_manifest <- data.frame(
  variant_role = c("primary", "same_seed_repeat"),
  stage = c("final_primary_thorough", "final_repeat_same_seed"),
  config = c(primary_config, repeat_config),
  stringsAsFactors = FALSE
)
if (length(sensitivity_configs)) {
  sensitivity_stages <- sub("\\.cfg$", "", basename(sensitivity_configs))
  stage_manifest <- rbind(
    stage_manifest,
    data.frame(
      variant_role = rep("sensitivity", length(sensitivity_configs)),
      stage = sensitivity_stages,
      config = sensitivity_configs,
      stringsAsFactors = FALSE
    )
  )
}

expected_primary_output <- file.path(output_dir, paste0(scenario, "_dated.tre"))
expected_repeat_output <- file.path(
  output_dir, paste0(scenario, "_dated_repeat_same_seed.tre")
)
stage_manifest$output <- c(
  expected_primary_output,
  expected_repeat_output,
  if (nrow(stage_manifest) > 2L) {
    file.path(
      output_dir,
      paste0(
        scenario, "_dated_sensitivity_",
        sub("^final_sensitivity_", "", stage_manifest$stage[-c(1L, 2L)]),
        ".tre"
      )
    )
  } else character()
)

config_info <- lapply(seq_len(nrow(stage_manifest)), function(i) {
  config_path <- canonical(stage_manifest$config[[i]], must_work = TRUE)
  parsed <- read_active_config(config_path)
  calibration_lines <- parsed$active[parsed$key %in% c("mrca", "min", "max")]
  optimization_lines <- parsed$active[parsed$key %in% expected_prime_keys]
  treefile <- canonical(one_config_value(parsed, "treefile", config_path), TRUE)
  outfile <- canonical(one_config_value(parsed, "outfile", config_path), FALSE)
  smoothing <- suppressWarnings(as.numeric(one_config_value(parsed, "smooth", config_path)))
  seed <- one_config_value(parsed, "seed", config_path)
  prohibited <- c(
    "prime", "cv", "randomcv", "cvstart", "cvstop", "cvmultstep",
    "cviter", "cvoutfile", "log_pen"
  )
  checks <- c(
    treefile_exact = identical(treefile, input_tree_path),
    outfile_exact = identical(outfile, canonical(stage_manifest$output[[i]], FALSE)),
    numsites_exact = identical(one_config_value(parsed, "numsites", config_path), EXPECTED_NUMSITES),
    single_scenario_calibration_exact = identical(calibration_lines, scenario_calibration$lines),
    exactly_one_mrca_min_max = sum(parsed$key == "mrca") == 1L &&
      sum(parsed$key == "min") == 1L && sum(parsed$key == "max") == 1L,
    scenario_prime_exact = identical(optimization_lines, prime$active),
    positive_smoothing = length(smoothing) == 1L && is.finite(smoothing) && smoothing > 0,
    seed_exact = identical(seed, EXPECTED_FINAL_SEED),
    thorough_exactly_once = sum(parsed$active == "thorough") == 1L,
    prohibited_modes_absent = !any(parsed$key %in% prohibited)
  )
  for (name in names(checks)) {
    add_check(
      paste0("config_", name), stage_manifest$stage[[i]], checks[[name]],
      checks[[name]], TRUE
    )
  }
  list(parsed = parsed, smoothing = smoothing, config_path = config_path)
})

manifest_smoothing <- vapply(config_info, function(x) x$smoothing, numeric(1))
stage_manifest$smoothing <- manifest_smoothing
observed_sensitivity_values <- manifest_smoothing[stage_manifest$variant_role == "sensitivity"]
add_check(
  "primary_smoothing_equals_cv_selection", "config_vs_cv",
  numeric_set_equal(manifest_smoothing[[1L]], selected_smoothing),
  manifest_smoothing[[1L]], selected_smoothing
)
add_check(
  "repeat_smoothing_equals_cv_selection", "config_vs_cv",
  numeric_set_equal(manifest_smoothing[[2L]], selected_smoothing),
  manifest_smoothing[[2L]], selected_smoothing
)
add_check(
  "sensitivity_config_set_exactly_matches_decision", "config_vs_cv",
  numeric_set_equal(observed_sensitivity_values, sensitivity_values) &&
    anyDuplicated(observed_sensitivity_values) == 0L,
  paste(sort(observed_sensitivity_values), collapse = ","),
  paste(sort(sensitivity_values), collapse = ",")
)
add_check(
  "sensitivity_presence_matches_required_flag", "config_vs_cv",
  (identical(sensitivity_required, TRUE) && length(sensitivity_configs) > 0L) ||
    (identical(sensitivity_required, FALSE) && length(sensitivity_configs) == 0L),
  paste0("required=", sensitivity_required, ";configs=", length(sensitivity_configs)),
  "TRUE with >=1 config, or FALSE with 0 configs"
)
add_check(
  "primary_excluded_from_sensitivity_set", "config_vs_cv",
  !any(abs(observed_sensitivity_values - selected_smoothing) <=
         1e-12 * pmax(abs(observed_sensitivity_values), abs(selected_smoothing))),
  paste(sort(observed_sensitivity_values), collapse = ","),
  paste0("must exclude ", selected_smoothing)
)

input_tree <- read_tree_strict(input_tree_path, "fixed input tree")
input_context <- descendant_context(input_tree)
input_nodes <- vapply(seq_len(nrow(node_specs)), function(i) {
  safe_mrca(input_tree, node_specs$anchor_tip_1[[i]], node_specs$anchor_tip_2[[i]])
}, integer(1))
input_desc_counts <- vapply(
  input_nodes,
  function(node) length(input_context$descendant_ids(node)),
  integer(1)
)
add_check(
  "input_focal_anchor_descendant_counts_exact", "input_tree",
  identical(input_desc_counts, node_specs$expected_descendant_tip_count),
  paste(input_desc_counts, collapse = ","),
  paste(node_specs$expected_descendant_tip_count, collapse = ",")
)
add_check(
  "input_root_anchor_maps_to_root", "input_tree",
  identical(input_nodes[[1L]], input_context$root), input_nodes[[1L]],
  input_context$root
)

variant_checks <- list()
variant_node_ages <- list()
provenance_rows <- list()
add_provenance <- function(role, stage, path) {
  provenance_rows[[length(provenance_rows) + 1L]] <<- data.frame(
    role = role,
    stage = stage,
    path = canonical(path, must_work = TRUE),
    size_bytes = as.numeric(file.info(path)$size),
    sha256 = sha256_file(path),
    hash_timing = paste0(
      "computed by 03_final_qa_leave_one_out.R after stage completion; ",
      "not an execution-time attestation"
    ),
    stringsAsFactors = FALSE
  )
}

add_provenance("fixed_input_tree", "global", input_tree_path)
add_provenance("patched_treepl_binary", "global", binary_path)
add_provenance("cv_decision", "global", decision_path)
add_provenance("scenario_prime_recommendations", "global", prime_path)
command_args_all <- commandArgs(trailingOnly = FALSE)
script_args <- grep("^--file=", command_args_all, value = TRUE)
if (length(script_args) != 1L) {
  stop("Could not identify exactly one executing QA script path", call. = FALSE)
}
script_path <- canonical(sub("^--file=", "", script_args[[1L]]), must_work = TRUE)
add_provenance("qa_script", "global", script_path)

for (i in seq_len(nrow(stage_manifest))) {
  role <- stage_manifest$variant_role[[i]]
  stage <- stage_manifest$stage[[i]]
  config_path <- config_info[[i]]$config_path
  output_path <- canonical(stage_manifest$output[[i]], must_work = TRUE)
  metadata_path <- file.path(log_dir, paste0(stage, ".metadata.tsv"))
  stdout_path <- file.path(log_dir, paste0(stage, ".stdout"))
  stderr_path <- file.path(log_dir, paste0(stage, ".stderr"))
  resources_path <- file.path(log_dir, paste0(stage, ".resources.txt"))
  metadata <- read_metadata(metadata_path)
  stdout_lines <- read_nonempty_lines(stdout_path, "stdout")
  stderr_exists_empty <- file.exists(stderr_path) && file.info(stderr_path)$size == 0
  resources_exit_zero <- validate_resource_log(resources_path)
  metadata_required <- c(
    "stage", "pid", "start_time", "binary", "config", "end_time",
    "exit_code", "elapsed_s"
  )
  metadata_counts <- vapply(
    metadata_required, function(key) sum(metadata$table$key == key), integer(1)
  )
  metadata_stage_exact <- identical(metadata$one("stage"), stage)
  metadata_config_exact <- identical(
    canonical(metadata$one("config"), TRUE), config_path
  )
  metadata_binary_exact <- identical(
    canonical(metadata$one("binary"), TRUE), binary_path
  )
  metadata_exit_zero <- identical(metadata$one("exit_code"), "0")
  failure_re <- paste(
    c(
      "complete failure", "segmentation",
      "(^|[^[:alpha:]])nan([^[:alpha:]]|$)",
      "(^|[^[:alpha:]])inf([^[:alpha:]]|$)",
      "(^|[^[:alpha:]])error([^[:alpha:]]|$)"
    ),
    collapse = "|"
  )
  stdout_no_failure_token <- !any(grepl(
    failure_re, stdout_lines, ignore.case = TRUE, perl = TRUE
  ))

  tree <- read_tree_strict(output_path, paste0(stage, " output tree"))
  context <- descendant_context(tree)
  same_tip_set <- setequal(tree$tip.label, input_tree$tip.label)
  rf <- if (same_tip_set) {
    suppressWarnings(as.numeric(dist.topo(input_tree, tree, method = "PH85")))
  } else NA_real_
  rooted_clades_exact <- same_tip_set &&
    identical(context$signatures, input_context$signatures)
  finite_edges <- all(is.finite(tree$edge.length))
  negative_count <- sum(tree$edge.length < 0, na.rm = TRUE)
  zero_count <- sum(tree$edge.length == 0, na.rm = TRUE)
  tip_depth <- context$depth[seq_len(Ntip(tree))]
  root_to_tip_range <- diff(range(tip_depth))

  node_rows <- do.call(rbind, lapply(seq_len(nrow(node_specs)), function(j) {
    node_age_record(
      tree, context, node_specs[j, , drop = FALSE], role, stage,
      manifest_smoothing[[i]]
    )
  }))
  rownames(node_rows) <- NULL
  focal_counts_exact <- all(node_rows$descendant_count_exact)
  root_anchor_exact <- node_rows$mrca_node[node_rows$node_id == "CROWN_PAPILIONOIDEA"] ==
    context$root
  calibrated_rows <- node_rows[node_rows$calibrated_in_scenario, , drop = FALSE]
  one_calibrated_node <- nrow(calibrated_rows) == 1L &&
    identical(calibrated_rows$node_id, scenario_calibration$node_id)
  calibrated_within_bounds <- one_calibrated_node &&
    all(calibrated_rows$calibrated_node_within_bounds_1e_4)

  stage_gate_values <- c(
    config_exists_nonempty = file.exists(config_path) && file.info(config_path)$size > 0,
    metadata_required_keys_once = all(metadata_counts == 1L),
    metadata_stage_exact = metadata_stage_exact,
    metadata_config_exact = metadata_config_exact,
    metadata_binary_exact = metadata_binary_exact,
    metadata_exit_zero = metadata_exit_zero,
    stderr_exists_empty = stderr_exists_empty,
    stdout_exists_nonempty = file.exists(stdout_path) && file.info(stdout_path)$size > 0,
    stdout_no_failure_token = stdout_no_failure_token,
    resources_exists_nonempty = file.exists(resources_path) && file.info(resources_path)$size > 0,
    resources_exit_zero = resources_exit_zero,
    output_exists_nonempty = file.exists(output_path) && file.info(output_path)$size > 0,
    tip_count_495 = Ntip(tree) == 495L,
    unique_tip_labels = !anyDuplicated(tree$tip.label),
    rooted = is.rooted(tree),
    fully_binary = is.binary.tree(tree) && tree$Nnode == 494L,
    tip_set_exact = same_tip_set,
    rooted_rf_zero = length(rf) == 1L && is.finite(rf) && rf == 0,
    rooted_clades_exact = rooted_clades_exact,
    finite_branch_lengths = finite_edges,
    negative_branch_count_zero = negative_count == 0L,
    zero_branch_count_zero = zero_count == 0L,
    root_to_tip_range_below_1e_4 = is.finite(root_to_tip_range) &&
      root_to_tip_range < ULTRAMETRIC_TOL_MA,
    focal_anchor_descendant_counts_exact = focal_counts_exact,
    root_anchor_maps_to_root = isTRUE(root_anchor_exact),
    exactly_scenario_calibration_mapped = one_calibrated_node,
    scenario_calibrated_node_within_bounds_1e_4 = calibrated_within_bounds
  )
  for (name in names(stage_gate_values)) {
    add_check(
      paste0("variant_", name), stage, stage_gate_values[[name]],
      stage_gate_values[[name]], TRUE
    )
  }

  variant_checks[[length(variant_checks) + 1L]] <- data.frame(
    scenario = scenario,
    variant_role = role,
    stage = stage,
    smoothing = manifest_smoothing[[i]],
    config_path = config_path,
    output_path = output_path,
    output_sha256 = sha256_file(output_path),
    tip_count = Ntip(tree),
    unique_tip_count = length(unique(tree$tip.label)),
    rooted = is.rooted(tree),
    fully_binary = is.binary.tree(tree) && tree$Nnode == 494L,
    tip_set_exact = same_tip_set,
    rooted_rf = rf,
    rooted_clades_exact = rooted_clades_exact,
    finite_branch_lengths = finite_edges,
    negative_branch_count = negative_count,
    zero_branch_count = zero_count,
    root_to_tip_min_ma = min(tip_depth),
    root_to_tip_max_ma = max(tip_depth),
    root_to_tip_range_ma = root_to_tip_range,
    root_to_tip_range_below_1e_4 = root_to_tip_range < ULTRAMETRIC_TOL_MA,
    focal_anchor_descendant_counts_exact = focal_counts_exact,
    exactly_one_scenario_calibration = one_calibrated_node,
    calibrated_node_within_bounds_1e_4 = calibrated_within_bounds,
    metadata_stage_config_binary_exit0 = metadata_stage_exact &&
      metadata_config_exact && metadata_binary_exact && metadata_exit_zero,
    stderr_empty = stderr_exists_empty,
    stdout_nonempty_no_failure_token =
      file.info(stdout_path)$size > 0 && stdout_no_failure_token,
    resources_exit_zero = resources_exit_zero,
    all_hard_gates_pass = all(stage_gate_values),
    stringsAsFactors = FALSE
  )
  variant_node_ages[[length(variant_node_ages) + 1L]] <- node_rows

  add_provenance("config", stage, config_path)
  add_provenance("metadata", stage, metadata_path)
  add_provenance("stdout", stage, stdout_path)
  add_provenance("stderr", stage, stderr_path)
  add_provenance("resources", stage, resources_path)
  add_provenance("dated_tree", stage, output_path)
}

variant_checks <- do.call(rbind, variant_checks)
variant_node_ages <- do.call(rbind, variant_node_ages)
provenance <- do.call(rbind, provenance_rows)
rownames(variant_checks) <- NULL
rownames(variant_node_ages) <- NULL
rownames(provenance) <- NULL

primary_hash <- variant_checks$output_sha256[variant_checks$variant_role == "primary"]
repeat_hash <- variant_checks$output_sha256[
  variant_checks$variant_role == "same_seed_repeat"
]
same_seed_byte_identical <- length(primary_hash) == 1L &&
  length(repeat_hash) == 1L && identical(primary_hash, repeat_hash)
add_check(
  "same_seed_repeat_raw_byte_identical", "repeatability",
  same_seed_byte_identical, paste(primary_hash, repeat_hash, sep = " vs "),
  "identical SHA-256"
)

checks <- do.call(rbind, checks_list)
rownames(checks) <- NULL
critical_failures <- checks$status == "FAIL"
overall_status <- if (any(critical_failures)) "FAIL" else "PASS"

primary_node_ages <- variant_node_ages[
  variant_node_ages$variant_role == "primary", , drop = FALSE
]

write.table(
  primary_node_ages, output_paths[[1L]], sep = "\t", quote = FALSE,
  row.names = FALSE, na = "NA"
)
write.table(
  checks, output_paths[[2L]], sep = "\t", quote = FALSE,
  row.names = FALSE, na = "NA"
)
write.table(
  variant_node_ages, output_paths[[4L]], sep = "\t", quote = FALSE,
  row.names = FALSE, na = "NA"
)
write.table(
  variant_checks, output_paths[[5L]], sep = "\t", quote = FALSE,
  row.names = FALSE, na = "NA"
)
write.table(
  provenance, output_paths[[6L]], sep = "\t", quote = FALSE,
  row.names = FALSE, na = "NA"
)

fmt <- function(x) {
  if (length(x) != 1L || is.na(x) || !is.finite(x)) return("NA")
  sprintf("%.12g", x)
}
md_escape <- function(x) gsub("\\|", "\\\\|", as.character(x))
node_table <- c(
  "| Node | Calibrated | Age (Ma) | Reference interval (Ma) | Boundary status | Descendants |",
  "|---|---:|---:|---:|---|---:|",
  vapply(seq_len(nrow(primary_node_ages)), function(i) {
    paste0(
      "| ", primary_node_ages$node_id[[i]],
      " | ", primary_node_ages$calibrated_in_scenario[[i]],
      " | ", fmt(primary_node_ages$node_age_ma[[i]]),
      " | ", fmt(primary_node_ages$reference_lower_bound_ma[[i]]), "-",
      fmt(primary_node_ages$reference_upper_bound_ma[[i]]),
      " | ", primary_node_ages$boundary_status_1e_4[[i]],
      " | ", primary_node_ages$descendant_tip_count[[i]], " |"
    )
  }, character(1))
)
variant_table <- c(
  "| Variant | Smooth | RF | Root-to-tip range (Ma) | Negative | Zero | Calibration pass | All hard gates |",
  "|---|---:|---:|---:|---:|---:|---:|---:|",
  vapply(seq_len(nrow(variant_checks)), function(i) {
    paste0(
      "| ", variant_checks$stage[[i]],
      " | ", fmt(variant_checks$smoothing[[i]]),
      " | ", fmt(variant_checks$rooted_rf[[i]]),
      " | ", fmt(variant_checks$root_to_tip_range_ma[[i]]),
      " | ", variant_checks$negative_branch_count[[i]],
      " | ", variant_checks$zero_branch_count[[i]],
      " | ", variant_checks$calibrated_node_within_bounds_1e_4[[i]],
      " | ", variant_checks$all_hard_gates_pass[[i]], " |"
    )
  }, character(1))
)
failed_ids <- checks$check_id[critical_failures]
report <- c(
  paste0("# Leave-one-out final QA: ", scenario),
  "",
  paste0("**Overall status:** ", overall_status),
  "",
  paste0("- Fixed input tree: `", input_tree_path, "`"),
  paste0("- Fixed input SHA-256: `", input_hash, "`"),
  paste0("- Patched treePL binary: `", binary_path, "`"),
  paste0("- Patched binary SHA-256: `", binary_hash, "`"),
  paste0("- CV decision: `", decision_path, "`"),
  paste0("- Selected smoothing: `", fmt(selected_smoothing), "`"),
  paste0("- Final seed: `", EXPECTED_FINAL_SEED, "`"),
  paste0("- Scenario calibration: `", scenario_calibration$node_id, "` only"),
  paste0("- Final variants validated: `", nrow(variant_checks), "`"),
  paste0("- Same-seed repeat raw-byte identical: `", same_seed_byte_identical, "`"),
  paste0("- Ultrametric and calibration tolerance: `", ULTRAMETRIC_TOL_MA, " Ma`"),
  "",
  "The three focal node ages below were independently extracted from each dated tree using fixed anchor pairs. Only the scenario's single calibration is treated as constrained; the other two intervals are reference values and are marked `not_calibrated`.",
  "",
  "## Primary focal node ages",
  "",
  node_table,
  "",
  "## Final-variant hard gates",
  "",
  variant_table,
  "",
  "Every primary, repeat, and sensitivity tree is required to have 495 unique tips, be rooted and fully binary, retain RF=0 and the exact rooted clade set, have finite strictly positive branch lengths, have root-to-tip range <1e-4 Ma, and satisfy the scenario calibration within 1e-4 Ma.",
  "",
  "## Critical failures",
  "",
  if (length(failed_ids)) {
    paste0("- `", md_escape(failed_ids), "`")
  } else {
    "None."
  },
  "",
  "## Provenance note",
  "",
  "`provenance_sha256.tsv` contains hashes computed after run completion. The stage metadata must independently identify the exact config and patched binary and record exit code 0; post-run hashes are not represented as launch-time attestations."
)
writeLines(report, output_paths[[3L]], useBytes = TRUE)

if (any(critical_failures)) {
  stop(
    "Leave-one-out final QA hard gate failed: ",
    paste(failed_ids, collapse = ", "),
    ". See ", output_paths[[2L]], call. = FALSE
  )
}

cat(
  "Leave-one-out final QA PASS: ", nrow(checks), " hard checks; ",
  nrow(variant_checks), " final variants\n", sep = ""
)
cat("Primary node ages: ", output_paths[[1L]], "\n", sep = "")
cat("QA report: ", output_paths[[3L]], "\n", sep = "")

