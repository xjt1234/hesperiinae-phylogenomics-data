#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop(
    "Usage: 02_summarize_leave_one_out_cv.R <scenario_root> ",
    "<C_root_only|D_papilionidae_only>",
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
if (!requireNamespace("ape", quietly = TRUE)) {
  stop("R package 'ape' is required", call. = FALSE)
}

config_dir <- file.path(scenario_root, "configs")
log_dir <- file.path(scenario_root, "logs")
output_dir <- file.path(scenario_root, "output")
qa_dir <- file.path(scenario_root, "qa")
if (!all(dir.exists(c(config_dir, log_dir, output_dir, qa_dir)))) {
  stop("Scenario is missing one or more required directories", call. = FALSE)
}

out_long <- file.path(qa_dir, "cv_scores_long.tsv")
out_summary <- file.path(qa_dir, "cv_summary.tsv")
out_decision <- file.path(qa_dir, "cv_decision.tsv")
out_report <- file.path(qa_dir, "cv_report.md")
output_paths <- c(out_long, out_summary, out_decision, out_report)
existing_outputs <- output_paths[file.exists(output_paths)]
if (length(existing_outputs)) {
  stop(
    "Refusing to overwrite existing output(s): ",
    paste(existing_outputs, collapse = ", "), call. = FALSE
  )
}

input_tree_expected <- paste0(
  "/home/data/t200301/xjt/Hesperiinae_review/",
  "treePL_R2_3_T25_deep_only_20260904_000707/",
  "input/verified_495_tip_input.treefile"
)
input_tree_sha256_expected <-
  "4d35d553ed5bbf08022daef2361ffa0c342f8694a53aaf08d5350692094439e5"
binary_expected <- paste0(
  "/home/data/t200301/xjt/Hesperiinae_review/",
  "treePL_R2_3_T25_deep_only_20260904_000707/",
  "software/treePL_randomcv_fix_src/treePL"
)
binary_sha256_expected <-
  "5d6bfac8d18d835f24ebb79f0373b070a948b58a179f437f0aacb8fed74ccc5d"

stage_ids <- sprintf("cv_initial_run_%02d", 1:3)
expected_seeds <- as.character(2026090510 + 1:3)
formal_grid <- c(
  1e5, 1e4, 1e3, 1e2, 1e1, 1, 1e-1, 1e-2, 1e-3, 1e-4,
  1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13
)
formal_labels <- c(
  "100000", "10000", "1000", "100", "10", "1", "0.1", "0.01",
  "0.001", "0.0001", "0.00001", "0.000001", "0.0000001",
  "0.00000001", "0.000000001", "0.0000000001",
  "0.00000000001", "0.000000000001", "0.0000000000001"
)
tie_tolerance <- 1e-10
stdout_score_tolerance <- 5e-5

calibration_lines_by_scenario <- list(
  C_root_only = c(
    paste(
      "mrca = CROWN_PAPILIONOIDEA",
      "Parnassius_apollo_ncbi Aeromachus_catocyanea_mydata"
    ),
    "min = CROWN_PAPILIONOIDEA 91.5046",
    "max = CROWN_PAPILIONOIDEA 100.8925"
  ),
  D_papilionidae_only = c(
    paste(
      "mrca = CROWN_PAPILIONIDAE_NONBARONIA",
      "Parnassius_apollo_ncbi Graphium_cloanthus_mydata"
    ),
    "min = CROWN_PAPILIONIDAE_NONBARONIA 44.1968",
    "max = CROWN_PAPILIONIDAE_NONBARONIA 52.9473"
  )
)
expected_calibration_lines <- calibration_lines_by_scenario[[scenario]]

canonical <- function(path, must_work = FALSE) {
  normalizePath(path, mustWork = must_work)
}

relative_difference <- function(a, b) {
  abs(a - b) / pmax(abs(a), abs(b), .Machine$double.xmin)
}

numeric_equal <- function(a, b, tolerance = 1e-12) {
  if (length(a) != length(b) || any(!is.finite(a)) || any(!is.finite(b))) {
    return(FALSE)
  }
  scale <- pmax(abs(a), abs(b), .Machine$double.xmin)
  all(abs(a - b) <= tolerance * scale)
}

read_active_config <- function(path) {
  if (!file.exists(path)) stop("Missing config: ", path, call. = FALSE)
  lines <- readLines(path, warn = FALSE)
  active <- trimws(sub("#.*$", "", lines))
  active[nzchar(active)]
}

get_one_value <- function(active, key, path) {
  pattern <- paste0("^", key, "[[:space:]]*=[[:space:]]*(.+?)[[:space:]]*$")
  hits <- grep(pattern, active, perl = TRUE)
  if (length(hits) != 1L) {
    stop("Expected exactly one ", key, " assignment in ", path, call. = FALSE)
  }
  sub(pattern, "\\1", active[hits], perl = TRUE)
}

sha256_file <- function(path) {
  if (!file.exists(path) || file.info(path)$isdir) {
    stop("Missing file for SHA-256: ", path, call. = FALSE)
  }
  sha256sum <- Sys.which("sha256sum")
  if (!nzchar(sha256sum)) stop("sha256sum executable not found", call. = FALSE)
  result <- system2(sha256sum, path, stdout = TRUE, stderr = TRUE)
  status <- attr(result, "status")
  if (!is.null(status) && status != 0L) {
    stop("sha256sum failed for ", path, call. = FALSE)
  }
  fields <- strsplit(result[[1L]], "[[:space:]]+", perl = TRUE)[[1L]]
  hash <- fields[nzchar(fields)][[1L]]
  if (!grepl("^[0-9a-f]{64}$", hash)) {
    stop("Malformed SHA-256 result for ", path, call. = FALSE)
  }
  hash
}

validate_fixed_inputs <- function() {
  observed_tree <- canonical(input_tree_expected, must_work = TRUE)
  observed_binary <- canonical(binary_expected, must_work = TRUE)
  if (!identical(sha256_file(observed_tree), input_tree_sha256_expected)) {
    stop("Fixed input-tree SHA-256 mismatch", call. = FALSE)
  }
  if (!identical(sha256_file(observed_binary), binary_sha256_expected)) {
    stop("Fixed treePL-binary SHA-256 mismatch", call. = FALSE)
  }
  invisible(TRUE)
}

read_prime_recommendations <- function() {
  cfg_path <- file.path(qa_dir, "prime_recommendations.cfg")
  tsv_path <- file.path(qa_dir, "prime_recommendations.tsv")
  cfg_lines <- read_active_config(cfg_path)
  if (length(cfg_lines) != 6L) {
    stop("prime_recommendations.cfg must contain exactly six active lines")
  }
  expected_keys <- c(
    "opt", "moredetail", "optad", "moredetailad", "optcvad", "moredetailcvad"
  )
  observed_keys <- ifelse(
    grepl("=", cfg_lines, fixed = TRUE),
    trimws(sub("=.*$", "", cfg_lines)),
    cfg_lines
  )
  if (!identical(observed_keys, expected_keys)) {
    stop("Prime recommendation directives are not the exact required six-line set")
  }
  integer_rows <- c(1L, 3L, 5L)
  detail_rows <- c(2L, 4L, 6L)
  integer_text <- vapply(
    cfg_lines[integer_rows],
    function(x) trimws(sub("^[^=]+=[[:space:]]*", "", x)),
    character(1)
  )
  integer_values <- suppressWarnings(as.numeric(integer_text))
  if (any(!grepl("^[1-9][0-9]*$", integer_text)) ||
      any(!is.finite(integer_values)) || any(integer_values < 1) ||
      any(integer_values > .Machine$integer.max)) {
    stop("Prime integer recommendations must be positive machine-range integers")
  }
  if (!identical(cfg_lines[detail_rows], expected_keys[detail_rows])) {
    stop("Prime detail recommendations must be enabled as bare directives")
  }

  if (!file.exists(tsv_path) || file.info(tsv_path)$size <= 0) {
    stop("Missing/empty prime recommendations TSV: ", tsv_path, call. = FALSE)
  }
  tab <- read.delim(
    tsv_path, stringsAsFactors = FALSE, check.names = FALSE,
    colClasses = "character"
  )
  if (!identical(names(tab), c("parameter", "recommendation")) || nrow(tab) != 6L ||
      !identical(tab$parameter, expected_keys)) {
    stop("prime_recommendations.tsv must contain the exact six ordered rows")
  }
  reconstructed <- ifelse(
    tab$parameter %in% expected_keys[detail_rows],
    ifelse(tab$recommendation == "enabled", tab$parameter, NA_character_),
    paste(tab$parameter, "=", tab$recommendation)
  )
  if (anyNA(reconstructed) || !identical(reconstructed, cfg_lines)) {
    stop("Prime CFG and TSV recommendations disagree", call. = FALSE)
  }
  cfg_lines
}

validate_config <- function(stage, replicate, prime_recommendations) {
  config_path <- file.path(config_dir, paste0(stage, ".cfg"))
  active <- read_active_config(config_path)
  score_path <- file.path(qa_dir, paste0(stage, "_scores.txt"))
  dated_path <- file.path(output_dir, paste0(stage, "_dated.tre"))

  calibration_lines <- active[grepl("^(mrca|min|max)[[:space:]]*=", active)]
  if (!identical(calibration_lines, expected_calibration_lines)) {
    stop(
      "Calibration block is not exactly the single prescribed calibration in ",
      config_path, call. = FALSE
    )
  }
  if (sum(grepl("^mrca[[:space:]]*=", active)) != 1L ||
      sum(grepl("^min[[:space:]]*=", active)) != 1L ||
      sum(grepl("^max[[:space:]]*=", active)) != 1L) {
    stop("Config must contain exactly one mrca/min/max block: ", config_path)
  }

  optimization_lines <- active[grepl("^(opt|moredetail)", active)]
  if (!identical(optimization_lines, prime_recommendations)) {
    stop("Config does not exactly use the scenario prime recommendations: ", config_path)
  }
  if (sum(active == "randomcv") != 1L ||
      any(active %in% c("cv", "prime", "thorough")) ||
      any(grepl("^smooth[[:space:]]*=", active))) {
    stop("Invalid or incompatible treePL mode directive in ", config_path)
  }

  expected_text <- list(
    treefile = input_tree_expected,
    numsites = "182682",
    outfile = dated_path,
    cvoutfile = score_path,
    cvstart = "100000",
    cvstop = "0.0000000000001",
    cvmultstep = "0.1",
    cviter = "3",
    nthreads = "10",
    seed = expected_seeds[[replicate]]
  )
  path_keys <- c("treefile", "outfile", "cvoutfile")
  for (key in names(expected_text)) {
    observed <- get_one_value(active, key, config_path)
    expected <- expected_text[[key]]
    if (key %in% path_keys) {
      observed <- canonical(observed, must_work = identical(key, "treefile"))
      expected <- canonical(expected, must_work = identical(key, "treefile"))
    }
    if (!identical(observed, expected)) {
      stop(
        key, " mismatch in ", config_path, ": observed=", observed,
        " expected=", expected, call. = FALSE
      )
    }
  }
  list(config = config_path, score = score_path, dated = dated_path)
}

read_metadata <- function(stage, config_path) {
  path <- file.path(log_dir, paste0(stage, ".metadata.tsv"))
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
      stop("Metadata key is missing or non-unique: ", key, " in ", path)
    }
    value
  }
  if (!identical(one("stage"), stage)) stop("Metadata stage mismatch: ", path)
  if (!identical(one("exit_code"), "0")) stop("Nonzero treePL exit: ", path)
  if (!identical(canonical(one("binary"), TRUE), canonical(binary_expected, TRUE))) {
    stop("Metadata binary mismatch: ", path)
  }
  if (!identical(canonical(one("config"), TRUE), canonical(config_path, TRUE))) {
    stop("Metadata config mismatch: ", path)
  }
  invisible(TRUE)
}

parse_scores <- function(path, require_all_nonblank = FALSE) {
  if (!file.exists(path) || file.info(path)$size <= 0) {
    stop("Missing/empty score-bearing file: ", path, call. = FALSE)
  }
  lines <- readLines(path, warn = FALSE)
  matches <- regexec(
    "^chisq:[[:space:]]*\\(([^)]+)\\)[[:space:]]+([^[:space:]]+)[[:space:]]*$",
    lines, perl = TRUE
  )
  parts <- regmatches(lines, matches)
  parts <- parts[lengths(parts) == 3L]
  if (require_all_nonblank && length(parts) != sum(nzchar(trimws(lines)))) {
    stop("Unexpected non-score content in ", path, call. = FALSE)
  }
  if (length(parts) != length(formal_grid)) {
    stop(
      "Expected exactly 19 score rows in ", path,
      "; observed ", length(parts), call. = FALSE
    )
  }
  smoothing_text <- vapply(parts, function(x) x[[2L]], character(1))
  score_text <- vapply(parts, function(x) x[[3L]], character(1))
  smoothing <- suppressWarnings(as.numeric(smoothing_text))
  raw_score <- suppressWarnings(as.numeric(score_text))
  if (!numeric_equal(smoothing, formal_grid) || any(!is.finite(raw_score))) {
    stop("Malformed, non-finite, duplicated, or out-of-order score grid in ", path)
  }
  data.frame(
    smoothing = smoothing,
    smoothing_text = smoothing_text,
    raw_score = raw_score,
    stringsAsFactors = FALSE
  )
}

descendant_clade_signatures <- function(tree) {
  n_tip <- length(tree$tip.label)
  children <- split(tree$edge[, 2L], tree$edge[, 1L])
  memo <- new.env(parent = emptyenv())
  descendants <- function(node) {
    key <- as.character(node)
    if (exists(key, envir = memo, inherits = FALSE)) {
      return(get(key, envir = memo, inherits = FALSE))
    }
    result <- if (node <= n_tip) {
      tree$tip.label[[node]]
    } else {
      kids <- children[[key]]
      if (is.null(kids) || length(kids) < 2L) {
        stop("Malformed internal node in phylogeny", call. = FALSE)
      }
      sort(unique(unlist(lapply(kids, descendants), use.names = FALSE)))
    }
    assign(key, result, envir = memo)
    result
  }
  internal_nodes <- seq.int(n_tip + 1L, n_tip + tree$Nnode)
  sort(vapply(
    internal_nodes,
    function(node) paste(descendants(node), collapse = "\037"),
    character(1)
  ))
}

read_strict_tree <- function(path, label) {
  if (!file.exists(path) || file.info(path)$size <= 0) {
    stop("Missing/empty ", label, ": ", path, call. = FALSE)
  }
  tree <- tryCatch(
    ape::read.tree(path),
    error = function(e) stop("Cannot parse ", label, ": ", conditionMessage(e))
  )
  if (!inherits(tree, "phylo") || length(tree$tip.label) != 495L ||
      anyNA(tree$tip.label) || any(!nzchar(tree$tip.label)) ||
      anyDuplicated(tree$tip.label)) {
    stop(label, " must contain exactly 495 unique, nonempty tip labels")
  }
  if (!ape::is.rooted(tree) || !ape::is.binary.tree(tree) || tree$Nnode != 494L) {
    stop(label, " must be a rooted, fully binary 495-tip tree")
  }
  tree
}

validate_dated_topology <- function(path, input_tree, input_signatures) {
  dated <- read_strict_tree(path, "dated output tree")
  if (!setequal(dated$tip.label, input_tree$tip.label)) {
    stop("Dated output tip-label set differs from the fixed input: ", path)
  }
  if (!identical(descendant_clade_signatures(dated), input_signatures)) {
    stop("Dated output rooted topology differs from the fixed input: ", path)
  }
  invisible(TRUE)
}

validate_resource_log <- function(path) {
  if (!file.exists(path) || file.info(path)$size <= 0) {
    stop("Missing/empty resource log: ", path, call. = FALSE)
  }
  lines <- readLines(path, warn = FALSE)
  exit_lines <- grep("^[[:space:]]*Exit status:[[:space:]]*", lines, value = TRUE)
  if (length(exit_lines) != 1L ||
      !grepl("^[[:space:]]*Exit status:[[:space:]]*0[[:space:]]*$", exit_lines)) {
    stop("Resource log does not contain exactly one Exit status: 0: ", path)
  }
  invisible(TRUE)
}

validate_stage <- function(stage, replicate, prime_recommendations,
                           input_tree, input_signatures) {
  paths <- validate_config(stage, replicate, prime_recommendations)
  read_metadata(stage, paths$config)

  stderr_path <- file.path(log_dir, paste0(stage, ".stderr"))
  stdout_path <- file.path(log_dir, paste0(stage, ".stdout"))
  resource_path <- file.path(log_dir, paste0(stage, ".resources.txt"))
  if (!file.exists(stderr_path) || file.info(stderr_path)$size != 0) {
    stop("Missing or non-empty stderr: ", stderr_path, call. = FALSE)
  }
  if (!file.exists(stdout_path) || file.info(stdout_path)$size <= 0) {
    stop("Missing/empty stdout: ", stdout_path, call. = FALSE)
  }
  validate_resource_log(resource_path)
  validate_dated_topology(paths$dated, input_tree, input_signatures)

  stdout_lines <- readLines(stdout_path, warn = FALSE)
  failure_re <- paste(
    c(
      "problem", "complete failure", "segmentation",
      "(^|[^[:alpha:]])nan([^[:alpha:]]|$)",
      "(^|[^[:alpha:]])inf([^[:alpha:]]|$)",
      "(^|[^[:alpha:]])error([^[:alpha:]]|$)"
    ),
    collapse = "|"
  )
  if (any(grepl(failure_re, stdout_lines, ignore.case = TRUE, perl = TRUE))) {
    stop("Matched treePL failure token in ", stdout_path, call. = FALSE)
  }

  score_file <- parse_scores(paths$score, require_all_nonblank = TRUE)
  stdout_scores <- parse_scores(stdout_path, require_all_nonblank = FALSE)
  if (!numeric_equal(score_file$smoothing, stdout_scores$smoothing) ||
      any(relative_difference(score_file$raw_score, stdout_scores$raw_score) >
          stdout_score_tolerance)) {
    stop("Score-file/stdout mismatch in ", stage, call. = FALSE)
  }

  score_file$scenario <- scenario
  score_file$stage <- stage
  score_file$replicate <- replicate
  score_file$seed <- expected_seeds[[replicate]]
  score_file
}

validate_fixed_inputs()
prime_recommendations <- read_prime_recommendations()
input_tree <- read_strict_tree(input_tree_expected, "fixed input tree")
input_signatures <- descendant_clade_signatures(input_tree)

formal <- do.call(rbind, lapply(seq_along(stage_ids), function(i) {
  validate_stage(
    stage_ids[[i]], i, prime_recommendations, input_tree, input_signatures
  )
}))
rownames(formal) <- NULL
if (nrow(formal) != 57L) stop("Formal score table must contain exactly 57 rows")
for (i in seq_along(stage_ids)) {
  rows <- formal[formal$replicate == i, , drop = FALSE]
  if (nrow(rows) != 19L || !numeric_equal(rows$smoothing, formal_grid)) {
    stop("Replicate does not contain the exact ordered 19-point grid: ", i)
  }
}

formal$smoothing_label <- formal_labels[match(formal$smoothing, formal_grid)]
formal$replicate_rank <- ave(
  formal$raw_score, formal$replicate,
  FUN = function(x) rank(x, ties.method = "min")
)
formal$replicate_winner <- ave(
  formal$raw_score, formal$replicate,
  FUN = function(x) relative_difference(x, min(x)) <= tie_tolerance
) == 1

summary_tab <- do.call(rbind, lapply(seq_along(formal_grid), function(i) {
  rows <- formal[formal$smoothing == formal_grid[[i]], , drop = FALSE]
  data.frame(
    scenario = scenario,
    smoothing = formal_grid[[i]],
    smoothing_label = formal_labels[[i]],
    median_score = median(rows$raw_score),
    mean_score = mean(rows$raw_score),
    min_score = min(rows$raw_score),
    max_score = max(rows$raw_score),
    score_range = diff(range(rows$raw_score)),
    replicate_winner_count = sum(rows$replicate_winner),
    stringsAsFactors = FALSE
  )
}))
summary_tab$aggregate_rank <- rank(summary_tab$median_score, ties.method = "min")
aggregate_min <- min(summary_tab$median_score)
summary_tab$aggregate_winner <-
  relative_difference(summary_tab$median_score, aggregate_min) <= tie_tolerance

aggregate_winners <- summary_tab$smoothing[summary_tab$aggregate_winner]
if (!length(aggregate_winners)) stop("No aggregate CV winner", call. = FALSE)
lower_boundary_hit <- min(formal_grid) %in% aggregate_winners
upper_boundary_hit <- max(formal_grid) %in% aggregate_winners
boundary_hit <- lower_boundary_hit || upper_boundary_hit
extension_direction <- if (lower_boundary_hit && upper_boundary_hit) {
  "lower_and_upper_smoothing"
} else if (lower_boundary_hit) {
  "lower_smoothing"
} else if (upper_boundary_hit) {
  "upper_smoothing"
} else {
  "none"
}

replicate_winner_sets <- lapply(seq_along(stage_ids), function(i) {
  rows <- formal[formal$replicate == i, , drop = FALSE]
  rows$smoothing[rows$replicate_winner]
})
distinct_replicate_winners <- formal_grid[
  formal_grid %in% unique(unlist(replicate_winner_sets, use.names = FALSE))
]
selected <- if (boundary_hit) NA_real_ else max(aggregate_winners)
replicate_consistent <- !boundary_hit && all(vapply(
  replicate_winner_sets,
  function(values) {
    length(values) == 1L && numeric_equal(values, selected)
  },
  logical(1)
))
aggregate_tie <- length(aggregate_winners) > 1L

decision_status <- if (boundary_hit) {
  "REQUIRES_GRID_EXTENSION"
} else if (aggregate_tie || !replicate_consistent) {
  "ACCEPTED_WITH_SENSITIVITIES"
} else {
  "ACCEPTED"
}

label_for <- function(values) {
  if (!length(values)) return(character())
  formal_labels[match(values, formal_grid)]
}
format_values <- function(values) {
  labels <- label_for(values)
  if (!length(labels)) "none" else paste(labels, collapse = ",")
}

other_aggregate_winners <- if (is.na(selected)) {
  aggregate_winners
} else {
  aggregate_winners[aggregate_winners != selected]
}
sensitivity_values <- if (is.na(selected)) {
  unique(c(aggregate_winners, distinct_replicate_winners))
} else {
  unique(c(other_aggregate_winners, distinct_replicate_winners[
    distinct_replicate_winners != selected
  ]))
}
sensitivity_values <- formal_grid[formal_grid %in% sensitivity_values]

formal$aggregate_rank <- summary_tab$aggregate_rank[
  match(formal$smoothing, summary_tab$smoothing)
]
formal$aggregate_winner <- summary_tab$aggregate_winner[
  match(formal$smoothing, summary_tab$smoothing)
]
formal <- formal[, c(
  "scenario", "stage", "replicate", "seed", "smoothing",
  "smoothing_label", "raw_score", "replicate_rank", "replicate_winner",
  "aggregate_rank", "aggregate_winner"
)]

replicate_winner_text <- vapply(
  replicate_winner_sets, format_values, character(1)
)
decision <- data.frame(
  scenario = scenario,
  status = decision_status,
  decision_basis = "lowest_median_raw_chisq_across_three_19_point_randomCV_runs",
  selected_smoothing = if (is.na(selected)) NA_character_ else format_values(selected),
  selected_median_score = if (is.na(selected)) NA_real_ else
    summary_tab$median_score[summary_tab$smoothing == selected],
  grid_extension_direction = extension_direction,
  aggregate_winner_at_boundary = boundary_hit,
  aggregate_winner_at_lower_boundary = lower_boundary_hit,
  aggregate_winner_at_upper_boundary = upper_boundary_hit,
  aggregate_winner_count = length(aggregate_winners),
  aggregate_winners = format_values(aggregate_winners),
  other_aggregate_winners = format_values(other_aggregate_winners),
  aggregate_tie = aggregate_tie,
  replicate_winners_run01 = replicate_winner_text[[1L]],
  replicate_winners_run02 = replicate_winner_text[[2L]],
  replicate_winners_run03 = replicate_winner_text[[3L]],
  distinct_replicate_winners = format_values(distinct_replicate_winners),
  replicate_winners_consistent_with_selection = replicate_consistent,
  sensitivity_required = !identical(decision_status, "ACCEPTED"),
  sensitivity_smoothing_values = format_values(sensitivity_values),
  formal_runs = paste(stage_ids, collapse = ","),
  formal_grid_points_per_run = length(formal_grid),
  formal_rows = nrow(formal),
  formal_seed_values = paste(expected_seeds, collapse = ","),
  tie_relative_tolerance = tie_tolerance,
  stdout_score_relative_tolerance = stdout_score_tolerance,
  config_validation_passed = TRUE,
  prime_six_line_validation_passed = TRUE,
  metadata_validation_passed = TRUE,
  stderr_validation_passed = TRUE,
  stdout_and_score_validation_passed = TRUE,
  resource_exit_zero_validation_passed = TRUE,
  dated_tree_495_tip_rooted_topology_validation_passed = TRUE,
  input_tree_sha256 = input_tree_sha256_expected,
  treepl_binary_sha256 = binary_sha256_expected,
  stringsAsFactors = FALSE
)

report_lines <- c(
  paste0("# ", scenario, " initial random-CV decision"),
  "",
  paste0("- Status: `", decision_status, "`"),
  paste0("- Aggregate tied winner(s): `", format_values(aggregate_winners), "`"),
  paste0(
    "- Selected smoothing: `",
    if (is.na(selected)) "NA" else format_values(selected), "`"
  ),
  paste0("- Required grid-extension direction: `", extension_direction, "`"),
  paste0(
    "- Replicate winner(s), runs 01/02/03: `",
    paste(replicate_winner_text, collapse = " ; "), "`"
  ),
  paste0("- Other aggregate tied winner(s): `", format_values(other_aggregate_winners), "`"),
  paste0("- Distinct replicate winner(s): `", format_values(distinct_replicate_winners), "`"),
  paste0("- Sensitivity smoothing value(s): `", format_values(sensitivity_values), "`"),
  "",
  "The decision uses the median raw chi-square score at each smoothing value",
  "across three complete 19-point random-CV runs. Lower scores are better;",
  "relative differences no greater than `1e-10` are treated as ties.",
  "",
  if (boundary_hit) {
    paste0(
      "Because the tied aggregate winner set includes a grid boundary, no ",
      "smoothing value is selected from this grid. Extend once toward `",
      extension_direction, "` as prescribed in `analysis_design.md`."
    )
  } else if (identical(decision_status, "ACCEPTED_WITH_SENSITIVITIES")) {
    paste0(
      "The largest interior aggregate tied winner is selected. The listed ",
      "alternative aggregate/replicate winners must be retained as sensitivity values."
    )
  } else {
    "The unique interior aggregate winner is accepted and all replicates agree."
  },
  "",
  "Before scoring, the script verified the exact single calibration block, fixed",
  "input tree and `numsites`, exact CV grid/settings/seeds, six scenario-specific",
  "prime recommendations, metadata/config/binary/exit status, empty stderr,",
  "stdout/score agreement, resource exit status, and 495-tip rooted topology of",
  "all three dated outputs."
)

tmp_dir <- tempfile(pattern = ".cv_summary_", tmpdir = qa_dir)
if (!dir.create(tmp_dir, mode = "0700")) {
  stop("Could not create temporary output directory: ", tmp_dir)
}
on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)
tmp_paths <- file.path(tmp_dir, basename(output_paths))
write.table(
  formal, tmp_paths[[1L]], sep = "\t", quote = FALSE, row.names = FALSE,
  na = "NA"
)
write.table(
  summary_tab, tmp_paths[[2L]], sep = "\t", quote = FALSE, row.names = FALSE,
  na = "NA"
)
write.table(
  decision, tmp_paths[[3L]], sep = "\t", quote = FALSE, row.names = FALSE,
  na = "NA"
)
writeLines(report_lines, tmp_paths[[4L]], useBytes = TRUE)

if (any(file.exists(output_paths))) {
  stop("Refusing to overwrite outputs created during validation", call. = FALSE)
}
for (i in seq_along(output_paths)) {
  if (!file.rename(tmp_paths[[i]], output_paths[[i]])) {
    stop("Could not atomically install output: ", output_paths[[i]], call. = FALSE)
  }
}

cat("Scenario:", scenario, "\n")
cat("Formal random-CV rows:", nrow(formal), "\n")
cat("Aggregate winner(s):", format_values(aggregate_winners), "\n")
cat("Selected smoothing:", if (is.na(selected)) "NA" else format_values(selected), "\n")
cat("Grid-extension direction:", extension_direction, "\n")
cat("Status:", decision_status, "\n")

if (identical(decision_status, "REQUIRES_GRID_EXTENSION")) {
  quit(status = 2L)
}
