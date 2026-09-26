#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: 06b_summarize_cv_23point.R <run_root>", call. = FALSE)
}
if (!requireNamespace("ape", quietly = TRUE)) {
  stop("R package 'ape' is required for strict dated-tree validation", call. = FALSE)
}

run_root <- normalizePath(args[[1]], mustWork = TRUE)
qa_dir <- file.path(run_root, "qa")
config_dir <- file.path(run_root, "configs")
log_dir <- file.path(run_root, "logs")
output_dir <- file.path(run_root, "output")
binary_expected <- file.path(
  run_root, "software", "treePL_randomcv_fix_src", "treePL"
)
tree_expected <- file.path(run_root, "input", "verified_495_tip_input.treefile")

out_long <- file.path(qa_dir, "cv_23point_scores_long.tsv")
out_summary <- file.path(qa_dir, "cv_23point_summary.tsv")
out_decision <- file.path(qa_dir, "cv_23point_decision.tsv")
for (path in c(out_long, out_summary, out_decision)) {
  if (file.exists(path)) {
    stop("Refusing to overwrite existing output: ", path, call. = FALSE)
  }
}

declarations <- file.path(
  qa_dir,
  c(
    "cv_selection_rule.md",
    "cv_grid_extension_decision.md",
    "cv_continuous_grid_decision.md",
    "cv_19point_boundary_extension_decision.md",
    "cv_tie_plateau_rule_declared_before_low_extension_scores.md"
  )
)
missing_declarations <- declarations[!file.exists(declarations)]
if (length(missing_declarations)) {
  stop(
    "Missing prespecified decision document(s): ",
    paste(basename(missing_declarations), collapse = ", "),
    call. = FALSE
  )
}
empty_declarations <- declarations[file.info(declarations)$size <= 0]
if (length(empty_declarations)) {
  stop(
    "Empty prespecified decision document(s): ",
    paste(basename(empty_declarations), collapse = ", "),
    call. = FALSE
  )
}

formal_grid <- 10^(5:-17)
formal_labels <- c(
  "100000", "10000", "1000", "100", "10", "1", "0.1", "0.01",
  "0.001", "0.0001", "0.00001", "0.000001", "0.0000001",
  "0.00000001", "0.000000001", "0.0000000001",
  "0.00000000001", "0.000000000001", "0.0000000000001",
  "0.00000000000001", "0.000000000000001",
  "0.0000000000000001", "0.00000000000000001"
)
diagnostic_grid <- formal_grid[1:19]
formal_ids <- sprintf("cv_extendedgrid_run_%02d", 1:3)
diagnostic_ids <- sprintf("cv_finalgrid_run_%02d", 1:3)
expected_seeds <- as.character(2026090440 + 1:3)
tie_tolerance <- 1e-10
stdout_score_tolerance <- 5e-5
ultrametric_tolerance_ma <- 1e-4

if (length(formal_grid) != 23L || length(formal_labels) != 23L) {
  stop("Internal 23-point grid definition error", call. = FALSE)
}
if (!identical(formal_grid, 10^(5:-17))) {
  stop("Formal grid is not the exact continuous 1e5-through-1e-17 grid")
}

canonical <- function(path) {
  normalizePath(path, mustWork = FALSE)
}

relative_difference <- function(a, b) {
  abs(a - b) / pmax(1, abs(a), abs(b))
}

numeric_equal <- function(a, b, tolerance = 1e-12) {
  if (length(a) != length(b) || any(!is.finite(a)) || any(!is.finite(b))) {
    return(FALSE)
  }
  scale <- pmax(abs(a), abs(b), .Machine$double.xmin)
  all(abs(a - b) <= tolerance * scale)
}

read_active_config <- function(path) {
  lines <- readLines(path, warn = FALSE)
  active <- trimws(sub("#.*$", "", lines))
  active[nzchar(active)]
}

prime_recommendations_path <- file.path(qa_dir, "prime_recommendations.cfg")
if (!file.exists(prime_recommendations_path)) {
  stop("Missing prime recommendations: ", prime_recommendations_path, call. = FALSE)
}
prime_recommendations <- read_active_config(prime_recommendations_path)
expected_prime_recommendations <- c(
  "opt = 4",
  "moredetail",
  "optad = 4",
  "moredetailad",
  "optcvad = 1",
  "moredetailcvad"
)
if (!identical(prime_recommendations, expected_prime_recommendations)) {
  stop(
    "Prime recommendations are not the exact three-calibration recommendations",
    call. = FALSE
  )
}

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

reference_tree <- tryCatch(
  ape::read.tree(tree_expected),
  error = function(e) {
    stop("Cannot parse reference tree: ", conditionMessage(e), call. = FALSE)
  }
)
if (!inherits(reference_tree, "phylo") ||
    ape::Ntip(reference_tree) != 495L ||
    anyDuplicated(reference_tree$tip.label)) {
  stop("Reference tree is not the expected unique 495-tip phylogeny", call. = FALSE)
}
reference_topology <- reference_tree
reference_topology$edge.length <- NULL
reference_topology$node.label <- NULL

validate_source_19point_decision <- function() {
  path <- file.path(qa_dir, "cv_decision.tsv")
  if (!file.exists(path) || file.info(path)$size <= 0) {
    stop("Missing/empty source 19-point decision: ", path, call. = FALSE)
  }
  tab <- read.delim(
    path, stringsAsFactors = FALSE, check.names = FALSE,
    na.strings = c("NA", "")
  )
  required <- c(
    "status", "decision_basis", "selected_smoothing",
    "aggregate_winner_at_grid_boundary", "aggregate_winner_count",
    "formal_runs", "formal_grid_points_per_run", "formal_rows",
    "formal_seed_values"
  )
  if (nrow(tab) != 1L || any(!required %in% names(tab))) {
    stop("Malformed source 19-point decision: ", path, call. = FALSE)
  }
  value <- function(key) tab[[key]][[1L]]
  is_true <- function(x) !is.na(x) && identical(toupper(as.character(x)), "TRUE")
  selected_missing <- is.na(value("selected_smoothing")) ||
    !nzchar(as.character(value("selected_smoothing")))
  expected_runs <- paste(diagnostic_ids, collapse = ",")
  expected_seed_text <- paste(expected_seeds, collapse = ",")
  if (!identical(as.character(value("status")), "REQUIRES_GRID_EXTENSION") ||
      !identical(
        as.character(value("decision_basis")),
        "lowest_median_raw_chisq_across_three_19_point_randomCV_runs"
      ) ||
      !selected_missing ||
      !is_true(value("aggregate_winner_at_grid_boundary")) ||
      !identical(as.character(value("formal_runs")), expected_runs) ||
      suppressWarnings(as.integer(value("formal_grid_points_per_run"))) != 19L ||
      suppressWarnings(as.integer(value("formal_rows"))) != 57L ||
      !identical(as.character(value("formal_seed_values")), expected_seed_text)) {
    stop(
      "Source 19-point decision does not document the required boundary extension",
      call. = FALSE
    )
  }

  summary_path <- file.path(qa_dir, "cv_summary.tsv")
  if (!file.exists(summary_path) || file.info(summary_path)$size <= 0) {
    stop("Missing/empty source 19-point summary: ", summary_path, call. = FALSE)
  }
  source_summary <- read.delim(
    summary_path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    colClasses = "character"
  )
  required_summary <- c(
    "smoothing", "smoothing_label", "median_score", "aggregate_winner"
  )
  if (nrow(source_summary) != 19L ||
      any(!required_summary %in% names(source_summary))) {
    stop("Malformed source 19-point summary: ", summary_path, call. = FALSE)
  }
  source_smoothing <- suppressWarnings(as.numeric(source_summary$smoothing))
  source_median <- suppressWarnings(as.numeric(source_summary$median_score))
  winner_flag_text <- toupper(as.character(source_summary$aggregate_winner))
  if (!numeric_equal(source_smoothing, diagnostic_grid) ||
      !identical(
        as.character(source_summary$smoothing_label), formal_labels[1:19]
      ) ||
      any(!is.finite(source_median)) ||
      any(!winner_flag_text %in% c("TRUE", "FALSE"))) {
    stop("Invalid source 19-point summary grid/scores/flags", call. = FALSE)
  }
  source_winner_flags <- winner_flag_text == "TRUE"
  recomputed_winner_flags <- relative_difference(
    source_median, min(source_median)
  ) <= tie_tolerance
  if (!identical(source_winner_flags, recomputed_winner_flags)) {
    stop(
      "Source 19-point aggregate-winner flags disagree with recomputed medians",
      call. = FALSE
    )
  }
  source_winners <- source_smoothing[recomputed_winner_flags]
  source_winner_labels <- formal_labels[1:19][recomputed_winner_flags]
  lower_boundary_in_winners <- diagnostic_grid[[19L]] %in% source_winners
  upper_boundary_is_only_winner <- length(source_winners) == 1L &&
    diagnostic_grid[[1L]] %in% source_winners
  if (!lower_boundary_in_winners || upper_boundary_is_only_winner ||
      suppressWarnings(as.integer(value("aggregate_winner_count"))) !=
        length(source_winners)) {
    stop(
      "Source 19-point evidence does not prove a lower-boundary aggregate winner",
      call. = FALSE
    )
  }

  list(
    lower_boundary_in_aggregate_winner_set = TRUE,
    aggregate_winner_values = source_winner_labels
  )
}

validate_lint <- function(stage, expected_stop_text) {
  path <- file.path(qa_dir, paste0("lint_", stage, ".strict.tsv"))
  if (!file.exists(path) || file.info(path)$size <= 0) {
    stop("Missing/empty config lint: ", path, call. = FALSE)
  }
  lint <- read.delim(
    path, stringsAsFactors = FALSE, check.names = FALSE,
    quote = "", comment.char = ""
  )
  if (!all(c("check", "pass", "observed", "expected") %in% names(lint)) ||
      nrow(lint) < 20L ||
      anyNA(lint$check) ||
      anyDuplicated(lint$check) ||
      !all(toupper(as.character(lint$pass)) == "TRUE")) {
    stop("Config lint is malformed or contains a failure: ", path, call. = FALSE)
  }
  required_checks <- c(
    "treefile_exact", "numsites_182682",
    "exactly_three_mrca", "exactly_three_min", "exactly_three_max",
    "root_mrca_exact", "root_bounds_exact",
    "papilionidae_mrca_exact", "papilionidae_bounds_exact",
    "hesperiinae_mrca_exact", "hesperiinae_bounds_exact",
    "log_pen_absent", "randomcv_enabled",
    "prime_disabled", "thorough_disabled_for_cv",
    "cv_grid_start", "cv_grid_stop", "cv_grid_multiplier",
    "cviter_3", "nthreads_10",
    "cvoutfile_within_run", "outfile_within_run"
  )
  if (any(!required_checks %in% lint$check)) {
    stop(
      "Config lint lacks required checks for ", stage, ": ",
      paste(required_checks[!required_checks %in% lint$check], collapse = ", "),
      call. = FALSE
    )
  }
  check_exact <- function(check, observed, expected) {
    row <- lint[lint$check == check, , drop = FALSE]
    if (nrow(row) != 1L ||
        !identical(as.character(row$observed[[1L]]), observed) ||
        !identical(as.character(row$expected[[1L]]), expected)) {
      stop("Unexpected lint evidence for ", check, " in ", path, call. = FALSE)
    }
  }
  check_exact("cv_grid_start", "100000", "100000")
  check_exact("cv_grid_stop", expected_stop_text, expected_stop_text)
  check_exact("cv_grid_multiplier", "0.1", "0.1")
  check_exact("cviter_3", "3", "3")
  check_exact("nthreads_10", "10", "10")
  invisible(path)
}

validate_config <- function(
    stage, expected_seed, expected_score, expected_dated_tree,
    expected_stop_text
) {
  path <- file.path(config_dir, paste0(stage, ".cfg"))
  if (!file.exists(path)) {
    stop("Missing config: ", path, call. = FALSE)
  }
  active <- read_active_config(path)
  expected_active <- c(
    paste0("treefile = ", tree_expected),
    "numsites = 182682",
    paste0("outfile = ", expected_dated_tree),
    expected_calibration_lines,
    "randomcv",
    "cvstart = 100000",
    paste0("cvstop = ", expected_stop_text),
    "cvmultstep = 0.1",
    "cviter = 3",
    paste0("cvoutfile = ", expected_score),
    "nthreads = 10",
    paste0("seed = ", expected_seed),
    expected_prime_recommendations
  )
  if (!identical(active, expected_active)) {
    stop(
      "Config is not the exact prescribed random-CV contract: ", path,
      call. = FALSE
    )
  }
  validate_lint(stage, expected_stop_text)
  invisible(path)
}

read_metadata <- function(stage, config_path) {
  path <- file.path(log_dir, paste0(stage, ".metadata.tsv"))
  if (!file.exists(path) || file.info(path)$size <= 0) {
    stop("Missing/empty metadata: ", path, call. = FALSE)
  }
  tab <- read.delim(
    path, header = FALSE, sep = "	", quote = "", comment.char = "",
    col.names = c("key", "value"), stringsAsFactors = FALSE
  )
  if (anyNA(tab$key) || anyDuplicated(tab$key)) {
    stop("Malformed/duplicated metadata keys: ", path, call. = FALSE)
  }
  one <- function(key) {
    value <- tab$value[tab$key == key]
    if (length(value) != 1L) {
      stop("Metadata key not unique: ", key, " in ", path, call. = FALSE)
    }
    value
  }
  required_keys <- c(
    "stage", "pid", "start_time", "binary", "config",
    "end_time", "exit_code", "elapsed_s"
  )
  invisible(lapply(required_keys, one))
  if (!identical(one("stage"), stage) ||
      !identical(one("exit_code"), "0") ||
      !identical(canonical(one("binary")), canonical(binary_expected)) ||
      !identical(canonical(one("config")), canonical(config_path))) {
    stop("Metadata contract mismatch: ", path, call. = FALSE)
  }
  pid <- suppressWarnings(as.numeric(one("pid")))
  elapsed <- suppressWarnings(as.numeric(one("elapsed_s")))
  if (!is.finite(pid) || pid <= 0 ||
      !is.finite(elapsed) || elapsed < 0 ||
      !nzchar(one("start_time")) || !nzchar(one("end_time"))) {
    stop("Invalid runtime fields in metadata: ", path, call. = FALSE)
  }
  if (!file.exists(binary_expected) ||
      file.info(binary_expected)$size <= 0 ||
      file.access(binary_expected, mode = 1) != 0) {
    stop("Expected treePL binary is missing, empty, or non-executable")
  }
  invisible(path)
}

parse_scores <- function(path, expected_grid, require_all_nonblank = FALSE) {
  if (!file.exists(path) || file.info(path)$size <= 0) {
    stop("Missing/empty score file: ", path, call. = FALSE)
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
  if (length(parts) != length(expected_grid)) {
    stop(
      "Expected ", length(expected_grid), " score rows in ", path,
      ", observed ", length(parts), call. = FALSE
    )
  }
  smoothing_text <- vapply(parts, function(x) x[[2L]], character(1))
  score_text <- vapply(parts, function(x) x[[3L]], character(1))
  smoothing <- suppressWarnings(as.numeric(smoothing_text))
  score <- suppressWarnings(as.numeric(score_text))
  if (!numeric_equal(smoothing, expected_grid) || any(!is.finite(score))) {
    stop(
      "Malformed, non-finite, duplicated, missing, or out-of-order score grid in ",
      path, call. = FALSE
    )
  }
  data.frame(
    smoothing = smoothing,
    smoothing_text = smoothing_text,
    raw_score = score,
    stringsAsFactors = FALSE
  )
}

validate_dated_tree <- function(path) {
  if (!file.exists(path) || file.info(path)$size <= 0) {
    stop("Missing/empty dated tree: ", path, call. = FALSE)
  }
  lines <- trimws(readLines(path, warn = FALSE))
  lines <- lines[nzchar(lines)]
  tree_text <- paste(lines, collapse = "")
  if (!length(lines) ||
      !endsWith(tree_text, ";") ||
      !grepl(":", tree_text, fixed = TRUE) ||
      grepl(
        "(^|[^[:alpha:]])(nan|inf)([^[:alpha:]]|$)",
        tree_text, ignore.case = TRUE, perl = TRUE
      )) {
    stop("Dated tree is not finite branch-length Newick: ", path, call. = FALSE)
  }
  dated <- tryCatch(
    ape::read.tree(path),
    error = function(e) {
      stop("Cannot parse dated tree ", path, ": ", conditionMessage(e), call. = FALSE)
    }
  )
  if (!inherits(dated, "phylo") ||
      ape::Ntip(dated) != ape::Ntip(reference_tree) ||
      anyDuplicated(dated$tip.label) ||
      !setequal(dated$tip.label, reference_tree$tip.label) ||
      nrow(dated$edge) != nrow(reference_tree$edge) ||
      is.null(dated$edge.length) ||
      length(dated$edge.length) != nrow(dated$edge) ||
      any(!is.finite(dated$edge.length)) ||
      any(dated$edge.length < 0)) {
    stop("Dated tree failed tip/edge/branch-length validation: ", path, call. = FALSE)
  }
  dated_topology <- dated
  dated_topology$edge.length <- NULL
  dated_topology$node.label <- NULL
  if (!isTRUE(ape::all.equal.phylo(
    reference_topology, dated_topology, use.edge.length = FALSE
  ))) {
    stop("Dated-tree topology differs from the input topology: ", path, call. = FALSE)
  }
  depths <- ape::node.depth.edgelength(dated)[seq_len(ape::Ntip(dated))]
  if (any(!is.finite(depths)) ||
      diff(range(depths)) > ultrametric_tolerance_ma) {
    stop(
      "Dated tree is not ultrametric within ",
      format(ultrametric_tolerance_ma, scientific = TRUE),
      " Ma: ", path, call. = FALSE
    )
  }
  invisible(TRUE)
}

validate_stage <- function(
    stage, replicate, expected_seed, expected_grid, expected_stop_text
) {
  score_path <- file.path(qa_dir, paste0(stage, "_scores.txt"))
  dated_path <- file.path(output_dir, paste0(stage, "_dated.tre"))
  config_path <- validate_config(
    stage, expected_seed, score_path, dated_path, expected_stop_text
  )
  read_metadata(stage, config_path)

  stderr_path <- file.path(log_dir, paste0(stage, ".stderr"))
  stdout_path <- file.path(log_dir, paste0(stage, ".stdout"))
  resource_path <- file.path(log_dir, paste0(stage, ".resources.txt"))
  for (path in c(stderr_path, stdout_path, resource_path, dated_path)) {
    if (!file.exists(path)) {
      stop("Missing stage artifact: ", path, call. = FALSE)
    }
  }
  if (file.info(stderr_path)$size != 0) {
    stop("Non-empty stderr: ", stderr_path, call. = FALSE)
  }
  if (file.info(stdout_path)$size <= 0 ||
      file.info(resource_path)$size <= 0 ||
      file.info(dated_path)$size <= 0) {
    stop("Empty required stage artifact for ", stage, call. = FALSE)
  }

  resource_lines <- readLines(resource_path, warn = FALSE)
  if (!any(grepl("Exit status: 0", resource_lines, fixed = TRUE)) ||
      !any(grepl(binary_expected, resource_lines, fixed = TRUE)) ||
      !any(grepl(config_path, resource_lines, fixed = TRUE))) {
    stop("Resource log does not bind a successful command for ", stage, call. = FALSE)
  }

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
    stop("Matched treePL failure string in ", stdout_path, call. = FALSE)
  }

  validate_dated_tree(dated_path)
  scores <- parse_scores(score_path, expected_grid, require_all_nonblank = TRUE)
  stdout_scores <- parse_scores(stdout_path, expected_grid)
  if (any(relative_difference(scores$raw_score, stdout_scores$raw_score) >
          stdout_score_tolerance)) {
    stop("cvoutfile/stdout score mismatch in ", stage, call. = FALSE)
  }

  scores$stage <- stage
  scores$replicate <- replicate
  scores$seed <- expected_seed
  scores
}

source_19point_decision <- validate_source_19point_decision()

formal <- do.call(
  rbind,
  lapply(seq_along(formal_ids), function(i) {
    validate_stage(
      formal_ids[[i]], i, expected_seeds[[i]], formal_grid,
      "0.00000000000000001"
    )
  })
)
rownames(formal) <- NULL
if (nrow(formal) != 69L) {
  stop("Formal 23-point score table must contain exactly 69 rows")
}
for (i in 1:3) {
  observed <- formal$smoothing[formal$replicate == i]
  if (!numeric_equal(observed, formal_grid)) {
    stop("Formal replicate lacks the exact ordered 23-point grid: ", i)
  }
}

formal$smoothing_label <- formal_labels[match(formal$smoothing, formal_grid)]
if (anyNA(formal$smoothing_label)) {
  stop("Cannot map one or more formal smoothing labels")
}
formal$replicate_rank <- ave(
  formal$raw_score, formal$replicate,
  FUN = function(x) rank(x, ties.method = "min")
)
formal$replicate_winner <- ave(
  formal$raw_score, formal$replicate,
  FUN = function(x) as.integer(relative_difference(x, min(x)) <= tie_tolerance)
) == 1L

summary_rows <- lapply(seq_along(formal_grid), function(i) {
  rows <- formal[formal$smoothing == formal_grid[[i]], , drop = FALSE]
  if (nrow(rows) != 3L) {
    stop("Formal smoothing value does not have exactly three replicates")
  }
  data.frame(
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
})
summary_tab <- do.call(rbind, summary_rows)
summary_tab$aggregate_rank <- rank(summary_tab$median_score, ties.method = "min")
aggregate_min <- min(summary_tab$median_score)
summary_tab$aggregate_winner <- relative_difference(
  summary_tab$median_score, aggregate_min
) <= tie_tolerance

aggregate_winners <- summary_tab$smoothing[summary_tab$aggregate_winner]
if (length(aggregate_winners) < 1L) {
  stop("No aggregate winner")
}
lower_boundary <- formal_grid[[length(formal_grid)]]
lower_boundary_score <- summary_tab$median_score[
  summary_tab$smoothing == lower_boundary
]
interior_scores <- summary_tab$median_score[
  summary_tab$smoothing %in% formal_grid[2:(length(formal_grid) - 1L)]
]
lower_boundary_in_aggregate_winner_set <- lower_boundary %in% aggregate_winners
lower_boundary_strictly_below_all_interior <-
  length(lower_boundary_score) == 1L &&
  length(interior_scores) == length(formal_grid) - 2L &&
  all(lower_boundary_score < interior_scores) &&
  all(relative_difference(lower_boundary_score, interior_scores) >
        tie_tolerance)
lowest_three_grid_values_all_aggregate_winners <- all(
  formal_grid[(length(formal_grid) - 2L):length(formal_grid)] %in%
    aggregate_winners
)
aggregate_winner_set_contains_interior <- any(
  aggregate_winners %in% formal_grid[2:(length(formal_grid) - 1L)]
)
finite_precision_plateau <-
  lower_boundary_in_aggregate_winner_set &&
  lowest_three_grid_values_all_aggregate_winners &&
  aggregate_winner_set_contains_interior
aggregate_boundary <- any(
  aggregate_winners %in% formal_grid[c(1L, length(formal_grid))]
)
boundary_requires_extension <- aggregate_boundary && !finite_precision_plateau
selected <- if (finite_precision_plateau) {
  max(aggregate_winners)
} else if (length(aggregate_winners) == 1L && !aggregate_boundary) {
  aggregate_winners[[1L]]
} else {
  NA_real_
}
selected_smoothing_rule <- if (finite_precision_plateau) {
  "largest_smoothing_in_aggregate_tied_minimum_plateau"
} else if (length(aggregate_winners) == 1L && !aggregate_boundary) {
  "unique_nonboundary_aggregate_minimum"
} else {
  "none_requires_further_cv_decision"
}

individual <- lapply(1:3, function(i) {
  rows <- formal[formal$replicate == i, , drop = FALSE]
  rows$smoothing[as.logical(rows$replicate_winner)]
})
individual_winner_counts <- lengths(individual)
individual_unique_winner <- individual_winner_counts == 1L
individual_flat <- unique(unlist(individual))
aggregate_support <- if (length(aggregate_winners) == 1L) {
  sum(vapply(individual, function(x) {
    length(x) == 1L &&
      relative_difference(x[[1L]], aggregate_winners[[1L]]) <= tie_tolerance
  }, logical(1)))
} else {
  0L
}
within_one_step <- length(aggregate_winners) == 1L &&
  all(abs(log10(individual_flat) - log10(aggregate_winners[[1L]])) <=
        1 + 1e-12)
stability <- if (length(aggregate_winners) == 1L && aggregate_support >= 2L) {
  "strong"
} else if (within_one_step) {
  "moderate"
} else {
  "weak"
}
discordant <- identical(stability, "weak")

sensitivity_values <- unique(c(aggregate_winners, individual_flat))
if (finite_precision_plateau) {
  sensitivity_values <- sensitivity_values[!sensitivity_values %in% selected]
}
sensitivity_values <- formal_grid[formal_grid %in% sensitivity_values]
label_for <- function(x) formal_labels[match(x, formal_grid)]
individual_text <- paste(vapply(
  1:3,
  function(i) paste(label_for(individual[[i]]), collapse = ","),
  character(1)
), collapse = ";")

formal$aggregate_rank <- summary_tab$aggregate_rank[
  match(formal$smoothing, summary_tab$smoothing)
]
formal$aggregate_winner <- summary_tab$aggregate_winner[
  match(formal$smoothing, summary_tab$smoothing)
]
formal <- formal[, c(
  "stage", "replicate", "seed", "smoothing", "smoothing_label", "raw_score",
  "replicate_rank", "replicate_winner", "aggregate_rank", "aggregate_winner"
)]

diagnostic_differences <- unlist(lapply(1:3, function(i) {
  diagnostic <- validate_stage(
    diagnostic_ids[[i]], i, expected_seeds[[i]], diagnostic_grid,
    "0.0000000000001"
  )
  formal_prefix <- formal[
    formal$replicate == i & formal$smoothing %in% diagnostic_grid,
    c("smoothing", "raw_score")
  ]
  if (nrow(formal_prefix) != 19L ||
      !numeric_equal(diagnostic$smoothing, formal_prefix$smoothing)) {
    stop("Paired 19-point prefix mismatch for replicate ", i)
  }
  relative_difference(diagnostic$raw_score, formal_prefix$raw_score)
}))
diagnostic_all_prefix_points_present <- length(diagnostic_differences) == 57L
if (!diagnostic_all_prefix_points_present ||
    any(!is.finite(diagnostic_differences))) {
  stop("Incomplete/non-finite paired 19-point prefix diagnostic")
}
diagnostic_max_relative_difference <- max(diagnostic_differences)

decision_status <- if (boundary_requires_extension) {
  "REQUIRES_GRID_EXTENSION"
} else if (finite_precision_plateau) {
  "ACCEPTED_WITH_SENSITIVITIES"
} else if (length(aggregate_winners) > 1L) {
  "AGGREGATE_TIE_REQUIRES_SENSITIVITIES"
} else if (discordant) {
  "ACCEPTED_WITH_SENSITIVITIES"
} else {
  "ACCEPTED"
}

decision <- data.frame(
  status = decision_status,
  decision_basis =
    "lowest_median_raw_chisq_across_three_23_point_randomCV_runs",
  selected_smoothing_rule = selected_smoothing_rule,
  selected_smoothing =
    if (is.na(selected)) NA_character_ else label_for(selected),
  selected_median_score = if (is.na(selected)) NA_real_ else
    summary_tab$median_score[summary_tab$smoothing == selected],
  aggregate_winner_at_grid_boundary = aggregate_boundary,
  boundary_requires_extension = boundary_requires_extension,
  aggregate_winner_count = length(aggregate_winners),
  aggregate_winner_values =
    paste(label_for(aggregate_winners), collapse = ","),
  lower_boundary_in_aggregate_winner_set =
    lower_boundary_in_aggregate_winner_set,
  lower_boundary_strictly_below_all_interior =
    lower_boundary_strictly_below_all_interior,
  lowest_three_grid_values_all_aggregate_winners =
    lowest_three_grid_values_all_aggregate_winners,
  aggregate_winner_set_contains_interior =
    aggregate_winner_set_contains_interior,
  finite_precision_plateau = finite_precision_plateau,
  stability = stability,
  aggregate_replicate_winner_support = aggregate_support,
  aggregate_unique_replicate_winner_support = aggregate_support,
  replicate_winner_count_run01 = individual_winner_counts[[1L]],
  replicate_winner_count_run02 = individual_winner_counts[[2L]],
  replicate_winner_count_run03 = individual_winner_counts[[3L]],
  replicate_has_unique_winner_run01 = individual_unique_winner[[1L]],
  replicate_has_unique_winner_run02 = individual_unique_winner[[2L]],
  replicate_has_unique_winner_run03 = individual_unique_winner[[3L]],
  replicate_winners_run01_run02_run03 = individual_text,
  sensitivity_required = finite_precision_plateau || discordant ||
    length(aggregate_winners) > 1L,
  sensitivity_smoothing_values =
    paste(label_for(sensitivity_values), collapse = ","),
  formal_runs = paste(formal_ids, collapse = ","),
  formal_grid_points_per_run = length(formal_grid),
  formal_rows = nrow(formal),
  formal_seed_values = paste(expected_seeds, collapse = ","),
  formal_exit_codes_all_zero = TRUE,
  formal_stderr_all_empty = TRUE,
  formal_dated_trees_all_nonempty = TRUE,
  formal_dated_trees_valid_branch_length_newick = TRUE,
  formal_dated_tree_tip_sets_match_input = TRUE,
  formal_dated_tree_topologies_match_input = TRUE,
  formal_dated_trees_ultrametric_within_tolerance = TRUE,
  formal_dated_trees_ultrametric_tolerance_ma = ultrametric_tolerance_ma,
  formal_config_contracts_validated = TRUE,
  formal_config_lints_all_pass = TRUE,
  formal_metadata_validated = TRUE,
  formal_stdout_score_sets_validated = TRUE,
  formal_resource_logs_validated = TRUE,
  source_19point_decision_validated = TRUE,
  diagnostic_19point_lower_boundary_in_aggregate_winner_set =
    source_19point_decision$lower_boundary_in_aggregate_winner_set,
  diagnostic_19point_aggregate_winner_values =
    paste(source_19point_decision$aggregate_winner_values, collapse = ","),
  diagnostic_19point_full_stage_validation_passed = TRUE,
  diagnostic_19point_runs = paste(diagnostic_ids, collapse = ","),
  diagnostic_19point_grid_points_per_run = length(diagnostic_grid),
  diagnostic_19point_rows = 3L * length(diagnostic_grid),
  diagnostic_19point_seed_values = paste(expected_seeds, collapse = ","),
  diagnostic_19point_runs_used_for_selection = FALSE,
  diagnostic_19point_prefix_comparisons = length(diagnostic_differences),
  diagnostic_19point_all_prefix_points_present =
    diagnostic_all_prefix_points_present,
  diagnostic_19point_max_relative_difference =
    diagnostic_max_relative_difference,
  tie_relative_tolerance = tie_tolerance,
  stringsAsFactors = FALSE
)

write.table(
  formal, out_long, sep = "	", quote = FALSE, row.names = FALSE, na = "NA"
)
write.table(
  summary_tab, out_summary, sep = "	", quote = FALSE,
  row.names = FALSE, na = "NA"
)
write.table(
  decision, out_decision, sep = "	", quote = FALSE,
  row.names = FALSE, na = "NA"
)

cat("Formal 23-point random-CV rows:", nrow(formal), "\n")
cat("Aggregate winner:", paste(label_for(aggregate_winners), collapse = ","), "\n")
cat("Boundary:", aggregate_boundary, "\n")
cat("Finite-precision plateau:", finite_precision_plateau, "\n")
cat("Stability:", stability, "\n")
cat("Status:", decision_status, "\n")
cat("Sensitivity values:", paste(label_for(sensitivity_values), collapse = ","), "\n")
cat(
  "19-point paired-prefix max relative difference:",
  format(diagnostic_max_relative_difference, digits = 12), "\n"
)

if (decision_status %in% c(
  "REQUIRES_GRID_EXTENSION", "AGGREGATE_TIE_REQUIRES_SENSITIVITIES"
)) {
  quit(status = 2L)
}
