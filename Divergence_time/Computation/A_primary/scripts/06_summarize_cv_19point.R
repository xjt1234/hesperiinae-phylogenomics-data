#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: 06_summarize_cv_19point.R <run_root>", call. = FALSE)
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

out_long <- file.path(qa_dir, "cv_scores_long.tsv")
out_summary <- file.path(qa_dir, "cv_summary.tsv")
out_decision <- file.path(qa_dir, "cv_decision.tsv")
for (p in c(out_long, out_summary, out_decision)) {
  if (file.exists(p)) stop("Refusing to overwrite existing output: ", p, call. = FALSE)
}

declarations <- file.path(
  qa_dir,
  c(
    "cv_selection_rule.md",
    "cv_grid_extension_decision.md",
    "cv_continuous_grid_decision.md",
    "cv_15point_boundary_decision.md"
  )
)
if (any(!file.exists(declarations))) {
  stop(
    "Missing prespecified decision document(s): ",
    paste(basename(declarations[!file.exists(declarations)]), collapse = ", "),
    call. = FALSE
  )
}

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
formal_ids <- sprintf("cv_finalgrid_run_%02d", 1:3)
diagnostic_ids <- sprintf("cv_full_run_%02d", 1:3)
expected_seeds <- as.character(2026090410 + 1:3)
tie_tolerance <- 1e-10
stdout_score_tolerance <- 5e-5

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

get_one_value <- function(active, key, path) {
  pattern <- paste0("^", key, "[[:space:]]*=[[:space:]]*(.+?)[[:space:]]*$")
  hits <- grep(pattern, active, perl = TRUE)
  if (length(hits) != 1L) {
    stop("Expected exactly one ", key, " assignment in ", path, call. = FALSE)
  }
  sub(pattern, "\\1", active[hits], perl = TRUE)
}

validate_config <- function(stage, expected_seed, expected_score, expected_tree) {
  path <- file.path(config_dir, paste0(stage, ".cfg"))
  if (!file.exists(path)) stop("Missing config: ", path, call. = FALSE)
  active <- read_active_config(path)

  required_exact <- c(
    "randomcv",
    "mrca = CROWN_PAPILIONOIDEA Parnassius_apollo_ncbi Aeromachus_catocyanea_mydata",
    "min = CROWN_PAPILIONOIDEA 91.5046",
    "max = CROWN_PAPILIONOIDEA 100.8925",
    "mrca = CROWN_PAPILIONIDAE_NONBARONIA Parnassius_apollo_ncbi Graphium_cloanthus_mydata",
    "min = CROWN_PAPILIONIDAE_NONBARONIA 44.1968",
    "max = CROWN_PAPILIONIDAE_NONBARONIA 52.9473",
    "opt = 2",
    "optad = 2",
    "moredetailad",
    "optcvad = 1",
    "moredetailcvad"
  )
  bad_counts <- required_exact[vapply(
    required_exact, function(x) sum(active == x) != 1L, logical(1)
  )]
  if (length(bad_counts)) {
    stop(
      "Required exact config line missing/non-unique in ", path, ": ",
      paste(bad_counts, collapse = "; "), call. = FALSE
    )
  }
  if (any(active %in% c("prime", "cv", "thorough", "moredetail"))) {
    stop("Incompatible mode/detail directive in ", path, call. = FALSE)
  }
  calibration_lines <- active[
    grepl("^(mrca|min|max)[[:space:]]*=", active)
  ]
  if (any(grepl("HESPERIINAE|101\\.4|36\\.2|40\\.7",
                calibration_lines, ignore.case = TRUE)) ||
      any(grepl("^log_pen([[:space:]]*=|$)", active, ignore.case = TRUE))) {
    stop("Forbidden calibration/penalty content in ", path, call. = FALSE)
  }

  expected_text <- list(
    treefile = tree_expected,
    numsites = "182682",
    outfile = expected_tree,
    cvoutfile = expected_score,
    cviter = "3",
    nthreads = "10",
    seed = expected_seed
  )
  for (key in names(expected_text)) {
    observed <- get_one_value(active, key, path)
    expected <- expected_text[[key]]
    if (key %in% c("treefile", "outfile", "cvoutfile")) {
      observed <- canonical(observed)
      expected <- canonical(expected)
    }
    if (!identical(observed, expected)) {
      stop(
        key, " mismatch in ", path, ": observed=", observed,
        " expected=", expected, call. = FALSE
      )
    }
  }

  expected_numeric <- c(cvstart = 1e5, cvstop = 1e-13, cvmultstep = 0.1)
  for (key in names(expected_numeric)) {
    observed_text <- get_one_value(active, key, path)
    observed <- suppressWarnings(as.numeric(observed_text))
    if (!numeric_equal(observed, expected_numeric[[key]])) {
      stop(key, " mismatch in ", path, ": ", observed_text, call. = FALSE)
    }
  }

  lint_path <- file.path(qa_dir, paste0("config_lint_", stage, ".tsv"))
  if (!file.exists(lint_path)) stop("Missing config lint: ", lint_path, call. = FALSE)
  lint <- read.delim(lint_path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("check", "pass") %in% names(lint)) || !all(lint$pass)) {
    stop("Config lint is absent or contains a failure: ", lint_path, call. = FALSE)
  }
  invisible(path)
}

read_metadata <- function(stage, config_path) {
  path <- file.path(log_dir, paste0(stage, ".metadata.tsv"))
  if (!file.exists(path)) stop("Missing metadata: ", path, call. = FALSE)
  tab <- read.delim(
    path, header = FALSE, sep = "\t", quote = "", comment.char = "",
    col.names = c("key", "value"), stringsAsFactors = FALSE
  )
  one <- function(key) {
    value <- tab$value[tab$key == key]
    if (length(value) != 1L) stop("Metadata key not unique: ", key, " in ", path)
    value
  }
  if (!identical(one("stage"), stage)) stop("Metadata stage mismatch: ", path)
  if (!identical(one("exit_code"), "0")) stop("Nonzero treePL exit in ", path)
  if (!identical(canonical(one("binary")), canonical(binary_expected))) {
    stop("Metadata binary mismatch: ", path)
  }
  if (!identical(canonical(one("config")), canonical(config_path))) {
    stop("Metadata config mismatch: ", path)
  }
  invisible(tab)
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
    stop("Malformed, non-finite, or out-of-order score grid in ", path, call. = FALSE)
  }
  data.frame(
    smoothing = smoothing,
    smoothing_text = smoothing_text,
    raw_score = score,
    stringsAsFactors = FALSE
  )
}

validate_stage <- function(stage, replicate, expected_seed) {
  score_path <- file.path(qa_dir, paste0(stage, "_scores.txt"))
  dated_path <- file.path(output_dir, paste0(stage, "_dated.tre"))
  config_path <- validate_config(stage, expected_seed, score_path, dated_path)
  read_metadata(stage, config_path)

  stderr_path <- file.path(log_dir, paste0(stage, ".stderr"))
  stdout_path <- file.path(log_dir, paste0(stage, ".stdout"))
  resource_path <- file.path(log_dir, paste0(stage, ".resources.txt"))
  for (path in c(stderr_path, stdout_path, resource_path, dated_path)) {
    if (!file.exists(path)) stop("Missing stage artifact: ", path, call. = FALSE)
  }
  if (file.info(stderr_path)$size != 0) stop("Non-empty stderr: ", stderr_path)
  if (file.info(stdout_path)$size <= 0 || file.info(resource_path)$size <= 0 ||
      file.info(dated_path)$size <= 0) {
    stop("Empty required stage artifact for ", stage, call. = FALSE)
  }
  stdout_lines <- readLines(stdout_path, warn = FALSE)
  failure_re <- paste(
    c(
      "problem", "complete failure", "segmentation", "(^|[^[:alpha:]])nan([^[:alpha:]]|$)",
      "(^|[^[:alpha:]])inf([^[:alpha:]]|$)",
      "(^|[^[:alpha:]])error([^[:alpha:]]|$)"
    ),
    collapse = "|"
  )
  if (any(grepl(failure_re, stdout_lines, ignore.case = TRUE, perl = TRUE))) {
    stop("Matched treePL failure string in ", stdout_path, call. = FALSE)
  }

  scores <- parse_scores(score_path, formal_grid, require_all_nonblank = TRUE)
  stdout_scores <- parse_scores(stdout_path, formal_grid)
  if (any(relative_difference(scores$raw_score, stdout_scores$raw_score) >
          stdout_score_tolerance)) {
    stop("cvoutfile/stdout score mismatch in ", stage, call. = FALSE)
  }

  scores$stage <- stage
  scores$replicate <- replicate
  scores$seed <- expected_seed
  scores
}

formal <- do.call(
  rbind,
  lapply(seq_along(formal_ids), function(i) {
    validate_stage(formal_ids[[i]], i, expected_seeds[[i]])
  })
)
rownames(formal) <- NULL

if (nrow(formal) != 57L) stop("Formal score table must contain 57 rows")
for (i in 1:3) {
  observed <- formal$smoothing[formal$replicate == i]
  if (!numeric_equal(observed, formal_grid)) {
    stop("Formal replicate does not contain exact ordered 19-point grid: ", i)
  }
}

formal$smoothing_label <- formal_labels[
  match(formal$smoothing, formal_grid)
]
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
if (length(aggregate_winners) < 1L) stop("No aggregate winner")
aggregate_boundary <- any(aggregate_winners %in% formal_grid[c(1L, length(formal_grid))])
selected <- if (length(aggregate_winners) == 1L && !aggregate_boundary) {
  aggregate_winners[[1]]
} else {
  NA_real_
}

individual <- lapply(1:3, function(i) {
  rows <- formal[formal$replicate == i, , drop = FALSE]
  rows$smoothing[as.logical(rows$replicate_winner)]
})
individual_flat <- unique(unlist(individual))
aggregate_support <- if (length(aggregate_winners) == 1L) {
  sum(vapply(individual, function(x) aggregate_winners[[1]] %in% x, logical(1)))
} else {
  0L
}
within_one_step <- length(aggregate_winners) == 1L &&
  all(abs(log10(individual_flat) - log10(aggregate_winners[[1]])) <= 1 + 1e-12)
stability <- if (length(aggregate_winners) == 1L && aggregate_support >= 2L) {
  "strong"
} else if (within_one_step) {
  "moderate"
} else {
  "weak"
}
discordant <- identical(stability, "weak")

sensitivity_values <- unique(c(aggregate_winners, individual_flat))
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

diagnostic_values <- list()
for (i in 1:3) {
  path <- file.path(qa_dir, paste0(diagnostic_ids[[i]], "_scores.txt"))
  tab <- parse_scores(path, formal_grid[1:15], require_all_nonblank = TRUE)
  formal_prefix <- formal[
    formal$replicate == i & formal$smoothing %in% formal_grid[1:15],
    c("smoothing", "raw_score")
  ]
  diagnostic_values[[i]] <- relative_difference(
    tab$raw_score, formal_prefix$raw_score
  )
}
diagnostic_differences <- unlist(diagnostic_values)
diagnostic_max_relative_difference <- max(diagnostic_differences)
diagnostic_all_prefix_points_present <- length(diagnostic_differences) == 45L

decision_status <- if (aggregate_boundary) {
  "REQUIRES_GRID_EXTENSION"
} else if (length(aggregate_winners) > 1L) {
  "AGGREGATE_TIE_REQUIRES_SENSITIVITIES"
} else if (discordant) {
  "ACCEPTED_WITH_SENSITIVITIES"
} else {
  "ACCEPTED"
}

decision <- data.frame(
  status = decision_status,
  decision_basis = "lowest_median_raw_chisq_across_three_19_point_randomCV_runs",
  selected_smoothing = if (is.na(selected)) NA_character_ else label_for(selected),
  selected_median_score = if (is.na(selected)) NA_real_ else
    summary_tab$median_score[summary_tab$smoothing == selected],
  aggregate_winner_at_grid_boundary = aggregate_boundary,
  aggregate_winner_count = length(aggregate_winners),
  stability = stability,
  aggregate_replicate_winner_support = aggregate_support,
  replicate_winners_run01_run02_run03 = individual_text,
  sensitivity_required = discordant || length(aggregate_winners) > 1L,
  sensitivity_smoothing_values = paste(label_for(sensitivity_values), collapse = ","),
  formal_runs = paste(formal_ids, collapse = ","),
  formal_grid_points_per_run = length(formal_grid),
  formal_rows = nrow(formal),
  formal_seed_values = paste(expected_seeds, collapse = ","),
  formal_exit_codes_all_zero = TRUE,
  formal_stderr_all_empty = TRUE,
  formal_dated_trees_all_nonempty = TRUE,
  diagnostic_15point_runs_used_for_selection = FALSE,
  diagnostic_15point_prefix_comparisons = length(diagnostic_differences),
  diagnostic_15point_all_prefix_points_present =
    diagnostic_all_prefix_points_present,
  diagnostic_15point_max_relative_difference =
    diagnostic_max_relative_difference,
  tie_relative_tolerance = tie_tolerance,
  stringsAsFactors = FALSE
)

write.table(
  formal, out_long, sep = "\t", quote = FALSE, row.names = FALSE,
  na = "NA"
)
write.table(
  summary_tab, out_summary, sep = "\t", quote = FALSE, row.names = FALSE,
  na = "NA"
)
write.table(
  decision, out_decision, sep = "\t", quote = FALSE, row.names = FALSE,
  na = "NA"
)

cat("Formal random-CV rows:", nrow(formal), "\n")
cat("Aggregate winner:", paste(label_for(aggregate_winners), collapse = ","), "\n")
cat("Boundary:", aggregate_boundary, "\n")
cat("Stability:", stability, "\n")
cat("Status:", decision_status, "\n")
cat("Sensitivity values:", paste(label_for(sensitivity_values), collapse = ","), "\n")

if (aggregate_boundary || length(aggregate_winners) > 1L) quit(status = 2L)
