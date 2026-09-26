#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, digits = 15)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: 13b_generate_final_configs_from_23point.R <run_root>", call. = FALSE)
}

run_root <- normalizePath(args[[1]], mustWork = TRUE)
config_dir <- file.path(run_root, "configs")
input_dir <- file.path(run_root, "input")
log_dir <- file.path(run_root, "logs")
output_dir <- file.path(run_root, "output")
qa_dir <- file.path(run_root, "qa")
script_dir <- file.path(run_root, "scripts")
required_dirs <- c(config_dir, input_dir, log_dir, output_dir, qa_dir, script_dir)
if (any(!dir.exists(required_dirs))) {
  stop(
    "Missing run directory/directories: ",
    paste(required_dirs[!dir.exists(required_dirs)], collapse = ", "),
    call. = FALSE
  )
}
if (grepl("[[:space:][:cntrl:]]", run_root)) {
  stop("run_root may not contain whitespace/control characters", call. = FALSE)
}

decision_path <- file.path(qa_dir, "cv_23point_decision.tsv")
prime_recommendations_path <- file.path(qa_dir, "prime_recommendations.cfg")
plateau_rule_path <- file.path(
  qa_dir, "cv_tie_plateau_rule_declared_before_low_extension_scores.md"
)
treefile <- file.path(input_dir, "verified_495_tip_input.treefile")
binary <- file.path(run_root, "software", "treePL_randomcv_fix_src", "treePL")
lint_script <- file.path(script_dir, "04_lint_treepl_config.R")
wrapper <- file.path(script_dir, "run_treepl_stage.sh")
required_inputs <- c(
  decision_path, prime_recommendations_path, plateau_rule_path,
  treefile, binary, lint_script, wrapper
)
if (any(!file.exists(required_inputs))) {
  stop(
    "Missing required input(s): ",
    paste(required_inputs[!file.exists(required_inputs)], collapse = ", "),
    call. = FALSE
  )
}
if (isTRUE(file.info(treefile)$isdir) || file.info(treefile)$size <= 0) {
  stop("Fixed input tree is not a non-empty regular file: ", treefile, call. = FALSE)
}
if (isTRUE(file.info(plateau_rule_path)$isdir) ||
    file.info(plateau_rule_path)$size <= 0) {
  stop(
    "Finite-precision plateau rule is empty/invalid: ", plateau_rule_path,
    call. = FALSE
  )
}

read_single_row_tsv <- function(path) {
  tab <- read.delim(
    path,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    colClasses = "character",
    check.names = FALSE,
    na.strings = character(),
    strip.white = FALSE
  )
  if (nrow(tab) != 1L) {
    stop("Expected exactly one decision row in ", path, call. = FALSE)
  }
  if (anyDuplicated(names(tab))) {
    stop("Duplicate column name(s) in ", path, call. = FALSE)
  }
  tab
}

decision <- read_single_row_tsv(decision_path)
required_decision_columns <- c(
  "status",
  "decision_basis",
  "selected_smoothing",
  "selected_median_score",
  "aggregate_winner_at_grid_boundary",
  "boundary_requires_extension",
  "aggregate_winner_count",
  "aggregate_winner_values",
  "lower_boundary_in_aggregate_winner_set",
  "lower_boundary_strictly_below_all_interior",
  "lowest_three_grid_values_all_aggregate_winners",
  "aggregate_winner_set_contains_interior",
  "finite_precision_plateau",
  "selected_smoothing_rule",
  "stability",
  "replicate_winners_run01_run02_run03",
  "sensitivity_required",
  "sensitivity_smoothing_values",
  "formal_runs",
  "formal_grid_points_per_run",
  "formal_rows",
  "formal_seed_values",
  "formal_exit_codes_all_zero",
  "formal_stderr_all_empty",
  "formal_dated_trees_all_nonempty",
  "formal_dated_trees_valid_branch_length_newick",
  "diagnostic_19point_full_stage_validation_passed",
  "diagnostic_19point_runs_used_for_selection",
  "diagnostic_19point_all_prefix_points_present"
)
missing_decision_columns <- setdiff(required_decision_columns, names(decision))
if (length(missing_decision_columns)) {
  stop(
    "CV decision is missing required column(s): ",
    paste(missing_decision_columns, collapse = ", "),
    call. = FALSE
  )
}

decision_value <- function(key) {
  value <- decision[[key]][[1L]]
  if (is.na(value)) value <- ""
  if (!identical(value, trimws(value))) {
    stop("CV decision field has surrounding whitespace: ", key, call. = FALSE)
  }
  value
}

exact_bool <- function(key) {
  value <- decision_value(key)
  if (!value %in% c("TRUE", "FALSE")) {
    stop("CV decision field must be exact TRUE/FALSE: ", key, call. = FALSE)
  }
  identical(value, "TRUE")
}

exact_uint <- function(key) {
  value <- decision_value(key)
  if (!grepl("^(0|[1-9][0-9]*)$", value)) {
    stop("CV decision field must be an unsigned integer: ", key, call. = FALSE)
  }
  parsed <- suppressWarnings(as.numeric(value))
  if (!is.finite(parsed) || parsed > 2^53) {
    stop("CV decision integer is outside the exact numeric range: ", key, call. = FALSE)
  }
  parsed
}

parse_positive_number <- function(value, field) {
  numeric_pattern <- paste0(
    "^(?:[0-9]+(?:\\.[0-9]*)?|\\.[0-9]+)",
    "(?:[eE][+-]?[0-9]+)?$"
  )
  if (!nzchar(value) || !grepl(numeric_pattern, value, perl = TRUE)) {
    stop(field, " must be one positive numeric token: ", value, call. = FALSE)
  }
  parsed <- suppressWarnings(as.numeric(value))
  if (length(parsed) != 1L || !is.finite(parsed) || parsed <= 0) {
    stop(field, " must be finite and positive: ", value, call. = FALSE)
  }
  parsed
}

same_number <- function(a, b, tolerance = 1e-12) {
  abs(a - b) <= tolerance * max(abs(a), abs(b), .Machine$double.xmin)
}

formal_grid <- 10^(5:-17)

parse_number_tokens <- function(text, field) {
  if (!nzchar(text)) stop(field, " may not be empty", call. = FALSE)
  raw_tokens <- strsplit(text, ",", fixed = TRUE)[[1L]]
  tokens <- trimws(raw_tokens)
  if (any(!nzchar(tokens)) || !identical(tokens, raw_tokens)) {
    stop(field, " contains an empty token or whitespace", call. = FALSE)
  }
  values <- vapply(
    seq_along(tokens),
    function(i) parse_positive_number(
      tokens[[i]], paste0(field, "[", i, "]")
    ),
    numeric(1)
  )
  list(text = tokens, numeric = values)
}

formal_grid_indices <- function(values, field) {
  vapply(seq_along(values), function(i) {
    hits <- which(vapply(
      formal_grid, same_number, logical(1), b = values[[i]]
    ))
    if (length(hits) != 1L) {
      stop(
        field, " contains a value outside the formal 23-point grid: ",
        values[[i]], call. = FALSE
      )
    }
    hits[[1L]]
  }, integer(1))
}

accepted_statuses <- c("ACCEPTED", "ACCEPTED_WITH_SENSITIVITIES")
status <- decision_value("status")
if (!status %in% accepted_statuses) {
  stop("CV decision is not qualified for final runs: status=", status, call. = FALSE)
}
if (!identical(
  decision_value("decision_basis"),
  "lowest_median_raw_chisq_across_three_23_point_randomCV_runs"
)) {
  stop("Unexpected CV decision basis", call. = FALSE)
}
if (exact_bool("boundary_requires_extension")) {
  stop(
    "Refusing final configs: unresolved boundary requires another CV extension",
    call. = FALSE
  )
}
aggregate_winner_count <- exact_uint("aggregate_winner_count")
aggregate_boundary <- exact_bool("aggregate_winner_at_grid_boundary")
finite_precision_plateau <- exact_bool("finite_precision_plateau")
if (aggregate_winner_count < 1) {
  stop("Refusing final configs: aggregate CV winner set is empty", call. = FALSE)
}
if (finite_precision_plateau) {
  if (!identical(status, "ACCEPTED_WITH_SENSITIVITIES") ||
      !aggregate_boundary || aggregate_winner_count < 3) {
    stop("Malformed finite-precision plateau decision", call. = FALSE)
  }
} else {
  if (aggregate_boundary) {
    stop(
      "Refusing final configs: ordinary boundary winner requires more CV",
      call. = FALSE
    )
  }
  if (!identical(aggregate_winner_count, 1)) {
    stop(
      "Refusing final configs: non-plateau aggregate CV winner is not unique",
      call. = FALSE
    )
  }
}

qualification_true_fields <- c(
  "formal_exit_codes_all_zero",
  "formal_stderr_all_empty",
  "formal_dated_trees_all_nonempty",
  "formal_dated_trees_valid_branch_length_newick",
  "diagnostic_19point_full_stage_validation_passed",
  "diagnostic_19point_all_prefix_points_present"
)
failed_qualification <- qualification_true_fields[!vapply(
  qualification_true_fields, exact_bool, logical(1)
)]
if (length(failed_qualification)) {
  stop(
    "CV decision failed qualification field(s): ",
    paste(failed_qualification, collapse = ", "),
    call. = FALSE
  )
}
if (exact_bool("diagnostic_19point_runs_used_for_selection")) {
  stop("Diagnostic 19-point runs may not have been used for selection", call. = FALSE)
}
if (!identical(exact_uint("formal_grid_points_per_run"), 23) ||
    !identical(exact_uint("formal_rows"), 69)) {
  stop("Formal CV decision does not describe three complete 23-point runs", call. = FALSE)
}
if (!identical(
  decision_value("formal_runs"),
  "cv_extendedgrid_run_01,cv_extendedgrid_run_02,cv_extendedgrid_run_03"
) || !identical(
  decision_value("formal_seed_values"),
  "2026090441,2026090442,2026090443"
)) {
  stop("Formal CV stage/seed provenance is not the prespecified design", call. = FALSE)
}

selected_smooth_text <- decision_value("selected_smoothing")
selected_smooth <- parse_positive_number(
  selected_smooth_text, "selected_smoothing"
)
selected_grid_index <- formal_grid_indices(
  selected_smooth, "selected_smoothing"
)
aggregate_winners <- parse_number_tokens(
  decision_value("aggregate_winner_values"), "aggregate_winner_values"
)
aggregate_grid_indices <- formal_grid_indices(
  aggregate_winners$numeric, "aggregate_winner_values"
)
if (length(aggregate_grid_indices) != aggregate_winner_count ||
    anyDuplicated(aggregate_grid_indices) ||
    (length(aggregate_grid_indices) > 1L &&
     any(diff(aggregate_grid_indices) <= 0L))) {
  stop(
    "aggregate_winner_values must contain exactly aggregate_winner_count ",
    "unique values in descending-smoothing grid order",
    call. = FALSE
  )
}
if (!selected_grid_index %in% aggregate_grid_indices) {
  stop("selected_smoothing is absent from aggregate_winner_values", call. = FALSE)
}

actual_boundary <- any(
  aggregate_grid_indices %in% c(1L, length(formal_grid))
)
actual_lower_boundary <- length(formal_grid) %in% aggregate_grid_indices
actual_lowest_three <- all(
  tail(seq_along(formal_grid), 3L) %in% aggregate_grid_indices
)
actual_contains_interior <- any(
  aggregate_grid_indices %in% seq.int(2L, length(formal_grid) - 1L)
)
if (!identical(aggregate_boundary, actual_boundary) ||
    !identical(
      exact_bool("lower_boundary_in_aggregate_winner_set"),
      actual_lower_boundary
    ) ||
    !identical(
      exact_bool("lowest_three_grid_values_all_aggregate_winners"),
      actual_lowest_three
    ) ||
    !identical(
      exact_bool("aggregate_winner_set_contains_interior"),
      actual_contains_interior
    )) {
  stop("Aggregate-winner flags disagree with aggregate_winner_values", call. = FALSE)
}

if (finite_precision_plateau) {
  if (!actual_lower_boundary || !actual_lowest_three ||
      !actual_contains_interior ||
      exact_bool("lower_boundary_strictly_below_all_interior") ||
      !identical(
        decision_value("selected_smoothing_rule"),
        "largest_smoothing_in_aggregate_tied_minimum_plateau"
      ) ||
      selected_grid_index != min(aggregate_grid_indices)) {
    stop(
      "Finite-precision plateau does not satisfy the predeclared conservative ",
      "tie-break rule",
      call. = FALSE
    )
  }
} else {
  if (exact_bool("lower_boundary_strictly_below_all_interior") ||
      !identical(
        decision_value("selected_smoothing_rule"),
        "unique_nonboundary_aggregate_minimum"
      ) ||
      !identical(selected_grid_index, aggregate_grid_indices[[1L]])) {
    stop("Malformed unique non-boundary aggregate decision", call. = FALSE)
  }
}
selected_median_score <- suppressWarnings(as.numeric(
  decision_value("selected_median_score")
))
if (!is.finite(selected_median_score) || selected_median_score < 0) {
  stop("selected_median_score must be finite and nonnegative", call. = FALSE)
}
if (!decision_value("stability") %in% c("strong", "moderate", "weak")) {
  stop("Unknown CV stability classification", call. = FALSE)
}

sensitivity_required <- exact_bool("sensitivity_required")
if (!identical(
  status,
  if (sensitivity_required) "ACCEPTED_WITH_SENSITIVITIES" else "ACCEPTED"
)) {
  stop("CV status and sensitivity_required are inconsistent", call. = FALSE)
}
if (identical(status, "ACCEPTED_WITH_SENSITIVITIES") &&
    !identical(decision_value("stability"), "weak")) {
  stop("Sensitivity-required decision must have stability=weak", call. = FALSE)
}


sensitivity_text <- character()
sensitivity_numeric <- numeric()
if (sensitivity_required) {
  raw_sensitivity <- decision_value("sensitivity_smoothing_values")
  if (!nzchar(raw_sensitivity)) {
    stop("Sensitivity runs are required but no smoothing values were supplied", call. = FALSE)
  }
  candidates <- parse_number_tokens(
    raw_sensitivity, "sensitivity_smoothing_values"
  )
  candidate_grid_indices <- formal_grid_indices(
    candidates$numeric, "sensitivity_smoothing_values"
  )
  if (anyDuplicated(candidate_grid_indices) ||
      (length(candidate_grid_indices) > 1L &&
       any(diff(candidate_grid_indices) <= 0L))) {
    stop(
      "sensitivity_smoothing_values must be deduplicated and in ",
      "descending-smoothing grid order",
      call. = FALSE
    )
  }

  replicate_groups <- strsplit(
    decision_value("replicate_winners_run01_run02_run03"),
    ";", fixed = TRUE
  )[[1L]]
  if (length(replicate_groups) != 3L || any(!nzchar(replicate_groups))) {
    stop("Expected three non-empty replicate-winner groups", call. = FALSE)
  }
  replicate_winners <- lapply(seq_along(replicate_groups), function(i) {
    parse_number_tokens(
      replicate_groups[[i]],
      paste0("replicate_winners_run", sprintf("%02d", i))
    )
  })
  replicate_grid_indices <- unlist(lapply(
    seq_along(replicate_winners),
    function(i) formal_grid_indices(
      replicate_winners[[i]]$numeric,
      paste0("replicate_winners_run", sprintf("%02d", i))
    )
  ), use.names = FALSE)

  expected_sensitivity_indices <- sort(unique(c(
    aggregate_grid_indices[aggregate_grid_indices != selected_grid_index],
    replicate_grid_indices[replicate_grid_indices != selected_grid_index]
  )))
  observed_sensitivity_indices <- candidate_grid_indices[
    candidate_grid_indices != selected_grid_index
  ]
  if (finite_precision_plateau &&
      any(candidate_grid_indices == selected_grid_index)) {
    stop("Plateau sensitivity list must exclude the primary smoothing", call. = FALSE)
  }
  if (!identical(
    observed_sensitivity_indices,
    expected_sensitivity_indices
  )) {
    stop(
      "Sensitivity list must contain every other aggregate tie and every ",
      "distinct replicate winner, after deduplication and primary exclusion",
      call. = FALSE
    )
  }

  keep <- candidate_grid_indices != selected_grid_index
  sensitivity_text <- candidates$text[keep]
  sensitivity_numeric <- candidates$numeric[keep]
  if (!length(sensitivity_text)) {
    stop(
      "Sensitivity runs are required but no distinct non-primary smoothing value remains",
      call. = FALSE
    )
  }
}

prime_raw <- readLines(prime_recommendations_path, warn = FALSE)
prime_active <- trimws(sub("#.*$", "", prime_raw))
prime_active <- prime_active[nzchar(prime_active)]
if (!length(prime_active)) {
  stop("Prime recommendations file is empty", call. = FALSE)
}
assignment_pattern <- "^(opt|optad|optcvad)[[:space:]]*=[[:space:]]*[1-9][0-9]*$"
detail_names <- c("moredetail", "moredetailad", "moredetailcvad")
valid_prime_line <- grepl(assignment_pattern, prime_active, perl = TRUE) |
  prime_active %in% detail_names
if (any(!valid_prime_line)) {
  stop(
    "Unexpected prime recommendation line(s): ",
    paste(prime_active[!valid_prime_line], collapse = "; "),
    call. = FALSE
  )
}
assignment_keys <- sub(assignment_pattern, "\\1", prime_active[
  grepl(assignment_pattern, prime_active, perl = TRUE)
], perl = TRUE)
for (key in c("opt", "optad", "optcvad")) {
  if (sum(assignment_keys == key) != 1L) {
    stop("Prime recommendations require exactly one assignment for ", key, call. = FALSE)
  }
}
for (key in detail_names) {
  if (sum(prime_active == key) > 1L) {
    stop("Duplicate prime recommendation directive: ", key, call. = FALSE)
  }
}

calibration_lines <- c(
  "mrca = CROWN_PAPILIONOIDEA Parnassius_apollo_ncbi Aeromachus_catocyanea_mydata",
  "min = CROWN_PAPILIONOIDEA 91.5046",
  "max = CROWN_PAPILIONOIDEA 100.8925",
  "",
  "mrca = CROWN_PAPILIONIDAE_NONBARONIA Parnassius_apollo_ncbi Graphium_cloanthus_mydata",
  "min = CROWN_PAPILIONIDAE_NONBARONIA 44.1968",
  "max = CROWN_PAPILIONIDAE_NONBARONIA 52.9473",
  "",
  "mrca = CROWN_HESPERIINAE Aeromachus_catocyanea_mydata Acada_biseriata_kawahara2023",
  "min = CROWN_HESPERIINAE 36.210861",
  "max = CROWN_HESPERIINAE 40.662537"
)
required_calibration_active <- calibration_lines[nzchar(calibration_lines)]
final_seed <- "2026090488"

canonical_smooth_id <- function(value) {
  scientific <- formatC(value, format = "e", digits = 14)
  scientific <- sub("e\\+", "e", scientific)
  scientific <- sub("e(-?)0+", "e\\1", scientific, perl = TRUE)
  scientific <- sub("\\.?0+e", "e", scientific)
  identifier <- gsub("\\.", "p", scientific)
  if (!grepl("^[A-Za-z0-9_.-]+$", identifier)) {
    stop("Could not make a safe smoothing identifier: ", identifier, call. = FALSE)
  }
  identifier
}

variants <- data.frame(
  variant_id = c(
    "final_primary",
    "final_baseline_no_thorough",
    "final_repeat_same_seed"
  ),
  variant_role = c(
    "primary",
    "optimization_diagnostic",
    "same_seed_repeat"
  ),
  smooth = rep(selected_smooth_text, 3L),
  thorough = c("TRUE", "FALSE", "TRUE"),
  seed = rep(final_seed, 3L),
  config_path = c(
    "configs/final.cfg",
    "configs/final_baseline_no_thorough.cfg",
    "configs/final_repeat_same_seed.cfg"
  ),
  tree_path = c(
    "output/T25_three_calibrations_dated.tre",
    "output/T25_three_calibrations_dated_baseline_no_thorough.tre",
    "output/T25_three_calibrations_dated_repeat_same_seed.tre"
  ),
  stage = c(
    "final_primary_thorough",
    "final_baseline_no_thorough",
    "final_repeat_same_seed"
  ),
  stringsAsFactors = FALSE
)

if (length(sensitivity_text)) {
  sensitivity_ids <- vapply(
    sensitivity_numeric, canonical_smooth_id, character(1)
  )
  sensitivity_stages <- paste0("final_sensitivity_", sensitivity_ids)
  sensitivity_rows <- data.frame(
    variant_id = sensitivity_stages,
    variant_role = rep("smoothing_sensitivity", length(sensitivity_text)),
    smooth = sensitivity_text,
    thorough = rep("TRUE", length(sensitivity_text)),
    seed = rep(final_seed, length(sensitivity_text)),
    config_path = paste0("configs/", sensitivity_stages, ".cfg"),
    tree_path = paste0(
      "output/T25_three_calibrations_dated_sensitivity_",
      sensitivity_ids,
      ".tre"
    ),
    stage = sensitivity_stages,
    stringsAsFactors = FALSE
  )
  variants <- rbind(variants, sensitivity_rows)
}

required_manifest_columns <- c(
  "variant_id", "variant_role", "smooth", "thorough", "seed",
  "config_path", "tree_path", "stage"
)
variants <- variants[, required_manifest_columns, drop = FALSE]
allowed_roles <- c(
  "primary", "optimization_diagnostic", "same_seed_repeat",
  "smoothing_sensitivity"
)
if (!all(variants$variant_role %in% allowed_roles) ||
    sum(variants$variant_role == "primary") != 1L ||
    sum(variants$variant_role == "optimization_diagnostic") != 1L ||
    sum(variants$variant_role == "same_seed_repeat") < 1L) {
  stop("Internal error: final variant roles violate the manifest contract", call. = FALSE)
}
if (anyDuplicated(variants$variant_id) || anyDuplicated(variants$stage) ||
    anyDuplicated(variants$config_path) || anyDuplicated(variants$tree_path)) {
  stop("Internal error: final variant identifiers/paths are not unique", call. = FALSE)
}
if (!all(grepl("^[A-Za-z0-9_.-]+$", variants$variant_id)) ||
    !all(grepl("^[A-Za-z0-9_.-]+$", variants$stage))) {
  stop("Internal error: unsafe final variant identifier/stage", call. = FALSE)
}
if (!all(variants$seed == final_seed)) {
  stop("Internal error: final variants do not share the fixed final seed", call. = FALSE)
}
primary <- variants[variants$variant_role == "primary", , drop = FALSE]
repeat_rows <- variants[variants$variant_role == "same_seed_repeat", , drop = FALSE]
if (!all(repeat_rows$smooth == primary$smooth) ||
    !all(repeat_rows$thorough == primary$thorough) ||
    !all(repeat_rows$seed == primary$seed)) {
  stop("Internal error: same-seed repeat does not match primary", call. = FALSE)
}

make_config <- function(outfile, smooth_text, thorough) {
  lines <- c(
    paste0("treefile = ", treefile),
    "numsites = 182682",
    paste0("outfile = ", outfile),
    "",
    calibration_lines,
    "",
    paste0("smooth = ", smooth_text),
    paste0("seed = ", final_seed),
    "",
    prime_active
  )
  if (isTRUE(thorough)) lines <- c(lines, "", "thorough")
  lines
}

config_contents <- vector("list", nrow(variants))
names(config_contents) <- variants$config_path
for (i in seq_len(nrow(variants))) {
  outfile <- file.path(run_root, variants$tree_path[[i]])
  lines <- make_config(
    outfile,
    variants$smooth[[i]],
    identical(variants$thorough[[i]], "TRUE")
  )
  active <- lines[nzchar(lines)]
  if (sum(grepl("^mrca[[:space:]]*=", active)) != 3L ||
      sum(grepl("^min[[:space:]]*=", active)) != 3L ||
      sum(grepl("^max[[:space:]]*=", active)) != 3L ||
      !all(vapply(
        required_calibration_active,
        function(x) sum(active == x) == 1L,
        logical(1)
      ))) {
    stop("Internal error: generated config calibration block is not exact", call. = FALSE)
  }
  if (sum(active == "thorough") != as.integer(variants$thorough[[i]] == "TRUE")) {
    stop("Internal error: generated thorough state is inconsistent", call. = FALSE)
  }
  if (any(grepl(
    "^(prime|cv|randomcv|cvstart|cvstop|cvmultstep|cvoutfile)([[:space:]]*=|$)",
    active
  ))) {
    stop("Internal error: generated final config contains a CV/prime directive", call. = FALSE)
  }
  config_contents[[i]] <- lines
}

manifest_path <- file.path(qa_dir, "final_variant_manifest.tsv")
commands_path <- file.path(qa_dir, "run_final_variants_commands.sh")
config_targets <- file.path(run_root, variants$config_path)
tree_targets <- file.path(run_root, variants$tree_path)
lint_targets <- file.path(
  qa_dir,
  paste0("config_lint_", sub("\\.cfg$", "", basename(config_targets)), ".tsv")
)
log_targets <- unlist(lapply(
  variants$stage,
  function(stage) file.path(
    log_dir,
    paste0(stage, c(".stdout", ".stderr", ".metadata.tsv", ".resources.txt"))
  )
), use.names = FALSE)
all_collision_targets <- c(
  config_targets, manifest_path, commands_path, tree_targets, lint_targets,
  log_targets
)
existing_targets <- all_collision_targets[file.exists(all_collision_targets)]
if (length(existing_targets)) {
  stop(
    "Refusing to overwrite/reuse existing final-stage artifact(s): ",
    paste(existing_targets, collapse = ", "),
    call. = FALSE
  )
}

manifest_lines <- c(
  paste(required_manifest_columns, collapse = "\t"),
  vapply(
    seq_len(nrow(variants)),
    function(i) paste(
      unlist(variants[i, required_manifest_columns], use.names = FALSE),
      collapse = "\t"
    ),
    character(1)
  )
)
if (any(grepl("[\t\r\n]", unlist(variants), perl = TRUE))) {
  stop("Internal error: tab/newline in final-variant manifest field", call. = FALSE)
}

shell_quote <- function(x) shQuote(x, type = "sh")
command_lines <- c(
  "#!/usr/bin/env bash",
  "set -euo pipefail",
  "",
  "# Generated only as a command list; this generator never starts treePL."
)
for (i in seq_len(nrow(variants))) {
  phase <- if (identical(variants$variant_role[[i]], "optimization_diagnostic")) {
    "final_baseline"
  } else {
    "final"
  }
  command_lines <- c(
    command_lines,
    paste(
      "Rscript",
      shell_quote(lint_script),
      shell_quote(config_targets[[i]]),
      phase,
      shell_quote(run_root),
      shell_quote(lint_targets[[i]])
    )
  )
}
command_lines <- c(command_lines, "")
for (i in seq_len(nrow(variants))) {
  command_lines <- c(
    command_lines,
    paste(
      shell_quote(wrapper),
      shell_quote(run_root),
      shell_quote(binary),
      shell_quote(config_targets[[i]]),
      shell_quote(variants$stage[[i]])
    )
  )
}

write_exclusive <- function(lines, path) {
  if (file.exists(path)) {
    stop("Refusing to overwrite existing output: ", path, call. = FALSE)
  }
  con <- file(path, open = "wx", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(lines, con, useBytes = TRUE)
}

for (i in seq_len(nrow(variants))) {
  write_exclusive(config_contents[[i]], config_targets[[i]])
}
write_exclusive(manifest_lines, manifest_path)
write_exclusive(command_lines, commands_path)

cat("Generated final treePL configs:", nrow(variants), "\n")
cat("Selected smoothing:", selected_smooth_text, "\n")
cat("Sensitivity configs:", length(sensitivity_text), "\n")
cat("Manifest:", manifest_path, "\n")
cat("Command list (not executed):", commands_path, "\n")
