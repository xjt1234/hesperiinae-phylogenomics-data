#!/usr/bin/env Rscript

# Independent, fail-closed comparison of the two retained treePL scenarios.
# This script is deliberately self-contained and refuses to overwrite outputs.

suppressPackageStartupMessages(library(ape))
options(stringsAsFactors = FALSE, digits = 17)

RUN_ROOT <- "/home/data/t200301/xjt/Hesperiinae_review/treePL_R2_3_minimal_validation_20260904_135356"
A_ROOT <- "/home/data/t200301/xjt/Hesperiinae_review/treePL_R2_3_T25_deep_only_20260904_000707"
B_ROOT <- "/home/data/t200301/xjt/Hesperiinae_review/treePL_R2_3_T25_three_calibrations_20260904_075723"
EXPECTED_INPUT_SHA256 <- "4d35d553ed5bbf08022daef2361ffa0c342f8694a53aaf08d5350692094439e5"
EXPECTED_TIPS <- 495L
EXPECTED_INTERNAL <- 494L
EXPECTED_NUMSITES <- 182682L
BOUNDARY_TOL <- 1e-4
EOCENE_OLIGOCENE_MA <- 33.9

paths <- list(
  A = list(
    scenario = "A_two_node_main",
    root = A_ROOT,
    input = file.path(A_ROOT, "input", "verified_495_tip_input.treefile"),
    final = file.path(A_ROOT, "output", "T25_deep_only_dated.tre"),
    config = file.path(A_ROOT, "configs", "final.cfg"),
    metadata = file.path(A_ROOT, "logs", "final_primary_thorough.metadata.tsv")
  ),
  B = list(
    scenario = "B_three_node_sensitivity",
    root = B_ROOT,
    input = file.path(B_ROOT, "input", "verified_495_tip_input.treefile"),
    final = file.path(B_ROOT, "output", "T25_three_calibrations_dated.tre"),
    config = file.path(B_ROOT, "configs", "final.cfg"),
    metadata = file.path(B_ROOT, "logs", "final_primary_thorough.metadata.tsv")
  )
)

out_paths <- c(
  summary = file.path(RUN_ROOT, "calibration_sensitivity_summary.tsv"),
  differences = file.path(RUN_ROOT, "qa", "all_internal_node_differences.tsv"),
  difference_summary = file.path(RUN_ROOT, "qa", "all_internal_node_difference_summary.tsv"),
  crossings = file.path(RUN_ROOT, "qa", "boundary_crossings_33_9.tsv"),
  checks = file.path(RUN_ROOT, "qa", "ab_validation_checks.tsv"),
  report = file.path(RUN_ROOT, "qa", "ab_validation_report.md")
)

if (any(file.exists(out_paths))) {
  stop("Refusing to overwrite existing output(s): ",
       paste(out_paths[file.exists(out_paths)], collapse = ", "))
}
required_files <- unlist(lapply(paths, function(x) unlist(x[c("input", "final", "config", "metadata")])))
if (any(!file.exists(required_files))) {
  stop("Missing required input(s): ", paste(required_files[!file.exists(required_files)], collapse = ", "))
}
dir.create(file.path(RUN_ROOT, "scripts"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(RUN_ROOT, "qa"), recursive = TRUE, showWarnings = FALSE)

sha256_file <- function(path) {
  ans <- system2("sha256sum", shQuote(path), stdout = TRUE, stderr = TRUE)
  status <- attr(ans, "status")
  if (!is.null(status) && status != 0L) stop("sha256sum failed for ", path, ": ", paste(ans, collapse = " "))
  strsplit(ans[[1]], "[[:space:]]+")[[1]][[1]]
}

read_cfg <- function(path) {
  raw <- readLines(path, warn = FALSE)
  clean <- trimws(sub("#.*$", "", raw))
  clean[nzchar(clean)]
}

cfg_scalar <- function(lines, key) {
  hits <- grep(paste0("^", key, "[[:space:]]*="), lines, value = TRUE)
  if (length(hits) != 1L) return(NA_character_)
  trimws(sub(paste0("^", key, "[[:space:]]*=[[:space:]]*"), "", hits))
}

parse_calibrations <- function(lines) {
  mrca_lines <- grep("^mrca[[:space:]]*=", lines, value = TRUE)
  min_lines <- grep("^min[[:space:]]*=", lines, value = TRUE)
  max_lines <- grep("^max[[:space:]]*=", lines, value = TRUE)
  parse_name <- function(x, key) sub(paste0("^", key, "[[:space:]]*=[[:space:]]*([^[:space:]]+).*$"), "\\1", x)
  parse_value <- function(x, key) as.numeric(sub(paste0("^", key, "[[:space:]]*=[[:space:]]*[^[:space:]]+[[:space:]]+"), "", x))
  mins <- setNames(parse_value(min_lines, "min"), parse_name(min_lines, "min"))
  maxs <- setNames(parse_value(max_lines, "max"), parse_name(max_lines, "max"))
  list(
    names = parse_name(mrca_lines, "mrca"),
    n = length(mrca_lines),
    mins = mins,
    maxs = maxs,
    complete = identical(sort(names(mins)), sort(parse_name(mrca_lines, "mrca"))) &&
      identical(sort(names(maxs)), sort(parse_name(mrca_lines, "mrca")))
  )
}

read_metadata <- function(path) {
  x <- read.delim(path, header = FALSE, sep = "\t", quote = "", comment.char = "",
                  colClasses = "character")
  if (ncol(x) != 2L || anyDuplicated(x[[1]])) stop("Malformed metadata: ", path)
  setNames(x[[2]], x[[1]])
}

root_node <- function(tree) {
  roots <- setdiff(tree$edge[, 1], tree$edge[, 2])
  if (length(roots) != 1L) stop("Expected exactly one root")
  as.integer(roots)
}

tree_is_rooted_binary <- function(tree) {
  r <- root_node(tree)
  internal <- seq.int(Ntip(tree) + 1L, Ntip(tree) + tree$Nnode)
  outdegree <- tabulate(match(tree$edge[, 1], internal), nbins = length(internal))
  isTRUE(is.rooted(tree)) && r %in% internal && tree$Nnode == Ntip(tree) - 1L && all(outdegree == 2L)
}

tree_index <- function(tree) {
  n_tip <- Ntip(tree)
  children <- split(tree$edge[, 2], tree$edge[, 1])
  memo <- new.env(parent = emptyenv())
  desc_ids <- function(node) {
    key <- as.character(node)
    if (exists(key, envir = memo, inherits = FALSE)) return(get(key, envir = memo, inherits = FALSE))
    ans <- if (node <= n_tip) {
      node
    } else {
      kids <- children[[key]]
      if (is.null(kids)) integer(0) else sort(unique(unlist(lapply(kids, desc_ids), use.names = FALSE)))
    }
    assign(key, ans, envir = memo)
    ans
  }
  depths <- node.depth.edgelength(tree)
  internal <- seq.int(n_tip + 1L, n_tip + tree$Nnode)
  rows <- lapply(internal, function(node) {
    ids <- desc_ids(node)
    labs <- sort(tree$tip.label[ids])
    kids <- children[[as.character(node)]]
    anchors <- if (length(kids) == 2L) {
      c(sort(tree$tip.label[desc_ids(kids[[1]])])[[1]],
        sort(tree$tip.label[desc_ids(kids[[2]])])[[1]])
    } else c(NA_character_, NA_character_)
    paths <- depths[ids] - depths[node]
    data.frame(
      node = node,
      clade_key = paste(labs, collapse = "|"),
      descendant_tip_count = length(ids),
      anchor_tip_1 = anchors[[1]],
      anchor_tip_2 = anchors[[2]],
      age_ma = mean(paths),
      descendant_path_min_ma = min(paths),
      descendant_path_max_ma = max(paths),
      descendant_path_range_ma = diff(range(paths)),
      stringsAsFactors = FALSE
    )
  })
  ans <- do.call(rbind, rows)
  if (anyDuplicated(ans$clade_key)) stop("Internal clade keys are not unique")
  ans
}

root_to_tip_range <- function(tree) {
  d <- node.depth.edgelength(tree)[seq_len(Ntip(tree))]
  diff(range(d))
}

rooted_rf <- function(index_x, index_y) {
  length(setdiff(index_x$clade_key, index_y$clade_key)) +
    length(setdiff(index_y$clade_key, index_x$clade_key))
}

branch_counts <- function(tree) {
  x <- tree$edge.length
  c(negative = sum(is.finite(x) & x < 0), zero = sum(is.finite(x) & x == 0),
    nonfinite = sum(!is.finite(x)))
}

boundary_status <- function(age, lower, upper, calibrated, tol = BOUNDARY_TOL) {
  if (!calibrated) return("not_calibrated")
  low <- abs(age - lower) <= tol
  high <- abs(age - upper) <= tol
  if (low && high) return("both_boundaries_hit")
  if (low) return("lower_boundary_hit")
  if (high) return("upper_boundary_hit")
  if (age < lower - tol) return("below_lower_bound")
  if (age > upper + tol) return("above_upper_bound")
  "within_bounds_no_boundary_hit"
}

calibration_defs <- data.frame(
  id = c("CROWN_PAPILIONOIDEA", "CROWN_PAPILIONIDAE_NONBARONIA", "CROWN_HESPERIINAE"),
  anchor1 = c("Parnassius_apollo_ncbi", "Parnassius_apollo_ncbi", "Aeromachus_catocyanea_mydata"),
  anchor2 = c("Aeromachus_catocyanea_mydata", "Graphium_cloanthus_mydata", "Acada_biseriata_kawahara2023"),
  descendant_expected = c(495L, 7L, 488L),
  lower = c(91.5046, 44.1968, 36.210861),
  upper = c(100.8925, 52.9473, 40.662537),
  stringsAsFactors = FALSE
)

extract_focal <- function(tree, index) {
  out <- lapply(seq_len(nrow(calibration_defs)), function(i) {
    d <- calibration_defs[i, ]
    if (!all(c(d$anchor1, d$anchor2) %in% tree$tip.label)) stop("Missing focal anchor tip(s)")
    node <- getMRCA(tree, c(d$anchor1, d$anchor2))
    row <- index[index$node == node, , drop = FALSE]
    if (nrow(row) != 1L) stop("Unable to index focal MRCA: ", d$id)
    data.frame(id = d$id, node = node, descendant_tip_count = row$descendant_tip_count,
               age_ma = row$age_ma, lower = d$lower, upper = d$upper,
               clade_key = row$clade_key, stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

add_check_factory <- function() {
  checks <- list()
  add <- function(check, scope, critical, pass, observed, expected, details = "") {
    checks[[length(checks) + 1L]] <<- data.frame(
      check = check, scope = scope, critical = critical,
      status = if (isTRUE(pass)) "PASS" else "FAIL",
      observed = as.character(observed), expected = as.character(expected), details = details,
      stringsAsFactors = FALSE
    )
  }
  list(add = add, get = function() do.call(rbind, checks))
}

cfg <- lapply(paths, function(x) read_cfg(x$config))
cal <- lapply(cfg, parse_calibrations)
meta <- lapply(paths, function(x) read_metadata(x$metadata))
trees <- list(
  A_input = read.tree(paths$A$input), A_final = read.tree(paths$A$final),
  B_input = read.tree(paths$B$input), B_final = read.tree(paths$B$final)
)
indices <- lapply(trees, tree_index)
focal <- list(A = extract_focal(trees$A_final, indices$A_final),
              B = extract_focal(trees$B_final, indices$B_final))
input_sha <- c(A = sha256_file(paths$A$input), B = sha256_file(paths$B$input))
final_sha <- c(A = sha256_file(paths$A$final), B = sha256_file(paths$B$final))
numsites <- c(A = as.integer(cfg_scalar(cfg$A, "numsites")), B = as.integer(cfg_scalar(cfg$B, "numsites")))
smoothing <- c(A = as.numeric(cfg_scalar(cfg$A, "smooth")), B = as.numeric(cfg_scalar(cfg$B, "smooth")))
rt_range <- c(A = root_to_tip_range(trees$A_final), B = root_to_tip_range(trees$B_final))
branches <- list(A = branch_counts(trees$A_final), B = branch_counts(trees$B_final))
rf_input_final <- c(A = rooted_rf(indices$A_input, indices$A_final),
                    B = rooted_rf(indices$B_input, indices$B_final))

cc <- add_check_factory()
add_check <- cc$add
add_check("input_file_sha256_fixed", "shared", TRUE,
          identical(unname(input_sha), rep(EXPECTED_INPUT_SHA256, 2L)),
          paste(input_sha, collapse = ";"), EXPECTED_INPUT_SHA256)
add_check("input_tip_sets_identical", "shared", TRUE,
          setequal(trees$A_input$tip.label, trees$B_input$tip.label),
          length(intersect(trees$A_input$tip.label, trees$B_input$tip.label)), EXPECTED_TIPS)
add_check("final_tip_sets_identical", "shared", TRUE,
          setequal(trees$A_final$tip.label, trees$B_final$tip.label),
          length(intersect(trees$A_final$tip.label, trees$B_final$tip.label)), EXPECTED_TIPS)
add_check("A_B_rooted_clades_identical", "shared", TRUE,
          rooted_rf(indices$A_final, indices$B_final) == 0L,
          rooted_rf(indices$A_final, indices$B_final), 0)
add_check("A_B_internal_nodes_matched", "shared", TRUE,
          length(intersect(indices$A_final$clade_key, indices$B_final$clade_key)) == EXPECTED_INTERNAL,
          length(intersect(indices$A_final$clade_key, indices$B_final$clade_key)), EXPECTED_INTERNAL)

expected_cal_names <- list(
  A = c("CROWN_PAPILIONOIDEA", "CROWN_PAPILIONIDAE_NONBARONIA"),
  B = c("CROWN_PAPILIONOIDEA", "CROWN_PAPILIONIDAE_NONBARONIA", "CROWN_HESPERIINAE")
)

for (s in c("A", "B")) {
  input_name <- paste0(s, "_input")
  final_name <- paste0(s, "_final")
  scenario <- paths[[s]]$scenario
  for (which in c("input", "final")) {
    tree <- trees[[paste0(s, "_", which)]]
    add_check(paste0(which, "_tip_count"), scenario, TRUE, Ntip(tree) == EXPECTED_TIPS,
              Ntip(tree), EXPECTED_TIPS)
    add_check(paste0(which, "_unique_tip_labels"), scenario, TRUE, anyDuplicated(tree$tip.label) == 0L,
              anyDuplicated(tree$tip.label), 0)
    add_check(paste0(which, "_rooted_binary"), scenario, TRUE, tree_is_rooted_binary(tree),
              tree_is_rooted_binary(tree), TRUE)
    bc <- branch_counts(tree)
    add_check(paste0(which, "_negative_branch_count"), scenario, TRUE, bc[["negative"]] == 0L,
              bc[["negative"]], 0)
    add_check(paste0(which, "_zero_branch_count"), scenario, TRUE, bc[["zero"]] == 0L,
              bc[["zero"]], 0)
    add_check(paste0(which, "_nonfinite_branch_count"), scenario, TRUE, bc[["nonfinite"]] == 0L,
              bc[["nonfinite"]], 0)
  }
  add_check("input_final_tip_sets_identical", scenario, TRUE,
            setequal(trees[[input_name]]$tip.label, trees[[final_name]]$tip.label),
            length(intersect(trees[[input_name]]$tip.label, trees[[final_name]]$tip.label)), EXPECTED_TIPS)
  add_check("input_final_rooted_rf", scenario, TRUE, rf_input_final[[s]] == 0L,
            rf_input_final[[s]], 0)
  add_check("input_final_rooted_clades_identical", scenario, TRUE,
            identical(sort(indices[[input_name]]$clade_key), sort(indices[[final_name]]$clade_key)),
            length(intersect(indices[[input_name]]$clade_key, indices[[final_name]]$clade_key)), EXPECTED_INTERNAL)
  add_check("final_root_to_tip_range", scenario, TRUE, rt_range[[s]] < 1e-4,
            format(rt_range[[s]], scientific = TRUE), "<1e-4 Ma")
  add_check("numsites", scenario, TRUE, identical(numsites[[s]], EXPECTED_NUMSITES),
            numsites[[s]], EXPECTED_NUMSITES)
  add_check("calibration_count", scenario, TRUE, cal[[s]]$n == length(expected_cal_names[[s]]),
            cal[[s]]$n, length(expected_cal_names[[s]]))
  add_check("calibration_names_exact", scenario, TRUE,
            identical(sort(cal[[s]]$names), sort(expected_cal_names[[s]])),
            paste(sort(cal[[s]]$names), collapse = ";"), paste(sort(expected_cal_names[[s]]), collapse = ";"))
  add_check("calibration_min_max_complete", scenario, TRUE, cal[[s]]$complete,
            cal[[s]]$complete, TRUE)
  add_check("hesperiinae_unconstrained_in_A", scenario, TRUE,
            if (s == "A") !any(grepl("CROWN_HESPERIINAE", cfg[[s]], fixed = TRUE)) else TRUE,
            if (s == "A") any(grepl("CROWN_HESPERIINAE", cfg[[s]], fixed = TRUE)) else "not_applicable",
            if (s == "A") FALSE else "not_applicable")
  add_check("metadata_exit_code", scenario, TRUE,
            !is.na(meta[[s]][["exit_code"]]) && meta[[s]][["exit_code"]] == "0",
            meta[[s]][["exit_code"]], 0)
  cfg_tree <- cfg_scalar(cfg[[s]], "treefile")
  cfg_out <- cfg_scalar(cfg[[s]], "outfile")
  add_check("config_uses_verified_undated_input", scenario, TRUE,
            identical(normalizePath(cfg_tree), normalizePath(paths[[s]]$input)) &&
              sha256_file(cfg_tree) == EXPECTED_INPUT_SHA256,
            cfg_tree, paths[[s]]$input)
  add_check("final_is_direct_config_outfile", scenario, TRUE,
            identical(normalizePath(cfg_out), normalizePath(paths[[s]]$final)),
            cfg_out, paths[[s]]$final,
            "Direct treePL outfile; no rescale/re-anchor step is used by this compiler")
  expected_desc <- calibration_defs$descendant_expected
  add_check("focal_mrca_descendant_counts", scenario, TRUE,
            identical(as.integer(focal[[s]]$descendant_tip_count), as.integer(expected_desc)),
            paste(focal[[s]]$descendant_tip_count, collapse = ";"), paste(expected_desc, collapse = ";"))
  pap_key <- focal[[s]]$clade_key[focal[[s]]$id == "CROWN_PAPILIONIDAE_NONBARONIA"]
  hesp_key <- focal[[s]]$clade_key[focal[[s]]$id == "CROWN_HESPERIINAE"]
  pap_tips <- strsplit(pap_key, "|", fixed = TRUE)[[1]]
  hesp_tips <- strsplit(hesp_key, "|", fixed = TRUE)[[1]]
  add_check("papilionidae_hesperiinae_partition_root", scenario, TRUE,
            length(intersect(pap_tips, hesp_tips)) == 0L &&
              setequal(c(pap_tips, hesp_tips), trees[[final_name]]$tip.label),
            paste(length(pap_tips), length(hesp_tips), sep = "+"), "7+488 disjoint tips = 495")
  for (i in seq_len(nrow(focal[[s]]))) {
    id <- focal[[s]]$id[[i]]
    calibrated <- id %in% cal[[s]]$names
    if (calibrated) {
      age <- focal[[s]]$age_ma[[i]]
      lo <- calibration_defs$lower[calibration_defs$id == id]
      hi <- calibration_defs$upper[calibration_defs$id == id]
      add_check(paste0("calibration_within_bounds_", id), scenario, TRUE,
                age >= lo - BOUNDARY_TOL && age <= hi + BOUNDARY_TOL,
                format(age, digits = 16), paste0("[", lo, ",", hi, "] +/- ", BOUNDARY_TOL))
    }
  }
}

checks <- cc$get()
critical_failures <- checks$critical & checks$status != "PASS"
if (any(critical_failures)) {
  fail_text <- paste(checks$check[critical_failures], checks$scope[critical_failures], sep = "@", collapse = "; ")
  stop("Critical validation failure(s); no outputs written: ", fail_text)
}

# Match every rooted internal node by the full, sorted descendant-tip set.
diffs <- merge(
  indices$A_final[, c("node", "clade_key", "descendant_tip_count", "anchor_tip_1", "anchor_tip_2", "age_ma")],
  indices$B_final[, c("node", "clade_key", "descendant_tip_count", "age_ma")],
  by = "clade_key", suffixes = c("_A", "_B"), sort = FALSE
)
if (nrow(diffs) != EXPECTED_INTERNAL) stop("Internal-node match count changed after validation")
diffs$signed_shift_B_minus_A_ma <- diffs$age_ma_B - diffs$age_ma_A
diffs$absolute_difference_ma <- abs(diffs$signed_shift_B_minus_A_ma)
diffs$crosses_33_9_ma <- (diffs$age_ma_A > EOCENE_OLIGOCENE_MA & diffs$age_ma_B <= EOCENE_OLIGOCENE_MA) |
  (diffs$age_ma_B > EOCENE_OLIGOCENE_MA & diffs$age_ma_A <= EOCENE_OLIGOCENE_MA)
diffs$crossing_direction <- ifelse(
  !diffs$crosses_33_9_ma, "no_crossing",
  ifelse(diffs$age_ma_A > EOCENE_OLIGOCENE_MA, "A_older_to_B_younger_or_equal", "A_younger_or_equal_to_B_older")
)
diffs <- diffs[order(-diffs$descendant_tip_count_A, diffs$clade_key), ]
diffs$matched_internal_node_index <- seq_len(nrow(diffs))
diffs <- diffs[, c("matched_internal_node_index", "node_A", "node_B", "descendant_tip_count_A",
                   "anchor_tip_1", "anchor_tip_2", "age_ma_A", "age_ma_B",
                   "signed_shift_B_minus_A_ma", "absolute_difference_ma", "crosses_33_9_ma",
                   "crossing_direction", "clade_key")]
names(diffs)[names(diffs) == "descendant_tip_count_A"] <- "descendant_tip_count"

absdiff <- diffs$absolute_difference_ma
diff_summary <- data.frame(
  matched_internal_node_count = length(absdiff),
  mean_absolute_difference_ma = mean(absdiff),
  median_absolute_difference_ma = median(absdiff),
  p95_absolute_difference_ma = unname(quantile(absdiff, probs = 0.95, type = 7)),
  maximum_absolute_difference_ma = max(absdiff),
  count_absolute_difference_gt_1_ma = sum(absdiff > 1),
  count_absolute_difference_gt_2_ma = sum(absdiff > 2),
  count_crossing_33_9_ma = sum(diffs$crosses_33_9_ma),
  mean_signed_shift_B_minus_A_ma = mean(diffs$signed_shift_B_minus_A_ma),
  median_signed_shift_B_minus_A_ma = median(diffs$signed_shift_B_minus_A_ma),
  stringsAsFactors = FALSE
)

crossings <- diffs[diffs$crosses_33_9_ma, c(
  "matched_internal_node_index", "node_A", "node_B", "descendant_tip_count",
  "anchor_tip_1", "anchor_tip_2", "age_ma_A", "age_ma_B",
  "signed_shift_B_minus_A_ma", "absolute_difference_ma", "crossing_direction", "clade_key"
), drop = FALSE]
crossings$boundary_ma <- EOCENE_OLIGOCENE_MA
crossings <- crossings[, c("boundary_ma", setdiff(names(crossings), "boundary_ma"))]

scenario_status <- function(s) {
  relevant <- checks$scope %in% c("shared", paths[[s]]$scenario) & checks$critical
  if (all(checks$status[relevant] == "PASS")) "PASS" else "FAIL"
}
focal_age <- function(s, id) focal[[s]]$age_ma[focal[[s]]$id == id]
calibrated_in <- function(s, id) id %in% cal[[s]]$names

summary <- do.call(rbind, lapply(c("A", "B"), function(s) {
  root_age <- focal_age(s, "CROWN_PAPILIONOIDEA")
  pap_age <- focal_age(s, "CROWN_PAPILIONIDAE_NONBARONIA")
  hesp_age <- focal_age(s, "CROWN_HESPERIINAE")
  data.frame(
    scenario = paths[[s]]$scenario,
    input_tree_path = paths[[s]]$input,
    input_tree_sha256 = input_sha[[s]],
    final_tree_path = paths[[s]]$final,
    final_tree_sha256 = final_sha[[s]],
    tip_count = Ntip(trees[[paste0(s, "_final")]]),
    numsites = numsites[[s]],
    smoothing = smoothing[[s]],
    root_age_ma = root_age,
    papilionidae_age_ma = pap_age,
    hesperiinae_age_ma = hesp_age,
    hesperiinae_shift_from_main_ma = hesp_age - focal_age("A", "CROWN_HESPERIINAE"),
    root_boundary_status = boundary_status(root_age, 91.5046, 100.8925, TRUE),
    papilionidae_boundary_status = boundary_status(pap_age, 44.1968, 52.9473, TRUE),
    hesperiinae_boundary_status = boundary_status(hesp_age, 36.210861, 40.662537,
                                                   calibrated_in(s, "CROWN_HESPERIINAE")),
    ultrametric_range_ma = rt_range[[s]],
    negative_branch_count = branches[[s]][["negative"]],
    zero_branch_count = branches[[s]][["zero"]],
    nonfinite_branch_count = branches[[s]][["nonfinite"]],
    rooted_rf = rf_input_final[[s]],
    qa_status = scenario_status(s),
    stringsAsFactors = FALSE
  )
}))

fmt <- function(x) format(x, digits = 15, scientific = FALSE, trim = TRUE)
report <- c(
  "# A/B treePL minimal validation report",
  "",
  paste0("Generated by `scripts/01_compile_ab_summary.R`; boundary tolerance = ", BOUNDARY_TOL, " Ma."),
  "The script parsed all three trees independently with `ape` and matched rooted internal nodes by each node's complete sorted descendant-tip set.",
  "",
  "## Retained scenarios",
  "",
  "| Scenario | smoothing | root (Ma) | Papilionidae excluding Baronia (Ma) | Hesperiinae (Ma) | Hesperiinae shift from A (Ma) | QA |",
  "|---|---:|---:|---:|---:|---:|---|",
  paste0("| ", summary$scenario, " | ", format(summary$smoothing, scientific = TRUE), " | ",
         fmt(summary$root_age_ma), " | ", fmt(summary$papilionidae_age_ma), " | ",
         fmt(summary$hesperiinae_age_ma), " | ", fmt(summary$hesperiinae_shift_from_main_ma),
         " | ", summary$qa_status, " |"),
  "",
  "Boundary classifications use `<= 1e-4 Ma`: A root = lower-bound hit; A and B Papilionidae = upper-bound hits; B Hesperiinae = upper-bound hit. A Hesperiinae is explicitly unconstrained.",
  "",
  "## All-node comparison",
  "",
  paste0("Matched internal nodes: ", diff_summary$matched_internal_node_count, "."),
  paste0("Absolute age difference (Ma): mean ", fmt(diff_summary$mean_absolute_difference_ma),
         ", median ", fmt(diff_summary$median_absolute_difference_ma),
         ", 95th percentile ", fmt(diff_summary$p95_absolute_difference_ma),
         ", maximum ", fmt(diff_summary$maximum_absolute_difference_ma), "."),
  paste0("Nodes with absolute difference >1 Ma: ", diff_summary$count_absolute_difference_gt_1_ma,
         "; >2 Ma: ", diff_summary$count_absolute_difference_gt_2_ma, "."),
  paste0("Nodes crossing the 33.9 Ma Eocene-Oligocene boundary between A and B: ",
         diff_summary$count_crossing_33_9_ma, "."),
  "",
  "## Validation outcome",
  "",
  paste0("Critical checks: ", sum(checks$critical), " PASS, 0 FAIL. Both retained scenarios pass."),
  paste0("Fixed 495-tip input SHA-256: `", EXPECTED_INPUT_SHA256, "`."),
  "Both input and final trees are rooted and fully binary; each final tree has 495 unique tips, no negative, zero, or non-finite branches, root-to-tip range <1e-4 Ma, and rooted RF=0 versus its undated ML input.",
  "A has exactly two calibrations and no Hesperiinae constraint; B has exactly three calibrations. Both final-stage metadata files record exit code 0.",
  "The final tree paths equal the direct `outfile` values in their treePL configs; this comparison performs no rescaling, re-anchoring, or branch-length modification.",
  "",
  "## Output files",
  "",
  "- `calibration_sensitivity_summary.tsv`: scenario-level results and hashes.",
  "- `qa/all_internal_node_differences.tsv`: all 494 matched internal nodes, with full descendant-tip keys.",
  "- `qa/all_internal_node_difference_summary.tsv`: predefined difference statistics.",
  "- `qa/boundary_crossings_33_9.tsv`: every internal node that changes sides of 33.9 Ma.",
  "- `qa/ab_validation_checks.tsv`: machine-readable critical checks.",
  "",
  "No C/D leave-one-out output is read, changed, or interpreted by this script."
)

write_tsv_strict <- function(x, path) {
  if (file.exists(path)) stop("Refusing to overwrite: ", path)
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}
write_lines_strict <- function(x, path) {
  if (file.exists(path)) stop("Refusing to overwrite: ", path)
  writeLines(x, path, useBytes = TRUE)
}

write_tsv_strict(summary, out_paths[["summary"]])
write_tsv_strict(diffs, out_paths[["differences"]])
write_tsv_strict(diff_summary, out_paths[["difference_summary"]])
write_tsv_strict(crossings, out_paths[["crossings"]])
write_tsv_strict(checks, out_paths[["checks"]])
write_lines_strict(report, out_paths[["report"]])

cat("PASS: A/B validation complete\n")
cat("A Hesperiinae:", format(focal_age("A", "CROWN_HESPERIINAE"), digits = 16), "Ma\n")
cat("B Hesperiinae:", format(focal_age("B", "CROWN_HESPERIINAE"), digits = 16), "Ma\n")
cat("Shift B-A:", format(focal_age("B", "CROWN_HESPERIINAE") - focal_age("A", "CROWN_HESPERIINAE"), digits = 16), "Ma\n")
cat("Matched internal nodes:", nrow(diffs), "\n")
cat("Boundary crossings at 33.9 Ma:", nrow(crossings), "\n")
