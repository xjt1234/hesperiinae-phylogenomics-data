#!/usr/bin/env Rscript

# Rebuild the exact environmental series used by the historical local RPANDA
# analysis, extending it with observed CENOGRID points (never extrapolation) to
# the automatically measured age of the Scenario A 417-tip chronogram.

suppressPackageStartupMessages(library(ape))
options(stringsAsFactors = FALSE, digits = 17, scipen = 999)

args_full <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", grep("^--file=", args_full, value = TRUE))
if (length(script_arg) != 1L) stop("Cannot determine script path")
run_dir <- normalizePath(file.path(dirname(script_arg), ".."), mustWork = TRUE)

tree_file <- file.path(run_dir, "01_tree", "T25_deep_only_dated.Hesperiinae_417sp_RPANDA.tre")
raw_file <- file.path(run_dir, "02_environment", "input", "CENOGRID_S34.tsv")
old_file <- file.path(run_dir, "02_environment", "input", "temperature_surface.tsv")
out_file <- file.path(run_dir, "02_environment", "temperature_surface_to_tree_age.tsv")
qa_file <- file.path(run_dir, "02_environment", "temperature_overlap_QA.tsv")
report_file <- file.path(run_dir, "02_environment", "temperature_QA_report.md")

needed <- c(tree_file, raw_file, old_file)
if (any(!file.exists(needed))) {
  stop("Missing required input(s):\n", paste(needed[!file.exists(needed)], collapse = "\n"))
}

tree <- read.tree(tree_file)
tip_depths <- node.depth.edgelength(tree)[seq_len(Ntip(tree))]
total_time <- max(tip_depths)

raw <- read.delim(raw_file, skip = 22L, check.names = FALSE,
                  quote = "", comment.char = "", na.strings = c("", "NA"))
time_col <- "Tuned time [Ma]"
d18o_col <- "Foram benth δ18O [‰ PDB] (ISOBENd18oLOESSsmoothLongTerm)"
if (!all(c(time_col, d18o_col) %in% names(raw))) {
  stop("Expected CENOGRID columns not found. Available columns:\n",
       paste(names(raw), collapse = "\n"))
}

time <- as.numeric(raw[[time_col]])
d18o <- as.numeric(raw[[d18o_col]])
if (any(!is.finite(time)) || any(!is.finite(d18o))) stop("Non-finite CENOGRID time or d18O")
if (is.unsorted(time, strictly = TRUE)) stop("CENOGRID time must be strictly increasing")
if (min(time) > 0 || max(time) < total_time) stop("Raw CENOGRID does not cover [0, tree age]")

# Exact transformation recovered from every shared point of the old local
# temperature files. The source column is already the published long-term
# LOESS series; df=80 smoothing is applied later inside fit_env_bd().
deep_temperature <- ifelse(d18o > 1.75, 8.5 - 2 * d18o, 12 - 4 * d18o)
present_idx <- which(time == 0)
if (length(present_idx) != 1L) stop("Expected exactly one CENOGRID time-zero row")
deep_temperature_present <- deep_temperature[present_idx]
surface_temperature <- ifelse(
  time <= 2.58,
  2 * deep_temperature + 12.25,
  ifelse(
    time <= 5.33,
    2.5 * deep_temperature + 12.15,
    14.15 + (deep_temperature - deep_temperature_present)
  )
)

# Keep the first observed point at or beyond total_time so the fitted spline
# covers the complete tree age. This is an observation, not an added endpoint.
last_idx <- which(time >= total_time)[1]
keep_idx <- seq_len(last_idx)
env <- data.frame(
  time = time[keep_idx],
  temp = surface_temperature[keep_idx],
  deep_temperature = deep_temperature[keep_idx],
  benthic_d18O_longterm_LOESS = d18o[keep_idx],
  stringsAsFactors = FALSE
)

if (is.unsorted(env$time, strictly = TRUE)) stop("Output time is not strictly increasing")
if (any(!is.finite(as.matrix(env)))) stop("Output contains NA/NaN/Inf")
if (min(env$time) != 0 || max(env$time) < total_time) stop("Output does not cover full tree age")

old <- read.delim(old_file, check.names = FALSE, stringsAsFactors = FALSE)
if (ncol(old) < 2L) stop("Historical temperature table must have at least two columns")
old <- data.frame(time = as.numeric(old[[1]]), old_temp = as.numeric(old[[2]]))
shared <- merge(old, env[c("time", "temp")], by = "time", all = FALSE, sort = TRUE)
if (nrow(shared) == 0L) stop("No shared time points with historical temperature table")
diff <- shared$temp - shared$old_temp

qa <- data.frame(
  metric = c(
    "tree_total_time_Ma", "raw_min_time_Ma", "raw_max_time_Ma",
    "output_last_observed_time_Ma", "output_points", "shared_points",
    "maximum_absolute_difference_C", "RMSE_C", "time_zero_definition",
    "time_direction", "RPANDA_environment_spline_df", "constant_extrapolation_used"
  ),
  value = c(
    format(total_time, digits = 17), format(min(time), digits = 17),
    format(max(time), digits = 17), format(max(env$time), digits = 17),
    nrow(env), nrow(shared), format(max(abs(diff)), digits = 17),
    format(sqrt(mean(diff^2)), digits = 17), "present", "increases_into_past",
    "80", "FALSE"
  ),
  stringsAsFactors = FALSE
)

write.table(env, out_file, sep = "\t", quote = FALSE, row.names = FALSE,
            na = "NA", fileEncoding = "UTF-8")
write.table(qa, qa_file, sep = "\t", quote = FALSE, row.names = FALSE,
            na = "NA", fileEncoding = "UTF-8")

sha <- function(path) unname(tools::md5sum(path))
report <- c(
  "# Temperature reconstruction QA",
  "",
  sprintf("- Scenario A 417-tip tree age measured from the tree: **%.15f Ma**.", total_time),
  sprintf("- Raw CENOGRID coverage: %.3f–%.3f Ma; output ends at the first observed point at/after the tree age: **%.3f Ma**.", min(time), max(time), max(env$time)),
  sprintf("- Output rows: **%d**; time is strictly increasing, unique, and finite.", nrow(env)),
  sprintf("- Historical overlap: **%d shared points**, maximum absolute difference **%.17g °C**, RMSE **%.17g °C**.", nrow(shared), max(abs(diff)), sqrt(mean(diff^2))),
  "- `t = 0` is the present and time increases into the past.",
  "- No constant endpoint, `approx(..., rule=2)`, or other extrapolation was used.",
  "- The selected CENOGRID variable is column 12, `ISOBENd18oLOESSsmoothLongTerm` (benthic δ18O, ‰ PDB). This input proxy is already long-term LOESS-smoothed.",
  "- Deep-ocean temperature was reconstructed as `8.5 - 2*d18O` where d18O > 1.75 and `12 - 4*d18O` otherwise.",
  "- Surface-air temperature was reconstructed as `2*Tdeep + 12.25` through 2.58 Ma, `2.5*Tdeep + 12.15` from >2.58 through 5.33 Ma, and `14.15 + (Tdeep - Tdeep_at_0)` before 5.33 Ma.",
  "- The discontinuities at the 2.58 and 5.33 Ma piecewise boundaries are retained because they are present in the historical local series.",
  "- Environmental models use the additional `pspline::sm.spline(..., df=80)` layer inside the archived Toussaint/RPANDA fitting function; no data-driven df selection is performed.",
  "- Provenance limitation: the original local conversion script was not retained. The equations above were recovered by an exhaustive column/formula audit and reproduce every shared historical value to floating-point precision; this script is the new explicit provenance record.",
  "",
  "## Input checksums (MD5; SHA-256 is recorded in the run manifest)",
  "",
  sprintf("- CENOGRID: `%s`", sha(raw_file)),
  sprintf("- Historical surface-temperature table: `%s`", sha(old_file)),
  sprintf("- 417-tip tree: `%s`", sha(tree_file))
)
writeLines(report, report_file, useBytes = TRUE)

cat(sprintf("Temperature reconstruction complete: %d rows, %.3f Ma coverage; max overlap error %.3g C\n",
            nrow(env), max(env$time), max(abs(diff))))
