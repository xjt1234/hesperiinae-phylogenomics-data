#!/usr/bin/env Rscript

# Final checksum-pinned wrapper for the boundary-aware v4 finalizer.
# Corrects BOM handling for R 4.5/PCRE2 and tightens terminology so the nine
# AICc-eligible fits are not broadly called "regular" when several extinction
# estimates carry an explicit mu0_near_zero boundary flag.

options(stringsAsFactors = FALSE, digits = 17, scipen = 999)

args_full <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", grep("^--file=", args_full, value = TRUE))
if (length(script_arg) != 1L) stop("Cannot determine script path")
script_file <- normalizePath(script_arg, mustWork = TRUE)
code_dir <- dirname(script_file)
fixed1_file <- file.path(code_dir, "03e_finalize_boundary_aware_v4_fixed.R")
if (!file.exists(fixed1_file)) stop("Missing checksum-pinned fixed1 wrapper")

sha256 <- function(path) {
  answer <- system2("sha256sum", path, stdout = TRUE, stderr = TRUE)
  status <- attr(answer, "status")
  if (!length(answer) || (!is.null(status) && status != 0L)) stop("sha256sum failed: ", path)
  strsplit(answer[[1L]], "[[:space:]]+")[[1L]][[1L]]
}

expected_fixed1_sha256 <- "8da5a30fa056359872b213aad3156c3f63911274a614311c4de3912d00c7a299"
if (!identical(sha256(fixed1_file), expected_fixed1_sha256)) {
  stop("fixed1 wrapper SHA-256 mismatch; refusing an unaudited transformation")
}

wrapper <- readLines(fixed1_file, warn = FALSE, encoding = "UTF-8")
parse_anchor <- which(wrapper == "parsed <- tryCatch(")
if (length(parse_anchor) != 1L) stop("Could not locate unique fixed1 parse anchor")

final_corrections <- c(
  "# R 4.5 uses PCRE2, which rejects a literal escaped \\u token in a PCRE pattern.",
  "bom_at <- which(grepl(\"y <- sub\", code, fixed = TRUE) & grepl(\"ufeff\", code, fixed = TRUE))",
  "if (length(bom_at) != 1L) stop(\"Could not locate unique BOM-cleaning line\")",
  "code[bom_at] <- '  if (length(y)) y <- sub(paste0(\"^\", intToUtf8(0xfeff)), \"\", y)'",
  "",
  "# Use finite/stable/reproducible and AICc-eligible terminology. Several of",
  "# these nine estimates retain mu0_near_zero flags and are not described",
  "# generically as regular interior birth-death estimates.",
  "code <- gsub('\"finite_MLE\"', '\"finite_stable_reproducible_estimate\"', code, fixed = TRUE)",
  "code <- gsub(\"finite regular MLEs\", \"finite, stable, reproducible estimates\", code, fixed = TRUE)",
  "code <- gsub(\"regular finite-MLE\", \"AICc-eligible finite-estimate\", code, fixed = TRUE)",
  "code <- gsub(\"regular finite MLE\", \"AICc-eligible finite estimate\", code, fixed = TRUE)",
  "code <- gsub(\"regular-MLE candidates\", \"AICc-eligible finite-estimate candidates\", code, fixed = TRUE)",
  "code <- gsub(\"ordinary regular-MLE assumptions\", \"ordinary finite-interior-MLE assumptions\", code, fixed = TRUE)",
  "code <- gsub(\"nine regular\", \"nine AICc-eligible\", code, fixed = TRUE)",
  "code <- gsub(\"eligible finite MLE\", \"AICc-eligible finite, stable, reproducible estimate\", code, fixed = TRUE)",
  "code <- gsub(\"finite_regular_MLE_models\", \"finite_stable_reproducible_estimate_models\", code, fixed = TRUE)",
  "code <- gsub(\"finite_regular_MLE_but_mixed_sensitivity_table\", \"finite_stable_reproducible_estimate_but_mixed_sensitivity_table\", code, fixed = TRUE)",
  "",
  "near_zero_code <- c(",
  "  'near_zero_models <- model_table$model[model_table$eligible_for_AICc &',",
  "  '  !is.na(model_table$boundary_flag) & grepl(\"mu0_near_zero\", model_table$boundary_flag, fixed = TRUE)]',",
  "  'near_zero_statement <- if (length(near_zero_models)) {',",
  "  '  paste0(\"Extinction-rate caution: \" , length(near_zero_models),',",
  "  '         \" AICc-eligible fits retain `mu0_near_zero` boundary flags: `\",',",
  "  '         paste(near_zero_models, collapse = \"`, `\"),',",
  "  '         \"`. Their extinction components are weakly resolved and should not be interpreted as robust interior estimates.\")',",
  "  '} else {',",
  "  '  \"No AICc-eligible fit carries a near-zero extinction boundary flag.\"',",
  "  '}',",
  "  ''",
  ")",
  "code <- insert_before(code, 'report <- c(', near_zero_code, 'near-zero extinction statement')",
  "code <- insert_before(",
  "  code, '  \"## Evidence files\",',",
  "  c('  \"## Extinction-rate boundary caution\",', '  \"\",',",
  "    '  near_zero_statement,', '  \"\",'),",
  "  'near-zero extinction report section'",
  ")",
  ""
)

wrapper <- append(wrapper, final_corrections, after = parse_anchor - 1L)
parsed_wrapper <- tryCatch(
  parse(text = wrapper, srcfile = "03e_finalize_boundary_aware_v4_fixed.R [PCRE2/terminology correction]"),
  error = function(e) stop("Corrected finalizer wrapper did not parse: ", conditionMessage(e))
)
eval(parsed_wrapper, envir = .GlobalEnv)
