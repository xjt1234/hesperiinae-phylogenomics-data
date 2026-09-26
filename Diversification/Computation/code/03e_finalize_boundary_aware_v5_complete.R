#!/usr/bin/env Rscript

# Checksum-pinned, non-overwriting v5 finalizer.
#
# This wrapper retains the formal boundary-aware inference from v4, writes to
# new final_v5_boundary_aware_complete directories, and closes the provenance
# gap in the optimization ledger by retaining every original BCSTDTempVar run
# that triggered the subsequent boundary analysis. It also uses deliberately
# neutral terminology for the nine AICc-eligible finite, stable, reproducible
# estimates because several carry explicit mu0_near_zero flags.

options(stringsAsFactors = FALSE, digits = 17, scipen = 999)

args_full <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", grep("^--file=", args_full, value = TRUE))
if (length(script_arg) != 1L) stop("Cannot determine script path")
script_file <- normalizePath(script_arg, mustWork = TRUE)
code_dir <- dirname(script_file)
fixed1_file <- file.path(code_dir, "03e_finalize_boundary_aware_v4_fixed.R")
fixed2_file <- file.path(code_dir, "03e_finalize_boundary_aware_v4_fixed2.R")
if (!file.exists(fixed1_file) || !file.exists(fixed2_file)) {
  stop("Missing checksum-pinned v4 wrapper(s)")
}

sha256 <- function(path) {
  answer <- system2("sha256sum", path, stdout = TRUE, stderr = TRUE)
  status <- attr(answer, "status")
  if (!length(answer) || (!is.null(status) && status != 0L)) {
    stop("sha256sum failed: ", path)
  }
  strsplit(answer[[1L]], "[[:space:]]+")[[1L]][[1L]]
}

expected_fixed1_sha256 <- "8da5a30fa056359872b213aad3156c3f63911274a614311c4de3912d00c7a299"
expected_fixed2_sha256 <- "9fd9b4abfef81fbe9864b7bed9140312a06d505d10ce62c196b11791dbdc7b76"
if (!identical(sha256(fixed1_file), expected_fixed1_sha256) ||
    !identical(sha256(fixed2_file), expected_fixed2_sha256)) {
  stop("v4 wrapper SHA-256 mismatch; refusing an unaudited transformation")
}

# Begin with the fixed1 transformer and inject both the already-audited fixed2
# corrections and the v5 complete-ledger changes before it parses the base.
wrapper <- readLines(fixed1_file, warn = FALSE, encoding = "UTF-8")
parse_anchor <- which(wrapper == "parsed <- tryCatch(")
if (length(parse_anchor) != 1L) stop("Could not locate unique fixed1 parse anchor")

final_corrections <- c(
  "# R 4.5 uses PCRE2, which rejects a literal escaped \\u token in a PCRE pattern.",
  "bom_at <- which(grepl(\"y <- sub\", code, fixed = TRUE) & grepl(\"ufeff\", code, fixed = TRUE))",
  "if (length(bom_at) != 1L) stop(\"Could not locate unique BOM-cleaning line\")",
  "code[bom_at] <- '  if (length(y)) y <- sub(paste0(\"^\", intToUtf8(0xfeff)), \"\", y)'",
  "",
  "# Use finite/stable/reproducible and AICc-eligible terminology.",
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
  "code <- gsub(\"regular_finite_MLE_models\", \"aicc_eligible_finite_estimate_models\", code, fixed = TRUE)",
  "code <- gsub(\"regular_AICc\", \"eligible_AICc\", code, fixed = TRUE)",
  "code <- gsub(\"regular_selected\", \"eligible_selected\", code, fixed = TRUE)",
  "code <- gsub(\"reproducible finite, stable, reproducible estimates\", \"finite, stable, reproducible estimates\", code, fixed = TRUE)",
  "code <- gsub(\"nine eligible AICc-eligible finite-estimate candidates\", \"nine AICc-eligible models with finite, stable, reproducible estimates\", code, fixed = TRUE)",
  "code <- gsub(\"a AICc-eligible finite estimate\", \"an AICc-eligible finite, stable, reproducible estimate\", code, fixed = TRUE)",
  "",
  "# Write a new immutable final product; retain v4 as an audit intermediate.",
  "code <- gsub(\"final_v4_boundary_aware\", \"final_v5_boundary_aware_complete\", code, fixed = TRUE)",
  "code <- gsub(\"analysis_signature_v4_boundary_aware_final\", \"analysis_signature_v5_boundary_aware_complete_final\", code, fixed = TRUE)",
  "code <- gsub(\"03e_finalize_boundary_aware_v4_signature_2\", \"03e_finalize_boundary_aware_v5_complete_signature_1\", code, fixed = TRUE)",
  "code <- gsub(\"final-v4\", \"final-v5\", code, fixed = TRUE)",
  "code <- gsub(\"Boundary-aware v4 finalization\", \"Boundary-aware v5 complete finalization\", code, fixed = TRUE)",
  "code <- gsub(\"9 finite-MLE eligible\", \"9 AICc-eligible finite/stable/reproducible estimates\", code, fixed = TRUE)",
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
  "",
  "# Retain whether a row was selected by its source analysis separately from",
  "# whether it is selected for the formal boundary-aware comparison.",
  "code <- insert_before(",
  "  code, '    selected_for_formal_MLE = selected, diagnostic_only = FALSE,',",
  "  '    selected_in_source_analysis = selected,', 'source-selection ledger field'",
  ")",
  "code <- insert_before(",
  "  code, '    selected_for_formal_MLE = FALSE, diagnostic_only = TRUE,',",
  "  '    selected_in_source_analysis = FALSE,', 'profile source-selection ledger field'",
  ")",
  "",
  "boundary_ledger_code <- c(",
  "  '# Preserve every original BCSTDTempVar run that triggered boundary follow-up.',",
  "  'boundary_round1_result <- round1$all_results[[boundary_model]]',",
  "  'for (i in seq_along(boundary_round1_result$runs)) {',",
  "  '  x <- boundary_round1_result$runs[[i]]',",
  "  '  row <- regular_run_row(x, boundary_round1_result$spec, i == boundary_round1_result$best_i)',",
  "  '  row$model_status <- \"boundary_nonidentifiable\"',",
  "  '  row$eligible_for_AICc <- FALSE',",
  "  '  row$run_source <- \"first_round_v2_boundary_trigger\"',",
  "  '  row$diagnostic_naive_AICc <- if (is.null(x$fit)) NA_real_ else x$fit$aicc',",
  "  '  row$selected_for_formal_MLE <- FALSE',",
  "  '  row$diagnostic_only <- TRUE',",
  "  '  row$diagnostic_sequence_best <- FALSE',",
  "  '  run_rows[[length(run_rows) + 1L]] <- row',",
  "  '}',",
  "  ''",
  ")",
  "code <- insert_before(",
  "  code, 'for (run in full_runs) run_rows[[length(run_rows) + 1L]] <- profile_run_row(run)',",
  "  boundary_ledger_code, 'complete first-round boundary ledger'",
  ")",
  "",
  "ledger_audit_code <- c(",
  "  'profile_runs_table <- read_tsv_strict(profile_runs_file)',",
  "  'fixed_beta_profile_reproduction_min <- min(profile_summary$best_start_count)',",
  "  'fixed_beta_no_higher_finite_nonconverged <- all(vapply(seq_len(nrow(profile_summary)), function(i) {',",
  "  '  candidate <- profile_runs_table[',",
  "  '    abs(profile_runs_table$fixed_beta - profile_summary$beta[i]) <= 1e-12 &',",
  "  '      is.finite(profile_runs_table$logLik) & profile_runs_table$convergence != 0L,',",
  "  '    , drop = FALSE',",
  "  '  ]',",
  "  '  !nrow(candidate) || all(candidate$logLik <= profile_summary$logLik[i] + 1e-8)',",
  "  '}, logical(1)))',",
  "  'expected_round1_ledger_rows <- sum(vapply(round1$all_results, function(x) length(x$runs), integer(1)))',",
  "  'expected_complete_ledger_rows <- expected_round1_ledger_rows + length(full_runs)',",
  "  'boundary_round1_ledger <- optimization_runs$model == boundary_model &',",
  "  '  optimization_runs$run_source == \"first_round_v2_boundary_trigger\"',",
  "  'boundary_profile_ledger <- optimization_runs$model == boundary_model &',",
  "  '  optimization_runs$run_source == \"03c_free_logcenter_boundary_sequence\"',",
  "  ''",
  ")",
  "code <- insert_after(",
  "  code, 'optimization_runs <- clean_frame(do.call(rbind, run_rows))',",
  "  ledger_audit_code, 'complete-ledger audit values'",
  ")",
  "",
  "signature_ledger_code <- c(",
  "  'signature_components <- rbind(signature_components, data.frame(',",
  "  '  component = c(\"optimization_ledger_schema\", \"optimization_ledger_rows\",',",
  "  '                \"round1_boundary_trigger_runs\", \"profile_free_boundary_runs\"),',",
  "  '  value = c(\"complete_round1_plus_03c_free_logcenter\",',",
  "  '            as.character(nrow(optimization_runs)),',",
  "  '            as.character(sum(boundary_round1_ledger)),',",
  "  '            as.character(sum(boundary_profile_ledger))),',",
  "  '  stringsAsFactors = FALSE',",
  "  '))',",
  "  ''",
  ")",
  "code <- insert_before(",
  "  code, 'signature_components <- clean_frame(signature_components)',",
  "  signature_ledger_code, 'ledger signature components'",
  ")",
  "",
  "qa_ledger_code <- c(",
  "  'qa <- rbind(qa, data.frame(',",
  "  '  check = c(',",
  "  '    \"complete_optimization_ledger_rows\",',",
  "  '    \"all_round1_runs_retained_in_ledger\",',",
  "  '    \"BCSTDTempVar_round1_boundary_trigger_runs_retained\",',",
  "  '    \"BCSTDTempVar_round1_runs_diagnostic_only_ineligible_unselected\",',",
  "  '    \"BCSTDTempVar_source_selection_preserved_but_not_formal\",',",
  "  '    \"BCSTDTempVar_free_logcenter_runs_retained\",',",
  "  '    \"fixed_beta_each_point_reproduced_twice\",',",
  "  '    \"fixed_beta_no_higher_finite_nonconverged_run\"',",
  "  '  ),',",
  "  '  observed = c(',",
  "  '    nrow(optimization_runs),',",
  "  '    sum(optimization_runs$run_source %in% c(\"first_round_v2\", \"first_round_v2_boundary_trigger\")),',",
  "  '    sum(boundary_round1_ledger),',",
  "  '    all(optimization_runs$diagnostic_only[boundary_round1_ledger] &',",
  "  '        !optimization_runs$eligible_for_AICc[boundary_round1_ledger] &',",
  "  '        !optimization_runs$selected_for_formal_MLE[boundary_round1_ledger]),',",
  "  '    paste(sum(optimization_runs$selected_in_source_analysis[boundary_round1_ledger]),',",
  "  '          sum(optimization_runs$selected_for_formal_MLE[boundary_round1_ledger]), sep = \"/\"),',",
  "  '    sum(boundary_profile_ledger),',",
  "  '    fixed_beta_profile_reproduction_min,',",
  "  '    fixed_beta_no_higher_finite_nonconverged',",
  "  '  ),',",
  "  '  expected = c(',",
  "  '    as.character(expected_complete_ledger_rows),',",
  "  '    as.character(expected_round1_ledger_rows),',",
  "  '    as.character(length(boundary_round1_result$runs)),',",
  "  '    \"TRUE\", \"1 source-selected / 0 formal-selected\",',",
  "  '    as.character(length(full_runs)), \">=2\", \"TRUE\"',",
  "  '  ),',",
  "  '  pass = c(',",
  "  '    nrow(optimization_runs) == expected_complete_ledger_rows,',",
  "  '    sum(optimization_runs$run_source %in% c(\"first_round_v2\", \"first_round_v2_boundary_trigger\")) ==',",
  "  '      expected_round1_ledger_rows,',",
  "  '    sum(boundary_round1_ledger) == length(boundary_round1_result$runs),',",
  "  '    all(optimization_runs$diagnostic_only[boundary_round1_ledger] &',",
  "  '        !optimization_runs$eligible_for_AICc[boundary_round1_ledger] &',",
  "  '        !optimization_runs$selected_for_formal_MLE[boundary_round1_ledger]),',",
  "  '    sum(optimization_runs$selected_in_source_analysis[boundary_round1_ledger]) == 1L &&',",
  "  '      sum(optimization_runs$selected_for_formal_MLE[boundary_round1_ledger]) == 0L,',",
  "  '    sum(boundary_profile_ledger) == length(full_runs),',",
  "  '    fixed_beta_profile_reproduction_min >= 2L,',",
  "  '    fixed_beta_no_higher_finite_nonconverged',",
  "  '  ),',",
  "  '  stringsAsFactors = FALSE',",
  "  '))',",
  "  ''",
  ")",
  "code <- insert_before(code, 'qa <- clean_frame(qa)', qa_ledger_code, 'complete-ledger QA')",
  "",
  "code <- insert_after(",
  "  code, 'all_results <- all_results[model_names]',",
  "  c('all_results[[boundary_model]]$round1_boundary_trigger <-',",
  "    '  round1$all_results[[boundary_model]]'),",
  "  'round-one boundary evidence in final RDS'",
  ")",
  "",
  "code <- gsub(",
  "  \"nine-model original multi-start runs plus ten free log-centred boundary-diagnostic runs\",",
  "  \"all first-round multi-start runs (including the BCSTDTempVar boundary-trigger runs) plus ten free log-centred boundary-diagnostic runs\",",
  "  code, fixed = TRUE",
  ")",
  ""
)

wrapper <- append(wrapper, final_corrections, after = parse_anchor - 1L)
parsed_wrapper <- tryCatch(
  parse(text = wrapper, srcfile = "03e_finalize_boundary_aware_v4_fixed.R [v5 complete-ledger correction]"),
  error = function(e) stop("Corrected v5 finalizer wrapper did not parse: ", conditionMessage(e))
)
eval(parsed_wrapper, envir = .GlobalEnv)
