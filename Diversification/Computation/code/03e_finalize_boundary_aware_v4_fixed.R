#!/usr/bin/env Rscript

# Checksum-pinned correction/extension of 03e_finalize_boundary_aware_v4.R.
# It uses the preferred "diagnostic numerical boundary platform" terminology
# and adds a segregated, explicitly non-inferential ten-row naive-score
# sensitivity table. The formal nine-eligible-model AICc weights are unchanged.

options(stringsAsFactors = FALSE, digits = 17, scipen = 999)

args_full <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", grep("^--file=", args_full, value = TRUE))
if (length(script_arg) != 1L) stop("Cannot determine script path")
script_file <- normalizePath(script_arg, mustWork = TRUE)
code_dir <- dirname(script_file)
base_file <- file.path(code_dir, "03e_finalize_boundary_aware_v4.R")
if (!file.exists(base_file)) stop("Missing immutable base finalizer: ", base_file)

sha256 <- function(path) {
  answer <- system2("sha256sum", path, stdout = TRUE, stderr = TRUE)
  status <- attr(answer, "status")
  if (!length(answer) || (!is.null(status) && status != 0L)) stop("sha256sum failed: ", path)
  strsplit(answer[[1L]], "[[:space:]]+")[[1L]][[1L]]
}

expected_base_sha256 <- "5f7d6bfd76ae28e3fd20629cffc2fd31879bf140823179d01be78924fe1c2c00"
if (!identical(sha256(base_file), expected_base_sha256)) {
  stop("Base finalizer SHA-256 mismatch; refusing an unaudited transformation")
}

code <- readLines(base_file, warn = FALSE, encoding = "UTF-8")

replace_fixed <- function(x, old, new, label, expected = 1L) {
  hits <- sum(grepl(old, x, fixed = TRUE))
  if (hits != expected) stop("Expected ", expected, " replacement(s) for ", label, "; observed ", hits)
  sub(old, new, x, fixed = TRUE)
}

insert_before <- function(x, anchor, addition, label) {
  at <- which(x == anchor)
  if (length(at) != 1L) stop("Expected one insertion anchor for ", label, "; observed ", length(at))
  append(x, addition, after = at - 1L)
}

insert_after <- function(x, anchor, addition, label) {
  at <- which(x == anchor)
  if (length(at) != 1L) stop("Expected one insertion anchor for ", label, "; observed ", length(at))
  append(x, addition, after = at)
}

# The evaluated extreme is a numerical boundary platform, not a demonstrated
# analytic likelihood supremum.
code <- gsub("diagnostic_supremum", "diagnostic_boundary_platform", code, fixed = TRUE)
code <- gsub("Diagnostic supremum", "Diagnostic numerical boundary platform", code, fixed = TRUE)
code <- gsub("diagnostic supremum", "diagnostic numerical boundary platform", code, fixed = TRUE)
code <- gsub("likelihood supremum", "numerical boundary platform", code, fixed = TRUE)
code <- gsub("03e_finalize_boundary_aware_v4_signature_1",
             "03e_finalize_boundary_aware_v4_signature_2", code, fixed = TRUE)

naive_score_code <- c(
  "# Segregated sensitivity analysis only: deliberately combine the nine regular",
  "# MLE log likelihoods with the largest evaluated BCSTDTempVar numerical boundary",
  "# platform. Every row is invalid for inference because this is not a regular",
  "# ten-model comparison and the boundary platform is not a finite MLE.",
  "diagnostic_naive_scores <- data.frame(",
  "  model = model_table$model, model_status = model_table$model_status,",
  "  score_source = ifelse(",
  "    model_table$eligible_for_AICc, \"finite_regular_MLE_but_mixed_sensitivity_table\",",
  "    \"approximate_numerical_boundary_platform_not_MLE\"",
  "  ),",
  "  k = model_table$k,",
  "  diagnostic_naive_logLik = ifelse(",
  "    model_table$eligible_for_AICc, model_table$logLik, diagnostic_boundary_platform",
  "  ),",
  "  stringsAsFactors = FALSE",
  ")",
  "diagnostic_naive_scores$diagnostic_naive_AIC <-",
  "  -2 * diagnostic_naive_scores$diagnostic_naive_logLik + 2 * diagnostic_naive_scores$k",
  "diagnostic_naive_scores$diagnostic_naive_AICc <-",
  "  diagnostic_naive_scores$diagnostic_naive_AIC +",
  "  2 * diagnostic_naive_scores$k * (diagnostic_naive_scores$k + 1) /",
  "  (n_obs - diagnostic_naive_scores$k - 1)",
  "diagnostic_naive_scores$diagnostic_naive_deltaAICc <-",
  "  diagnostic_naive_scores$diagnostic_naive_AICc -",
  "  min(diagnostic_naive_scores$diagnostic_naive_AICc)",
  "naive_relative <- exp(-0.5 * diagnostic_naive_scores$diagnostic_naive_deltaAICc)",
  "diagnostic_naive_scores$diagnostic_naive_Akaike_weight <- naive_relative / sum(naive_relative)",
  "diagnostic_naive_scores$valid_for_inference <- FALSE",
  "diagnostic_naive_scores$reason <- ifelse(",
  "  diagnostic_naive_scores$model == boundary_model,",
  "  paste(",
  "    \"BCSTDTempVar has no identified finite MLE; this row substitutes only the\",",
  "    \"largest evaluated approximate numerical boundary platform and its ordinary\",",
  "    \"AICc/weight are invalid.\"",
  "  ),",
  "  paste(",
  "    \"Sensitivity-only row: although this model has a regular finite MLE, the\",",
  "    \"ten-row normalization includes a non-MLE BCSTDTempVar boundary platform;\",",
  "    \"use the formal boundary-aware table instead.\"",
  "  )",
  ")",
  "diagnostic_naive_scores$boundary_platform_beta <- ifelse(",
  "  diagnostic_naive_scores$model == boundary_model,",
  "  as.numeric(diagnostic_best$native_final[\"beta\"]), NA_real_",
  ")",
  "diagnostic_naive_scores$boundary_platform_max_mu_per_Ma <- ifelse(",
  "  diagnostic_naive_scores$model == boundary_model, diagnostic_best$rate_qa$max_mu, NA_real_",
  ")",
  "diagnostic_naive_scores <- clean_frame(diagnostic_naive_scores)",
  ""
)
code <- insert_before(
  code, "# Build one clean optimization ledger: all original starts for the nine regular",
  naive_score_code, "non-inferential naive-score sensitivity table"
)

code <- insert_after(
  code, "optimization_runs$analysis_signature <- analysis_signature",
  c("diagnostic_naive_scores$analysis_signature <- analysis_signature", ""),
  "sensitivity signature"
)

qa_extension <- c(
  "qa <- rbind(qa, data.frame(",
  "  check = c(",
  "    \"diagnostic_naive_10model_sensitivity_rows\",",
  "    \"diagnostic_naive_scores_all_invalid_for_inference\",",
  "    \"diagnostic_naive_scores_carry_explicit_reason\",",
  "    \"diagnostic_naive_weights_sum_only_as_arithmetic_sensitivity\",",
  "    \"formal_ranking_excludes_boundary_platform\",",
  "    \"diagnostic_sensitivity_file_not_formal_ranking_source\"",
  "  ),",
  "  observed = c(",
  "    nrow(diagnostic_naive_scores),",
  "    all(!diagnostic_naive_scores$valid_for_inference),",
  "    all(nzchar(diagnostic_naive_scores$reason)),",
  "    format(sum(diagnostic_naive_scores$diagnostic_naive_Akaike_weight), digits = 17),",
  "    all(is.na(model_table$Akaike_weight[model_table$model == boundary_model])),",
  "    all(!diagnostic_naive_scores$valid_for_inference) &&",
  "      all(model_table$eligible_for_AICc == (model_table$model != boundary_model))",
  "  ),",
  "  expected = c(",
  "    \"10\", \"TRUE\", \"TRUE\", \"1 (arithmetic only; not inferential)\",",
  "    \"TRUE\", \"TRUE\"",
  "  ),",
  "  pass = c(",
  "    nrow(diagnostic_naive_scores) == 10L,",
  "    all(!diagnostic_naive_scores$valid_for_inference),",
  "    all(nzchar(diagnostic_naive_scores$reason)),",
  "    abs(sum(diagnostic_naive_scores$diagnostic_naive_Akaike_weight) - 1) < 1e-12,",
  "    all(is.na(model_table$Akaike_weight[model_table$model == boundary_model])) &&",
  "      abs(sum(model_table$Akaike_weight[model_table$eligible_for_AICc]) - 1) < 1e-12,",
  "    all(!diagnostic_naive_scores$valid_for_inference) &&",
  "      all(model_table$eligible_for_AICc == (model_table$model != boundary_model))",
  "  ),",
  "  stringsAsFactors = FALSE",
  "))",
  ""
)
code <- insert_before(code, "qa <- clean_frame(qa)", qa_extension, "sensitivity QA")

code <- insert_before(
  code,
  'write_clean_tsv(optimization_runs, file.path(final_model_dir, "optimization_runs.tsv"))',
  c(
    'write_clean_tsv(diagnostic_naive_scores, file.path(final_table_dir, "diagnostic_naive_10model_scores_NOT_FOR_INFERENCE.tsv"))',
    ""
  ),
  "sensitivity TSV output"
)

code <- insert_after(
  code,
  '  model_selection_QA = qa,',
  '  diagnostic_naive_10model_scores_NOT_FOR_INFERENCE = diagnostic_naive_scores,',
  "sensitivity RDS payload"
)

evidence_anchor <- '  "- `model_selection_QA.tsv`: strict expected-failure-aware QA; `ten_models_have_finite_MLE = FALSE` is an expected, passing fact.",'
code <- insert_after(
  code, evidence_anchor,
  '  "- `diagnostic_naive_10model_scores_NOT_FOR_INFERENCE.tsv`: transparent arithmetic sensitivity only; every row is marked `valid_for_inference = FALSE`, and this file is excluded from the formal ranking and figures.",',
  "sensitivity report evidence"
)

parsed <- tryCatch(
  parse(text = code, srcfile = "03e_finalize_boundary_aware_v4.R [fixed boundary-platform extension]"),
  error = function(e) stop("Transformed fixed finalizer did not parse: ", conditionMessage(e))
)
eval(parsed, envir = .GlobalEnv)
