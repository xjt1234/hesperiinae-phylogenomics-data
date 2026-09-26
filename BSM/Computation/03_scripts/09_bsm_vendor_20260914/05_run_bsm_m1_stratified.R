#!/usr/bin/env Rscript

# Reproducible, fail-closed BSM runner for the accepted Scenario-A / M1 DEC fit.
#
# This script NEVER searches for a result object.  The fit, postfit, independent
# model-scope FINAL_QA.json, and its validator-owned acceptance record must be
# supplied explicitly.  It defaults to a read-only preflight; stochastic maps
# are generated only with --action run.
#
# Acceptance TSV schema (exactly one data row):
#   status  model  fit_sha256
# where status=ACCEPTED, model=M1, and fit_sha256 is the SHA-256 of --fit-rds.
# The independent post-fit validator, not this script, is responsible for
# creating that acceptance record and the sibling FINAL_QA.json.  Both must be
# in 08_qa/<model-validation-dir>/ and identify M1_R4 PASS/COMPLETE.
#
# First run (100 successful maps):
#   bash 03_scripts/r44_env.sh \
#     --file=03_scripts/05_run_bsm_m1_stratified.R --args \
#     --action run --fit-rds /ABS/PATH/fit_result.rds \
#     --postfit-rds /ABS/PATH/postfit/postfit_recalculated_ancestral_states.rds \
#     --final-qa-json /ABS/PATH/08_qa/M1_SCOPE/FINAL_QA.json \
#     --acceptance-tsv /ABS/PATH/08_qa/M1_SCOPE/M1_R4_final_acceptance.tsv \
#     --outdir /ABS/PATH/05_bsm/M1_BSM_<timestamp> --target-maps 100
#
# Deterministic extension (same directory and seed schedule):
#   bash 03_scripts/r44_env.sh \
#     --file=03_scripts/05_run_bsm_m1_stratified.R --args \
#     --action run ... --outdir SAME --target-maps 500 --resume true
#
# target-maps must be 100, 200, 300, 400, or 500.  A checkpoint summary is
# made after each complete block of 100 maps.  With --stop-when-stable true,
# an extension stops at the first >=200 checkpoint passing the pre-specified
# precision, cumulative-change, and spatial-ranking criteria below.

suppressPackageStartupMessages({
  library(ape)
  library(BioGeoBEARS)
  library(methods)
  library(digest)
  library(jsonlite)
})

parse_args <- function(x) {
  out <- list()
  i <- 1L
  while (i <= length(x)) {
    if (!startsWith(x[[i]], "--") || i == length(x)) {
      stop("Arguments must be --key value pairs")
    }
    out[[sub("^--", "", x[[i]])]] <- x[[i + 1L]]
    i <- i + 2L
  }
  out
}

as_bool <- function(x, name) {
  y <- tolower(x)
  if (!(y %in% c("true", "false"))) stop(name, " must be true or false")
  identical(y, "true")
}

required <- function(opts, keys) {
  absent <- setdiff(keys, names(opts))
  if (length(absent) > 0L) stop("Missing arguments: ", paste(absent, collapse=", "))
}

opts <- parse_args(commandArgs(trailingOnly=TRUE))
required(opts, c("action", "fit-rds", "postfit-rds", "final-qa-json",
                 "acceptance-tsv", "outdir"))

action <- tolower(opts$action)
if (!(action %in% c("preflight", "run", "summarize"))) {
  stop("action must be preflight, run, or summarize")
}
target_maps <- as.integer(if (is.null(opts[["target-maps"]])) "100" else opts[["target-maps"]])
if (!(target_maps %in% seq(100L, 500L, 100L))) {
  stop("target-maps must be one of 100, 200, 300, 400, 500")
}
maxtries <- as.integer(if (is.null(opts$maxtries)) "40000" else opts$maxtries)
if (!is.finite(maxtries) || maxtries < 1L) stop("maxtries must be a positive integer")
seed_base <- as.numeric(if (is.null(opts[["seed-base"]])) "202609050" else opts[["seed-base"]])
if (!is.finite(seed_base) || seed_base < 1 || seed_base > 2e9) stop("seed-base must be in [1, 2e9]")
max_attempts <- as.integer(if (is.null(opts[["max-attempts"]])) as.character(2L * target_maps) else opts[["max-attempts"]])
if (!is.finite(max_attempts) || max_attempts < target_maps) stop("max-attempts must be >= target-maps")
resume <- as_bool(if (is.null(opts$resume)) "false" else opts$resume, "resume")
stop_when_stable <- as_bool(if (is.null(opts[["stop-when-stable"]])) "true" else opts[["stop-when-stable"]], "stop-when-stable")

script_arg <- grep("^--file=", commandArgs(FALSE), value=TRUE)
if (length(script_arg) != 1L) stop("Cannot resolve script path")
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork=TRUE)
job_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork=TRUE)
postfit_rds <- normalizePath(opts[["postfit-rds"]], mustWork=TRUE)
fit_rds <- normalizePath(opts[["fit-rds"]], mustWork=TRUE)
final_qa_json <- normalizePath(opts[["final-qa-json"]], mustWork=TRUE)
acceptance_tsv <- normalizePath(opts[["acceptance-tsv"]], mustWork=TRUE)
outdir <- normalizePath(opts$outdir, mustWork=FALSE)

inside <- function(path, root) startsWith(path, paste0(root, "/"))
if (!inside(fit_rds, job_root)) stop("fit-rds must be inside this rerun job")
if (!inside(postfit_rds, job_root)) stop("postfit-rds must be inside this rerun job")
if (!inside(final_qa_json, job_root)) stop("final-qa-json must be inside this rerun job")
if (!inside(acceptance_tsv, job_root)) stop("acceptance-tsv must be inside this rerun job")
expected_postfit_parent <- normalizePath(file.path(dirname(fit_rds), "postfit"), mustWork=TRUE)
if (!identical(dirname(postfit_rds), expected_postfit_parent) ||
    !identical(basename(postfit_rds), "postfit_recalculated_ancestral_states.rds")) {
  stop("postfit-rds must be the named postfit object directly under the accepted fit directory")
}
qa_root <- normalizePath(file.path(job_root, "08_qa"), mustWork=TRUE)
qa_scope_dir <- dirname(final_qa_json)
if (!identical(basename(final_qa_json), "FINAL_QA.json") ||
    !identical(dirname(qa_scope_dir), qa_root)) {
  stop("final-qa-json must be canonical 08_qa/<model-validation-dir>/FINAL_QA.json")
}
if (!identical(dirname(acceptance_tsv), qa_scope_dir) ||
    !identical(basename(acceptance_tsv), "M1_R4_final_acceptance.tsv")) {
  stop("acceptance-tsv must be the sibling validator-owned M1_R4_final_acceptance.tsv")
}
allowed_bsm_root <- normalizePath(file.path(job_root, "05_bsm"), mustWork=TRUE)
if (!inside(outdir, allowed_bsm_root)) stop("outdir must be a child of ", allowed_bsm_root)

timestamp <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
stamp_id <- function() format(Sys.time(), "%Y%m%dT%H%M%S%z")
sha256_file <- function(path) digest(file=path, algo="sha256", serialize=FALSE)
sha256_object <- function(x) digest(x, algo="sha256", serialize=TRUE)
clean_text <- function(x) gsub("[\t\r\n]+", " ", as.character(x))

require_json_keys <- function(x, keys, label) {
  if (!is.list(x)) stop(label, " must be a JSON object")
  missing <- setdiff(keys, names(x))
  if (length(missing) > 0L) stop(label, " lacks fields: ", paste(missing, collapse=", "))
}

json_scalar_character <- function(x, label) {
  if (length(x) != 1L || is.null(x) || is.na(x)) stop(label, " must be one JSON scalar")
  as.character(x)
}

require_sha256 <- function(x, label) {
  value <- tolower(json_scalar_character(x, label))
  if (!grepl("^[0-9a-f]{64}$", value)) stop(label, " is not a SHA-256")
  value
}

resolve_qa_artifact <- function(report, role, expected_path=NULL) {
  artifacts <- report$artifacts
  if (!is.list(artifacts) || !is.list(artifacts[[role]])) {
    stop("FINAL_QA artifacts lacks ", role)
  }
  entry <- artifacts[[role]]
  require_json_keys(entry, c("path", "sha256"), paste0("FINAL_QA ", role))
  path <- normalizePath(json_scalar_character(entry$path, paste0(role, " path")), mustWork=TRUE)
  if (!inside(path, job_root)) stop("FINAL_QA ", role, " artifact must be inside this rerun job")
  digest_expected <- require_sha256(entry$sha256, paste0(role, " sha256"))
  if (!identical(sha256_file(path), digest_expected)) {
    stop("FINAL_QA ", role, " artifact hash is stale")
  }
  if (!is.null(expected_path) && !identical(path, normalizePath(expected_path, mustWork=TRUE))) {
    stop("FINAL_QA ", role, " is not the current required implementation")
  }
  list(path=path, sha256=digest_expected)
}

atomic_save_rds_new <- function(x, path) {
  if (file.exists(path)) stop("Refusing to overwrite: ", path)
  tmp <- paste0(path, ".tmp.", Sys.getpid())
  saveRDS(x, tmp, version=3)
  if (!file.rename(tmp, path)) stop("Atomic rename failed for ", path)
  invisible(path)
}

atomic_write_tsv_new <- function(x, path, row.names=FALSE) {
  if (file.exists(path)) stop("Refusing to overwrite: ", path)
  tmp <- paste0(path, ".tmp.", Sys.getpid())
  write.table(x, tmp, sep="\t", quote=FALSE, row.names=row.names,
              col.names=if (row.names) NA else TRUE, na="NA")
  if (!file.rename(tmp, path)) stop("Atomic rename failed for ", path)
  invisible(path)
}

atomic_write_json_new <- function(x, path) {
  if (file.exists(path)) stop("Refusing to overwrite: ", path)
  tmp <- paste0(path, ".tmp.", Sys.getpid())
  jsonlite::write_json(x, tmp, auto_unbox=TRUE, pretty=TRUE, digits=NA,
                       null="null", na="null")
  if (!file.rename(tmp, path)) stop("Atomic rename failed for ", path)
  invisible(path)
}

write_tsv <- function(x, path) {
  write.table(x, path, sep="\t", quote=FALSE, row.names=FALSE,
              col.names=TRUE, na="NA")
}

bgb_fun <- function(name) {
  ns <- asNamespace("BioGeoBEARS")
  if (!exists(name, envir=ns, inherits=FALSE)) stop("Missing BioGeoBEARS function: ", name)
  get(name, envir=ns, inherits=FALSE)
}

expected_paths <- list(
  tree=file.path(job_root, "01_inputs", "frozen", "tree_scenarioA417.tre"),
  geog=file.path(job_root, "01_inputs", "frozen", "geog_scenarioA417_analysis_order.LagrangePHYLIP"),
  times=file.path(job_root, "02_config", "timeperiods_5epochs.txt"),
  multipliers=file.path(job_root, "02_config", "M1_conservative_dispersal_multipliers.txt"),
  area_order=file.path(job_root, "02_config", "area_order.tsv")
)
expected_hashes <- c(
  tree="06725d5a0c1ef29007aeb0f7657175c2d2329c75cb4d8fb2170e39b5fdaf6962",
  geog="89c5b776bdce88d9c0074b3a6117d25e50c8604a69dba158d39fa5544527f1c0",
  times="90dadc7755385b44d9333f22e2acaeac163bd0ce4507289e1720671eccc4e51f",
  multipliers="34cc3ee7f0789d978948e35085fb9ed89de09f2600b50384bb074bd11db7029c",
  area_order="7daa151006fbe79d00ed92a83938aca4b04dbc89bbaa814710bf244dcb1a0804"
)
for (nm in names(expected_paths)) {
  p <- normalizePath(expected_paths[[nm]], mustWork=TRUE)
  if (!identical(sha256_file(p), unname(expected_hashes[[nm]]))) {
    stop("Frozen ", nm, " hash mismatch")
  }
  expected_paths[[nm]] <- p
}

area_order <- read.delim(expected_paths$area_order, stringsAsFactors=FALSE,
                         check.names=FALSE)
needed_area_cols <- c("analysis_bit_index", "analysis_internal_code", "abbrev", "full_name")
if (!all(needed_area_cols %in% names(area_order))) stop("area_order.tsv schema mismatch")
area_order <- area_order[order(area_order$analysis_bit_index), ]
internal_codes <- LETTERS[1:11]
actual_names <- c("AF", "AUS", "CAM", "ENA", "EPA", "IND", "MDG", "ORI", "SAM", "WNA", "WPA")
if (!identical(area_order$analysis_internal_code, internal_codes) ||
    !identical(area_order$abbrev, actual_names)) stop("Area-order semantic gate failed")
code_to_area <- setNames(area_order$abbrev, area_order$analysis_internal_code)

fit_hash <- sha256_file(fit_rds)
postfit_hash <- sha256_file(postfit_rds)
final_qa_hash <- sha256_file(final_qa_json)

final_qa <- jsonlite::read_json(final_qa_json, simplifyVector=FALSE)
require_json_keys(
  final_qa,
  c("schema_version", "status", "scope", "model_id", "job_root", "checks", "models", "artifacts"),
  "FINAL_QA"
)
qa_job_root <- normalizePath(
  json_scalar_character(final_qa$job_root, "FINAL_QA job_root"), mustWork=TRUE
)
if (!identical(json_scalar_character(final_qa$schema_version, "FINAL_QA schema_version"), "1.0") ||
    !identical(json_scalar_character(final_qa$status, "FINAL_QA status"), "PASS") ||
    !identical(json_scalar_character(final_qa$scope, "FINAL_QA scope"), "model") ||
    !identical(json_scalar_character(final_qa$model_id, "FINAL_QA model_id"), "M1_R4") ||
    !identical(qa_job_root, job_root)) {
  stop("FINAL_QA identity/status/scope/job mismatch; require model-scope M1_R4 PASS")
}
if (!is.list(final_qa$checks) || length(final_qa$checks) == 0L ||
    any(!vapply(final_qa$checks, is.list, logical(1))) ||
    any(!vapply(final_qa$checks, function(x) {
      length(x$status) == 1L && identical(as.character(x$status), "PASS")
    }, logical(1)))) {
  stop("FINAL_QA checks must be a non-empty all-PASS list")
}
has_check <- function(check_id) {
  sum(vapply(final_qa$checks, function(x) {
    length(x$check_id) == 1L && identical(as.character(x$check_id), check_id) &&
      length(x$status) == 1L && identical(as.character(x$status), "PASS")
  }, logical(1))) == 1L
}
if (!has_check("model.M1_R4.postfit") || !has_check("acceptance.emission")) {
  stop("FINAL_QA lacks unique PASS model.M1_R4.postfit and acceptance.emission checks")
}
if (!is.list(final_qa$models) || !identical(names(final_qa$models), "M1_R4")) {
  stop("FINAL_QA must contain exactly its own M1_R4 model result")
}
qa_result <- final_qa$models$M1_R4
require_json_keys(qa_result, c("status", "outcome", "inspection"), "FINAL_QA M1_R4 result")
if (!identical(json_scalar_character(qa_result$status, "M1_R4 status"), "PASS") ||
    !identical(json_scalar_character(qa_result$outcome, "M1_R4 outcome"), "COMPLETE")) {
  stop("FINAL_QA M1_R4 result is not PASS/COMPLETE")
}
qa_inspection <- qa_result$inspection
require_json_keys(qa_inspection, c("status", "model", "max_range", "fit", "postfit"),
                  "FINAL_QA M1_R4 inspection")
if (!identical(json_scalar_character(qa_inspection$status, "inspection status"), "PASS") ||
    !identical(json_scalar_character(qa_inspection$model, "inspection model"), "M1") ||
    length(qa_inspection$max_range) != 1L ||
    !isTRUE(as.numeric(qa_inspection$max_range) == 4)) {
  stop("FINAL_QA inspection must be PASS for M1 max-range 4")
}
require_json_keys(qa_inspection$fit, c("path", "sha256", "j"), "FINAL_QA inspected fit")
require_json_keys(qa_inspection$postfit, c("path", "sha256"), "FINAL_QA inspected postfit")
qa_fit_path <- normalizePath(json_scalar_character(qa_inspection$fit$path, "inspected fit path"), mustWork=TRUE)
qa_postfit_path <- normalizePath(json_scalar_character(qa_inspection$postfit$path, "inspected postfit path"), mustWork=TRUE)
qa_fit_hash <- require_sha256(qa_inspection$fit$sha256, "inspected fit sha256")
qa_postfit_hash <- require_sha256(qa_inspection$postfit$sha256, "inspected postfit sha256")
if (!identical(qa_fit_path, fit_rds) || !identical(qa_postfit_path, postfit_rds) ||
    !identical(qa_fit_hash, fit_hash) || !identical(qa_postfit_hash, postfit_hash) ||
    length(qa_inspection$fit$j) != 1L || !isTRUE(as.numeric(qa_inspection$fit$j) == 0)) {
  stop("Current fit/postfit paths, hashes, or j=0 differ from FINAL_QA inspection")
}

qa_bound_artifacts <- list(
  validator=resolve_qa_artifact(
    final_qa, "validator", file.path(job_root, "03_scripts", "08_validate_postfit_and_final.py")
  ),
  rds_inspector=resolve_qa_artifact(
    final_qa, "rds_inspector", file.path(job_root, "03_scripts", "08_inspect_postfit_rds.R")
  ),
  spec=resolve_qa_artifact(
    final_qa, "spec", file.path(job_root, "03_scripts", "final_qa_spec_v1.json")
  ),
  manifest=resolve_qa_artifact(final_qa, "manifest")
)

acceptance <- read.delim(acceptance_tsv, stringsAsFactors=FALSE, check.names=FALSE,
                         quote="", comment.char="")
needed_acceptance <- c("status", "model", "fit_sha256")
if (nrow(acceptance) != 1L || !identical(names(acceptance), needed_acceptance)) {
  stop("Acceptance TSV must have exactly one row and exact columns status, model, fit_sha256")
}
if (!identical(acceptance$status[[1]], "ACCEPTED") ||
    !identical(acceptance$model[[1]], "M1") ||
    !identical(acceptance$fit_sha256[[1]], fit_hash)) {
  stop("Independent acceptance content does not match the current accepted M1 fit")
}
qa_acceptances <- final_qa$artifacts$acceptances
if (!is.list(qa_acceptances) ||
    !identical(names(qa_acceptances), "M1_R4_final_acceptance.tsv") ||
    !is.list(qa_acceptances[["M1_R4_final_acceptance.tsv"]])) {
  stop("FINAL_QA must own exactly the M1_R4 acceptance artifact")
}
qa_acceptance <- qa_acceptances[["M1_R4_final_acceptance.tsv"]]
require_json_keys(qa_acceptance, c("path", "sha256"), "FINAL_QA acceptance artifact")
qa_acceptance_path <- normalizePath(
  json_scalar_character(qa_acceptance$path, "FINAL_QA acceptance path"), mustWork=TRUE
)
qa_acceptance_hash <- require_sha256(qa_acceptance$sha256, "FINAL_QA acceptance sha256")
if (!identical(qa_acceptance_path, acceptance_tsv) ||
    !identical(qa_acceptance_hash, sha256_file(acceptance_tsv))) {
  stop("Acceptance path/hash differs from its validator-owned FINAL_QA record")
}

if (!identical(as.character(packageVersion("BioGeoBEARS")), "1.1.3")) {
  stop("This runner was audited only for BioGeoBEARS 1.1.3")
}
remote_sha <- unname(packageDescription("BioGeoBEARS")$RemoteSha)
expected_remote_sha <- "1672cc0c171b1a05effad69fa426427b3d9ef4e1"
if (length(remote_sha) != 1L || !identical(remote_sha, expected_remote_sha)) {
  stop("BioGeoBEARS RemoteSha mismatch")
}

fit_optimizer <- readRDS(fit_rds)
if (!is.list(fit_optimizer) || is.null(fit_optimizer$inputs) || is.null(fit_optimizer$outputs)) {
  stop("Invalid optimizer result structure")
}
if (!is.finite(as.numeric(fit_optimizer$total_loglikelihood))) stop("Optimizer fit has non-finite total_loglikelihood")
if (!isS4(fit_optimizer$outputs) || !("params_table" %in% slotNames(fit_optimizer$outputs))) {
  stop("Fit lacks outputs@params_table")
}
pt_fit <- fit_optimizer$outputs@params_table
if (!all(c("d", "e", "j", "w") %in% rownames(pt_fit))) stop("Fit parameter table is incomplete")
if (!identical(rownames(pt_fit)[pt_fit[, "type"] == "free"], c("d", "e"))) stop("M1 must have only d and e free")
if (!identical(as.character(pt_fit["j", "type"]), "fixed") ||
    !isTRUE(as.numeric(pt_fit["j", "init"]) == 0) || !isTRUE(as.numeric(pt_fit["j", "est"]) == 0)) {
  stop("j must be exactly fixed at zero")
}
if (!identical(as.character(pt_fit["w", "type"]), "fixed") ||
    !isTRUE(as.numeric(pt_fit["w", "init"]) == 1) || !isTRUE(as.numeric(pt_fit["w", "est"]) == 1)) {
  stop("w must be exactly fixed at one")
}
prov <- fit_optimizer$inputs$provenance
if (is.null(prov) || !identical(prov$model, "M1") || !isTRUE(prov$final_fit) ||
    !identical(as.integer(prov$max_range_size), 4L)) {
  stop("Result provenance does not identify a final M1 max-range-4 fit")
}
if (!isTRUE(fit_optimizer$inputs$max_range_size == 4L) || !isTRUE(fit_optimizer$inputs$include_null_range) ||
    !identical(fit_optimizer$inputs$force_sparse, FALSE) || !identical(fit_optimizer$inputs$speedup, FALSE)) {
  stop("M1 fit settings fail the final-fit gate")
}

# BSM consumes the independently recalculated, fixed-parameter ancestral-state
# object, while retaining the accepted optimizer result as the inferential
# source of d/e and lnL.
res <- readRDS(postfit_rds)
if (!is.list(res) || is.null(res$inputs) || is.null(res$outputs)) stop("Invalid post-fit result structure")
if (!is.finite(as.numeric(res$total_loglikelihood))) stop("Post-fit object has non-finite total_loglikelihood")
if (!isS4(res$outputs) || !("params_table" %in% slotNames(res$outputs))) {
  stop("Post-fit object lacks outputs@params_table")
}
pt_postfit <- res$outputs@params_table
if (!all(c("d", "e", "j", "w") %in% rownames(pt_postfit)) ||
    !isTRUE(all.equal(as.numeric(pt_postfit[c("d", "e"), "est"]),
                      as.numeric(pt_fit[c("d", "e"), "est"]), tolerance=1e-12)) ||
    !isTRUE(as.numeric(pt_postfit["j", "est"]) == 0) ||
    !isTRUE(as.numeric(pt_postfit["w", "est"]) == 1)) {
  stop("Post-fit parameters do not match the accepted optimizer fit")
}
ll_optimizer <- as.numeric(fit_optimizer$total_loglikelihood)
ll_postfit <- as.numeric(res$total_loglikelihood)
if (abs(ll_postfit - ll_optimizer) > 1e-6) {
  stop(sprintf("Post-fit likelihood differs from accepted optimizer: %.12g", ll_postfit - ll_optimizer))
}
postfit_prov <- res$inputs$postfit_provenance
if (is.null(postfit_prov) ||
    !identical(postfit_prov$real_optimizer_fit_sha256, fit_hash) ||
    is.null(res$postfit_completion) ||
    !identical(res$postfit_completion$real_optimizer_fit_sha256, fit_hash)) {
  stop("Post-fit object does not preserve the accepted optimizer hash")
}
if (!isTRUE(res$inputs$return_condlikes_table) || !isTRUE(res$inputs$calc_ancprobs)) {
  stop("Post-fit object was not recalculated with conditional tables and ancestral probabilities")
}

result_input_paths <- list(
  tree=normalizePath(res$inputs$trfn, mustWork=TRUE),
  geog=normalizePath(res$inputs$geogfn, mustWork=TRUE),
  times=normalizePath(res$inputs$timesfn, mustWork=TRUE),
  multipliers=normalizePath(res$inputs$dispersal_multipliers_fn, mustWork=TRUE)
)
optimizer_input_paths <- list(
  tree=normalizePath(fit_optimizer$inputs$trfn, mustWork=TRUE),
  geog=normalizePath(fit_optimizer$inputs$geogfn, mustWork=TRUE),
  times=normalizePath(fit_optimizer$inputs$timesfn, mustWork=TRUE),
  multipliers=normalizePath(fit_optimizer$inputs$dispersal_multipliers_fn, mustWork=TRUE)
)
for (nm in names(result_input_paths)) {
  if (!identical(result_input_paths[[nm]], expected_paths[[nm]]) ||
      !identical(optimizer_input_paths[[nm]], expected_paths[[nm]]) ||
      !identical(sha256_file(result_input_paths[[nm]]), unname(expected_hashes[[nm]])) ||
      !identical(sha256_file(optimizer_input_paths[[nm]]), unname(expected_hashes[[nm]]))) {
    stop("Optimizer/post-fit input mismatch for ", nm)
  }
}
expected_times <- c(3.3, 13.9, 23.03, 33.9, 45)
if (!isTRUE(all.equal(as.numeric(res$inputs$timeperiods), expected_times, tolerance=0))) {
  stop("Result does not contain the five frozen M1 time periods")
}

tr <- ape::read.tree(expected_paths$tree)
geog <- BioGeoBEARS::getranges_from_LagrangePHYLIP(expected_paths$geog)
if (length(tr$tip.label) != 417L || tr$Nnode != 416L || nrow(geog@df) != 417L ||
    !setequal(tr$tip.label, rownames(geog@df)) ||
    !identical(colnames(geog@df), internal_codes)) {
  stop("Tree/geography/area-code gate failed")
}

# Time-stratified BSM requires these post-fit arrays.  A fit lacking them must
# be regenerated; fabricating them after acceptance would break provenance.
required_result_fields <- c(
  "condlikes_table",
  "relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS_TABLE",
  "ML_marginal_prob_each_state_at_branch_top_AT_node"
)
for (nm in required_result_fields) {
  if (is.null(res[[nm]]) || length(res[[nm]]) == 0L) stop("BSM-required result field is absent: ", nm)
}
for (nm in c("master_table", "tree_sections_list")) {
  if (is.null(res$inputs[[nm]]) || length(res$inputs[[nm]]) == 0L) stop("BSM-required input field is absent: ", nm)
}

required_bgb_functions <- c(
  "get_inputs_for_stochastic_mapping_stratified",
  "stochastic_mapping_on_stratified",
  "get_dmat_times_from_res",
  "simulate_source_area_ana",
  "simulate_source_area_clado",
  "get_huge_events_tables_from_clado_ana_events_tables",
  "count_events_huge_tables",
  "count_ana_dispersal_events"
)
function_hashes <- setNames(vapply(required_bgb_functions, function(nm) {
  sha256_object(list(formals=formals(bgb_fun(nm)), body=body(bgb_fun(nm))))
}, character(1)), required_bgb_functions)

# Confirm the exact installed 1.1.3 extinction-counter defect.  This is an
# executable compatibility test, not a replacement for the independent count.
synthetic <- data.frame(
  abs_event_time=c(1, 2), event_type=c("d", "e"),
  current_rangetxt=c("A", "AB"), new_rangetxt=c("AB", "A"),
  ana_dispersal_from=c("A", ""), dispersal_to=c("B", "-"),
  extirpation_from=c("-", "B"), stringsAsFactors=FALSE
)
native_bug_probe <- bgb_fun("count_ana_dispersal_events")(
  ana_events_table=synthetic, areanames=c("A", "B"), actual_names=c("x", "y")
)
native_probe_e <- sum(as.numeric(unlist(native_bug_probe$e_counts_df)), na.rm=TRUE)
if (!identical(native_probe_e, 0)) {
  stop("The known count_ana_dispersal_events bug is no longer reproduced; re-audit before running")
}

preflight <- data.frame(
  check=c("model_scope_final_QA", "acceptance", "fit_and_postfit_hashes", "M1_final_recipe", "tree_geog", "five_strata",
          "BSM_arrays", "area_order", "BioGeoBEARS_version", "native_e_bug_probe"),
  status="PASS",
  detail=c(
    paste0("M1_R4 PASS/COMPLETE;qa=", final_qa_hash,
           ";validator=", qa_bound_artifacts$validator$sha256,
           ";inspector=", qa_bound_artifacts$rds_inspector$sha256,
           ";spec=", qa_bound_artifacts$spec$sha256,
           ";manifest=", qa_bound_artifacts$manifest$sha256),
    qa_acceptance_hash, paste0("optimizer=", fit_hash, ";postfit=", postfit_hash),
    sprintf("optimizer_lnL=%.12f;postfit_lnL=%.12f;delta=%.3g;d=%.12g;e=%.12g;j=0;w=1",
            ll_optimizer, ll_postfit, ll_postfit - ll_optimizer,
            as.numeric(pt_fit["d", "est"]), as.numeric(pt_fit["e", "est"])),
    "417 tips;416 internal;11 A:K areas", paste(expected_times, collapse=","),
    paste(required_result_fields, collapse=","), paste(paste(internal_codes, actual_names, sep="="), collapse=","),
    paste0("1.1.3@", remote_sha), "synthetic true e=1; native observed e=0; independent counter required"
  ), stringsAsFactors=FALSE
)

cat("BSM_PREFLIGHT_PASS\n")
print(preflight, row.names=FALSE)
if (identical(action, "preflight")) quit(status=0L)

epoch_bounds <- c(0, expected_times)
epoch_labels <- c("0-3.3", "3.3-13.9", "13.9-23.03", "23.03-33.9", "33.9-45")
epoch_for_age <- function(age) {
  idx <- findInterval(age, epoch_bounds, rightmost.closed=FALSE, all.inside=FALSE)
  if (length(idx) != 1L || is.na(idx) || !(idx %in% seq_along(epoch_labels))) {
    stop("Event age is outside the frozen epochs: ", age)
  }
  idx
}

tree_exposure <- function(tree) {
  depths <- ape::node.depth.edgelength(tree)
  root_depth <- max(depths[seq_along(tree$tip.label)])
  ages <- root_depth - depths
  if (any(ages[tree$edge[, 1]] + 1e-8 < ages[tree$edge[, 2]])) stop("Tree age direction is inconsistent")
  out <- lapply(seq_along(epoch_labels), function(i) {
    lo <- epoch_bounds[[i]]
    hi <- epoch_bounds[[i + 1L]]
    older <- ages[tree$edge[, 1]]
    younger <- ages[tree$edge[, 2]]
    overlap <- pmax(0, pmin(older, hi) - pmax(younger, lo))
    data.frame(epoch=epoch_labels[[i]], young_ma=lo, old_ma=hi,
               lineage_myr_exposure=sum(overlap), stringsAsFactors=FALSE)
  })
  z <- do.call(rbind, out)
  if (any(!is.finite(z$lineage_myr_exposure)) || any(z$lineage_myr_exposure <= 0)) {
    stop("Invalid lineage-time exposure")
  }
  z
}
exposure <- tree_exposure(tr)

parse_range <- function(x) {
  if (length(x) != 1L || is.na(x) || x %in% c("", "-", "_")) return(character())
  chars <- strsplit(as.character(x), "", fixed=TRUE)[[1]]
  if (anyDuplicated(chars) || any(!(chars %in% internal_codes))) stop("Invalid range text: ", x)
  chars
}

standardize_anagenetic <- function(tab, map_id) {
  empty <- data.frame(
    map_id=integer(), event_id=integer(), age_ma=numeric(), epoch=character(),
    event_type=character(), current_range=character(), new_range=character(),
    from_code=character(), to_code=character(), affected_code=character(),
    from_area=character(), to_area=character(), affected_area=character(),
    stringsAsFactors=FALSE
  )
  if (is.null(tab) || (length(tab) == 1L && is.na(tab))) return(empty)
  if (!is.data.frame(tab)) stop("Anagenetic events object is neither data.frame nor NA")
  if (nrow(tab) == 0L) return(empty)
  needed <- c("abs_event_time", "event_type", "current_rangetxt", "new_rangetxt",
              "ana_dispersal_from", "dispersal_to", "extirpation_from")
  if (!all(needed %in% names(tab))) stop("Anagenetic event schema mismatch")
  rows <- vector("list", nrow(tab))
  for (i in seq_len(nrow(tab))) {
    typ <- as.character(tab$event_type[[i]])
    if (!(typ %in% c("d", "e"))) stop("Unexpected M1 anagenetic event_type: ", typ)
    cur <- parse_range(tab$current_rangetxt[[i]])
    new <- parse_range(tab$new_rangetxt[[i]])
    gained <- setdiff(new, cur)
    lost <- setdiff(cur, new)
    age <- as.numeric(tab$abs_event_time[[i]])
    if (!is.finite(age) || age < 0 || age >= 45) stop("Invalid anagenetic event age")
    ep <- epoch_for_age(age)
    if (all(c("time_top", "time_bot") %in% names(tab))) {
      top <- as.numeric(tab$time_top[[i]])
      bot <- as.numeric(tab$time_bot[[i]])
      if (is.finite(top) && is.finite(bot) && !(age >= top - 1e-8 && age < bot + 1e-8)) {
        stop("Anagenetic event falls outside its source stratum")
      }
    }
    from <- to <- affected <- ""
    if (typ == "d") {
      if (length(gained) != 1L || length(lost) != 0L) stop("d event is not a one-area gain")
      to <- gained[[1]]
      if (!identical(as.character(tab$dispersal_to[[i]]), to)) stop("dispersal_to disagrees with range delta")
      from <- as.character(tab$ana_dispersal_from[[i]])
      if (length(from) != 1L || !(from %in% cur)) stop("Imputed source is not in the ancestral range")
      affected <- to
    } else {
      if (length(lost) != 1L || length(gained) != 0L) stop("e event is not a one-area loss")
      affected <- lost[[1]]
      if (!identical(as.character(tab$extirpation_from[[i]]), affected)) stop("extirpation_from disagrees with range delta")
    }
    rows[[i]] <- data.frame(
      map_id=map_id, event_id=i, age_ma=age, epoch=epoch_labels[[ep]], event_type=typ,
      current_range=paste(cur, collapse=""), new_range=paste(new, collapse=""),
      from_code=from, to_code=to, affected_code=affected,
      from_area=if (nzchar(from)) unname(code_to_area[[from]]) else "",
      to_area=if (nzchar(to)) unname(code_to_area[[to]]) else "",
      affected_area=unname(code_to_area[[affected]]), stringsAsFactors=FALSE
    )
  }
  do.call(rbind, rows)
}

standardize_cladogenetic <- function(tab, map_id) {
  if (!is.data.frame(tab)) stop("Cladogenetic events object must be a data.frame")
  needed <- c("node", "time_bp", "clado_event_type")
  if (!all(needed %in% names(tab))) stop("Cladogenetic event schema mismatch")
  typ <- trimws(as.character(tab$clado_event_type))
  keep <- !is.na(tab$clado_event_type) & nzchar(typ)
  z <- tab[keep, , drop=FALSE]
  typ <- typ[keep]
  expected_nodes <- (length(tr$tip.label) + 1L):(length(tr$tip.label) + tr$Nnode)
  nodes <- as.integer(z$node)
  if (nrow(z) != tr$Nnode || anyDuplicated(nodes) || !setequal(nodes, expected_nodes)) {
    stop("Each successful map must contain exactly 416 unique cladogenetic node events")
  }
  ages <- as.numeric(z$time_bp)
  if (any(!is.finite(ages)) || any(ages < 0) || any(ages >= 45)) stop("Invalid cladogenetic event age")
  eps <- vapply(ages, epoch_for_age, integer(1))
  out <- data.frame(map_id=map_id, event_id=seq_len(nrow(z)), node=nodes,
                    age_ma=ages, epoch=epoch_labels[eps], event_type=typ,
                    event_text=if ("clado_event_txt" %in% names(z)) as.character(z$clado_event_txt) else "",
                    stringsAsFactors=FALSE)
  if (any(grepl("\\(j\\)", out$event_type))) stop("Founder event observed although j is fixed at zero")
  out
}

native_total_counts <- function(clado, ana) {
  huge <- bgb_fun("get_huge_events_tables_from_clado_ana_events_tables")(
    clado_events_tables=list(clado), ana_events_tables=list(ana), model_name="M1"
  )
  x <- bgb_fun("count_events_huge_tables")(huge_tables=huge, BSM_i=1)
  c(d=as.numeric(x$d[[1]]), e=as.numeric(x$e[[1]]), j=as.numeric(x$j[[1]]))
}

make_seed <- function(base, offset) {
  as.integer(((as.double(base) + as.double(offset) - 1) %% 2147483646) + 1)
}

build_one_map <- function(map_id, attempt_id, bsm_inputs) {
  map_seed <- make_seed(seed_base, attempt_id)
  source_seed <- make_seed(seed_base, 1000000 + attempt_id)
  sm <- bgb_fun("stochastic_mapping_on_stratified")(
    res=res, stochastic_mapping_inputs_list=bsm_inputs, maxtries=maxtries,
    seedval=map_seed, master_nodenum_toPrint=0
  )
  if (!is.list(sm) || !is.data.frame(sm$master_table_cladogenetic_events)) {
    stop("Unrecognized stratified stochastic-map output")
  }
  clado_raw <- sm$master_table_cladogenetic_events
  ana_raw <- sm$table_w_anagenetic_events
  if (is.null(ana_raw)) ana_raw <- NA

  dmat_times <- bgb_fun("get_dmat_times_from_res")(res=res, numstates=NULL)
  set.seed(source_seed)
  clado_src <- bgb_fun("simulate_source_area_clado")(
    clado_events_table=clado_raw, areanames=internal_codes,
    dmat=dmat_times$dmat, times=dmat_times$times
  )
  ana_src <- if (is.data.frame(ana_raw)) {
    bgb_fun("simulate_source_area_ana")(
      ana_events_table=ana_raw, areanames=internal_codes,
      dmat=dmat_times$dmat, times=dmat_times$times
    )
  } else ana_raw

  ana_std <- standardize_anagenetic(ana_src, map_id)
  clado_std <- standardize_cladogenetic(clado_src, map_id)
  indep <- c(
    d=sum(ana_std$event_type == "d"),
    e=sum(ana_std$event_type == "e"),
    j=sum(grepl("\\(j\\)", clado_std$event_type))
  )
  native <- native_total_counts(clado_src, ana_src)
  if (!identical(as.numeric(indep), as.numeric(native))) {
    stop("Independent d/e/j totals disagree with BioGeoBEARS count_events_huge_tables")
  }
  if (!isTRUE(indep[["j"]] == 0)) stop("j invariant failed")

  list(
    schema_version="1.0",
    metadata=list(
      map_id=map_id, attempt_id=attempt_id, map_seed=map_seed,
      source_seed=source_seed, maxtries=maxtries,
      accepted_optimizer_fit_sha256=fit_hash, postfit_ancestral_states_sha256=postfit_hash,
      final_qa_sha256=final_qa_hash, acceptance_sha256=qa_acceptance_hash,
      created=timestamp(),
      area_order_sha256=expected_hashes[["area_order"]],
      independent_counts=as.list(indep), native_counts=as.list(native)
    ),
    clado_events_raw=clado_raw,
    ana_events_raw=ana_raw,
    clado_events_with_source=clado_src,
    ana_events_with_source=ana_src,
    independent_cladogenetic_events=clado_std,
    independent_anagenetic_events=ana_std
  )
}

map_files <- function() {
  if (!dir.exists(file.path(outdir, "maps"))) return(character())
  sort(list.files(file.path(outdir, "maps"), pattern="^map_[0-9]{4}\\.rds$", full.names=TRUE))
}

load_maps <- function(n=NULL) {
  fs <- map_files()
  if (!is.null(n)) {
    if (length(fs) < n) stop("Only ", length(fs), " maps exist; requested ", n)
    fs <- fs[seq_len(n)]
  }
  maps <- lapply(fs, readRDS)
  ids <- vapply(maps, function(x) as.integer(x$metadata$map_id), integer(1))
  if (!identical(ids, seq_along(maps))) stop("Map IDs are not unique and contiguous")
  attempt_ids <- vapply(maps, function(x) as.integer(x$metadata$attempt_id), integer(1))
  if (anyNA(attempt_ids) || anyDuplicated(attempt_ids) || is.unsorted(attempt_ids, strictly=TRUE)) {
    stop("Successful-map attempt IDs are absent, duplicated, or not strictly increasing")
  }
  map_seeds <- vapply(maps, function(x) as.integer(x$metadata$map_seed), integer(1))
  source_seeds <- vapply(maps, function(x) as.integer(x$metadata$source_seed), integer(1))
  recorded_maxtries <- vapply(maps, function(x) as.integer(x$metadata$maxtries), integer(1))
  if (!identical(map_seeds, vapply(attempt_ids, function(i) make_seed(seed_base, i), integer(1))) ||
      !identical(source_seeds, vapply(attempt_ids, function(i) make_seed(seed_base, 1000000 + i), integer(1)))) {
    stop("A checkpoint map violates the deterministic map/source seed schedule")
  }
  if (anyNA(recorded_maxtries) || any(recorded_maxtries != maxtries)) {
    stop("A checkpoint map was generated with a different maxtries value")
  }
  if (any(vapply(maps, function(x) {
    !identical(x$metadata$accepted_optimizer_fit_sha256, fit_hash) ||
      !identical(x$metadata$postfit_ancestral_states_sha256, postfit_hash) ||
      !identical(x$metadata$final_qa_sha256, final_qa_hash) ||
      !identical(x$metadata$acceptance_sha256, qa_acceptance_hash)
  }, logical(1)))) {
    stop("A checkpoint map belongs to another optimizer/post-fit/FINAL_QA/acceptance bundle")
  }
  maps
}

summarize_vector <- function(x) {
  data.frame(mean=mean(x), sd=if (length(x) > 1L) sd(x) else NA_real_,
             q025=as.numeric(quantile(x, 0.025, names=FALSE)),
             q975=as.numeric(quantile(x, 0.975, names=FALSE)),
             mcse=if (length(x) > 1L) sd(x) / sqrt(length(x)) else NA_real_,
             interval_definition="empirical map-to-map quantiles; not a CI for the mean",
             stringsAsFactors=FALSE)
}

precision_metric <- function(mean, mcse) {
  if (!is.finite(mean) || !is.finite(mcse)) return(Inf)
  if (mean >= 1) mcse / mean else mcse
}

precision_ok <- function(mean, mcse) precision_metric(mean, mcse) <= 0.05

change_metric <- function(old, new) {
  if (!is.finite(old) || !is.finite(new)) return(Inf)
  if (max(abs(old), abs(new)) >= 1) abs(new - old) / max(abs(old), 1e-12) else abs(new - old)
}

change_ok <- function(old, new) change_metric(old, new) <= 0.05

max_or_zero <- function(x) if (length(x) == 0L) 0 else max(x)

top_keys <- function(df, key, value, n=3L) {
  z <- df[df[[value]] > 0, c(key, value), drop=FALSE]
  if (nrow(z) == 0L) return(character())
  z <- z[order(-z[[value]], z[[key]]), , drop=FALSE]
  sort(head(z[[key]], n))
}

summary_dir_for <- function(n) file.path(outdir, "summaries", sprintf("n_%04d", n))

read_checkpoint_means <- function(n) {
  sd <- summary_dir_for(n)
  list(
    totals=read.delim(file.path(sd, "event_totals_summary.tsv"), stringsAsFactors=FALSE),
    routes=read.delim(file.path(sd, "d_route_summary.tsv"), stringsAsFactors=FALSE),
    ext=read.delim(file.path(sd, "e_area_summary.tsv"), stringsAsFactors=FALSE)
  )
}

summarize_checkpoint <- function(n) {
  final_dir <- summary_dir_for(n)
  if (dir.exists(final_dir) && file.exists(file.path(final_dir, "summary_status.tsv"))) {
    return(read.delim(file.path(final_dir, "summary_status.tsv"), stringsAsFactors=FALSE))
  }
  tmp_dir <- paste0(final_dir, ".tmp.", Sys.getpid())
  if (!dir.create(tmp_dir, recursive=TRUE, showWarnings=FALSE)) stop("Cannot create summary temp directory")
  maps <- load_maps(n)
  attempt_audit <- validate_attempt_records(maps)

  per_map <- do.call(rbind, lapply(maps, function(x) {
    a <- x$independent_anagenetic_events
    ctab <- x$independent_cladogenetic_events
    data.frame(
      map_id=as.integer(x$metadata$map_id), attempt_id=as.integer(x$metadata$attempt_id),
      map_seed=as.integer(x$metadata$map_seed), source_seed=as.integer(x$metadata$source_seed),
      maxtries=as.integer(x$metadata$maxtries),
      d=sum(a$event_type == "d"), e=sum(a$event_type == "e"),
      j=sum(grepl("\\(j\\)", ctab$event_type)), cladogenetic_events=nrow(ctab),
      stringsAsFactors=FALSE
    )
  }))

  ana_all <- do.call(rbind, lapply(maps, `[[`, "independent_anagenetic_events"))
  clado_all <- do.call(rbind, lapply(maps, `[[`, "independent_cladogenetic_events"))

  totals_summary <- do.call(rbind, lapply(c("d", "e", "j"), function(ev) {
    cbind(data.frame(event=ev, n_maps=n, stringsAsFactors=FALSE), summarize_vector(per_map[[ev]]))
  }))

  route_grid <- expand.grid(map_id=seq_len(n), from_code=internal_codes,
                            to_code=internal_codes, stringsAsFactors=FALSE)
  route_grid <- route_grid[route_grid$from_code != route_grid$to_code, ]
  route_obs <- if (nrow(ana_all) > 0L) {
    z <- ana_all[ana_all$event_type == "d", ]
    as.data.frame(table(map_id=z$map_id, from_code=z$from_code, to_code=z$to_code),
                  stringsAsFactors=FALSE)
  } else data.frame(map_id=integer(), from_code=character(), to_code=character(), Freq=integer())
  route_obs$map_id <- as.integer(as.character(route_obs$map_id))
  routes <- merge(route_grid, route_obs, by=c("map_id", "from_code", "to_code"), all.x=TRUE, sort=FALSE)
  routes$count <- ifelse(is.na(routes$Freq), 0L, as.integer(routes$Freq))
  routes$Freq <- NULL
  routes$from_area <- unname(code_to_area[routes$from_code])
  routes$to_area <- unname(code_to_area[routes$to_code])
  routes <- routes[, c("map_id", "from_code", "to_code", "from_area", "to_area", "count")]

  ext_grid <- expand.grid(map_id=seq_len(n), area_code=internal_codes, stringsAsFactors=FALSE)
  ext_obs <- if (nrow(ana_all) > 0L) {
    z <- ana_all[ana_all$event_type == "e", ]
    as.data.frame(table(map_id=z$map_id, area_code=z$affected_code), stringsAsFactors=FALSE)
  } else data.frame(map_id=integer(), area_code=character(), Freq=integer())
  ext_obs$map_id <- as.integer(as.character(ext_obs$map_id))
  ext <- merge(ext_grid, ext_obs, by=c("map_id", "area_code"), all.x=TRUE, sort=FALSE)
  ext$count <- ifelse(is.na(ext$Freq), 0L, as.integer(ext$Freq))
  ext$Freq <- NULL
  ext$area <- unname(code_to_area[ext$area_code])
  ext <- ext[, c("map_id", "area_code", "area", "count")]

  period_grid <- expand.grid(map_id=seq_len(n), epoch=epoch_labels,
                             event=c("d", "e"), stringsAsFactors=FALSE)
  period_obs <- if (nrow(ana_all) > 0L) {
    as.data.frame(table(map_id=ana_all$map_id, epoch=ana_all$epoch,
                        event=ana_all$event_type), stringsAsFactors=FALSE)
  } else data.frame(map_id=integer(), epoch=character(), event=character(), Freq=integer())
  period_obs$map_id <- as.integer(as.character(period_obs$map_id))
  period <- merge(period_grid, period_obs, by=c("map_id", "epoch", "event"), all.x=TRUE, sort=FALSE)
  period$count <- ifelse(is.na(period$Freq), 0L, as.integer(period$Freq))
  period$Freq <- NULL
  period <- merge(period, exposure, by="epoch", all.x=TRUE, sort=FALSE)
  period$rate_per_lineage_myr <- period$count / period$lineage_myr_exposure
  period <- period[, c("map_id", "epoch", "young_ma", "old_ma", "event", "count",
                       "lineage_myr_exposure", "rate_per_lineage_myr")]

  route_period_grid <- expand.grid(
    map_id=seq_len(n), epoch=epoch_labels, from_code=internal_codes,
    to_code=internal_codes, stringsAsFactors=FALSE
  )
  route_period_grid <- route_period_grid[route_period_grid$from_code != route_period_grid$to_code, ]
  route_period_obs <- if (nrow(ana_all) > 0L) {
    z <- ana_all[ana_all$event_type == "d", ]
    as.data.frame(table(map_id=z$map_id, epoch=z$epoch,
                        from_code=z$from_code, to_code=z$to_code), stringsAsFactors=FALSE)
  } else data.frame(map_id=integer(), epoch=character(), from_code=character(),
                    to_code=character(), Freq=integer())
  route_period_obs$map_id <- as.integer(as.character(route_period_obs$map_id))
  routes_period <- merge(route_period_grid, route_period_obs,
                         by=c("map_id", "epoch", "from_code", "to_code"),
                         all.x=TRUE, sort=FALSE)
  routes_period$count <- ifelse(is.na(routes_period$Freq), 0L, as.integer(routes_period$Freq))
  routes_period$Freq <- NULL
  routes_period$from_area <- unname(code_to_area[routes_period$from_code])
  routes_period$to_area <- unname(code_to_area[routes_period$to_code])
  routes_period <- merge(routes_period, exposure[, c("epoch", "young_ma", "old_ma")],
                         by="epoch", all.x=TRUE, sort=FALSE)
  routes_period <- routes_period[, c("map_id", "epoch", "young_ma", "old_ma",
                                     "from_code", "to_code", "from_area", "to_area", "count")]

  ext_period_grid <- expand.grid(map_id=seq_len(n), epoch=epoch_labels,
                                 area_code=internal_codes, stringsAsFactors=FALSE)
  ext_period_obs <- if (nrow(ana_all) > 0L) {
    z <- ana_all[ana_all$event_type == "e", ]
    as.data.frame(table(map_id=z$map_id, epoch=z$epoch,
                        area_code=z$affected_code), stringsAsFactors=FALSE)
  } else data.frame(map_id=integer(), epoch=character(), area_code=character(), Freq=integer())
  ext_period_obs$map_id <- as.integer(as.character(ext_period_obs$map_id))
  ext_period <- merge(ext_period_grid, ext_period_obs,
                      by=c("map_id", "epoch", "area_code"), all.x=TRUE, sort=FALSE)
  ext_period$count <- ifelse(is.na(ext_period$Freq), 0L, as.integer(ext_period$Freq))
  ext_period$Freq <- NULL
  ext_period$area <- unname(code_to_area[ext_period$area_code])
  ext_period <- merge(ext_period, exposure[, c("epoch", "young_ma", "old_ma")],
                      by="epoch", all.x=TRUE, sort=FALSE)
  ext_period <- ext_period[, c("map_id", "epoch", "young_ma", "old_ma",
                               "area_code", "area", "count")]

  route_summary <- do.call(rbind, lapply(split(routes, interaction(routes$from_code, routes$to_code, drop=TRUE)), function(z) {
    cbind(z[1, c("from_code", "to_code", "from_area", "to_area")], n_maps=n,
          summarize_vector(z$count))
  }))
  rownames(route_summary) <- NULL
  route_summary$route <- paste(route_summary$from_area, route_summary$to_area, sep="->")
  ext_summary <- do.call(rbind, lapply(split(ext, ext$area_code), function(z) {
    cbind(z[1, c("area_code", "area")], n_maps=n, summarize_vector(z$count))
  }))
  rownames(ext_summary) <- NULL
  temporal_summary <- do.call(rbind, lapply(split(period, interaction(period$epoch, period$event, drop=TRUE)), function(z) {
    s1 <- summarize_vector(z$count)
    s2 <- summarize_vector(z$rate_per_lineage_myr)
    data.frame(epoch=z$epoch[[1]], young_ma=z$young_ma[[1]], old_ma=z$old_ma[[1]],
               event=z$event[[1]], n_maps=n, lineage_myr_exposure=z$lineage_myr_exposure[[1]],
               mean_count=s1$mean, sd_count=s1$sd, q025_count=s1$q025, q975_count=s1$q975,
               mcse_count=s1$mcse, mean_rate=s2$mean, sd_rate=s2$sd,
               q025_rate=s2$q025, q975_rate=s2$q975, mcse_rate=s2$mcse,
               interval_definition=s1$interval_definition, stringsAsFactors=FALSE)
  }))
  rownames(temporal_summary) <- NULL

  native_x <- do.call(rbind, lapply(maps, function(x) {
    data.frame(map_id=x$metadata$map_id,
               independent_d=x$metadata$independent_counts$d,
               native_d=x$metadata$native_counts$d,
               independent_e=x$metadata$independent_counts$e,
               native_e=x$metadata$native_counts$e,
               independent_j=x$metadata$independent_counts$j,
               native_j=x$metadata$native_counts$j,
               totals_match=identical(as.numeric(unlist(x$metadata$independent_counts)),
                                      as.numeric(unlist(x$metadata$native_counts))),
               native_area_e_counter_used=FALSE,
               note="count_ana_dispersal_events e-by-area bypassed: installed 1.1.3 references events_df2$extirpation_from",
               stringsAsFactors=FALSE)
  }))

  conservation <- data.frame(
    map_id=per_map$map_id,
    d_route_sum=vapply(seq_len(n), function(i) sum(routes$count[routes$map_id == i]), numeric(1)),
    d_total=per_map$d,
    e_area_sum=vapply(seq_len(n), function(i) sum(ext$count[ext$map_id == i]), numeric(1)),
    e_total=per_map$e,
    j_total=per_map$j,
    cladogenetic_events=per_map$cladogenetic_events,
    native_totals_match=native_x$totals_match,
    stringsAsFactors=FALSE
  )
  period_conservation <- merge(
    aggregate(count ~ map_id + epoch, routes_period, sum),
    period[period$event == "d", c("map_id", "epoch", "count")],
    by=c("map_id", "epoch"), suffixes=c("_route", "_d"), all=TRUE
  )
  period_conservation <- merge(
    period_conservation,
    aggregate(count ~ map_id + epoch, ext_period, sum),
    by=c("map_id", "epoch"), all=TRUE
  )
  names(period_conservation)[names(period_conservation) == "count"] <- "count_extinction_area"
  period_e <- period[period$event == "e", c("map_id", "epoch", "count")]
  names(period_e)[names(period_e) == "count"] <- "count_e"
  period_conservation <- merge(period_conservation, period_e,
                               by=c("map_id", "epoch"), all=TRUE)
  period_conservation$pass <- with(period_conservation,
    count_route == count_d & count_extinction_area == count_e)
  period_pass_by_map <- tapply(period_conservation$pass, period_conservation$map_id, all)
  conservation$period_conservation_pass <- as.logical(period_pass_by_map[as.character(conservation$map_id)])
  conservation$pass <- with(conservation,
    d_route_sum == d_total & e_area_sum == e_total & j_total == 0 &
      cladogenetic_events == 416 & native_totals_match & period_conservation_pass)
  if (!all(conservation$pass)) stop("Event conservation/native cross-check failed")

  totals_de <- totals_summary[totals_summary$event %in% c("d", "e"), ]
  total_precision_metrics <- vapply(seq_len(nrow(totals_de)), function(i) {
    z <- totals_de[i, ]
    precision_metric(z$mean, z$mcse)
  }, numeric(1))
  total_precision <- all(total_precision_metrics <= 0.05)
  temporal_precision_metrics <- vapply(seq_len(nrow(temporal_summary)), function(i) {
    precision_metric(temporal_summary$mean_count[[i]], temporal_summary$mcse_count[[i]])
  }, numeric(1))
  temporal_precision <- all(temporal_precision_metrics <= 0.05)
  top_route_now <- top_keys(route_summary, "route", "mean", 3L)
  top_ext_now <- top_keys(ext_summary, "area", "mean", 3L)
  route_precision_metrics <- vapply(which(route_summary$route %in% top_route_now), function(i) {
    precision_metric(route_summary$mean[[i]], route_summary$mcse[[i]])
  }, numeric(1))
  route_top_precision <- all(route_precision_metrics <= 0.05)
  ext_precision_metrics <- vapply(which(ext_summary$area %in% top_ext_now), function(i) {
    precision_metric(ext_summary$mean[[i]], ext_summary$mcse[[i]])
  }, numeric(1))
  ext_top_precision <- all(ext_precision_metrics <= 0.05)

  vs_prev_totals <- vs_prev_routes <- vs_prev_ext <- FALSE
  vs_100_totals <- vs_100_routes <- vs_100_ext <- FALSE
  max_change_vs_prev <- max_change_vs_100 <- NA_real_
  if (n >= 200L) {
    prev <- read_checkpoint_means(n - 100L)
    base <- read_checkpoint_means(100L)
    curr_de <- totals_summary[totals_summary$event %in% c("d", "e"), ]
    prev_de <- prev$totals[prev$totals$event %in% c("d", "e"), ]
    base_de <- base$totals[base$totals$event %in% c("d", "e"), ]
    curr_de <- curr_de[match(c("d", "e"), curr_de$event), ]
    prev_de <- prev_de[match(c("d", "e"), prev_de$event), ]
    base_de <- base_de[match(c("d", "e"), base_de$event), ]
    changes_prev <- mapply(change_metric, prev_de$mean, curr_de$mean)
    changes_100 <- mapply(change_metric, base_de$mean, curr_de$mean)
    max_change_vs_prev <- max(changes_prev)
    max_change_vs_100 <- max(changes_100)
    vs_prev_totals <- all(changes_prev <= 0.05)
    vs_100_totals <- all(changes_100 <= 0.05)
    vs_prev_routes <- identical(top_route_now, top_keys(prev$routes, "route", "mean", 3L))
    vs_prev_ext <- identical(top_ext_now, top_keys(prev$ext, "area", "mean", 3L))
    vs_100_routes <- identical(top_route_now, top_keys(base$routes, "route", "mean", 3L))
    vs_100_ext <- identical(top_ext_now, top_keys(base$ext, "area", "mean", 3L))
  }
  overall_stable <- n >= 200L && total_precision && temporal_precision &&
    route_top_precision && ext_top_precision && vs_prev_totals && vs_prev_routes &&
    vs_prev_ext && vs_100_totals && vs_100_routes && vs_100_ext

  convergence_current <- data.frame(
    n_maps=n, total_precision_pass=total_precision,
    max_total_precision_metric=max_or_zero(total_precision_metrics),
    temporal_precision_pass=temporal_precision,
    max_temporal_precision_metric=max_or_zero(temporal_precision_metrics),
    top3_route_precision_pass=route_top_precision,
    max_top3_route_precision_metric=max_or_zero(route_precision_metrics),
    top3_extinction_precision_pass=ext_top_precision,
    max_top3_extinction_precision_metric=max_or_zero(ext_precision_metrics),
    totals_change_vs_previous_100_pass=vs_prev_totals,
    max_total_change_metric_vs_previous_100=max_change_vs_prev,
    top3_routes_vs_previous_100_stable=vs_prev_routes,
    top3_extinction_vs_previous_100_stable=vs_prev_ext,
    totals_change_vs_100_pass=vs_100_totals,
    max_total_change_metric_vs_100=max_change_vs_100,
    top3_routes_vs_100_stable=vs_100_routes,
    top3_extinction_vs_100_stable=vs_100_ext,
    top3_routes=paste(top_route_now, collapse=";"),
    top3_extinction_areas=paste(top_ext_now, collapse=";"),
    overall_stable=overall_stable,
    metric_threshold=0.05,
    checkpoint_increment=100L,
    maximum_maps=500L,
    criteria="MCSE/mean<=0.05 if mean>=1 else MCSE<=0.05; cumulative mean change<=5% (or abs<=0.05 below 1); identical top-3 sets; comparisons to prior 100 and initial 100",
    stringsAsFactors=FALSE
  )
  if (n >= 200L) {
    prior_convergence_path <- file.path(summary_dir_for(n - 100L), "bsm_convergence.tsv")
    if (!file.exists(prior_convergence_path)) stop("Previous checkpoint lacks bsm_convergence.tsv")
    prior_convergence <- read.delim(prior_convergence_path, stringsAsFactors=FALSE,
                                    check.names=FALSE)
    if (!identical(names(prior_convergence), names(convergence_current))) {
      stop("Previous convergence table schema differs from the current runner")
    }
    convergence <- rbind(prior_convergence, convergence_current)
  } else {
    convergence <- convergence_current
  }
  if (!identical(as.integer(convergence$n_maps), seq(100L, n, 100L))) {
    stop("Cumulative convergence checkpoints are not exactly 100,200,...,n")
  }
  status <- data.frame(
    n_maps=n,
    status=if (n == 100L) "NEEDS_EXTENSION_FOR_CHANGE_CHECK" else if (overall_stable) "STABLE" else if (n < 500L) "EXTEND" else "INCOMPLETE_AT_500",
    overall_stable=overall_stable,
    conservation_pass=all(conservation$pass), native_totals_crosscheck_pass=all(native_x$totals_match),
    attempted=attempt_audit$summary$attempted,
    successful=attempt_audit$summary$successful,
    failed=attempt_audit$summary$failed,
    seed_schedule_pass=attempt_audit$summary$seed_schedule_pass,
    maxtries_pass=attempt_audit$summary$maxtries_pass,
    attempt_ledger_pass=attempt_audit$summary$all_pass,
    native_area_e_counter_used=FALSE,
    interval_definition="q025/q975 are empirical map-to-map quantiles, not confidence intervals for the mean",
    generated=timestamp(), stringsAsFactors=FALSE
  )

  write_tsv(per_map, file.path(tmp_dir, "per_map_counts.tsv"))
  write_tsv(period, file.path(tmp_dir, "per_map_period_counts.tsv"))
  write_tsv(routes, file.path(tmp_dir, "dispersal_routes_by_map.tsv"))
  write_tsv(ext, file.path(tmp_dir, "extinction_by_area_by_map.tsv"))
  write_tsv(routes_period, file.path(tmp_dir, "dispersal_routes_by_map_period.tsv"))
  write_tsv(ext_period, file.path(tmp_dir, "extinction_by_area_by_map_period.tsv"))
  write_tsv(totals_summary, file.path(tmp_dir, "event_totals_summary.tsv"))
  write_tsv(temporal_summary, file.path(tmp_dir, "temporal_event_rate_summary.tsv"))
  write_tsv(route_summary, file.path(tmp_dir, "d_route_summary.tsv"))
  write_tsv(ext_summary, file.path(tmp_dir, "e_area_summary.tsv"))
  write_tsv(native_x, file.path(tmp_dir, "native_count_crosscheck.tsv"))
  write_tsv(conservation, file.path(tmp_dir, "conservation_checks.tsv"))
  write_tsv(period_conservation, file.path(tmp_dir, "period_conservation_checks.tsv"))
  write_tsv(attempt_audit$summary, file.path(tmp_dir, "attempt_counts_audit.tsv"))
  write_tsv(attempt_audit$records, file.path(tmp_dir, "attempt_seed_maxtries_audit.tsv"))
  write_tsv(convergence, file.path(tmp_dir, "bsm_convergence.tsv"))
  write_tsv(exposure, file.path(tmp_dir, "epoch_lineage_exposure.tsv"))
  write_tsv(ana_all, file.path(tmp_dir, "anagenetic_events_long.tsv"))
  write_tsv(clado_all, file.path(tmp_dir, "cladogenetic_events_long.tsv"))
  write_tsv(status, file.path(tmp_dir, "summary_status.tsv"))
  if (!dir.create(dirname(final_dir), recursive=TRUE, showWarnings=FALSE) && !dir.exists(dirname(final_dir))) {
    stop("Cannot create summaries directory")
  }
  if (!file.rename(tmp_dir, final_dir)) stop("Atomic summary-directory rename failed")
  status
}

static_manifest <- list(
  schema_version="1.0",
  analysis="ScenarioA417_M1_time_stratified_DEC_BSM",
  accepted_optimizer_fit_rds=fit_rds,
  accepted_optimizer_fit_sha256=fit_hash,
  postfit_ancestral_states_rds=postfit_rds,
  postfit_ancestral_states_sha256=postfit_hash,
  optimizer_lnL=ll_optimizer,
  postfit_lnL=ll_postfit,
  postfit_minus_optimizer_lnL=ll_postfit - ll_optimizer,
  final_qa_json=final_qa_json,
  final_qa_sha256=final_qa_hash,
  final_qa_identity=list(schema_version="1.0", status="PASS", scope="model",
                         model_id="M1_R4", model_status="PASS", outcome="COMPLETE"),
  validator_owned_artifacts=qa_bound_artifacts,
  acceptance_tsv=acceptance_tsv,
  acceptance_sha256=qa_acceptance_hash,
  runner=script_path,
  runner_sha256=sha256_file(script_path),
  seed_base=seed_base,
  map_seed_rule="wrap_int32(seed_base + attempt_id)",
  source_seed_rule="wrap_int32(seed_base + 1000000 + attempt_id)",
  maxtries_per_branch=maxtries,
  initial_successful_maps=100L,
  maximum_successful_maps=500L,
  checkpoint_every_successful_maps=100L,
  warning_policy="any warning rejects that attempted map",
  source_direction="one source area sampled from ancestral range using epoch-specific dmat weights and a recorded independent source_seed",
  temporal_rate_denominator="fixed-tree lineage-million-years intersecting each [young,old) epoch",
  native_total_crosscheck="BioGeoBEARS::count_events_huge_tables",
  native_extinction_by_area="BYPASSED_KNOWN_BUG",
  files=lapply(seq_along(expected_paths), function(i) list(
    role=names(expected_paths)[[i]], path=expected_paths[[i]],
    sha256=expected_hashes[[names(expected_paths)[[i]]]]
  )),
  area_codes=as.list(setNames(actual_names, internal_codes)),
  epochs=lapply(seq_along(epoch_labels), function(i) list(
    label=epoch_labels[[i]], young_ma=epoch_bounds[[i]], old_ma=epoch_bounds[[i + 1L]]
  )),
  package=list(BioGeoBEARS="1.1.3", RemoteSha=remote_sha,
               ape=as.character(packageVersion("ape")), digest=as.character(packageVersion("digest")),
               jsonlite=as.character(packageVersion("jsonlite"))),
  BSM_function_hashes=as.list(function_hashes),
  native_e_bug_probe=list(true_e=1L, observed_native_e=native_probe_e,
                          function_name="count_ana_dispersal_events", action="bypass for e-by-area"),
  stability=list(
    precision="MCSE/mean <= 0.05 when mean >=1; otherwise MCSE <=0.05",
    cumulative_change="<=5% when either mean >=1; otherwise absolute change <=0.05",
    ranks="top-3 directed d routes and top-3 extinction areas must be identical as sets",
    comparisons="current vs previous 100-map checkpoint and current vs initial 100",
    at_500="failure is reported as INCOMPLETE_AT_500, never silently accepted"
  )
)

prepare_outdir <- function() {
  if (dir.exists(outdir) && length(list.files(outdir, all.files=TRUE, no..=TRUE)) > 0L) {
    if (!resume) stop("Refusing non-empty outdir without --resume true")
    mf <- file.path(outdir, "bsm_run_manifest.json")
    if (!file.exists(mf)) stop("Resume directory lacks bsm_run_manifest.json")
    old <- jsonlite::read_json(mf, simplifyVector=TRUE)
    for (nm in c("accepted_optimizer_fit_sha256", "postfit_ancestral_states_sha256",
                 "final_qa_sha256", "acceptance_sha256", "runner_sha256",
                 "seed_base", "maxtries_per_branch")) {
      if (!identical(as.character(old[[nm]]), as.character(static_manifest[[nm]]))) {
        stop("Resume manifest mismatch: ", nm)
      }
    }
  } else {
    dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
    for (d in c("maps", "attempts", "logs", "summaries", "locks", "invocations")) {
      dir.create(file.path(outdir, d), showWarnings=FALSE)
    }
    atomic_write_json_new(static_manifest, file.path(outdir, "bsm_run_manifest.json"))
    manifest_tsv <- data.frame(
      key=c("accepted_optimizer_fit_sha256", "postfit_ancestral_states_sha256",
            "final_qa_sha256", "acceptance_sha256", "runner_sha256",
            "seed_base", "maxtries_per_branch",
            "initial_successful_maps", "maximum_successful_maps", "area_order_sha256"),
      value=c(fit_hash, postfit_hash, final_qa_hash, qa_acceptance_hash,
              sha256_file(script_path), seed_base, maxtries,
              100, 500, expected_hashes[["area_order"]]), stringsAsFactors=FALSE
    )
    atomic_write_tsv_new(manifest_tsv, file.path(outdir, "bsm_run_manifest.tsv"))
    atomic_write_tsv_new(preflight, file.path(outdir, "preflight.tsv"))
    atomic_write_tsv_new(exposure, file.path(outdir, "epoch_lineage_exposure.tsv"))
  }
}

attempt_ids_existing <- function(maps=NULL) {
  ids <- integer()
  if (is.null(maps)) maps <- load_maps()
  if (length(maps) > 0L) ids <- c(ids, vapply(maps, function(x) as.integer(x$metadata$attempt_id), integer(1)))
  ff <- list.files(file.path(outdir, "attempts"), pattern="^attempt_[0-9]{6}_FAILED\\.tsv$", full.names=TRUE)
  if (length(ff) > 0L) ids <- c(ids, as.integer(sub("^attempt_([0-9]{6})_FAILED\\.tsv$", "\\1", basename(ff))))
  if (anyDuplicated(ids)) stop("Duplicate attempt IDs in checkpoints")
  ids
}

validate_attempt_records <- function(maps=NULL, through_attempt=NULL) {
  if (is.null(maps)) maps <- load_maps()
  attempt_dir <- file.path(outdir, "attempts")
  fs <- sort(list.files(
    attempt_dir,
    pattern="^attempt_[0-9]{6}_(SUCCESS|FAILED)\\.tsv$",
    full.names=TRUE
  ))
  if (length(fs) == 0L) {
    if (length(maps) > 0L) stop("Successful maps exist without attempt records")
    empty_records <- data.frame(
      attempt_id=integer(), status=character(), map_id=integer(), map_seed=integer(),
      source_seed=integer(), maxtries=integer(), warnings=integer(), error=character(),
      log=character(), completed=character(), record_file=character(),
      filename_attempt_id=integer(), filename_status=character(),
      expected_map_seed=integer(), expected_source_seed=integer(),
      filename_pass=logical(), seed_pass=logical(), maxtries_pass=logical(),
      log_exists=logical(), map_record_pass=logical(), stringsAsFactors=FALSE
    )
    empty_summary <- data.frame(
      through_attempt=0L, attempted=0L, successful=0L, failed=0L,
      configured_maxtries=maxtries, attempt_ids_contiguous=TRUE,
      count_balance_pass=TRUE,
      map_ids_contiguous=TRUE, filename_pass=TRUE, seed_schedule_pass=TRUE,
      maxtries_pass=TRUE, logs_exist=TRUE, map_record_pass=TRUE,
      all_pass=TRUE, stringsAsFactors=FALSE
    )
    return(list(summary=empty_summary, records=empty_records))
  }

  records <- lapply(fs, function(f) {
    z <- read.delim(f, stringsAsFactors=FALSE, check.names=FALSE,
                    quote="", comment.char="", fill=TRUE)
    if (nrow(z) != 1L) stop("Attempt record must contain exactly one row: ", f)
    z$record_file <- normalizePath(f, mustWork=TRUE)
    z
  })
  required_attempt_cols <- c(
    "attempt_id", "status", "map_id", "map_seed", "source_seed", "maxtries",
    "warnings", "error", "log", "completed", "record_file"
  )
  if (any(!vapply(records, function(z) all(required_attempt_cols %in% names(z)), logical(1)))) {
    stop("Attempt record schema mismatch")
  }
  records <- do.call(rbind, lapply(records, function(z) z[, required_attempt_cols, drop=FALSE]))
  records$attempt_id <- as.integer(records$attempt_id)
  records$map_id <- suppressWarnings(as.integer(records$map_id))
  records$map_seed <- as.integer(records$map_seed)
  records$source_seed <- as.integer(records$source_seed)
  records$maxtries <- as.integer(records$maxtries)
  records$warnings <- as.integer(records$warnings)
  bn <- basename(records$record_file)
  records$filename_attempt_id <- as.integer(sub(
    "^attempt_([0-9]{6})_(SUCCESS|FAILED)\\.tsv$", "\\1", bn
  ))
  records$filename_status <- sub(
    "^attempt_([0-9]{6})_(SUCCESS|FAILED)\\.tsv$", "\\2", bn
  )

  if (is.null(through_attempt)) {
    through_attempt <- if (length(maps) > 0L) {
      max(vapply(maps, function(x) as.integer(x$metadata$attempt_id), integer(1)))
    } else {
      max(records$attempt_id)
    }
  }
  through_attempt <- as.integer(through_attempt)
  records <- records[records$attempt_id <= through_attempt, , drop=FALSE]
  records <- records[order(records$attempt_id), , drop=FALSE]
  if (nrow(records) == 0L && through_attempt > 0L) stop("No attempt records through requested cutoff")
  if (anyNA(records$attempt_id) || anyDuplicated(records$attempt_id) ||
      !identical(records$attempt_id, seq_len(through_attempt))) {
    stop("Attempt ledger IDs must be unique and exactly contiguous from 1 through cutoff")
  }
  if (any(!(records$status %in% c("SUCCESS", "FAILED")))) stop("Unknown attempt status")

  records$expected_map_seed <- vapply(
    records$attempt_id, function(i) make_seed(seed_base, i), integer(1)
  )
  records$expected_source_seed <- vapply(
    records$attempt_id, function(i) make_seed(seed_base, 1000000 + i), integer(1)
  )
  records$filename_pass <- records$filename_attempt_id == records$attempt_id &
    records$filename_status == records$status
  records$seed_pass <- records$map_seed == records$expected_map_seed &
    records$source_seed == records$expected_source_seed
  records$maxtries_pass <- !is.na(records$maxtries) & records$maxtries == maxtries
  records$log_exists <- !is.na(records$log) & file.exists(records$log)
  records$map_record_pass <- FALSE

  map_attempt_ids <- if (length(maps) == 0L) integer() else {
    vapply(maps, function(x) as.integer(x$metadata$attempt_id), integer(1))
  }
  success_rows <- which(records$status == "SUCCESS")
  failure_rows <- which(records$status == "FAILED")
  if (length(success_rows) != length(maps) || !setequal(records$attempt_id[success_rows], map_attempt_ids)) {
    stop("SUCCESS records and saved maps are not one-to-one through the audit cutoff")
  }
  if (length(success_rows) > 0L) {
    for (i in success_rows) {
      map_pos <- match(records$attempt_id[[i]], map_attempt_ids)
      x <- maps[[map_pos]]
      records$map_record_pass[[i]] <-
        !is.na(records$map_id[[i]]) && records$map_id[[i]] == as.integer(x$metadata$map_id) &&
        records$map_seed[[i]] == as.integer(x$metadata$map_seed) &&
        records$source_seed[[i]] == as.integer(x$metadata$source_seed) &&
        records$maxtries[[i]] == as.integer(x$metadata$maxtries) &&
        records$warnings[[i]] == 0L
    }
  }
  if (length(failure_rows) > 0L) {
    records$map_record_pass[failure_rows] <- is.na(records$map_id[failure_rows])
  }
  map_ids_contiguous <- identical(records$map_id[success_rows], seq_along(maps))
  count_balance_pass <- length(success_rows) + length(failure_rows) == nrow(records) &&
    nrow(records) == through_attempt
  checks <- c(
    count_balance_pass=count_balance_pass,
    filename_pass=all(records$filename_pass),
    seed_schedule_pass=all(records$seed_pass),
    maxtries_pass=all(records$maxtries_pass),
    logs_exist=all(records$log_exists),
    map_record_pass=all(records$map_record_pass),
    map_ids_contiguous=map_ids_contiguous
  )
  if (!all(checks)) {
    stop("Attempt ledger validation failed: ", paste(names(checks)[!checks], collapse=", "))
  }
  summary <- data.frame(
    through_attempt=through_attempt,
    attempted=nrow(records),
    successful=length(success_rows),
    failed=length(failure_rows),
    configured_maxtries=maxtries,
    attempt_ids_contiguous=TRUE,
    count_balance_pass=count_balance_pass,
    map_ids_contiguous=map_ids_contiguous,
    filename_pass=checks[["filename_pass"]],
    seed_schedule_pass=checks[["seed_schedule_pass"]],
    maxtries_pass=checks[["maxtries_pass"]],
    logs_exist=checks[["logs_exist"]],
    map_record_pass=checks[["map_record_pass"]],
    all_pass=all(checks), stringsAsFactors=FALSE
  )
  list(summary=summary, records=records)
}

attempt_with_log <- function(map_id, attempt_id, bsm_inputs) {
  log_final <- file.path(outdir, "logs", sprintf("attempt_%06d.log", attempt_id))
  if (file.exists(log_final)) stop("Attempt log already exists: ", log_final)
  log_tmp <- paste0(log_final, ".tmp.", Sys.getpid())
  con <- file(log_tmp, open="wt")
  sink(con, type="output")
  sink(con, type="message")
  warnings_seen <- character()
  value <- NULL
  error <- NULL
  tryCatch({
    value <- withCallingHandlers(
      build_one_map(map_id=map_id, attempt_id=attempt_id, bsm_inputs=bsm_inputs),
      warning=function(w) {
        warnings_seen <<- c(warnings_seen, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
  }, error=function(e) error <<- conditionMessage(e))
  sink(type="message")
  sink(type="output")
  close(con)
  if (!file.rename(log_tmp, log_final)) stop("Could not finalize attempt log")
  if (length(warnings_seen) > 0L && is.null(error)) {
    error <- paste("Warnings reject map:", paste(unique(warnings_seen), collapse=" | "))
  }
  list(value=value, error=error, warnings=warnings_seen, log=log_final)
}

with_lock <- function(expr) {
  active <- file.path(outdir, "locks", "ACTIVE")
  if (!dir.create(active, showWarnings=FALSE)) stop("An ACTIVE lock exists; refuse concurrent or ambiguous resume")
  lock_info <- data.frame(pid=Sys.getpid(), host=Sys.info()[["nodename"]], started=timestamp(),
                          action=action, target_maps=target_maps, stringsAsFactors=FALSE)
  write_tsv(lock_info, file.path(active, "lock.tsv"))
  on.exit({
    released <- file.path(outdir, "locks", paste0("RELEASED_", stamp_id(), "_", Sys.getpid()))
    file.rename(active, released)
  }, add=TRUE)
  force(expr)
}

prepare_outdir()

if (identical(action, "summarize")) {
  all_maps <- load_maps()
  attempt_audit <- validate_attempt_records(all_maps)
  nmaps <- length(all_maps)
  if (nmaps < 100L) stop("At least 100 successful maps are needed for a formal checkpoint summary")
  for (n in seq(100L, min(500L, (nmaps %/% 100L) * 100L), 100L)) summarize_checkpoint(n)
  cat("BSM_SUMMARY_PASS maps=", nmaps,
      " attempts=", attempt_audit$summary$attempted,
      " failed=", attempt_audit$summary$failed, "\n", sep="")
  quit(status=0L)
}

with_lock({
  invocation <- data.frame(
    pid=Sys.getpid(), started=timestamp(), action=action, target_maps=target_maps,
    max_attempts=max_attempts, maxtries=maxtries, seed_base=seed_base,
    resume=resume, stop_when_stable=stop_when_stable, stringsAsFactors=FALSE
  )
  inv_path <- file.path(outdir, "invocations", paste0("invocation_", stamp_id(), "_", Sys.getpid(), ".tsv"))
  atomic_write_tsv_new(invocation, inv_path)

  bsm_inputs_fn <- file.path(outdir, "stochastic_mapping_inputs.rds")
  if (file.exists(bsm_inputs_fn)) {
    bsm_inputs <- readRDS(bsm_inputs_fn)
  } else {
    # Explicit NULL and one core avoid the package's FALSE-as-cluster defect.
    res$inputs$num_cores_to_use <- 1L
    options(mc.cores=1L)
    bsm_inputs <- bgb_fun("get_inputs_for_stochastic_mapping_stratified")(
      res=res, cluster_already_open=NULL, rootedge=FALSE,
      statenum_bottom_root_branch_1based=NULL, printlevel=1,
      min_branchlength=1e-6
    )
    if (!is.list(bsm_inputs) || length(bsm_inputs) != length(expected_times) ||
        any(!vapply(bsm_inputs, function(x) identical(x$stratified, TRUE), logical(1)))) {
      stop("Time-stratified BSM input builder did not return exactly five strata")
    }
    atomic_save_rds_new(bsm_inputs, bsm_inputs_fn)
  }

  maps <- load_maps()
  invisible(validate_attempt_records(maps))
  ids <- attempt_ids_existing(maps)
  next_attempt <- if (length(ids) == 0L) 1L else max(ids) + 1L
  if (length(maps) > target_maps) stop("Existing maps exceed requested target; use a target >= existing count")
  resumed_stable <- FALSE
  if (length(maps) >= 100L) {
    for (completed_n in seq(100L, (length(maps) %/% 100L) * 100L, 100L)) {
      completed_status <- summarize_checkpoint(completed_n)
      resumed_stable <- resumed_stable || isTRUE(completed_status$overall_stable[[1]])
    }
  }

  while (length(maps) < target_maps && next_attempt <= max_attempts &&
         !(stop_when_stable && resumed_stable)) {
    map_id <- length(maps) + 1L
    attempt_id <- next_attempt
    cat(sprintf("BSM_ATTEMPT attempt=%d map_goal=%d/%d seed=%d\n",
                attempt_id, map_id, target_maps, make_seed(seed_base, attempt_id)))
    ans <- attempt_with_log(map_id, attempt_id, bsm_inputs)
    if (is.null(ans$error)) {
      map_path <- file.path(outdir, "maps", sprintf("map_%04d.rds", map_id))
      atomic_save_rds_new(ans$value, map_path)
      ok <- data.frame(
        attempt_id=attempt_id, status="SUCCESS", map_id=map_id,
        map_seed=ans$value$metadata$map_seed, source_seed=ans$value$metadata$source_seed,
        maxtries=maxtries, warnings=0L, error="", log=ans$log, completed=timestamp(),
        stringsAsFactors=FALSE
      )
      atomic_write_tsv_new(ok, file.path(outdir, "attempts", sprintf("attempt_%06d_SUCCESS.tsv", attempt_id)))
      maps[[map_id]] <- ans$value
      cat(sprintf("BSM_SUCCESS map=%d attempt=%d d=%d e=%d j=%d\n", map_id, attempt_id,
                  ans$value$metadata$independent_counts$d,
                  ans$value$metadata$independent_counts$e,
                  ans$value$metadata$independent_counts$j))
      if ((map_id %% 100L) == 0L) {
        checkpoint_status <- summarize_checkpoint(map_id)
        cat("BSM_CHECKPOINT n=", map_id, " status=", checkpoint_status$status[[1]], "\n", sep="")
        if (stop_when_stable && map_id >= 200L && isTRUE(checkpoint_status$overall_stable[[1]])) break
      }
    } else {
      fail <- data.frame(
        attempt_id=attempt_id, status="FAILED", map_id=NA_integer_,
        map_seed=make_seed(seed_base, attempt_id),
        source_seed=make_seed(seed_base, 1000000 + attempt_id),
        maxtries=maxtries, warnings=length(ans$warnings), error=clean_text(ans$error),
        log=ans$log, completed=timestamp(), stringsAsFactors=FALSE
      )
      atomic_write_tsv_new(fail, file.path(outdir, "attempts", sprintf("attempt_%06d_FAILED.tsv", attempt_id)))
      cat("BSM_FAILURE attempt=", attempt_id, " error=", clean_text(ans$error), "\n", sep="")
    }
    next_attempt <- next_attempt + 1L
  }

  maps <- load_maps()
  attempt_audit <- validate_attempt_records(maps)
  success_n <- as.integer(attempt_audit$summary$successful)
  failed_n <- as.integer(attempt_audit$summary$failed)
  attempted_n <- as.integer(attempt_audit$summary$attempted)
  ledger <- attempt_audit$records
  ledger_path <- file.path(outdir, "invocations", paste0("attempt_ledger_through_", sprintf("%06d", attempted_n), "_", stamp_id(), ".tsv"))
  atomic_write_tsv_new(ledger, ledger_path)

  if (success_n >= 100L) {
    last_checkpoint <- (success_n %/% 100L) * 100L
    checkpoint_status <- summarize_checkpoint(last_checkpoint)
  } else {
    checkpoint_status <- data.frame(status="INSUFFICIENT_SUCCESSFUL_MAPS", overall_stable=FALSE)
  }
  terminal_status <- if (isTRUE(checkpoint_status$overall_stable[[1]])) {
    "RUN_STABLE"
  } else if (success_n >= 500L) {
    "INCOMPLETE_AT_500"
  } else if (success_n >= target_maps) {
    "RUN_TARGET_REACHED_NEEDS_EXTENSION"
  } else {
    "FAILED_TO_REACH_TARGET"
  }
  final_status <- data.frame(
    status=terminal_status,
    attempted=attempted_n, successful=success_n, failed=failed_n,
    requested_target=target_maps, max_attempts=max_attempts,
    attempt_ids_contiguous=attempt_audit$summary$attempt_ids_contiguous,
    count_balance_pass=attempt_audit$summary$count_balance_pass,
    seed_schedule_pass=attempt_audit$summary$seed_schedule_pass,
    maxtries_pass=attempt_audit$summary$maxtries_pass,
    attempt_ledger_pass=attempt_audit$summary$all_pass,
    latest_checkpoint_status=checkpoint_status$status[[1]],
    latest_checkpoint_stable=isTRUE(checkpoint_status$overall_stable[[1]]),
    finished=timestamp(), stringsAsFactors=FALSE
  )
  status_path <- file.path(outdir, "invocations", paste0("run_status_", stamp_id(), "_", Sys.getpid(), ".tsv"))
  atomic_write_tsv_new(final_status, status_path)
  print(final_status, row.names=FALSE)
  if (final_status$status[[1]] %in% c("FAILED_TO_REACH_TARGET", "INCOMPLETE_AT_500")) {
    stop(if (identical(final_status$status[[1]], "INCOMPLETE_AT_500")) {
      "Reached 500 successful maps without satisfying the frozen Monte Carlo stability rules"
    } else {
      "Failed to obtain the requested number of accepted maps within max-attempts"
    })
  }
})
