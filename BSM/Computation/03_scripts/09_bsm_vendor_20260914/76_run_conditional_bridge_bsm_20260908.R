#!/usr/bin/env Rscript
# Versioned overlay of the frozen 69 runner. No installed namespace is mutated.
suppressPackageStartupMessages({library(digest); library(jsonlite)})
overlay_self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE))
stopifnot(length(overlay_self)==1L)
overlay_self <- normalizePath(overlay_self, mustWork=TRUE)
overlay_job <- dirname(dirname(overlay_self))
overlay_base <- file.path(overlay_job,"03_scripts/69_run_conditional_native_bsm_20260908.R")
overlay_sha <- function(path) digest(file=path,algo="sha256",serialize=FALSE)
stopifnot(overlay_sha(overlay_base)=="30f69a8f322bebdfca1cc18d45750641de434eb51f502100b98c62ea96f7512e")
overlay_injected <- character()

overlay_setup <- function() {
  stopifnot(identical(contract$bsm_path_sampler,"UNIFORMIZATION_NATIVE_EFFECTIVE_Q"),
    identical(contract$r5_execution_mode,"NOT_PERFORMED_BY_USER_DECISION"),
    isTRUE(contract$bridge_relative_tail_tolerance==1e-12),
    isTRUE(contract$bridge_expected_branch_calls==1473))
  helper <- file.path(job_root,"03_scripts/74_ctmc_uniformization_bridge_20260908.R")
  expected <- contract$frozen_helpers[["03_scripts/74_ctmc_uniformization_bridge_20260908.R"]]
  if(is.null(expected) || !identical(overlay_sha(helper),expected)) stop("Bridge helper pin mismatch")
  for(item in contract$bridge_required_evidence) {
    path <- file.path(job_root,item$path)
    if(!identical(overlay_sha(path),item$sha256)) stop("Bridge evidence pin mismatch: ",item$path)
  }
  source(helper,local=.GlobalEnv)
  bridge_context <<- new_bsm_bridge_context(rel_tol=1e-12,max_vector_bytes=2147483648)
  bridge_native_bgb_fun <<- bgb_fun
  bridge_private <<- new.env(parent=asNamespace("BioGeoBEARS"))
  for(n in c("stochastic_mapping_on_stratified","stochastic_map_given_inputs")) {
    fun <- bridge_native_bgb_fun(n)
    environment(fun) <- bridge_private
    assign(n,fun,envir=bridge_private)
  }
  bridge_private$stochastic_map_branch <- bridge_context$branch
  bgb_fun <<- function(name) {
    if(exists(name,envir=bridge_private,inherits=FALSE)) get(name,envir=bridge_private,inherits=FALSE)
    else bridge_native_bgb_fun(name)
  }
  # Substitute the old function body inside a block: its bgb_fun lookup still
  # resolves through run_attempt's capture environment, preserving raw maps.
  modified_builder <- build_one_map
  body(modified_builder) <- substitute({
    bridge_context$reset_audit()
    on.exit({
      bridge_audit_partial <- bridge_context$get_audit()
      if(nrow(bridge_audit_partial)>0L)
        atomic_write_tsv_new(bridge_audit_partial,"bridge_branch_audit.tsv")
    },add=TRUE)
    bridge_payload <- ORIGINAL_BODY
    bridge_audit <- bridge_context$get_audit()
    if(nrow(bridge_audit)!=contract$bridge_expected_branch_calls)
      stop("Bridge branch-call count mismatch: ",nrow(bridge_audit))
    if(any(!is.finite(bridge_audit$relative_tail_bound)) ||
       any(bridge_audit$relative_tail_bound>contract$bridge_relative_tail_tolerance))
      stop("Bridge numerical tail gate failed")
    bridge_payload$bridge_branch_audit <- bridge_audit
    bridge_payload$metadata$branch_sampler <- "UNIFORMIZATION_NATIVE_EFFECTIVE_Q"
    bridge_payload$metadata$bridge_helper_sha256 <- contract$frozen_helpers[["03_scripts/74_ctmc_uniformization_bridge_20260908.R"]]
    bridge_payload$metadata$bridge_tail_union_bound <- sum(bridge_audit$relative_tail_bound)
    bridge_payload$metadata$RNGkind <- RNGkind()
    bridge_payload
  },list(ORIGINAL_BODY=body(build_one_map)))
  build_one_map <<- modified_builder
  helper_hashes <<- setNames(vapply(helper_names,function(n)
    sha256_object(list(formals=formals(get(n)),body=body(get(n)))),character(1)),helper_names)
}

overlay_reuse_preparation <- function() {
  if(action %in% c("preflight","selftest")) return(invisible(NULL))
  cache <- file.path(outdir,"stochastic_mapping_inputs.rds")
  identity <- file.path(outdir,"stochastic_mapping_inputs_identity.json")
  if(file.exists(cache)) return(invisible(NULL)) # 69 independently validates resume identity.
  if(resume) stop("Bridge resume cannot silently create a missing preparation cache")
  src <- file.path(job_root,contract$bridge_native_cache$path)
  src_id <- file.path(job_root,contract$bridge_native_cache$identity_path)
  stopifnot(identical(overlay_sha(src),contract$bridge_native_cache$sha256),
            identical(overlay_sha(src_id),contract$bridge_native_cache$identity_sha256))
  old_id <- read_json(src_id,simplifyVector=TRUE)
  stopifnot(identical(old_id$sha256,contract$bridge_native_cache$sha256),
            identical(old_id$technical_validation_sha256,contract$bridge_native_cache$original_technical_validation_sha256))
  # Copy only deterministic native preparation, never failed candidate histories.
  if(!file.copy(src,cache,overwrite=FALSE)) stop("Cannot copy pinned native preparation")
  if(!identical(overlay_sha(cache),old_id$sha256)) stop("Copied cache digest mismatch")
  atomic_write_json_new(list(sha256=old_id$sha256,
    technical_validation_sha256=technical_validation_hash,
    cache_origin="BYTE_IDENTICAL_NATIVE_PREPARATION_REUSED_NO_ENDPOINT_RECALCULATION",
    source_cache=src,source_identity_sha256=overlay_sha(src_id),
    source_technical_validation_sha256=old_id$technical_validation_sha256,
    source_preparation_elapsed_seconds=old_id$elapsed_seconds,
    object_bytes=old_id$object_bytes),identity)
}

for(overlay_expr in parse(overlay_base)) {
  overlay_lhs <- if(is.call(overlay_expr) && identical(overlay_expr[[1]],as.name("<-")) &&
    is.symbol(overlay_expr[[2]])) as.character(overlay_expr[[2]]) else ""
  if(identical(overlay_lhs,"active")) {
    overlay_reuse_preparation()
    overlay_injected <- c(overlay_injected,"cache")
  }
  eval(overlay_expr,envir=.GlobalEnv)
  if(identical(overlay_lhs,"contract")) {
    overlay_setup()
    overlay_injected <- c(overlay_injected,"setup")
  }
  if(identical(overlay_lhs,"technical_validation")) {
    technical_validation$branch_sampler <- contract$bsm_path_sampler
    technical_validation$api_hashes_scope <- "Runtime bindings: private branch kernel replaced; native dispatcher bodies and other functions unchanged"
    technical_validation$bridge_helper_sha256 <- contract$frozen_helpers[["03_scripts/74_ctmc_uniformization_bridge_20260908.R"]]
    technical_validation$generator_interpretation <- "Qeff exactly represents native waiting rates and normalized outgoing probabilities; Qraw and endpoint preparation remain unchanged"
    overlay_injected <- c(overlay_injected,"technical")
  }
  if(identical(overlay_lhs,"static_manifest")) {
    static_manifest$branch_sampler <- contract$bsm_path_sampler
    static_manifest$api_hashes_scope <- "Runtime bindings, not an assertion that every binding is the unmodified installed function"
    static_manifest$branch_sampler_reference <- "https://doi.org/10.1214/09-AOAS247"
    static_manifest$bridge_relative_tail_tolerance <- contract$bridge_relative_tail_tolerance
    static_manifest$bridge_expected_branch_calls <- contract$bridge_expected_branch_calls
    static_manifest$bridge_helper_sha256 <- contract$frozen_helpers[["03_scripts/74_ctmc_uniformization_bridge_20260908.R"]]
    static_manifest$base_runner_sha256 <- overlay_sha(overlay_base)
    static_manifest$native_preparation_reuse <- contract$bridge_native_cache
    static_manifest$maxtries_interpretation <- "Legacy ledger field retained; uniformization has no rejection loop or manual fallback; see bridge term and tail diagnostics"
    static_manifest$forced_history_policy <- "No manual fallback exists in the replacement branch kernel; original rejection markers remain rejection gates"
    static_manifest$local_bindings <- c("stochastic_mapping_on_stratified","stochastic_map_given_inputs","stochastic_map_branch")
    overlay_injected <- c(overlay_injected,"manifest")
    stopifnot(setequal(overlay_injected,c("setup","technical","manifest")))
  }
}
