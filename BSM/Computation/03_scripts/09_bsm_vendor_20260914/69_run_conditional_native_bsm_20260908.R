#!/usr/bin/env Rscript
# Conditional native M1 BSM continuation. Never emits formal model acceptance.
# Reviewed counting/stability routines are imported by exact source hash only;
# the old runner's authorization gates and top-level execution are not executed.
suppressPackageStartupMessages({
  library(ape); library(BioGeoBEARS); library(methods); library(digest); library(jsonlite)
})
invocation_started <- format(Sys.time(),"%Y-%m-%dT%H:%M:%OS6Z",tz="UTC")

args_parse <- function(x) {
  z <- list(); i <- 1L
  while (i <= length(x)) {
    if (!startsWith(x[[i]], "--") || i == length(x)) stop("Use --key value pairs")
    key <- sub("^--", "", x[[i]])
    if (key %in% names(z)) stop("Duplicate argument: ", key)
    z[[key]] <- x[[i+1L]]; i <- i+2L
  }
  z
}
opts <- args_parse(commandArgs(TRUE))
required <- c("action", "contract-json", "contract-sha", "outdir")
if (!all(required %in% names(opts))) stop("Required: ", paste(required, collapse=", "))
if (length(setdiff(names(opts), c(required,"target-maps","max-attempts","resume")))) stop("Unknown argument")
action <- opts$action
if (!(action %in% c("preflight","selftest","pilot","run","summarize"))) stop("Invalid action")
target_maps <- as.integer(if (is.null(opts[["target-maps"]])) if(action=="pilot") 1L else 100L else opts[["target-maps"]])
if (!(target_maps %in% c(1L,seq(100L,500L,100L)))) stop("Invalid target")
if (action=="pilot" && target_maps!=1L) stop("Pilot target must equal one")
if (action=="run" && target_maps==1L) stop("Use action=pilot for one map")
max_attempts <- as.integer(if (is.null(opts[["max-attempts"]])) if(action=="pilot") 10L else min(1000L,2L*target_maps) else opts[["max-attempts"]])
if (!is.finite(max_attempts) || max_attempts<target_maps || max_attempts>1000L) stop("Global attempt ceiling must be target..1000")
resume <- if(is.null(opts$resume)) FALSE else identical(tolower(opts$resume),"true")
if (!is.null(opts$resume) && !(tolower(opts$resume)%in%c("true","false"))) stop("Invalid resume flag")
stop_when_stable <- TRUE
maxtries <- 40000L
seed_base <- 202609080
no_success_failure_stop <- 10L
forced_history_in_log <- function(lines) {
  markers <- c("Error_in_stochastic_simulation; no success on branch","manually devising a history to force-fit",
               "program is in the manually-sorted events section","Attempted manual history sampled disallowed state")
  any(vapply(markers,function(x)any(grepl(x,lines,fixed=TRUE)),logical(1)))
}
trailing_failure_count <- function(statuses) {
  if(anyNA(statuses) || any(!statuses %in% c("SUCCESS","FAILED"))) stop("Invalid attempt status")
  if(!length(statuses)) return(0L)
  as.integer(sum(cumprod(as.integer(rev(statuses)=="FAILED"))))
}

script_arg <- grep("^--file=",commandArgs(FALSE),value=TRUE)
if(length(script_arg)!=1L) stop("Cannot resolve script")
script_path <- normalizePath(sub("^--file=","",script_arg),mustWork=TRUE)
job_root <- normalizePath(file.path(dirname(script_path),".."),mustWork=TRUE)
original_runner <- file.path(job_root,"03_scripts/05_run_bsm_m1_stratified.R")
hash_file <- function(p) digest(file=p,algo="sha256",serialize=FALSE)
if(!identical(hash_file(original_runner),"cf03d8019906827ca3c6fcb49874d210acd4f6983682956c684d7c52a100abc5")) stop("Reviewed helper source changed")

# Rename legacy metadata fields, not algorithmic operations or thresholds.
renames <- c(accepted_optimizer_fit_sha256="source_optimizer_fit_sha256",
             final_qa_sha256="conditional_contract_sha256", acceptance_sha256="technical_validation_sha256",
             final_qa_hash="conditional_contract_hash", qa_acceptance_hash="technical_validation_hash")
rewrite_metadata <- function(x) {
  if(is.symbol(x)) { k<-as.character(x); return(if(k%in%names(renames))as.name(renames[[k]]) else x) }
  if(is.character(x)) { hit<-x%in%names(renames); x[hit]<-unname(renames[x[hit]]); return(x) }
  if(is.call(x) || is.pairlist(x) || is.expression(x)) {
    y<-lapply(as.list(x),rewrite_metadata); n<-names(y)
    if(!is.null(n)) { hit<-n%in%names(renames); n[hit]<-unname(renames[n[hit]]); names(y)<-n }
    return(if(is.call(x))as.call(y) else if(is.pairlist(x))as.pairlist(y) else as.expression(y))
  }
  x
}
helper_names <- c("timestamp","stamp_id","sha256_file","sha256_object","clean_text",
  "atomic_save_rds_new","atomic_write_tsv_new","atomic_write_json_new","write_tsv","bgb_fun",
  "epoch_for_age","tree_exposure","parse_range","standardize_anagenetic","standardize_cladogenetic",
  "native_total_counts","make_seed","build_one_map","map_files","load_maps","summarize_vector",
  "precision_metric","precision_ok","change_metric","change_ok","max_or_zero","top_keys",
  "summary_dir_for","read_checkpoint_means","summarize_checkpoint","validate_attempt_records")
found <- character()
for(expr in parse(original_runner)) {
  if(is.call(expr) && identical(expr[[1]],as.name("<-")) && is.symbol(expr[[2]]) &&
     as.character(expr[[2]])%in%helper_names && is.call(expr[[3]]) && identical(expr[[3]][[1]],as.name("function"))) {
    eval(rewrite_metadata(expr),envir=.GlobalEnv); found<-c(found,as.character(expr[[2]]))
  }
}
if(!setequal(found,helper_names) || anyDuplicated(found)) stop("Helper import mismatch")
helper_hashes <- setNames(vapply(helper_names,function(n)sha256_object(list(formals=formals(get(n)),body=body(get(n)))),character(1)),helper_names)

contract_path <- normalizePath(opts[["contract-json"]],mustWork=TRUE)
conditional_contract_hash <- hash_file(contract_path)
if(!identical(conditional_contract_hash,opts[["contract-sha"]])) stop("Conditional contract hash mismatch")
contract <- read_json(contract_path,simplifyVector=TRUE)
require_contract <- list(schema_version="1.0",status="APPROVED_CONDITIONAL_CONTINUATION",backend="native",
                         scientific_acceptance="NONE",model="M1",max_successful_maps=500,max_attempts=1000,KKT1=FALSE,
                         seed_base=seed_base,maxtries=maxtries,no_success_failure_stop=no_success_failure_stop,
                         first_stability_checkpoint=200,minimum_successful_maps=100,map_batch_size=100,pilot_successful_maps=1)
for(k in names(require_contract)) if(!identical(as.character(contract[[k]]),as.character(require_contract[[k]]))) stop("Conditional contract field mismatch: ",k)
fit_rds <- file.path(job_root,"04_runs/M1_R4_final_interruption_r2_20260907_6c/fit_result.rds")
postfit_rds <- file.path(dirname(fit_rds),"postfit/postfit_recalculated_ancestral_states.rds")
fit_hash <- hash_file(fit_rds); postfit_hash <- hash_file(postfit_rds)
if(!identical(fit_hash,"0d66ac670ef971687cee55a8f4af2a8661eebb99dc4a158d51f3e3598a800ca7") ||
   !identical(postfit_hash,"d635c71472d93424b8825b9b8f613e9b1c42497f762421ed5262b6648a8b1504") ||
   !identical(contract$fit_sha256,fit_hash) || !identical(contract$postfit_sha256,postfit_hash)) stop("Pinned native input identity mismatch")
allowed_bsm_root <- normalizePath(file.path(job_root,"05_bsm"),mustWork=TRUE)
outdir <- file.path(normalizePath(dirname(opts$outdir),mustWork=TRUE),basename(opts$outdir))
if(!identical(dirname(outdir),allowed_bsm_root) || !grepl("^M1_BSM_CONDITIONAL_[A-Za-z0-9_]+$",basename(outdir))) stop("Output must be a new direct conditional BSM child")

expected_paths <- list(tree=file.path(job_root,"01_inputs/frozen/tree_scenarioA417.tre"),
  geog=file.path(job_root,"01_inputs/frozen/geog_scenarioA417_analysis_order.LagrangePHYLIP"),
  times=file.path(job_root,"02_config/timeperiods_5epochs.txt"),
  multipliers=file.path(job_root,"02_config/M1_conservative_dispersal_multipliers.txt"),
  area_order=file.path(job_root,"02_config/area_order.tsv"))
expected_hashes <- c(tree="06725d5a0c1ef29007aeb0f7657175c2d2329c75cb4d8fb2170e39b5fdaf6962",
 geog="89c5b776bdce88d9c0074b3a6117d25e50c8604a69dba158d39fa5544527f1c0",
 times="90dadc7755385b44d9333f22e2acaeac163bd0ce4507289e1720671eccc4e51f",
 multipliers="34cc3ee7f0789d978948e35085fb9ed89de09f2600b50384bb074bd11db7029c",
 area_order="7daa151006fbe79d00ed92a83938aca4b04dbc89bbaa814710bf244dcb1a0804")
for(k in names(expected_paths)) if(!identical(hash_file(expected_paths[[k]]),expected_hashes[[k]])) stop("Frozen input changed: ",k)
if(!identical(as.character(packageVersion("BioGeoBEARS")),"1.1.3") ||
   !identical(unname(packageDescription("BioGeoBEARS")$RemoteSha),"1672cc0c171b1a05effad69fa426427b3d9ef4e1")) stop("Native package identity changed")
api_names <- c("get_inputs_for_stochastic_mapping_stratified","get_inputs_for_stochastic_mapping_from_results_object",
  "stochastic_mapping_on_stratified","stochastic_map_given_inputs","stochastic_map_branch","get_dmat_times_from_res",
  "simulate_source_area_ana","simulate_source_area_clado","get_huge_events_tables_from_clado_ana_events_tables",
  "count_events_huge_tables","count_ana_dispersal_events")
api_hashes <- setNames(vapply(api_names,function(n)sha256_object(list(formals=formals(bgb_fun(n)),body=body(bgb_fun(n)))),character(1)),api_names)
if(!identical(api_hashes[["get_inputs_for_stochastic_mapping_stratified"]],"df3161e2a22ddf1ea330d95691e7a1c529f56f61b2b109f281c1bc4de5b6e97c") ||
   !identical(api_hashes[["stochastic_mapping_on_stratified"]],"54e87de693f0f15b0a256b299e2c7bd85afed0c06b09672975800295e190e428")) stop("Native BSM dispatcher body changed")
fit <- readRDS(fit_rds); res <- readRDS(postfit_rds)
prov <- fit$inputs$provenance
if(!identical(prov$model,"M1") || !isTRUE(prov$final_fit) || as.integer(prov$max_range_size)!=4L ||
   !identical(fit$inputs$speedup,FALSE) || !identical(fit$inputs$force_sparse,FALSE) || !isTRUE(fit$inputs$include_null_range)) stop("Source final recipe mismatch")
pt <- fit$outputs@params_table; pt2 <- res$outputs@params_table
if(!identical(rownames(pt)[pt[,"type"]=="free"],c("d","e")) || !isTRUE(as.numeric(pt["j","est"])==0) ||
   !isTRUE(as.numeric(pt["w","est"])==1) || !isTRUE(all.equal(as.numeric(pt[c("d","e","j","w"),"est"]),as.numeric(pt2[c("d","e","j","w"),"est"]),tolerance=1e-12))) stop("Postfit/source parameter mismatch")
ll_optimizer <- as.numeric(fit$total_loglikelihood); ll_postfit <- as.numeric(res$total_loglikelihood)
if(length(ll_optimizer)!=1L || length(ll_postfit)!=1L || !is.finite(ll_optimizer) || !is.finite(ll_postfit) || abs(ll_optimizer-ll_postfit)>1e-6) stop("Postfit likelihood consistency failed")
if(!identical(res$inputs$postfit_provenance$real_optimizer_fit_sha256,fit_hash) ||
   !identical(res$postfit_completion$real_optimizer_fit_sha256,fit_hash)) stop("Postfit optimizer binding failed")
expected_times <- c(3.3,13.9,23.03,33.9,45)
epoch_bounds <- c(0,expected_times)
epoch_labels <- c("0-3.3","3.3-13.9","13.9-23.03","23.03-33.9","33.9-45")
internal_codes <- LETTERS[1:11]
actual_names <- c("AF","AUS","CAM","ENA","EPA","IND","MDG","ORI","SAM","WNA","WPA")
ao <- read.delim(expected_paths$area_order,stringsAsFactors=FALSE); ao<-ao[order(ao$analysis_bit_index),]
if(!identical(ao$analysis_internal_code,internal_codes) || !identical(ao$abbrev,actual_names)) stop("Area dictionary mismatch")
code_to_area <- setNames(actual_names,internal_codes)
tr <- read.tree(expected_paths$tree)
geog <- getranges_from_LagrangePHYLIP(expected_paths$geog)
if(length(tr$tip.label)!=417L || tr$Nnode!=416L || !setequal(tr$tip.label,rownames(geog@df)) || !identical(colnames(geog@df),internal_codes)) stop("Frozen tree/geography mismatch")
if(!identical(as.numeric(res$inputs$timeperiods),expected_times) || length(res$inputs$tree_sections_list)!=5L ||
   !identical(dim(res$condlikes_table),c(1891L,562L)) ||
   !identical(dim(res$relative_probs_of_each_state_at_branch_bottom_below_node_DOWNPASS_TABLE),c(1891L,562L))) stop("Stratified BSM arrays absent/inconsistent")
posterior <- res$ML_marginal_prob_each_state_at_branch_top_AT_node
if(!identical(dim(posterior),c(833L,562L)) || any(!is.finite(posterior)) || any(posterior<0) || max(abs(rowSums(posterior)-1))>1e-10) stop("Native node probability gate failed")
for(r in list(fit,res)) {
  pp <- list(tree=r$inputs$trfn,geog=r$inputs$geogfn,times=r$inputs$timesfn,multipliers=r$inputs$dispersal_multipliers_fn)
  for(k in names(pp)) if(!identical(normalizePath(pp[[k]],mustWork=TRUE),normalizePath(expected_paths[[k]],mustWork=TRUE))) stop("Native runtime input path mismatch")
}
exposure <- tree_exposure(tr)
technical_validation <- list(status="TECHNICAL_PASS_CONDITIONAL_ONLY",scientific_acceptance="NONE",model="M1",
  fit_sha256=fit_hash,postfit_sha256=postfit_hash,lnL_delta=ll_postfit-ll_optimizer,
  node_probability_row_error=max(abs(rowSums(posterior)-1)),tips=417L,internal_nodes=416L,states=562L,
  time_strata=5L,immutable_source_KKT1=FALSE,native_api_hashes=as.list(api_hashes),frozen_hashes=as.list(expected_hashes))
technical_validation_hash <- sha256_object(technical_validation)
static_manifest <- list(schema_version="conditional-bsm-1.0",scientific_acceptance="NONE",analysis="M1_NATIVE_CONDITIONAL_BSM",
  source_optimizer_fit_sha256=fit_hash,postfit_ancestral_states_sha256=postfit_hash,
  conditional_contract_sha256=conditional_contract_hash,technical_validation_sha256=technical_validation_hash,
  runner_sha256=hash_file(script_path),helper_source_sha256=hash_file(original_runner),helper_hashes=as.list(helper_hashes),
  native_api_hashes=as.list(api_hashes),seed_base=seed_base,maxtries_per_branch=maxtries,maximum_successful_maps=500L,
  global_maximum_attempts=1000L,checkpoint_every_successful_maps=100L,warning_policy="reject every warning",
  no_success_failure_stop=no_success_failure_stop,failure_stop_definition="10 consecutive failed attempts, including across resume",
  forced_history_policy="reject printed maxtries/manual force-fit histories; never alter native namespace",
  temporal_rate_denominator="fixed-tree lineage-million-years in half-open geological epochs",
  source_direction="epoch-specific weighted source imputation with separate seed; not directly identified source",
  probability_status="conditional existing native postfit; native KKT1 remains FALSE",
  precision_rule="original <=0.05 MCSE/change/top-three set criteria; minimum200 maps for stability")
cat("CONDITIONAL_BSM_PREFLIGHT_PASS\n")
cat(toJSON(technical_validation,auto_unbox=TRUE,digits=NA,pretty=TRUE),"\n")
if(action=="preflight") quit(status=0)
if(action=="selftest") {
  rejected <- function(expr) inherits(tryCatch({force(expr); NULL},error=identity),"error")
  global_names<-unique(unlist(lapply(helper_names,function(n)codetools::findGlobals(get(n),merge=TRUE))))
  missing_globals<-global_names[!vapply(global_names,exists,logical(1),envir=.GlobalEnv,inherits=TRUE)]
  # Original summarize_checkpoint uses with(data.frame, ...) for these columns.
  nse_columns<-c("cladogenetic_events","count_d","count_e","count_extinction_area","count_route",
    "d_route_sum","d_total","e_area_sum","e_total","j_total","native_totals_match","period_conservation_pass")
  missing_globals<-setdiff(missing_globals,nse_columns)
  if(length(missing_globals)) stop("Imported helper unresolved globals: ",paste(missing_globals,collapse=", "))
  checks <- c(epoch_boundary=vapply(epoch_bounds[-length(epoch_bounds)],epoch_for_age,integer(1))==1:5,
    upper_epoch_reject=rejected(epoch_for_age(45)),negative_epoch_reject=rejected(epoch_for_age(-1)),
    empty_failures=trailing_failure_count(character())==0L,
    trailing_failures=trailing_failure_count(c("FAILED","SUCCESS","FAILED","FAILED"))==2L,
    reset_failures=trailing_failure_count(c("FAILED","SUCCESS"))==0L,
    cap_failures=trailing_failure_count(rep("FAILED",10))==no_success_failure_stop,
    invalid_failure_reject=rejected(trailing_failure_count("INVALID")),
    clean_log=!forced_history_in_log(c("normal native output","success")),
    manual_log=forced_history_in_log("manually devising a history to force-fit"),
    maxtries_log=forced_history_in_log("Error_in_stochastic_simulation; no success on branch"),
    source_seed_separate=make_seed(seed_base,1)!=make_seed(seed_base,1000001),
    metadata_renames=!any(grepl("accepted_optimizer_fit_sha256|final_qa_hash|qa_acceptance_hash",deparse(body(build_one_map)))),
    precision_pass=precision_ok(100,5),precision_fail=!precision_ok(100,5.01),
    precision_nonfinite_reject=!precision_ok(100,NaN),change_nonfinite_reject=!change_ok(Inf,100))
  fixture <- data.frame(abs_event_time=c(3.3,13.9),event_type=c("d","e"),current_rangetxt=c("A","AB"),
    new_rangetxt=c("AB","A"),ana_dispersal_from=c("A",""),dispersal_to=c("B",""),extirpation_from=c("","B"))
  z<-standardize_anagenetic(fixture,1L)
  checks<-c(checks,event_epoch=identical(z$epoch,epoch_labels[c(2,3)]),area_semantics=identical(z$affected_area,c("AUS","AUS")))
  fixture$new_rangetxt[1]<-"ABC"
  checks<-c(checks,multiple_gain_reject=rejected(standardize_anagenetic(fixture,1L)))
  if(!all(checks)) stop("No-map selftest failed: ",paste(names(checks)[!checks],collapse=", "))
  cat("CONDITIONAL_BSM_SELFTEST_PASS checks=",length(checks)," native_maps_called=0 filesystem_writes=0\n",sep="")
  quit(status=0)
}

threads <- c("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","VECLIB_MAXIMUM_THREADS","NUMEXPR_NUM_THREADS","BLIS_NUM_THREADS")
if(any(Sys.getenv(threads)!="1")) stop("All nested thread limits must equal one")
res$inputs$num_cores_to_use<-1L
options(mc.cores=1L)
if(dir.exists(outdir)) {
  if(!resume) stop("Existing conditional namespace requires explicit resume")
  old <- read_json(file.path(outdir,"bsm_run_manifest.json"),simplifyVector=FALSE)
  for(k in c("source_optimizer_fit_sha256","postfit_ancestral_states_sha256","conditional_contract_sha256",
             "technical_validation_sha256","runner_sha256","seed_base","maxtries_per_branch"))
    if(!identical(as.character(old[[k]]),as.character(static_manifest[[k]]))) stop("Resume manifest mismatch: ",k)
} else {
  dir.create(outdir)
  for(d in c("maps","attempts","commits","work","summaries","locks","invocations","preparation")) dir.create(file.path(outdir,d))
  atomic_write_json_new(static_manifest,file.path(outdir,"bsm_run_manifest.json"))
  atomic_write_json_new(technical_validation,file.path(outdir,"technical_validation.json"))
  atomic_write_tsv_new(exposure,file.path(outdir,"epoch_lineage_exposure.tsv"))
}

copy_view <- function(src,dst) {
  if(file.exists(dst)) { if(!identical(hash_file(src),hash_file(dst))) stop("Immutable view mismatch: ",dst) }
  else if(!file.copy(src,dst,overwrite=FALSE)) stop("Cannot create derived checkpoint view: ",dst)
}
commit_dirs <- function() sort(list.dirs(file.path(outdir,"commits"),full.names=TRUE,recursive=FALSE))
verify_commit <- function(path) {
  seal <- readRDS(file.path(path,"commit_seal.rds"))
  if(!identical(seal$conditional_contract_sha256,conditional_contract_hash) || !identical(seal$runner_sha256,hash_file(script_path))) stop("Commit ownership mismatch")
  for(k in names(seal$files)) if(!identical(hash_file(file.path(path,k)),seal$files[[k]])) stop("Commit payload digest mismatch: ",k)
  record <- read.delim(file.path(path,seal$record_file),stringsAsFactors=FALSE,quote="",comment.char="")
  if(nrow(record)!=1L || as.integer(record$attempt_id)!=seal$attempt_id) stop("Commit record identity mismatch")
  list(seal=seal,record=record)
}
rebuild_views <- function() {
  dirs <- commit_dirs()
  ids <- unname(vapply(dirs,function(p)verify_commit(p)$seal$attempt_id,integer(1)))
  if(!identical(ids,seq_along(ids))) stop("Commit attempt sequence is not contiguous")
  for(p in dirs) {
    z<-verify_commit(p); record<-z$record; id<-as.integer(record$attempt_id)
    copy_view(file.path(p,z$seal$record_file),file.path(outdir,"attempts",sprintf("attempt_%06d_%s.tsv",id,record$status)))
    if(record$status=="SUCCESS") copy_view(file.path(p,"map.rds"),file.path(outdir,"maps",sprintf("map_%04d.rds",as.integer(record$map_id))))
  }
  ids
}
seal_workspace <- function(workspace,record_file,attempt_id) {
  files <- list.files(workspace,full.names=FALSE,recursive=FALSE)
  files <- files[file.info(file.path(workspace,files))$isdir==FALSE & files!="commit_seal.rds"]
  seal <- list(attempt_id=as.integer(attempt_id),conditional_contract_sha256=conditional_contract_hash,
               runner_sha256=hash_file(script_path),record_file=record_file,
               files=setNames(vapply(files,function(f)hash_file(file.path(workspace,f)),character(1)),files))
  atomic_save_rds_new(seal,file.path(workspace,"commit_seal.rds"))
}
commit_workspace <- function(workspace,attempt_id) {
  target <- file.path(outdir,"commits",sprintf("attempt_%06d",attempt_id))
  if(dir.exists(target)) stop("Commit already exists")
  invisible(verify_commit(workspace))
  if(!file.rename(workspace,target)) stop("Atomic attempt commit failed")
}
recover_orphans <- function() {
  workdirs <- sort(list.dirs(file.path(outdir,"work"),full.names=TRUE,recursive=FALSE))
  for(w in workdirs) {
    if(!grepl("^attempt_[0-9]{6}$",basename(w))) stop("Unknown workspace: ",w)
    id <- as.integer(sub("attempt_","",basename(w)))
    if(file.exists(file.path(w,"commit_seal.rds"))) { commit_workspace(w,id); next }
    # An incomplete attempt is retained as FAILED, never silently reused or accepted.
    log <- file.path(w,"stdout.log")
    if(!file.exists(log)) writeLines("Interrupted before stdout log was created",log)
    record <- data.frame(attempt_id=id,status="FAILED",map_id=NA_integer_,map_seed=make_seed(seed_base,id),
      source_seed=make_seed(seed_base,1000000+id),maxtries=maxtries,warnings=0L,
      error="INTERRUPTED_BEFORE_ATOMIC_COMMIT: original workspace retained; no map accepted",
      log=file.path(outdir,"commits",basename(w),"stdout.log"),completed=timestamp(),stringsAsFactors=FALSE)
    atomic_write_tsv_new(record,file.path(w,"recovery_record.tsv"))
    seal_workspace(w,"recovery_record.tsv",id); commit_workspace(w,id)
  }
}

run_attempt <- function(map_id,attempt_id,bsm_inputs) {
  workspace <- file.path(outdir,"work",sprintf("attempt_%06d",attempt_id))
  if(!dir.create(workspace,showWarnings=FALSE)) stop("Attempt workspace already exists")
  wd <- getwd(); setwd(workspace); on.exit(setwd(wd),add=TRUE)
  started <- Sys.time(); warnings_seen<-character(); value<-NULL; error<-NULL
  # Retain the raw native history even if later source-imputation/count checks fail.
  map_env <- new.env(parent=.GlobalEnv)
  map_env$bgb_fun <- function(name) {
    f <- bgb_fun(name)
    if(name!="stochastic_mapping_on_stratified") return(f)
    function(...) { ans<-f(...); atomic_save_rds_new(ans,"raw_native_map.rds"); ans }
  }
  map_builder <- build_one_map; environment(map_builder)<-map_env
  con <- file("stdout.log",open="wt"); sink(con,type="output"); sink(con,type="message")
  tryCatch({value <- withCallingHandlers(map_builder(map_id,attempt_id,bsm_inputs),warning=function(w) {
    warnings_seen<<-c(warnings_seen,conditionMessage(w)); invokeRestart("muffleWarning")
  })},error=function(e)error<<-conditionMessage(e))
  sink(type="message"); sink(type="output"); close(con)
  logs <- readLines("stdout.log",warn=FALSE)
  forced <- forced_history_in_log(logs)
  if(is.null(value) && is.null(error)) error<-"Native map builder returned NULL"
  if(length(warnings_seen)>0L && is.null(error)) error<-paste("Warnings reject map:",paste(unique(warnings_seen),collapse=" | "))
  if(forced) error<-paste("FORCED_MANUAL_HISTORY_REJECTED",if(is.null(error))"" else error)
  if(!is.null(value)) {
    value$metadata$scientific_acceptance <- "NONE"
    value$metadata$conditional_inference <- TRUE
    value$metadata$elapsed_seconds <- as.numeric(difftime(Sys.time(),started,units="secs"))
    atomic_save_rds_new(value,if(is.null(error))"map.rds" else "rejected_candidate_map.rds")
  }
  record <- data.frame(attempt_id=attempt_id,status=if(is.null(error))"SUCCESS" else "FAILED",
    map_id=if(is.null(error))map_id else NA_integer_,map_seed=make_seed(seed_base,attempt_id),
    source_seed=make_seed(seed_base,1000000+attempt_id),maxtries=maxtries,warnings=length(warnings_seen),
    error=if(is.null(error))"" else clean_text(error),
    log=file.path(outdir,"commits",basename(workspace),"stdout.log"),completed=timestamp(),stringsAsFactors=FALSE)
  atomic_write_tsv_new(record,"record.tsv")
  atomic_write_json_new(list(attempt_id=attempt_id,elapsed_seconds=as.numeric(difftime(Sys.time(),started,units="secs")),
                            forced_history_detected=forced,warnings=as.list(warnings_seen)),"timing_and_warnings.json")
  setwd(wd)
  seal_workspace(workspace,"record.tsv",attempt_id); commit_workspace(workspace,attempt_id)
  rebuild_views()
  record
}

active <- file.path(outdir,"locks/ACTIVE")
if(!dir.create(active,showWarnings=FALSE)) stop("ACTIVE lock exists: verify old process externally; never delete or bypass blindly")
write_tsv(data.frame(pid=Sys.getpid(),host=Sys.info()[["nodename"]],started=timestamp()),file.path(active,"lock.tsv"))
on.exit_release <- function() file.rename(active,file.path(outdir,"locks",paste0("RELEASED_",stamp_id(),"_",Sys.getpid())))
run_body <- function() {
  on.exit(on.exit_release(),add=TRUE)
  original_wd<-getwd(); on.exit(setwd(original_wd),add=TRUE)
  recover_orphans(); ids<-rebuild_views(); maps<-load_maps()
  initial_audit<-validate_attempt_records(maps,through_attempt=length(ids))
  consecutive_failures<-trailing_failure_count(initial_audit$records$status)
  if(length(maps)>target_maps && action!="summarize") stop("Existing maps exceed target")
  cache <- file.path(outdir,"stochastic_mapping_inputs.rds")
  cache_id <- file.path(outdir,"stochastic_mapping_inputs_identity.json")
  if(file.exists(cache)) {
    if(!file.exists(cache_id)) stop("Preparation cache exists without identity: review before reuse")
    identity <- read_json(cache_id,simplifyVector=TRUE)
    if(!identical(identity$sha256,hash_file(cache)) || !identical(identity$technical_validation_sha256,technical_validation_hash)) stop("Cache digest/ownership mismatch")
    bsm_inputs<-readRDS(cache)
  } else {
    if(action=="summarize") stop("Summarize cannot create preparation cache")
    wd<-getwd(); setwd(file.path(outdir,"preparation")); t0<-Sys.time()
    res$inputs$num_cores_to_use<-1L; options(mc.cores=1L)
    bsm_inputs<-withCallingHandlers(bgb_fun("get_inputs_for_stochastic_mapping_stratified")(
      res=res,cluster_already_open=NULL,rootedge=FALSE,statenum_bottom_root_branch_1based=NULL,
      printlevel=1,min_branchlength=1e-6),warning=function(w)stop("Preparation warning rejected: ",conditionMessage(w)))
    setwd(wd)
    if(length(bsm_inputs)!=5L || any(!vapply(bsm_inputs,function(x)identical(x$stratified,TRUE),logical(1)))) stop("Prepared strata invalid")
    atomic_save_rds_new(bsm_inputs,cache)
    atomic_write_json_new(list(sha256=hash_file(cache),technical_validation_sha256=technical_validation_hash,
      elapsed_seconds=as.numeric(difftime(Sys.time(),t0,units="secs")),object_bytes=as.numeric(object.size(bsm_inputs))),cache_id)
  }
  if(length(maps)>=100L) for(n in seq(100L,(length(maps)%/%100L)*100L,100L)) summarize_checkpoint(n)
  stable<-FALSE
  if(length(maps)>=200L) stable<-isTRUE(summarize_checkpoint((length(maps)%/%100L)*100L)$overall_stable[[1]])
  if(action!="summarize") while(length(maps)<target_maps && length(ids)<max_attempts && !stable && consecutive_failures<no_success_failure_stop) {
    attempt_id<-length(ids)+1L; map_id<-length(maps)+1L
    cat(sprintf("CONDITIONAL_BSM_ATTEMPT attempt=%d map_goal=%d/%d\n",attempt_id,map_id,target_maps)); flush.console()
    record<-run_attempt(map_id,attempt_id,bsm_inputs)
    consecutive_failures<-if(record$status=="SUCCESS")0L else consecutive_failures+1L
    cat("CONDITIONAL_BSM_",record$status," attempt=",attempt_id," error=",record$error,"\n",sep=""); flush.console()
    ids<-rebuild_views(); maps<-load_maps()
    if(record$status=="SUCCESS" && length(maps)%%100L==0L) {
      status<-summarize_checkpoint(length(maps)); stable<-isTRUE(status$overall_stable[[1]])
      cat("CONDITIONAL_BSM_CHECKPOINT maps=",length(maps)," status=",status$status,"\n",sep="")
    }
  }
  audit<-validate_attempt_records(maps,through_attempt=length(ids))
  status<-if(stable)"MONTE_CARLO_STABLE_CONDITIONAL_ONLY" else if(consecutive_failures>=no_success_failure_stop)"NO_SUCCESS_FAILURE_STOP_REVIEW_REQUIRED" else if(action=="pilot" && length(maps)==1L)"PILOT_COMPLETE_CONDITIONAL_ONLY" else if(length(maps)>=500L)"INCOMPLETE_AT_500" else if(action=="summarize")"SUMMARY_COMPLETE_CONDITIONAL_ONLY" else if(length(maps)>=target_maps)"TARGET_REACHED_NEEDS_EXTENSION" else "ATTEMPT_CAP_REACHED_REVIEW_REQUIRED"
  result<-list(status=status,scientific_acceptance="NONE",successful=length(maps),attempted=length(ids),
    invocation_pid=Sys.getpid(),action=action,resume=resume,outdir=outdir,invocation_started=invocation_started,
    failed=as.integer(audit$summary$failed),requested_target=target_maps,invocation_attempt_ceiling=max_attempts,
    maximum_successful_maps=500L,global_maximum_attempts=1000L,ledger_pass=isTRUE(audit$summary$all_pass),
    consecutive_failures=consecutive_failures,no_success_failure_stop=no_success_failure_stop,
    conditional_contract_sha256=conditional_contract_hash,runner_sha256=hash_file(script_path),
    fit_sha256=fit_hash,postfit_sha256=postfit_hash,
    source_optimizer_fit_sha256=fit_hash,postfit_ancestral_states_sha256=postfit_hash,
    technical_validation_sha256=technical_validation_hash,finished=timestamp(),downstream_actions_launched=list())
  atomic_write_json_new(result,file.path(outdir,"invocations",paste0("status_",stamp_id(),"_",Sys.getpid(),".json")))
  cat(toJSON(result,auto_unbox=TRUE,pretty=TRUE),"\n")
  if(status%in%c("INCOMPLETE_AT_500","ATTEMPT_CAP_REACHED_REVIEW_REQUIRED","NO_SUCCESS_FAILURE_STOP_REVIEW_REQUIRED")) stop(status)
}
run_body()
