#!/usr/bin/env Rscript
# Zero-likelihood, zero-optimizer verification of new, independently built inputs.
suppressPackageStartupMessages({library(ape); library(BioGeoBEARS); library(digest); library(jsonlite)})
script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)), mustWork=TRUE)
job <- normalizePath(file.path(dirname(script), ".."), mustWork=TRUE)
sha <- function(p) digest(file=p, algo="sha256", serialize=FALSE)
mf <- file.path(job,"01_inputs/frozen/input_manifest.json")
stopifnot(sha(mf)=="f7e121a3ce4aff304767eef6df39eba6f9ef473242aad48411a3b2c204d4a306")
manifest <- read_json(mf, simplifyVector=FALSE)
stopifnot(manifest$schema_version=="geography_pelopidas_legacy_exception_v1", manifest$job_root==job)
for (item in manifest$files) stopifnot(sha(file.path(job,item$relative_path))==item$sha256)
treefn <- file.path(job,"01_inputs/frozen/tree_scenarioA417.tre")
geogfn <- file.path(job,"01_inputs/frozen/geog_scenarioA417_analysis_order.LagrangePHYLIP")
tr <- read.tree(treefn)
geog <- getranges_from_LagrangePHYLIP(geogfn)
g <- as.matrix(geog@df); storage.mode(g)<-"numeric"
stopifnot(length(tr$tip.label)==417L, tr$Nnode==416L, is.rooted(tr), is.binary(tr),
          !anyDuplicated(tr$tip.label), !anyDuplicated(rownames(g)),
          setequal(tr$tip.label,rownames(g)), identical(dim(g),c(417L,11L)),
          identical(colnames(g),LETTERS[1:11]), all(g %in% c(0,1)),
          all(rowSums(g)>0), max(rowSums(g))==4, all(tr$edge.length>0))
depth <- node.depth.edgelength(tr)[seq_len(417)]
stopifnot(diff(range(depth))<=5.1e-6, abs(max(depth)-42.0830895)<1e-4)
states <- cladoRcpp::rcpp_areas_list_to_states_list(areas=LETTERS[1:11],maxareas=4L,include_null_range=TRUE)
if(length(states[[1]])==1L && !is.na(states[[1]]) && states[[1]]=="_") states[[1]]<-NA
tips <- tipranges_to_tip_condlikes_of_data_on_each_state(geog,tr,states_list=states,maxareas=4L,include_null_range=TRUE)
stopifnot(length(states)==562L, identical(dim(tips),c(417L,562L)), all(rowSums(tips)==1), all(tips %in% c(0,1)))
stopifnot(paste(g["Pelopidas_mathias",],collapse="")=="10001101000")
runner <- file.path(job,"03_scripts/04_fit_dec_model_R4_v1.R")
results <- list()
write_new <- function(x,p) {
  if(file.exists(p) || dir.exists(p) || (!is.na(Sys.readlink(p)) && nzchar(Sys.readlink(p)))) stop("Refuse existing QA output: ",p)
  write_json(x,p,auto_unbox=TRUE,pretty=TRUE,digits=NA)
}
for(model in c("M1","M2","M0")) {
  d <- file.path(job,"04_runs",paste0(model,"_R4_prepare_v1"))
  status <- read.delim(file.path(d,"STATUS.tsv"),stringsAsFactors=FALSE)
  stopifnot(nrow(status)==1L, status$status=="VALIDATED_NOT_FITTED", status$cores==2L)
  prefit <- file.path(d,"prefit_run.rds"); run <- readRDS(prefit)
  stopifnot(run$provenance$model==model, run$provenance$mode=="prepare", run$max_range_size==4,
            run$include_null_range, !run$speedup, !run$force_sparse, run$num_cores_to_use==2,
            run$provenance$input_manifest_sha256==sha(mf),
            normalizePath(run$trfn)==treefn, normalizePath(run$geogfn)==geogfn,
            is.null(run$states_list), isTRUE(check_BioGeoBEARS_run(run)))
  p <- run$BioGeoBEARS_model_object@params_table
  stopifnot(identical(rownames(p)[p[,"type"]=="free"],c("d","e")),
            as.numeric(p["j","est"])==0, as.numeric(p["w","est"])==1)
  if(model!="M0") {
    stopifnot(identical(as.numeric(run$timeperiods),c(3.3,13.9,23.03,33.9,45)),length(run$list_of_dispersal_multipliers_mats)==5)
    for(raw_matrix in run$list_of_dispersal_multipliers_mats) {
      m <- as.matrix(raw_matrix)
      stopifnot(identical(dim(m),c(11L,11L)),
        max(abs(m-t(m)))<1e-15, all(diag(m)==1),all(m>0))
    }
  } else stopifnot(length(run$timeperiods)==0, length(run$list_of_dispersal_multipliers_mats)==0)
  qa <- list(status="PASS_NATIVE_INPUT_QA", model=model,max_range=4L,n_tips=417L,n_internal=416L,
             n_states=562L, input_manifest_sha256=sha(mf), fit_runner_sha256=sha(runner),
             exact_tip_state_matches=417L,likelihood_calls=0L,optimizer_calls=0L,
             R6_inputs_or_processes_modified=FALSE,scientific_acceptance="NOT_ESTABLISHED_INPUT_QA_ONLY")
  write_new(qa,file.path(d,"native_input_QA.json"))
  write_new(list(schema_version="r4_exception_native_preparation_v1",model=model,
       input_manifest_sha256=sha(mf),fit_runner_sha256=sha(runner),
       prefit_run_sha256=sha(prefit),native_input_QA_sha256=sha(file.path(d,"native_input_QA.json")),
       QA_script_sha256=sha(script),R_version=as.character(getRversion()),
       BioGeoBEARS_version=as.character(packageVersion("BioGeoBEARS")),
       BioGeoBEARS_remote_sha=unname(packageDescription("BioGeoBEARS")$RemoteSha),
       likelihood_calls=0L,optimizer_calls=0L),file.path(d,"preflight_manifest.json"))
  results[[model]]<-qa
}
write_new(list(status="PASS_NATIVE_INPUT_QA",models=results,tree_sha256=sha(treefn),
  geography_sha256=sha(geogfn),root_age_ma=max(depth),root_to_tip_range_ma=diff(range(depth)),
  likelihood_calls=0L,optimizer_calls=0L),file.path(job,"08_qa/native_input_QA.json"))
cat("PASS_NATIVE_INPUT_QA models=M0,M1,M2 tips=417 states=562 likelihood_calls=0 optimizer_calls=0\n")
