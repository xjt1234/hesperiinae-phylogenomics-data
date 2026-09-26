#!/usr/bin/env Rscript
# Bounded tests only: no whole-tree BSM, optimization, or cache preparation.
args<-commandArgs(TRUE)
if(length(args)!=2L || args[1L]!="--outdir") stop("Use --outdir NEW_ABSOLUTE_QA_DIRECTORY")
script_arg<-grep("^--file=",commandArgs(FALSE),value=TRUE)
script_path<-normalizePath(sub("^--file=","",script_arg),mustWork=TRUE)
job<-normalizePath(file.path(dirname(script_path),"../.."),mustWork=TRUE)
qa<-normalizePath(file.path(job,"08_qa/bsm_bridge_preflight_20260914"),mustWork=TRUE)
outdir<-file.path(normalizePath(dirname(args[2L]),mustWork=TRUE),basename(args[2L]))
if(dirname(outdir)!=qa || dir.exists(outdir)) stop("Output must be new direct child of diagnosis QA directory")
source(file.path(dirname(script_path),"74_ctmc_uniformization_bridge_20260908.R"))
suppressPackageStartupMessages(library(Matrix))
checks<-list()
check<-function(name,pass,detail="") {
  checks[[length(checks)+1L]]<<-data.frame(test=name,pass=isTRUE(pass),detail=as.character(detail),stringsAsFactors=FALSE)
  if(!isTRUE(pass)) stop("SELFTEST_FAIL: ",name," ",detail)
}
fails<-function(expr) inherits(tryCatch({force(expr);NULL},error=identity),"error")
ctx<-new_bsm_bridge_context()
set.seed(202609141L)
a<-.7; b<-.4; t<-.8; Q<-matrix(c(-a,a,b,-b),2,byrow=TRUE)
p11<-b/(a+b)+a/(a+b)*exp(-(a+b)*t)
p12<-a/(a+b)*(1-exp(-(a+b)*t))
z11<-ctx$kernel(Q,t,1L,1L); z12<-ctx$kernel(Q,t,1L,2L)
check("two_state_P11",abs(z11$diagnostics$endpoint_probability-p11)<1e-12)
check("two_state_P12",abs(z12$diagnostics$endpoint_probability-p12)<1e-12)
check("two_state_tail",max(z11$diagnostics$relative_tail_bound,z12$diagnostics$relative_tail_bound)<=1e-12)
draws11<-replicate(1000L,ctx$kernel(Q,t,1L,1L),simplify=FALSE)
draws12<-replicate(300L,ctx$kernel(Q,t,1L,2L),simplify=FALSE)
counts11<-vapply(draws11,function(x)length(x$times),integer(1)); counts12<-vapply(draws12,function(x)length(x$times),integer(1))
check("two_state_jump_parity",all(counts11%%2L==0L) && all(counts12%%2L==1L))
p_nojump<-exp(-a*t)/p11; observed<-mean(counts11==0L)
check("two_state_zero_jump_frequency",abs(observed-p_nojump)<6*sqrt(p_nojump*(1-p_nojump)/1000)+.002,
      sprintf("observed=%.8f expected=%.8f",observed,p_nojump))
reject_once<-function(Q,t,a,b) {
  repeat {
    time<-0; cur<-a; states<-a; times<-numeric()
    repeat {
      rate<--Q[cur,cur]
      if(rate==0) break
      time<-time+rexp(1L,rate)
      if(time>=t) break
      pr<-Q[cur,];pr[cur]<-0
      cur<-sample.int(nrow(Q),1L,prob=pr); states<-c(states,cur);times<-c(times,time)
    }
    if(cur==b)return(list(states=states,times=times))
  }
}
set.seed(202609142L)
rejection<-replicate(1000L,reject_once(Q,t,1L,1L),simplify=FALSE)
counts_rej<-vapply(rejection,function(x)length(x$times),integer(1))
se_count<-sqrt(var(counts11)/length(counts11)+var(counts_rej)/length(counts_rej))
check("rejection_mean_jump_agreement",abs(mean(counts11)-mean(counts_rej))<6*se_count+.002)
dwell1<-function(x)sum(diff(c(0,x$times,t))[x$states==1L])
dwell_bridge<-vapply(draws11,dwell1,numeric(1));dwell_rej<-vapply(rejection,dwell1,numeric(1))
se_dwell<-sqrt(var(dwell_bridge)/length(dwell_bridge)+var(dwell_rej)/length(dwell_rej))
check("rejection_dwell_agreement",abs(mean(dwell_bridge)-mean(dwell_rej))<6*se_dwell+.002)
eps<-1e-12; Qtiny<-matrix(c(-eps,eps,0,0,-eps,eps,0,0,0),3,byrow=TRUE)
ztiny<-ctx$kernel(Qtiny,1,1L,3L); tiny_expected<-exp(-eps)*eps^2/2
check("tiny_probability_relative_accuracy",abs(ztiny$diagnostics$endpoint_probability/tiny_expected-1)<1e-10)
check("tiny_probability_tail",ztiny$diagnostics$relative_tail_bound<=1e-12 && ztiny$diagnostics$last_series_term>=2L)
check("tiny_probability_endpoints",identical(ztiny$states,c(1L,2L,3L)))
Qsub<-matrix(c(-1,.7,.4,-.9),2,byrow=TRUE)
Qcem<-matrix(c(-1,.7,.3,.4,-.9,.5,0,0,0),3,byrow=TRUE)
zsub<-ctx$kernel(Qsub,.8,1L,2L); zcem<-ctx$kernel(Qcem,.8,1L,2L)
p_expm<-as.matrix(expm::expm(Qsub*.8))[1,2]
check("subgenerator_expm",abs(zsub$diagnostics$endpoint_probability-p_expm)<1e-12)
check("subgenerator_cemetery_equivalence",abs(zsub$diagnostics$endpoint_probability-zcem$diagnostics$endpoint_probability)<1e-12)
check("native_converter_rejects_material_killing",fails(bgb_native_effective_generator(Qsub)))
check("unreachable_rejected",fails(ctx$kernel(matrix(0,2,2),1,1L,2L)))
zero<-ctx$kernel(matrix(0,2,2),1,1L,1L)
check("zero_generator_no_events",!length(zero$times) && identical(zero$states,1L))
check("invalid_negative_rate_rejected",fails(ctx$kernel(matrix(c(-1,1,-.1,.1),2,byrow=TRUE),1,1L,2L)))
set.seed(777L);repro1<-ctx$kernel(Q,t,1L,2L)
set.seed(777L);repro2<-ctx$kernel(Q,t,1L,2L)
check("seed_reproducibility",identical(repro1$states,repro2$states) && identical(repro1$times,repro2$times))
check("event_time_bounds",all(repro1$times>0 & repro1$times<t) && !is.unsorted(repro1$times))
Qrounded<-Q;Qrounded[1L,2L]<-Qrounded[1L,2L]+1e-8
eff<-bgb_native_effective_generator(Qrounded)
check("native_effective_generator_row_sum",max(abs(rowSums(eff$Q)))<1e-14)
check("native_effective_generator_rate_preserved",identical(diag(eff$Q),diag(Qrounded)))
trtable<-data.frame(edge.length=t,time_bp=0,sampled_states_AT_brbots=1L,sampled_states_AT_nodes=2L,node=2L)
ctx$reset_audit();set.seed(55L)
events<-ctx$branch(1L,trtable,Q,list(0L,c(0L,1L)),c("A","AB"),c("A","B"))
cols<-c("nodenum_at_top_of_branch","trynum","brlen","current_rangenum_1based","new_rangenum_1based","current_rangetxt","new_rangetxt","abs_event_time","event_time","event_type","event_txt","new_area_num_1based","lost_area_num_1based","dispersal_to","extirpation_from")
check("native_event_schema",is.data.frame(events) && ncol(events)==15L && identical(names(events),cols))
check("native_event_endpoint",events$current_rangetxt[1L]=="A" && tail(events$new_rangetxt,1L)=="AB")
check("native_event_age_transform",max(abs(as.numeric(events$abs_event_time)+as.numeric(events$event_time)-t))<1e-14)
check("native_audit_separate",nrow(ctx$get_audit())==1L && ctx$get_audit()$real_jumps==nrow(events))


# Synthetic guards only: old fitted-Q or failed-branch fixtures are NOT imported.
bad_rounding<-Q;bad_rounding[1L,2L]<-bad_rounding[1L,2L]+1e-5
check("native_rounding_budget_rejects_large_discrepancy",fails(bgb_native_effective_generator(bad_rounding)))
result<-do.call(rbind,checks)
dir.create(outdir,recursive=FALSE)
write.table(result,file.path(outdir,"checks.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
saveRDS(list(checks=result,analytic=list(P11=p11,P12=p12,no_jump=p_nojump),tiny=ztiny,
  subgenerator=zsub,native_events=events,native_audit=ctx$get_audit(),cache=ctx$cache_stats(),
  session=capture.output(sessionInfo())),file.path(outdir,"selftest.rds"))
jsonlite::write_json(list(status="PASS",all_checks_pass=all(result$pass),n_checks=nrow(result),
  scope="GENERIC_SYNTHETIC_ALGORITHM_SELFTEST_ONLY",whole_tree_maps=0L,current_fit_tested=FALSE,
  historical_scientific_fixtures_imported=FALSE,
  helper_sha256=digest::digest(file=file.path(dirname(script_path),"74_ctmc_uniformization_bridge_20260908.R"),algo="sha256",serialize=FALSE),
  selftest_script_sha256=digest::digest(file=script_path,algo="sha256",serialize=FALSE),
  next_gate="New R4 native-input preflight followed by freshly prepared one-map pilot",
  created=format(Sys.time(),"%Y-%m-%dT%H:%M:%OS6Z",tz="UTC")),
  file.path(outdir,"selftest_summary.json"),auto_unbox=TRUE,pretty=TRUE)
cat("CTMC_BRIDGE_GENERIC_SELFTEST_PASS tests=",nrow(result)," whole_tree_maps=0 current_fit_tested=FALSE\n",sep="")
