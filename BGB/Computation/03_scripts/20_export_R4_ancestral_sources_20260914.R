#!/usr/bin/env Rscript
# Export existing matrices only: no fitting, likelihood recalculation or graphics.
options(stringsAsFactors=FALSE,digits=17)
Sys.setenv(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
suppressPackageStartupMessages({library(ape);library(digest);library(jsonlite)})
self <- normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)))
job <- dirname(dirname(self)); setwd(job)
out <- file.path(job,"07_tables/ancestral_R4_20260914_v1")
stopifnot(!dir.exists(out))
sha <- function(p) digest(file=p,algo="sha256",serialize=FALSE)
rd <- function(p) read.delim(p,check.names=FALSE)
models <- c("M0","M1","M2")
fitpins <- c(M0="0e425cc6ccc55f810e87e3f6eb531a0b3e2e500ffec581748f82eb303d81725e",M1="bb855fe8bfc2919d5eed94224b177463d517dc78fbe958cc34faf44fcca1bff4",M2="f57721b5d11857211edecd905ba1c9279475c46a884b996f0efad10de4d2a775")
postpins <- c(M0="fe80e8e5791d501536856f1b4c22e6eafcaa66f3fa75b8d9ab232e698bc0e542",M1="566b86c007d123be0fccf1064abf076ad26e6c082ceac8133cbed447ef52e08d",M2="3e23818be5d98b3150e3be09fbfdaa1e1cf1c178c34415cc8349e80db1ad49e2")
stopifnot(sha("01_inputs/frozen/input_manifest.json")=="f7e121a3ce4aff304767eef6df39eba6f9ef473242aad48411a3b2c204d4a306",
sha("01_inputs/frozen/tree_scenarioA417.tre")=="06725d5a0c1ef29007aeb0f7657175c2d2329c75cb4d8fb2170e39b5fdaf6962",
sha("01_inputs/frozen/geog_scenarioA417_analysis_order.LagrangePHYLIP")=="3aa227e4600629eb77178ab68a30bb7cf408189bbcb9f8989ae1f184ce834ca4")
sourcepaths <- c("01_inputs/frozen/input_manifest.json","01_inputs/frozen/tree_scenarioA417.tre","01_inputs/frozen/geog_scenarioA417_analysis_order.LagrangePHYLIP","01_inputs/frozen/tip_tribe_metadata_scenarioA417.tsv","02_config/area_order.tsv")
probs <- fits <- list()
for(m in models) {
  fr <- sprintf("04_runs/%s_R4_final_v1",m); pr <- sprintf("04_runs/%s_R4_postfit_v1",m)
  stopifnot(sha(file.path(fr,"fit_result.rds"))==fitpins[m],sha(file.path(pr,"postfit_recalculated_ancestral_states.rds"))==postpins[m])
  entry <- fromJSON(file.path(pr,"POSTFIT_ENTERED.json"))
  fits[[m]] <- readRDS(file.path(fr,"fit_result.rds"))
  probs[[m]] <- readRDS(file.path(pr,"ancestral_probability_matrices.rds"))
  f <- fits[[m]]; p <- probs[[m]]
  stopifnot(identical(dim(p$top),c(833L,562L)),all(is.finite(p$top)),min(p$top)>=0,max(p$top)<=1,max(abs(rowSums(p$top)-1))<1e-12,f$inputs$max_range_size==4L)
  sourcepaths <- c(sourcepaths,file.path(fr,c("fit_result.rds","fit_summary.tsv")),file.path(pr,c("POSTFIT_ENTERED.json","postfit_recalculated_ancestral_states.rds","ancestral_probability_matrices.rds","postfit_QA.tsv")))
}
dict <- probs$M1$state_dictionary
stopifnot(identical(dict,probs$M0$state_dictionary),identical(dict,probs$M2$state_dictionary),nrow(dict)==562L,identical(as.integer(dict$state_index),1:562))
tr <- read.tree("01_inputs/frozen/tree_scenarioA417.tre"); nt <- Ntip(tr); nn <- nt+tr$Nnode
stopifnot(nt==417L,tr$Nnode==416L,nrow(tr$edge)==832L,!anyDuplicated(tr$tip.label),all(tr$edge.length>=0))
area <- rd("02_config/area_order.tsv");area <- area[order(area$analysis_bit_index),]
stopifnot(identical(area$abbrev,c("AF","AUS","CAM","ENA","EPA","IND","MDG","ORI","SAM","WNA","WPA")))
lines <- readLines("01_inputs/frozen/geog_scenarioA417_analysis_order.LagrangePHYLIP")
geo <- read.table(text=lines[-1],col.names=c("tip_label","geography_bits"),colClasses="character")
stopifnot(nrow(geo)==417L,!anyDuplicated(geo$tip_label),setequal(geo$tip_label,tr$tip.label))
geo <- geo[match(tr$tip.label,geo$tip_label),]
observed <- vapply(strsplit(geo$geography_bits,""),function(b) paste(area$abbrev[b=="1"],collapse="+"),character(1))
tipstate <- match(observed,dict$semantic_state);stopifnot(!anyNA(tipstate))
for(m in models) stopifnot(all(probs[[m]]$top[cbind(1:417,tipstate)]==1))
meta <- rd("01_inputs/frozen/tip_tribe_metadata_scenarioA417.tsv")
stopifnot(nrow(meta)==417L,!anyDuplicated(meta$target_species),setequal(meta$target_species,tr$tip.label))
meta <- meta[match(tr$tip.label,meta$target_species),];stopifnot(!anyNA(meta$tribe))
parent <- rep(NA_integer_,nn);parent[tr$edge[,2]] <- tr$edge[,1]
children <- split(tr$edge[,2],tr$edge[,1]);descendants <- vector("list",nn)
desc <- function(n) {if(!is.null(descendants[[n]])) return(descendants[[n]]);ans <- if(n<=nt)n else unlist(lapply(children[[as.character(n)]],desc),use.names=FALSE);descendants[[n]] <<- ans;ans}
stopifnot(identical(which(is.na(parent)),418L));invisible(desc(418L))
depth <- node.depth.edgelength(tr)
age <- vapply(1:nn,function(n)max(depth[descendants[[n]]])-depth[n],numeric(1))
tiphash <- vapply(descendants,function(t)digest(paste0(paste(sort(enc2utf8(tr$tip.label[t]),method="radix"),collapse="\n"),"\n"),algo="sha256",serialize=FALSE),character(1))
ref <- file.path(dirname(job),"R2_4_ScenarioA417_timestrat_DEC_20260904_232443/02_config")
apath <- file.path(ref,"fig3_node_anchors.tsv");dpath <- file.path(ref,"fig3_display_anchors.tsv")
stopifnot(sha(apath)=="a77b4629036f9a34615b8b06759e3d1d1507254c7b0833d61d76510db4c659ab",sha(dpath)=="d7484b5cb52b6bc2c583552008e78c06de24742a27190e34a4e55e27a341509e")
a <- rd(apath);display <- rd(dpath);a <- a[match(display$anchor_id,a$anchor_id),]
an <- match(a$descendant_tipset_sha256,tiphash)
stopifnot(!anyNA(an),all(an==a$derived_ape_node_number),all(lengths(descendants)[an]==a$descendant_tip_count),max(abs(age[an]-a$node_age_Ma))<1e-8)
sourcepaths <- c(sourcepaths,apath,dpath)
inc <- t(vapply(strsplit(dict$semantic_state,"+",fixed=TRUE),function(s)as.integer(area$abbrev %in% s),integer(11)))
stopifnot(all(rowSums(inc)==dict$range_size))
tv <- list(M0_M1=rowSums(abs(probs$M0$top-probs$M1$top))/2,M1_M2=rowSums(abs(probs$M1$top-probs$M2$top))/2,M0_M2=rowSums(abs(probs$M0$top-probs$M2$top))/2)
summarize_node <- function(m,n) {p<-probs[[m]]$top[n,];ii<-order(-p,seq_along(p))[1:3];data.frame(model=m,ape_node=n,node_age_Ma=age[n],top1_state=dict$semantic_state[ii[1]],top1_prob=p[ii[1]],top2_state=dict$semantic_state[ii[2]],top2_prob=p[ii[2]],top3_state=dict$semantic_state[ii[3]],top3_prob=p[ii[3]],other_prob=1-sum(p[ii]),entropy_nats=-sum(p[p>0]*log(p[p>0])),p_range4=sum(p[dict$range_size==4]),TV_M0_M1=tv$M0_M1[n],TV_M1_M2=tv$M1_M2[n],TV_M0_M2=tv$M0_M2[n])}
summary <- do.call(rbind,lapply(models,function(m)do.call(rbind,lapply(seq_along(an),function(i)cbind(data.frame(anchor_number=display$anchor_number[i],anchor_name=display$display_label[i],anchor_kind=a$anchor_kind[i]),summarize_node(m,an[i]))))))
node_summary <- do.call(rbind,lapply(models,function(m)do.call(rbind,lapply(1:nn,function(n)summarize_node(m,n)))))
dir.create(out,recursive=TRUE)
write_tsv <- function(x,name) {for(j in seq_along(x))if(is.double(x[[j]]))x[[j]]<-ifelse(is.na(x[[j]]),NA_character_,sprintf("%.17g",x[[j]]));write.table(x,file.path(out,name),sep="\t",quote=FALSE,row.names=FALSE,na="NA")}
write_tsv(data.frame(ape_node=1:nn,is_tip=1:nn<=nt,tip_label=c(tr$tip.label,rep("",tr$Nnode)),parent=parent,node_age_Ma=age,root_distance=depth,desc_tipcount=lengths(descendants),tipset_sha256=tiphash),"nodes.tsv")
write_tsv(data.frame(parent=tr$edge[,1],child=tr$edge[,2],branch_length_Ma=tr$edge.length),"edges.tsv")
write_tsv(data.frame(state_index_1based=dict$state_index,semantic_state=dict$semantic_state,range_size=dict$range_size),"state_dictionary.tsv")
write_tsv(data.frame(tip_label=tr$tip.label,tribe=meta$tribe,geography_bits=geo$geography_bits,observed_state=observed),"tip_metadata.tsv")
write_tsv(summary,"anchors_summary.tsv");write_tsv(node_summary,"node_summaries.tsv")
write_tsv(display,"display_anchors.tsv");write_tsv(a,"anchor_definitions.tsv");write_tsv(area,"area_order.tsv")
for(m in models) {wide<-as.data.frame(probs[[m]]$top);names(wide)<-sprintf("p%04d",1:562);write_tsv(cbind(data.frame(ape_node=1:833),wide),paste0(m,"_node_posteriors.tsv"));mar<-as.data.frame(probs[[m]]$top %*% inc);names(mar)<-area$abbrev;write_tsv(cbind(data.frame(ape_node=1:833),mar),paste0(m,"_area_inclusion.tsv"))}
model_table <- do.call(rbind,lapply(models,function(m)rd(sprintf("04_runs/%s_R4_final_v1/fit_summary.tsv",m))))
model_table$delta_AIC <- model_table$AIC-min(model_table$AIC)
model_table$relative_AIC_weight_R4_only <- exp(-model_table$delta_AIC/2)/sum(exp(-model_table$delta_AIC/2))
model_table$prespecified_role <- c("static_baseline","conservative_main","permissive_sensitivity")
write_tsv(model_table,"model_comparison.tsv")
outputpaths <- list.files(out,full.names=TRUE)
manifest <- list(status="CONDITIONAL_NEW_R4_NATIVE_POSTFIT",exported_at=format(Sys.time(),"%Y-%m-%dT%H:%M:%S%z"),new_fits_or_postfits=0,scientific_acceptance="AUTHOR_REVIEW_REQUIRED_KKT1_FALSE",models=models,scope="Pelopidas mathias legacy R4; other 416 audited new rows; no R6 pooling",exporter=list(path=self,sha256=sha(self)),sources=lapply(sourcepaths,function(p)list(path=normalizePath(p),sha256=sha(p))),outputs=lapply(outputpaths,function(p)list(file=basename(p),sha256=sha(p))),validation=list(tips=417L,nodes=833L,states=562L,anchors=15L,tree_geog_metadata_exact_join=TRUE,tip_probabilities_match_observed_ranges=TRUE,anchor_descendant_hash_match=TRUE,state_order_equal=TRUE,max_row_sum_error=lapply(probs,function(p)max(abs(rowSums(p$top)-1)))),definitions=list(probabilities="Existing independent postfit node-top marginal matrices, all 562 states retained",area_inclusion="Sum of full-state probabilities over states containing each area; not a distribution over areas",TV="0.5 sum absolute differences between complete state distributions; not uncertainty interval",AIC="Conditional descriptive comparison within identical new R4 inputs; no model averaging or global optimum assertion",node_age="Maximum descendant depth minus node depth; plotting uses original root distances"))
write_json(manifest,file.path(out,"source_manifest.json"),auto_unbox=TRUE,pretty=TRUE,digits=17,na="null")
cat("R4_SOURCE_EXPORT_COMPLETE ",out,"\n",sep="")
