#!/usr/bin/env Rscript
# New exact-input downstream driver. It never changes or discovers biological inputs.
options(stringsAsFactors=FALSE, warn=1, digits=7)
required_packages <- c(paco='0.4.2',ape='5.8.1',vegan='2.7.2',jsonlite='2.0.0',digest='0.6.39')
need <- function(ok,message) if (!isTRUE(ok)) stop(message,call.=FALSE)
sha <- function(path) digest::digest(file=path,algo='sha256',serialize=FALSE)
jsave <- function(x,path) {need(!file.exists(path),paste('Refusing overwrite:',path));jsonlite::write_json(x,path,auto_unbox=TRUE,pretty=TRUE,digits=17,null='null')}
tsave <- function(x,path) {need(!file.exists(path),paste('Refusing overwrite:',path));write.table(x,path,sep='\t',quote=FALSE,row.names=FALSE,na='NA',fileEncoding='UTF-8')}
read_tsv <- function(path) read.delim(path,check.names=FALSE,stringsAsFactors=FALSE,na.strings=character(),quote='',comment.char='')
record_path <- function(rec,what) {need(is.list(rec)&&is.character(rec$path)&&length(rec$path)==1&&startsWith(rec$path,'/')&&file.exists(rec$path),paste('Missing absolute input:',what));need(identical(sha(rec$path),rec$sha256),paste('SHA mismatch:',what));normalizePath(rec$path)}
checkpoint_outputs <- function(cp) {
  records<-cp$outputs
  if(!is.null(cp$output_sha256)){need(is.list(cp$output_sha256)&&length(names(cp$output_sha256))==length(cp$output_sha256),'Invalid checkpoint output_sha256 mapping');records<-c(records,lapply(names(cp$output_sha256),function(p)list(path=p,sha256=cp$output_sha256[[p]])))}
  need(is.list(records)&&length(records)>0,'Completed checkpoint has no path/hash outputs')
  records
}
scalar_text <- function(x) is.character(x)&&length(x)==1&&!is.na(x)&&nzchar(trimws(x))
validate_ids <- function(ids,what) {need(is.character(ids)&&length(ids)>0&&!anyNA(ids)&&all(nzchar(ids))&&!any(grepl('[[:space:]]',ids))&&!anyDuplicated(ids),paste('Missing/duplicate/whitespace exact IDs:',what))}
mc_result <- function(null,observed,tail) {
  need(length(null)>0&&all(is.finite(null))&&is.finite(observed),'Invalid null/observed statistic')
  tolerance <- 64*.Machine$double.eps*max(1,abs(observed))
  extreme <- if(tail=='lower')null<=observed+tolerance else null>=observed-tolerance
  list(p=(1+sum(extreme))/(length(null)+1),extreme=extreme,count=sum(extreme),tolerance=tolerance,tail=tail,rule='(1 + extreme_count) / (nperm + 1); inclusive ties within declared numerical tolerance')
}
pc_vectors <- function(pc) {v<-if(as.integer(pc$correction[2])==1)pc$vectors else pc$vectors.cor;need(is.matrix(v)&&ncol(v)>0&&all(is.finite(v)),'Invalid or zero-rank PCo coordinates');v}
fit_paco <- function(Hcoords,Wcoords,HP,keep=FALSE) {
  pos<-which(HP>0,arr.ind=TRUE)
  X<-Hcoords[pos[,1],,drop=FALSE];Y<-Wcoords[pos[,2],,drop=FALSE]
  proc<-suppressWarnings(vegan::procrustes(X,Y,symmetric=FALSE))
  need(is.finite(proc$ss),'Nonfinite PACo fit (possibly degenerate null association); no silent replacement/redraw')
  if(keep)list(proc=proc,pos=pos) else proc$ss
}
fit_parafit <- function(Hcoords,Wcoords,HP) sum((t(Hcoords)%*%HP%*%Wcoords)^2)
read_scenario <- function(s,synthetic) {
  need(length(s$scenario_index)==1&&!is.na(as.integer(s$scenario_index))&&as.integer(s$scenario_index)>=0&&as.double(s$scenario_index)==as.integer(s$scenario_index),'Explicit nonnegative scenario_index is required')
  need(scalar_text(s$id)&&grepl('^[A-Za-z0-9][A-Za-z0-9_.-]*$',s$id),'Unsafe/missing scenario id')
  for(k in c('host_tree_basis','wol_tree_basis','association_basis','host_distance_unit','wol_distance_unit'))need(scalar_text(s[[k]]),paste('Missing scenario semantics:',k))
  paths<-lapply(c('host_tree','wol_tree','links'),function(k)record_path(s$inputs[[k]],paste(s$id,k)));names(paths)<-c('host_tree','wol_tree','links')
  links<-read_tsv(paths$links);need(all(c('host_tip','wol_tip')%in%names(links)),'Links require exact host_tip and wol_tip')
  validate_ids(links$host_tip,'host_tip');validate_ids(links$wol_tip,'wol_tip')
  need(nrow(links)==as.integer(s$n_links)&&nrow(links)>=3,'Scenario n_links mismatch or fewer than three pairs')
  need(synthetic||nrow(links)==25,'This frozen driver permits real 25-pair scenarios only')
  if(all(c('group','link_group')%in%names(links)))need(identical(links$group,links$link_group),'Conflicting group aliases')
  for(role in c('host_tip','wol_tip'))need(!anyDuplicated(links[[role]]),'Not a bijection')
  host<-ape::read.tree(paths$host_tree);wol<-ape::read.tree(paths$wol_tree)
  for(item in list(list(tree=host,ids=links$host_tip,role='host'),list(tree=wol,ids=links$wol_tip,role='wol'))){
    tr<-item$tree;need(inherits(tr,'phylo'),'Expected one Newick tree');validate_ids(tr$tip.label,paste(item$role,'tree'))
    need(setequal(tr$tip.label,item$ids),paste('Exact tree/link tips mismatch:',item$role))
    need(!is.null(tr$edge.length)&&all(is.finite(tr$edge.length))&&all(tr$edge.length>=0),'Missing/nonfinite/negative tree branch length')
  }
  hids<-sort(host$tip.label,method='radix');wids<-sort(wol$tip.label,method='radix')
  HP<-matrix(0L,length(hids),length(wids),dimnames=list(hids,wids));HP[cbind(match(links$host_tip,hids),match(links$wol_tip,wids))]<-1L
  need(all(rowSums(HP)==1L)&&all(colSums(HP)==1L),'Observed links must be a strict one-to-one bijection')
  H<-ape::cophenetic.phylo(host)[hids,hids,drop=FALSE];W<-ape::cophenetic.phylo(wol)[wids,wids,drop=FALSE]
  for(d in list(H,W))need(all(is.finite(d))&&all(d>=0)&&max(abs(d-t(d)))<1e-10&&all(diag(d)==0)&&any(d>0),'Invalid phylogenetic distance matrix')
  list(id=s$id,spec=s,paths=paths,links=links,H=H,W=W,HP=HP)
}
run_scenario <- function(d,out,nperm,seed) {
  dir.create(out,recursive=FALSE);cat(format(Sys.time(),tz='UTC'),d$id,'observed fits\n')
  H<-d$H;W<-d$W;HP<-d$HP;n<-nrow(HP)
  coord<-get('coordpcoa',asNamespace('paco'))
  hc<-coord(H,correction='cailliez');wc<-coord(W,correction='cailliez');HX<-pc_vectors(hc);WY<-pc_vectors(wc)
  ah<-ape::pcoa(H,correction='cailliez');aw<-ape::pcoa(W,correction='cailliez');AH<-pc_vectors(ah);AW<-pc_vectors(aw)
  observed<-fit_paco(HX,WY,HP,TRUE);ss<-observed$proc$ss;pf<-fit_parafit(AH,AW,HP)
  referenceD<-paco::add_pcoord(paco::prepare_paco_data(H,W,HP),correction='cailliez')
  referenceProc<-suppressWarnings(vegan::procrustes(referenceD$H_PCo,referenceD$P_PCo,symmetric=FALSE))
  need(isTRUE(all.equal(ss,referenceProc$ss,tolerance=1e-10)),'Cached PACo differs from public add_pcoord path')
  referencePF<-ape::parafit(H,W,HP,nperm=0,test.links=FALSE,seed=seed,correction='cailliez',silent=TRUE)
  need(isTRUE(all.equal(pf,referencePF$ParaFitGlobal,tolerance=1e-10)),'ParaFit statistic differs from ape::parafit')
  raw<-unname(paco::residuals_paco(observed$proc,type='interaction'));pos<-observed$pos
  need(length(raw)==n&&all(is.finite(raw))&&all(raw>=0)&&abs(sum(raw^2)-ss)<=1e-7*max(1,abs(ss)),'Raw interaction residuals do not reconstruct PACo SS')
  rh<-rownames(H)[pos[,1]];rw<-rownames(W)[pos[,2]]
  need(identical(rownames(observed$proc$X),rh)&&identical(rownames(observed$proc$Yrot),rw),'Procrustes row labels changed; cannot map residuals safely')
  indices<-match(rh,d$links$host_tip);need(identical(d$links$wol_tip[indices],rw),'Residual host/wol exact pairing mismatch')
  residuals<-d$links[indices,,drop=FALSE];residuals$residual<-raw;residuals$residual_sq<-raw^2
  tsave(residuals,file.path(out,'paco_link_residuals.tsv'))
  tsave(data.frame(host_tip=rownames(HP),HP,check.names=FALSE),file.path(out,'HP_associations.tsv'))
  tsave(data.frame(host_tip=rownames(H),H,check.names=FALSE),file.path(out,'host_distances.tsv'))
  tsave(data.frame(wol_tip=rownames(W),W,check.names=FALSE),file.path(out,'wol_distances.tsv'))
  saveRDS(list(paco_host=hc,paco_wol=wc,parafit_host=ah,parafit_wol=aw,observed_proc=observed$proc),file.path(out,'observed_coordinates_and_proc.rds'),version=3)
  RNGkind('Mersenne-Twister','Inversion','Rejection');set.seed(seed)
  permutations<-t(replicate(nperm,sample.int(n)));need(nrow(permutations)==nperm&&ncol(permutations)==n,'Bad permutation dimensions')
  primarySS<-primaryPF<-numeric(nperm)
  for(i in seq_len(nperm)){
    permHP<-HP[permutations[i,],,drop=FALSE];dimnames(permHP)<-dimnames(HP)
    need(all(rowSums(permHP)==1L)&&all(colSums(permHP)==1L),'Primary permutation violated bijection')
    primarySS[i]<-fit_paco(HX,WY,permHP);primaryPF[i]<-fit_parafit(AH,AW,permHP)
    if(i%%1000==0)cat(format(Sys.time(),tz='UTC'),d$id,'bijection',i,'/',nperm,'\n')
  }
  pmc<-mc_result(primarySS,ss,'lower');fmc<-mc_result(primaryPF,pf,'upper')
  tsave(data.frame(permutation=seq_len(nperm),paco_ss=primarySS,parafit_global=primaryPF,paco_extreme=pmc$extreme,parafit_extreme=fmc$extreme,bijection_verified=TRUE),file.path(out,'bijection_null_distribution.tsv'))
  colnames(permutations)<-paste0('host_slot_',seq_len(n));tsave(data.frame(permutation=seq_len(nperm),permutations,check.names=FALSE),file.path(out,'bijection_host_label_permutations.tsv'))
  tsave(data.frame(slot=seq_len(n),host_tip=rownames(HP),observed_wol_tip=colnames(HP)[max.col(HP,ties.method='first')]),file.path(out,'permutation_slot_key.tsv'))
  RNGkind('Mersenne-Twister','Inversion','Rejection');set.seed(seed)
  rands<-stats::simulate(vegan::nullmodel(HP,'r0'),nsim=nperm)
  rSS<-numeric(nperm);rassign<-matrix(NA_integer_,nperm,n);rWcount<-integer(nperm)
  for(i in seq_len(nperm)){
    permHP<-rands[,,i];dimnames(permHP)<-dimnames(HP)
    need(all(rowSums(permHP)==1L)&&all(permHP%in%c(0L,1L)),'r0 draw does not preserve observed row sums')
    rSS[i]<-fit_paco(HX,WY,permHP);rassign[i,]<-max.col(permHP,ties.method='first');rWcount[i]<-sum(colSums(permHP)>0)
    if(i%%1000==0)cat(format(Sys.time(),tz='UTC'),d$id,'r0 sensitivity',i,'/',nperm,'\n')
  }
  rmc<-mc_result(rSS,ss,'lower')
  tsave(data.frame(permutation=seq_len(nperm),paco_ss=rSS,paco_extreme=rmc$extreme,occupied_wol_taxa=rWcount,all_host_degree_one=TRUE,both_margins_preserved=rWcount==n),file.path(out,'paco_r0_null_distribution.tsv'))
  colnames(rassign)<-paste0('host_slot_',seq_len(n));tsave(data.frame(permutation=seq_len(nperm),rassign,check.names=FALSE),file.path(out,'r0_assigned_wol_slots.tsv'))
  tsave(data.frame(slot=seq_len(n),wol_tip=colnames(HP)),file.path(out,'wol_slot_key.tsv'))
  hashes<-list(host_tree_sha256=sha(d$paths$host_tree),wol_tree_sha256=sha(d$paths$wol_tree),links_sha256=sha(d$paths$links))
  base<-data.frame(scenario_id=d$id,host_distance_unit=d$spec$host_distance_unit,wol_distance_unit=d$spec$wol_distance_unit,n_links=n,nperm=nperm,seed=seed,correction='cailliez',symmetric=FALSE,host_tree_sha256=hashes$host_tree_sha256,wol_tree_sha256=hashes$wol_tree_sha256,links_sha256=hashes$links_sha256)
  tsave(cbind(base,paco_ss=ss,paco_p=pmc$p,null_model='bijection_host_labels',p_rule=pmc$rule,tail='lower',extreme_count=pmc$count,comparison_tolerance=pmc$tolerance,p_resolution=1/(nperm+1),MCSE_approx=sqrt(pmc$p*(1-pmc$p)/(nperm+1))),file.path(out,'paco_global.tsv'))
  tsave(cbind(base,ParaFitGlobal=pf,p_global=fmc$p,null_model='bijection_host_labels',p_rule=fmc$rule,tail='upper',extreme_count=fmc$count,comparison_tolerance=fmc$tolerance,p_resolution=1/(nperm+1),MCSE_approx=sqrt(fmc$p*(1-fmc$p)/(nperm+1))),file.path(out,'parafit_global.tsv'))
  tsave(cbind(base,paco_ss=ss,paco_p=rmc$p,legacy_package_p=sum(rSS<=ss)/nperm,null_model='vegan_r0',p_rule=rmc$rule,tail='lower',extreme_count=rmc$count,comparison_tolerance=rmc$tolerance,p_resolution=1/(nperm+1),MCSE_approx=sqrt(rmc$p*(1-rmc$p)/(nperm+1))),file.path(out,'paco_r0_global.tsv'))
  result<-list(status='complete',scenario_id=d$id,scenario_index=d$spec$scenario_index,host_distance_unit=d$spec$host_distance_unit,wol_distance_unit=d$spec$wol_distance_unit,n_links=n,nperm=nperm,seed=seed,p_resolution=1/(nperm+1),MCSE_approx=list(paco_bijection=sqrt(pmc$p*(1-pmc$p)/(nperm+1)),parafit_bijection=sqrt(fmc$p*(1-fmc$p)/(nperm+1)),paco_r0=sqrt(rmc$p*(1-rmc$p)/(nperm+1))),input_hashes=hashes,paco_ss=ss,raw_residual_squared_sum=sum(raw^2),paco_bijection_p=pmc$p,parafit_global=pf,parafit_bijection_p=fmc$p,paco_r0_p=rmc$p,paco_r0_legacy_package_p=sum(rSS<=ss)/nperm,paco_r0_occupied_wol_range=range(rWcount),paco_residual_kind='paco_interaction_raw',correction='cailliez',symmetric=FALSE,primary_null='uniform host-label permutations preserving both margins; unrestricted sample.int including occasional observed mapping',sensitivity_null='vegan r0; observed host row sum one is preserved, Wol column sums are free',parafit_default_null_used=FALSE,link_significance_computed=FALSE,paco_coordinate_notes=list(host=hc$note,wol=wc$note),parafit_coordinate_notes=list(host=ah$note,wol=aw$note))
  jsave(result,file.path(out,'summary.json'));result
}
main <- function(){
  args<-commandArgs(trailingOnly=TRUE);values<-list();flags<-c('--synthetic','--check-only')
  i<-1L;while(i<=length(args)){key<-args[i];if(key%in%flags){values[[key]]<-TRUE;i<-i+1L}else{need(key%in%c('--inputs-json','--out-dir')&&i<length(args),'Use --inputs-json MANIFEST --out-dir NEW [--synthetic] [--check-only]');values[[key]]<-args[i+1L];i<-i+2L}}
  for(p in names(required_packages)){need(requireNamespace(p,quietly=TRUE),paste('Missing installed package:',p));need(as.character(packageVersion(p))==required_packages[p],paste('Package version differs from audited implementation:',p))}
  need(!is.null(values[['--inputs-json']])&&!is.null(values[['--out-dir']]),'Both --inputs-json and --out-dir are required')
  mp<-normalizePath(values[['--inputs-json']],mustWork=TRUE);manifest<-jsonlite::read_json(mp,simplifyVector=FALSE);synthetic<-isTRUE(values[['--synthetic']])
  need(identical(manifest$status,'complete')&&identical(as.integer(manifest$schema),1L),'Input manifest must be schema 1 status complete')
  need(identical(manifest$data_kind,if(synthetic)'synthetic' else 'real'),'Synthetic flag/data_kind mismatch')
  need(identical(manifest$correction,'cailliez')&&identical(manifest$symmetric,FALSE),'Frozen analysis requires cailliez and symmetric=false')
  nperm<-as.integer(manifest$nperm);seed<-as.integer(manifest$seed)
  need(length(nperm)==1&&!is.na(nperm)&&nperm>0&&nperm<=9999&&(synthetic||nperm==9999),'Real analyses require exactly 9999 permutations')
  need(length(seed)==1&&!is.na(seed)&&seed>0,'A positive fixed RNG seed is required')
  need(scalar_text(manifest$association_status)&&scalar_text(manifest$association_statement),'Explicit source-confirmation status/statement required')
  need(synthetic||manifest$association_status=='author_confirmed_same_reads','Real analysis requires declared author-confirmed same-read provenance; species name alone is insufficient')
  evidence<-record_path(manifest$association_evidence,'source confirmation');need(length(manifest$scenarios)>0,'No scenarios declared')
  ids<-vapply(manifest$scenarios,function(s)s$id,character(1));need(!anyDuplicated(ids),'Duplicate scenario id')
  indices<-vapply(manifest$scenarios,function(s)as.integer(s$scenario_index),integer(1));need(!anyDuplicated(indices)&&all(indices>=0)&&all(as.double(seed)+indices<.Machine$integer.max),'Duplicate/invalid scenario indices or seed overflow')
  self_arg<-commandArgs()[startsWith(commandArgs(),'--file=')];need(length(self_arg)==1,'Unable to bind executing script path');self<-normalizePath(sub('^--file=','',self_arg),mustWork=TRUE)
  frozen<-list(list(path=mp,sha256=sha(mp)),list(path=evidence,sha256=sha(evidence)),list(path=self,sha256=sha(self)))
  covered<-character();checkpoints<-manifest$analysis_checkpoints
  need(synthetic||length(checkpoints)>0,'Real inputs require completed checkpoints')
  for(rec in checkpoints){p<-record_path(rec,'input completion checkpoint');cp<-jsonlite::read_json(p,simplifyVector=FALSE);need(identical(cp$status,'complete'),'Input checkpoint incomplete');frozen[[length(frozen)+1L]]<-list(path=p,sha256=sha(p));for(rec2 in checkpoint_outputs(cp)){f<-record_path(rec2,'completed input');covered<-c(covered,f);frozen[[length(frozen)+1L]]<-list(path=f,sha256=sha(f))}}
  scenarios<-lapply(manifest$scenarios,read_scenario,synthetic=synthetic)
  for(d in scenarios)for(p in unlist(d$paths)){need(synthetic||p%in%covered,paste('Input absent from complete checkpoint outputs:',p));frozen[[length(frozen)+1L]]<-list(path=p,sha256=sha(p))}
  need(synthetic||evidence%in%covered,'Source confirmation must be covered by completed input checkpoint')
  out<-values[['--out-dir']];need(startsWith(out,'/')&&!file.exists(out)&&!dir.exists(out),'Output directory must be an unused absolute path')
  if(isTRUE(values[['--check-only']])){cat(jsonlite::toJSON(list(status='ready',data_kind=manifest$data_kind,scenarios=ids,nperm=nperm,seed=seed),auto_unbox=TRUE),'\n');return(invisible(NULL))}
  need(dir.create(out,recursive=TRUE),'Cannot create new output directory')
  jsave(manifest,file.path(out,'frozen_input_manifest.json'));writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
  functions<-list(paco_prepare=paco::prepare_paco_data,paco_add_pcoord=paco::add_pcoord,paco_PACo=paco::PACo,paco_residuals=paco::residuals_paco,paco_coordpcoa=get('coordpcoa',asNamespace('paco')),vegan_procrustes=vegan::procrustes,ape_parafit=ape::parafit)
  source<-unlist(lapply(names(functions),function(k)c(paste('FUNCTION',k),deparse(functions[[k]]),'')));writeLines(source,file.path(out,'audited_package_function_sources.txt'))
  definitions<-list(status='frozen_before_fit',packages=as.list(required_packages),R_version=R.version.string,script_sha256=sha(self),R_options_digits=7L,nperm=nperm,base_seed=seed,scenario_seed_rule='base_seed + explicit scenario_index, independent of scheduling or manifest subset',thread_environment=as.list(Sys.getenv(c('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'))),RNGkind=c('Mersenne-Twister','Inversion','Rejection'),paco_tail='lower ss',parafit_tail='upper global',primary='one-to-one host-label permutation of observed association matrix, preserving both host and Wol margins; same draws for PACo and ParaFit',paco_sensitivity='vegan r0 preserving row sums only; plus-one p and legacy package count/n both retained',correction='cailliez',symmetric=FALSE,link_tests=FALSE,source_identity_claim=manifest$association_statement,source_hashes=frozen)
  jsave(definitions,file.path(out,'method_frozen_before_fit.json'))
  results<-lapply(scenarios,function(d)run_scenario(d,file.path(out,d$id),nperm,seed+as.integer(d$spec$scenario_index)))
  for(rec in frozen)need(identical(sha(rec$path),rec$sha256),paste('Input changed during analysis:',rec$path))
  outputs<-list.files(out,recursive=TRUE,full.names=TRUE);outputs<-outputs[!dir.exists(outputs)]
  cp<-list(status='complete',data_kind=manifest$data_kind,completed_at_UTC=format(Sys.time(),tz='UTC',usetz=TRUE),script_sha256=sha(self),scenario_count=length(results),scenarios=results,source_hashes=frozen,outputs=lapply(outputs,function(p)list(path=normalizePath(p),sha256=sha(p))))
  jsave(cp,file.path(out,'completed.json'));cat(jsonlite::toJSON(list(status='complete',data_kind=manifest$data_kind,out_dir=out,scenarios=ids),auto_unbox=TRUE),'\n')
}
if(sys.nframe()==0L)tryCatch(main(),error=function(e){cat('COPHYLOGENY_GATE:',conditionMessage(e),'\n',file=stderr());quit(status=2L)})
