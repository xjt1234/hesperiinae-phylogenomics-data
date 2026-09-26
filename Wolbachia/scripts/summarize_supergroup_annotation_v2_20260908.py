#!/usr/bin/env python3
"""Apply eight reviewed literature labels to frozen local HSPs; evidence only, never assign samples."""
import argparse,collections,csv,datetime,hashlib,json
from decimal import Decimal
from pathlib import Path
PROPOSAL_SHA='d7221d0ce0bed0220e9c1ec0812e2063515e9b6c263d23d0540d04cc8d016417'
SOURCE_CP_SHA='2927e08a27573879773f57a81c8490f698cfb1d946ad6edfa75e4e26091b634d'
APPROVED={'NZ_CP041215.1':'A','CP003884.1':'A','NZ_CP054598.1':'A','NZ_CP054557.1':'A','CP042904.1':'A','NZ_CP042446.1':'A','NZ_CP011148.1':'A','NC_010981.1':'B'}
CATEGORIES=['A_only','B_only','AB_tie','unknown_only','A_unknown_tie','B_unknown_tie','AB_unknown_tie','E_only','other_explicit_tie','other_explicit_unknown_tie','no_qualifying_HSP']
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
def tsv(p):
 with Path(p).open() as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rows,fields=None):
 with p.open('x',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields or list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def jwrite(p,obj):
 with p.open('x') as f:json.dump(obj,f,ensure_ascii=False,indent=2);f.write('\n')
def category(gs):
 s=frozenset(gs)
 return {frozenset():'no_qualifying_HSP',frozenset(['A']):'A_only',frozenset(['B']):'B_only',frozenset(['A','B']):'AB_tie',frozenset(['unknown']):'unknown_only',frozenset(['A','unknown']):'A_unknown_tie',frozenset(['B','unknown']):'B_unknown_tie',frozenset(['A','B','unknown']):'AB_unknown_tie',frozenset(['E']):'E_only'}.get(s,'other_explicit_unknown_tie' if 'unknown' in s else 'other_explicit_tie')
def counts(rows,col='evidence_category_v2'):
 c=collections.Counter(r[col] for r in rows)
 return {k:c[k] for k in CATEGORIES}
def signature(r):
 return tuple(r[k] for k in ['gene_id','cluster_id','taxon','reference_accession','bitscore','raw_score','pident','query_coverage','qstart','qend','sstart','send'])
def sample_summaries(gene_rows,mainset):
 samples=collections.defaultdict(list)
 for r in gene_rows:samples[r['taxon']].append(r)
 out=[]
 for taxon,rs in sorted(samples.items()):
  c=counts(rs);families={r['cluster_id'] for r in rs}
  if len(families)!=len(rs):raise ValueError('More than one gene per sample/family')
  A={r['cluster_id'] for r in rs if r['evidence_category_v2']=='A_only'};B={r['cluster_id'] for r in rs if r['evidence_category_v2']=='B_only'}
  unknown=[r for r in rs if 'unknown' in r['top_tie_groups_v2'].split(';')];ab=[r for r in rs if {'A','B'}<=set(r['top_tie_groups_v2'].split(';'))]
  status='A_only_and_B_only_evidence_in_distinct_families' if A and B else 'A_B_tied_reference_evidence_present' if ab else 'residual_unknown_reference_evidence_present' if unknown else 'no_observed_A_B_or_unknown_top_ties'
  out.append({'taxon':taxon,'in_main25':int(taxon in mainset),'strict_genes':len(rs),'strict_families':len(families),**{k+'_gene_count':c[k] for k in CATEGORIES},'unknown_including_gene_count':len(unknown),'A_B_tie_including_unknown_gene_count':len(ab),'A_only_family_count':len(A),'B_only_family_count':len(B),'A_only_family_ids':';'.join(sorted(A)),'B_only_family_ids':';'.join(sorted(B)),'both_A_only_and_B_only_distinct_families':int(bool(A and B)),'evidence_pattern':status,'final_sample_supergroup':'NOT_ASSIGNED','interpretation_scope':'highest reported bitscore reference evidence; no majority label or infection/transfer inference'})
 return out
def family_summaries(rows):
 by=collections.defaultdict(list)
 for r in rows:by[r['cluster_id']].append(r)
 result=[]
 for family,rs in sorted(by.items()):
  c=counts(rs)
  result.append({'cluster_id':family,'genes':len(rs),'taxa':len({r['taxon'] for r in rs}),**{k+'_gene_count':c[k] for k in CATEGORIES},'A_only_taxa':';'.join(sorted(r['taxon'] for r in rs if r['evidence_category_v2']=='A_only')),'B_only_taxa':';'.join(sorted(r['taxon'] for r in rs if r['evidence_category_v2']=='B_only')),'contains_A_only_and_B_only_across_taxa':int(c['A_only']>0 and c['B_only']>0),'not_independent_sample_replicates':True})
 return result
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--run-dir',type=Path,required=True);ap.add_argument('--stage',choices=['all'],default='all');args=ap.parse_args();run=args.run_dir.resolve()
 out=run/'03_MATRICES_AND_TREES/supergroup_annotation_review_20260908';cp=run/'checkpoints/supergroup_annotation_review_20260908.complete.json'
 if out.exists() or cp.exists():raise FileExistsError('Refusing existing output')
 proposal=run/'02_QC/upload_archive_recovery_20260908/reference_annotation_sources.tsv';public=proposal.parent/'reference_public_metadata';manifest=public/'raw_response_manifest.tsv';oldroot=run/'03_MATRICES_AND_TREES/supergroup_reference_review_20260908';sourcecp=run/'checkpoints/homology_supergroup_reference_review.complete.json'
 if sha(proposal)!=PROPOSAL_SHA or sha(sourcecp)!=SOURCE_CP_SHA:raise ValueError('Frozen proposal/source checkpoint mismatch')
 checkpoint=json.loads(sourcecp.read_text())
 if checkpoint['status']!='complete':raise ValueError('Incomplete HSP stage')
 # Every prior output is rehashed, including the raw search output and immutable query input.
 oldoutputs={r['path']:r['sha256'] for r in checkpoint['outputs']}
 for path,h in oldoutputs.items():
  if sha(path)!=h:raise ValueError('Source output changed '+path)
 refpath=oldroot/'inputs/reference_annotations.tsv';hsppath=oldroot/'evidence/all_reported_hsp_evidence.tsv';toppath=oldroot/'evidence/all_top_bitscore_ties.tsv';genepath=oldroot/'evidence/gene_supergroup_evidence.tsv';memberpath=run/'03_MATRICES_AND_TREES/coding_orthology_qc/strict_candidate_members.tsv';mainplan=run/'01_PROVENANCE/frozen_tree_inputs.json';gap=proposal.parent/'reference_archive_unknown_label_gaps.tsv';assembly=proposal.parent/'reference_archive_local_assembly_metadata.tsv'
 sources=[proposal,sourcecp,manifest,public/'summary.json',refpath,hsppath,toppath,genepath,memberpath,mainplan,gap,assembly,Path(__file__).resolve()]
 rawrecords=tsv(manifest);rawsha={}
 for r in rawrecords:
  p=Path(r['raw_path']);h=sha(p)
  if h!=r['raw_sha256'] or p.stat().st_size!=int(r['raw_bytes']):raise ValueError('Raw primary/metadata mismatch')
  rawsha[str(p)]=h;sources.append(p)
  req=p.with_name(p.name.replace('.response.raw','.request_result.json'))
  if not req.exists():raise ValueError('No saved requestmetadata')
  sources.append(req)
 sources=list(dict.fromkeys(sources));inputhash={str(p):sha(p) for p in sources}
 prop=tsv(proposal);proposed={r['local_reference_accession']:r for r in prop}
 if len(proposed)!=31:raise ValueError('Expected priority31 proposals')
 accepted={a:r['proposed_supergroup'] for a,r in proposed.items() if r['proposed_supergroup']!='unknown'}
 if accepted!=APPROVED:raise ValueError('Approved exact accession/group set changed')
 refs={r['accession']:r for r in tsv(refpath)}
 if len(refs)!=244:raise ValueError('Expected244references')
 for a,r in proposed.items():
  if refs[a]['sequence_sha256']!=r['local_sequence_sha256'] or refs[a]['explicit_group']!=r['previous_label']:raise ValueError('Proposal reference identity mismatch')
  if r['source_gap_table_sha256']!=inputhash[str(gap)] or r['source_assembly_table_sha256']!=inputhash[str(assembly)]:raise ValueError('Proposal accession/assembly source changed')
  if rawsha[r['metadata_raw_path']]!=r['metadata_raw_sha256']:raise ValueError('Metadata binding mismatch')
  if a in APPROVED:
   if r['review_status']!='SOURCE_SUPPORTED_PROPOSAL_AWAITING_ROOT_REVIEW' or rawsha[r['primary_raw_path']]!=r['primary_raw_sha256']:raise ValueError('Approved source changed')
   if refs[a]['explicit_group']!='unknown':raise ValueError('Attempt to override an existing explicit title label')
 members={r['gene_id']:r for r in tsv(memberpath)};oldgenes={r['gene_id']:r for r in tsv(genepath)};oldtop=tsv(toppath)
 if len(members)!=889 or set(oldgenes)!=set(members) or len(oldtop)!=10914:raise ValueError('Frozen gene/top counts differ')
 mainset=set(next(r['expected_taxa'] for r in json.loads(mainplan.read_text())['scenarios'] if r['label']=='main'))
 if len(mainset)!=25:raise ValueError('Expectedmain25')
 groupv2={a:APPROVED.get(a,r['explicit_group']) for a,r in refs.items()}
 out.mkdir()
 freeze={'schema':1,'status':'frozen_before_annotation_summary','created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'proposal_path':str(proposal),'proposal_sha256':PROPOSAL_SHA,'approved_exact_accession_groups':APPROVED,'review_record':{'reviewer':'root analysis agent','review_scope':'parent explicitly reviewed proposal and five primary XML originals before authorizing this version','approval_basis':'exact accession/strain bridge plus cited primary source, not same host name, strain suffix, or majority hits'},'source_hsp_checkpoint':str(sourcecp),'source_hsp_checkpoint_sha256':SOURCE_CP_SHA,'primary_and_metadata_raw_manifest':str(manifest),'all_input_sha256':inputhash,'raw_primary_and_metadata_sha256':rawsha,'old_output_hash_assertions':len(oldoutputs),'search_rerun':False,'original_panel_or_results_modified':False}
 jwrite(out/'frozen_annotation_inputs.json',freeze)
 reference_rows=[]
 for a,r in sorted(refs.items()):
  p=proposed.get(a,{});updated=a in APPROVED
  reference_rows.append({'accession':a,'title':r['title'],'sequence_length':r['sequence_length'],'sequence_sha256':r['sequence_sha256'],'supergroup_v1_title_only':r['explicit_group'],'supergroup_v2':groupv2[a],'annotation_basis_v2':'reviewed_primary_literature_exact_accession_or_documented_strain_bridge' if updated else 'retained_original_explicit_title' if r['explicit_group']!='unknown' else 'unknown_no_approved_explicit_source','changed_from_v1':int(updated),'proposal_path':str(proposal) if a in proposed else '','proposal_sha256':PROPOSAL_SHA if a in proposed else '','primary_url':p.get('primary_url','') if updated else '','primary_locator':p.get('primary_locator','') if updated else '','evidence_paraphrase':p.get('evidence_text','') if updated else '','accession_strain_bridge':p.get('accession_strain_bridge','') if updated else '','primary_raw_path':p.get('primary_raw_path','') if updated else '','primary_raw_sha256':p.get('primary_raw_sha256','') if updated else '','metadata_raw_path':p.get('metadata_raw_path',''),'metadata_raw_sha256':p.get('metadata_raw_sha256',''),'source_review_status_v2':'accepted_by_root_review_for_reference_only' if updated else 'unchanged','sample_group_assigned':False})
 write(out/'reference_annotations_v2.tsv',reference_rows)
 hsps=tsv(hsppath);bygene=collections.defaultdict(list);allgroupqual=collections.defaultdict(lambda:collections.defaultdict(list));hspchanged=0
 if len(hsps)!=103957:raise ValueError('Expected103957HSPs')
 for r in hsps:
  g=r['gene_id'];m=members[g];a=r['reference_accession']
  if (r['cluster_id'],r['taxon'])!=(m['cluster_id'],m['taxon']) or r['explicit_reference_group']!=refs[a]['explicit_group'] or r['exact_query_and_reference_tracebacks']!='1':raise ValueError('HSP immutable provenance differs')
  qual=Decimal(r['pident'])>=75 and Decimal(r['query_coverage'])>=Decimal('.5') and int(r['length'])>=120 and Decimal(r['evalue'])<=Decimal('1e-20')
  if qual!=bool(int(r['qualifying_hsp'])):raise ValueError('HSP qualifying status differs')
  r['reference_group_v1']=r['explicit_reference_group'];r['reference_group_v2']=groupv2[a];r['annotation_changed_v2']=int(a in APPROVED);hspchanged+=r['annotation_changed_v2']
  if qual:bygene[g].append(r);allgroupqual[g][groupv2[a]].append(r)
 write(out/'all_reported_hsp_evidence_v2.tsv',hsps)
 top=[];genelist=[]
 for g,m in sorted(members.items()):
  hs=bygene[g];best=max((Decimal(r['bitscore']) for r in hs),default=None);tie=[r for r in hs if Decimal(r['bitscore'])==best]
  oldset=sorted({r['reference_group_v1'] for r in tie});newset=sorted({r['reference_group_v2'] for r in tie})
  if ';'.join(oldset)!=oldgenes[g]['top_tie_explicit_groups'] or best!=Decimal(oldgenes[g]['top_reported_bitscore']):raise ValueError('Recomputed original tops differ')
  for r in tie:top.append({k:r[k] for k in ['gene_id','cluster_id','taxon','reference_accession','reference_group_v1','reference_group_v2','annotation_changed_v2','bitscore','raw_score','pident','query_coverage','qstart','qend','sstart','send','BLAST_line']})
  row={'gene_id':g,'cluster_id':m['cluster_id'],'taxon':m['taxon'],'in_main25':int(m['taxon'] in mainset),'qualifying_HSPs':len(hs),'top_reported_bitscore':str(best) if best is not None else '','top_tied_HSP_count':len(tie),'top_tie_reference_accessions':';'.join(sorted({r['reference_accession'] for r in tie})),'top_tie_groups_v1':';'.join(oldset),'top_tie_groups_v2':';'.join(newset),'evidence_category_v1':category(oldset),'evidence_category_v2':category(newset),'top_group_set_changed':int(oldset!=newset),'category_changed':int(category(oldset)!=category(newset)),'remaining_unknown_top_reference_accessions':';'.join(sorted({r['reference_accession'] for r in tie if r['reference_group_v2']=='unknown'})),'A_and_B_tied_within_gene':int('A' in newset and 'B' in newset),'unknown_in_top_ties':int('unknown' in newset),'final_sample_supergroup':'NOT_ASSIGNED'}
  for label in ['A','B','E','unknown']:
   arr=allgroupqual[g][label];row[label+'_best_qualifying_reported_bitscore']=str(max((Decimal(x['bitscore']) for x in arr),default=''));row[label+'_qualifying_reference_count']=len({x['reference_accession'] for x in arr})
  genelist.append(row)
 if collections.Counter(signature(r) for r in top)!=collections.Counter(signature(r) for r in oldtop):raise ValueError('Original top HSP membership or scores changed')
 write(out/'all_top_bitscore_ties_v2.tsv',top)
 write(out/'gene_supergroup_evidence_v2.tsv',genelist)
 maingenes=[r for r in genelist if r['in_main25']];write(out/'main25_gene_supergroup_evidence_v2.tsv',maingenes)
 sample=sample_summaries(genelist,mainset);mainsample=[r for r in sample if r['in_main25']]
 write(out/'sample_evidence_counts_v2.tsv',sample);write(out/'main25_sample_evidence_counts_v2.tsv',mainsample)
 fam=family_summaries(genelist);mainfam=family_summaries(maingenes);write(out/'family_evidence_summary_v2.tsv',fam);write(out/'main25_family_evidence_summary_v2.tsv',mainfam)
 conflict={r['taxon'] for r in sample if r['both_A_only_and_B_only_distinct_families']};conflictevidence=[r for r in genelist if r['taxon'] in conflict and r['evidence_category_v2'] in ['A_only','B_only']]
 write(out/'sample_cross_family_A_B_evidence.tsv',conflictevidence,list(genelist[0]))
 changed=[r for r in genelist if r['top_group_set_changed']];write(out/'gene_annotation_changes_v1_to_v2.tsv',changed,list(genelist[0]))
 for p,h in inputhash.items():
  if sha(p)!=h:raise ValueError('Frozen input changed during stage')
 for p,h in oldoutputs.items():
  if sha(p)!=h:raise ValueError('Old output changed during stage')
 if len(sample)!=37 or len(fam)!=52 or len(mainsample)!=25 or len(mainfam)!=52 or len(maingenes)!=773:raise ValueError('Expected889/52/37/main25/773')
 summary={'schema':1,'status':'complete_reference_annotation_v2_evidence_only','finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'reference_records':244,'references_reannotated':8,'reference_groups_v1':dict(collections.Counter(r['explicit_group'] for r in refs.values())),'reference_groups_v2':dict(collections.Counter(groupv2.values())),'genes':len(genelist),'families':len(fam),'taxa':len(sample),'main25_genes':len(maingenes),'main25_families':len(mainfam),'main25_taxa':len(mainsample),'reported_HSPs':len(hsps),'qualifying_HSPs':sum(len(x) for x in bygene.values()),'annotated_HSP_rows_changed':hspchanged,'top_HSPs':len(top),'top_HSP_membership_and_scores_unchanged':True,'all_genes_top_bitscore_and_reference_set_unchanged':True,'genes_with_changed_top_group_sets':len(changed),'gene_categories_v1':counts(genelist,'evidence_category_v1'),'gene_categories_v2':counts(genelist),'main25_gene_categories_v1':counts(maingenes,'evidence_category_v1'),'main25_gene_categories_v2':counts(maingenes),'genes_with_unknown_top_ties_v1':sum('unknown' in r['top_tie_groups_v1'].split(';') for r in genelist),'genes_with_unknown_top_ties_v2':sum(r['unknown_in_top_ties'] for r in genelist),'samples_with_Aonly_and_Bonly_distinct_family_evidence':sorted(conflict),'main25_samples_with_Aonly_and_Bonly_distinct_family_evidence':sorted(conflict&mainset),'samples_with_within_gene_AB_ties':sorted({r['taxon'] for r in genelist if r['A_and_B_tied_within_gene']}),'sample_final_assignments':0,'majority_vote_used':False,'BLAST_rerun':False,'original_panel_modified':False,'statistical_significance_tests_performed':False,'scope_limit':'Reported BLAST similarity evidence with reviewed reference labels, not strain identification, orthology proof, mixed infection proof or transfer inference','input_sha256':inputhash}
 jwrite(out/'summary.json',summary)
 txt=['# 参考超群注释v2与逐基因证据汇总（2026-09-08）','','已按主线程对5份一手原文及accession桥接的审查，将8条参考的unknown标签补为7条A、1条B。原244条参考序列、原BLAST输出及旧证据表均未修改；本阶段不联网、不重新比对。','',f"新参考注释为 {summary['reference_groups_v2']}。全部103,957条HSP、103,942条合格HSP及10,914条最高bitscore并列HSP都保留原分数和成员身份。889基因、52家族、37样本完整；主25含773基因和52家族。",'',f"旧版含unknown最高分并列的基因 {summary['genes_with_unknown_top_ties_v1']}，新版本为 {summary['genes_with_unknown_top_ties_v2']}。标签改进并不自动变成样本分组：每一基因仍保留全部最高分参考的标签集合。",'', '| 逐基因最高分参考类别 | 原37样本集v1 | 原37样本集v2 | 主25 v2 |','|---|---:|---:|---:|']
 for c in CATEGORIES:txt.append(f"| {c} | {summary['gene_categories_v1'][c]} | {summary['gene_categories_v2'][c]} | {summary['main25_gene_categories_v2'][c]} |")
 txt+=['','A_only/B_only表示该基因最高分并列参考只含已注释A/只含已注释B；AB_tie表示同一基因最高分并列包含A和B；A_unknown/B_unknown/AB_unknown分别保留未知参考参与的并列。unknown不是Other，也不是某个确定超群。','',f"37样本中有 {len(conflict)} 个样本同时在不同严格家族看到A_only与B_only证据；其中主25有 {len(conflict&mainset)} 个。样本名：{', '.join(sorted(conflict)) or '无'}。",'',f"同一基因出现A/B最高分并列的样本：{', '.join(summary['samples_with_within_gene_AB_ties']) or '无'}。",'', '不同家族的A/B证据不一致，与同一基因的A/B并列，是两类需区分的现象。它们可受保守同源区段、参考取样、重组/多拷贝等因素影响；本表仅记录现象，不足以判定混合感染、确定菌株、宿主转换或污染。没有按多数票给样本强制A/B；原树高支持也不能替代这些检查。','','每个sample表同时列出完整互斥类别计数、含unknown总数、跨家族A_only/B_only家族ID；这些家族计数是描述性证据，不作为独立重复做显著性检验。52家族表显示各家族跨样本的证据分布，sample_cross_family_A_B_evidence.tsv保留具体基因以便下一步定位。','', '应用范围和原文SHA在frozen_annotation_inputs.json；reference_annotations_v2.tsv逐条记录原标签、新标签、证据类型、原文位置与accession桥接。完成checkpoint和QA包含输入/输出哈希，不覆盖旧输出或修改原始FASTA。']
 with (out/'interpretation_zh.md').open('x') as f:f.write('\n'.join(txt)+'\n')
 qa={'status':'PASS_machine_annotation_v2_evidence_only','checked_upstream_output_hashes':len(oldoutputs),'checked_frozen_inputs':len(inputhash),'primary_and_metadata_raw_responses':len(rawsha),'primary_publications':len({proposed[a]['primary_raw_path'] for a in APPROVED}),'accepted_exact_accessions':APPROVED,'no_previously_known_label_overridden':True,'all244_sequence_hashes_unchanged':True,'all103957_HSP_provenance_and_qualification_verified':True,'all10914_top_HSP_multiset_unchanged':True,'counts889_52_37_and_main25_773_verified':True,'all_gene_category_counts_reconcile':sum(summary['gene_categories_v2'].values())==889,'unknown_not_collapsed_to_Other':True,'final_sample_assignments':0,'mixed_infection_not_inferred':True,'output_sha256':{str(p):sha(p) for p in sorted(out.iterdir()) if p.is_file()}}
 jwrite(out/'QA.json',qa)
 outputhash={str(p):sha(p) for p in sorted(out.iterdir()) if p.is_file()}
 jwrite(cp,{'schema':1,'status':'complete','stage':'supergroup_annotation_review_20260908','completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'summary':summary,'input_sha256':inputhash,'output_sha256':outputhash,'outputs':[{'path':p,'sha256':h} for p,h in outputhash.items()],'final_sample_assignments':0})
 print(json.dumps({k:v for k,v in summary.items() if k!='input_sha256'},ensure_ascii=False,indent=2));print('checkpoint',cp);print('checkpoint_sha256',sha(cp))
if __name__=='__main__':main()
