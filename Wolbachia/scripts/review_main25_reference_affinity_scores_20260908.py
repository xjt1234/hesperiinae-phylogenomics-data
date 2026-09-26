#!/usr/bin/env python3
"""Describe main25 reference-affinity score gaps, without voting, new searches, or gene removal."""
import collections,csv,datetime,hashlib,json,itertools,statistics
from decimal import Decimal
from pathlib import Path
RUN=Path('/home/data/t200301/xjt/Hesperiinae_wolbachia/Wolbachia_R2_6_QC_rerun_20260907_112808')
BASE=RUN/'03_MATRICES_AND_TREES'
V2=BASE/'supergroup_annotation_review_20260908'
OUT=RUN/'02_QC/target_figure_recovery_20260908/main25_reference_affinity_scores_20260908'
FIELDS='qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore qseq sseq score stitle'.split()
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
def read(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rs,fields=None):
 with p.open('x',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields or list(rs[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rs)
def jwrite(p,o):
 with p.open('x') as f:json.dump(o,f,ensure_ascii=False,indent=2);f.write('\n')
def d(v):return Decimal(str(v))
def ds(v):return str(v) if v is not None else ''
def summary_group(hs):
 if not hs:return {'qualifying_HSPs':0,'qualifying_references':0,'best_bitscore':'','best_raw_score':'','best_bitscore_tie_raw_score_min':'','best_bitscore_tie_raw_score_max':'','best_bitscore_pident_min':'','best_bitscore_pident_max':'','best_bitscore_qcov_min':'','best_bitscore_qcov_max':'','best_bitscore_query_intervals':'','best_bitscore_reference_accessions':'','best_bitscore_HSPs':[]}
 best=max(d(r['bitscore']) for r in hs);ties=[r for r in hs if d(r['bitscore'])==best]
 return {'qualifying_HSPs':len(hs),'qualifying_references':len({r['reference_accession'] for r in hs}),'best_bitscore':ds(best),'best_raw_score':max(int(r['raw_score']) for r in hs),'best_bitscore_tie_raw_score_min':min(int(r['raw_score']) for r in ties),'best_bitscore_tie_raw_score_max':max(int(r['raw_score']) for r in ties),'best_bitscore_pident_min':ds(min(d(r['pident']) for r in ties)),'best_bitscore_pident_max':ds(max(d(r['pident']) for r in ties)),'best_bitscore_qcov_min':ds(min(d(r['query_coverage_exact']) for r in ties)),'best_bitscore_qcov_max':ds(max(d(r['query_coverage_exact']) for r in ties)),'best_bitscore_query_intervals':';'.join(sorted({r['qstart']+'-'+r['qend'] for r in ties})),'best_bitscore_reference_accessions':';'.join(sorted({r['reference_accession'] for r in ties})),'best_bitscore_HSPs':ties}
def main():
 if OUT.exists():raise FileExistsError('Refusing existing audit')
 cp=RUN/'checkpoints/supergroup_annotation_review_20260908.complete.json';hsp=BASE/'supergroup_reference_review_20260908/search/strict889_vs_all244.tsv';oldcp=RUN/'checkpoints/homology_supergroup_reference_review.complete.json';planpath=RUN/'01_PROVENANCE/frozen_tree_inputs.json';predpath=BASE/'coding_loci_rebuild/prediction/gene_provenance_and_qc.tsv'
 inputs=[cp,oldcp,hsp,planpath,predpath,V2/'reference_annotations_v2.tsv',V2/'main25_gene_supergroup_evidence_v2.tsv',V2/'main25_sample_evidence_counts_v2.tsv',Path(__file__).resolve()];hashes={str(p):sha(p) for p in inputs}
 cj=json.loads(cp.read_text());oj=json.loads(oldcp.read_text())
 assert cj['status']=='complete' and hashes[str(cp)]=='7fe0ddec216e19aa169a863c8a72317a84fe1a5e32648312ded129376e7a0d40'
 for p in [V2/'reference_annotations_v2.tsv',V2/'main25_gene_supergroup_evidence_v2.tsv',V2/'main25_sample_evidence_counts_v2.tsv']:assert hashes[str(p)]==cj['output_sha256'][str(p)]
 assert hashes[str(hsp)]==next(r['sha256'] for r in oj['outputs'] if r['path']==str(hsp))
 refs={r['accession']:r for r in read(V2/'reference_annotations_v2.tsv')};genes={r['gene_id']:r for r in read(V2/'main25_gene_supergroup_evidence_v2.tsv')};samples=read(V2/'main25_sample_evidence_counts_v2.tsv');plan=json.loads(planpath.read_text());repeat=set(plan['reference_repeat_families_excluded_in_sensitivity']);mainset={r['taxon'] for r in samples};pred={r['gene_id']:r for r in read(predpath) if r['gene_id'] in genes}
 bgenes={g for g,r in genes.items() if r['evidence_category_v2']=='B_only'};conflicts={r['taxon'] for r in samples if r['both_A_only_and_B_only_distinct_families']=='1'};targets=bgenes|{g for g,r in genes.items() if r['taxon'] in conflicts}
 assert len(genes)==773 and len(mainset)==25 and len(bgenes)==21 and len(conflicts)==5
 by=collections.defaultdict(lambda:collections.defaultdict(list));total=0;qualified=0
 with hsp.open() as f:
  for lineno,line in enumerate(f,1):
   total+=1;v=line.rstrip('\n').split('\t')
   if len(v)!=16:raise ValueError('RawHSPfields')
   raw=dict(zip(FIELDS,v));g=raw['qseqid']
   if g not in genes:continue
   r=genes[g];a=raw['sseqid'];qcov=d(abs(int(raw['qend'])-int(raw['qstart']))+1)/d(raw['qlen'])
   qual=d(raw['pident'])>=75 and qcov>=d('.5') and int(raw['length'])>=120 and d(raw['evalue'])<=d('1e-20')
   if not qual:continue
   qualified+=1
   if raw['stitle']!=refs[a]['title'] or int(raw['slen'])!=int(refs[a]['sequence_length']):raise ValueError('Exactreferencebinding')
   # Keep score and percent-identity text exactly as emitted by BLAST; do not turn a rounded bit-score tie into extra precision.
   rec={'gene_id':g,'cluster_id':r['cluster_id'],'taxon':r['taxon'],'reference_group_v2':refs[a]['supergroup_v2'],'reference_accession':a,'reference_title':raw['stitle'],'BLAST_line':lineno,'bitscore':raw['bitscore'],'raw_score':raw['score'],'pident':raw['pident'],'query_coverage_exact':str(qcov),'query_covered_bp':abs(int(raw['qend'])-int(raw['qstart']))+1,'query_length':raw['qlen'],'alignment_columns':raw['length'],'evalue':raw['evalue'],'qstart':raw['qstart'],'qend':raw['qend'],'sstart':raw['sstart'],'send':raw['send'],'query_hsp_sha256':hashlib.sha256(raw['qseq'].encode()).hexdigest(),'subject_hsp_sha256':hashlib.sha256(raw['sseq'].encode()).hexdigest(),'ungapped_query_hsp_sha256':hashlib.sha256(raw['qseq'].replace('-','').encode()).hexdigest(),'in_frozen_repeat_sensitivity_family':int(r['cluster_id'] in repeat)}
   by[g][refs[a]['supergroup_v2']].append(rec)
 if total!=103957:raise ValueError('Raw rowcount')
 OUT.mkdir();allrows=[];bestrows=[]
 for g,r in sorted(genes.items(),key=lambda kv:(kv[1]['taxon'],kv[1]['cluster_id'])):
  summaries={k:summary_group(by[g][k]) for k in ['A','B','E','unknown']}
  maxscore=max(d(s['best_bitscore']) for s in summaries.values() if s['best_bitscore']);globalgroups={k for k,s in summaries.items() if s['best_bitscore'] and d(s['best_bitscore'])==maxscore}
  if globalgroups!=set(r['top_tie_groups_v2'].split(';')) or maxscore!=d(r['top_reported_bitscore']):raise ValueError('Raw HSP reconfirmed tops differ')
  a=summaries['A'];b=summaries['B'];comparable=bool(a['best_bitscore'] and b['best_bitscore'])
  bitdiff=d(b['best_bitscore'])-d(a['best_bitscore']) if comparable else None;rawdiff=b['best_raw_score']-a['best_raw_score'] if comparable else None
  row={'gene_id':g,'cluster_id':r['cluster_id'],'taxon':r['taxon'],'evidence_category_v2':r['evidence_category_v2'],'top_tie_groups_v2':r['top_tie_groups_v2'],'main25_B_only_gene':int(g in bgenes),'sample_with_Aonly_Bonly_distinct_families':int(r['taxon'] in conflicts),'targeted_detail_gene':int(g in targets),'in_frozen_repeat_sensitivity_family':int(r['cluster_id'] in repeat),'full_CDS_sequence_sha256':pred[g]['cds_sha256'],'B_minus_A_best_bitscore':ds(bitdiff),'B_minus_A_best_raw_score':rawdiff if comparable else '','bitscore_difference_divided_by_query_nt':ds(bitdiff/d(pred[g]['cds_length'])) if comparable else '','A_B_best_query_interval_sets_equal':int(a['best_bitscore_query_intervals']==b['best_bitscore_query_intervals']) if comparable else '','no_A_qualifying_HSP_under_frozen_search':int(not a['best_bitscore']),'no_B_qualifying_HSP_under_frozen_search':int(not b['best_bitscore']),'significance_or_classification_threshold_applied':False,'final_sample_supergroup':'NOT_ASSIGNED'}
  for k,s in summaries.items():
   for key,value in s.items():
    if key!='best_bitscore_HSPs':row[k+'_'+key]=value
   if g in targets:
    for z in s['best_bitscore_HSPs']:bestrows.append({**z,'is_B_only_gene':int(g in bgenes),'score_selection':'all ties at this group maximum reported bitscore','not_jointly_maximizing_identity_or_coverage':True})
  allrows.append(row)
 write(OUT/'main25_all773_gene_group_best_scores.tsv',allrows)
 br=[r for r in allrows if r['main25_B_only_gene']];write(OUT/'main25_B_only21_gene_score_review.tsv',br)
 write(OUT/'main25_conflict_samples_all_gene_score_review.tsv',[r for r in allrows if r['sample_with_Aonly_Bonly_distinct_families']])
 write(OUT/'target_gene_group_best_all_tied_HSPs.tsv',bestrows)
 # Exact sequence duplicates among B-only genes: retain identities, never treat them as independent confirmation.
 seqs=collections.defaultdict(list)
 for r in br:seqs[r['full_CDS_sequence_sha256']].append(r)
 duplicate_rows=[{'CDS_sha256':h,'genes':';'.join(x['gene_id'] for x in rs),'taxa':';'.join(x['taxon'] for x in rs),'families':';'.join(sorted({x['cluster_id'] for x in rs})),'count':len(rs)} for h,rs in seqs.items() if len(rs)>1]
 write(OUT/'B_only_exact_CDS_duplicate_groups.tsv',duplicate_rows,['CDS_sha256','genes','taxa','families','count'])
 classes=[]
 for s in samples:
  taxon=s['taxon'];rs=[r for r in genes.values() if r['taxon']==taxon];nA=int(s['A_only_gene_count']);nB=int(s['B_only_gene_count']);ntie=int(s['A_B_tie_including_unknown_gene_count']);groups=set().union(*(set(r['top_tie_groups_v2'].split(';')) for r in rs))
  unambig='A_only_reference_affinity_in_unambiguous_genes' if nA and not nB else 'B_only_reference_affinity_in_unambiguous_genes' if nB and not nA else 'cross_family_A_B_affinity_disagreement' if nA and nB else 'no_unambiguous_A_or_B_affinity'
  # Conservative optional figure categories are descriptive evidence patterns, not final supergroup diagnoses.
  if nA and nB:color='conflicting_reference_affinity';reason='不同家族分别出现A-only与B-only；保留两方，不多数票定组'
  elif nA and 'B' not in groups:color='A_affinity_no_observed_B';reason='有A-only证据，全部最高分集合未见明确B；仍保留unknown基因，不等于确诊A菌株'
  elif nB and 'A' not in groups:color='B_affinity_no_observed_A';reason='有B-only证据，全部最高分集合未见明确A；仍保留unknown基因，不等于确诊B菌株'
  else:color='unclassified';reason='存在A/B并列或其他未解决证据；不把未并列基因的一侧标签推广为样本确定超群'
  classes.append({'taxon':taxon,'strict_genes':s['strict_genes'],'A_only':nA,'B_only':nB,'AB_tie_without_unknown':int(s['AB_tie_gene_count']),'unknown_only':int(s['unknown_only_gene_count']),'A_unknown_tie':int(s['A_unknown_tie_gene_count']),'B_unknown_tie':int(s['B_unknown_tie_gene_count']),'AB_unknown_tie':int(s['AB_unknown_tie_gene_count']),'AB_tie_including_unknown':ntie,'unknown_including_all_patterns':int(s['unknown_including_gene_count']),'unambiguous_gene_pattern':unambig,'conservative_figure_evidence_category_suggestion':color,'suggestion_reason_zh':reason,'unknown_retained':True,'final_strain_or_sample_supergroup':'UNCLASSIFIED_NOT_ASSIGNED','use_legend':'Reference-affinity evidence pattern; not strain supergroup diagnosis'})
 write(OUT/'main25_sample_descriptive_affinity_categories.tsv',classes)
 family=[]
 for c in sorted({r['cluster_id'] for r in br}):
  rs=[r for r in br if r['cluster_id']==c]
  family.append({'cluster_id':c,'B_only_genes':len(rs),'taxa':';'.join(r['taxon'] for r in rs),'in_frozen_repeat_sensitivity_family':int(c in repeat),'B_minus_A_bitscore_values':';'.join(sorted(set(r['B_minus_A_best_bitscore'] for r in rs),key=d)),'B_minus_A_raw_score_values':';'.join(map(str,sorted(set(r['B_minus_A_best_raw_score'] for r in rs))))})
 write(OUT/'main25_B_only_family_summary.tsv',family)
 for p,h in hashes.items():
  if sha(Path(p))!=h:raise ValueError('Input changed')
 bitdiffs=[d(r['B_minus_A_best_bitscore']) for r in br];rawdiffs=[r['B_minus_A_best_raw_score'] for r in br];rep_b=[r for r in br if r['in_frozen_repeat_sensitivity_family']]
 persist=[]
 for t in conflicts:
  rs=[r for r in allrows if r['taxon']==t and not r['in_frozen_repeat_sensitivity_family']]
  if any(r['evidence_category_v2']=='A_only' for r in rs) and any(r['evidence_category_v2']=='B_only' for r in rs):persist.append(t)
 summary={'status':'complete_descriptive_affinity_score_review','created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'main25_genes':len(genes),'targeted_detail_genes':len(targets),'main25_B_only_genes':len(br),'main25_B_only_families':len(family),'main25_B_only_samples':len({r['taxon'] for r in br}),'raw_HSP_rows_read':total,'main25_qualifying_HSPs':qualified,'group_best_tied_HSP_detail_rows':len(bestrows),'B_minus_A_bitscore_range':[str(min(bitdiffs)),str(max(bitdiffs))],'B_minus_A_bitscore_median':str(statistics.median(bitdiffs)),'B_minus_A_raw_score_range':[min(rawdiffs),max(rawdiffs)],'B_minus_A_bitscore_observed_values_counts':dict(collections.Counter(str(x) for x in bitdiffs)),'B_only_genes_with_A_and_B_reported_qualifying_hits':sum(bool(r['A_best_bitscore'] and r['B_best_bitscore']) for r in br),'B_only_genes_in_frozen_repeat_excluded_families':len(rep_b),'B_only_families_in_frozen_repeat_excluded_set':sorted({r['cluster_id'] for r in rep_b}),'five_cross_family_disagreements_still_present_after_existing_repeat_sensitivity_filter':sorted(persist),'exact_duplicate_CDS_groups_among_B_only_genes':len(duplicate_rows),'unambiguous_gene_pattern_counts':dict(collections.Counter(r['unambiguous_gene_pattern'] for r in classes)),'conservative_optional_figure_category_counts':dict(collections.Counter(r['conservative_figure_evidence_category_suggestion'] for r in classes)),'figure_categories_are_reference_affinity_patterns_not_strain_supergroups':True,'classification_score_threshold_set':False,'genes_dropped':0,'sample_supergroup_assignments':0,'new_BLAST_or_tree_runs':0,'input_sha256':hashes}
 jwrite(OUT/'summary.json',summary)
 txt=['# 主25 B亲和基因分差与着色解释（2026-09-08）','','当前主25共有21个B-only基因，分布于9个家族和8个样本；每个均有满足冻结条件的A和B参考命中。B-only在此只表示“该基因的最高报告bitscore参考标签集合为B”，不代表样本/菌株确诊。','',f"这些基因B−A的最高报告bitscore差为 {min(bitdiffs)}–{max(bitdiffs)}（中位数{statistics.median(bitdiffs)}），原始score差为{min(rawdiffs)}–{max(rawdiffs)}。不能把全部差异概括为非常小或全部稳健；本次逐值报告，未设置按结果选择的分差阈值。",'','| B−A bitscore | 基因数 |','|---:|---:|']
 for val,n in sorted(collections.Counter(bitdiffs).items()):txt.append(f'| {val} | {n} |')
 txt+=['',f"21个B-only中仅{len(rep_b)}个基因属于已冻结repeat敏感性排除家族（{';'.join(summary['B_only_families_in_frozen_repeat_excluded_set'])}）；剩余不在该7家族名单。五个跨家族A/B分歧样本在按既有名单暂看非repeat家族证据后仍全部存在分歧，因此不能仅以repeat解释全部现象。没有新增排除操作。",'','主要集中模式：Acerbas、Aeromachus nana、Sovia subflava、Thoressa kuata的P1b45d62aeeff39a57fb5f62e最高B−A约84 bits；Aeromachus nana/Thoressa kuata的P23e7a20ecc5930faafdadf7b约140 bits；Polytremis eltola两个B-only家族的差为9和32 bits。全部精确accession、rawscore、identity、qcov与query区段保存在21基因表和group-best HSP表，未选择一条并列参考代替全体。','', '各组的pident和qcov列来自该组最高bitscore的全部并列HSP，使用min/max范围；它们不是分别从不同命中挑最大identity、最大coverage后拼成一条不存在的优良命中。raw score另报全组最大值和最高bitscore并列中的范围，保留bitscore输出舍入现象；A/B最佳命中覆盖区段不同时明确标记。','', '同一样本多家族分歧只说明参考亲和证据不一致，无法单凭此判断混合感染、重组、污染、整合或某个确定菌株。B-only跨样本也有完全相同CDS，表中保留SHA重复组；这些相同序列不是独立生物学确认。','', '建议将最终图的颜色明确改为“参考亲和证据类别”，不要写作已经确定的超群分型。为避免以少数未并列基因推广全样本，给出保守的四类描述性建议：','- A_affinity_no_observed_B：至少一个A-only家族，所有最高分参考集合未见明确B；unknown继续保留。','- B_affinity_no_observed_A：对称定义；本主25没有满足者。','- conflicting_reference_affinity：不同家族同时有A-only和B-only，保留冲突。','- unclassified：存在A/B并列等未解决证据，不强制选一方。','',f"这套可选着色建议计数为{summary['conservative_optional_figure_category_counts']}。它不写入最终样本supergroup，也不按分差设阈值，仍需要图例解释来源和未解决部分。",'','另附较宽松的“仅看未并列基因”模式（17侧A、3侧B、5跨家族分歧），只供证据梳理，不建议将其直接当最终组别：Isoteinon只有1个B-only，其余14个含A/B并列；Potanthus trachala和Sovia lucasii各6个B-only、9个含A/B并列。因此本报告的保守图示将这3个保持unclassified。','', 'PACo/ParaFit全局分析本身不依赖以上标签；D/E解释也不能从参考亲和类别直接推出宿主转换。全部773个主25基因的A/B/E/unknown最优HSP分数已汇总，21个B-only与5个分歧样本全部基因另有明确子表。未新增BLAST、未建树、未更改任何基因纳入规则。']
 with (OUT/'interpretation_zh.md').open('x') as f:f.write('\n'.join(txt)+'\n')
 qa={'status':'PASS_machine_raw_score_rejoin','checked_input_hashes':len(hashes),'raw_rows103957_verified':True,'all773_gene_top_group_sets_and_bitscores_reproduced_from_raw_HSPs':True,'B_only21_gene_coverage_complete':True,'all5_conflict_samples_all_genes_included':True,'all_group_best_ties_preserved':True,'no_gene_exclusion_or_sample_diagnosis':True,'output_sha256':{str(p):sha(p) for p in sorted(OUT.iterdir()) if p.is_file()}}
 jwrite(OUT/'QA.json',qa)
 print(json.dumps({k:v for k,v in summary.items() if k!='input_sha256'},ensure_ascii=False,indent=2))
if __name__=='__main__':main()
