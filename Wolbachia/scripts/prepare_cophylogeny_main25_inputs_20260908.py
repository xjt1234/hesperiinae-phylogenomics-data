#!/usr/bin/env python3
"""Prepare source-confirmed AHE25 and explicitly historical USCO25 inputs, with exact Decimal distances."""
import argparse,collections,copy,csv,datetime,hashlib,itertools,json,re
from decimal import Decimal,getcontext
from pathlib import Path
getcontext().prec=50
BASE=Path('/home/data/t200301/xjt')
A_TREE=BASE/'Hesperiinae_review/treePL_R2_3_T25_deep_only_20260904_000707/output/T25_deep_only_dated.tre'
A_SHA='28d4460885cc7425c9c703b79aff43c24a2ab931c86110c5ffee8d6533309eba'
AUTHOR_SHA='d2bb0a2bd2c9fea217e2e43f4894fcd93f07cf061a877bc168a00c0ba984bf77'
OLD_AHE=BASE/'AHE_Hesperiinae/AHE_NT_supermatrix/Hesperiinae_NT_supermatrix.fasta'
FORMAL_AHE=BASE/'deep_Confilcts/Hesperiinae_NT_supermatrix.fasta'
LEG=BASE/'Wolchbia/cophylogeny_analysis_20260226_paco_v13'
ZERO=Decimal(0)
class Node:
 def __init__(self,children=None,label='',length=None):self.children=children or [];self.label=label;self.length=length;self.source_id='';self.path_ids=[]
def parse_newick(text):
 # This stage accepts plain Newick labels used by frozen inputs; unknown syntax fails closed.
 if any(x in text for x in "[]\"'"):raise ValueError('Unsupported comment/quoted Newick syntax')
 tok=re.findall(r'\(|\)|,|:|;|[^\s(),:;]+',text);i=0
 def node():
  nonlocal i
  children=[]
  if tok[i]=='(':
   i+=1;children.append(node())
   while tok[i]==',':i+=1;children.append(node())
   if tok[i]!=')':raise ValueError('Expected close parenthesis')
   i+=1
  label='';length=None
  if i<len(tok) and tok[i] not in [')',',',':',';']:label=tok[i];i+=1
  if i<len(tok) and tok[i]==':':
   i+=1;length=Decimal(tok[i]);i+=1
   if not length.is_finite() or length<0:raise ValueError('Invalid tree branch')
  if not children and not label:raise ValueError('Unnamed tip')
  return Node(children,label,length)
 root=node()
 if tok[i:]!=[';']:raise ValueError('Expected one complete tree')
 for k,n in enumerate(preorder(root),1):n.source_id=f'N{k:04d}';n.path_ids=[n.source_id]
 tips=[n.label for n in preorder(root) if not n.children]
 if len(tips)!=len(set(tips)):raise ValueError('Duplicate tree tips')
 return root
def preorder(root):
 yield root
 for c in root.children:yield from preorder(c)
def descend(n):
 return frozenset([n.label]) if not n.children else frozenset().union(*(descend(c) for c in n.children))
def decstr(x):return format(x,'f')
def newick(n):
 s=('('+','.join(newick(c) for c in n.children)+')' if n.children else '')+n.label
 return s+(':'+decstr(n.length) if n.length is not None else '')
def prune(root,keep):
 if not keep<=descend(root):raise ValueError('Target absent')
 def rec(n,isroot=False):
  if not n.children:return copy.deepcopy(n) if n.label in keep else None
  children=[q for c in n.children if (q:=rec(c)) is not None]
  if not children:return None
  q=copy.copy(n);q.children=children;q.path_ids=[n.source_id]
  if len(children)==1 and not isroot:
   child=children[0];child.length=(child.length or ZERO)+(n.length or ZERO);child.path_ids=[n.source_id]+child.path_ids
   return child
  return q
 return rec(root,True)
def maps(root):
 parent={root.source_id:None};depth={root.source_id:ZERO};nodes={};tips={}
 def rec(n):
  nodes[n.source_id]=n
  if not n.children:tips[n.label]=n
  for c in n.children:parent[c.source_id]=n.source_id;depth[c.source_id]=depth[n.source_id]+(c.length or ZERO);rec(c)
 rec(root)
 return nodes,parent,depth,tips
def lca(a,b,parent):
 ancestors=set()
 while a is not None:ancestors.add(a);a=parent[a]
 while b not in ancestors:b=parent[b]
 return b
def pair_distance(a,b,parent,depth,tips):
 m=lca(tips[a].source_id,tips[b].source_id,parent)
 return depth[tips[a].source_id]+depth[tips[b].source_id]-2*depth[m],m
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
def table(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def write_tsv(p,rows,fields=None):
 with p.open('x',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields or list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def write_json(p,obj):
 with p.open('x') as f:json.dump(obj,f,ensure_ascii=False,indent=2);f.write('\n')
def fasta_selected(p,needed):
 found={};name=None;seq=[];allids=set()
 def add():
  if name in needed:
   s=''.join(seq).upper();found[name]={'length':len(s),'sha256':hashlib.sha256(s.encode()).hexdigest(),'sequence':s}
 with p.open() as f:
  for line in f:
   if line.startswith('>'):
    add();name=line[1:].strip().split()[0]
    if name in allids:raise ValueError('Duplicate FASTA ID')
    allids.add(name);seq=[]
   elif name in needed:seq.append(line.strip())
 add()
 if set(found)!=needed:raise ValueError('Missing exact FASTA row')
 return found,len(allids)
def tree_artifacts(out,prefix,root,keep,dated=False):
 pruned=prune(root,keep)
 if pruned.source_id!=root.source_id:raise ValueError('Root lost')
 fn,fp,fd,ft=maps(root);pn,pp,pd,pt=maps(pruned)
 anchor=max(fd[n.source_id] for n in ft.values()) if dated else None
 min_height=min(fd[n.source_id] for n in ft.values())
 max_height=max(fd[n.source_id] for n in ft.values())
 pairrows=[]
 for a,b in itertools.combinations(sorted(keep),2):
  d,m=pair_distance(a,b,fp,fd,ft);d2,m2=pair_distance(a,b,pp,pd,pt)
  if d!=d2 or m!=m2:raise ValueError('Pruned distance or pair MRCA changed')
  pairrows.append({'host_tip_a':a,'host_tip_b':b,'full_distance_decimal':decstr(d),'pruned_distance_decimal':decstr(d2),'exact_decimal_equal':1,'full_MRCA_node_id':m,'pruned_MRCA_original_node_id':m2,'MRCA_age_from_full_tree_anchor':decstr(anchor-fd[m]) if dated else ''})
 if len(pairrows)!=300:raise ValueError('Expected300pairs')
 # Mapping includes excluded nodes and collapsed paths; no old internal support label is reinterpreted.
 path_target={sid:n.source_id for n in pn.values() if n.source_id!=pruned.source_id for sid in n.path_ids}
 node_rows=[]
 for sid,n in fn.items():
  full_desc=descend(n);kept=full_desc&keep
  retained=sid in pn
  target=sid if retained else path_target.get(sid,'')
  state='retained_original_root' if sid==root.source_id else 'retained_node' if retained else 'collapsed_into_retained_path' if target else 'excluded_no_selected_descendants'
  node_rows.append({'full_node_id':sid,'full_node_label':n.label,'full_is_tip':int(not n.children),'full_parent_node_id':fp[sid] or '', 'full_incoming_branch_length':decstr(n.length) if n.length is not None else '', 'full_root_distance':decstr(fd[sid]),'full_descendant_count':len(full_desc),'full_descendant_tips':';'.join(sorted(full_desc)),'selected_descendant_count':len(kept),'selected_descendant_tips':';'.join(sorted(kept)),'pruning_status':state,'pruned_target_original_node_id':target,'pruned_incoming_branch_length':decstr(pn[sid].length) if retained and pn[sid].length is not None else '', 'pruned_incoming_full_path_node_ids':';'.join(pn[sid].path_ids) if retained else '', 'full_age_from_full_tree_anchor':decstr(anchor-fd[sid]) if dated else '', 'pruned_age_from_same_anchor':decstr(anchor-pd[sid]) if dated and retained else ''})
  if retained and fd[sid]!=pd[sid]:raise ValueError('Root depth changed')
 # Validate serialization via both exact parser and Bio.Phylo's independent Newick parser.
 treefile=out/(prefix+'.root_preserved.nwk')
 with treefile.open('x') as f:f.write(newick(pruned)+';\n')
 reparsed=parse_newick(treefile.read_text());rn,rp,rd,rt=maps(reparsed)
 for a,b in itertools.combinations(sorted(keep),2):
  if pair_distance(a,b,fp,fd,ft)[0]!=pair_distance(a,b,rp,rd,rt)[0]:raise ValueError('Serialization changed distance')
 from Bio import Phylo
 bio=Phylo.read(treefile,'newick')
 if set(t.name for t in bio.get_terminals())!=keep:raise ValueError('Bio tips differ')
 bio_max=max(abs(float(r['full_distance_decimal'])-bio.distance(r['host_tip_a'],r['host_tip_b'])) for r in pairrows)
 if bio_max>1e-9:raise ValueError('Independent float parser distances differ')
 write_tsv(out/(prefix+'.all300_pair_distances.tsv'),pairrows)
 write_tsv(out/(prefix+'.full_to_pruned_node_map.tsv'),node_rows)
 matrix=[dict(host_tip=a,**{b:decstr(pair_distance(a,b,pp,pd,pt)[0]) for b in sorted(keep)}) for a in sorted(keep)]
 write_tsv(out/(prefix+'.patristic_distance_matrix.tsv'),matrix)
 return {'tree':str(treefile),'tree_sha256':sha(treefile),'source_tips':len(ft),'selected_tips':len(pt),'source_root_original_node_id':root.source_id,'source_root_label':root.label,'pruned_root_label':pruned.label,'original_root_preserved':True,'source_root_child_count':len(root.children),'pruned_root_child_count':len(pruned.children),'single_child_root_preserved':len(pruned.children)==1,'all300_pair_distances_exact_decimal_equal':True,'all300_pair_MRCA_original_node_ids_identical':True,'all_retained_node_depths_and_ages_exact':True,'BioPhylo_independent_max_abs_float_distance_error':bio_max,'full_root_to_tip_min':decstr(min_height),'full_root_to_tip_max':decstr(max_height),'full_root_to_tip_range':decstr(max_height-min_height),'node_age_anchor_full_tree_max_root_tip':decstr(anchor) if dated else None,'age_unit':'Ma as declared by dated tree source' if dated else 'not dated; no age interpretation','rounding_treatment':'No re-ultrametricization; retain original decimals; dated ages use fixed full-tree maximum root-tip anchor and retain residual tip ages' if dated else 'No rescaling or dating','expected_taxa':sorted(keep)}
def selftest():
 r=parse_newick('((A:0.1,X:0.2)L:0.3,(B:0.4,C:0.5)R:0.6)ROOT;')
 p=prune(r,{'A','B'});fn,fp,fd,ft=maps(r);pn,pp,pd,pt=maps(p)
 assert pair_distance('A','B',fp,fd,ft)[0]==Decimal('1.4')==pair_distance('A','B',pp,pd,pt)[0]
 q=prune(r,{'B','C'});assert len(q.children)==1 and q.source_id==r.source_id
 assert q.children[0].length==Decimal('0.6')
 assert str(pt['A'].length)=='0.4'
 for bad in ['(A:-1,B:1);','(A:1,A:2);','(A:NaN,B:1);']:
  try:parse_newick(bad)
  except ValueError:pass
  else:raise AssertionError('Accepted invalid fixture')
 return {'status':'PASS','exact_decimal_contraction':True,'single_child_original_root_retained':True,'negative_nonfinite_duplicate_inputs_rejected':True}
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--run-dir',type=Path,required=True);args=ap.parse_args();run=args.run_dir.resolve()
 out=run/'03_MATRICES_AND_TREES/cophylogeny_inputs_main25_20260908';cp=run/'checkpoints/cophylogeny_inputs_main25_20260908.complete.json'
 if out.exists() or cp.exists():raise FileExistsError('Refusing existing stage output')
 authorpath=run/'01_PROVENANCE/host_source_author_confirmation_20260908.json';review=run/'02_QC/target_figure_recovery_20260908/main25_identity_review.tsv';legacycross=run/'02_QC/host_sample_crosswalk.tsv';planpath=run/'01_PROVENANCE/frozen_tree_inputs.json'
 sources=[A_TREE,authorpath,review,legacycross,planpath,OLD_AHE,FORMAL_AHE,LEG/'fna_60_ABS60_75_123_ASTRAL_IV.fixed.treefile',LEG/'cophylo_main/DATA/host.pruned_for_cophylogeny.treefile',LEG/'host_wolbachia_43_associations.tsv',LEG/'cophylo_main/DATA/host_wol_association.used.tsv',LEG/'cophylo_main/run.log',Path(__file__).resolve()]
 hashes={str(p):sha(p) for p in sources}
 if hashes[str(authorpath)]!=AUTHOR_SHA or hashes[str(A_TREE)]!=A_SHA:raise ValueError('Frozen author/tree SHA changed')
 author=json.loads(authorpath.read_text());rev=table(review);plan=json.loads(planpath.read_text());main=next(s for s in plan['scenarios'] if s['label']=='main')
 if sha(review)!=author['main25_source_review_sha256'] or len(rev)!=25:raise ValueError('Review not authorized exact25')
 expected=set(main['expected_taxa']);bywol={r['wolbachia_sequence_id']:r for r in rev};confirmed={r['wolbachia_sequence_id']:r for r in author['approved_context_pairs']}
 if len(bywol)!=25 or set(bywol)!=expected or set(confirmed)!=expected:raise ValueError('Main25 mismatch')
 if author['independent_machine_specimen_verification'] is not False or author['evidence_type']!='author_provided_source_confirmation':raise ValueError('Author evidence type changed')
 old,nold=fasta_selected(OLD_AHE,{r['AHE_original_row_ID'] for r in rev});formal,nformal=fasta_selected(FORMAL_AHE,{r['A495_exact_ID_candidate'] for r in rev})
 if nformal!=495:raise ValueError('FormalAHE count')
 cross=[]
 for w,r in sorted(bywol.items()):
  host=r['A495_exact_ID_candidate'];oid=r['AHE_original_row_ID'];a=confirmed[w]
  if any(a[k]!=r[k] for k in ['AHE_original_row_ID','A495_exact_ID_candidate','AHE_row_sha256']):raise ValueError('Author reviewed pair differs')
  if old[oid]['sequence']!=formal[host]['sequence'] or old[oid]['sha256']!=r['AHE_row_sha256'] or old[oid]['length']!=int(r['AHE_row_length']):raise ValueError('AHE continuity failed')
  summary=Path(r['Wol_summary_evidence_path']);hashes[str(summary)]=sha(summary)
  sample=re.search(r'^Sample:\s*(.+)$',summary.read_text(),re.M)
  if not sample or sample.group(1).strip()!=r['Wol_original_sample_ID']:raise ValueError('Wol sample ID differs')
  cross.append({'host_tip':host,'wol_tip':w,'AHE_original_row_ID':oid,'AHE_row_length':old[oid]['length'],'AHE_original_row_sha256':old[oid]['sha256'],'AHE_formal_row_sha256':formal[host]['sha256'],'machine_AHE_sequence_continuity':'EXACT_BODY_IDENTITY','machine_specimen_or_raw_read_identity_independently_verified':False,'author_same_sequencing_files_confirmation':True,'author_evidence_type':author['evidence_type'],'author_confirmation_path':str(authorpath),'author_confirmation_sha256':AUTHOR_SHA,'author_recorded_at':author['recorded_at'],'author_scope':'current_main25_only','source_link_status':'SOURCE_CONFIRMED_BY_AUTHOR_WITH_EXACT_AHE_ROW_CONTINUITY','Wol_original_sample_ID':r['Wol_original_sample_ID'],'Wol_contig_source':r['Wol_contig_source'],'Wol_summary_evidence_path':str(summary),'Wol_summary_sha256':hashes[str(summary)],'special_label_case':r['specific_conflict_or_rename_case'],'final_sample_supergroup':'NOT_ASSIGNED'})
 if len({r['host_tip'] for r in cross})!=25:raise ValueError('Duplicatehost')
 root=parse_newick(A_TREE.read_text())
 if len(descend(root))!=495:raise ValueError('A495 count')
 fixture=selftest();out.mkdir()
 write_tsv(out/'main25_source_confirmed_crosswalk.tsv',cross)
 write_tsv(out/'AHE_main25_associations.tsv',[{'host_tip':r['host_tip'],'wol_tip':r['wol_tip']} for r in cross])
 ahe=tree_artifacts(out,'AHE_A495_main25',root,{r['host_tip'] for r in cross},dated=True)
 ahe.update({'scenario_label':'AHE_A495_main25_author_confirmed','source_tree':str(A_TREE),'source_tree_sha256':A_SHA,'association_table':str(out/'AHE_main25_associations.tsv'),'association_sha256':sha(out/'AHE_main25_associations.tsv'),'source_confirmed_crosswalk':str(out/'main25_source_confirmed_crosswalk.tsv'),'source_confirmed_crosswalk_sha256':sha(out/'main25_source_confirmed_crosswalk.tsv'),'source_confirmation_status':'author_provided_plus_machine_AHE_row_continuity','independent_machine_specimen_verification':False})
 # Old USCO source is explicitly identified by its actual run log and used association table.
 usco={};ulog=(LEG/'cophylo_main/run.log').read_text();upath=LEG/'fna_60_ABS60_75_123_ASTRAL_IV.fixed.treefile';usedpath=LEG/'cophylo_main/DATA/host.pruned_for_cophylogeny.treefile'
 if '[OK] host_tree = '+str(upath) not in ulog:raise ValueError('USCO source not in old run log')
 oldassoc=table(LEG/'host_wolbachia_43_associations.tsv');usedassoc=table(LEG/'cophylo_main/DATA/host_wol_association.used.tsv')
 direct={(r['host_taxon'],r['wolbachia_taxon']) for r in oldassoc};used={(r['host_tip'],r['wol_tip']) for r in usedassoc}
 if direct!=used:raise ValueError('Original/used historicalassociation differs')
 uscomap=[{'host_tip':h,'wol_tip':w,'association_evidence':'exact_pair_in_historical_input_and_actual_used_table','current_author_confirmation_scope':'AHE25 only; this USCO association remains historical sensitivity','source_association_table':str(LEG/'cophylo_main/DATA/host_wol_association.used.tsv')} for h,w in sorted(used) if w in expected]
 if len(uscomap)!=25 or {r['wol_tip'] for r in uscomap}!=expected or len({r['host_tip'] for r in uscomap})!=25:raise ValueError('USCO exact25 not complete')
 uraw=parse_newick(upath.read_text());uused=parse_newick(usedpath.read_text());un,up,ud,ut=maps(uraw);vn,vp,vd,vt=maps(uused);ukeep={r['host_tip'] for r in uscomap}
 uscopairs=[]
 for a,b in itertools.combinations(sorted(ukeep),2):
  d,_=pair_distance(a,b,up,ud,ut);d2,_=pair_distance(a,b,vp,vd,vt)
  if abs(d-d2)>Decimal('0.000000001'):raise ValueError('RawUSCO vs usedUSCO distance differs')
  uscopairs.append({'host_tip_a':a,'host_tip_b':b,'USCO_raw107_distance':decstr(d),'USCO_actual_used43_distance':decstr(d2),'absolute_difference':decstr(abs(d-d2)),'within_1e_minus9':1})
 write_tsv(out/'USCO_main25_historical_associations.tsv',uscomap)
 write_tsv(out/'USCO_raw107_vs_actual_used43_all300_distances.tsv',uscopairs)
 usco=tree_artifacts(out,'USCO_actual_used43_main25_legacy_sensitivity',uused,ukeep,dated=False)
 usco.update({'scenario_label':'USCO_legacy_association_main25_sensitivity','source_tree':str(usedpath),'source_tree_sha256':hashes[str(usedpath)],'raw_backbone_tree':str(upath),'raw_backbone_sha256':hashes[str(upath)],'association_table':str(out/'USCO_main25_historical_associations.tsv'),'association_sha256':sha(out/'USCO_main25_historical_associations.tsv'),'association_status':'historical_exact_pair_sensitivity_not_new_AHE_source_confirmation','raw107_vs_used43_pair_distance_max_abs_diff':decstr(max(Decimal(r['absolute_difference']) for r in uscopairs)),'not_a_dated_tree':True,'no_fuzzy_matching':True})
 for p,h in hashes.items():
  if sha(Path(p))!=h:raise ValueError('Input changed '+p)
 manifest={'schema':1,'status':'complete_source_confirmed_main25_inputs','created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'main25_wol_taxa':sorted(expected),'association_count':25,'AHE_original_alignment_rows':nold,'AHE_formal_alignment_rows':nformal,'AHE_exact_row_continuity_count':25,'author_evidence':{'path':str(authorpath),'sha256':AUTHOR_SHA,'type':author['evidence_type'],'scope':author['scope_applied'],'independent_machine_specimen_verification':False},'source_confirmed_crosswalk':str(out/'main25_source_confirmed_crosswalk.tsv'),'AHE':ahe,'USCO_legacy_sensitivity':usco,'synthetic_algorithm_checks':fixture,'old43_crosswalk_unchanged_sha256':hashes[str(legacycross)],'input_sha256':hashes,'Wolbachia_tree_selected_here':False,'statistics_run_here':False}
 write_json(out/'input_schema.json',{'schema':1,'association_columns':{'host_tip':'exact tree label; never strip suffix','wol_tip':'frozen original Wolbachia sample label'},'source_confirmed_crosswalk_fields_separate':['author_same_sequencing_files_confirmation','machine_AHE_sequence_continuity','machine_specimen_or_raw_read_identity_independently_verified'],'tree_root_policy':'preserve source root, including unary root if pruning removes its other children','age_policy':'use AHE full-tree age anchor and original node map; do not infer a new origin from the pruned MRCA or force tip ages to zero','distance_validation':'all300 original vs pruned distances equal as Decimal parsed from frozen Newick; no rounding during path contraction','USCO_policy':'historical explicit association sensitivity, not dated and not independently source-confirmed by current AHE author statement','manifest':'input_manifest.json','completion_checkpoint':str(cp)})
 write_json(out/'input_manifest.json',manifest)
 with (out/'README_zh.md').open('x') as f:f.write('\n'.join(['# 主25共系统发育输入（2026-09-08）','','25对AHE宿主—Wolbachia配对已依据作者明确同源测序文件声明与25/25旧→正式AHE行序列完全一致性建立。作者声明不冒称机器独立标本鉴定；两个证据字段分开保留。旧43行UNRESOLVED快照不改。','',f"AHE宿主树来自9月4日正式A495定年树。保留原根（剪枝后子节点数{ahe['pruned_root_child_count']}）和全部路径枝长；300对patristic距离以及共同祖先节点ID逐对精确相等。只将被删除侧枝造成的单子内部路径长度相加，不重定根或重定年。",'',f"原树根到tip范围为{ahe['full_root_to_tip_min']}–{ahe['full_root_to_tip_max']} Ma，差异{ahe['full_root_to_tip_range']}来自保存小数的非严格超度量性。节点年龄采用原完整树最大根到tip距离{ahe['node_age_anchor_full_tree_max_root_tip']}为固定基准，保留原始残余tip年龄；新子集不能自行把共同祖先当成原根或用子集最大深度替换基准。",'','USCO敏感性树来自旧真实PACo运行保存的host.pruned_for_cophylogeny.treefile，并已与其日志指定的107-tip USCO原骨架核对全部300对距离。25对关联直接取历史输入表与实际使用表中一致的精确pair，不按同种名/后缀模糊转换。这一输入保留为历史对应敏感性，不是定年树，也不把当前AHE来源声明外推成USCO独立来源确认。','','本阶段未选择Wolbachia树、未运行PACo/ParaFit/eMPress、未分配最终超群。下游必须按host_tip/wol_tip关联，读取input_manifest.json和完成checkpoint的哈希；AHE分支尺度与USCO原报告枝长不能当作同一时间尺度。','']))
 outputs={str(p):sha(p) for p in sorted(out.iterdir()) if p.is_file()}
 completion={'schema':1,'status':'complete','stage':'cophylogeny_inputs_main25_20260908','completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'manifest':str(out/'input_manifest.json'),'manifest_sha256':sha(out/'input_manifest.json'),'output_sha256':outputs,'input_sha256':hashes,'AHE_gate':'source_confirmed_by_author_and_exact_AHE_sequence_continuity','machine_independent_specimen_verification':False,'paired_associations':25,'AHE_all300_exact_distances':True,'USCO_status':'legacy_exact_association_sensitivity','source_confirmed_crosswalk':str(out/'main25_source_confirmed_crosswalk.tsv'),'AHE_tree':ahe['tree'],'AHE_associations':ahe['association_table'],'USCO_tree':usco['tree'],'USCO_associations':usco['association_table']}
 write_json(cp,completion)
 print(json.dumps({'checkpoint':str(cp),'checkpoint_sha256':sha(cp),'AHE':ahe,'USCO':usco,'crosswalk_sha256':sha(out/'main25_source_confirmed_crosswalk.tsv'),'outputs':len(outputs),'fixture':fixture},ensure_ascii=False,indent=2))
if __name__=='__main__':main()
