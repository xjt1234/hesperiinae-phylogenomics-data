#!/usr/bin/env python3
"""Freeze 100 preselected raw UFBoot lines and AHE-only manifests; never fit statistics."""
import argparse,collections,copy,csv,datetime,hashlib,importlib.util,io,itertools,json,random,sys
from pathlib import Path
sys.dont_write_bytecode=True
SELECTION_SHA='d04ab1a2fb60f6f34950de672f67a15ba304bef2370d4dea68b6eead0fabb943'
FORMAL_PLAN_SHA='96fc5349d49f7cb73d758a6c094d8c815f34181098e40ac60a0e369cc5ca6094'
HELPER_SHA='6d5b02a881541c4d68298dc5ae0325e59afd1e9efa7d8b21421a8175c5a2bb46'
R_SHA='df207f4e1cf052f2a2e86490f272aa73689f2b68eb92e80f103b5c7578a7c982'
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
def bsha(b):return hashlib.sha256(b).hexdigest()
def tsv(p):
 with Path(p).open() as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rows):
 with p.open('x',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def jwrite(p,o):
 with p.open('x') as f:json.dump(o,f,ensure_ascii=False,indent=2);f.write('\n')
def cp_records(o):
 d={r['path']:r['sha256'] for r in o.get('outputs',[])}
 for p,h in o.get('output_sha256',{}).items():
  if p in d and d[p]!=h:raise ValueError('Conflicting checkpoint output hashes')
  d[p]=h
 return d
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--run-dir',type=Path,required=True);args=ap.parse_args();run=args.run_dir.resolve()
 out=run/'03_MATRICES_AND_TREES/cophylogeny_bootstrap100_inputs_20260908';inputcp=run/'checkpoints/cophylogeny_bootstrap100_inputs.complete.json';completecp=run/'checkpoints/cophylogeny_bootstrap100_preparation.complete.json';planpath=run/'01_PROVENANCE/cophylogeny_bootstrap100_100scenarios_20260908.json'
 if any(p.exists() for p in [out,inputcp,completecp,planpath]):raise FileExistsError('Refusing existing preparation')
 selectionpath=run/'01_PROVENANCE/cophylogeny_bootstrap100_selection_20260908.json';formalpath=run/'01_PROVENANCE/cophylogeny_main25_16scenarios_20260908.json';helperpath=run/'scripts/prepare_cophylogeny_main25_inputs_20260908.py';rpath=run/'scripts/frozen/cophylogeny_exact_inputs_R26_df207f4e1cf052f2.R'
 if sha(selectionpath)!=SELECTION_SHA or sha(formalpath)!=FORMAL_PLAN_SHA or sha(helperpath)!=HELPER_SHA or sha(rpath)!=R_SHA:raise ValueError('Frozen selection/method/helper changed')
 selection=json.loads(selectionpath.read_text());formal=json.loads(formalpath.read_text());sel=selection['one_based_tree_indices'];expected=sorted(random.Random(selection['selection_seed']).sample(range(1,1001),100))
 if selection['status']!='frozen_before_bootstrap_cophylogeny' or sel!=expected or len(sel)!=100 or selection['nperm_each']!=9999 or selection['host_scenario']!='AHE_A495_main25_author_confirmed':raise ValueError('Selection/analysis settings differ')
 primary=next(s for s in formal['scenarios'] if s['id']==formal['primary_scenario']);assert primary['id']=='AHE_A495_dated__main_mfp' and formal['seed']==20260908 and formal['nperm']==9999 and formal['correction']=='cailliez' and formal['symmetric'] is False
 hostrec=primary['inputs']['host_tree'];linksrec=primary['inputs']['links'];hostcp_rec=next(r for r in formal['analysis_checkpoints'] if 'covered_evidence' in r['path']);hostcp=json.loads(Path(hostcp_rec['path']).read_text())
 if sha(hostcp_rec['path'])!=hostcp_rec['sha256'] or hostcp['status']!='complete':raise ValueError('Author-covered host checkpoint differs')
 hostcoverage=cp_records(hostcp)
 for p,h in hostcoverage.items():
  if sha(p)!=h:raise ValueError('Host checkpoint output differs')
 for r in [hostrec,linksrec,formal['association_evidence']]:
  if hostcoverage.get(r['path'])!=r['sha256'] or sha(r['path'])!=r['sha256']:raise ValueError('Author/host/links not covered')
 linkrows=tsv(linksrec['path']);woltips={r['wol_tip'] for r in linkrows};hosttips={r['host_tip'] for r in linkrows}
 if len(linkrows)!=len(woltips)!=25:raise ValueError('Association count')
 if len(linkrows)!=25 or len(woltips)!=25 or len(hosttips)!=25:raise ValueError('Nonbijective25')
 treeplanpath=run/'01_PROVENANCE/frozen_tree_inputs.json';treeplan=json.loads(treeplanpath.read_text());ms=next(s for s in treeplan['scenarios'] if s['label']=='main');matrix=Path(ms['matrix'])
 if set(ms['expected_taxa'])!=woltips or sha(matrix)!=ms['matrix_sha256']:raise ValueError('Source matrix main25 differs')
 with matrix.open() as f:matrixids=[line[1:].strip().split()[0] for line in f if line.startswith('>')]
 if len(matrixids)!=25 or set(matrixids)!=woltips:raise ValueError('Matrix exact IDs differ')
 source=Path(selection['source_ufboot']);sourcecp=run/'checkpoints/homology_coding_tree_main.complete.json';sco=json.loads(sourcecp.read_text());scoverage=cp_records(sco)
 if sco['status']!='complete' or scoverage.get(str(source))!=selection['source_ufboot_sha256'] or sha(source)!=selection['source_ufboot_sha256']:raise ValueError('Source UFBoot checkpoint binding failed')
 for p,h in scoverage.items():
  if sha(p)!=h:raise ValueError('Main-tree completed output changed')
 blob=source.read_bytes();lines=blob.splitlines(keepends=True)
 if len(lines)!=1000 or any(not x.strip() for x in lines):raise ValueError('One tree per source line required')
 spec=importlib.util.spec_from_file_location('exact_host_input_helpers',helperpath);h=importlib.util.module_from_spec(spec);spec.loader.exec_module(h)
 from Bio import Phylo
 inputpaths=[selectionpath,formalpath,helperpath,rpath,source,sourcecp,Path(hostcp_rec['path']),Path(hostrec['path']),Path(linksrec['path']),Path(formal['association_evidence']['path']),treeplanpath,matrix,Path(__file__).resolve()];inputhash={str(p):sha(p) for p in inputpaths}
 out.mkdir();(out/'trees').mkdir();inventory=[];distrows=[];rootrows=[];allbranchrows=[];offset=0;selectedbytes=[];negative=0;maxfloaterror=0
 for i,raw in enumerate(lines,1):
  start=offset;offset+=len(raw)
  if i not in sel:continue
  root=h.parse_newick(raw.decode('ascii').strip());nodes,parent,depth,tips=h.maps(root)
  if set(tips)!=woltips or len(tips)!=25:raise ValueError('Tree exact25 mismatch at '+str(i))
  branch=[n for n in nodes.values() if n.source_id!=root.source_id]
  if any(n.length is None or not n.length.is_finite() or n.length<0 for n in branch):raise ValueError('Missing or invalid branch')
  # Root is IQ-TREE's storage node. Its child groups are audited, never biologically rerooted.
  bio=Phylo.read(io.StringIO(raw.decode('ascii')),'newick')
  if len(bio.get_terminals())!=25 or set(t.name for t in bio.get_terminals())!=woltips:raise ValueError('Independent Bio parser mismatch')
  rootparts=sorted(';'.join(sorted(h.descend(c))) for c in root.children)
  treepath=out/'trees'/f'main_mfp_ufboot_{i:04d}.nwk'
  with treepath.open('xb') as f:f.write(raw)
  if treepath.read_bytes()!=raw:raise ValueError('Tree line copy not exact bytes')
  selectedbytes.append(raw);pairvalues=[];floating=[]
  for a,b in itertools.combinations(sorted(woltips),2):
   distance,mrca=h.pair_distance(a,b,parent,depth,tips)
   if not distance.is_finite() or distance<0:raise ValueError('Invalid pairwise distance')
   error=abs(float(distance)-bio.distance(a,b));floating.append(error);pairvalues.append(distance)
   if error>1e-12:raise ValueError('Independent cophenetic mismatch')
   distrows.append({'source_tree_1based_index':i,'taxon_a':a,'taxon_b':b,'patristic_distance_decimal':h.decstr(distance),'BioPhylo_distance':format(bio.distance(a,b),'.17g'),'independent_abs_float_error':format(error,'.17g'),'MRCA_storage_node_id':mrca,'source_line_bytes_sha256':bsha(raw),'no_rerooting_or_branch_changes':1})
  if len(pairvalues)!=300 or not any(x>0 for x in pairvalues):raise ValueError('Degenerate distances')
  maxfloaterror=max(maxfloaterror,max(floating));newhash=sha(treepath)
  inventory.append({'execution_order_1based':len(inventory)+1,'source_tree_1based_index':i,'source_ufboot_path':str(source),'source_ufboot_sha256':selection['source_ufboot_sha256'],'source_line_byte_offset_0based':start,'source_line_bytes':len(raw),'source_line_bytes_sha256':bsha(raw),'source_line_newick_without_lineending_sha256':bsha(raw.rstrip(b'\r\n')),'new_tree_path':str(treepath),'new_tree_sha256':newhash,'exact_source_line_bytes_preserved':int(newhash==bsha(raw)),'tip_count':25,'branch_count':len(branch),'minimum_branch_length':h.decstr(min(n.length for n in branch)),'maximum_branch_length':h.decstr(max(n.length for n in branch)),'zero_branch_count':sum(n.length==0 for n in branch),'tree_total_branch_length':h.decstr(sum((n.length for n in branch),h.ZERO)),'root_child_count':len(root.children),'root_label':root.label,'root_branch_length':h.decstr(root.length) if root.length is not None else '', 'root_interpretation':'unrooted_IQTREE_storage_node_not_biological_root','all300_distances_nonnegative_finite':1,'zero_pairwise_distances':sum(x==0 for x in pairvalues),'minimum_pairwise_distance':h.decstr(min(pairvalues)),'maximum_pairwise_distance':h.decstr(max(pairvalues)),'scenario_index':1000+i,'R_seed':formal['seed']+1000+i})
  rootrows.append({'source_tree_1based_index':i,'root_label':root.label,'root_child_count':len(root.children),'root_child_descendant_groups_json':json.dumps(rootparts,separators=(',',':')),'internal_outdegree_counts_json':json.dumps(dict(collections.Counter(len(n.children) for n in nodes.values() if n.children)),sort_keys=True),'unary_internal_nodes':sum(len(n.children)==1 for n in nodes.values() if n.children),'all_exact_tips':';'.join(sorted(tips)),'no_biological_root_or_dates_inferred':1})
  for sid,n in nodes.items():allbranchrows.append({'source_tree_1based_index':i,'storage_node_id':sid,'storage_parent_id':parent[sid] or '', 'label':n.label,'is_tip':int(not n.children),'branch_length_original_decimal':h.decstr(n.length) if n.length is not None else '', 'storage_root_distance_decimal':h.decstr(depth[sid]),'descendant_tips':';'.join(sorted(h.descend(n)))})
 if len(inventory)!=100 or len(distrows)!=30000:raise ValueError('Count mismatch')
 with (out/'selected100_exact_source_lines.ufboot').open('xb') as f:f.write(b''.join(selectedbytes))
 write(out/'tree_inventory.tsv',inventory);write(out/'all100_trees_300pairs.tsv',distrows);write(out/'tree_root_and_tip_review.tsv',rootrows);write(out/'tree_node_branch_inventory.tsv',allbranchrows)
 freeze={'schema':1,'status':'frozen_selected_tree_inputs','source_selection':{'path':str(selectionpath),'sha256':SELECTION_SHA},'source_ufboot':{'path':str(source),'sha256':selection['source_ufboot_sha256']},'source_tree_checkpoint':{'path':str(sourcecp),'sha256':inputhash[str(sourcecp)]},'host_checkpoint':hostcp_rec,'method_source_plan':{'path':str(formalpath),'sha256':FORMAL_PLAN_SHA},'R_entrypoint':{'path':str(rpath),'sha256':R_SHA},'tree_indices_1based':sel,'selection_random_draw_reproduced':True,'sequence_taxa':sorted(woltips),'host_taxa':sorted(hosttips),'original_source_lines_copied_verbatim':True,'statistics_run':False,'new_tree_inference':False,'all_input_sha256':inputhash}
 jwrite(out/'frozen_input_sources.json',freeze)
 qa={'status':'PASS_machine_bootstrap100_tree_inputs','trees':100,'source_trees':1000,'all25_exact_tips_each':True,'all_source_line_byte_hashes_equal_output_tree_hashes':True,'all_nonroot_edges_explicit_nonnegative_finite':True,'patristic_pairs_checked':30000,'all_distance_matrices_nondegenerate':True,'BioPhylo_independent_max_abs_float_error':maxfloaterror,'root_child_count_distribution':dict(collections.Counter(r['root_child_count'] for r in inventory)),'biological_root_inferred':False,'statistics_run':False,'selection_reproduced_from_predeclared_seed':True}
 jwrite(out/'tree_input_QA.json',qa)
 for p,s in inputhash.items():
  if sha(p)!=s:raise ValueError('Source changed while preparing')
 records={str(p):sha(p) for p in sorted(out.rglob('*')) if p.is_file()}
 # Tree input checkpoint precedes plans, avoiding cyclic plan/checkpoint hashing.
 jwrite(inputcp,{'schema':1,'status':'complete','stage':'cophylogeny_bootstrap100_tree_inputs','completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'summary':qa,'source_checkpoints':[{'path':str(sourcecp),'sha256':inputhash[str(sourcecp)]},hostcp_rec],'selection_path':str(selectionpath),'selection_sha256':SELECTION_SHA,'input_sha256':inputhash,'outputs':[{'path':p,'sha256':s} for p,s in records.items()]})
 scenarios=[]
 for inv in inventory:
  i=inv['source_tree_1based_index'];sid=f'AHE_A495_dated__main_mfp_ufboot_{i:04d}'
  scenarios.append({'id':sid,'scenario_index':1000+i,'execution_order_1based':inv['execution_order_1based'],'role':'bootstrap_sensitivity','n_links':25,'host_tree_basis':primary['host_tree_basis'],'wol_tree_basis':f'Preselected source main MFP UFBoot tree {i}/1000, original line bytes and all branch lengths retained; unrooted storage representation','association_basis':primary['association_basis'],'host_distance_unit':primary['host_distance_unit'],'wol_distance_unit':primary['wol_distance_unit'],'bootstrap_source_tree_1based_index':i,'bootstrap_source_ufboot':str(source),'bootstrap_source_ufboot_sha256':selection['source_ufboot_sha256'],'bootstrap_source_line_bytes_sha256':inv['source_line_bytes_sha256'],'expected_R_seed':inv['R_seed'],'inputs':{'host_tree':copy.deepcopy(hostrec),'wol_tree':{'path':inv['new_tree_path'],'sha256':inv['new_tree_sha256']},'links':copy.deepcopy(linksrec)}})
 analysis={'schema':1,'status':'complete','data_kind':'real','created_before_fitting':datetime.datetime.now(datetime.timezone.utc).isoformat(),'analysis_label':'R2.6 AHE main25 sensitivity to 100 preselected original main UFBoot trees','nperm':9999,'seed':formal['seed'],'correction':formal['correction'],'symmetric':formal['symmetric'],'association_status':formal['association_status'],'association_statement':'User explicitly confirmed that frozen main25 Wolbachia and AHE sequences came from the same sequencing files. Exact old-to-formal AHE row continuity is separately machine-verified. This is author-provided source metadata, not independent machine specimen validation. This plan contains AHE host inputs only.','association_evidence':copy.deepcopy(formal['association_evidence']),'analysis_checkpoints':[copy.deepcopy(hostcp_rec),{'path':str(inputcp),'sha256':sha(inputcp)}],'source_selection':{'path':str(selectionpath),'sha256':SELECTION_SHA},'reference_primary_scenario':formal['primary_scenario'],'reference_primary_plan':{'path':str(formalpath),'sha256':FORMAL_PLAN_SHA},'first_execution_scenario':scenarios[0]['id'],'scenario_order_fixed_before_statistics':True,'all_cases_are_bootstrap_sensitivities':True,'R_seed_rule':'base seed 20260908 + scenario_index (1000 + 1-based source UFBoot index)','frozen_R_entrypoint':{'path':str(rpath),'sha256':R_SHA},'required_R_packages':{'paco':'0.4.2','ape':'5.8.1','vegan':'2.7.2','jsonlite':'2.0.0','digest':'0.6.39'},'interpretation_rules':formal['interpretation_rules'][:3]+['These are the 100 source UFBoot line indices selected before cophylogeny statistics; no selection by fit or p value.','All host trees and exact25 associations are fixed to the AHE primary input.','Sensitivity includes each UFBoot topology and its reported branch lengths; tree storage roots have no ancestral meaning.','UFBoot samples are not a Bayesian posterior, and resulting p-value variability is not transfer or divergence-time probability.','Original main ML fit remains the primary; no UFBoot replicate is promoted to a new primary based on results.','No per-link hypothesis tests or automatic supergroup assignments.'],'scenarios':scenarios}
 jwrite(planpath,analysis);cases=out/'case_manifests';cases.mkdir()
 for s in scenarios:jwrite(cases/(s['id']+'.json'),dict(analysis,scenarios=[s]))
 planqa={'status':'PASS_machine_manifest_preparation_no_R_invoked','scenario_count':100,'unique_scenario_ids':len({s['id'] for s in scenarios}),'exact_selection_order_preserved':[s['bootstrap_source_tree_1based_index'] for s in scenarios]==sel,'scenario_indices_unique_nonnegative':len({s['scenario_index'] for s in scenarios})==100 and all(s['scenario_index']>=0 for s in scenarios),'all_R_seeds_match_frozen_selection_rule':all(analysis['seed']+s['scenario_index']==20260908+1000+s['bootstrap_source_tree_1based_index']==s['expected_R_seed'] for s in scenarios),'nperm_each':9999,'method_source_plan_sha256':FORMAL_PLAN_SHA,'R_entrypoint_unchanged_sha256':R_SHA,'host_author_evidence_covered':True,'all_tree_inputs_covered_by_complete_checkpoints':True,'R_check_only_executed':False,'R_statistics_executed':False,'R_version_gate_not_bypassed':True,'package_compatibility_note':'R df207 gate must use previously audited versions, including vegan2.7.2; no R execution attempted while root restores the audited library','plan':str(planpath),'plan_sha256':sha(planpath)}
 # Prove every case's three inputs and author evidence is covered by declared complete records.
 covered=dict(hostcoverage);covered.update(records)
 for s in scenarios:
  for r in s['inputs'].values():
   if covered.get(r['path'])!=r['sha256']:raise ValueError('Case input not covered')
 if covered.get(formal['association_evidence']['path'])!=formal['association_evidence']['sha256']:raise ValueError('Author evidence not covered')
 jwrite(out/'analysis_plan_QA.json',planqa)
 with (out/'README_zh.md').open('x') as f:f.write('\n'.join(['# 100棵预选UFBoot树的AHE主25敏感性输入','','按已冻结selection的100个1-based行号提取原main MFP的UFBoot树；同一seed重新抽样验证索引完全一致。每个新树文件保留源行全部字节（包含原换行），表中记录完整源文件SHA、源行字节偏移、行字节数、带换行及去换行SHA。未重新推断树、未改任何枝长或tip。','','100棵树逐棵均为同样25个Wol样本，全部非根边有有限非负长度；30,000对patristic距离使用Decimal计算并由独立Bio.Phylo读取核对。IQ-TREE的存储根不是生物学根，也不据此解释传播方向或时间。AHE树保留原单子根与完整树年龄基准，宿主exactID—菌ID仍按作者来源确认和机器行连续性分别记录。','','100个情景只使用固定AHE宿主；R驱动与正式16方案同为df207版本、9999置换、Cailliez、非对称PACo、共享双边际一一对应置换主零模型和PACo r0敏感性。保持base_seed=20260908，scenario_index=1000+源树行号，实现预声明的R种子规则。没有按拟合值或p值选树；原主ML情景仍是正式主分析。','','计划文件：'+str(planpath),'树输入checkpoint：'+str(inputcp),'完成准备checkpoint：'+str(completecp),'','树输入checkpoint先生成，计划再引用其SHA，避免循环哈希；100份single-case manifest已新建，供根线程wrapper使用。所有情景role=bootstrap_sensitivity；没有将第一棵或“最好”的bootstrap树当作primary。','','本阶段没有调用R（包括check-only），未规避R包版本门槛；正式执行应使用已核验库，vegan版本仍要求2.7.2。没有启动置换统计。','']))
 finalpaths=[p for p in sorted(out.rglob('*')) if p.is_file()]+[inputcp,planpath]
 finalhash={str(p):sha(p) for p in finalpaths}
 jwrite(completecp,{'schema':1,'status':'complete','stage':'cophylogeny_bootstrap100_preparation','completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'tree_count':100,'scenario_count':100,'statistics_executed':False,'tree_input_checkpoint':{'path':str(inputcp),'sha256':sha(inputcp)},'plan':{'path':str(planpath),'sha256':sha(planpath)},'selection':{'path':str(selectionpath),'sha256':SELECTION_SHA},'input_sha256':inputhash,'outputs':[{'path':p,'sha256':s} for p,s in finalhash.items()]})
 print(json.dumps({'tree_input_checkpoint':str(inputcp),'tree_input_checkpoint_sha256':sha(inputcp),'plan':str(planpath),'plan_sha256':sha(planpath),'preparation_checkpoint':str(completecp),'preparation_checkpoint_sha256':sha(completecp),'tree_QA':qa,'manifest_QA':planqa,'output_files':len(finalhash)},ensure_ascii=False,indent=2))
if __name__=='__main__':main()
