#!/usr/bin/env python3
"""Compare two frozen full495 calibrations for the same 24 local branch pairs."""
import argparse, collections, csv, datetime, hashlib, itertools, json, math, platform, sys, time
from decimal import Decimal, getcontext
from pathlib import Path
import Bio
from Bio import Phylo
getcontext().prec = 50
TOL = Decimal('1e-12')
A_SHA = '28d4460885cc7425c9c703b79aff43c24a2ab931c86110c5ffee8d6533309eba'
B_SHA = '0e1e17d539019d7c5c869941992838d65f371d77e219d90fb9033f9507d839a9'
ALGORITHM_SHA = '2ebe337e2e39291d315f04050dfdf800ef943850929451237cf9c58a76dd21b4'
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def readrows(p):
    with Path(p).open() as f: return list(csv.DictReader(f, delimiter='\t'))
def save(path, value):
    with path.open('x', encoding='utf-8') as f: json.dump(value, f, ensure_ascii=False, indent=2); f.write('\n')
def tsv(path, data):
    assert data
    with path.open('x', encoding='utf-8', newline='') as f:
        w=csv.DictWriter(f, list(data[0]), delimiter='\t'); w.writeheader(); w.writerows(data)
def textdec(x): return '' if x is None else str(x)
def signature(tips): return sha_bytes(json.dumps(sorted(tips), separators=(',',':')).encode())
def sha_bytes(b): return hashlib.sha256(b).hexdigest()
def tree_data(path):
    tr=Phylo.read(path, 'newick'); nodes=list(tr.find_clades(order='preorder'))
    tips=[n.name for n in tr.get_terminals()]
    assert len(tips)==len(set(tips))==495 and all(tips), 'Expected 495 exact unique host IDs'
    assert len(nodes)==989 and all(len(n.clades) in (0,2) for n in nodes), 'Expected rooted full495 binary tree'
    parent={c:n for n in nodes for c in n.clades}; depth={tr.root:Decimal(0)}
    for n in nodes:
        for c in n.clades:
            assert c.branch_length is not None and math.isfinite(c.branch_length) and c.branch_length>=0
            depth[c]=depth[n]+Decimal(str(c.branch_length))
    desc={}
    for n in reversed(nodes):
        desc[n]=frozenset([n.name]) if not n.clades else frozenset().union(*(desc[c] for c in n.clades))
    lookup={v:k for k,v in desc.items()}; assert len(lookup)==989, 'Ambiguous full descendant clade'
    anchor=max(depth[n] for n in tr.get_terminals()); age={n:anchor-depth[n] for n in nodes}
    assert all(v>=0 and v.is_finite() for v in age.values())
    return dict(tree=tr,nodes=nodes,tips=set(tips),parent=parent,depth=depth,desc=desc,lookup=lookup,anchor=anchor,age=age)
def intersect(ca,pa,cb,pb):
    if any(x is None for x in (ca,pa,cb,pb)):return (None,None,'not_available')
    assert 0<=ca<=pa and 0<=cb<=pb
    low,high=max(ca,cb),min(pa,pb)
    return low,high,('feasible' if low<high else 'boundary_only' if low==high else 'infeasible')
def main():
    started=time.monotonic()
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--run-dir',type=Path,required=True)
    p.add_argument('--out-dir',type=Path); args=p.parse_args(); r=args.run_dir.resolve()
    out=args.out_dir or r/'05_RECONCILIATION_AND_TIMING/calibration_AB_main25_20260908'
    assert not out.exists(), 'Refuse overwrite existing output directory'
    paths={
      'tree_A':Path('/home/data/t200301/xjt/Hesperiinae_review/treePL_R2_3_T25_deep_only_20260904_000707/output/T25_deep_only_dated.tre'),
      'tree_B':Path('/home/data/t200301/xjt/Hesperiinae_review/treePL_R2_3_T25_three_calibrations_20260904_075723/output/T25_three_calibrations_dated.tre'),
      'algorithm':r/'01_PROVENANCE/calibration_AB_main25_algorithm_before_run_20260908.md',
      'candidate':r/'05_RECONCILIATION_AND_TIMING/empress_root_topology_20260908/candidate_host_pairs.tsv',
      'branchmap':r/'05_RECONCILIATION_AND_TIMING/empress_root_topology_20260908/host_branch_map.tsv',
      'fullmap':r/'03_MATRICES_AND_TREES/cophylogeny_inputs_main25_20260908/AHE_A495_main25.full_to_pruned_node_map.tsv',
      'associations':r/'03_MATRICES_AND_TREES/cophylogeny_inputs_main25_20260908/AHE_main25_associations.tsv',
      'windows_A':r/'05_RECONCILIATION_AND_TIMING/conditional_windows_main25_20260908/conditional_timewindows.tsv',
      'age_audit_A':r/'05_RECONCILIATION_AND_TIMING/conditional_windows_main25_20260908/conditional_branch_age_audit.tsv',
      'windows_cp':r/'05_RECONCILIATION_AND_TIMING/conditional_windows_main25_20260908/completed.json',
      'export_cp':r/'checkpoints/empress_figure_export_20260908.complete.json',
      'host_cp':r/'checkpoints/cophylogeny_inputs_main25_20260908.covered_evidence.complete.json',
      'script':Path(__file__).resolve()}
    sources={str(q.resolve()):sha(q) for q in paths.values()}
    assert sha(paths['tree_A'])==A_SHA and sha(paths['tree_B'])==B_SHA, 'Frozen tree hash changed'
    assert sha(paths['algorithm'])==ALGORITHM_SHA
    assert sha(paths['candidate'])=='bfa48f9cee7937686ac88a9caeb703325b89bf5d3c935d84d67f8d991cb6ed12'
    assert sha(paths['branchmap'])=='7f3b21b125cf6536749d1529a7255c3b6ca7e723349d22bf1ff64cd627389173'
    assert sha(paths['windows_A'])=='1908c1abf62a5a786f98c8dddd64ccbf368e1df972e5d68b65d4d208e85555fc'
    for role in ['windows_cp','export_cp','host_cp']:
        cp=json.loads(paths[role].read_text()); assert cp['status']=='complete'
        covered=dict(cp.get('output_sha256',{}))
        covered.update({x['path']:x['sha256'] for x in cp.get('outputs',[])})
        for q,h in covered.items(): assert sha(q)==h; sources[str(Path(q).resolve())]=h
        for q,h in cp.get('source_sha256',{}).items(): assert sha(q)==h; sources[str(Path(q).resolve())]=h
    A,B=tree_data(paths['tree_A']),tree_data(paths['tree_B'])
    assert A['tips']==B['tips'], 'STOP: full495 tip identities differ'
    assert set(A['lookup'])==set(B['lookup']), 'STOP: full rooted topology clade sets differ'
    for tips,an in A['lookup'].items():
        bn=B['lookup'][tips]; ap=A['parent'].get(an); bp=B['parent'].get(bn)
        assert (A['desc'][ap] if ap is not None else None)==(B['desc'][bp] if bp is not None else None), 'STOP: full parent identities differ'
    err=[]; exact=0
    def compare(x, raw):
        nonlocal exact
        if x is None: assert raw=='', 'Missing parent must remain empty'; return
        delta=abs(x-Decimal(raw)); err.append(delta); exact+=int(delta==0)
        assert delta<=TOL, 'A baseline age/distance mismatch: '+str(delta)
    fm=readrows(paths['fullmap']); assert len(fm)==989
    aid={}; aname={}
    for row in fm:
        fulltips=frozenset(row['full_descendant_tips'].split(';')); n=A['lookup'][fulltips]
        assert row['full_node_id'] not in aid and n not in aname
        aid[row['full_node_id']]=n; aname[n]=row['full_node_id']
        compare(A['age'][n],row['full_age_from_full_tree_anchor']); compare(A['depth'][n],row['full_root_distance'])
        if n in A['parent']:compare(Decimal(str(n.branch_length)),row['full_incoming_branch_length'])
        else:assert row['full_incoming_branch_length']==''
    assert len(aid)==989
    for row in fm:
        n=aid[row['full_node_id']]; pa=A['parent'].get(n)
        assert row['full_parent_node_id']==(aname[pa] if pa is not None else '')
    bname={n:'B_preorder_'+str(i).zfill(4) for i,n in enumerate(B['nodes'],1)}
    assoc=readrows(paths['associations']); selected={x['host_tip'] for x in assoc}
    assert len(assoc)==len(selected)==len({x['wol_tip'] for x in assoc})==25 and selected<=A['tips']
    display={x['host_tip']:x['wol_tip'].replace('_',' ') for x in assoc}
    br=readrows(paths['branchmap']); assert len(br)==49
    brby={x['host_branch_id']:x for x in br}; sets={k:frozenset(json.loads(v['selected_descendant_tips_json'])) for k,v in brby.items()}
    assert len(brby)==len(set(sets.values()))==49 and all(v and v<=selected for v in sets.values())
    assert sum(len(v)==1 for v in sets.values())==25
    retained={}
    for k,tips in sets.items():
        h=signature(tips); assert brby[k]['clade_hash']==h
        assert k==(next(iter(tips)) if len(tips)==1 else 'H_'+h)
        supers=[j for j,tt in sets.items() if tips<tt]
        retained[k]=min(supers,key=lambda j:len(sets[j])) if supers else None
        assert brby[k]['retained_parent_node_id']==(retained[k] or '')
    assert sum(v is None for v in retained.values())==1
    assert all(sum(v==k for v in retained.values())==(0 if len(sets[k])==1 else 2) for k in sets)
    branches={}; branch_rows=[]
    for k,row in brby.items():
        an=aid[row['full_target_node_id']]; fulltips=A['desc'][an]; bn=B['lookup'][fulltips]
        assert fulltips&selected==sets[k]
        assert A['tree'].common_ancestor(sorted(sets[k])) is an
        assert B['tree'].common_ancestor(sorted(sets[k])) is bn
        pa=A['parent'].get(an); pb=B['parent'].get(bn)
        assert row['full_immediate_parent_node_id']==(aname[pa] if pa is not None else '')
        rp=retained[k]; ra=aid[brby[rp]['full_target_node_id']] if rp else None
        rb=B['lookup'][A['desc'][ra]] if ra is not None else None
        assert row['retained_parent_full_node_id']==(aname[ra] if ra is not None else '')
        ages={}
        for label,t,n,op,pp in [('A',A,an,pa,ra),('B',B,bn,pb,rb)]:
            vals={'child':t['age'][n],'original_parent':t['age'][op] if op is not None else None,'retained_parent':t['age'][pp] if pp is not None else None}
            for parent_age in [vals['original_parent'],vals['retained_parent']]:
                if parent_age is not None:assert parent_age>=vals['child']
            if all(vals[z] is not None for z in ['original_parent','retained_parent']): assert vals['retained_parent']>=vals['original_parent']
            ages[label]=vals
        compare(ages['A']['child'],row['node_age_Ma']);compare(ages['A']['original_parent'],row['original_immediate_parent_age_Ma']);compare(ages['A']['retained_parent'],row['retained_parent_age_Ma'])
        assert (row['is_sampled_root']=='1')==(rp is None)
        rec={'host_branch_id':k,'selected_descendant_tips_json':json.dumps(sorted(sets[k])),'selected_tip_count':len(sets[k]),'full_clade_sha256':signature(fulltips),'full_descendant_tips_json':json.dumps(sorted(fulltips)),'full_tip_count':len(fulltips),'A_full_node_id':aname[an],'B_node_debug_id_not_identity':bname[bn],'is_sampled_root':rp is None,'retained_parent_branch_id':rp or '', 'A_original_parent_id':aname[pa] if pa is not None else '', 'B_original_parent_debug_id':bname[pb] if pb is not None else ''}
        for kind in ['child','retained_parent','original_parent']:
            av,bv=ages['A'][kind],ages['B'][kind]
            rec['A_'+kind+'_age_Ma']=textdec(av);rec['B_'+kind+'_age_Ma']=textdec(bv)
            rec['B_minus_A_'+kind+'_age_Ma']=textdec(bv-av) if av is not None and bv is not None else ''
        branch_rows.append(rec); branches[k]=ages
    candidates=readrows(paths['candidate']); awin={x['candidate_id']:x for x in readrows(paths['windows_A'])}
    assert len(candidates)==len(awin)==24 and {x['candidate_id'] for x in candidates}==set(awin)
    assert len({tuple(sorted([x['host_branch_a'],x['host_branch_b']])) for x in candidates})==24
    ageaudit={(x['candidate_id'],x['side']):x for x in readrows(paths['age_audit_A'])};assert len(ageaudit)==48
    windows=[];counts={s:{ch:collections.Counter() for ch in ['clade','terminal']} for s in ['A','B']}
    for c in candidates:
        cid=c['candidate_id'];ka,kb=c['host_branch_a'],c['host_branch_b']; assert sets[ka].isdisjoint(sets[kb])
        for side,k in [('a',ka),('b',kb)]:
            assert sets[k]==frozenset(json.loads(c['host_'+side+'_tips_json']))==frozenset(json.loads(awin[cid]['lineage_'+side+'_tips_json']))
            old=ageaudit[(cid,side)];assert old['host_branch_id']==k
            for kind in ['child','retained_parent','original_parent']: compare(branches[k]['A'][kind],old[kind+'_age'])
        row=dict(c);row['lineage_a_label']=awin[cid]['lineage_a_label'];row['lineage_b_label']=awin[cid]['lineage_b_label']
        for ch,pk in [('clade','retained_parent'),('terminal','original_parent')]:
            ws={}
            for s in ['A','B']:
                aa,bb=branches[ka][s],branches[kb][s];lo,hi,st=intersect(aa['child'],aa[pk],bb['child'],bb[pk]);ws[s]=(lo,hi,st);counts[s][ch][st]+=1
                row[s+'_'+ch+'_lower_Ma']=textdec(lo);row[s+'_'+ch+'_upper_Ma']=textdec(hi);row[s+'_'+ch+'_status']=st
                if s=='A':
                    compare(lo,awin[cid][ch+'_lower_Ma']);compare(hi,awin[cid][ch+'_upper_Ma']);assert st==awin[cid][ch+'_status']
            for ix,key in [(0,'lower'),(1,'upper')]:row['B_minus_A_'+ch+'_'+key+'_Ma']=textdec(ws['B'][ix]-ws['A'][ix]) if ws['B'][ix] is not None and ws['A'][ix] is not None else ''
            row[ch+'_status_changed']=ws['A'][2]!=ws['B'][2]
        for side,k in [('a',ka),('b',kb)]:
            for s in ['A','B']:row[s+'_lineage_'+side+'_original_stem_parent_Ma']=textdec(branches[k][s]['original_parent'])
            compare(branches[k]['A']['original_parent'],awin[cid]['lineage_'+side+'_stem_age_Ma'])
        windows.append(row)
    assert counts['A']['clade']=={'feasible':16,'infeasible':8} and counts['A']['terminal']=={'feasible':14,'infeasible':10}
    fullrows=[]
    for row in fm:
        an=aid[row['full_node_id']];tips=A['desc'][an];bn=B['lookup'][tips]
        fullrows.append({'A_full_node_id':aname[an],'B_node_debug_id_not_identity':bname[bn],'full_clade_sha256':signature(tips),'full_tip_count':len(tips),'full_descendant_tips_json':json.dumps(sorted(tips)),'rooted_parent_clade_identical':True,'A_age_Ma':str(A['age'][an]),'B_age_Ma':str(B['age'][bn]),'B_minus_A_age_Ma':str(B['age'][bn]-A['age'][an])})
    delta_summary={}
    for ch in ['clade','terminal']:
        delta_summary[ch]={'status_changed_count':sum(x[ch+'_status_changed'] for x in windows),'status_unchanged_count':sum(not x[ch+'_status_changed'] for x in windows),'transitions':dict(collections.Counter(x['A_'+ch+'_status']+'->'+x['B_'+ch+'_status'] for x in windows))}
        for bound in ['lower','upper']:
            vals=[Decimal(x['B_minus_A_'+ch+'_'+bound+'_Ma']) for x in windows if x['B_minus_A_'+ch+'_'+bound+'_Ma']!='']
            delta_summary[ch][bound+'_delta_B_minus_A_Ma']={'min':str(min(vals)),'max':str(max(vals)),'max_abs':str(max(map(abs,vals)))}
    methods={'schema':1,'analysis_kind':'same_candidates_calibration_sensitivity','A':'T25_deep_only','B':'T25_three_calibrations','n_candidates':24,'host_tips_selected':25,'full_host_tips':495,'full_nodes':989,'retained_nodes':49,'mapping_identity':'exact full descendant tip set; IDs used only as local labels','age_anchor':'each original full-tree maximum root-to-tip distance; original tip rounding retained','decimal_precision':50,'A_reproduction_absolute_tolerance_Ma':str(TOL),'boundary_rule':'exact Decimal comparison, tolerance never used for classifying overlap','channels':{'clade':'overlap of incoming edges after pruning to same25 and contraction','terminal':'overlap of original full-tree immediate parent stems'},'interpretation':'B is calibration sensitivity, not a confidence/HPD interval, dating posterior, globally time-consistent reconciliation, or a transfer date','no_new_candidates_or_reconciliation':True,'sources':sources,'command':[sys.executable,*sys.argv],'software':{'Python':platform.python_version(),'Biopython':Bio.__version__}}
    out.mkdir(parents=True,exist_ok=False)
    tsv(out/'full495_clade_identity_AB.tsv',fullrows);tsv(out/'all49_branch_ages_AB.tsv',branch_rows);tsv(out/'all24_candidate_windows_AB.tsv',windows)
    summary={'status':'complete','data_kind':'real','full495_tip_and_rooted_clade_identity_exact':True,'complete_full_node_clades_compared':989,'retained_nodes':49,'retained_nonroot_edges':48,'candidate_count':24,'A_numeric_comparisons':len(err),'A_exact_comparisons':exact,'A_max_abs_reproduction_error_Ma':str(max(err)),'A_root_anchor_Ma':str(A['anchor']),'B_root_anchor_Ma':str(B['anchor']),'A_tip_age_range_Ma':[str(min(A['age'][n] for n in A['tree'].get_terminals())),str(max(A['age'][n] for n in A['tree'].get_terminals()))],'B_tip_age_range_Ma':[str(min(B['age'][n] for n in B['tree'].get_terminals())),str(max(B['age'][n] for n in B['tree'].get_terminals()))],'interval_status_counts':{s:{ch:dict(v) for ch,v in cc.items()} for s,cc in counts.items()},'deltas':delta_summary,'scientific_limitations':methods['interpretation']}
    save(out/'methods.json',methods);save(out/'summary.json',summary)
    report=['同24候选的校准A/B年龄敏感性已完成。两个完整495样本树的989个有根clade及父子身份完全一致；49个主25保留节点按完整后代集合映射，没有按N编号猜测。','A的全部既有节点/窗口数值核对：'+str(len(err))+'项，最大绝对误差 '+str(max(err))+' Ma。','']
    for ch,lab in [('clade','裁剪后保留分支交集'),('terminal','原完整树紧邻stem交集')]:
        report.extend([lab+'：A '+json.dumps(dict(counts['A'][ch]),ensure_ascii=False)+'；B '+json.dumps(dict(counts['B'][ch]),ensure_ascii=False)+'。','状态变化 '+str(delta_summary[ch]['status_changed_count'])+'/24；上下界B−A差值见 summary.json 与完整24行表。'])
    report.extend(['','B是校准设定敏感性，不是CI、HPD或后验区间；局部交集仍不是已证实转移日期或全协调时间一致性。原tip微小年龄保留，空交集上下界保留，没有截断到0。候选集合、主A图、原树、PACo及DTL分析均未修改。'])
    with (out/'review_zh.md').open('x',encoding='utf-8') as f:f.write('\n\n'.join(report)+'\n')
    for q,h in sources.items():assert sha(q)==h, 'Input changed during execution: '+q
    done={'schema':1,'status':'complete','data_kind':'real','stage':'calibration_AB_main25','finished_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'exit_code':0,'elapsed_seconds':time.monotonic()-started,'candidate_count':24,'host_node_count':49,'full_host_node_count':989,'command':methods['command'],'software':methods['software'],'source_sha256':sources,'outputs':[{'path':str(q.resolve()),'sha256':sha(q)} for q in sorted(out.iterdir()) if q.is_file()]}
    save(out/'completed.json',done)
    print(json.dumps({'status':'complete','out_dir':str(out),'completed_sha256':sha(out/'completed.json'),'script_sha256':sha(Path(__file__)),'elapsed_seconds':done['elapsed_seconds'],'A_max_abs_error_Ma':str(max(err)),'counts':summary['interval_status_counts'],'changes':delta_summary},ensure_ascii=False))
if __name__=='__main__': main()
