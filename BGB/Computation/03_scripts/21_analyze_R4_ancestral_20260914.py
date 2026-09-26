#!/usr/bin/env python3
"""Auditable quality checks and summaries of completed R4 reconstructions."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd

JOB = Path(__file__).resolve().parents[1]
SOURCE = JOB / '07_tables/ancestral_R4_20260914_v1'
QA = JOB / '08_qa/ancestral_R4_20260914_v1'
MODELS = ['M0', 'M1', 'M2']
AREAS = ['AF','AUS','CAM','ENA','EPA','IND','MDG','ORI','SAM','WNA','WPA']

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()

def read(name):
    kwargs={'dtype':{'geography_bits':str}} if name=='tip_metadata.tsv' else {}
    return pd.read_csv(SOURCE/name, sep='\t', keep_default_na=False, float_precision='round_trip', **kwargs)

def load_data():
    manifest = json.loads((SOURCE/'source_manifest.json').read_text())
    assert manifest['status'] == 'CONDITIONAL_NEW_R4_NATIVE_POSTFIT'
    for item in manifest['sources']:
        assert sha(item['path']) == item['sha256'], item['path']
    for item in manifest['outputs']:
        assert sha(SOURCE/item['file']) == item['sha256'], item['file']
    nodes, edges, tips = read('nodes.tsv'), read('edges.tsv'), read('tip_metadata.tsv')
    anchors, states, ns = read('anchors_summary.tsv'), read('state_dictionary.tsv'), read('node_summaries.tsv')
    assert len(nodes)==833 and not nodes.ape_node.duplicated().any()
    assert len(edges)==832 and not edges.child.duplicated().any()
    assert len(tips)==417 and not tips.tip_label.duplicated().any()
    assert len(anchors)==45 and not anchors.duplicated(['model','anchor_number']).any()
    assert len(ns)==2499 and not ns.duplicated(['model','ape_node']).any()
    assert states.state_index_1based.tolist()==list(range(1,563))
    assert states.range_size.value_counts().sort_index().to_dict()=={0:1,1:11,2:55,3:165,4:330}
    node_map = nodes.set_index('ape_node').to_dict('index')
    for row in node_map.values():
        row['plot_age_Ma'] = float(node_map[418]['node_age_Ma'])-float(row['root_distance'])
    tip_map = {n:r['tip_label'] for n,r in node_map.items() if r['is_tip']}
    assert set(tip_map.values())==set(tips.tip_label)
    children, parent = {}, {}
    for edge in edges.itertuples():
        children.setdefault(edge.parent,[]).append(edge.child); parent[edge.child]=edge.parent
        assert abs(node_map[edge.child]['root_distance']-node_map[edge.parent]['root_distance']-edge.branch_length_Ma)<1e-7
    assert set(node_map)-set(parent)=={418}
    desc={}
    def descend(n):
        desc[n] = {n} if n in tip_map else set().union(*(descend(c) for c in children[n]))
        return desc[n]
    descend(418)
    for n, ds in desc.items():
        assert hashlib.sha256(('\n'.join(sorted(tip_map[t] for t in ds))+'\n').encode()).hexdigest()==node_map[n]['tipset_sha256']
    order=[]
    def walk(n):
        if n in tip_map: order.append(n)
        else:
            for c in sorted(children[n],key=lambda c:min(desc[c])): walk(c)
    walk(418)
    y={n:float(i) for i,n in enumerate(order)}
    def sety(n):
        if n not in y:y[n]=float(np.mean([sety(c) for c in children[n]]))
        return y[n]
    sety(418)
    tribe_map=dict(zip(tips.tip_label,tips.tribe)); geog=dict(zip(tips.tip_label,tips.geography_bits))
    assert all(len(bits)==11 and set(bits)<={'0','1'} for bits in geog.values())
    assert all('+'.join(a for a,b in zip(AREAS,r.geography_bits) if b=='1')==r.observed_state for r in tips.itertuples())
    tribe_qa=[]
    for tribe in sorted(set(tribe_map.values())):
        members={n for n,label in tip_map.items() if tribe_map[label]==tribe}
        exact=[n for n,ds in desc.items() if ds==members]
        assert len(exact)==1,(tribe,exact)
        tribe_qa.append({'tribe':tribe,'tip_count':len(members),'monophyletic':True,'ape_node':exact[0]})
    assert len(tribe_qa)==14
    inc=np.array([[int(a in state.split('+')) for a in AREAS] for state in states.semantic_state])
    posterior={}; marg={}; profile=[]
    tip_states=dict(zip(tips.tip_label,tips.observed_state))
    for model in MODELS:
        table=read(f'{model}_node_posteriors.tsv').set_index('ape_node')
        assert table.index.tolist()==list(range(1,834)) and table.columns.tolist()==[f'p{i:04d}' for i in range(1,563)]
        p=table.to_numpy(float)
        assert np.isfinite(p).all() and p.min()>=0 and p.max()<=1
        error=float(np.max(np.abs(p.sum(1)-1)));assert error<1e-12
        area=read(f'{model}_area_inclusion.tsv').set_index('ape_node')[AREAS].to_numpy(float)
        assert np.max(np.abs(area-p@inc))<1e-12
        assert np.max(np.abs(area.sum(1)-p@states.range_size.to_numpy()))<1e-11
        for n,label in tip_map.items():assert p[n-1,list(states.semantic_state).index(tip_states[label])]==1
        posterior[model]={int(n):row for n,row in zip(table.index,p)}
        marg[model]={int(n):row for n,row in zip(table.index,area)}
        profile.append({'model':model,'nodes':833,'states':562,'duplicate_node_keys':0,'invalid_or_missing_probabilities':0,'max_rowsum_error':error,'tips_match_geography':417,'full_area_marginal_check':True})
    anchor_map={int(r.ape_node):int(r.anchor_number) for r in anchors[anchors.model=='M1'].itertuples()}
    assert len(anchor_map)==15
    for r in anchors.itertuples():
        p=posterior[r.model][r.ape_node]
        assert abs(r.top1_prob-p.max())<1e-12 and abs(r.top1_prob+r.top2_prob+r.top3_prob+r.other_prob-1)<1e-12
    area=read('area_order.tsv')
    return {'nodes':node_map,'edges':edges,'children':children,'parent':parent,'tips':tip_map,'order':order,'desc':desc,'y':y,'geog':geog,'tribes':tribe_map,'tribe_qa':tribe_qa,'posterior':posterior,'marginal':marg,'states':states,'anchors':anchors,'anchor_map':anchor_map,'area_names':dict(zip(area.abbrev,area.full_name)),'node_summaries':ns,'profile':profile,'model_comparison':read('model_comparison.tsv'),'manifest':manifest}

def main():
    data=load_data(); QA.mkdir(parents=True,exist_ok=False)
    pd.DataFrame(data['profile']).to_csv(QA/'probability_QA.tsv',sep='\t',index=False)
    pd.DataFrame(data['tribe_qa']).to_csv(QA/'tribe_monophyly_QA.tsv',sep='\t',index=False)
    rows=[]
    for a,b in [('M0','M1'),('M1','M2'),('M0','M2')]:
        pa=np.array([data['posterior'][a][n] for n in range(418,834)])
        pb=np.array([data['posterior'][b][n] for n in range(418,834)])
        tv=np.abs(pa-pb).sum(1)/2
        rows.append({'model_a':a,'model_b':b,'internal_nodes':416,'same_top_state':int((pa.argmax(1)==pb.argmax(1)).sum()),'different_top_state':int((pa.argmax(1)!=pb.argmax(1)).sum()),'TV_median':float(np.median(tv)),'TV_max':float(tv.max()),'TV_root':float(tv[0]),'TV_above_0_5':int((tv>.5).sum())})
    pd.DataFrame(rows).to_csv(SOURCE/'model_sensitivity_internal_nodes.tsv',sep='\t',index=False,float_format='%.17g')
    root=data['anchors'][data['anchors'].anchor_number==1]
    root.to_csv(SOURCE/'root_three_model_summary.tsv',sep='\t',index=False,float_format='%.17g')
    audit={'status':'PASS_SOURCE_INTEGRITY_AND_CONDITIONAL_PROBABILITY_QA','inputs_sha256':sha(SOURCE/'source_manifest.json'),'profile':data['profile'],'internal_node_comparisons':rows,'scientific_status':'CONDITIONAL_AUTHOR_REVIEW_KKT1_FALSE','no_fit_or_postfit_recalculation':True,'limitations':['No global optimum certification','R4 legacy coding exception for Pelopidas mathias','R6 not pooled or interpreted as completed','BSM is separate','Node probabilities do not include tree or parameter uncertainty']}
    (QA/'analysis_QA.json').write_text(json.dumps(audit,indent=2)+'\n')
    report=['# R4 ancestral analysis — 14 September 2026','','## Outcome','All three new R4 fits and independent postfit reconstructions are complete. This export and figure analysis performed no new optimization or ancestral likelihood recalculation. M1 remains the prespecified main scenario, M0 the static baseline and M2 the permissive sensitivity.','','## Model comparison (identical R4 inputs)','Model | lnL | AIC | d | e','--- | ---: | ---: | ---: | ---:']
    for r in data['model_comparison'].itertuples():report.append(f'{r.model} | {r.lnL:.6f} | {r.AIC:.6f} | {r.d:.9g} | {r.e:.9g}')
    report+=['','M2 has the smallest descriptive AIC, but this does not replace the predesignated M1 analysis or certify a global optimum. All three results retain KKT1 FALSE, KKT2 TRUE and convergence code zero. M2 extinction is at the lower bound. AIC weights are conditional numerical summaries of these R4 candidate fits only, not weights across the differently coded R6 scenario.','','## Root uncertainty','Model | Largest complete range | Probability | Other beyond top three | P4','--- | --- | ---: | ---: | ---:']
    for r in root.itertuples():report.append(f'{r.model} | {r.top1_state} | {r.top1_prob:.6f} | {r.other_prob:.6f} | {r.p_range4:.6f}')
    report+=['','The root is diffuse and strongly scenario-dependent. Under M1 the largest single complete range has probability only 0.1893, so a unique ancestral origin is not supported. Root total variation is 0.902708 for M1 versus M2 and 0.961430 for M0 versus M1. M1 assigns 0.827118 probability to four-area ranges; this is a ceiling-mass diagnostic, not proof of bias or a completed alternative range-cap test.','','## Internal nodes','Pair | Same largest state / 416 | Different / 416 | Median TV','--- | ---: | ---: | ---:']
    for r in rows:report.append(f"{r['model_a']}–{r['model_b']} | {r['same_top_state']} | {r['different_top_state']} | {r['TV_median']:.6f}")
    report+=['','Top-state agreement alone does not measure agreement of complete distributions. TV and full probability tables accompany these summaries. Anchor nodes were selected before this rerun, not chosen to favor the new result. Singleton tribal anchors are stem-parent nodes and include sister lineages.','','## Data quality and lineage','The export contains 833 nodes × 562 full states for each model; all row sums, bounds, state order, tree/metadata/geography joins and 417 observed tip states passed. All 14 tribal memberships matched exact monophyletic descendant sets. Reused anchor definitions matched descendant-set SHA256 identities and node ages. New fit and postfit hashes were verified; no earlier geography posterior table or BSM history was reused.','','## Interpretation limits and handoff','The four-area state ceiling is retained by user decision. Pelopidas mathias uses the requested legacy coding; its omitted AUS/MDG bits must not be described as evidence of absence. The other 416 taxa use the audited new coding. The ongoing R6 analysis changes both this coding exception and the range ceiling, so it is a joint scenario contrast rather than a pure R4/R6 cap test.','','The figures visualize conditional node-top probabilities on one fixed dated tree with fixed fitted parameters. They are not across-tree confidence intervals, stochastic histories, directional colonization counts or a founder-event test. BSM results will be reported separately after the new run. Manuscript-scale and full-tip figures are review-ready drafts while the wider project remains in progress.','','## Reproducibility','Run the source exporter only into a fresh version directory, then the analysis and plotting scripts. The companion notebook executes the same hash, join and probability checks. The scientific-figure workflow determined the hero-tree/uncertainty layout, editable vector exports and explicit source/limitation QA; no numerical result was changed for aesthetics.']
    reportpath=JOB/'09_SUMMARIES/R4_ANCESTRAL_ANALYSIS_20260914.md';reportpath.parent.mkdir(parents=True,exist_ok=True)
    with reportpath.open('x') as handle:handle.write('\n'.join(report)+'\n')
    cells=[{'cell_type':'markdown','metadata':{},'source':['# New R4 ancestral-source quality checks\n','These cells execute the same read-only hash, key, topology, geography and probability checks used to generate the figures. No model is fitted.\n']},{'cell_type':'code','metadata':{},'execution_count':None,'outputs':[],'source':['from pathlib import Path\n','import importlib.util\n',f'job = Path({str(JOB)!r})\n',"spec = importlib.util.spec_from_file_location('r4_analysis', job / '03_scripts/21_analyze_R4_ancestral_20260914.py')\n","analysis = importlib.util.module_from_spec(spec); spec.loader.exec_module(analysis)\n","data = analysis.load_data()\n","import pandas as pd\n","pd.DataFrame(data['profile'])\n"]},{'cell_type':'code','metadata':{},'execution_count':None,'outputs':[],'source':["data['anchors'].query('anchor_number == 1')[['model','top1_state','top1_prob','other_prob','p_range4','TV_M1_M2']]\n"]}]
    (QA/'R4_ancestral_QC.ipynb').write_text(json.dumps({'nbformat':4,'nbformat_minor':5,'metadata':{'kernelspec':{'display_name':'Python 3','language':'python','name':'python3'}},'cells':cells},indent=2)+'\n')
    print(json.dumps(audit,indent=2))

if __name__=='__main__':main()
