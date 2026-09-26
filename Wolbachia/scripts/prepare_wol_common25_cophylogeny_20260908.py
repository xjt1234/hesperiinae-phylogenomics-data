#!/usr/bin/env python3
"""Project eight already completed bacterial trees onto the frozen main25.
This prepares distance inputs only, without re-inferring or biologically rooting.
"""
import argparse, csv, datetime, hashlib, itertools, json, math, shutil, sys
from pathlib import Path
sys.dont_write_bytecode = True
from Bio import Phylo

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1048576), b''): h.update(block)
    return h.hexdigest()

def dump(path, data):
    with path.open('x') as f: json.dump(data, f, ensure_ascii=False, indent=2)

def check_tree(tree, expected):
    names = [x.name for x in tree.get_terminals()]
    assert len(names) == len(set(names)) and set(names) == set(expected)
    for node in tree.find_clades():
        assert node.branch_length is None or (math.isfinite(node.branch_length) and node.branch_length >= 0)

def main(run):
    confirmation = run/'01_PROVENANCE/host_source_author_confirmation_20260908.json'
    assert sha(confirmation) == 'd2bb0a2bd2c9fea217e2e43f4894fcd93f07cf061a877bc168a00c0ba984bf77'
    statement = json.loads(confirmation.read_text())
    taxa = sorted(row['wolbachia_sequence_id'] for row in statement['approved_context_pairs'])
    assert len(taxa) == len(set(taxa)) == 25
    frozen_input = run/'01_PROVENANCE/frozen_tree_inputs.json'
    assert set(json.loads(frozen_input.read_text())['scenarios'][0]['expected_taxa']) == set(taxa)
    spec = [
        ('main_mfp', 'coding_trees/main', 'coding_tree_main', 25, 'nucleotide substitutions/site', 'primary'),
        ('sample_expanded_mfp_project25', 'coding_trees/sample_expanded', 'coding_tree_sample_expanded', 31, 'nucleotide substitutions/site', 'sensitivity'),
        ('reference_repeat_excluded_mfp', 'coding_trees/reference_repeat_excluded', 'coding_tree_reference_repeat_excluded', 25, 'nucleotide substitutions/site', 'sensitivity'),
        ('main_gtr_gamma', 'coding_tree_diagnostics_20260908/main_gtr_gamma', 'coding_diagnostic_main_gtr_gamma', 25, 'nucleotide substitutions/site', 'sensitivity'),
        ('sample_expanded_gtr_gamma_project25', 'coding_tree_diagnostics_20260908/sample_expanded_gtr_gamma', 'coding_diagnostic_sample_expanded_gtr_gamma', 31, 'nucleotide substitutions/site', 'sensitivity'),
        ('reference_repeat_excluded_gtr_gamma', 'coding_tree_diagnostics_20260908/reference_repeat_excluded_gtr_gamma', 'coding_diagnostic_reference_repeat_excluded_gtr_gamma', 25, 'nucleotide substitutions/site', 'sensitivity'),
        ('main_aa_gamma', 'coding_tree_diagnostics_20260908/main_aa_gamma', 'coding_diagnostic_main_aa_gamma', 25, 'amino acid substitutions/site', 'sensitivity'),
        ('main_nt12_gamma', 'coding_tree_diagnostics_20260908/main_nt12_gamma', 'coding_diagnostic_main_nt12_gamma', 25, 'nucleotide substitutions/site at codon positions 1+2', 'sensitivity'),
    ]
    source_hashes = {str(confirmation):sha(confirmation), str(frozen_input):sha(frozen_input)}
    for label, directory, stage, count, unit, role in spec:
        source = run/'03_MATRICES_AND_TREES'/directory/'tree.treefile'
        checkpoint = run/'checkpoints'/('homology_'+stage+'.complete.json')
        completed = json.loads(checkpoint.read_text()); assert completed['status'] == 'complete'
        verified = {r['path']:r['sha256'] for r in completed['outputs']}
        assert str(source) in verified
        for p, expected in verified.items():
            assert sha(p) == expected, p
            source_hashes[p] = expected
        source_hashes[str(checkpoint)] = sha(checkpoint)
    out = run/'03_MATRICES_AND_TREES/cophylogeny_wol_common25_20260908'
    out.mkdir(exist_ok=False)
    scenarios=[]
    for label, directory, stage, count, unit, role in spec:
        source=run/'03_MATRICES_AND_TREES'/directory/'tree.treefile'
        original=Phylo.read(source,'newick'); names=[x.name for x in original.get_terminals()]
        assert len(names)==count and set(taxa)<=set(names)
        check_tree(original,names)
        target=out/label;target.mkdir();treepath=target/'wol_main25.newick'
        if count==25: shutil.copyfile(source,treepath)
        else:
            projected=Phylo.read(source,'newick')
            for tip in sorted(set(names)-set(taxa)): projected.prune(tip)
            Phylo.write(projected,treepath,'newick',format_branch_length='%1.17g',format_confidence='%1.17g')
        projected=Phylo.read(treepath,'newick'); check_tree(projected,taxa)
        pair_rows=[]
        for a,b in itertools.combinations(taxa,2):
            full=original.distance(a,b);new=projected.distance(a,b)
            assert math.isfinite(new) and new>=0 and math.isclose(full,new,rel_tol=1e-12,abs_tol=1e-12)
            pair_rows.append({'taxon_a':a,'taxon_b':b,'source_distance':full,'projected_distance':new,'absolute_difference':abs(full-new)})
        with (target/'projection_distance_QA.tsv').open('x',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(pair_rows[0]),delimiter='\t');w.writeheader();w.writerows(pair_rows)
        with (target/'patristic_distances.tsv').open('x',newline='') as f:
            w=csv.writer(f,delimiter='\t');w.writerow(['wol_tip']+taxa)
            for a in taxa:w.writerow([a]+[format(projected.distance(a,b),'.17g') for b in taxa])
        entry={'label':label,'role':role,'source_tree':str(source),'source_tree_sha256':sha(source),'tree':str(treepath),'tree_sha256':sha(treepath),'source_inference_taxa':count,'associated_taxa':25,'projection_pairs_validated':300,'max_pairwise_distance_error':max(r['absolute_difference'] for r in pair_rows),'distance_unit':unit,'biological_rerooting_performed':False,'support_labels_on_projected_tree':'not recomputed; no support conclusions derived from these labels'}
        dump(target/'projection.json',entry);scenarios.append(entry)
    manifest={'schema':1,'status':'complete','stage':'cophylogeny_wol_common25_inputs','main_analysis':'main_mfp','pre_fit_selection':'Keep the originally frozen main25/52 MFP+MERGE tree as primary; all other completed trees are sensitivity scenarios, selected before new cophylogeny statistics.','exact_wol_taxa':taxa,'scenarios':scenarios,'source_hashes':source_hashes,'script_sha256':sha(__file__),'date':datetime.datetime.now(datetime.timezone.utc).isoformat()}
    dump(out/'manifest.json',manifest)
    assert all(sha(p)==expected for p,expected in source_hashes.items())
    outputs=[{'path':str(p),'sha256':sha(p)} for p in sorted(out.rglob('*')) if p.is_file()]
    checkpoint={'status':'complete','stage':'cophylogeny_wol_common25_inputs','summary':{'scenarios':8,'associated_taxa_each':25,'pairwise_projection_checks':2400,'new_tree_inference':False,'new_statistics_fitted':False},'outputs':outputs,'source_hashes':source_hashes,'script_sha256':sha(__file__)}
    dump(run/'checkpoints/cophylogeny_wol_common25_inputs.complete.json',checkpoint)
    print(json.dumps(checkpoint['summary']))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run-dir',type=Path,required=True);args=p.parse_args();main(args.run_dir.resolve())
