#!/usr/bin/env python3
"""Two explicit local age intersections; never estimates transfer dates."""
import argparse
import csv
import datetime
import hashlib
import json
from decimal import Decimal
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rows(path):
    with Path(path).open() as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


WINDOW_FIELDS = ['candidate_id','lineage_a_label','lineage_b_label','lineage_a_tips_json','lineage_b_tips_json','clade_lower_Ma','clade_upper_Ma','clade_status','terminal_lower_Ma','terminal_upper_Ma','terminal_status','lineage_a_stem_age_Ma','lineage_b_stem_age_Ma','evidence_note','support_definition','plot_order']
AUDIT_FIELDS = ['candidate_id','side','host_branch_id','host_node_id','tips','full_node_id','full_parent_node_id','child_age','original_parent_age','retained_parent_id','retained_parent_age']

def write_tsv(path, data, fields):
    with Path(path).open('x') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(data)


def intersection(child_a, parent_a, child_b, parent_b):
    if any(x is None or x == '' for x in (child_a, parent_a, child_b, parent_b)):
        return '', '', 'not_available'
    values = [Decimal(str(x)) for x in (child_a, parent_a, child_b, parent_b)]
    assert all(x.is_finite() and x >= 0 for x in values)
    ca, pa, cb, pb = values
    assert pa >= ca and pb >= cb, 'A branch itself has reversed ages'
    low, high = max(ca, cb), min(pa, pb)
    state = 'feasible' if low < high else ('boundary_only' if low == high else 'infeasible')
    return str(low), str(high), state


def full_ages(full_rows):
    by_id = {r['full_node_id']: r for r in full_rows}
    assert full_rows and len(by_id) == len(full_rows)
    roots = [r for r in full_rows if not r['full_parent_node_id']]
    assert len(roots) == 1
    anchor = Decimal(roots[0]['full_age_from_full_tree_anchor'])
    for row in full_rows:
        age = Decimal(row['full_age_from_full_tree_anchor'])
        assert age.is_finite() and age >= 0
        assert age + Decimal(row['full_root_distance']) == anchor
        parent = row['full_parent_node_id']
        if parent:
            parent_age = Decimal(by_id[parent]['full_age_from_full_tree_anchor'])
            assert parent_age >= age
            assert parent_age - age == Decimal(row['full_incoming_branch_length'])
            assert Decimal(row['full_incoming_branch_length']).is_finite()
    return by_id


def resolve_branches(branch_rows, full_rows, associations):
    full = full_ages(full_rows)
    selected = {r['host_tip'] for r in associations}
    assert len(selected) == len(associations) == len({r['wol_tip'] for r in associations}) == 25
    for r in full_rows:
        tips = set(r['full_descendant_tips'].split(';')) & selected
        assert tips == set(filter(None,r['selected_descendant_tips'].split(';')))
        assert len(tips) == int(r['selected_descendant_count'])
    branch_source = {r['host_branch_id']:r for r in branch_rows}
    branch = {}
    for row in branch_rows:
        key = row['host_branch_id']
        assert key not in branch
        tips = frozenset(json.loads(row['selected_descendant_tips_json']))
        assert tips and tips <= selected
        original = full[row['full_target_node_id']]
        assert tips == set(original['selected_descendant_tips'].split(';'))
        assert original['pruning_status'] == 'retained_node', 'Target must be the retained full-tree node, not another node sharing selected descendants'
        assert row['host_node_id'] == key
        clade_hash = hashlib.sha256(json.dumps(sorted(tips),separators=(',',':')).encode()).hexdigest()
        assert row['clade_hash'] == clade_hash
        assert key == (next(iter(tips)) if len(tips)==1 else 'H_'+clade_hash)
        assert row['is_tip'] == str(int(len(tips)==1))
        assert row['full_immediate_parent_node_id'] == original['full_parent_node_id']
        assert Decimal(row['node_age_Ma']) == Decimal(original['full_age_from_full_tree_anchor'])
        path = json.loads(row['full_path_node_ids_json'])
        assert path == original['pruned_incoming_full_path_node_ids'].split(';')
        assert path and path[-1] == original['full_node_id']
        assert Decimal(row['sampled_incoming_branch_length_Ma']) == Decimal(original['pruned_incoming_branch_length'])
        segments = json.loads(row['full_path_segments_json'], parse_float=Decimal)
        assert len(segments) == len(path)
        for fid,segment in zip(path,segments):
            fr=full[fid];pid=fr['full_parent_node_id'];assert pid
            assert segment['full_child_node_id']==fid and segment['full_parent_node_id']==pid
            assert Decimal(str(segment['child_age_Ma']))==Decimal(fr['full_age_from_full_tree_anchor'])
            assert Decimal(str(segment['parent_age_Ma']))==Decimal(full[pid]['full_age_from_full_tree_anchor'])
            assert Decimal(str(segment['branch_length_Ma']))==Decimal(fr['full_incoming_branch_length'])
        for above,below in zip(path,path[1:]):assert full[below]['full_parent_node_id']==above
        assert sum(Decimal(full[fid]['full_incoming_branch_length']) for fid in path)==Decimal(row['sampled_incoming_branch_length_Ma'])
        original_parent = original['full_parent_node_id']
        assert original_parent and Decimal(row['original_immediate_parent_age_Ma'])==Decimal(full[original_parent]['full_age_from_full_tree_anchor'])
        branch[key] = {
            'host_branch_id': key,
            'host_node_id': row['host_node_id'],
            'tips': tips,
            'full_node_id': original['full_node_id'],
            'full_parent_node_id': original_parent,
            'child_age': Decimal(original['full_age_from_full_tree_anchor']),
            'original_parent_age': Decimal(full[original_parent]['full_age_from_full_tree_anchor']) if original_parent else None,
        }
    assert len(branch) == 49, 'Expected rooted binary main25 host tree, with unary full root contracted'
    assert len({v['tips'] for v in branch.values()}) == len(branch)
    for value in branch.values():
        supersets = [x for x in branch.values() if value['tips'] < x['tips']]
        if supersets:
            parent = min(supersets, key=lambda x: len(x['tips']))
            assert sum(len(x['tips']) == len(parent['tips']) for x in supersets) == 1
            value['retained_parent_id'] = parent['host_branch_id']
            value['retained_parent_age'] = parent['child_age']
            assert value['retained_parent_age'] >= value['original_parent_age'] >= value['child_age']
            row=branch_source[value['host_branch_id']]
            assert row['retained_parent_node_id']==parent['host_node_id']
            assert row['retained_parent_full_node_id']==parent['full_node_id']
            assert Decimal(row['retained_parent_age_Ma'])==value['retained_parent_age']
            assert row['is_sampled_root']=='0' and row['root_stem_status']=='not_root'
            path=json.loads(row['full_path_node_ids_json'])
            assert full[path[0]]['full_parent_node_id']==parent['full_node_id']
        else:
            assert value['tips'] == selected
            value['retained_parent_id'] = ''
            value['retained_parent_age'] = None
            row=branch_source[value['host_branch_id']]
            assert not row['retained_parent_node_id'] and not row['retained_parent_full_node_id'] and not row['retained_parent_age_Ma']
            assert row['is_sampled_root']=='1' and row['root_stem_status']=='known_finite_source_stem_saved_but_no_sampled_parent'
    assert sum(not v['retained_parent_id'] for v in branch.values())==1
    return branch


def make_windows(candidate_rows, branch, associations):
    names = {r['host_tip']: r['wol_tip'].replace('_', ' ') for r in associations}
    windows = []
    branch_audit = []
    seen = set()
    seen_pairs = set()
    for row in candidate_rows:
        cid = row['candidate_id']
        assert cid and cid not in seen
        seen.add(cid)
        a, b = branch[row['host_branch_a']], branch[row['host_branch_b']]
        pair=tuple(sorted([a['host_branch_id'],b['host_branch_id']]))
        assert pair not in seen_pairs, 'Duplicate branch pair under different candidate IDs'
        seen_pairs.add(pair)
        for side,value in [('a',a),('b',b)]:
            if 'host_'+side+'_tips_json' in row:assert set(json.loads(row['host_'+side+'_tips_json']))==value['tips']
        assert a['tips'].isdisjoint(b['tips']), 'Transfer landing cannot be an ancestor/descendant branch'
        gray = intersection(a['child_age'], a['retained_parent_age'], b['child_age'], b['retained_parent_age'])
        blue = intersection(a['child_age'], a['original_parent_age'], b['child_age'], b['original_parent_age'])
        if gray[2] != 'not_available' and blue[2] != 'not_available':
            assert Decimal(gray[0]) == Decimal(blue[0])
            assert Decimal(gray[1]) >= Decimal(blue[1])
        def label(x):
            return names[next(iter(x['tips']))] if len(x['tips']) == 1 else f"Clade {x['full_node_id']} (n={len(x['tips'])})"
        data = {
            'candidate_id': cid,
            'lineage_a_label': label(a), 'lineage_b_label': label(b),
            'lineage_a_tips_json': json.dumps(sorted(a['tips'])),
            'lineage_b_tips_json': json.dumps(sorted(b['tips'])),
            'clade_lower_Ma': gray[0], 'clade_upper_Ma': gray[1], 'clade_status': gray[2],
            'terminal_lower_Ma': blue[0], 'terminal_upper_Ma': blue[1], 'terminal_status': blue[2],
            'lineage_a_stem_age_Ma': str(a['original_parent_age']) if a['original_parent_age'] is not None else '',
            'lineage_b_stem_age_Ma': str(b['original_parent_age']) if b['original_parent_age'] is not None else '',
            'evidence_note': 'Gray: overlap of incoming edges in the actual pruned reconciliation host tree. Blue: overlap restricted to immediate stem segments above the corresponding nodes in the original A495 tree. Blue is a narrower conditional scenario, not a correction of gray. Both are local necessary conditions, not a globally time-consistent reconciliation or transfer dates. Tip rounding is retained, lower bounds are computed and never replaced with zero. Root with no defined incoming edge is unavailable.',
            'support_definition': 'Candidate set is the union of every transfer branch pair in three deterministic optimal medians at the predeclared display root of main MFP. Exact conditional at-least-once MPR counts and denominators by model, edge root and cost are separate source data; roots/models are sensitivity conditions, not probability samples.',
            'plot_order': len(windows) + 1,
        }
        windows.append(data)
        for side, value in [('a', a), ('b', b)]:
            branch_audit.append({'candidate_id': cid, 'side': side, **{
                k: (json.dumps(sorted(v)) if k == 'tips' else (str(v) if isinstance(v, Decimal) else v))
                for k, v in value.items()
            }})
    return windows, branch_audit


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir', type=Path, required=True)
    parser.add_argument('--reconciliation-checkpoint', type=Path, required=True)
    parser.add_argument('--host-branch-map', type=Path, required=True)
    parser.add_argument('--candidate-tsv', type=Path, required=True)
    parser.add_argument('--out-dir', type=Path, required=True)
    args = parser.parse_args()
    r = args.run_dir.resolve()
    map_path = r/'03_MATRICES_AND_TREES/cophylogeny_inputs_main25_20260908/AHE_A495_main25.full_to_pruned_node_map.tsv'
    assoc_path = r/'03_MATRICES_AND_TREES/cophylogeny_inputs_main25_20260908/AHE_main25_associations.tsv'
    host_cp = r/'checkpoints/cophylogeny_inputs_main25_20260908.covered_evidence.complete.json'
    cp = json.loads(args.reconciliation_checkpoint.read_text())
    assert cp['status'] == 'complete'
    coverage = {x['path']: x['sha256'] for x in cp['outputs']}
    for path, digest in coverage.items():
        assert sha(path) == digest
    for path in [args.host_branch_map.resolve(), args.candidate_tsv.resolve()]:
        assert coverage[str(path)] == sha(path)
    hcp = json.loads(host_cp.read_text())
    assert hcp['status'] == 'complete'
    for path, digest in hcp['output_sha256'].items():
        assert sha(path) == digest
    assert hcp['output_sha256'][str(map_path)] == sha(map_path)
    assert hcp['output_sha256'][str(assoc_path)] == sha(assoc_path)
    source = {str(x.resolve()): sha(x) for x in [Path(__file__), args.reconciliation_checkpoint, args.host_branch_map, args.candidate_tsv, map_path, assoc_path, host_cp]}
    branches = resolve_branches(rows(args.host_branch_map), rows(map_path), rows(assoc_path))
    windows, branch_audit = make_windows(rows(args.candidate_tsv), branches, rows(assoc_path))
    args.out_dir.mkdir(parents=True, exist_ok=False)
    write_tsv(args.out_dir/'conditional_timewindows.tsv', windows, WINDOW_FIELDS)
    write_tsv(args.out_dir/'conditional_branch_age_audit.tsv', branch_audit, AUDIT_FIELDS)
    (args.out_dir/'methods.json').write_text(json.dumps({'schema': 1, 'window_channel_labels': {'clade': 'Pruned-branch overlap', 'terminal': 'Original-stem overlap'}, 'local_necessary_conditions_only': True, 'original_tip_rounding_retained': True, 'biological_direction_inferred': False, 'age_anchor': 'Original full495 maximum root-to-tip distance; preserved source node table, no recalibration', 'sources': source}, indent=2)+'\n')
    for path, digest in source.items():
        assert sha(path) == digest
    done = {'status': 'complete', 'finished_at': datetime.datetime.now().astimezone().isoformat(), 'candidate_count': len(windows), 'candidate_set_status': 'nonempty' if windows else 'empty_completed_inference_no_candidates_invented', 'host_branches_verified_against_full_A495': len(branches), 'source_sha256': source, 'interval_status_counts': {prefix: {state: sum(x[prefix+'_status'] == state for x in windows) for state in ('feasible', 'boundary_only', 'infeasible', 'not_available')} for prefix in ('clade', 'terminal')}, 'outputs': [{'path': str(p.resolve()), 'sha256': sha(p)} for p in sorted(args.out_dir.iterdir()) if p.is_file()]}
    with (args.out_dir/'completed.json').open('x') as stream:
        json.dump(done, stream, indent=2)
    print(json.dumps(done))


if __name__ == '__main__':
    main()
