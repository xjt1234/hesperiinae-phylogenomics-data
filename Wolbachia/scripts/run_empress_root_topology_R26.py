#!/usr/bin/env python3
"""Frozen all-edge-root DTL sensitivity; no biological outgroup or event dating.

--freeze only builds verified input/root catalogs and a protocol. The ordinary
--stage all entry is intended for the existing bounded background runner.
"""
import argparse
import copy
import csv
import datetime
import hashlib
import importlib.util
import io
import itertools
import json
import math
import os
from pathlib import Path
import resource
import signal
import sys
import time
import traceback

EMPRESS_PYTHON = '/home/data/t200301/.local/share/virtualenvs/empress-rWZu1KWm/bin/python'
EMPRESS_SOURCE = Path('/home/data/t200301/software/empress')
PLAN_REL = '01_PROVENANCE/empress_root_topology_protocol_20260908.json'
OUT_REL = '05_RECONCILIATION_AND_TIMING/empress_root_topology_20260908'
SINK = (None, None)
COSTS = [(2, 2, 1), (2, 3, 1), (2, 4, 1)]


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def canonical(obj):
    return json.dumps(obj, ensure_ascii=True, separators=(',', ':'))


def digest(obj):
    return hashlib.sha256(canonical(obj).encode()).hexdigest()


def save_json(path, obj):
    with Path(path).open('x') as f:
        json.dump(obj, f, indent=2, ensure_ascii=False)
        f.write('\n')


def save_tsv(path, rows, fields):
    with Path(path).open('x', newline='') as f:
        w = csv.DictWriter(f, fields, delimiter='\t', extrasaction='raise')
        w.writeheader()
        w.writerows(rows)


def verified(spec):
    if sha(spec['path']) != spec['sha256']:
        raise ValueError('Frozen source changed: ' + spec['path'])
    return Path(spec['path'])


def file_spec(path):
    return {'path': str(Path(path).resolve()), 'sha256': sha(path)}


def load_adapter(spec):
    path = verified(spec)
    ms = importlib.util.spec_from_file_location('empress_graph_adapter_frozen', path)
    module = importlib.util.module_from_spec(ms)
    ms.loader.exec_module(module)
    return module


def newick(tree):
    h = io.StringIO()
    Phylo.write(tree, h, 'newick', format_branch_length='%1.17g')
    return h.getvalue().strip() + '\n'


def tipset(clade):
    return sorted(t.name for t in clade.get_terminals())


def node_names(tree, prefix):
    rows = []
    for c in tree.find_clades(order='postorder'):
        tips = tipset(c)
        if not c.is_terminal():
            c.name = prefix + digest(tips)
        c.confidence = None
        rows.append({'node_id': c.name, 'clade_hash': digest(tips), 'tips_json': canonical(tips),
                     'is_tip': int(c.is_terminal())})
    names = [x['node_id'] for x in rows]
    if len(names) != len(set(names)):
        raise ValueError('Duplicate node names')
    return rows


def require_binary(tree, expected_tips):
    names = [t.name for t in tree.get_terminals()]
    if len(names) != len(set(names)) or set(names) != set(expected_tips):
        raise ValueError('Exact tree tip set mismatch')
    for c in tree.find_clades():
        if not c.is_terminal() and len(c.clades) != 2:
            raise ValueError('eMPRess requires strictly binary rooted trees')


def unrooted_edges(tree):
    """Suppress degree-two display root only, preserve every topology edge."""
    adj = {c: {} for c in tree.find_clades()}
    for parent in tree.find_clades():
        for child in parent.clades:
            length = child.branch_length
            if length is None or not math.isfinite(length) or length < 0:
                raise ValueError('Missing/negative/nonfinite Wol edge length')
            adj[parent][child] = length
            adj[child][parent] = length
    r = tree.root
    suppressed = len(adj[r]) == 2
    if suppressed:
        a, b = list(adj[r])
        length = adj[r][a] + adj[r][b]
        del adj[a][r]
        del adj[b][r]
        del adj[r]
        adj[a][b] = adj[b][a] = length
    tips = {c for c in adj if c.is_terminal()}
    if any(len(adj[c]) != (1 if c in tips else 3) for c in adj):
        raise ValueError('Wol tree is not fully bifurcating unrooted topology')

    def side(a, b):
        stack, seen, labels = [a], {b}, []
        while stack:
            node = stack.pop()
            seen.add(node)
            if node in tips:
                labels.append(node.name)
            stack.extend(x for x in adj[node] if x not in seen)
        return sorted(labels)

    catalog, seen = [], set()
    for a in adj:
        for b, length in adj[a].items():
            if frozenset((a, b)) in seen:
                continue
            seen.add(frozenset((a, b)))
            sides = sorted([side(a, b), side(b, a)])
            key = canonical(sides)
            catalog.append({'a': a, 'b': b, 'length': length, 'sides': sides,
                            'canonical_split': key, 'root_edge_hash': digest(sides)})
    catalog.sort(key=lambda x: x['canonical_split'])
    if len(catalog) != 2 * len(tips) - 3:
        raise ValueError('Unexpected root edge count')
    return adj, catalog, suppressed


def rooted_edge(adj, edge):
    def branch(node, parent, length):
        children = [branch(c, node, value) for c, value in adj[node].items() if c != parent]
        children.sort(key=lambda c: canonical(tipset(c)))
        return Clade(branch_length=length, name=node.name if not children else None, clades=children)
    r = Clade(branch_length=0, clades=[branch(edge['a'], edge['b'], edge['length'] / 2),
                                      branch(edge['b'], edge['a'], edge['length'] / 2)])
    r.clades.sort(key=lambda c: canonical(tipset(c)))
    result = Tree(root=r, rooted=True)
    node_names(result, 'P_')
    return result


def pairwise(tree):
    tips = sorted(t.name for t in tree.get_terminals())
    return {(a, b): tree.distance(a, b) for a, b in itertools.combinations(tips, 2)}


def prepare_host(host_path, map_path, expected_tips, out):
    original = Phylo.read(host_path, 'newick')
    source_dist = pairwise(original)
    tree = copy.deepcopy(original)
    removed = []
    while len(tree.root.clades) == 1:
        removed.append({'source_root_name': tree.root.name,
                        'source_root_branch_length': tree.root.branch_length,
                        'removed_child_stem_length': tree.root.clades[0].branch_length})
        tree.root = tree.root.clades[0]
    if len(removed) != 1:
        raise ValueError('Expected exactly one known AHE source unary root')
    require_binary(tree, expected_tips)
    rows = node_names(tree, 'H_')
    parent = {c: p for p in tree.find_clades() for c in p.clades}
    with Path(map_path).open() as f:
        full_rows = list(csv.DictReader(f, delimiter='\t'))
    by_full = {r['full_node_id']: r for r in full_rows}
    retained = {tuple(r['selected_descendant_tips'].split(';')): r for r in full_rows
                if r['pruning_status'] == 'retained_node'}
    if len(retained) != 49:
        raise ValueError('Expected exactly 49 retained binary host nodes')
    outrows = []
    for clade in tree.find_clades():
        tips = tipset(clade)
        row = retained[tuple(tips)]
        full_path = row['pruned_incoming_full_path_node_ids'].split(';')
        segments = []
        for fid in full_path:
            fr = by_full[fid]
            pid = fr['full_parent_node_id']
            if not pid:
                raise ValueError('Missing finite original incoming branch')
            segments.append({'full_child_node_id': fid, 'full_parent_node_id': pid,
                             'child_age_Ma': float(fr['full_age_from_full_tree_anchor']),
                             'parent_age_Ma': float(by_full[pid]['full_age_from_full_tree_anchor']),
                             'branch_length_Ma': float(fr['full_incoming_branch_length'])})
        for previous, following in zip(segments, segments[1:]):
            if previous['full_child_node_id'] != following['full_parent_node_id']:
                raise ValueError('Full path is not contiguous')
        sampled_parent = parent.get(clade)
        parent_full = retained[tuple(tipset(sampled_parent))] if sampled_parent else None
        crow = {'host_node_id': clade.name, 'host_branch_id': clade.name,
                'selected_descendant_tips_json': canonical(tips), 'clade_hash': digest(tips),
                'is_tip': int(clade.is_terminal()), 'is_sampled_root': int(sampled_parent is None),
                'retained_parent_node_id': sampled_parent.name if sampled_parent else '',
                'full_target_node_id': row['full_node_id'],
                'full_immediate_parent_node_id': row['full_parent_node_id'],
                'retained_parent_full_node_id': parent_full['full_node_id'] if parent_full else '',
                'node_age_Ma': row['full_age_from_full_tree_anchor'],
                'retained_parent_age_Ma': parent_full['full_age_from_full_tree_anchor'] if parent_full else '',
                'original_immediate_parent_age_Ma': by_full[row['full_parent_node_id']]['full_age_from_full_tree_anchor'],
                'full_path_node_ids_json': canonical(full_path), 'full_path_segments_json': canonical(segments),
                'sampled_incoming_branch_length_Ma': row['pruned_incoming_branch_length'],
                'root_stem_status': 'known_finite_source_stem_saved_but_no_sampled_parent' if sampled_parent is None else 'not_root'}
        outrows.append(crow)
    tree.root.branch_length = 0
    if max(abs(v - pairwise(tree)[k]) for k, v in source_dist.items()) > 1e-8:
        raise ValueError('Host format transformation changed pairwise distances')
    host_out = out / 'host_binary_format.nwk'
    host_out.write_text(newick(tree))
    save_tsv(out / 'host_branch_map.tsv', outrows, list(outrows[0]))
    save_json(out / 'host_format_conversion.json', {'removed_unary_roots': removed,
              'distance_max_abs_error': max(abs(v - pairwise(tree)[k]) for k, v in source_dist.items()),
              'source_map': file_spec(map_path), 'source_host': file_spec(host_path),
              'dated_semantics': 'Use full A495 map; never eMPRess host_distances. Root virtual edge has no invented time upper bound.'})
    return tree, host_out


def freeze(run):
    manifest_path = run / '01_PROVENANCE/cophylogeny_main25_16scenarios_20260908.json'
    manifest = json.loads(manifest_path.read_text())
    scenarios = [s for s in manifest['scenarios'] if s['id'].startswith('AHE_A495_dated__')]
    if len(scenarios) != 8:
        raise ValueError('Expected eight frozen AHE / Wol common25 scenarios')
    srcs = [file_spec(manifest_path)]
    for s in scenarios:
        for key in ('host_tree', 'wol_tree', 'links'):
            verified(s['inputs'][key])
            srcs.append(s['inputs'][key])
    host_spec = scenarios[0]['inputs']['host_tree']
    link_spec = scenarios[0]['inputs']['links']
    if any(s['inputs']['host_tree'] != host_spec or s['inputs']['links'] != link_spec for s in scenarios):
        raise ValueError('Sensitivity scenarios do not share exact host and links')
    with verified(link_spec).open() as f:
        links = list(csv.DictReader(f, delimiter='\t'))
    htip = [x['host_tip'] for x in links]
    ptip = [x['wol_tip'] for x in links]
    if len(links) != 25 or len(set(htip)) != 25 or len(set(ptip)) != 25:
        raise ValueError('Expected complete 25-pair bijection in different ID domains')
    map_path = Path(host_spec['path']).with_name('AHE_A495_main25.full_to_pruned_node_map.tsv')
    srcs.append(file_spec(map_path))
    full_date = Path('/home/data/t200301/xjt/Hesperiinae_review/treePL_R2_3_T25_deep_only_20260904_000707/output/T25_deep_only_dated.tre')
    srcs.append(file_spec(full_date))
    out = run / OUT_REL
    out.mkdir(exist_ok=False)
    inputs = out / 'inputs'
    inputs.mkdir()
    host_tree, host_out = prepare_host(verified(host_spec), map_path, htip, out)
    mapping_path = out / 'exact_associations.mapping'
    mapping_path.write_text(''.join(p + ':' + h + '\n' for p, h in sorted(zip(ptip, htip))))
    prepared = [file_spec(host_out), file_spec(out / 'host_branch_map.tsv'), file_spec(mapping_path)]
    catalog_rows, models = [], []
    for s in scenarios:
        model = s['id'].split('__', 1)[1]
        tree = Phylo.read(verified(s['inputs']['wol_tree']), 'newick')
        if set(tipset(tree.root)) != set(ptip):
            raise ValueError('Wol exact tip mismatch')
        original_distances = pairwise(tree)
        adj, edges, suppressed = unrooted_edges(tree)
        if len(edges) != 47:
            raise ValueError('Expected 47 roots per 25-tip model')
        mdir = inputs / model
        mdir.mkdir()
        root_specs = []
        for rank, e in enumerate(edges, 1):
            rooted = rooted_edge(adj, e)
            require_binary(rooted, ptip)
            current_dist = pairwise(rooted)
            max_error = max(abs(v - current_dist[k]) for k, v in original_distances.items())
            if max_error > 1e-8:
                raise ValueError('Rerooting changed distances')
            path = mdir / ('root_' + e['root_edge_hash'] + '.nwk')
            path.write_text(newick(rooted))
            nodes_path = path.with_suffix('.nodes.json')
            save_json(nodes_path, node_names(rooted, 'P_'))
            spec = {'rank': rank, 'root_edge_hash': e['root_edge_hash'],
                    'tree': file_spec(path), 'nodes': file_spec(nodes_path)}
            root_specs.append(spec)
            prepared.extend([spec['tree'], spec['nodes']])
            catalog_rows.append({'model': model, 'root_rank': rank, 'root_edge_hash': e['root_edge_hash'],
                                 'canonical_split': e['canonical_split'], 'source_edge_length': e['length'],
                                 'zero_length_edge': int(e['length'] == 0),
                                 'display_root_degree_two_suppressed': int(suppressed),
                                 'rooting_pairwise_distance_max_abs_error': max_error,
                                 'rooted_tree_path': str(path), 'rooted_tree_sha256': spec['tree']['sha256']})
        models.append({'model': model, 'source': s['inputs']['wol_tree'], 'roots': root_specs})
    save_tsv(out / 'root_edge_catalog.tsv', catalog_rows, list(catalog_rows[0]))
    prepared.append(file_spec(out / 'root_edge_catalog.tsv'))
    helper = run / 'scripts/empress_graph_adapter_R26.py'
    frozen_helper = run / 'scripts/frozen' / ('empress_graph_adapter_R26_' + sha(helper)[:16] + '.py')
    with frozen_helper.open('xb') as f:
        f.write(helper.read_bytes())
    source_files = [file_spec(p) for p in sorted((EMPRESS_SOURCE / 'empress').rglob('*.py'))]
    plan = {'schema': 1, 'data_kind': 'real', 'status': 'frozen_before_real_reconciliation', 'created_at': now(),
            'stage': 'empress_root_topology', 'stage_scope': 'Undated DTL reconciliation under every topology edge root; conditional event evidence, not demonstrated horizontal transfer or event dates',
            'output_dir': str(out), 'driver_sha256': sha(Path(__file__)),
            'adapter': file_spec(frozen_helper), 'python': EMPRESS_PYTHON,
            'empress_source': str(EMPRESS_SOURCE), 'empress_source_files': source_files,
            'inputs': list({x['path']: x for x in srcs}.values()), 'prepared_inputs': prepared,
            'host_tree': file_spec(host_out), 'mapping': file_spec(mapping_path),
            'n_tips': 25, 'costs_D_T_L': COSTS, 'models': models, 'states_expected': 1128,
            'cpu_workers': 1, 'thread_limit': 1, 'memory_limit_gib': 4,
            'state_timeout_seconds': 60, 'total_timeout_seconds': 1800,
            'primary_model': 'main_mfp', 'illustration_root_rule': 'First lexicographically sorted canonical split JSON; no biological outgroup claim',
            'median_rule': 'Minimize expected parent-specific event symmetric difference across all optimal reconstructions; exact integer score; lexicographic root and per-node full event tie breaks',
            'candidate_rule': 'Union of all unordered actual donor/recipient host incoming-branch pairs in the three primary-model root-rank-1 deterministic medians; retain every transfer instance and moving parasite clade; no fixed number or significance/appearance filter',
            'support_rule': 'Per model/root/cost, exact MPR presence is total minus count avoiding every matching event. Host pair and host pair + moving parasite clade are separate keys. Root/model proportions are descriptive sensitivity, not root posterior or biological probability.',
            'time_semantics': 'No dates from eMPRess; original A495 node/path map is preserved for separate local interval analysis. Pairwise overlap is not complete reconciliation chronological feasibility.',
            'zero_length_rule': 'Retain all topology edges, annotate exact zeros; do not collapse or silently drop duplicate topology scenarios',
            'association_basis': manifest['association_status']}
    save_json(run / PLAN_REL, plan)
    print(json.dumps({'status': 'frozen_no_real_reconciliations', 'plan': str(run / PLAN_REL),
                      'plan_sha256': sha(run / PLAN_REL), 'states': 1128, 'rooted_inputs': 376}, indent=2))


def deserialize_graph(rows):
    graph = {}
    for row in rows:
        node = tuple(row['mapping'])
        e = row['event']
        graph.setdefault(node, []).append((e[0], tuple(e[1]), tuple(e[2])))
    return graph


def transfer_info(node, event, pnodes, hosttips):
    if event[0] != 'T':
        raise ValueError('Expected transfer')
    donor = node[1]
    resident, moving = event[1:]
    if resident[1] != donor or moving[1] == donor:
        raise ValueError('Installed T tuple violates resident/landing convention')
    if set(hosttips[donor]) <= set(hosttips[moving[1]]) or set(hosttips[moving[1]]) <= set(hosttips[donor]):
        raise ValueError('Transfer host edges are ancestrally comparable')
    pair = tuple(sorted([donor, moving[1]]))
    return {'host_pair': pair, 'donor_host_branch': donor, 'recipient_host_branch': moving[1],
            'parasite_parent_node': node[0], 'moving_parasite_node': moving[0],
            'moving_parasite_clade_hash': pnodes[moving[0]]['clade_hash'],
            'moving_parasite_tips_json': pnodes[moving[0]]['tips_json'],
            'parent_event_tuple_json': canonical([node, event])}


def run_all(run, synthetic=False):
    plan_path = run / PLAN_REL
    plan = json.loads(plan_path.read_text())
    if sha(Path(__file__)) != plan['driver_sha256']:
        raise ValueError('Driver snapshot differs from frozen protocol')
    if plan['costs_D_T_L'] != [list(c) for c in COSTS]:
        raise ValueError('Unexpected frozen scope')
    calculated_states = sum(len(m['roots']) for m in plan['models']) * 3
    if plan['states_expected'] != calculated_states:
        raise ValueError('State count differs from models/root catalog')
    if not synthetic and (plan.get('data_kind') != 'real' or calculated_states != 1128
                          or len(plan['models']) != 8 or any(len(m['roots']) != 47 for m in plan['models'])
                          or plan['n_tips'] != 25):
        raise ValueError('Unexpected real scope')
    if synthetic and (plan.get('data_kind') != 'synthetic' or '08_QA' not in run.parts):
        raise ValueError('Synthetic execution requires an explicit fixture under 08_QA')
    for spec in plan['inputs'] + plan['prepared_inputs'] + plan['empress_source_files']:
        verified(spec)
    adapter = load_adapter(plan['adapter'])
    out = Path(plan['output_dir'])
    states_dir = out / 'states'
    states_dir.mkdir(exist_ok=False)
    cap = plan['memory_limit_gib'] * 1024 ** 3
    resource.setrlimit(resource.RLIMIT_AS, (cap, cap))
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    os.environ.update({k: '1' for k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS']})
    with (out / 'host_branch_map.tsv').open() as f:
        hrows = list(csv.DictReader(f, delimiter='\t'))
    hosttips = {x['host_node_id']: json.loads(x['selected_descendant_tips_json']) for x in hrows}
    mapping = {}
    for line in verified(plan['mapping']).read_text().splitlines():
        p, h = line.split(':')
        if p in mapping:
            raise ValueError('Duplicate exact association')
        mapping[p] = h
    if len(mapping) != plan['n_tips'] or len(set(mapping.values())) != plan['n_tips']:
        raise ValueError('Association not a 25-pair bijection')
    if synthetic and any(not x.startswith('SYN_') for x in list(mapping) + list(mapping.values())):
        raise ValueError('Synthetic tips must be explicitly synthetic')
    all_summaries, primary_transfers = [], []
    started = time.monotonic()

    def timeout(sig, frame):
        raise TimeoutError('Bounded eMPRess state exceeded time limit')

    signal.signal(signal.SIGALRM, timeout)
    for model in plan['models']:
        for root in model['roots']:
            pnodes = {x['node_id']: x for x in json.loads(verified(root['nodes']).read_text())}
            for costs in COSTS:
                if time.monotonic() - started > plan['total_timeout_seconds']:
                    raise TimeoutError('Batch exceeded total wall time')
                signal.alarm(plan['state_timeout_seconds'])
                state_start = time.monotonic()
                label = '%s__r%02d__D%dT%dL%d' % (model['model'], root['rank'], *costs)
                state_out = states_dir / label
                state_out.mkdir(exist_ok=False)
                inp = empress.ReconInputWrapper.from_files(str(verified(plan['host_tree'])),
                        str(verified(root['tree'])), str(verified(plan['mapping'])))
                # Topology/association DP only: never read input_reader host_distances.
                graph, optimum, count, graph_roots = recongraph_tools.DP(inp, *costs)
                stats = adapter.exact_graph_stats(graph, graph_roots, count, costs, optimum)
                c_events = {(node[0], node[1]): value for (node, event), value in stats['event_counts'].items() if event[0] == 'C'}
                if c_events != {(p, h): count for p, h in mapping.items()}:
                    raise ValueError('Contemporary events are not the exact 25 observed pairs at frequency 1')
                representative = adapter.deterministic_median(graph, graph_roots, stats)
                counts = {kind: representative['event_counts'].get(kind, 0) for kind in 'SDTLC'}
                if counts['C'] != plan['n_tips'] or sum(counts[k] * v for k, v in zip('DTL', costs)) != optimum:
                    raise ValueError('Representative count/cost check failed')
                full_rows = adapter.json_graph_rows(graph, stats)
                save_json(state_out / 'full_graph.json', {'mapping_roots': graph_roots,
                          'mpr_count': str(count), 'optimum': optimum, 'events': full_rows})
                rep_rows = adapter.json_graph_rows(representative['graph'], stats)
                save_json(state_out / 'representative_median.json', {'root': representative['root'],
                          'selection_rule': plan['median_rule'], 'integer_median_score': str(representative['integer_median_score']),
                          'event_counts': counts, 'events': rep_rows})
                event_rows = []
                rep_keys = {(n, e[0]) for n, e in representative['graph'].items()}
                for row in full_rows:
                    node = tuple(row['mapping'])
                    e = row['event']
                    is_rep = (node, e) in rep_keys
                    event_rows.append({'mapping_json': canonical(node), 'event_json': canonical(e),
                                       'type': e[0], 'is_representative': int(is_rep),
                                       'exact_numerator': row['mpr_presence_numerator'],
                                       'exact_denominator': row['mpr_count_denominator'],
                                       'conditional_frequency': row['conditional_mpr_frequency']})
                    if e[0] == 'T':
                        detail = transfer_info(node, e, pnodes, hosttips)
                        if is_rep and model['model'] == plan['primary_model'] and root['rank'] == 1:
                            primary_transfers.append({'state_id': label, 'cost_D': costs[0],
                               'cost_T': costs[1], 'cost_L': costs[2], **detail})
                save_tsv(state_out / 'events.tsv', event_rows, list(event_rows[0]))
                row = {'state_id': label, 'model': model['model'], 'root_rank': root['rank'],
                       'root_edge_hash': root['root_edge_hash'], 'D_cost': costs[0], 'T_cost': costs[1],
                       'L_cost': costs[2], 'optimum': optimum, 'mpr_count': str(count),
                       'graph_mapping_nodes': len(graph), 'graph_events': len(full_rows),
                       **{'representative_' + k: counts[k] for k in 'SDTLC'},
                       **{'expected_' + k: str(stats['expected_counts'][k]) for k in 'SDTLC'},
                       'full_graph_path': str(state_out / 'full_graph.json'),
                       'full_graph_sha256': sha(state_out / 'full_graph.json'),
                       'representative_path': str(state_out / 'representative_median.json'),
                       'representative_sha256': sha(state_out / 'representative_median.json'),
                       'wol_nodes_path': root['nodes']['path'], 'wol_nodes_sha256': root['nodes']['sha256'],
                       'seconds': time.monotonic() - state_start}
                save_json(state_out / 'complete.json', {'status': 'complete_undated_conditional_state', **row})
                all_summaries.append(row)
                signal.alarm(0)
                if len(all_summaries) % 47 == 0:
                    print(json.dumps({'completed_states': len(all_summaries), 'expected': plan['states_expected'],
                                      'elapsed_seconds': time.monotonic() - started, 'latest_state': label}), flush=True)
    if len(all_summaries) != plan['states_expected']:
        raise ValueError('Incomplete state count')
    remaining_seconds = plan['total_timeout_seconds'] - math.ceil(time.monotonic() - started)
    if remaining_seconds <= 0:
        raise TimeoutError('Batch exceeded total wall time before candidate aggregation')
    signal.alarm(remaining_seconds)
    save_tsv(out / 'state_summary.tsv', all_summaries, list(all_summaries[0]))
    pairs = sorted(set(tuple(x['host_pair']) for x in primary_transfers))
    candidates = []
    pair_ids = {}
    for pair in pairs:
        cid = 'candidate_' + digest(list(pair))
        pair_ids[pair] = cid
        candidates.append({'candidate_id': cid, 'host_branch_a': pair[0], 'host_branch_b': pair[1],
                           'host_a_tips_json': canonical(hosttips[pair[0]]),
                           'host_b_tips_json': canonical(hosttips[pair[1]]),
                           'candidate_basis': 'all_transfer_branch_pairs_from_three_primary_canonical_root_medians',
                           'direction_status': 'undirected_for_interpretation_actual_DTLevent_direction_retained_in_instances'})
    cfields = ['candidate_id', 'host_branch_a', 'host_branch_b', 'host_a_tips_json', 'host_b_tips_json', 'candidate_basis', 'direction_status']
    save_tsv(out / 'candidate_host_pairs.tsv', candidates, cfields)
    instances = []
    moving_keys = set()
    for instance in primary_transfers:
        pair = tuple(instance.pop('host_pair'))
        cid = pair_ids[pair]
        moving_keys.add((cid, instance['moving_parasite_clade_hash']))
        instances.append({'candidate_id': cid, **instance})
    ifields = ['candidate_id', 'state_id', 'cost_D', 'cost_T', 'cost_L', 'donor_host_branch', 'recipient_host_branch',
               'parasite_parent_node', 'moving_parasite_node', 'moving_parasite_clade_hash', 'moving_parasite_tips_json', 'parent_event_tuple_json']
    save_tsv(out / 'candidate_transfer_instances.tsv', instances, ifields)
    support_path = out / 'candidate_support.tsv'
    fields = ['candidate_id', 'evidence_key_type', 'moving_parasite_clade_hash', 'state_id', 'model',
              'root_rank', 'root_edge_hash', 'D_cost', 'T_cost', 'L_cost', 'exact_presence_numerator',
              'exact_mpr_denominator', 'conditional_mpr_frequency', 'has_any_qualifying_MPR']
    with support_path.open('x', newline='') as f:
        writer = csv.DictWriter(f, fields, delimiter='\t')
        writer.writeheader()
        for state in all_summaries:
            blob = json.loads(verified({'path': state['full_graph_path'], 'sha256': state['full_graph_sha256']}).read_text())
            graph = deserialize_graph(blob['events'])
            roots = [tuple(x) for x in blob['mapping_roots']]
            stats = adapter.exact_graph_stats(graph, roots, int(blob['mpr_count']))
            pnodes = {x['node_id']: x for x in json.loads(verified({'path': state['wol_nodes_path'], 'sha256': state['wol_nodes_sha256']}).read_text())}
            for candidate in candidates:
                cid = candidate['candidate_id']
                pair = tuple(sorted([candidate['host_branch_a'], candidate['host_branch_b']]))
                for moving in [None] + sorted(h for c, h in moving_keys if c == cid):
                    def predicate(node, event):
                        return (event[0] == 'T' and tuple(sorted([node[1], event[2][1]])) == pair
                                and (moving is None or pnodes[event[2][0]]['clade_hash'] == moving))
                    numerator = adapter.count_candidate_presence(graph, roots, stats, predicate)
                    writer.writerow({'candidate_id': cid, 'evidence_key_type': 'host_pair' if moving is None else 'host_pair_moving_clade',
                         'moving_parasite_clade_hash': moving or '', **{k: state[k] for k in ['state_id', 'model', 'root_rank', 'root_edge_hash', 'D_cost', 'T_cost', 'L_cost']},
                         'exact_presence_numerator': str(numerator), 'exact_mpr_denominator': str(stats['total']),
                         'conditional_mpr_frequency': numerator / stats['total'], 'has_any_qualifying_MPR': int(numerator > 0)})
    for spec in plan['inputs'] + plan['prepared_inputs'] + plan['empress_source_files'] + [plan['adapter']]:
        verified(spec)
    if sha(plan_path) != PLAN_SHA_AT_START:
        raise ValueError('Protocol changed during execution')
    summary = {'status': 'complete_undated_event_evidence_review_required', 'created_at': now(),
               'states_complete': len(all_summaries), 'models': len(plan['models']), 'roots_per_model': [len(m['roots']) for m in plan['models']],
               'costs_per_root': 3, 'primary_candidate_host_pairs': len(candidates),
               'primary_transfer_instances': len(instances), 'real_reconciliation_performed': not synthetic, 'data_kind': plan['data_kind'],
               'event_dating_performed': False, 'biological_outgroup_assumed': False,
               'cpu_workers': 1, 'thread_limit': 1, 'peak_rss_KiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
               'elapsed_seconds': time.monotonic() - started, 'protocol': file_spec(plan_path),
               'semantic_limit': plan['stage_scope'], 'candidate_rule': plan['candidate_rule'],
               'support_rule': plan['support_rule']}
    save_json(out / 'summary.json', summary)
    outputs = [file_spec(p) for p in sorted(out.rglob('*')) if p.is_file()]
    save_json(run / 'checkpoints/empress_root_topology.complete.json', {'status': 'complete',
              'stage': 'empress_root_topology', 'finished_at': now(), 'summary': summary, 'outputs': outputs})
    signal.alarm(0)
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == '__main__':
    a = argparse.ArgumentParser(description=__doc__)
    a.add_argument('--run-dir', required=True, type=Path)
    a.add_argument('--threads', type=int, default=1)
    a.add_argument('--stage', choices=['all'], default='all')
    a.add_argument('--freeze', action='store_true')
    a.add_argument('--synthetic', action='store_true')
    args = a.parse_args()
    if not 1 <= args.threads <= 4:
        a.error('threads must be 1..4; this worker always uses one')
    if Path(sys.prefix).absolute() != Path(EMPRESS_PYTHON).parent.parent.absolute():
        env = dict(os.environ)
        env.update({'PYTHONDONTWRITEBYTECODE': '1', 'MPLBACKEND': 'Agg', 'OMP_NUM_THREADS': '1',
                    'OPENBLAS_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1'})
        os.execve(EMPRESS_PYTHON, [EMPRESS_PYTHON, '-B', '-u', str(Path(__file__).resolve())] + sys.argv[1:], env)
    sys.path.insert(0, str(EMPRESS_SOURCE))
    from Bio import Phylo
    from Bio.Phylo.Newick import Tree, Clade
    import empress
    from empress.reconcile import recongraph_tools
    args.run_dir = args.run_dir.resolve()
    if args.freeze:
        if args.synthetic:
            a.error('--freeze cannot be combined with --synthetic')
        freeze(args.run_dir)
    else:
        PLAN_SHA_AT_START = sha(args.run_dir / PLAN_REL)
        # A direct reentry must not add a failure marker beside an old success.
        existing = [args.run_dir / 'checkpoints/empress_root_topology.complete.json',
                    args.run_dir / 'checkpoints/empress_root_topology.failed.json',
                    args.run_dir / OUT_REL / 'states']
        if any(path.exists() for path in existing):
            raise SystemExit('Existing stage output/checkpoint: refuse reentry without altering previous state')
        try:
            run_all(args.run_dir, args.synthetic)
        except BaseException as exc:
            failure = args.run_dir / 'checkpoints/empress_root_topology.failed.json'
            if not failure.exists():
                save_json(failure, {'status': 'failed', 'stage': 'empress_root_topology', 'finished_at': now(),
                          'error': str(exc), 'traceback': traceback.format_exc(), 'no_complete_checkpoint_written': True})
            raise
