#!/usr/bin/env python3
"""Exact post-processing of an unchanged eMPRess optimal reconciliation DAG.

No file I/O or reconciliation is performed on import. Events are identified by
(parent mapping, complete event tuple), including distinct contemporary tips.
"""
from collections import Counter, defaultdict
from fractions import Fraction
import json

SINK = (None, None)


def canonical(value):
    return json.dumps(value, ensure_ascii=True, separators=(",", ":"))


def graph_order(graph, roots):
    """Return reachable mapping nodes in children-before-parents order."""
    if not roots or len(set(roots)) != len(roots):
        raise ValueError("Missing or duplicate optimal mapping roots")
    state, order = {}, []

    def visit(node):
        if node == SINK:
            return
        if state.get(node) == 1:
            raise ValueError("Cycle in reconciliation graph")
        if state.get(node) == 2:
            return
        if node not in graph or not graph[node]:
            raise ValueError("Missing node or empty event list")
        state[node] = 1
        if len(set(graph[node])) != len(graph[node]):
            raise ValueError("Duplicate complete event tuple under same mapping")
        for event in graph[node]:
            if len(event) != 3 or event[0] not in "SDTLC":
                raise ValueError("Malformed event")
            children = event[1:]
            if event[0] == "C" and children != (SINK, SINK):
                raise ValueError("Contemporary event has descendants")
            if event[0] == "L" and (children[0] == SINK or children[1] != SINK):
                raise ValueError("Loss must have one real first child")
            if event[0] in "SDT" and (SINK in children or children[0] == children[1]):
                raise ValueError("Binary event needs two distinct real children")
            for child in children:
                visit(child)
        state[node] = 2
        order.append(node)

    for root in sorted(roots, key=canonical):
        visit(root)
    if set(order) != set(graph):
        raise ValueError("Graph contains unreachable mappings")
    return order


def exact_graph_stats(graph, roots, expected_n=None, costs=None, optimum=None):
    """Exact inside/outside counts and conditional presence frequencies.

    This assumes the eMPRess tree reconciliation DAG (a mapping cannot recur
    twice within one valid reconciliation), not an arbitrary shared-node DAG.
    """
    order = graph_order(graph, roots)
    inside = {SINK: 1}
    for node in order:
        inside[node] = sum(inside[e[1]] * inside[e[2]] for e in graph[node])
    total = sum(inside[root] for root in roots)
    if total <= 0 or (expected_n is not None and total != expected_n):
        raise ValueError("Recomputed MPR count disagrees with installed DP")
    outside = defaultdict(int)
    for root in roots:
        outside[root] += 1
    event_counts = {}
    for node in reversed(order):
        for event in graph[node]:
            a, b = event[1:]
            event_counts[node, event] = outside[node] * inside[a] * inside[b]
            if a != SINK:
                outside[a] += outside[node] * inside[b]
            if b != SINK:
                outside[b] += outside[node] * inside[a]
    node_counts = {node: outside[node] * inside[node] for node in order}
    if any(c < 0 or c > total for c in event_counts.values()):
        raise ValueError("Invalid event presence count; not a tree-reconciliation DAG")
    if any(c < 0 or c > total for c in node_counts.values()):
        raise ValueError("Invalid mapping presence count")
    expected_counts = {kind: sum((Fraction(c, total) for (n, e), c in event_counts.items()
                                  if e[0] == kind), Fraction(0)) for kind in "SDTLC"}
    if costs is not None and optimum is not None:
        expected_cost = sum(expected_counts[k] * v for k, v in zip("DTL", costs))
        if expected_cost != optimum:
            raise ValueError("Expected event cost differs from optimum")
    return {"order": order, "inside": inside, "outside": dict(outside),
            "total": total, "event_counts": event_counts, "node_counts": node_counts,
            "expected_counts": expected_counts}


def deterministic_median(graph, roots, stats):
    """One event-distance median; exact integer score, deterministic tie rule.

    Maximize sum(2 * event_presence_count - total_MPR_count). This minimizes
    mean symmetric event-set difference from all MPRs. First select the
    lexicographically smallest full root mapping among tied optimum roots,
    then the lexicographically smallest complete event at each tied node.
    """
    total = stats["total"]
    scores, chosen = {SINK: 0}, {}
    for node in stats["order"]:
        pairs = [(2 * stats["event_counts"][node, event] - total
                  + scores[event[1]] + scores[event[2]], event) for event in graph[node]]
        best = max(x[0] for x in pairs)
        scores[node] = best
        chosen[node] = min((e for score, e in pairs if score == best), key=canonical)
    best = max(scores[root] for root in roots)
    root = min((r for r in roots if scores[r] == best), key=canonical)
    representative = {}

    def walk(node):
        if node == SINK:
            return
        if node in representative:
            raise ValueError("Representative repeats a mapping node")
        event = chosen[node]
        representative[node] = [event]
        walk(event[1])
        walk(event[2])

    walk(root)
    return {"root": root, "graph": representative, "integer_median_score": best,
            "event_counts": dict(Counter(e[0][0] for e in representative.values()))}


def count_candidate_presence(graph, roots, stats, predicate):
    """Number of MPRs with >=1 matching event, never a sum of frequencies."""
    avoiding = {SINK: 1}
    for node in stats["order"]:
        avoiding[node] = sum(avoiding[e[1]] * avoiding[e[2]] for e in graph[node]
                             if not predicate(node, e))
    return stats["total"] - sum(avoiding[root] for root in roots)


def json_graph_rows(graph, stats):
    """Lossless tuples with exact integer counts serialized as decimal strings."""
    rows = []
    for node in sorted(graph, key=canonical):
        for event in sorted(graph[node], key=canonical):
            count = stats["event_counts"][node, event]
            rows.append({"mapping": node, "event": event,
                         "mpr_presence_numerator": str(count),
                         "mpr_count_denominator": str(stats["total"]),
                         "conditional_mpr_frequency": float(Fraction(count, stats["total"]))})
    return rows
