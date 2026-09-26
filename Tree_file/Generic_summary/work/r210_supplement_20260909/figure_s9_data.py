"""Prepare auditable display trees; no estimation, rerooting or tip renaming."""
import copy
import csv
import hashlib
import json
import re
from pathlib import Path

ROOT = Path(r"D:\博士课题\AHE\SE_2026年9月1日返修")
SOURCE_WORK = ROOT / "work/r210_revision_20260909"
WORK = Path(__file__).resolve().parent
OUT = ROOT / "outputs/R2.10_附图S9_20260909/figures"
DATA_PATH = SOURCE_WORK / "complex_panel_data_AHE_NT.json"
DATA = json.loads(DATA_PATH.read_text(encoding="utf-8"))
PALETTE = ["#006C9E", "#B34814", "#7450A4", "#147B68"]
TITLES = [
    "Pedesta / Thoressa / Onryza", "Creteus", "Astictopterus",
    "Caenides / Paronymus", "Ancistroides", "Amblyscirtes / Notamblyscirtes",
    "Polites / Pompeius / Appia / Wallengrenia", "Polytremis",
    "Potanthus / Taractrocera",
]
SUBTITLES = [
    "Three non-exclusive genera", "Single-species placement", "Non-exclusive genus",
    "Two non-exclusive genera", "Three separate genus-labelled components",
    "Two simius samples; identity unresolved",
    "Mixed circumscription and sampling cases", "Three genus-labelled components",
    "Two exclusive sampled genus clades",
]
NOTES = [
    "Archived Onryza labels retained; type species absent.",
    "Two C. cyrina samples; genus monophyly not tested.",
    "A. jama and A. punctulata occur in separate positions.",
    "Paronymus MRCA: 75.2 / 77; interpret cautiously.",
    "Three components, not two exclusive clades.",
    "simius samples are not sisters.\nAmblyscirtes type species absent.",
    "Appia: one species. Wallengrenia: exclusive clade.",
    "Historical Polytremis labels retained.",
    "Boundary comparison; neither genus is non-exclusive.",
]


def tipset_hash(tips):
    return hashlib.sha256(("\n".join(sorted(tips)) + "\n").encode()).hexdigest()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


with (SOURCE_WORK / "figure7_genus_status_16_S4.tsv").open(encoding="utf-8-sig") as f:
    STATUS = list(csv.DictReader(f, delimiter="\t"))
TYPE_NAMES = {}
for row in STATUS:
    match = re.match(r"Yes \(([^)]+)\)", row["type_sampled_as_S4"])
    if match:
        TYPE_NAMES[row["genus_as_archived"]] = match.group(1).replace(" ", "_")


def type_tips(node, targets):
    return [t for t in node["all_descendants"]
            if any(t.startswith(n + "_") for g, n in TYPE_NAMES.items() if g in targets)]


def children(node):
    return node.get("children", [])


def leaves(node):
    if not children(node):
        return [node]
    return [leaf for ch in children(node) for leaf in leaves(ch)]


def nodes(node):
    return [node] + [n for ch in children(node) for n in nodes(ch)]


def prepare(panel, letter, detailed=False):
    targets = panel["target_genera_as_labelled"]
    palette = dict(zip(targets, PALETTE))
    root = copy.deepcopy(panel["tree"])

    def walk(n):
        genera = {t.split("_")[0] for t in n["all_descendants"]}
        n["type_representative_tips"] = type_tips(n, targets)
        n["focal_color"] = palette.get(next(iter(genera)), "#555C62") if len(genera) == 1 else "#555C62"
        if not detailed and len(genera) == 1:
            n.pop("children", None)
        for ch in children(n):
            walk(ch)
    walk(root)
    for i, n in enumerate(leaves(root), 1):
        n["unit_id"] = f"{letter}{i:02d}"
        tips = n["all_descendants"]
        genus = {t.split("_")[0] for t in tips}
        binomials = {" ".join(t.split("_")[:2]) for t in tips}
        n["is_focal"] = any(t.split("_")[0] in targets for t in tips)
        n["is_collapsed"] = len(tips) > 1
        if len(binomials) == 1:
            label = next(iter(binomials))
        elif len(genus) == 1:
            label = next(iter(genus))
        else:
            label = f"Other clade {n['unit_id']}"
        if len(tips) > 1:
            label += f" [{len(tips)}]"
        if n["type_representative_tips"]:
            label += " ★"
        if detailed and len(tips) == 1:
            label = tips[0] + (" ★" if n["type_representative_tips"] else "")
        n["display_label"] = label
        n["source_label"] = ""
        if not detailed and any("_simius_" in t for t in tips):
            n["display_label"] = tips[0].split("_")[0]
            n["source_label"] = "simius [Kawahara 2023]" if tips[0].endswith("kawahara2023") else "simius [NCBI]"
    return root


def export_membership(trees, filename):
    rows = []
    for panel, tree in zip(DATA["panels"], trees):
        for node in leaves(tree):
            rows.append({
                "panel": node["unit_id"][0], "complex_id": panel["complex_id"],
                "unit_id": node["unit_id"], "display_label": node["display_label"],
                "source_label": node["source_label"], "source_node_id": node["audit_id"],
                "support_raw": node["support_raw"], "collapsed": node["is_collapsed"],
                "sample_count": len(node["all_descendants"]),
                "tipset_sha256": node["descendant_tipset_sha256"],
                "type_representative_tips": ";".join(node["type_representative_tips"]),
                "all_descendant_tips": ";".join(node["all_descendants"]),
            })
    with (OUT / filename).open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys(), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    return rows


def validate(trees):
    checks = []
    for source, rendered in zip(DATA["panels"], trees):
        source_nodes = {n["audit_id"]: n for n in nodes(source["tree"])}
        terminal_tips = [t for n in leaves(rendered) for t in n["all_descendants"]]
        assert len(terminal_tips) == len(set(terminal_tips))
        assert set(terminal_tips) == set(source["tree"]["all_descendants"])
        for n in nodes(rendered):
            s = source_nodes[n["audit_id"]]
            assert tipset_hash(n["all_descendants"]) == n["descendant_tipset_sha256"]
            assert n["all_descendants"] == s["all_descendants"]
            assert n["support_raw"] == s["support_raw"]
            if children(n):
                assert [c["audit_id"] for c in children(n)] == [c["audit_id"] for c in children(s)]
        checks.append({"complex_id": source["complex_id"], "status": "PASS",
                       "source_tips": len(terminal_tips), "display_leaves": len(leaves(rendered)),
                       "nodes_checked": len(nodes(rendered))})
    assert sha(DATA["source_tree"]) == DATA["source_sha256"]
    return checks
