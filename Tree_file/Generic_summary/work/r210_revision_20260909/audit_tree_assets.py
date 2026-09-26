"""Read-only genus topology and asset audit of the frozen 495-tip AHE NT tree."""
import csv
import hashlib
import json
import re
import subprocess
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from zipfile import ZipFile
from lxml import etree

ROOT = Path(r"D:\博士课题\AHE\SE_2026年9月1日返修")
OUT = ROOT / "work/r210_revision_20260909"
OUT.mkdir(parents=True, exist_ok=True)
PACKAGE = ROOT / "返修的分析/分歧时间/Hesperiinae_divergence_time_manuscript_package_20260904_145425"
TREE = PACKAGE / "02_INPUT_AND_CALIBRATIONS/verified_495_tip_input.treefile"
META = PACKAGE / "06_TREE_BEAUTIFICATION/tip_or_clade_colors.tsv"
FOCAL = "Pedesta Thoressa Onryza Creteus Astictopterus Caenides Paronymus Ancistroides Amblyscirtes Notamblyscirtes Polites Pompeius Appia Wallengrenia Polytremis Potanthus Taractrocera".split()


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def tsv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def save(name, rows):
    with (OUT / name).open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


@dataclass(eq=False)
class Node:
    name: str = ""
    length: float = 0
    children: list = field(default_factory=list)
    parent: object = None
    tips: set = field(default_factory=set)
    audit_id: str = ""


def parse(text):
    tokens = re.findall(r"[(),:;]|[^(),:;\s]+", text)
    pos = 0

    def sub():
        nonlocal pos
        node = Node()
        if tokens[pos] == "(":
            pos += 1
            while True:
                child = sub()
                child.parent = node
                node.children.append(child)
                node.tips.update(child.tips)
                marker = tokens[pos]
                pos += 1
                if marker == ")":
                    break
                assert marker == ","
            if tokens[pos] not in (":", ",", ")", ";"):
                node.name = tokens[pos]
                pos += 1
        else:
            node.name = tokens[pos]
            node.tips = {node.name}
            pos += 1
        if tokens[pos] == ":":
            pos += 1
            node.length = float(tokens[pos])
            pos += 1
        return node

    root = sub()
    nodes = []

    def visit(node):
        node.audit_id = f"N{len(nodes)+1:04d}"
        nodes.append(node)
        for child in node.children:
            visit(child)

    visit(root)
    return root, nodes


root, nodes = parse(TREE.read_text())
metadata = {r["tip_label"]: r for r in tsv(META)}
assert len(metadata) == 495 and set(metadata) == root.tips
assert sha(TREE) == "4d35d553ed5bbf08022daef2361ffa0c342f8694a53aaf08d5350692094439e5"
genus = {tip: metadata[tip]["biological_name"].split("_")[0] for tip in root.tips}
summary, components, memberships = [], [], []
for name in sorted(set(genus.values()) | set(FOCAL)):
    tips = {tip for tip, g in genus.items() if g == name}
    if not tips:
        summary.append({"genus": name, "focal": name in FOCAL, "sample_tips": 0, "named_binomials": 0,
                        "status": "NOT_SAMPLED_UNDER_THIS_NAME", "maximal_pure_components": 0,
                        "mrca_total_tips": 0, "mrca_support_raw": "", "mrca_other_genera": "",
                        "mrca_audit_node": "", "sampled_binomials": ""})
        continue
    mrca = min((n for n in nodes if tips <= n.tips), key=lambda n: len(n.tips))
    pure = [n for n in nodes if n.tips <= tips and (n.parent is None or not n.parent.tips <= tips)]
    binomials = {"_".join(metadata[t]["biological_name"].split("_")[:2]) for t in tips}
    mono = mrca.tips == tips
    status = "MONOPHYLETIC_IN_THIS_TREE" if mono else "NONMONOPHYLETIC_IN_THIS_TREE"
    if len(tips) == 1:
        status = "SINGLE_TIP_PLACEMENT_ONLY"
    elif len(binomials) == 1 and mono:
        status = "SINGLE_SPECIES_PLACEMENT_ONLY"
    summary.append({"genus": name, "focal": name in FOCAL, "sample_tips": len(tips), "named_binomials": len(binomials),
                    "status": status, "maximal_pure_components": len(pure), "mrca_total_tips": len(mrca.tips),
                    "mrca_support_raw": mrca.name if mrca.children else "",
                    "mrca_other_genera": ";".join(sorted({genus[t] for t in mrca.tips-tips})),
                    "mrca_audit_node": mrca.audit_id, "sampled_binomials": ";".join(sorted(binomials))})
    if name not in FOCAL:
        continue
    for idx, node in enumerate(pure, 1):
        sister = node.parent.tips-node.tips if node.parent else set()
        components.append({"genus": name, "component": idx, "audit_node": node.audit_id,
                           "sample_tips": len(node.tips), "support_raw": node.name if node.children else "",
                           "parent_support_raw": node.parent.name if node.parent else "",
                           "tips": ";".join(sorted(node.tips)), "sister_tips": ";".join(sorted(sister)),
                           "sister_genera": ";".join(sorted({genus[t] for t in sister}))})
        for tip in sorted(node.tips):
            memberships.append({"genus": name, "component": idx, "tip_label": tip,
                                "biological_name": metadata[tip]["biological_name"], "tribe": metadata[tip]["group"]})
save("all_genus_monophyly_AHE_NT.tsv", summary)
save("focal_genus_components_AHE_NT.tsv", components)
save("focal_tip_membership_AHE_NT.tsv", memberships)

def clade_hash(tips):
    return hashlib.sha256(("\n".join(sorted(tips)) + "\n").encode("utf-8")).hexdigest()


complex_definitions = {
    "Pedesta_Thoressa_Onryza": ["Pedesta", "Thoressa", "Onryza"],
    "Creteus_placement": ["Creteus"],
    "Astictopterus": ["Astictopterus"],
    "Caenides_Paronymus": ["Caenides", "Paronymus"],
    "Ancistroides": ["Ancistroides"],
    "Amblyscirtes_Notamblyscirtes": ["Amblyscirtes", "Notamblyscirtes"],
    "Polites_Pompeius_Appia_Wallengrenia": ["Polites", "Pompeius", "Appia", "Wallengrenia"],
    "Polytremis": ["Polytremis"],
    "Potanthus_Taractrocera": ["Potanthus", "Taractrocera"],
}
panels = []
node_rows = []
for complex_id, focal_genera in complex_definitions.items():
    target = {t for t in root.tips if genus[t] in focal_genera}
    mrca = min((n for n in nodes if target <= n.tips), key=lambda n: len(n.tips))
    panel_root = mrca.parent if complex_id == "Creteus_placement" else mrca

    def skeleton(node):
        focal_here = node.tips & target
        info = {"audit_id": node.audit_id, "descendant_tipset_sha256": clade_hash(node.tips),
                "support_raw": node.name if node.children else "", "branch_length": node.length,
                "sample_tips": len(node.tips), "target_samples": len(focal_here),
                "all_descendants": sorted(node.tips), "target_descendants": sorted(focal_here),
                "other_descendants": sorted(node.tips-target)}
        if not focal_here:
            info["display_role"] = "collapsible_background_clade"
            info["genera"] = sorted({genus[t] for t in node.tips})
        elif not node.children:
            info["display_role"] = "focal_tip"
            info["tip_label"] = node.name
            info["biological_name_as_recorded"] = metadata[node.name]["biological_name"]
            info["genus_as_labelled"] = genus[node.name]
        else:
            info["display_role"] = "retained_internal_node"
            info["children"] = [skeleton(child) for child in node.children]
        node_rows.append({"complex_id": complex_id, "audit_node": node.audit_id,
                          "descendant_tipset_sha256": info["descendant_tipset_sha256"],
                          "support_raw": info["support_raw"], "role": info["display_role"],
                          "all_descendants": ";".join(info["all_descendants"]),
                          "target_descendants": ";".join(info["target_descendants"]),
                          "other_descendants": ";".join(info["other_descendants"])})
        return info

    panels.append({"complex_id": complex_id, "target_genera_as_labelled": focal_genera,
                   "target_sample_count": len(target),
                   "distinct_label_binomials": len({"_".join(metadata[t]["biological_name"].split("_")[:2]) for t in target}),
                   "mrca_audit_node": mrca.audit_id, "mrca_tipset_sha256": clade_hash(mrca.tips),
                   "mrca_all_descendants": sorted(mrca.tips), "mrca_other_descendants": sorted(mrca.tips-target),
                   "context_rule": "parent_of_single_species_clade" if complex_id == "Creteus_placement" else "exact_focal_MRCA",
                   "tree": skeleton(panel_root)})
(OUT / "complex_panel_data_AHE_NT.json").write_text(json.dumps({"source_tree": str(TREE),
    "source_sha256": sha(TREE), "hash_convention": "SHA256 of sorted exact tip labels joined by LF with trailing LF, UTF-8",
    "taxon_policy": "Labels preserved. Historical/current names are not synonymized; simius labels must remain distinct samples.",
    "support_policy": "Raw paired numeric tokens preserved; manuscript convention is SH-aLRT/UFBoot, original ML log absent.",
    "panels": panels}, ensure_ascii=False, indent=2), encoding="utf-8")
save("complex_key_nodes_AHE_NT.tsv", node_rows)

files = subprocess.run(["rg", "--files", str(ROOT)], check=True, capture_output=True, text=True, encoding="utf-8").stdout.splitlines()
tree_paths = [Path(p) for p in files if Path(p).suffix.lower() in {".tre", ".nwk", ".treefile", ".contree", ".newick", ".nex"}]
groups = {}
for path in tree_paths:
    native = Path("\\\\?\\" + str(path))
    if native.stat().st_size > 1_000_000:
        continue
    text = native.read_text(encoding="utf-8", errors="replace")
    if "Creteus_cyrina" not in text or not ("Pedesta_" in text or "Thoressa_" in text):
        continue
    digest = sha(native)
    groups.setdefault(digest, {"sha256": digest, "representative_path": str(path), "copies": [], "text": text})["copies"].append(str(path))
assets = []
for group in groups.values():
    try:
        asset_root, asset_nodes = parse(group.pop("text"))
        labels = [n.name for n in asset_nodes if n.children and n.name]
        assets.append({**group, "tip_count": len(asset_root.tips), "internal_labels": len(labels),
                       "paired_numeric_support_labels": sum(bool(re.fullmatch(r"\d+(?:\.\d+)?/\d+(?:\.\d+)?", lab)) for lab in labels),
                       "exact_495_tip_set": asset_root.tips == root.tips})
    except Exception as exc:
        assets.append({**group, "parse_error": str(exc)})
(OUT / "tree_asset_inventory.json").write_text(json.dumps({"all_tree_files_scanned": len(tree_paths),
    "candidate_definition": "Files containing Creteus_cyrina and Pedesta or Thoressa; duplicate SHA grouped", "assets": assets}, ensure_ascii=False, indent=2), encoding="utf-8")

baseline = OUT / "baseline/manuscript.docx"
with ZipFile(baseline) as z:
    xml = etree.fromstring(z.read("word/document.xml"))
ns = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
quotes = []
for idx, p in enumerate(xml.xpath("/w:document/w:body/w:p", namespaces=ns)):
    text = "".join(p.xpath(".//w:t[not(ancestor::w:del)]/text()", namespaces=ns))
    if any(word.lower() in text.lower() for word in ["non-monophy", "Pedesta", "SH-aLRT", "UFBoot"]):
        quotes.append({"paragraph_index_zero_based": idx, "text": text})
(OUT / "audited_manuscript_quotes.json").write_text(json.dumps({"baseline_sha256": sha(baseline), "quotes": quotes}, ensure_ascii=False, indent=2), encoding="utf-8")
paired = [n.name for n in nodes if n.children and re.fullmatch(r"\d+(?:\.\d+)?/\d+(?:\.\d+)?", n.name)]
result = {"status": "PASS", "tree_sha256": sha(TREE), "sample_tips": len(root.tips), "metadata_exact_match": True,
          "paired_numeric_support_labels": len(paired), "internal_nodes": sum(bool(n.children) for n in nodes),
          "focal_summary": [r for r in summary if r["focal"]], "candidate_unique_trees": len(assets)}
(OUT / "audit_tree_assets_QA.json").write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
print(json.dumps(result, ensure_ascii=False, indent=2))
