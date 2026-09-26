"""Read-only Figure 7 genus, rooting, and sample provenance audit; no tree fitting."""
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

import openpyxl

ROOT = Path(r"D:\博士课题\AHE\SE_2026年9月1日返修")
OUT = ROOT / "work/r210_revision_20260909"
PARSER = OUT / "audit_tree_assets.py"
# Reuse only declarations, not the source script's audit/export/inventory body.
library = {}
exec(compile(PARSER.read_text(encoding="utf-8").split("root, nodes = parse")[0], str(PARSER), "exec"), library)
TREE, META = library["TREE"], library["META"]
root, nodes = library["parse"](TREE.read_text())
metadata = {r["tip_label"]: r for r in library["tsv"](META)}
summary = {r["genus"]: r for r in library["tsv"](OUT / "all_genus_monophyly_AHE_NT.tsv")}
WORKBOOK = ROOT / "outputs/r25_temperature/Hesperiinae_Supplementary_Tables_S1-S7_R2.5_温度模型修订.xlsx"
workbook = openpyxl.load_workbook(WORKBOOK, read_only=True, data_only=True)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def tiphash(tips):
    return hashlib.sha256(("\n".join(sorted(tips)) + "\n").encode("utf-8")).hexdigest()


def dump(name, value):
    (OUT / name).write_text(json.dumps(value, ensure_ascii=False, indent=2), encoding="utf-8")


def save(name, rows):
    with (OUT / name).open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def node_info(node):
    return {"audit_node": node.audit_id, "descendant_tipset_sha256": tiphash(node.tips),
            "support_raw": node.name if node.children else "", "sample_count": len(node.tips),
            "descendants": sorted(node.tips),
            "tribal_groups_from_metadata": dict(Counter(metadata[t]["group"] for t in node.tips))}


notes = {
    "Creteus": "Two samples of C. cyrina form an exclusive clade; this is single-species placement and does not test genus monophyly.",
    "Appia": "Only A. appia sampled; placement within the Pompeius-labelled complex does not test Appia monophyly.",
    "Potanthus": "Four sampled species form an exclusive clade; treat as a boundary-comparison case, not recovered non-monophyly.",
    "Taractrocera": "Two sampled species form an exclusive clade; treat as a boundary-comparison case, not recovered non-monophyly.",
    "Wallengrenia": "Three species represented by four samples form an exclusive clade within the broader Polites-labelled complex.",
    "Paronymus": "P. xanthias is sister to Ceratricula semilutea (100/100); their parent with P. ligora has only 75.2/77. Preserve this weaker support.",
    "Ancistroides": "Three maximal pure components, not two exclusive clades: gemmifer and nigrita successively precede the longicornis/Notocrypta/Udaspes complex.",
    "Amblyscirtes": "Five uncontested Amblyscirtes-labelled species form an exclusive clade (100/100). A. simius is separate and also represented under Notamblyscirtes by another non-sister sample; source/voucher validation is required before genus-level interpretation.",
    "Onryza": "Type species not sampled according to Table S4; historical labels are retained without a new taxonomic assignment.",
    "Notamblyscirtes": "Single sampled species; its simius sample and the separate Amblyscirtes simius sample require source/voucher validation. Associated context row, not a Table S4 genus row.",
}
s4 = {}
for rownum, row in enumerate(workbook["Table S4"].iter_rows(values_only=True), 1):
    if row[0] in summary:
        s4[row[0]] = {"row": rownum, "type_species": row[1], "sampled": row[2], "old_claim": row[3]}
assert len(s4) == 16
rows = []
for name in list(s4) + ["Notamblyscirtes"]:
    current = summary[name]
    raw = current["status"]
    if name == "Amblyscirtes":
        status = "historical-label-discrepancy"
    elif raw.startswith("SINGLE_"):
        status = "single-species not-testable"
    elif raw == "MONOPHYLETIC_IN_THIS_TREE":
        status = "multi-species exclusive clade"
    else:
        status = "nonexclusive"
    focal = {t for t, r in metadata.items() if r["biological_name"].split("_")[0] == name}
    mrca = min((n for n in nodes if focal <= n.tips), key=lambda n: len(n.tips))
    record = s4.get(name, {})
    rows.append({"genus_as_archived": name, "status": status, "sample_tips": len(focal),
                 "distinct_binomials_under_archived_genus_label": int(current["named_binomials"]),
                 "maximal_pure_components": current["maximal_pure_components"],
                 "mrca_sample_count": len(mrca.tips), "mrca_support_raw": current["mrca_support_raw"],
                 "mrca_audit_node": mrca.audit_id, "mrca_tipset_sha256": tiphash(mrca.tips),
                 "type_species_as_S4": record.get("type_species", "not a Table S4 row"),
                 "type_sampled_as_S4": record.get("sampled", "not assessed here"),
                 "S4_excel_row": record.get("row", ""),
                 "interpretation": notes.get(name, "Sampled archived genus labels do not form an exclusive clade in this tree; this is not a formal topology test or a taxonomic recircumscription."),
                 "support_limit": "Raw pairs follow manuscript SH-aLRT/UFBoot convention. MRCA support alone is not a test of genus monophyly; original ML inference log is absent.",
                 "tips": ";".join(sorted(focal))})
save("figure7_genus_status_16_S4.tsv", rows[:16])
save("figure7_genus_status_17_with_context.tsv", rows)

# Preserve the existing Newick root; do not use a parser's rooted flag as evidence.
ingroup = {t for t, m in metadata.items() if m["group"] != "Outgroup"}
outgroup = set(metadata) - ingroup
assert len(root.children) == 2
assert {frozenset(c.tips) for c in root.children} == {frozenset(ingroup), frozenset(outgroup)}
rooting = {"source_tree": str(TREE), "source_sha256": sha(TREE), "root_degree": len(root.children),
           "rooting_operation_performed": False, "root_support_raw": root.name,
           "outgroup_count": len(outgroup), "ingroup_count": len(ingroup),
           "root_children": [node_info(c) for c in root.children],
           "finding": "Actual root bipartition is 7 Papilionidae outgroup samples versus 488 Hesperiinae samples, not 5 versus 490.",
           "outgroup_tip_labels": sorted(outgroup)}
dump("rooting_and_outgroup_audit.json", rooting)

# Source metadata retains legacy terminal labels; match only the explicit two records.
headers = None
source_rows = {}
legacy_targets = {"Notamblyscirtes_simius_mydata": "Notamblyscirtes_simius_ncbi",
                  "Amblyscirtes_simius_NEE": "Amblyscirtes_simius_kawahara2023"}
for rownum, row in enumerate(workbook["Table S1"].iter_rows(values_only=True), 1):
    if row[0] == "Original terminal label":
        headers = list(row)
    if row[0] in legacy_targets:
        source_rows[legacy_targets[row[0]]] = {"excel_row": rownum, "record": dict(zip(headers, row))}
assert len(source_rows) == 2
details = []
for tip, record in source_rows.items():
    leaf = next(n for n in nodes if not n.children and n.name == tip)
    lineage = []
    ancestor = leaf.parent
    while ancestor is not None:
        lineage.append(node_info(ancestor))
        ancestor = ancestor.parent
    details.append({"current_tree_tip": tip, "current_metadata": metadata[tip], "table_S1": record,
                    "source_correspondence_basis": "Exact biological name and provenance category; S1 intentionally retains legacy suffixes. No input name was changed.",
                    "immediate_sister_tips": sorted(leaf.parent.tips - leaf.tips),
                    "immediate_parent": node_info(leaf.parent), "ancestor_path_to_root": lineage})
five = {t for t in metadata if t.startswith("Amblyscirtes_") and "simius" not in t}
five_root = min((n for n in nodes if five <= n.tips), key=lambda n: len(n.tips))
assert five_root.tips == five and len(five) == 5
pair = set(source_rows)
pair_root = min((n for n in nodes if pair <= n.tips), key=lambda n: len(n.tips))
dump("simius_label_provenance_audit.json", {
    "status": "label/provenance discrepancy requiring voucher/source validation",
    "source_workbook": str(WORKBOOK), "source_workbook_sha256": sha(WORKBOOK),
    "nominal_species_count_for_simius": 1, "separate_sample_records": 2,
    "policy": "Two archive labels refer to the same nominal species but have distinct source/voucher records and separated positions. Do not infer contamination, relabel, merge terminals, or count them as two nominal species.",
    "simius_samples": details, "two_simius_MRCA": node_info(pair_root),
    "five_other_Amblyscirtes_clade": node_info(five_root),
    "five_other_Amblyscirtes_sister_tips": sorted(five_root.parent.tips-five_root.tips),
    "five_other_Amblyscirtes_parent": node_info(five_root.parent),
    "Table_S4_Amblyscirtes_original_claim": s4["Amblyscirtes"],
    "limit": "Local source metadata audited. Voucher identity, sequence provenance, and current nomenclature were not independently revalidated against original repositories in this script."})

panel_path = OUT / "complex_panel_data_AHE_NT.json"
panel_data = json.loads(panel_path.read_text(encoding="utf-8"))
for panel in panel_data["panels"]:
    if panel["complex_id"] == "Amblyscirtes_Notamblyscirtes":
        panel["distinct_nominal_species_count_for_reporting"] = 6
        panel["nominal_count_note"] = "Seven samples represent six nominal species: Amblyscirtes simius and Notamblyscirtes simius are two sample labels for the same nominal species. Samples are not merged in the topology."
    else:
        panel["distinct_nominal_species_count_for_reporting"] = panel["distinct_label_binomials"]
        panel["nominal_count_note"] = "First two biological-name components, under archived labels; no general nomenclatural revision is implied."
dump("complex_panel_data_AHE_NT.json", panel_data)
qa = {"status": "PASS", "source_tree_sha256": sha(TREE), "source_metadata_sha256": sha(META),
      "source_workbook_sha256": sha(WORKBOOK), "source_parser_sha256": sha(PARSER),
      "status_rows_S4": 16, "status_rows_with_context": 17,
      "16_row_status_counts": dict(Counter(r["status"] for r in rows[:16])),
      "root_bipartition": [len(outgroup), len(ingroup)], "simius_sample_records": len(details),
      "two_simius_nominal_species": 1, "five_uncontested_Amblyscirtes_exclusive": True,
      "inputs_modified": False, "new_phylogeny_fitted": False,
      "tool_failure_resolved": "Bio.Phylo import unavailable; used the existing custom Newick parser without installing dependencies."}
dump("figure7_status_and_provenance_QA.json", qa)
print(json.dumps(qa, ensure_ascii=False, indent=2))
