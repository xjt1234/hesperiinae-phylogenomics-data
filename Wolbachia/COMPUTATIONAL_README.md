# Wolbachia analysis files

The primary host-association analysis uses the 25 matched samples in
`AHE_main25_associations.tsv`, the dated AHE host tree
`AHE_A495_main25.root_preserved.nwk`, and the bacterial tree in
`Wolbachia_tree/Tree/tree.treefile`. The 102-terminal screening table has a
different scope and is retained as `screening_status_102.tsv`.

## Files and analysis steps

| Step | Original source in `scripts/` | Inputs and results |
| --- | --- | --- |
| Host and bacterial distance inputs | `prepare_cophylogeny_main25_inputs_20260908.py`; `prepare_wol_common25_cophylogeny_20260908.py` | The root-level host tree, association table and crosswalk; `computational_inputs/AHE_A495_main25.full_to_pruned_node_map.tsv`; seven completed model-sensitivity trees under `computational_inputs/sensitivity_trees/`. |
| PACo and ParaFit | `run_cophylogeny_exact_inputs_R26.R`; `run_cophylogeny_main25_recovery_20260908.py` | Primary global results and link residuals are at the Wolbachia directory root. `all16_scenarios.tsv` includes the primary analysis, bacterial model sensitivities and explicitly labelled historical USCO-host comparisons. The original 16-scenario specification is under `computational_inputs/protocols/`. |
| Bootstrap sensitivity | `prepare_cophylogeny_bootstrap100_inputs_20260908.py`; `run_cophylogeny_bootstrap100_20260908.py` | Selection and analysis specifications are under `computational_inputs/protocols/`; `selected100_exact_source_lines.ufboot` contains the 100 selected concatenated-matrix bootstrap trees in the recorded order. Results are in the existing `all100_bootstrap_scenarios.tsv`. These are bootstrap trees, not per-gene trees. |
| DTL reconciliation and root sensitivity | `run_empress_root_topology_R26.py`; `empress_graph_adapter_R26.py` | The original protocol is under `computational_inputs/protocols/`. Existing `empress_runs/` tables record root edges, host branches, state summaries, transfer instances and conditional candidate support. |
| Conditional time windows and calibration comparison | `build_conditional_branch_windows_20260908.py`; `compare_candidate_windows_calibration_AB_20260908.py` | Existing `conditional_timewindows.tsv`, `candidate_endpoint_types.tsv` and `all24_candidate_windows_AB.tsv`; additional branch-age and full-tree clade correspondence tables are under `additional_results/conditional_timing/`. The A and B dated trees are supplied in the repository's `Divergence_time/` directory. |
| Screening and reference affinity | `rebuild_screening_table.py`; `summarize_supergroup_annotation_v2_20260908.py`; `review_main25_reference_affinity_scores_20260908.py` | Existing `screening_status_102.tsv`; additional 25-sample affinity categories, reference annotations, per-gene evidence and score tables under `additional_results/reference_affinity/`. |

## Numerical records

`additional_results/primary_cophylogeny/` contains the original host and bacterial
distance matrices, 9,999 primary bijective permutation draws and null statistics,
the corresponding PACo `r0` sensitivity records, slot keys, and the saved R
coordinates/Procrustes object. These files retain the precision of the original
exports; they have not been recalculated or rounded for this repository.

The primary PACo and ParaFit tests use the same bijective host-label permutations,
Cailliez correction and a fixed seed of 20260908 plus the explicit scenario index.
The supplementary `r0` analysis retains its separate null distribution. Bootstrap
tree selection uses Python seed 20260908; the exact selected source-tree indices
and per-scenario R seeds are recorded in the supplied specifications.

The Empress protocol covers eight bacterial models, all 47 candidate root edges
per model, and three `(D, T, L)` cost combinations, `(2, 2, 1)`, `(2, 3, 1)` and
`(2, 4, 1)`, giving 1,128 conditional states. The timing tables apply the recorded
host-branch correspondences to calibration scenarios A and B.

## Software and source layout

The computational scripts are unchanged copies of the scripts used for these
analyses. `SOURCE_FILES.tsv` records the archive member and SHA-256 of each copied
file. Their original directory names, absolute Linux paths, checkpoint checks,
frozen hashes and historical scenario identifiers remain in the source. The
`recovery` filename identifies the completed execution with the existing R
library; `USCO_legacy` identifies the historical host sensitivity rather than the
primary AHE analysis.

The archived execution used Python 3.10.19 and R 4.5.2. The R driver checks paco
0.4.2, ape 5.8-1, vegan 2.7-2, jsonlite 2.0.0 and digest 0.6.39. Python input
preparation and calibration comparison use Biopython. The Empress driver calls a
separately installed Empress source tree and uses Unix `resource`; the parallel R
wrappers use Unix `fcntl`. The adapter is included, while the Empress package
itself is an external dependency.

These files document the computations and provide the selected analysis inputs
and outputs. They are not a configured one-command workflow: execution on another
machine requires adapting paths and reconstructing the original checkpoint/input
layout while retaining the scientific settings and validating the file hashes.

## Scope of the archived evidence

The supplied reference-affinity tables preserve the main25 evidence and exported
score summaries. The complete raw `strict889_vs_all244.tsv` search table referenced
by the affinity scripts was not included in the source-data archive used to
assemble this repository supplement. The screening reconstruction also refers to
historical source tables and sample summary files beyond the selected screening
result supplied here. These scripts therefore provide the original computation
logic, not a claim that every upstream screening/search step can be rerun from
this selection alone.

Images, figure-generation scripts, inference-command wrappers, per-gene trees,
process logs and complete historical execution archives are outside this
supplement. The manuscript describes phylogenetic inference commands.
