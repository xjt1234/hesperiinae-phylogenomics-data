#!/usr/bin/env python3
"""Freeze a NEW R4 M1 bounded bridge-BSM recipe; never fit or import old data."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

JOB = Path(__file__).resolve().parents[1]
TARGET = JOB / "02_config/conditional_bridge_bsm_R4_20260914.json"
FIT = "04_runs/M1_R4_final_v1/fit_result.rds"
POSTFIT = "04_runs/M1_R4_postfit_v1/postfit_recalculated_ancestral_states.rds"
FIT_SHA = "bb855fe8bfc2919d5eed94224b177463d517dc78fbe958cc34faf44fcca1bff4"
POSTFIT_SHA = "566b86c007d123be0fccf1064abf076ad26e6c082ceac8133cbed447ef52e08d"
INPUT_MANIFEST_SHA = "f7e121a3ce4aff304767eef6df39eba6f9ef473242aad48411a3b2c204d4a306"
VENDOR = "03_scripts/09_bsm_vendor_20260914"
QA = "08_qa/bsm_bridge_preflight_20260914"
FROZEN_INPUTS = {
    "01_inputs/frozen/tree_scenarioA417.tre": "06725d5a0c1ef29007aeb0f7657175c2d2329c75cb4d8fb2170e39b5fdaf6962",
    "01_inputs/frozen/geog_scenarioA417_analysis_order.LagrangePHYLIP": "3aa227e4600629eb77178ab68a30bb7cf408189bbcb9f8989ae1f184ce834ca4",
    "01_inputs/frozen/input_manifest.json": INPUT_MANIFEST_SHA,
    "02_config/timeperiods_5epochs.txt": "90dadc7755385b44d9333f22e2acaeac163bd0ce4507289e1720671eccc4e51f",
    "02_config/M1_conservative_dispersal_multipliers.txt": "34cc3ee7f0789d978948e35085fb9ed89de09f2600b50384bb074bd11db7029c",
    "02_config/area_order.tsv": "7daa151006fbe79d00ed92a83938aca4b04dbc89bbaa814710bf244dcb1a0804",
}

def sha(path: Path) -> str:
    if path.is_symlink() or not path.is_file() or path.resolve(strict=True) != path:
        raise ValueError(f"Required regular, non-redirected project file: {path}")
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()

def single_tsv(rel: str) -> dict:
    with (JOB / rel).open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 1:
        raise ValueError(f"Expected one-row table: {rel}")
    return rows[0]

def pin(rel: str) -> dict:
    return {"path": rel, "sha256": sha(JOB / rel)}

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--selftest-summary", default=QA + "/generic_v1/selftest_summary.json")
    args = parser.parse_args()
    if TARGET.exists():
        raise ValueError("Refuse to overwrite an existing frozen BSM contract")
    for rel, expected in {**FROZEN_INPUTS, FIT: FIT_SHA, POSTFIT: POSTFIT_SHA}.items():
        if sha(JOB / rel) != expected:
            raise ValueError(f"New R4 input pin mismatch: {rel}")
    summary_path = Path(args.selftest_summary)
    if not summary_path.is_absolute():
        summary_path = JOB / summary_path
    if not summary_path.is_relative_to(JOB / QA):
        raise ValueError("Generic selftest must be inside new R4 QA directory")
    summary = json.loads(summary_path.read_text())
    helper = f"{VENDOR}/74_ctmc_uniformization_bridge_20260908.R"
    test = f"{VENDOR}/selftest_generic_bridge_R4_20260914.R"
    expected_test = {"status": "PASS", "all_checks_pass": True, "n_checks": 25,
                     "whole_tree_maps": 0, "current_fit_tested": False,
                     "historical_scientific_fixtures_imported": False,
                     "scope": "GENERIC_SYNTHETIC_ALGORITHM_SELFTEST_ONLY"}
    if any(summary.get(k) != v for k, v in expected_test.items()):
        raise ValueError("Generic selftest missing or incorrectly represented")
    if summary.get("helper_sha256") != sha(JOB / helper) or summary.get("selftest_script_sha256") != sha(JOB / test):
        raise ValueError("Generic selftest does not bind final script/helper bytes")

    fit = single_tsv("04_runs/M1_R4_final_v1/fit_summary.tsv")
    postfit = single_tsv("04_runs/M1_R4_postfit_v1/postfit_QA.tsv")
    if any(fit.get(k) != v for k, v in {"model": "M1", "mode": "final", "max_range": "4",
            "k": "2", "n_tips": "417", "convcode": "0", "kkt1": "FALSE", "kkt2": "TRUE"}.items()):
        raise ValueError("Completed M1 fit recipe/diagnostics differ")
    if any(postfit.get(k) != v for k, v in {"status": "TECHNICAL_PASS_CONDITIONAL_ONLY",
            "scientific_acceptance": "NOT_GRANTED", "model": "M1", "max_range": "4",
            "n_tips": "417", "n_internal": "416", "n_states": "562", "warning_count": "0",
            "optimizer_fit_rds_sha256": FIT_SHA, "postfit_rds_sha256": POSTFIT_SHA}.items()):
        raise ValueError("Actual new R4 fixed-parameter postfit QA differs")
    if abs(float(postfit["delta_lnL"])) > 1e-6 or float(postfit["top_max_rowsum_error"]) > 1e-10:
        raise ValueError("Actual new postfit probability/likelihood gate failed")
    d, e, j, w = (float(fit[k]) for k in ("d", "e", "j", "w"))
    if not (d > 0 and e >= 0 and j == 0 and w == 1):
        raise ValueError("M1 parameters outside declared DEC recipe")

    helper_files = [str(p.relative_to(JOB)) for p in sorted((JOB / VENDOR).iterdir()) if p.is_file()]
    helper_files += ["03_scripts/10_run_R4_bridge_bsm_20260914.R",
                     "03_scripts/11_prepare_R4_bridge_contract_20260914.py",
                     "03_scripts/12_supervise_R4_BSM_20260914.py", "03_scripts/r44_env.sh"]
    frozen_helpers = {rel: sha(JOB / rel) for rel in helper_files}
    if frozen_helpers[f"{VENDOR}/05_run_bsm_m1_stratified.R"] != "cf03d8019906827ca3c6fcb49874d210acd4f6983682956c684d7c52a100abc5":
        raise ValueError("Original counting/stability helper bytes changed")
    if frozen_helpers[f"{VENDOR}/70_supervise_conditional_bsm_20260908.py"] != "39a8bc0d3a96740453fac8956e088b3589aa613a170091f00606806514b94d76":
        raise ValueError("Original transaction supervisor bytes changed")
    note = f"{QA}/NUMERICAL_AND_SCOPE_NOTE.md"
    evidence = {
        "generic_selftest": pin(str(summary_path.relative_to(JOB))),
        "generic_selftest_checks": pin(str((summary_path.parent / "checks.tsv").relative_to(JOB))),
        "actual_new_R4_postfit_QA": pin("04_runs/M1_R4_postfit_v1/postfit_QA.tsv"),
        "numerical_scope_note": pin(note),
    }
    contract = dict(
        schema_version="1.0", status="APPROVED_CONDITIONAL_CONTINUATION",
        created_date="2026-09-14", created=datetime.now(timezone.utc).isoformat(),
        user_authorization="R4按照计划继续分析",
        purpose="Complete the finite preplanned M1 BSM under the NEW R4 geography exception scenario; no refitting, no old caches or histories.",
        backend="native", scientific_acceptance="NONE", model="M1", KKT1=False,
        fit_rds=FIT, fit_sha256=FIT_SHA, postfit_rds=POSTFIT, postfit_sha256=POSTFIT_SHA,
        d=d, e=e, j=j, w=w, n_tips=417, n_states=562, max_range_size=4, n_time_strata=5,
        input_manifest_sha256=INPUT_MANIFEST_SHA,
        area_order=["AF","AUS","CAM","ENA","EPA","IND","MDG","ORI","SAM","WNA","WPA"],
        old_world_areas=["AF","AUS","EPA","IND","MDG","ORI","WPA"],
        new_world_areas=["CAM","ENA","SAM","WNA"],
        null_range_aggregation="Separate NULL category; never assign it silently to either world.",
        geography_scenario="All corrected geographic rows except explicitly authorized Pelopidas_mathias legacy coding; joint coding/range comparison to R6, not pure range-cap sensitivity.",
        resource_limits=dict(unit="bgb-pelopidas-r4-bsm-20260914.service", cpu_quota_cores=4,
            tasks_max=32, memory_max_bytes=16000000000, memory_swap_max_bytes=0,
            bsm_bgb_cores=1, native_math_threads=1, overall_walltime_seconds=172800),
        r5_execution_mode="NOT_PERFORMED_BY_USER_DECISION",
        minimum_successful_maps=100, first_stability_checkpoint=200, max_successful_maps=500,
        map_batch_size=100, max_attempts=1000, maxtries=40000, seed_base=202609140,
        pilot_successful_maps=1, no_success_failure_stop=10,
        bsm_stability_policy="Unchanged vendored05 rules: 5% relative MCSE (absolute0.05 below mean1), cumulative d/e change and identical top-three route/extinction sets vs previous and initial100; first passing checkpoint>=200, cap500.",
        bsm_required_qa=["actual_new_fit_and_postfit_hashes","fresh_native_preparation_only",
            "562_states_five_epochs","1473_audited_branch_calls","416_unique_cladogenetic_nodes",
            "one_area_gains_and_losses","j_zero","independent_vs_native_counts",
            "total_area_epoch_conservation","zero_warnings","no_forced_histories",
            "separate_map_and_source_seeds","all_attempts_and_orphans_retained","cache_identity_before_resume"],
        preparation_policy="FRESH_NATIVE_PREPARATION_FROM_NEW_R4_POSTFIT_ONLY",
        bsm_path_sampler="UNIFORMIZATION_NATIVE_EFFECTIVE_Q",
        bridge_relative_tail_tolerance=1e-12, bridge_expected_branch_calls=1473,
        bridge_vector_cache_max_bytes=2147483648,
        bridge_native_q_absolute_rounding_tolerance=1e-7,
        bridge_native_q_relative_rounding_tolerance=4 * 2**-23,
        bridge_generator_interpretation="Qeff matches native waiting rates and normalized outgoing probabilities; freshly recomputed native endpoint selection uses unchanged Qraw. Single-rounded residuals are bounded and disclosed, not claimed algebraically identical.",
        bridge_reference="https://doi.org/10.1214/09-AOAS247",
        bridge_required_evidence=evidence,
        maxtries_interpretation="40000 is a legacy provenance field only; bridge has max10000series terms,lambda<=1000,relative Poisson-tail control and no rejection/manual fallback.",
        forced_history_policy="No manual fallback in replacement kernel. Reject any warning, force-fit marker, invalid event, or failed conservation; retain all attempts and rejected candidates.",
        completion_policy="COMPUTE_COMPLETE_CONDITIONAL_REVIEW_REQUIRED, INCOMPLETE, or FAILED only; never silently promote KKT1 or claim unconditional scientific acceptance.",
        no_automatic_tasks=["more_R4_optimization","KKT_fixed_point_probes","R5","M0_or_M2_BSM",
            "R6_stop_restart_or_edit","Scenario_B","posterior_tree_refits","six_model_revival"],
        frozen_inputs=FROZEN_INPUTS, frozen_helpers=frozen_helpers,
        numerical_audit={**pin(note),"interpretation":"New-fit conditional scope only; old geography-specific numerical claims are not inherited."},
        code_reuse_policy="Vendored source only; old scientific cache, raw histories, fitted-Q fixtures and old-figure data are not runtime inputs.",
        table_S2_policy="Three DEC variants within this R4 scenario only; do not pool R4/R6 AIC weights.",
    )
    with TARGET.open("x") as handle:
        json.dump(contract, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    print(json.dumps({"status":"FROZEN_COMPUTATIONAL_RECIPE_NOT_SCIENTIFIC_ACCEPTANCE",
        "path":str(TARGET),"sha256":sha(TARGET),"scientific_acceptance":"NONE"}))

if __name__ == "__main__":
    main()
