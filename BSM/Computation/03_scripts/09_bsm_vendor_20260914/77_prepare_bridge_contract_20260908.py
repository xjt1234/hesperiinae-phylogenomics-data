#!/usr/bin/env python3
"""Freeze the tested bridge continuation recipe; create a new file, never a fit."""
import argparse
import hashlib
import json
from pathlib import Path

JOB = Path(__file__).resolve().parents[1]
PARENT = JOB / "02_config/conditional_completion_contract_20260908.json"
TARGET = JOB / "02_config/conditional_bridge_bsm_contract_20260908.json"

def sha(path):
    if path.is_symlink() or not path.is_file():
        raise ValueError(f"Not a regular evidence file: {path}")
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()

def main():
    args = argparse.ArgumentParser()
    args.add_argument("--selftest-summary", required=True)
    args.add_argument("--selftest-sha", required=True)
    opt = args.parse_args()
    if TARGET.exists():
        raise ValueError("Refuse to overwrite frozen bridge contract")
    if sha(PARENT) != "f3336caf6d5d398062c128d832f9a520a937cd7cf1995746dac9b7a19a9ff714":
        raise ValueError("Parent numerical contract changed")
    summary_path = Path(opt.selftest_summary)
    if not summary_path.is_absolute() or not summary_path.resolve().is_relative_to(JOB / "08_qa/bsm_failure_diagnosis_20260908"):
        raise ValueError("Selftest evidence outside exact diagnosis root")
    if sha(summary_path) != opt.selftest_sha:
        raise ValueError("Selftest evidence hash mismatch")
    summary = json.loads(summary_path.read_text())
    if summary.get("status") != "PASS" or summary.get("all_checks_pass") is not True or summary.get("whole_tree_maps") != 0:
        raise ValueError("Finite bridge selftest has not passed")
    helper = JOB / "03_scripts/74_ctmc_uniformization_bridge_20260908.R"
    test = JOB / "03_scripts/75_selftest_ctmc_bridge_20260908.R"
    if summary.get("helper_sha256") != sha(helper) or summary.get("selftest_script_sha256") != sha(test):
        raise ValueError("Selftest does not bind final helper and test bytes")
    contract = json.loads(PARENT.read_text())
    contract["purpose"] = "BSM-only repair: replace inefficient branch rejection with finite-tolerance endpoint-conditioned uniformization of the native simulator's effective generator; preserve all fitted inputs and endpoint-selection formulas."
    contract["user_authorization"] = "那你就继续排查吧 指导正常运行; earlier: 不需要再跑R5了，并且继续跑BSM"
    contract["r5_execution_mode"] = "NOT_PERFORMED_BY_USER_DECISION"
    contract.pop("r5", None)
    contract["resource_limits"]["unit"] = "bgb-r24-bridge-bsm-20260908.service"
    contract["bsm_path_sampler"] = "UNIFORMIZATION_NATIVE_EFFECTIVE_Q"
    contract["bridge_relative_tail_tolerance"] = 1e-12
    contract["bridge_expected_branch_calls"] = 1473
    contract["bridge_vector_cache_max_bytes"] = 2147483648
    contract["bridge_native_q_absolute_rounding_tolerance"] = 1e-7
    contract["bridge_native_q_relative_rounding_tolerance"] = 4 * 2**-23
    contract["bridge_generator_interpretation"] = "Qeff reproduces native exit rate -diag(Qraw) and normalized outgoing state probabilities; endpoint probabilities still use unchanged native Qraw/cache. These are not algebraically identical because Qraw has single-rounded row residuals. The diagnosed differences are documented, not hidden."
    contract["bridge_reference"] = "https://doi.org/10.1214/09-AOAS247"
    contract["parent_contract"] = {"path": str(PARENT.relative_to(JOB)), "sha256": sha(PARENT)}
    for rel in ("03_scripts/69_run_conditional_native_bsm_20260908.R", "03_scripts/70_supervise_conditional_bsm_20260908.py",
                "03_scripts/74_ctmc_uniformization_bridge_20260908.R", "03_scripts/75_selftest_ctmc_bridge_20260908.R",
                "03_scripts/76_run_conditional_bridge_bsm_20260908.R"):
        contract["frozen_helpers"][rel] = sha(JOB / rel)
    qa = "08_qa/bsm_failure_diagnosis_20260908/"
    evidence = {
        "finite_selftest": str(summary_path.relative_to(JOB)),
        "finite_selftest_checks": str((summary_path.parent / "checks.tsv").relative_to(JOB)),
        "branch_probability_checks": qa + "branch_numeric/branch_probability_summary.tsv",
        "generator_checks": qa + "branch_numeric/generator_checks.tsv",
        "structure_and_time_checks": qa + "native_code/time_state_audit.json",
        "native_generator_equivalence": qa + "native_code/NATIVE_EFFECTIVE_GENERATOR_EQUIVALENCE_20260908.md",
    }
    # A JSON object, deliberately not an array: R simplifyVector must retain lists.
    contract["bridge_required_evidence"] = {key: {"path": rel, "sha256": sha(JOB / rel)} for key, rel in evidence.items()}
    cache_rel = "05_bsm/M1_BSM_CONDITIONAL_20260908/stochastic_mapping_inputs.rds"
    identity_rel = "05_bsm/M1_BSM_CONDITIONAL_20260908/stochastic_mapping_inputs_identity.json"
    cache_sha = sha(JOB / cache_rel)
    if cache_sha != "4a54324e3fff499221a5b9613c32f8249a444710abb0f802659c84207c1c344f":
        raise ValueError("Native preparation cache changed")
    identity = json.loads((JOB / identity_rel).read_text())
    if identity["sha256"] != cache_sha or identity["technical_validation_sha256"] != "e107ade365a1523c8c888fcb591f0bb4ac88f46b1e3f6b9f56297929076af644":
        raise ValueError("Native preparation identity mismatch")
    contract["bridge_native_cache"] = {"path": cache_rel, "sha256": cache_sha,
        "identity_path": identity_rel, "identity_sha256": sha(JOB / identity_rel),
        "original_technical_validation_sha256": identity["technical_validation_sha256"]}
    contract["forced_history_policy"] = "Bridge kernel has no manual force-fit fallback. Retain rejection markers, original failed attempts and strict one-area/event-count gates. Do not import any old candidate map into the new namespace."
    contract["maxtries_interpretation"] = "40000 retained solely as a legacy provenance field; the branch bridge uses a maximum 10000 series terms, lambda<=1000 and a relative Poisson-tail bound, not a rejection loop."
    with TARGET.open("x") as handle:
        json.dump(contract, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    print(json.dumps({"status": "FROZEN_COMPUTATIONAL_RECIPE_NOT_SCIENTIFIC_ACCEPTANCE", "path": str(TARGET), "sha256": sha(TARGET), "scientific_acceptance": "NONE"}))

if __name__ == "__main__":
    main()
