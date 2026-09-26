#!/usr/bin/env python3
"""Bounded fresh R4 BSM, followed by separately hash-pinned result processing.

The immutable 70 supervisor supplies the transactional finite-stage machinery.
This adapter rebinds every scientific path and limits resources alongside R6.
Default is a non-mutating preflight, not a launch or a restart.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys
from datetime import datetime, timezone

SELF = Path(__file__).resolve()
JOB = SELF.parents[1]
BASE = JOB / "03_scripts/09_bsm_vendor_20260914/70_supervise_conditional_bsm_20260908.py"
BASE_SHA = "39a8bc0d3a96740453fac8956e088b3589aa613a170091f00606806514b94d76"
CONTRACT = JOB / "02_config/conditional_bridge_bsm_R4_20260914.json"
RUNNER = JOB / "03_scripts/10_run_R4_bridge_bsm_20260914.R"
WRAPPER = JOB / "03_scripts/r44_env.sh"
OUT = JOB / "08_qa/conditional_bridge_R4_20260914"
BSM = JOB / "05_bsm/M1_R4_BSM_BRIDGE_20260914"
FINALIZER = JOB / "03_scripts/33_finalize_R4_BSM_20260914.py"
UNIT = "bgb-pelopidas-r4-bsm-20260914.service"
R6_UNITS = (
    "bgb-geog417-r6-live-three-models-20260910.service",
    "bgb-geog417-r6-queue-pathguard-v2-20260910.service",
)


def digest(path):
    if path.is_symlink() or not path.is_file() or path.resolve() != path:
        raise ValueError(f"Expected unredirected regular file: {path}")
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def check_hash(path, expected):
    if not re.fullmatch("[a-f0-9]{64}", expected) or digest(path) != expected:
        raise ValueError(f"SHA256 mismatch: {path}")


def now():
    return datetime.now(timezone.utc).isoformat()


def r6_reservations():
    """Read only. Never signal or reconfigure either R6 controller."""
    total_cpu = 0.0
    total_memory = 0
    observations = []
    for unit in R6_UNITS:
        cp = subprocess.run(
            ["systemctl", "--user", "show", unit, "-p", "ActiveState", "-p", "ControlGroup"],
            check=True, text=True, capture_output=True, timeout=15,
        )
        info = dict(line.split("=", 1) for line in cp.stdout.splitlines() if "=" in line)
        observation = {"unit": unit, **info}
        if info.get("ActiveState") in ("active", "activating", "deactivating"):
            group = info.get("ControlGroup", "")
            if not group or ".." in Path(group).parts:
                raise ValueError("Cannot establish active R6 resource ownership")
            cg = Path("/sys/fs/cgroup") / group.lstrip("/")
            quota, period = (cg / "cpu.max").read_text().split()
            mem = (cg / "memory.max").read_text().strip()
            if quota == "max" or mem == "max":
                raise ValueError("R6 resource reservation is unexpectedly uncapped")
            cpu = int(quota) / int(period)
            total_cpu += cpu
            total_memory += int(mem)
            observation.update(cpu_quota_cores=cpu, memory_max_bytes=int(mem))
        observations.append(observation)
    if total_cpu > 40 or total_memory > 180000000000:
        raise ValueError("R6 reservations exceed the preserved 40-core / 180-GB envelope")
    return {"observed_utc": now(), "units": observations,
            "r6_cpu_quota_cores": total_cpu, "r6_memory_max_bytes": total_memory,
            "combined_with_bsm_cpu_quota_cores": total_cpu + 4,
            "combined_with_bsm_memory_max_bytes": total_memory + 16000000000,
            "r6_mutations": []}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner-sha", required=True)
    parser.add_argument("--contract-sha", required=True)
    parser.add_argument("--finalizer-sha", required=True)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--run", action="store_true")
    mode.add_argument("--selftest", action="store_true")
    args = parser.parse_args()
    check_hash(BASE, BASE_SHA)
    spec = importlib.util.spec_from_file_location("frozen_bsm_supervisor_70", BASE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    for name, value in {
        "__file__": str(SELF), "JOB": JOB, "CONTRACT": CONTRACT,
        "CONTRACT_SHA": args.contract_sha, "RUNNER": RUNNER, "WRAPPER": WRAPPER,
        "OUT": OUT, "BSM": BSM, "UNIT": UNIT,
    }.items():
        setattr(module, name, value)

    def checked_inputs(runner_sha, fresh=True):
        check_hash(CONTRACT, args.contract_sha)
        check_hash(RUNNER, runner_sha)
        check_hash(FINALIZER, args.finalizer_sha)
        for path in (OUT, BSM):
            module.exact_path(path)
        c = json.loads(CONTRACT.read_text())
        expected = {
            "status": "APPROVED_CONDITIONAL_CONTINUATION", "scientific_acceptance": "NONE",
            "backend": "native", "model": "M1", "KKT1": False,
            "pilot_successful_maps": 1, "minimum_successful_maps": 100,
            "first_stability_checkpoint": 200, "max_successful_maps": 500,
            "max_attempts": 1000, "seed_base": 202609140, "no_success_failure_stop": 10,
            "maxtries": 40000, "map_batch_size": 100,
            "bsm_path_sampler": "UNIFORMIZATION_NATIVE_EFFECTIVE_Q",
            "bridge_relative_tail_tolerance": 1e-12,
        }
        for key, wanted in expected.items():
            if c.get(key) != wanted:
                raise ValueError(f"Unexpected current R4 contract field: {key}")
        limits = c["resource_limits"]
        for key, wanted in {"unit": UNIT, "cpu_quota_cores": 4,
                            "memory_max_bytes": 16000000000, "tasks_max": 32,
                            "memory_swap_max_bytes": 0, "overall_walltime_seconds": 172800,
                            "bsm_bgb_cores": 1, "native_math_threads": 1}.items():
            if limits.get(key) != wanted:
                raise ValueError(f"Unexpected BSM resource contract: {key}")
        if c.get("bridge_native_cache"):
            raise ValueError("This R4 run must not import any historical preparation cache")
        pins = {**c["frozen_inputs"], **c["frozen_helpers"],
                c["fit_rds"]: c["fit_sha256"], c["postfit_rds"]: c["postfit_sha256"],
                str(BASE.relative_to(JOB)): BASE_SHA,
                str(FINALIZER.relative_to(JOB)): args.finalizer_sha}
        if c.get("numerical_audit"):
            item = c["numerical_audit"]
            pins[item["path"]] = item["sha256"]
        for item in c.get("bridge_required_evidence", {}).values():
            pins[item["path"]] = item["sha256"]
        for rel, expected_sha in pins.items():
            path = JOB / rel
            if not path.resolve().is_relative_to(JOB) or path.resolve() != path:
                raise ValueError("Source binding escapes the current R4 root")
            check_hash(path, expected_sha)
        if fresh and (OUT.exists() or BSM.exists()):
            raise ValueError("Output already exists; do not duplicate or restart this finite run")
        return c, pins

    def resource_guard():
        lines = Path("/proc/self/cgroup").read_text().splitlines()
        groups = [line.split(":", 2)[2] for line in lines if line.startswith("0::")]
        if len(groups) != 1 or not groups[0].endswith("/" + UNIT):
            raise ValueError("Run only in the dedicated detached BSM user service")
        cg = Path("/sys/fs/cgroup") / groups[0].lstrip("/")
        limits = {key: (cg / key).read_text().strip()
                  for key in ("memory.max", "memory.swap.max", "pids.max", "cpu.max")}
        quota, period = limits["cpu.max"].split()
        if not (limits["memory.max"] != "max" and 0 < int(limits["memory.max"]) <= 16000000000
                and limits["memory.swap.max"] == "0" and limits["pids.max"] != "max"
                and 0 < int(limits["pids.max"]) <= 32 and quota != "max"
                and 0 < int(quota) / int(period) <= 4):
            raise ValueError(f"Unsafe BSM cgroup: {limits}")
        r6_reservations()
        return cg, limits

    module.checked_inputs = checked_inputs
    module.resource_guard = resource_guard
    c, pins = checked_inputs(args.runner_sha)
    if args.selftest:
        tests = module.selftest(c, args.runner_sha)
        tests.update(current_R4_contract=True, postprocess_sha256=args.finalizer_sha)
        print(json.dumps(tests, indent=2))
        return 0
    launch = [
        "systemd-run", "--user", "--unit=" + UNIT,
        "--property=CPUQuota=400%", "--property=TasksMax=32",
        "--property=MemoryMax=16000000000", "--property=MemorySwapMax=0",
        "--property=RuntimeMaxSec=48h", "--property=KillMode=control-group",
        "--property=Restart=no", "--property=RemainAfterExit=yes",
        *["--setenv=" + key + "=1" for key in module.THREADS],
        sys.executable, str(SELF), "--runner-sha", args.runner_sha,
        "--contract-sha", args.contract_sha, "--finalizer-sha", args.finalizer_sha, "--run",
    ]
    if not args.run:
        print(json.dumps({
            "status": "READ_ONLY_PREFLIGHT_PASS", "scientific_acceptance": "NONE",
            "r6": r6_reservations(), "source_bindings": len(pins),
            "stages": [module.command(stage) for stage in module.STAGES],
            "suggested_systemd_command_not_executed": shlex.join(launch),
            "launch_argv_not_executed": launch,
            "completion_hook": [sys.executable, str(FINALIZER), "--run"],
            "restarts": "forbidden; outputs are unique and append-only per attempt",
        }, indent=2))
        return 0
    # The old main constructs a launch suggestion but cannot execute it. Its --run
    # branch uses the new guard and all rebound paths above, never old resources.
    sys.argv = [str(SELF), "--runner-sha", args.runner_sha, "--run"]
    result = module.main()
    state = json.loads((OUT / "STATUS.json").read_text())
    if result != 0 or state.get("status") != "COMPUTE_COMPLETE_CONDITIONAL_REVIEW_REQUIRED":
        return result or 1
    checked_inputs(args.runner_sha, fresh=False)
    hook_log = OUT / "postprocess.log"
    started = now()
    post_state = {"status": "POSTPROCESS_RUNNING", "started_utc": started,
                  "finalizer_sha256": args.finalizer_sha,
                  "bsm_status_preserved": state["status"], "scientific_acceptance": "NONE"}
    module.write_json(OUT / "POSTPROCESS_STATUS.json", post_state)
    try:
        with hook_log.open("x") as handle:
            completed = subprocess.run([sys.executable, str(FINALIZER), "--run"],
                                       cwd=JOB, stdin=subprocess.DEVNULL, stdout=handle,
                                       stderr=subprocess.STDOUT, timeout=7200, check=False)
    except BaseException as error:
        module.write_json(OUT / "POSTPROCESS_STATUS.json", {
            **post_state, "status": "POSTPROCESS_FAILED_BSM_PRESERVED",
            "finished_utc": now(), "error": repr(error),
        }, replace=True)
        raise
    module.write_json(OUT / "POSTPROCESS_STATUS.json", {
        **post_state,
        "status": "POSTPROCESS_EXIT_ZERO_REVIEW_REQUIRED" if completed.returncode == 0
                  else "POSTPROCESS_FAILED_BSM_PRESERVED",
        "started_utc": started, "finished_utc": now(), "exit_code": completed.returncode,
        "finalizer_sha256": args.finalizer_sha, "log_sha256": digest(hook_log),
        "bsm_status_preserved": state["status"], "scientific_acceptance": "NONE",
    }, replace=True)
    return completed.returncode


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)
