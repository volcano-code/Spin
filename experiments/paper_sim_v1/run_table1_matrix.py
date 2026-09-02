#!/usr/bin/env python3
"""Run the five Table 1 cases without mixing evidence classes or profiles."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
from typing import Mapping

from make_it_stand.paper_evidence import capture_execution_fingerprint


CASE_ORDER = ("spheres", "dinosaur", "armadillo1", "horse", "toyrocking")


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build_table_plan(
    *,
    beta_root: Path,
    output: Path,
    profile: str,
    lambdas: Mapping[str, float],
    repository_root: Path,
) -> list[dict[str, object]]:
    if set(lambdas) != set(CASE_ORDER):
        missing = sorted(set(CASE_ORDER) - set(lambdas))
        extra = sorted(set(lambdas) - set(CASE_ORDER))
        raise ValueError(f"lambda map mismatch; missing={missing}, extra={extra}")
    script_root = repository_root / "experiments" / "paper_sim_v1"
    plan = []
    for case in CASE_ORDER:
        value = float(lambdas[case])
        if not math.isfinite(value) or value < 0.0:
            raise ValueError(f"{case}: lambda must be finite and nonnegative")
        driver = script_root / (
            "run_toy_reconstruction.py"
            if case == "toyrocking"
            else "run_official_case.py"
        )
        command = [
            sys.executable,
            str(driver),
            "--beta-root", str(beta_root),
            "--profile", profile,
            "--shape-scale-lambda", format(value, ".17g"),
            "--maximum-carving-blocks", "400",
            "--maximum-deformation-iterations", "400",
            "--time-limit-seconds", "3600",
            "--memory-limit-gib", "24",
            "--repeat", "2",
            "--output", str(output / case),
        ]
        if case != "toyrocking":
            command[command.index("--profile"):command.index("--profile")] = [
                "--case", case
            ]
        if case == "horse":
            command.extend((
                "--replacement-normalized-surface",
                str(
                    script_root
                    / "reconstructed_assets"
                    / "horse_100_occupancy_outer.obj"
                ),
            ))
        plan.append({
            "case": case,
            "driver": str(driver),
            "lambda": value,
            "command": command,
        })
    return plan


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--beta-root", required=True, type=Path)
    parser.add_argument("--profile", required=True, choices=("paper", "released-beta"))
    parser.add_argument("--lambda-map", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    repository_root = Path(__file__).resolve().parents[2]
    beta_root = arguments.beta_root.expanduser().resolve()
    output = arguments.output.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    lambda_payload = json.loads(
        arguments.lambda_map.read_text(encoding="utf-8")
    )
    if not isinstance(lambda_payload, dict):
        parser.error("--lambda-map must contain a JSON object keyed by case")
    plan = build_table_plan(
        beta_root=beta_root,
        output=output,
        profile=arguments.profile,
        lambdas=lambda_payload,
        repository_root=repository_root,
    )
    records = []
    all_pass = True
    for item in plan:
        completed = subprocess.run(item["command"], check=False)
        case_root = output / str(item["case"])
        manifest_path = case_root / "manifest.json"
        manifest = (
            json.loads(manifest_path.read_text(encoding="utf-8"))
            if manifest_path.is_file()
            else {}
        )
        accepted = bool(
            completed.returncode == 0
            and manifest.get("repeat_count_requested") == 2
            and manifest.get("deterministic_output_hashes") is True
            and len(manifest.get("runs", [])) == 2
            and all(
                run.get("completion") == "stable"
                and run.get("sealed_artifact_status") == "SEALED"
                and run.get("sealed_safe_region_state") == "SAFE"
                and run.get("sealed_outer_flattening_repair_status") in {
                    "NOT_REQUIRED",
                    "APPLIED",
                }
                and run.get("sealed_wall_thickness_status") == "PASS"
                and run.get("sealed_containment_status") == "PASS"
                for run in manifest.get("runs", [])
            )
        )
        all_pass = all_pass and accepted
        records.append({
            "case": item["case"],
            "lambda": item["lambda"],
            "process_returncode": completed.returncode,
            "accepted": accepted,
            "manifest": str(manifest_path),
            "manifest_sha256": _sha256(manifest_path) if manifest_path.is_file() else None,
            "classification": manifest.get("classification"),
        })
    fingerprint = capture_execution_fingerprint(
        repository_root,
        source_paths=[
            repository_root / "CMakeLists.txt",
            repository_root / "CMakePresets.json",
            repository_root / "bindings",
            repository_root / "include",
            repository_root / "src",
            repository_root / "python",
            Path(__file__).resolve(),
            Path(__file__).resolve().with_name("run_official_case.py"),
            Path(__file__).resolve().with_name("run_toy_reconstruction.py"),
        ],
        executable_paths=[Path(sys.executable).resolve()],
        configuration={
            "profile": arguments.profile,
            "lambdas": lambda_payload,
            "case_order": list(CASE_ORDER),
        },
    )
    report = {
        "schema": "make-it-stand-table1-matrix-v1",
        "schema_version": 1,
        "profile": arguments.profile,
        "freshness": {"status": "CURRENT", "mismatches": []},
        "execution_fingerprint": fingerprint,
        "all_pass": all_pass,
        "case_order": list(CASE_ORDER),
        "cases": records,
        "claim_boundary": (
            "All lambdas are reconstructed parameters; horse and rocking toy "
            "are additionally reconstructed assets."
        ),
    }
    (output / "table1-matrix.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return 0 if all_pass else 2


if __name__ == "__main__":
    raise SystemExit(main())
