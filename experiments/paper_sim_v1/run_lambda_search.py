#!/usr/bin/env python3
"""Bounded full-resolution λ reconstruction under an explicit claim downgrade."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import subprocess
import sys

from make_it_stand.paper_evidence import capture_execution_fingerprint


COARSE_LAMBDAS = (5.0, 10.0, 20.0, 40.0, 80.0)


def _token(value: float) -> str:
    return format(value, ".12g").replace(".", "p")


def _candidate_is_eligible(process_returncode: int, run: dict) -> bool:
    return bool(
        process_returncode == 0
        and run.get("completion") == "stable"
        and run.get("sealed_artifact_status") == "SEALED"
        and run.get("sealed_safe_region_state") == "SAFE"
    )


def _run_candidate(
    *,
    beta_root: Path,
    case: str,
    value: float,
    output: Path,
    replacement_surface: Path | None,
) -> dict:
    root = output / f"lambda-{_token(value)}"
    command = [
        sys.executable,
        str(Path(__file__).resolve().with_name("run_official_case.py")),
        "--beta-root", str(beta_root),
        "--case", case,
        "--profile", "paper",
        "--shape-scale-lambda", str(value),
        "--maximum-carving-blocks", "320",
        "--maximum-deformation-iterations", "320",
        "--time-limit-seconds", "900",
        "--memory-limit-gib", "24",
        "--repeat", "1",
        "--output", str(root),
    ]
    if replacement_surface is not None:
        command.extend((
            "--replacement-normalized-surface",
            str(replacement_surface.resolve()),
        ))
    completed = subprocess.run(command, check=False)
    manifest_path = root / "manifest.json"
    manifest = (
        json.loads(manifest_path.read_text(encoding="utf-8"))
        if manifest_path.is_file()
        else None
    )
    run = manifest["runs"][0] if manifest and manifest.get("runs") else {}
    diagnostics_path = root / "repeat-1/diagnostics.json"
    diagnostics = (
        json.loads(diagnostics_path.read_text(encoding="utf-8"))
        if diagnostics_path.is_file()
        else {}
    )
    metrics = diagnostics.get("r10_metrics", {})
    eligible = _candidate_is_eligible(completed.returncode, run)
    return {
        "lambda": value,
        "eligible": eligible,
        "process_returncode": completed.returncode,
        "completion": run.get("completion", "missing_manifest"),
        "sealed_artifact_status": run.get("sealed_artifact_status"),
        "sealed_safe_region_state": run.get("sealed_safe_region_state"),
        "sealed_safe_region_residual": run.get("sealed_safe_region_residual"),
        "normalized_rms_vertex_displacement": metrics.get(
            "normalized_rms_vertex_displacement"
        ),
        "final_projected_target_energy": metrics.get(
            "final_projected_target_energy"
        ),
        "result_directory": str(root.resolve()),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--beta-root", required=True, type=Path)
    parser.add_argument(
        "--case",
        required=True,
        choices=("spheres", "armadillo1", "dinosaur", "horse"),
    )
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--horse-replacement", type=Path)
    arguments = parser.parse_args()
    if arguments.case == "horse" and arguments.horse_replacement is None:
        parser.error(
            "Horse λ calibration requires --horse-replacement and remains "
            "classified as RECONSTRUCTED_ASSET"
        )
    arguments.output.mkdir(parents=True, exist_ok=True)
    evaluated = [
        _run_candidate(
            beta_root=arguments.beta_root.resolve(),
            case=arguments.case,
            value=value,
            output=arguments.output,
            replacement_surface=arguments.horse_replacement,
        )
        for value in COARSE_LAMBDAS
    ]
    eligible = [item for item in evaluated if item["eligible"]]
    if eligible:
        winner = min(
            eligible,
            key=lambda item: (
                float(item["normalized_rms_vertex_displacement"]),
                float(item["final_projected_target_energy"]),
                float(item["lambda"]),
            ),
        )
        index = COARSE_LAMBDAS.index(float(winner["lambda"]))
        neighbors = []
        if index > 0:
            neighbors.append(math.sqrt(COARSE_LAMBDAS[index - 1] * winner["lambda"]))
        if index + 1 < len(COARSE_LAMBDAS):
            neighbors.append(math.sqrt(COARSE_LAMBDAS[index + 1] * winner["lambda"]))
        evaluated.extend(
            _run_candidate(
                beta_root=arguments.beta_root.resolve(),
                case=arguments.case,
                value=value,
                output=arguments.output,
                replacement_surface=arguments.horse_replacement,
            )
            for value in neighbors
        )
        eligible = [item for item in evaluated if item["eligible"]]
    selected = (
        min(
            eligible,
            key=lambda item: (
                float(item["normalized_rms_vertex_displacement"]),
                float(item["final_projected_target_energy"]),
                float(item["lambda"]),
            ),
        )
        if eligible
        else None
    )
    repository_root = Path(__file__).resolve().parents[2]
    candidate_manifests = [
        Path(item["result_directory"]) / "manifest.json" for item in evaluated
    ]
    fingerprint = capture_execution_fingerprint(
        repository_root,
        source_paths=[
            repository_root / "CMakeLists.txt",
            repository_root / "bindings",
            repository_root / "include",
            repository_root / "src",
            repository_root / "python",
            Path(__file__).resolve(),
            Path(__file__).resolve().with_name("run_official_case.py"),
            *(path for path in candidate_manifests if path.is_file()),
        ],
        executable_paths=[Path(sys.executable).resolve()],
        configuration={
            "case": arguments.case,
            "coarse_lambdas": list(COARSE_LAMBDAS),
            "horse_replacement": (
                str(arguments.horse_replacement.resolve())
                if arguments.horse_replacement is not None
                else None
            ),
        },
    )
    payload = {
        "schema": "make-it-stand-lambda-reconstruction-v1",
        "schema_version": 1,
        "case": arguments.case,
        "classification": (
            "RECONSTRUCTED_ASSET"
            if arguments.horse_replacement is not None
            else "RECONSTRUCTED_PARAMETER"
        ),
        "freshness": {"status": "CURRENT", "mismatches": []},
        "execution_fingerprint": fingerprint,
        "full_resolution_only": True,
        "maximum_outer_iterations": 320,
        "wall_clock_limit_seconds_per_candidate": 900,
        "memory_limit_gib": 24,
        "selection_order": [
            "stable_and_sealed",
            "minimum_normalized_rms_vertex_displacement",
            "minimum_final_projected_target_energy",
            "minimum_lambda",
        ],
        "selected": selected,
        "candidates": evaluated,
    }
    (arguments.output / "lambda-search.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return 0 if selected is not None else 2


if __name__ == "__main__":
    raise SystemExit(main())
