#!/usr/bin/env python3
"""Capture reduced rank-deficient quadric source/Boolean differentials."""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

from make_it_stand.verification.run_evidence import (
    EXPLICIT_EFFECTIVE_SOURCE_POLICY,
    begin_run,
    check_freshness,
    finish_run,
)
from make_it_stand.verification.spin_common import canonical


ROOT = Path(__file__).resolve().parents[3]
MARKERS = (
    "[PASS] native source all coefficients, hyperboloid=1",
    "[PASS] native source all coefficients, hyperboloid=0",
    "[PASS] doubled radius, rank-deficient plane source and budget rejected",
    "[PASS] rank-deficient source explicitly rejected type=ELLIPTIC_CYLINDER",
    "[PASS] rank-deficient source explicitly rejected type=PARABOLOID",
    "[PASS] paraboloid source motion in reduced stratum coordinates",
    "[PASS] clipped paraboloid Boolean moments in reduced stratum coordinates",
    "[PASS] paraboloid eigenvector sign-pivot boundary rejected",
    "[PASS] same-sign eigenvector magnitude tie remains admissible",
    "[PASS] cylinder eigenvector sign-pivot boundary rejected",
)


def write(path: Path, value: object) -> None:
    path.write_text(canonical(value) + "\n", encoding="ascii")


def run(output: Path, executable: Path, dependency_roots: dict[str, Path]) -> dict[str, object]:
    output = output.absolute()
    output.mkdir(parents=False, exist_ok=False)
    executable = executable.resolve(strict=True)
    started = time.perf_counter()
    inputs = {
        "schema": "faster-source-stratum-differential-inputs-v3",
        "classification": "SYNTHETIC_REGRESSION",
        "executable": {
            "path": str(executable),
            "sha256": hashlib.sha256(executable.read_bytes()).hexdigest(),
        },
        "ambient_eq14_policy": "reject_rank_deficient_source_differential_without_fabricating_10D_jacobian",
        "reduced_strata": {
            "supported": {"ELLIPTIC_CYLINDER": 7, "PARABOLOID": 8},
            "unsupported": ["PLANE"],
            "result": "10_moment_rows_by_actual_stratum_dimension_columns",
        },
        "checks": [
            "paraboloid_source_grid_and_extrusion_vertex_motion_vs_central_finite_differences",
            "actual_cgal_boolean_material_moments_and_intersection_vertex_motion_vs_finite_differences",
            "fixed_topology_compared_by_birth_feature_provenance",
            "paraboloid_and_cylinder_eigenvector_sign_pivot_ties_fail_closed",
            "same_sign_equal_magnitude_pivots_remain_admissible",
        ],
    }
    input_path = output / "inputs.json"
    write(input_path, inputs)
    begin = begin_run(
        repository_root=ROOT,
        runner_files=[Path(__file__)],
        checker_files=[
            ROOT / "tests/quadric_source_motion_regression.cpp",
            ROOT / "src/spin/QuadricSourceMotion.cpp",
            ROOT / "src/spin/QuadricDiscreteMomentDifferential.cpp",
            ROOT / "src/spin/QuadricStratumTangent.cpp",
            ROOT / "include/spin/QuadricSourceMotion.h",
            ROOT / "include/spin/QuadricDiscreteMomentDifferential.h",
            ROOT / "include/spin/QuadricStratumTangent.h",
        ],
        asset_files=[input_path, executable],
        build_configuration_files=[ROOT / "CMakeLists.txt"],
        parameters=inputs,
        seed_protocol={"seed": 0, "protocol": "deterministic source-stratum guard"},
        dependency_roots=dependency_roots,
        native_modules=[],
        pure_python=True,
        source_policy=EXPLICIT_EFFECTIVE_SOURCE_POLICY,
    )
    write(output / "begin.json", begin.as_dict())
    completed = subprocess.run(
        [str(executable)], cwd=output, capture_output=True, text=True, check=False
    )
    transcript = completed.stdout + completed.stderr
    transcript_path = output / "source-stratum-differential.stdout"
    transcript_path.write_text(transcript, encoding="utf-8")
    missing = [marker for marker in MARKERS if marker not in transcript]
    if completed.returncode != 0 or missing:
        raise RuntimeError(
            f"source-stratum differential run failed: returncode={completed.returncode} missing={missing}"
        )
    scientific = {
        "schema": "faster-source-stratum-differential-scientific-v3",
        "classification": "SYNTHETIC_REGRESSION",
        "scope": "REDUCED_PARABOLOID_SOURCE_AND_FIXED_TOPOLOGY_BOOLEAN_MOMENT_DIFFERENTIAL_NOT_EQ14_SOLVER_OR_AUTHOR_CASE",
        "stdout_sha256": hashlib.sha256(transcript.encode()).hexdigest(),
        "required_markers": list(MARKERS),
    }
    write(output / "scientific.json", scientific)
    evidence = finish_run(
        begin,
        events=[
            {
                "iteration": 0,
                "status": "ok",
                "data": {"kind": "FASTER_REDUCED_PARABOLOID_BOOLEAN_DIFFERENTIAL"},
            }
        ],
        output_files=[transcript_path, output / "scientific.json"],
        measurements={"elapsed_seconds": time.perf_counter() - started},
    )
    write(output / "evidence.json", evidence.as_dict())
    freshness = check_freshness(evidence)
    summary = {
        "status": "SYNTHETIC_REDUCED_PARABOLOID_BOOLEAN_DIFFERENTIAL_VALIDATED",
        "classification": "SYNTHETIC_REGRESSION",
        "project_status": "PARTIAL_BLOCKED",
        "original_evidence_admission": "BLOCKED_MISSING_AUTHOR_ASSETS_AND_PARAMETERS",
        "freshness": {"status": freshness.status, "mismatches": list(freshness.mismatches)},
    }
    write(output / "summary.json", summary)
    return summary


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    parser.add_argument("--dependency-root", action="append", required=True)
    args = parser.parse_args()
    roots: dict[str, Path] = {}
    for item in args.dependency_root:
        name, separator, value = item.partition("=")
        if not name or not separator or name in roots:
            parser.error("invalid or duplicate dependency root")
        roots[name] = Path(value).resolve()
    try:
        result = run(args.output_directory, args.executable, roots)
    except Exception as error:
        print(str(error), file=sys.stderr)
        return 2
    print(canonical(result))
    return 0 if result["freshness"]["status"] == "CURRENT" else 1


if __name__ == "__main__":
    raise SystemExit(main())
