#!/usr/bin/env python3
"""Capture the strict native Faster discrete-moment regression as evidence.

The regression holds the quadric construction topology fixed, differentiates
the actual triangulated/Boolean material boundary, and compares the ten
moments against three central-difference steps and Richardson extrapolation.
It is deliberately a synthetic method fixture: it cannot verify an author
asset, a CAD export, or the paper's independent Boolean implementation.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib
from pathlib import Path
import subprocess
import sys
import time
import traceback

from make_it_stand.verification.run_evidence import (
    EXPLICIT_EFFECTIVE_SOURCE_POLICY, begin_run, check_freshness, finish_run,
)
from make_it_stand.verification.spin_common import canonical, parse


ROOT = Path(__file__).resolve().parents[3]
_BINARY = "quadric-discrete-moment-regression"
_MARKERS = (
    "historical-sphere supplemental_richardson_entry_relative=",
    "contained-positive supplemental_richardson_entry_relative=",
    "clipped-positive supplemental_richardson_entry_relative=",
    "contained-negative supplemental_richardson_entry_relative=",
    "clipped-negative supplemental_richardson_entry_relative=",
    "topology-transition below=",
    "topology-open-stratum label=below",
    "topology-open-stratum label=above",
    "hyperboloid cap-radius requested=",
    "hyperboloid discrete-map h=3.0000000000000001e-05 full_map_relative=",
    "hyperboloid discrete-map h=1.0000000000000001e-05 full_map_relative=",
    "hyperboloid discrete-map h=3.0000000000000001e-06 full_map_relative=",
    "[PASS] discrete quadric moment differential",
)


def _write(path: Path, value: object) -> None:
    with path.open("x", encoding="ascii") as stream:
        stream.write(canonical(value) + "\n")


def _native_from_build(build: Path):
    """Load exactly the extension colocated with the declared CMake build."""
    directory = (build / "python").resolve(strict=True)
    candidates = sorted(directory.glob("_make_it_stand*.so"))
    if len(candidates) != 1:
        raise ValueError("declared build does not contain exactly one native extension")
    sys.path.insert(0, str(directory))
    try:
        extension = importlib.import_module("_make_it_stand")
    finally:
        sys.path.pop(0)
    if Path(getattr(extension, "__file__", "")).resolve() != candidates[0].resolve():
        raise ValueError("loaded extension differs from declared build directory")
    return extension


def run(
    output_directory: Path, *, build_directory: Path, dependency_roots: dict[str, Path]
) -> dict[str, object]:
    output = output_directory.absolute()
    output.mkdir(parents=False, exist_ok=False)
    build = build_directory.resolve()
    executable = build / _BINARY
    if not executable.is_file() or not executable.stat().st_mode & 0o111:
        raise ValueError("required native Faster discrete-moment executable is unavailable")

    started, tick, begin, events = (
        datetime.now(timezone.utc).isoformat(), time.perf_counter(), None, []
    )
    try:
        inputs = parse(canonical({
            "schema": "faster-discrete-moment-inputs-v2",
            "classification": "SYNTHETIC_REGRESSION",
            "binary": {
                "path": str(executable),
                "sha256": hashlib.sha256(executable.read_bytes()).hexdigest(),
            },
            "claims": [
                "fixed_topology_discrete_ten_moment_jacobian_matches_finite_differences",
                "inverse_gauss_map_source_vertices_match_native_triangulation",
                "intersection_vertex_motion_is_exercised_for_clipped_fixtures",
                "hyperboloid_boolean_moment_map_matches_finite_differences_in_all_ten_directions",
                "hyperboloid_truncation_roots_and_cap_centroids_are_differentiated",
                "rank_deficient_quadric_derivatives_fail_closed_with_explicit_classification",
            ],
        }))
        input_path = output / "inputs.json"
        _write(input_path, inputs)
        extension = _native_from_build(build)
        begin = begin_run(
            repository_root=ROOT,
            runner_files=[Path(__file__)],
            checker_files=[
                ROOT / "tests/quadric_discrete_moment_regression.cpp",
                ROOT / "src/spin/QuadricDiscreteMomentDifferential.cpp",
                ROOT / "src/spin/QuadricInverseGaussMap.cpp",
                ROOT / "src/spin/QuadricSourceMotion.cpp",
                ROOT / "src/spin/QuadricMomentEvaluator.cpp",
                ROOT / "src/spin/QuadricTriangulator.cpp",
                ROOT / "include/spin/QuadricDiscreteMomentDifferential.h",
                ROOT / "include/spin/QuadricInverseGaussMap.h",
                ROOT / "include/spin/QuadricMomentEvaluator.h",
                ROOT / "include/spin/QuadricSourceMotion.h",
                ROOT / "include/spin/QuadricTriangulator.h",
                ROOT / "CMakeLists.txt",
            ],
            asset_files=[input_path, executable],
            build_configuration_files=[ROOT / "CMakeLists.txt", build / "CMakeCache.txt"],
            parameters=inputs,
            seed_protocol={"seed": 0, "protocol": "deterministic native discrete-moment regression"},
            dependency_roots=dependency_roots,
            native_modules=[extension],
            pure_python=False,
            source_policy=EXPLICIT_EFFECTIVE_SOURCE_POLICY,
        )
        _write(output / "begin.json", begin.as_dict())
        completed = subprocess.run(
            [str(executable)], cwd=build, capture_output=True, text=True, check=False
        )
        transcript = completed.stdout + completed.stderr
        transcript_path = output / "discrete-moments.stdout"
        with transcript_path.open("x", encoding="utf-8") as stream:
            stream.write(transcript)
        if completed.returncode != 0 or any(marker not in transcript for marker in _MARKERS):
            raise RuntimeError("native Faster discrete-moment regression failed")
        result = {
            "path": str(transcript_path),
            "sha256": hashlib.sha256(transcript.encode()).hexdigest(),
            "returncode": completed.returncode,
            "required_markers": list(_MARKERS),
        }
        scientific = {
            "schema": "faster-discrete-moment-scientific-v2",
            "classification": "SYNTHETIC_REGRESSION",
            "native_regression": result,
        }
        _write(output / "scientific.json", scientific)
        events.append({
            "iteration": 0,
            "status": "ok",
            "data": {"kind": "FASTER_DISCRETE_MOMENT_NATIVE_REGRESSION", **result},
        })
        evidence = finish_run(
            begin,
            events=events,
            output_files=sorted(
                path for path in output.rglob("*")
                if path.is_file() and path.name != "begin.json"
            ),
            measurements={"started_at": started, "elapsed_seconds": time.perf_counter() - tick},
        )
        _write(output / "evidence.json", evidence.as_dict())
        freshness = check_freshness(evidence)
        summary = {
            "status": "SYNTHETIC_DISCRETE_FULLRANK_QUADRIC_METHOD_VALIDATED",
            "classification": "SYNTHETIC_REGRESSION",
            "project_status": "PARTIAL_BLOCKED",
            "original_evidence_admission": "BLOCKED_UNKNOWN",
            "paper_case_admission": "BLOCKED_MISSING_AUTHOR_ASSETS_AND_PARAMETERS",
            "freshness": {"status": freshness.status, "mismatches": list(freshness.mismatches)},
        }
        _write(output / "summary.json", summary)
        return summary
    except Exception as error:
        if begin is not None and not (output / "evidence.json").exists():
            events.append({"iteration": len(events), "status": "failed", "data": {"error": str(error)}})
            _write(output / "evidence.json", finish_run(
                begin,
                events=events,
                output_files=sorted(
                    path for path in output.rglob("*")
                    if path.is_file() and path.name != "begin.json"
                ),
                measurements={"started_at": started, "elapsed_seconds": time.perf_counter() - tick},
            ).as_dict())
        _write(output / "failure.json", {
            "status": "FAILED", "error": str(error), "traceback": traceback.format_exc(),
        })
        raise


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-directory", type=Path, required=True)
    parser.add_argument("--build-directory", type=Path, required=True)
    parser.add_argument("--dependency-root", action="append", required=True, metavar="NAME=PATH")
    args = parser.parse_args()
    roots: dict[str, Path] = {}
    for item in args.dependency_root:
        name, separator, value = item.partition("=")
        if not name or not separator or name in roots:
            parser.error("invalid or duplicate dependency root")
        roots[name] = Path(value).resolve()
    result = run(args.output_directory, build_directory=args.build_directory, dependency_roots=roots)
    print(canonical(result))
    return 0 if result["freshness"]["status"] == "CURRENT" else 1


if __name__ == "__main__":
    raise SystemExit(main())
