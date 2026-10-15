#!/usr/bin/env python3
"""Run a preregistered G8 reconstructed-Heart Spin-It method experiment.

This is not an author asset or paper-case reproduction. It keeps the current
reconstructed Heart geometry and all prior solver parameters fixed; only the
source geometry depth and allowed optimization level are raised from 5/4 to
8/8. Every attempt gets a new output directory and records full evidence.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time
import traceback

import numpy as np

ROOT = Path(__file__).resolve().parents[3]


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write(path: Path, value: object) -> None:
    with path.open("x", encoding="ascii") as stream:
        stream.write(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                allow_nan=False) + "\n")


def _write_obj(path: Path, vertices: np.ndarray, faces: np.ndarray) -> None:
    with path.open("x", encoding="ascii") as stream:
        for vertex in vertices:
            stream.write("v %.17g %.17g %.17g\n" % tuple(map(float, vertex)))
        for face in faces:
            stream.write("f %d %d %d\n" % tuple(int(index) + 1 for index in face))


def build_protocol():
    """Build the frozen parent Heart input and G8-only level changes."""
    generator_dir = str(ROOT / "experiments/reconstructed_heart_top")
    if generator_dir not in sys.path:
        sys.path.insert(0, generator_dir)
    from heart_top import heart_mesh
    from make_it_stand.verification import source_spin_material as bridge
    from make_it_stand.verification.spin_problem import SpinArtifactProblem

    vertices, faces = heart_mesh()
    center = np.array([0.0, 0.0, 0.95])
    design_vertices = center + 0.72 * (vertices - center)
    shell = (np.vstack((vertices, design_vertices)),
             np.vstack((faces, faces[:, ::-1] + len(vertices))))

    source_settings = bridge.SourceGeometrySettings(
        root_origin=(-1.5, -1.5, -0.25), root_side=3.0,
        maximum_depth=8, maximum_nodes=5_000_000)
    solver_settings = bridge.SourceSolverSettings(
        objective="TOP_MASS_AND_LOWERING", formulation="PAPER_GIVENS",
        initial_level=3, maximum_optimization_level=8,
        maximum_adaptive_passes=10, maximum_leaf_count=4000,
        maximum_output_cuboids=1_000_000, initial_fill=0.5,
        initial_phi=0.2, laplacian_weight=0.001,
        slqp=bridge.AuditedSlqpSettings(maximum_iterations=160))
    boundary_settings = bridge.SourceBoundarySettings(maximum_input_cuboids=4000)
    problem = SpinArtifactProblem(
        original=(vertices, faces), fixed_regions=[("shell", shell, 1.0)],
        design_regions=[("design", (design_vertices, faces), 1.0)],
        length_unit="25mm", mass_unit="18.75g",
        target={"mode": "SPINNING_TOP", "axis": [0, 0, 1],
                "origin": [0, 0, 0], "support_points": []},
        minimum_outer_wall=0.015)
    return problem, source_settings, solver_settings, boundary_settings


def _observed_levels(native_run: dict) -> dict:
    counts: dict[int, int] = {}
    passes = native_run.get("solve", {}).get("passes", [])
    for pass_record in passes:
        for field in ("level_counts", "output_level_counts",
                      "fractional_level_counts", "output_fractional_level_counts"):
            values = pass_record.get(field, {})
            if isinstance(values, dict):
                entries = ((int(level), int(count))
                           for level, count in values.items())
            elif isinstance(values, list):
                entries = enumerate(values)
            else:
                continue
            for level, count in entries:
                if int(count) > 0:
                    counts[int(level)] = max(counts.get(int(level), 0), int(count))
    return {"maximum_observed_level": max(counts, default=None),
            "level_counts_union": {str(level): counts[level] for level in sorted(counts)}}


def declared_dependency_roots() -> dict[str, Path]:
    """Return actual, bounded build-dependency roots captured by run evidence.

    The loaded extension itself is separately hashed in ``native_modules``.
    pybind11 headers are the external build dependency used by this extension;
    this declaration does not claim to enumerate all transitive system/runtime
    dependencies.
    """
    try:
        import pybind11
    except ImportError as error:
        raise RuntimeError("pybind11 build dependency is unavailable") from error
    root = Path(pybind11.get_include()).resolve(strict=True)
    if not root.is_dir() or not (root / "pybind11" / "pybind11.h").is_file():
        raise RuntimeError("pybind11 include root is incomplete")
    return {"pybind11": root}


def evidence_parameters(problem_identity: dict, config: dict) -> dict:
    """Canonicalize evidence inputs through strict JSON before sealing them."""
    from make_it_stand.verification.spin_common import canonical, parse

    return parse(canonical({
        "schema": "reconstructed-heart-g8-input-v1",
        "classification": "RECONSTRUCTED_METHOD_VALIDATION",
        "problem_identity": problem_identity,
        "config": config,
    }))


def run(output_directory: Path, *, build_directory: Path) -> dict:
    from make_it_stand.verification import source_spin_material as bridge
    from make_it_stand.verification.run_evidence import (
        EXPLICIT_EFFECTIVE_SOURCE_POLICY, begin_run, check_freshness, finish_run)
    from make_it_stand.verification.spin_bundle import save_spin_artifact, load_spin_artifact
    from make_it_stand.verification.spin_common import native

    output = output_directory.absolute()
    output.mkdir(parents=False, exist_ok=False)
    build = build_directory.resolve(strict=True)
    extension = native()
    extension_path = Path(extension.__file__).resolve(strict=True)
    if extension_path.parents[1] != build:
        raise ValueError("declared build directory does not contain the loaded native extension")

    problem, source_settings, solver_settings, boundary_settings = build_protocol()
    started = datetime.now(timezone.utc).isoformat()
    tick = time.perf_counter()
    events = []
    begin = None
    try:
        vertices, faces = problem.original.arrays()
        design = problem.design_regions[0].mesh.arrays()
        shell = problem.fixed_regions[0].mesh.arrays()
        config = {
            "asset": "RECONSTRUCTED_SELF_BUILT_HEART_V1",
            "original_paper_case": False,
            "source_geometry": asdict(source_settings),
            "solver": asdict(solver_settings),
            "boundary": asdict(boundary_settings),
            "one_factor_change_from_parent": {
                "source_geometry_maximum_depth": [5, 8],
                "maximum_optimization_level": [4, 8],
            },
            "resource_limit": {"memory_max_bytes": 24 * 1024**3,
                               "memory_swap_max_bytes": 0,
                               "wall_time_limit_seconds": 3600},
        }
        _write(output / "config.json", config)
        _write(output / "input.json", problem.identity)
        _write_obj(output / "outer.obj", vertices, faces)
        _write_obj(output / "design.obj", *design)
        _write_obj(output / "fixed-shell.obj", *shell)

        inputs = evidence_parameters(problem.identity, config)
        paths = [Path(__file__), Path(__file__).parents[2] / "reconstructed_heart_top" / "heart_top.py",
                 extension_path, ROOT / "python/make_it_stand/verification/source_spin_material.py",
                 ROOT / "python/make_it_stand/verification/spin_artifact.py",
                 ROOT / "python/make_it_stand/verification/mass_moments.py",
                 ROOT / "CMakeLists.txt", build / "CMakeCache.txt"]
        source = {str(path.resolve()): _sha256(path.resolve(strict=True)) for path in paths}
        _write(output / "source.json", {
            "git_commit": subprocess.check_output(["git", "rev-parse", "HEAD"],
                                                   cwd=ROOT, text=True).strip(),
            "files": source,
            "scope": "LOADED_BINARY_AND_DECLARED_FILES; DOES_NOT_ATTEST_ALL_BUILD_INPUTS"})
        begin = begin_run(
            repository_root=ROOT,
            runner_files=[Path(__file__)],
            checker_files=[ROOT / "python/make_it_stand/verification/source_spin_material.py",
                           ROOT / "python/make_it_stand/verification/spin_artifact.py",
                           ROOT / "python/make_it_stand/verification/spin_problem.py",
                           ROOT / "python/make_it_stand/verification/mass_moments.py"],
            asset_files=[output / "config.json", output / "input.json", output / "outer.obj",
                         output / "design.obj", output / "fixed-shell.obj"],
            build_configuration_files=[ROOT / "CMakeLists.txt", build / "CMakeCache.txt"],
            parameters=inputs,
            seed_protocol={"seed": 0, "protocol": "frozen source-spin2014 G8 parent protocol"},
            dependency_roots=declared_dependency_roots(), native_modules=[extension], pure_python=False,
            source_policy=EXPLICIT_EFFECTIVE_SOURCE_POLICY)
        _write(output / "begin.json", begin.as_dict())

        result = bridge.solve_and_seal_source_spin2014(
            problem, source_settings=source_settings, solver_settings=solver_settings,
            boundary_settings=boundary_settings)
        _write(output / "native-run.json", result.identity["native_run"])
        _write(output / "bridge-observations.json", result.evidence)
        if result.material is not None and len(result.material[0]):
            _write_obj(output / "material.obj", *result.material)
        if result.cavity is not None and len(result.cavity[0]):
            _write_obj(output / "cavity.obj", *result.cavity)

        sealed = False
        if result.artifact is not None:
            save_spin_artifact(result.artifact, output / "sealed")
            reloaded = load_spin_artifact(output / "sealed",
                                          expected_problem_digest=problem.input_digest)
            if reloaded.scientific_digest != result.artifact.scientific_digest:
                raise ValueError("sealed Heart artifact changed on reload")
            sealed = True

        observed = _observed_levels(result.identity["native_run"])
        summary = {
            "status": result.status, "sealed": sealed,
            "seal_error": result.seal_error,
            "scientific_digest": result.scientific_digest,
            "problem_digest": problem.input_digest,
            "observed_levels": observed,
            "g8_optimization_observed": observed["maximum_observed_level"] == 8,
            "asset_class": "RECONSTRUCTED_METHOD_VALIDATION",
            "paper_case_credit": False,
            "project_status": "PARTIAL_BLOCKED",
            "elapsed_seconds": time.perf_counter() - tick,
        }
        _write(output / "summary.json", summary)
        events.append({"iteration": 0,
                       "status": "ok" if sealed else "failed",
                       "data": {"kind": "RECONSTRUCTED_HEART_G8_ACTUAL_SOLVE",
                                "status": result.status,
                                "sealed": sealed,
                                "observed_levels": observed}})
        evidence = finish_run(
            begin, events=events,
            output_files=sorted(path for path in output.rglob("*")
                                if path.is_file() and path.name != "begin.json"),
            measurements={"started_at": started,
                          "elapsed_seconds": time.perf_counter() - tick})
        _write(output / "evidence.json", evidence.as_dict())
        freshness = check_freshness(evidence)
        summary["freshness"] = {"status": freshness.status,
                                "mismatches": list(freshness.mismatches)}
        _write(output / "final-summary.json", summary)
        return summary
    except Exception as error:
        failure = {"status": "FAILED", "error": str(error),
                   "traceback": traceback.format_exc(),
                   "elapsed_seconds": time.perf_counter() - tick}
        _write(output / "failure.json", failure)
        if begin is not None and not (output / "evidence.json").exists():
            events.append({"iteration": len(events), "status": "failed",
                           "data": {"error": str(error)}})
            evidence = finish_run(
                begin, events=events,
                output_files=sorted(path for path in output.rglob("*")
                                    if path.is_file() and path.name != "begin.json"),
                measurements={"started_at": started,
                              "elapsed_seconds": time.perf_counter() - tick})
            _write(output / "evidence.json", evidence.as_dict())
        raise


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--build-directory", type=Path, required=True)
    args = parser.parse_args()
    result = run(args.output, build_directory=args.build_directory)
    print(json.dumps(result, sort_keys=True, allow_nan=False))
    return 0 if result.get("sealed") and result.get("g8_optimization_observed") else 2


if __name__ == "__main__":
    raise SystemExit(main())
