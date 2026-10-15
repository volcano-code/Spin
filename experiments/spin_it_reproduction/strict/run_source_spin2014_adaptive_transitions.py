#!/usr/bin/env python3
"""Record a source-surface adaptive split trajectory without promoting success.

This is a deliberately small ``SYNTHETIC_REGRESSION`` fixture.  Its intended
outcome is ``BINARY_VALIDATION_FAILED``: the run demonstrates source-derived
octree refinement and preserves the rounded constraint failure as evidence.
It is neither an author asset nor a paper-case result.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
import time
import traceback

import numpy as np

from make_it_stand.verification import source_spin_material as bridge
from make_it_stand.verification.run_evidence import (
    EXPLICIT_EFFECTIVE_SOURCE_POLICY, begin_run, check_freshness, finish_run,
)
from make_it_stand.verification.spin_common import canonical, native, parse
from make_it_stand.verification.spin_problem import SpinArtifactProblem


ROOT = Path(__file__).resolve().parents[3]


def _box(lo, hi):
    vertices = np.array(
        [[hi[axis] if index & (1 << axis) else lo[axis] for axis in range(3)]
         for index in range(8)], dtype="<f8"
    )
    faces = np.array(
        [[0, 2, 1], [1, 2, 3], [4, 5, 6], [5, 7, 6], [0, 1, 4], [1, 5, 4],
         [2, 6, 3], [3, 6, 7], [0, 4, 2], [2, 4, 6], [1, 3, 5], [3, 7, 5]],
        dtype="<i4",
    )
    return vertices, faces


def fixture():
    # This is the asymmetric source-surface fixture used by the native
    # source-adaptive regression: it forces fractional leaves to split, then
    # intentionally retains a binary constraint residual after paper stopping.
    offset = 1.0 / 8192.0
    outer = _box((-.5 + offset, -.5, 0.), (.5 + offset, .5, 1.))
    design = _box((.125, -.125, .625), (.375, .125, .875))
    shell = (
        np.vstack((outer[0], design[0])),
        np.vstack((outer[1], design[1][:, ::-1] + len(outer[0]))),
    )
    problem = SpinArtifactProblem(
        original=outer,
        fixed_regions=[("shell", shell, 1.0)],
        design_regions=[("design", design, 1.0)],
        length_unit="m", mass_unit="kg",
        target={"mode": "STATIC_BALANCE", "axis": [0, 0, 1],
                "origin": [0, 0, 0], "support_points": [[0, 0, 0]]},
        minimum_outer_wall=.125,
    )
    settings = {
        "source_settings": bridge.SourceGeometrySettings(
            root_origin=(-1., -1., 0.), root_side=2., maximum_depth=5),
        "solver_settings": bridge.SourceSolverSettings(
            initial_level=2, maximum_optimization_level=5,
            maximum_adaptive_passes=3, laplacian_weight=1.,
            termination_policy="PAPER_POST_TRANSITION_V1",
            slqp=bridge.AuditedSlqpSettings(
                maximum_iterations=300, equality_tolerance=1e-10,
                projected_gradient_tolerance=1e-8)),
        "boundary_settings": bridge.SourceBoundarySettings(),
    }
    return problem, settings


def _write(path: Path, value) -> None:
    with path.open("x", encoding="ascii") as stream:
        stream.write(canonical(value) + "\n")


def _expected_transition(native_run: dict) -> None:
    solve = native_run.get("solve", {})
    passes = solve.get("passes")
    if solve.get("status") != "BINARY_VALIDATION_FAILED" or \
            solve.get("termination_policy") != "PAPER_POST_TRANSITION_V1" or \
            not isinstance(passes, list) or len(passes) != 3:
        raise ValueError("fixture did not produce the declared paper-transition binary failure")
    if not all(item.get("status") == "CONVERGED" for item in passes) or \
            not any(item.get("split", {}).get("changes") for item in passes) or \
            any(item.get("merge", {}).get("changes") for item in passes) or \
            passes[-1].get("paper_resolution_complete") is not True:
        raise ValueError("fixture did not preserve the required source split trajectory")
    stats = native_run.get("source", {}).get("stats", {})
    if not isinstance(stats.get("boundary_nodes"), int) or stats["boundary_nodes"] <= 0 or \
            stats.get("maximum_visited_level") != 5:
        raise ValueError("fixture lacks source-surface boundary-subtree evidence")


def run(output_directory, *, build_directory, dependency_roots):
    output = Path(output_directory).absolute()
    output.mkdir(parents=False, exist_ok=False)
    started = datetime.now(timezone.utc).isoformat()
    tick = time.perf_counter()
    evidence_start = None
    events = []
    try:
        problem, settings = fixture()
        inputs = parse(canonical({
            "schema": "source-spin2014-adaptive-transition-inputs-v1",
            "classification": "SYNTHETIC_REGRESSION",
            "expected_outcome": "BINARY_VALIDATION_FAILED_AFTER_REAL_SPLIT",
            "problem_identity": problem.identity,
            "settings": {name: asdict(value) for name, value in settings.items()},
        }))
        assets = []
        for label, mesh in (("outer", problem.original),
                            ("fixed_shell", problem.fixed_regions[0].mesh),
                            ("design", problem.design_regions[0].mesh)):
            for kind, array in zip(("vertices", "faces"), mesh.arrays()):
                path = output / f"{label}_{kind}.npy"
                with path.open("xb") as stream:
                    np.lib.format.write_array(stream, array, version=(1, 0), allow_pickle=False)
                assets.append(path)
        input_path = output / "inputs.json"
        _write(input_path, inputs)
        assets.append(input_path)
        extension = native()
        build = Path(build_directory).resolve()
        if Path(extension.__file__).resolve().parents[1] != build:
            raise ValueError("declared build directory does not contain the actually loaded extension")
        evidence_start = begin_run(
            repository_root=ROOT, runner_files=[Path(__file__)],
            checker_files=[Path(bridge.__file__), ROOT / "python/make_it_stand/verification/spin_problem.py"],
            asset_files=assets,
            build_configuration_files=[ROOT / "CMakeLists.txt", build / "CMakeCache.txt"],
            parameters=inputs,
            seed_protocol={"seed": 0, "protocol": "native audited HiGHS simplex: threads=1, parallel=off, random_seed=0"},
            dependency_roots=dependency_roots, native_modules=[extension], pure_python=False,
            source_policy=EXPLICIT_EFFECTIVE_SOURCE_POLICY,
        )
        _write(output / "begin.json", evidence_start.as_dict())
        inspection = extension.inspect_source_spin2014_domain(
            *problem.original.arrays(), *problem.design_regions[0].mesh.arrays(),
            source_settings=asdict(settings["source_settings"]),
            initial_level=settings["solver_settings"].initial_level,
            maximum_optimization_level=settings["solver_settings"].maximum_optimization_level,
            maximum_leaf_count=settings["solver_settings"].maximum_leaf_count,
        )
        _write(output / "domain-inspection.json", inspection)
        result = bridge.solve_and_seal_source_spin2014(problem, **settings)
        native_run = result.identity["native_run"]
        _expected_transition(native_run)
        if inspection["source"]["source_hash"] != native_run["source"]["source_hash"] or \
                inspection["source"]["partition_hash"] != native_run["source"]["partition_hash"]:
            raise ValueError("read-only domain inspection differs from actual transition source")
        _write(output / "native-run.json", native_run)
        scientific = {
            "schema": "source-spin2014-adaptive-transition-scientific-v1",
            "classification": "SYNTHETIC_REGRESSION",
            "scientific_digest": result.scientific_digest,
            "identity": result.identity,
        }
        _write(output / "scientific.json", scientific)
        for pass_index, item in enumerate(native_run["solve"]["passes"]):
            for record in item["trace"]:
                events.append({"iteration": len(events), "status": "ok", "data": {
                    "kind": "SLQP_ITERATION_OBSERVATION", "pass_index": pass_index,
                    "record": record,
                }})
        events.append({"iteration": len(events), "status": "ok", "data": {
            "kind": "EXPECTED_BINARY_VALIDATION_FAILURE_AFTER_SPLIT",
            "numerical_status": result.status,
            "scientific_digest": result.scientific_digest,
        }})
        outputs = sorted(path for path in output.rglob("*")
                         if path.is_file() and path.name != "begin.json")
        evidence = finish_run(evidence_start, events=events, output_files=outputs,
                              measurements={"started_at": started,
                                            "elapsed_seconds": time.perf_counter() - tick})
        _write(output / "evidence.json", evidence.as_dict())
        freshness = check_freshness(evidence)
        summary = {
            "status": result.status, "sealed": result.artifact is not None,
            "classification": "SYNTHETIC_REGRESSION", "project_status": "PARTIAL_BLOCKED",
            "build_source_attestation": "NOT_VERIFIED", "physical_validation": "NOT_IN_SCOPE",
            "original_evidence_admission": "BLOCKED_UNKNOWN",
            "freshness": {"status": freshness.status, "mismatches": list(freshness.mismatches)},
        }
        _write(output / "summary.json", summary)
        return summary
    except Exception as error:
        failure = {"schema": "source-spin2014-adaptive-transition-failure-v1", "status": "FAILED",
                   "error": str(error), "exception_type": type(error).__name__,
                   "traceback": traceback.format_exc()}
        _write(output / "failure.json", failure)
        if evidence_start is not None and not (output / "evidence.json").exists():
            events.append({"iteration": len(events), "status": "failed", "data": failure})
            evidence = finish_run(evidence_start, events=events,
                                  output_files=sorted(path for path in output.rglob("*")
                                                      if path.is_file() and path.name != "begin.json"),
                                  measurements={"started_at": started,
                                                "elapsed_seconds": time.perf_counter() - tick})
            _write(output / "evidence.json", evidence.as_dict())
        raise


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-directory", type=Path, required=True)
    parser.add_argument("--build-directory", type=Path, required=True)
    parser.add_argument("--dependency-root", action="append", required=True, metavar="NAME=PATH")
    args = parser.parse_args()
    roots = {}
    for item in args.dependency_root:
        name, separator, path = item.partition("=")
        if not name or not separator or name in roots:
            parser.error("invalid or duplicate dependency root")
        roots[name] = Path(path).resolve()
    summary = run(args.output_directory, build_directory=args.build_directory, dependency_roots=roots)
    print(canonical(summary))
    return 0 if summary["freshness"]["status"] == "CURRENT" else 1


if __name__ == "__main__":
    raise SystemExit(main())
