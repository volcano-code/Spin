#!/usr/bin/env python3
"""Inspect a source-only Medit surface with the Spin-It source octree.

This read-only diagnostic does not optimize, alter geometry, or claim a paper
case. It is intended for audited public source models whose Spin-It pose and
design parameters are unknown.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
from pathlib import Path
import time

from make_it_stand.verification.medit_surface import read_medit_surface
from make_it_stand.verification.run_evidence import (
    EXPLICIT_EFFECTIVE_SOURCE_POLICY, begin_run, check_freshness, finish_run,
)
from make_it_stand.verification.source_spin_material import SourceGeometrySettings
from make_it_stand.verification.spin_common import canonical, native, parse


ROOT = Path(__file__).resolve().parents[3]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _write(path: Path, value) -> None:
    with path.open("x", encoding="ascii") as stream:
        stream.write(canonical(value) + "\n")


def _run_in_output(output: Path, *, mesh_path: Path, expected_sha256: str,
                   build_directory: Path, initial_level: int, maximum_depth: int,
                   maximum_leaf_count: int, dependency_roots: dict[str, Path]) -> dict:
    output = output.absolute()
    mesh_path = mesh_path.resolve(strict=True)
    build = build_directory.resolve(strict=True)
    observed_sha256 = _sha256(mesh_path)
    if observed_sha256 != expected_sha256.lower():
        raise ValueError("source model SHA-256 does not match the locked input")
    surface = read_medit_surface(mesh_path)
    vertices, faces = surface.vertices, surface.faces
    lo, hi = vertices.min(axis=0), vertices.max(axis=0)
    side = float(max(hi - lo))
    source_settings = SourceGeometrySettings(
        root_origin=tuple(float(value) for value in lo),
        root_side=float(__import__("math").nextafter(side, float("inf"))),
        maximum_depth=maximum_depth,
    )
    if not 0 <= initial_level <= maximum_depth <= 12:
        raise ValueError("diagnostic depth must satisfy 0 <= initial <= maximum <= 12")

    inputs = parse(canonical({
        "schema": "source-model-domain-inspection-inputs-v1",
        "classification": "SOURCE_MODEL_ONLY",
        "optimizer_invoked": False,
        "original_paper_case_credit": False,
        "asset": {"path": str(mesh_path), "sha256": observed_sha256,
                  "medit_vertex_count": surface.medit_vertex_count,
                  "surface_vertex_count": surface.surface_vertex_count,
                  "triangle_count": surface.triangle_count,
                  "excluded_unreferenced_vertex_count": surface.excluded_unreferenced_vertex_count},
        "surface_ingestion": "MEDIT_TRIANGLES_UNCHANGED_ORDER_AND_COORDINATES; UNUSED_VOLUME_VERTICES_EXCLUDED",
        "source_settings": asdict(source_settings),
        "initial_level": initial_level,
        "maximum_optimization_level": maximum_depth,
        "maximum_leaf_count": maximum_leaf_count,
        "diagnostic_limitation": "DESIGN_SURFACE_EQUALS_SOURCE_SURFACE; READ_ONLY_DOMAIN_INSPECTION_ONLY",
    }))
    input_path = output / "inputs.json"
    _write(input_path, inputs)

    extension = native()
    if Path(extension.__file__).resolve().parents[1] != build:
        raise ValueError("declared build directory differs from loaded native extension")
    started, tick, events = datetime.now(timezone.utc).isoformat(), time.perf_counter(), []
    evidence_start = begin_run(
        repository_root=ROOT,
        runner_files=[Path(__file__)],
        checker_files=[ROOT / "python/make_it_stand/verification/medit_surface.py",
                       ROOT / "python/make_it_stand/verification/source_spin_material.py",
                       ROOT / "bindings/SourceSpinMaterialBindings.cpp",
                       ROOT / "src/spin2014/SourceSurfaceOctree.cpp"],
        asset_files=[mesh_path, input_path],
        build_configuration_files=[ROOT / "CMakeLists.txt", build / "CMakeCache.txt"],
        parameters=inputs,
        seed_protocol={"seed": 0, "protocol": "deterministic read-only source-domain construction"},
        dependency_roots=dependency_roots, native_modules=[extension], pure_python=False,
        source_policy=EXPLICIT_EFFECTIVE_SOURCE_POLICY,
    )
    _write(output / "begin.json", evidence_start.as_dict())
    observation = extension.inspect_source_spin2014_domain(
        vertices, faces, vertices, faces,
        source_settings=asdict(source_settings), initial_level=initial_level,
        maximum_optimization_level=maximum_depth,
        maximum_leaf_count=maximum_leaf_count, include_details=False,
    )
    stats = observation["source"]["stats"]
    if stats.get("maximum_visited_level") != maximum_depth:
        raise ValueError("source octree did not visit the requested maximum depth")
    if stats.get("boundary_nodes", 0) <= 0 or observation.get("leaf_count", 0) <= 0:
        raise ValueError("source-domain inspection lacks boundary nodes or adaptive leaves")
    observation_path = output / "domain-inspection.json"
    _write(observation_path, observation)
    events.append({"iteration": 0, "status": "ok", "data": {
        "kind": "SOURCE_MODEL_ONLY_READ_ONLY_DOMAIN_INSPECTION",
        "source_hash": observation["source"]["source_hash"],
        "partition_hash": observation["partition_hash"],
        "maximum_visited_level": stats["maximum_visited_level"],
        "boundary_nodes": stats["boundary_nodes"],
        "leaf_count": observation["leaf_count"],
        "optimizer_invoked": False,
    }})
    evidence = finish_run(
        evidence_start, events=events,
        output_files=[input_path, observation_path],
        measurements={"started_at": started, "elapsed_seconds": time.perf_counter() - tick},
    )
    _write(output / "evidence.json", evidence.as_dict())
    freshness = check_freshness(evidence)
    summary = {
        "schema": "source-model-domain-inspection-summary-v1",
        "status": "SOURCE_DOMAIN_INSPECTED" if freshness.status == "CURRENT" else "STALE_RESULT",
        "classification": "SOURCE_MODEL_ONLY",
        "optimizer_invoked": False,
        "original_paper_case_credit": False,
        "project_status": "PARTIAL_BLOCKED",
        "asset_sha256": observed_sha256,
        "source_hash": observation["source"]["source_hash"],
        "partition_hash": observation["partition_hash"],
        "source_stats": stats,
        "adaptive_domain": {"leaf_count": observation["leaf_count"],
                            "edge_count": observation["edge_count"],
                            "leaf_level_counts": observation["leaf_level_counts"]},
        "freshness": {"status": freshness.status, "mismatches": list(freshness.mismatches)},
    }
    _write(output / "summary.json", summary)
    return summary


def run(output_directory: Path, *, mesh_path: Path, expected_sha256: str,
        build_directory: Path, initial_level: int, maximum_depth: int,
        maximum_leaf_count: int, dependency_roots: dict[str, Path]) -> dict:
    output = output_directory.absolute()
    output.mkdir(parents=False, exist_ok=False)
    try:
        return _run_in_output(
            output, mesh_path=mesh_path, expected_sha256=expected_sha256,
            build_directory=build_directory, initial_level=initial_level,
            maximum_depth=maximum_depth, maximum_leaf_count=maximum_leaf_count,
            dependency_roots=dependency_roots,
        )
    except Exception as error:
        failure_path = output / "failure.json"
        if not failure_path.exists():
            _write(failure_path, {
                "schema": "source-model-domain-inspection-failure-v1",
                "status": "FAILED", "stage": "RUNNER", "error": str(error),
            })
        raise


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--expected-sha256", required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    parser.add_argument("--build-directory", type=Path, required=True)
    parser.add_argument("--initial-level", type=int, default=1)
    parser.add_argument("--maximum-depth", type=int, required=True)
    parser.add_argument("--maximum-leaf-count", type=int, default=200000)
    parser.add_argument("--dependency-root", action="append", required=True, metavar="NAME=PATH")
    args = parser.parse_args()
    roots = {}
    for item in args.dependency_root:
        name, separator, value = item.partition("=")
        if not name or not separator or name in roots:
            parser.error("invalid or duplicate dependency root")
        roots[name] = Path(value).resolve()
    result = run(args.output_directory, mesh_path=args.mesh,
                 expected_sha256=args.expected_sha256, build_directory=args.build_directory,
                 initial_level=args.initial_level, maximum_depth=args.maximum_depth,
                 maximum_leaf_count=args.maximum_leaf_count, dependency_roots=roots)
    print(canonical(result))
    return 0 if result["freshness"]["status"] == "CURRENT" else 1


if __name__ == "__main__":
    raise SystemExit(main())
