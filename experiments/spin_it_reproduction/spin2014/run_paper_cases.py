#!/usr/bin/env python3
"""Run every public Spin-It 2014 comparison on reconstructed voxel assets.

The author meshes, interaction history and fabrication files are unavailable.
Consequently this runner measures the published equations and experiment roles
on deterministic structural counterparts and never labels them author-exact.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import time
from typing import Any, Iterable

import numpy as np


ROOT = Path(__file__).resolve().parent
DEFAULT_CASES = ROOT / "paper_cases.json"
TIMING_KEYS = {"elapsed_seconds"}


def _without_timings(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: _without_timings(item) for key, item in value.items() if key not in TIMING_KEYS}
    if isinstance(value, list):
        return [_without_timings(item) for item in value]
    return value


def scientific_hash(payload: Any) -> str:
    encoded = json.dumps(
        _without_timings(payload), sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode()
    return hashlib.sha256(encoded).hexdigest()


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def decode_rle(runs: Iterable[Iterable[int]], expected: int) -> tuple[int, ...]:
    values: list[int] = []
    for row in runs:
        pair = tuple(row)
        if len(pair) != 2:
            raise ValueError("RLE rows must be value/count pairs")
        value, count = pair
        if value not in (0, 1) or isinstance(count, bool) or not isinstance(count, int) or count <= 0:
            raise ValueError("RLE rows must contain a binary value and positive count")
        values.extend([int(value)] * count)
    if len(values) != expected:
        raise ValueError(f"RLE expands to {len(values)} values, expected {expected}")
    return tuple(values)


def count_components(mask: Iterable[int], dimensions: tuple[int, int, int]) -> int:
    values = tuple(int(value) for value in mask)
    nx, ny, nz = dimensions
    if len(values) != nx * ny * nz:
        raise ValueError("component mask size does not match dimensions")
    remaining = {index for index, value in enumerate(values) if value}
    count = 0
    while remaining:
        count += 1
        stack = [remaining.pop()]
        while stack:
            index = stack.pop()
            x = index % nx
            y = (index // nx) % ny
            z = index // (nx * ny)
            neighbors = []
            if x: neighbors.append(index - 1)
            if x + 1 < nx: neighbors.append(index + 1)
            if y: neighbors.append(index - nx)
            if y + 1 < ny: neighbors.append(index + nx)
            if z: neighbors.append(index - nx * ny)
            if z + 1 < nz: neighbors.append(index + nx * ny)
            for neighbor in neighbors:
                if neighbor in remaining:
                    remaining.remove(neighbor)
                    stack.append(neighbor)
    return count


def _central_inertia(raw: Iterable[float]) -> np.ndarray:
    r = np.asarray(tuple(raw), dtype=float)
    if r.shape != (10,) or not np.isfinite(r).all() or r[0] <= 0.0:
        raise ValueError("raw moments must contain ten finite values with positive mass")
    mass = r[0]
    first = r[1:4]
    second = np.array(((r[4], r[9], r[8]), (r[9], r[5], r[7]), (r[8], r[7], r[6])))
    inertia = np.trace(second) * np.eye(3) - second
    center = first / mass
    return inertia - mass * ((center @ center) * np.eye(3) - np.outer(center, center))


def paper_metrics(mode: str, raw: Iterable[float], center_weight: float) -> dict[str, float]:
    r = np.asarray(tuple(raw), dtype=float)
    inertia = _central_inertia(r)
    lateral = inertia[:2, :2]
    axial = float(inertia[2, 2])
    yoyo = float(np.sum(lateral * lateral) / (axial * axial))
    if mode == "YO_YO":
        objective = yoyo
    elif mode == "SPINNING_TOP":
        objective = yoyo + float(center_weight) * float(r[3] * r[3])
    elif mode == "STATIC_BALANCE":
        objective = float(center_weight) * float((r[3] / r[0]) ** 2)
    else:
        raise ValueError(f"unsupported Spin-It 2014 mode: {mode}")
    eigenvalues, eigenvectors = np.linalg.eigh(inertia)
    # The equality constraints make the requested z axis a principal axis.
    # Identify that eigenpair by direction, then separately audit whether its
    # eigenvalue is the largest (paper Section 3.2 requirements 3 and 4).
    axial_index = int(np.argmax(np.abs(eigenvectors[2, :])))
    axial_vector = eigenvectors[:, axial_index]
    cosine = min(1.0, max(0.0, abs(float(axial_vector[2]))))
    lateral_eigenvalues = np.linalg.eigvalsh(lateral)
    if mode == "STATIC_BALANCE":
        constraint_residual = max(abs(r[1]), abs(r[2]))
    elif mode == "YO_YO":
        constraint_residual = max(
            abs(r[1]), abs(r[2]), abs(r[3]),
            abs(inertia[0, 2]), abs(inertia[1, 2]),
        )
    else:
        constraint_residual = max(
            abs(r[1]), abs(r[2]),
            abs(inertia[0, 2]), abs(inertia[1, 2]),
        )
    return {
        "mass": float(r[0]),
        "center_height": float(r[3] / r[0]),
        "objective": objective,
        "axis_alignment_error": math.degrees(math.acos(cosine)),
        "axial_is_maximum_principal_moment": bool(
            eigenvalues[axial_index] >= float(np.max(eigenvalues)) - 1e-10
        ),
        "principal_moment_ratio": axial / max(float(np.max(lateral_eigenvalues)), 1e-30),
        "constraint_infinity_norm": float(constraint_residual),
    }


def load_asset(case: dict[str, Any]) -> dict[str, Any]:
    path = ROOT / case["occupancy"]["path"]
    if file_sha256(path) != case["occupancy"]["sha256"]:
        raise RuntimeError(f"{case['id']}: occupancy hash mismatch")
    asset = json.loads(path.read_text(encoding="utf-8"))
    dimensions = tuple(int(value) for value in asset["dimensions"])
    expected = math.prod(dimensions)
    return {
        "dimensions": dimensions,
        "minimum": tuple(float(value) for value in asset["minimum"]),
        "maximum": tuple(float(value) for value in asset["maximum"]),
        "material": decode_rle(asset["material_rle"], expected),
        "carveable": decode_rle(asset["carveable_rle"], expected),
        "sha256": case["occupancy"]["sha256"],
    }


def _mesh_component_count(vertices: np.ndarray, faces: np.ndarray) -> int:
    if len(faces) == 0:
        return 0
    parent = list(range(len(vertices)))
    def find(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index
    def join(first: int, second: int) -> None:
        a, b = find(first), find(second)
        if a != b:
            parent[b] = a
    used: set[int] = set()
    for face in np.asarray(faces, dtype=int):
        a, b, c = map(int, face)
        used.update((a, b, c))
        join(a, b); join(b, c)
    return len({find(index) for index in used})


def _incline_oracle(result: Any) -> dict[str, Any] | None:
    vertices = np.asarray(result.outer_vertices, dtype=float)
    center = np.asarray(result.diagnostics["center_of_mass_model"], dtype=float)
    if len(vertices) == 0 or not np.isfinite(center).all():
        return None
    zmin = float(np.min(vertices[:, 2]))
    height_scale = max(float(np.ptp(vertices[:, 2])), 1e-12)
    contact = vertices[vertices[:, 2] <= zmin + 1e-6 * height_scale]
    if len(contact) == 0 or center[2] <= zmin:
        return None
    margins = np.max(np.abs(contact[:, :2] - center[:2]), axis=0)
    critical = math.degrees(math.atan2(float(np.min(margins)), float(center[2] - zmin)))
    return {
        "protocol": "quasi_static_support_extent_oracle",
        "critical_angle_degrees": critical,
        "paper_reported_spin_it_degrees": 8.0,
        "paper_reported_prevost_degrees": 1.0,
        "author_geometry_comparison": False,
    }


def run_case(
    case: dict[str, Any],
    repeats: int = 2,
    maximum_level: int | None = None,
    bundle_root: Path | None = None,
) -> dict[str, Any]:
    from make_it_stand.spin2014 import (
        SpinIt2014Request,
        solve_spin_it_2014,
        solve_spin_it_2014_cage,
    )

    asset = load_asset(case)
    density = case.get("densities_g_cm3", {})
    effective_level = int(
        maximum_level if maximum_level is not None else case["maximum_octree_level"]
    )
    effective_binary_tolerance = float(
        case["paper_binary_constraint_tolerance"]
        if effective_level >= int(case["paper_maximum_octree_level"])
        else case["binary_constraint_tolerance"]
    )
    effective_laplacian_weight = float(
        case["paper_level_laplacian_weight"]
        if effective_level >= int(case["paper_maximum_octree_level"])
        else case["laplacian_weight"]
    )
    request = SpinIt2014Request(
        dimensions=asset["dimensions"],
        minimum=asset["minimum"],
        maximum=asset["maximum"],
        material_mask=asset["material"],
        carveable_mask=asset["carveable"],
        mode=case["mode"],
        evidence="RECONSTRUCTED_ASSET",
        initial_level=2,
        maximum_level=effective_level,
        maximum_adaptive_passes=max(8, effective_level - 2 + 2),
        maximum_leaf_count=1_000_000,
        initial_fill=float(case.get("initial_fill", 0.5)),
        laplacian_weight=effective_laplacian_weight,
        center_of_mass_weight=float(case["center_of_mass_weight"]),
        boundary_density=float(density.get("shell", 1.0)),
        interior_density=float(density.get("insert", 1.0)),
        binary_constraint_tolerance=effective_binary_tolerance,
        slqp_maximum_iterations=600 if effective_level >= 9 else 240,
    )
    rows: list[dict[str, Any]] = []
    for repeat in range(1, repeats + 1):
        begin = time.perf_counter()
        result = solve_spin_it_2014(request)
        elapsed = time.perf_counter() - begin
        diagnostics = result.diagnostics
        metrics = (
            paper_metrics(case["mode"], diagnostics["final_raw_moments"], case["center_of_mass_weight"])
            if diagnostics.get("mass_properties_available")
            else None
        )
        row = {
            "repeat": repeat,
            "status": result.status,
            "geometry_sealed": result.geometry_sealed,
            "input_hash": diagnostics["input_hash"],
            "trace_hash": diagnostics["trace_hash"],
            "geometry_hash": diagnostics["geometry_hash"],
            "output_hash": diagnostics["output_hash"],
            "binary_constraint_infinity_norm": diagnostics["binary_constraint_infinity_norm"],
            "cavity_components": _mesh_component_count(result.cavity_vertices, result.cavity_faces),
            "metrics": metrics,
            "elapsed_seconds": elapsed,
            "cage_output_hash": "",
        }
        if bundle_root is not None:
            paths = result.save_viewer_bundle(
                bundle_root / case["id"] / f"repeat-{repeat}"
            )
            row["viewer_bundle"] = {
                Path(path).name: file_sha256(Path(path)) for path in paths
            }
        if case.get("cage", {}).get("deformation_enabled", False):
            cage_begin = time.perf_counter()
            try:
                cage = solve_spin_it_2014_cage(
                    result,
                    request,
                    maximum_handle_count=int(case["cage"].get("maximum_handle_count", 12)),
                )
            except (RuntimeError, ValueError) as error:
                cage = {
                    "status": "EXECUTION_FAILED",
                    "diagnostic": str(error),
                    "geometry_valid": False,
                    "silhouette_changed": False,
                    "trace_hash": "",
                    "output_hash": "",
                    "binding_hash": "",
                }
            cage["elapsed_seconds"] = time.perf_counter() - cage_begin
            cage.pop("deformed_outer_vertices", None)
            cage.pop("deformed_cavity_vertices", None)
            row["cage_deformation"] = cage
            row["cage_output_hash"] = cage.get("output_hash", "")
        if case["id"] == "trex-static":
            row["incline"] = _incline_oracle(result)
        rows.append(row)
    deterministic_fields = ("status", "geometry_sealed", "input_hash", "trace_hash", "geometry_hash", "output_hash", "cavity_components", "metrics", "cage_output_hash")
    deterministic = all(rows[0][field] == row[field] for row in rows[1:] for field in deterministic_fields)
    passed = bool(rows) and deterministic and all(
        row["status"] == "CONVERGED"
        and row["geometry_sealed"]
        and row["binary_constraint_infinity_norm"] <= effective_binary_tolerance
        and row["metrics"] is not None
        and row["metrics"]["constraint_infinity_norm"] <= effective_binary_tolerance
        and (
            case["mode"] == "STATIC_BALANCE"
            or row["metrics"]["axial_is_maximum_principal_moment"] is True
        )
        and (
            not case.get("cage", {}).get("deformation_enabled", False)
            or (
                row.get("cage_deformation", {}).get("status") == "CONVERGED"
                and row.get("cage_deformation", {}).get("geometry_valid") is True
                and row.get("cage_deformation", {}).get("silhouette_changed") is True
            )
        )
        for row in rows
    )
    payload = {
        "id": case["id"],
        "paper_role": case["paper_role"],
        "mode": case["mode"],
        "asset_class": "RECONSTRUCTED_ASSET",
        "author_exact": False,
        "physical_validation": False,
        "maximum_octree_level": effective_level,
        "paper_maximum_octree_level": int(case["paper_maximum_octree_level"]),
        "paper_resolution_reproduced": effective_level >= int(case["paper_maximum_octree_level"]),
        "binary_constraint_tolerance": effective_binary_tolerance,
        "paper_binary_constraint_tolerance": float(case["paper_binary_constraint_tolerance"]),
        "paper_numerics_reproduced": effective_binary_tolerance <= float(case["paper_binary_constraint_tolerance"]),
        "laplacian_weight": effective_laplacian_weight,
        "laplacian_weight_evidence": case["laplacian_weight_evidence"],
        "initial_fill": float(case.get("initial_fill", 0.5)),
        "initial_fill_evidence": case.get("initial_fill_evidence", "PAPER_DEFAULT"),
        "occupancy_sha256": asset["sha256"],
        "paper_reported_void_components": case["expected_void_components"],
        "paper_reported_manual_cuts": case["manual_planar_cuts"],
        "automatic_cut_optimization": False,
        "mass_reduction_extension": bool(case.get("mass_reduction_extension", False)),
        "cage_deformation_requested": bool(case.get("cage", {}).get("deformation_enabled", False)),
        "cage_deformation_executed": bool(rows) and all(
            row.get("cage_deformation", {}).get("status") == "CONVERGED"
            for row in rows
        ) if bool(case.get("cage", {}).get("deformation_enabled", False)) else True,
        "deterministic": deterministic,
        "passed": passed,
        "repeats": rows,
    }
    payload["scientific_hash"] = scientific_hash(payload)
    return payload


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=DEFAULT_CASES)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--case", action="append", default=[])
    parser.add_argument("--repeats", type=int, default=2)
    parser.add_argument("--maximum-level", type=int)
    parser.add_argument("--bundle-root", type=Path)
    args = parser.parse_args()
    if args.repeats <= 0:
        raise SystemExit("repeats must be positive")
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    cases = [row for row in manifest["cases"] if not args.case or row["id"] in args.case]
    if args.case and len(cases) != len(set(args.case)):
        raise SystemExit("one or more requested case IDs are absent")
    results = [
        run_case(case, args.repeats, args.maximum_level, args.bundle_root)
        for case in cases
    ]
    payload = {
        "schema": "spin-it-2014-paper-results-v1",
        "schema_version": 1,
        "paper_doi": "10.1145/2601097.2601157",
        "asset_class": "RECONSTRUCTED_ASSET",
        "author_exact": False,
        "physical_validation": False,
        "paper_resolution_reproduced": bool(results) and all(row["paper_resolution_reproduced"] for row in results),
        "paper_numerics_reproduced": bool(results) and all(row["paper_numerics_reproduced"] for row in results),
        "cage_deformation_executed": bool(results) and all(
            row["cage_deformation_executed"] for row in results
        ),
        "all_public_computational_cases_passed": all(row["passed"] for row in results),
        "cases": results,
    }
    payload["scientific_hash"] = scientific_hash(payload)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0 if payload["all_public_computational_cases_passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
