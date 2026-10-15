#!/usr/bin/env python3
"""Audit a saved Faster Jacobian against finite differences of real Cherchi outputs.

This is a local fixed-topology diagnostic, not an optimizer or paper-case run.
Every perturbed point is rebuilt by the native source triangulator and the
pinned Cherchi Boolean executable; topology/operation changes are reported as
ineligible differences instead of being silently compared.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import subprocess
import sys

import numpy as np


def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _geometry_from_dual(q: np.ndarray) -> np.ndarray:
    to_dual = np.eye(10)
    to_dual[7:, 7:] *= 2.0
    s = to_dual @ np.asarray(q, dtype=np.float64)
    norm = float(np.linalg.norm(s))
    if not np.isfinite(norm) or norm <= 0.0:
        raise ValueError("quadric has invalid dual-monomial coordinates")
    return s / norm


def dual_tangent_basis(q: object) -> tuple[np.ndarray, np.ndarray]:
    """Return an orthonormal basis of the 9D dual-monomial unit-sphere tangent."""
    sphere = _geometry_from_dual(np.asarray(q, dtype=np.float64))
    _, _, right = np.linalg.svd(sphere.reshape(1, 10), full_matrices=True)
    basis = right[1:, :].T.copy()
    if basis.shape != (10, 9) or not np.all(np.isfinite(basis)):
        raise ValueError("could not construct the dual-sphere tangent basis")
    return basis, sphere


def geometry_from_dual_geodesic(q: object, tangent: object, step: float) -> np.ndarray:
    """Map a dual-sphere geodesic point back to normalized geometry coefficients."""
    if isinstance(step, bool) or not isinstance(step, (int, float, np.floating)) or \
            not np.isfinite(step):
        raise ValueError("geodesic step must be finite")
    s0 = _geometry_from_dual(np.asarray(q, dtype=np.float64))
    direction = np.asarray(tangent, dtype=np.float64)
    if direction.shape != (10,) or not np.all(np.isfinite(direction)):
        raise ValueError("dual tangent must contain ten finite coordinates")
    direction = direction - s0 * float(s0 @ direction)
    tangent_norm = float(np.linalg.norm(direction))
    if tangent_norm <= 0.0:
        raise ValueError("dual geodesic tangent must be nonzero")
    direction /= tangent_norm
    dual = math.cos(float(step)) * s0 + math.sin(float(step)) * direction
    from_dual = np.eye(10)
    from_dual[7:, 7:] *= 0.5
    geometry = from_dual @ dual
    norm = float(np.linalg.norm(geometry))
    if not np.isfinite(norm) or norm <= 0.0:
        raise ValueError("dual geodesic produced invalid geometry coefficients")
    return geometry / norm


def geometry_direction_from_dual_tangent(q: object, tangent: object) -> np.ndarray:
    """Analytic derivative of the geometry-coordinate map at a dual-sphere point."""
    s0 = _geometry_from_dual(np.asarray(q, dtype=np.float64))
    direction = np.asarray(tangent, dtype=np.float64)
    if direction.shape != (10,) or not np.all(np.isfinite(direction)):
        raise ValueError("dual tangent must contain ten finite coordinates")
    direction = direction - s0 * float(s0 @ direction)
    tangent_norm = float(np.linalg.norm(direction))
    if tangent_norm <= 0.0:
        raise ValueError("dual geodesic tangent must be nonzero")
    direction /= tangent_norm
    from_dual = np.eye(10)
    from_dual[7:, 7:] *= 0.5
    unnormalized = from_dual @ s0
    norm = float(np.linalg.norm(unnormalized))
    geometry = unnormalized / norm
    return (np.eye(10) - np.outer(geometry, geometry)) @ (from_dual @ direction) / norm


def _load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ValueError(f"cannot load diagnostic dependency: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def _append_jsonl(path: Path, record: dict[str, object]) -> None:
    with path.open("a", encoding="utf-8") as stream:
        stream.write(json.dumps(record, sort_keys=True, separators=(",", ":"), allow_nan=False))
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())


def same_differential_stratum(base_operation: str, base_topology: str, base_cap: str,
                              operation: str, topology: str, cap: str) -> bool:
    """Require Boolean, construction topology, and source-cap branches to match."""
    return bool(
        base_operation and base_topology and base_cap and
        operation == base_operation and topology == base_topology and cap == base_cap
    )


def cap_branch_identity(payload: object) -> str:
    """Hash cap-generation controls, excluding the moving triangulation geometry hash."""
    if not isinstance(payload, dict):
        raise ValueError("source cap identity payload must be an object")
    keys = (
        "surface_type", "truncation_radius_hex", "initial_icosphere_subdivisions",
        "effective_icosphere_subdivisions", "maximum_radius_doublings",
    )
    if any(key not in payload for key in keys):
        raise ValueError("source cap identity payload is incomplete")
    branch = {key: payload[key] for key in keys}
    if not isinstance(branch["surface_type"], str) or not branch["surface_type"] or \
            not isinstance(branch["truncation_radius_hex"], str) or not branch["truncation_radius_hex"]:
        raise ValueError("source cap type and radius must be explicit")
    for key in keys[2:]:
        value = branch[key]
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise ValueError("source cap discretization controls must be nonnegative integers")
    return hashlib.sha256(json.dumps(
        branch, sort_keys=True, separators=(",", ":"), allow_nan=False,
    ).encode("ascii")).hexdigest()


def _load_parent(run: Path) -> tuple[dict[str, object], dict[str, object], np.ndarray]:
    inputs = json.loads((run / "inputs.json").read_text(encoding="utf-8"))
    result = json.loads((run / "result.json").read_text(encoding="utf-8"))
    freshness = json.loads((run / "freshness.json").read_text(encoding="utf-8"))
    if inputs.get("classification") != "DIAGNOSTIC_REPLAY" or \
            inputs.get("initialization_mode") != "exact_line_search_replay":
        raise ValueError("parent run is not an explicit exact line-search replay")
    if freshness.get("status") != "CURRENT" or freshness.get("mismatches"):
        raise ValueError("parent run evidence is stale or mismatched")
    if result.get("evaluation_map_binding_audit", {}).get("status") != "CONSISTENT":
        raise ValueError("parent evaluation-to-Newton map binding is not consistent")
    if len(result.get("evaluations", [])) != 1:
        raise ValueError("parent replay must contain exactly one actual evaluation")
    evaluation = result["evaluations"][0]
    q = np.asarray(inputs.get("initial_coefficients"), dtype=np.float64)
    jacobian = np.asarray(evaluation.get("moment_map", {}).get("moment_jacobian"), dtype=np.float64)
    if q.shape != (10,) or not np.all(np.isfinite(q)) or abs(float(np.linalg.norm(q)) - 1.0) > 1e-12:
        raise ValueError("parent replay coefficients are invalid")
    if jacobian.shape != (10, 10) or not np.all(np.isfinite(jacobian)):
        raise ValueError("parent discrete moment Jacobian must be finite 10x10")
    return inputs, result, q


def run(*, run_directory: Path, build: Path, probe: Path, output: Path,
        step_sizes: tuple[float, ...] = (1e-5, 1e-6, 1e-7)) -> dict[str, object]:
    if len(step_sizes) < 3 or any(not np.isfinite(h) or h <= 0.0 or h >= 0.5 for h in step_sizes):
        raise ValueError("at least three finite finite-difference steps in (0, 0.5) are required")
    if len(set(step_sizes)) != len(step_sizes):
        raise ValueError("finite-difference step sizes must be distinct")
    run_directory = run_directory.resolve(strict=True)
    build = build.resolve(strict=True)
    probe = probe.resolve(strict=True)
    if not probe.is_file() or not probe.stat().st_mode & 0o111:
        raise ValueError("pinned Cherchi probe is unavailable or non-executable")
    inputs, parent_result, q0 = _load_parent(run_directory)
    if _sha(probe) != inputs.get("cherchi_probe_sha256"):
        raise ValueError("probe digest differs from the parent replay")
    extension_directory = (build / "python").resolve(strict=True)
    extensions = sorted(extension_directory.glob("_make_it_stand*.so"))
    if len(extensions) != 1 or _sha(extensions[0]) != inputs.get("native_extension_sha256"):
        raise ValueError("native extension is absent, ambiguous, or differs from the parent replay")
    output = output.absolute()
    output.mkdir(parents=False, exist_ok=False)
    runner_path = Path(__file__).with_name("run_faster_cherchi_newton_integration.py").resolve()
    sys.path.insert(0, str(runner_path.parent))
    runner = _load_module(runner_path, "faster_newton_runner_for_fd_audit")
    sys.path.insert(0, str(extension_directory))
    try:
        extension = __import__("_make_it_stand")
    finally:
        sys.path.pop(0)
    if Path(extension.__file__).resolve() != extensions[0]:
        raise ValueError("loaded native extension differs from declared build")
    import make_it_stand as spin

    design_path = run_directory / "design.obj"
    design_vertices, _ = runner._read_obj(design_path)
    triangulation, _ = runner._make_triangulation_configuration(
        extension,
        icosphere_subdivisions=int(inputs["quadric_icosphere_subdivisions"]),
        maximum_icosphere_subdivisions=int(inputs["maximum_icosphere_subdivisions"]),
    )
    base_evaluation = parent_result["evaluations"][0]
    base_operation = base_evaluation["operation"]
    base_topology = base_evaluation["construction_topology_hash"]
    base_cap_payload = base_evaluation.get("cap_identity_payload")
    base_cap_branch = cap_branch_identity(base_cap_payload)
    jacobian = np.asarray(base_evaluation["moment_map"]["moment_jacobian"], dtype=np.float64)
    basis, _ = dual_tangent_basis(q0)
    records_path = output / "finite-difference-evaluations.jsonl"
    records_path.touch(exist_ok=False)
    directions: list[dict[str, object]] = []

    def evaluate_point(coefficients: np.ndarray, direction_index: int, h: float, sign: int):
        source = extension.evaluate_spin_it_faster_source_motion(
            coefficients, design_vertices, triangulation,
        )
        if not source["defined"]:
            raise ValueError("source motion undefined: " + str(source["diagnostic"]))
        source_vertices = np.asarray(source["vertices"], dtype=np.float64)
        source_faces = np.asarray(source["faces"], dtype=np.int64)
        selector = runner._source_side_selection_details(source)
        operation = str(selector["operation"])
        source_path = output / f"source-d{direction_index:02d}-h{h:.0e}-s{sign:+d}.obj"
        material_path = output / f"material-d{direction_index:02d}-h{h:.0e}-s{sign:+d}.obj"
        runner._write_obj(source_path, source_vertices, source_faces)
        completed = subprocess.run(
            [str(probe), operation, str(design_path), str(source_path), str(material_path),
             "--dump-exact-rational"],
            cwd=output, capture_output=True, text=True,
            timeout=float(inputs["cherchi_boolean_timeout_seconds"]), check=False,
        )
        if completed.returncode != 0:
            raise ValueError("Cherchi Boolean failed: " + completed.stderr.strip())
        material_vertices, material_faces = runner._read_obj(material_path)
        oriented_faces, orientation_audit = runner.normalize_closed_component_orientation(
            material_vertices, material_faces,
        )
        moments = np.asarray(runner._independent_polyhedral_raw_moments(
            material_vertices, oriented_faces, 1.0,
        ), dtype=np.float64)
        construction = runner._summarize_vertex_construction_provenance(
            (design_path, source_path), material_path,
            Path(str(material_path) + ".construction.tsv"),
        )
        topology = runner._construction_topology_signature(construction, oriented_faces)
        cap_payload = {
            "surface_type": str(source["surface_type"]),
            "truncation_radius_hex": float(source["truncation_radius"]).hex(),
            "initial_icosphere_subdivisions": int(source["initial_icosphere_subdivisions"]),
            "effective_icosphere_subdivisions": int(source["effective_icosphere_subdivisions"]),
            "maximum_radius_doublings": int(triangulation.maximum_radius_doublings),
            "triangulation_hash": str(source["triangulation_hash"]),
        }
        cap_identity = hashlib.sha256(json.dumps(
            cap_payload, sort_keys=True, separators=(",", ":"), allow_nan=False,
        ).encode("ascii")).hexdigest()
        cap_branch = cap_branch_identity(cap_payload)
        row = {
            "direction_index": direction_index,
            "step": h,
            "sign": sign,
            "quadric_coefficients": coefficients.tolist(),
            "operation": operation,
            "selector_margin": selector["normalized_signed_margin"],
            "topology_hash": topology,
            "cap_identity": cap_identity,
            "cap_branch_identity": cap_branch,
            "cap_identity_payload": cap_payload,
            "face_count": int(len(oriented_faces)),
            "orientation_component_count": orientation_audit["component_count"],
            "raw_moments": moments.tolist(),
            "raw_moments_sha256": hashlib.sha256(np.asarray(moments, dtype="<f8").tobytes()).hexdigest(),
        }
        _append_jsonl(records_path, row)
        return moments, operation, topology, cap_identity, cap_branch

    for direction_index in range(9):
        dual_direction = basis[:, direction_index]
        geometry_direction = geometry_direction_from_dual_tangent(q0, dual_direction)
        analytic_derivative = jacobian @ geometry_direction
        step_records = []
        for h in step_sizes:
            pair = [evaluate_point(
                geometry_from_dual_geodesic(q0, dual_direction, sign * h),
                direction_index, h, sign,
            ) for sign in (-1, 1)]
            minus, plus = pair
            topology_stable = same_differential_stratum(
                base_operation, base_topology, base_cap_branch,
                minus[1], minus[2], minus[4],
            ) and same_differential_stratum(
                base_operation, base_topology, base_cap_branch,
                plus[1], plus[2], plus[4],
            )
            finite_difference = (plus[0] - minus[0]) / (2.0 * h)
            error = float(np.linalg.norm(finite_difference - analytic_derivative) /
                          max(1.0, float(np.linalg.norm(finite_difference))))
            step_records.append({
                "step": h,
                "operation_pair": [minus[1], plus[1]],
                "topology_pair_equal": minus[2] == plus[2],
                "topology_matches_base": minus[2] == plus[2] == base_topology,
                "cap_identity_pair_equal": minus[3] == plus[3],
                "cap_identity_matches_base": minus[3] == plus[3] == base_evaluation["cap_identity"],
                "cap_branch_pair_equal": minus[4] == plus[4],
                "cap_branch_matches_base": minus[4] == plus[4] == base_cap_branch,
                "eligible_fixed_topology_difference": topology_stable,
                "finite_difference": finite_difference.tolist(),
                "analytic": analytic_derivative.tolist(),
                "relative_error": error if topology_stable else None,
            })
        eligible = [item for item in step_records if item["eligible_fixed_topology_difference"]]
        accurate = [item for item in eligible if item["relative_error"] <= 1e-6]
        directions.append({
            "direction_index": direction_index,
            "dual_tangent": dual_direction.tolist(),
            "geometry_coefficient_velocity": geometry_direction.tolist(),
            "analytic_moment_velocity": analytic_derivative.tolist(),
            "steps": step_records,
            "eligible_step_count": len(eligible),
            "accurate_step_count": len(accurate),
            "best_eligible_relative_error": min(
                (float(item["relative_error"]) for item in eligible), default=None,
            ),
        })

    passed = all(
        item["eligible_step_count"] >= 2 and item["accurate_step_count"] >= 1
        for item in directions
    )
    summary = {
        "schema": "faster-actual-cherchi-jacobian-audit-v1",
        "classification": "METHOD_VALIDATION",
        "status": "PASSED" if passed else "INCOMPLETE_OR_FAILED",
        "scope": "actual Cherchi Boolean-to-ten-moment differential on one fixed topology stratum and a 9D dual-sphere tangent basis",
        "parent_run": str(run_directory),
        "parent_result_sha256": _sha(run_directory / "result.json"),
        "parent_evidence_sha256": _sha(run_directory / "evidence.json"),
        "native_extension_sha256": _sha(extensions[0]),
        "cherchi_probe_sha256": _sha(probe),
        "audit_script_sha256": _sha(Path(__file__).resolve()),
        "quadric_coefficients": q0.tolist(),
        "boolean_operation": base_operation,
        "construction_topology_hash": base_topology,
        "cap_identity": base_evaluation["cap_identity"],
        "cap_branch_identity": base_cap_branch,
        "cap_identity_payload": base_cap_payload,
        "finite_difference_step_sizes": list(step_sizes),
        "required_relative_error": 1e-6,
        "required_eligible_steps_per_direction": 2,
        "directions": directions,
        "evaluation_ledger": {
            "path": records_path.name,
            "rows": sum(len(item["steps"]) * 2 for item in directions),
            "sha256": _sha(records_path),
        },
    }
    summary_path = output / "summary.json"
    summary_path.write_text(
        json.dumps(summary, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return summary


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True, dest="run_directory")
    parser.add_argument("--build", type=Path, required=True)
    parser.add_argument("--probe", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--step-sizes", type=float, nargs="+", default=(1e-5, 1e-6, 1e-7))
    return parser


def main() -> None:
    args = _parser().parse_args()
    summary = run(
        run_directory=args.run_directory,
        build=args.build,
        probe=args.probe,
        output=args.output,
        step_sizes=tuple(args.step_sizes),
    )
    print(json.dumps({
        "status": summary["status"],
        "scope": summary["scope"],
        "directions": len(summary["directions"]),
        "evaluations": summary["evaluation_ledger"]["rows"],
    }, sort_keys=True))
    if summary["status"] != "PASSED":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
