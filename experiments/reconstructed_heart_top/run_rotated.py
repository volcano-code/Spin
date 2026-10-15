"""Run the reconstructed heart with a frozen non-principal initial orientation.

This is an algorithm validation case, not an author asset.  The outer surface,
fixed shell, design region, target axis, and rotation are all recorded before
the solver starts; no result-dependent input selection is performed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import time
from dataclasses import asdict
from pathlib import Path

import numpy as np
import trimesh

from heart_top import heart_mesh, integrate
from run import write


def _rotate_about_center(vertices: np.ndarray, center: np.ndarray, degrees: float) -> np.ndarray:
    angle = np.deg2rad(float(degrees))
    c, s = float(np.cos(angle)), float(np.sin(angle))
    rotation = np.array(((c, 0.0, s), (0.0, 1.0, 0.0), (-s, 0.0, c)))
    points = np.asarray(vertices, dtype=float)
    center = np.asarray(center, dtype=float)
    return (points - center) @ rotation.T + center


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--rotation-degrees", type=float, default=25.0)
    parser.add_argument("--initial-level", type=int, default=3)
    parser.add_argument("--maximum-level", type=int, default=4)
    args = parser.parse_args()
    output = args.output
    output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()

    from make_it_stand.verification import source_spin_material as bridge
    from make_it_stand.verification.spin_common import native
    from make_it_stand.verification.spin_problem import SpinArtifactProblem
    from make_it_stand.verification.spin_bundle import save_spin_artifact

    original_vertices, faces = heart_mesh()
    original_moments = integrate(original_vertices, faces)
    rotation_center = original_moments[1:4] / original_moments[0]
    original_vertices = _rotate_about_center(
        original_vertices, rotation_center, args.rotation_degrees
    ) - rotation_center
    # The source solver currently accepts the native target frame only:
    # target origin zero and target axis +Z. Translating by the input centroid
    # makes that fixed frame pass through the centroid while retaining the
    # deliberately rotated principal axes.
    center = np.zeros(3, dtype=float)
    design_vertices = center + 0.72 * (original_vertices - center)
    fixed_shell = (
        np.vstack((original_vertices, design_vertices)),
        np.vstack((faces, faces[:, ::-1] + len(original_vertices))),
    )
    for name, mesh in (
        ("outer", (original_vertices, faces)),
        ("design", (design_vertices, faces)),
        ("fixed-shell", fixed_shell),
    ):
        trimesh.Trimesh(*mesh, process=False).export(output / f"{name}.obj")

    geometry = bridge.SourceGeometrySettings(
        root_origin=(-2.0, -2.0, -2.0), root_side=4.0, maximum_depth=5
    )
    solver = bridge.SourceSolverSettings(
        objective="TOP_MASS_AND_LOWERING",
        initial_level=int(args.initial_level),
        maximum_optimization_level=int(args.maximum_level),
        maximum_adaptive_passes=10,
        maximum_leaf_count=4000,
        initial_fill=0.5,
        initial_phi=0.2,
        laplacian_weight=0.001,
        slqp=bridge.AuditedSlqpSettings(maximum_iterations=160),
    )
    boundary = bridge.SourceBoundarySettings(maximum_input_cuboids=4000)
    config = {
        "asset": "RECONSTRUCTED_SELF_BUILT_HEART_ROTATED_V1",
        "author_original": False,
        "rotation_about_y_degrees": float(args.rotation_degrees),
        "rotation_center_input_model": rotation_center.tolist(),
        "translation_to_target_frame": (-rotation_center).tolist(),
        "initial_level": int(args.initial_level),
        "maximum_optimization_level": int(args.maximum_level),
        "target_axis_model": [0.0, 0.0, 1.0],
        "source_center_model": center.tolist(),
        "source": asdict(geometry),
        "solver": asdict(solver),
        "boundary": asdict(boundary),
        "minimum_outer_wall": 0.015,
        "solver_density": 1.0,
        "meters_per_model_unit": 0.025,
        "physical_density_kg_per_m3": 1200.0,
        "protocol": "NON_PRINCIPAL_INITIAL_AXIS_ALGORITHM_VALIDATION",
    }
    write(output / "config.json", config)
    root = Path(__file__).resolve().parents[2]
    identity_paths = [
        Path(__file__), Path(__file__).with_name("heart_top.py"), Path(native().__file__),
        root / "python/make_it_stand/verification/source_spin_material.py",
    ]
    identities = {
        str(path.resolve()): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in identity_paths
    }
    write(
        output / "source.json",
        {
            "head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip(),
            "files": identities,
            "scope": "LOADED_BINARY_AND_DECLARED_FILES_NOT_FULL_BUILD_ATTESTATION",
        },
    )

    problem = SpinArtifactProblem(
        original=(original_vertices, faces),
        fixed_regions=[("shell", fixed_shell, 1.0)],
        design_regions=[("design", (design_vertices, faces), 1.0)],
        length_unit="25mm",
        mass_unit="18.75g",
        target={"mode": "SPINNING_TOP", "axis": [0, 0, 1], "origin": [0, 0, 0], "support_points": []},
        minimum_outer_wall=0.015,
    )
    write(output / "input.json", problem.identity)
    print("Running frozen non-principal Heart optimization", flush=True)
    result = bridge.solve_and_seal_source_spin2014(
        problem, source_settings=geometry, solver_settings=solver, boundary_settings=boundary
    )
    write(output / "native-run.json", result.identity)
    write(output / "evidence.json", result.evidence)
    for name, mesh in (("material", result.material), ("cavity", result.cavity)):
        if mesh is not None and len(mesh[0]):
            trimesh.Trimesh(*mesh, process=False).export(output / f"{name}.obj")
    if result.artifact is not None:
        save_spin_artifact(result.artifact, output / "sealed")
    unchanged = all(
        hashlib.sha256(Path(path).read_bytes()).hexdigest() == digest
        for path, digest in identities.items()
    )
    summary = {
        "status": result.status,
        "sealed": result.artifact is not None,
        "seal_error": result.seal_error,
        "elapsed_seconds": time.monotonic() - started,
        "source_unchanged": unchanged,
        "problem_digest": problem.input_digest,
        "outer_moments": integrate(original_vertices, faces).tolist(),
        "scientific_digest": result.scientific_digest,
        "initial_axis_is_nonprincipal_by_construction": True,
        "original_paper_case": False,
    }
    write(output / "summary.json", summary)
    print(json.dumps({k: summary[k] for k in ("status", "sealed", "source_unchanged", "problem_digest")}), flush=True)
    return 0 if result.artifact is not None and unchanged else 2


if __name__ == "__main__":
    raise SystemExit(main())
