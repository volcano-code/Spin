"""Run a frozen rotated reconstructed Heart with an interior cube design domain.

This is a reconstructed algorithm-validation asset.  The outer Heart, the
rotation, and the centered interior cube are fixed before the solver starts;
the cube is not selected from an optimization result.  It is deliberately
separate from the scaled-Heart design domain used by run_rotated.py.
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


def rotate_about_center(vertices: np.ndarray, center: np.ndarray, degrees: float) -> np.ndarray:
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
    parser.add_argument("--core-half-side", type=float, default=0.18)
    parser.add_argument("--initial-level", type=int, default=3)
    parser.add_argument("--maximum-level", type=int, default=4)
    args = parser.parse_args()
    output = args.output
    output.mkdir(parents=False, exist_ok=False)
    started = time.monotonic()

    from make_it_stand.verification import source_spin_material as bridge
    from make_it_stand.verification.spin_common import native
    from make_it_stand.verification.spin_problem import SpinArtifactProblem
    from make_it_stand.verification.spin_bundle import save_spin_artifact

    original_vertices, outer_faces = heart_mesh()
    moments = integrate(original_vertices, outer_faces)
    input_centroid = moments[1:4] / moments[0]
    outer_vertices = rotate_about_center(original_vertices, input_centroid, args.rotation_degrees) - input_centroid
    outer = (outer_vertices, outer_faces)

    core = trimesh.creation.box(extents=(2.0 * args.core_half_side,) * 3)
    core_vertices = np.asarray(core.vertices, dtype=float)
    core_faces = np.asarray(core.faces, dtype=int)
    core = (core_vertices, core_faces)
    # The inner boundary is wound oppositely in the fixed-shell component.
    fixed_shell = (
        np.vstack((outer_vertices, core_vertices)),
        np.vstack((outer_faces, core_faces[:, ::-1] + len(outer_vertices))),
    )

    for name, mesh in (("outer", outer), ("design", core), ("fixed-shell", fixed_shell)):
        trimesh.Trimesh(*mesh, process=False).export(output / f"{name}.obj")

    geometry = bridge.SourceGeometrySettings(root_origin=(-2.0, -2.0, -2.0), root_side=4.0, maximum_depth=5)
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
        "asset": "RECONSTRUCTED_SELF_BUILT_HEART_ROTATED_INTERIOR_CUBE_V1",
        "author_original": False,
        "rotation_about_y_degrees": float(args.rotation_degrees),
        "rotation_center_input_model": input_centroid.tolist(),
        "translation_to_target_frame": (-input_centroid).tolist(),
        "core_half_side": float(args.core_half_side),
        "initial_level": int(args.initial_level),
        "maximum_optimization_level": int(args.maximum_level),
        "target_axis_model": [0.0, 0.0, 1.0],
        "source_center_model": [0.0, 0.0, 0.0],
        "source": asdict(geometry),
        "solver": asdict(solver),
        "boundary": asdict(boundary),
        "minimum_outer_wall": 0.015,
        "solver_density": 1.0,
        "meters_per_model_unit": 0.025,
        "physical_density_kg_per_m3": 1200.0,
        "protocol": "NON_PRINCIPAL_INITIAL_AXIS_FIXED_INTERIOR_CUBE_V1",
        "design_domain_selection": "DECLARED_BEFORE_SOLVE_CENTERED_CUBE",
    }
    write(output / "config.json", config)
    root = Path(__file__).resolve().parents[2]
    identity_paths = [
        Path(__file__), Path(__file__).with_name("heart_top.py"), Path(native().__file__),
        root / "python/make_it_stand/verification/source_spin_material.py",
    ]
    identities = {str(path.resolve()): hashlib.sha256(path.read_bytes()).hexdigest() for path in identity_paths}
    write(output / "source.json", {
        "head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip(),
        "files": identities,
        "scope": "LOADED_BINARY_AND_DECLARED_FILES_NOT_FULL_BUILD_ATTESTATION",
    })

    problem = SpinArtifactProblem(
        original=outer,
        fixed_regions=[("shell", fixed_shell, 1.0)],
        design_regions=[("design", core, 1.0)],
        length_unit="25mm",
        mass_unit="18.75g",
        target={"mode": "SPINNING_TOP", "axis": [0, 0, 1], "origin": [0, 0, 0], "support_points": []},
        minimum_outer_wall=0.015,
    )
    write(output / "input.json", problem.identity)
    print("Running frozen rotated Heart with declared interior cube", flush=True)
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
    unchanged = all(hashlib.sha256(Path(path).read_bytes()).hexdigest() == digest for path, digest in identities.items())
    summary = {
        "status": result.status,
        "sealed": result.artifact is not None,
        "seal_error": result.seal_error,
        "elapsed_seconds": time.monotonic() - started,
        "source_unchanged": unchanged,
        "problem_digest": problem.input_digest,
        "scientific_digest": result.scientific_digest,
        "initial_axis_is_nonprincipal_by_construction": True,
        "original_paper_case": False,
    }
    write(output / "summary.json", summary)
    print(json.dumps({k: summary[k] for k in ("status", "sealed", "source_unchanged", "problem_digest")}), flush=True)
    return 0 if result.artifact is not None and unchanged else 2


if __name__ == "__main__":
    raise SystemExit(main())
