"""Run a frozen, reconstructed non-principal cuboid yo-yo case.

The outer cuboid is rotated about its centroid.  The design domain is a pair
of equal rotated boxes at opposite diagonal positions, so removing equal
amounts preserves the first moments while changing Q_xz.  This is a method
validation asset, never an author-case claim.
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


def write(path: Path, value):
    path.write_text(json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n")


def target_frame_for_axis(axis):
    """Return the deterministic model-to-target right-handed frame for an axis.

    The third *column* is the declared physical target axis.  The other two
    columns are selected from the least-aligned model basis direction, so the
    coordinate convention is frozen by the input axis rather than a solver
    outcome or an arbitrary runtime choice.
    """
    direction = np.asarray(axis, dtype=float)
    if direction.shape != (3,) or not np.all(np.isfinite(direction)):
        raise ValueError("target axis must contain three finite values")
    length = float(np.linalg.norm(direction))
    if not length > 0.0:
        raise ValueError("target axis must be nonzero")
    direction /= length
    reference = np.eye(3)[int(np.argmin(np.abs(direction)))]
    first = reference - float(reference @ direction) * direction
    first /= float(np.linalg.norm(first))
    second = np.cross(direction, first)
    return tuple(tuple(float(value) for value in row)
                 for row in np.column_stack((first, second, direction)))


def rotated_box(center, extents, angle):
    mesh = trimesh.creation.box(extents=np.asarray(extents, dtype=float))
    rotation = trimesh.transformations.rotation_matrix(float(angle), [0.0, 1.0, 0.0])
    transform = rotation.copy()
    transform[:3, 3] = np.asarray(center, dtype=float)
    mesh.apply_transform(transform)
    return mesh


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--initial-level", type=int, default=2)
    parser.add_argument("--maximum-level", type=int, default=4)
    parser.add_argument(
        "--design-layout",
        choices=("paired-boxes", "interior-core", "offaxis-pairs"),
        default="paired-boxes",
        help="pre-registered carveable-domain layout; does not depend on solve results",
    )
    parser.add_argument(
        "--outer-rotation-degrees",
        type=float,
        default=20.0,
        help="fixed outer cuboid rotation about y; recorded as part of the asset identity",
    )
    parser.add_argument(
        "--target-axis",
        type=float,
        nargs=3,
        metavar=("X", "Y", "Z"),
        default=(0.0, 0.0, 1.0),
        help="pre-registered physical yo-yo axis in model coordinates",
    )
    args = parser.parse_args()
    output = args.output
    output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()

    from make_it_stand.verification import source_spin_material as bridge
    from make_it_stand.verification.spin_common import native
    from make_it_stand.verification.spin_problem import SpinArtifactProblem
    from make_it_stand.verification.spin_bundle import save_spin_artifact

    angle_degrees = float(args.outer_rotation_degrees)
    if not np.isfinite(angle_degrees) or abs(angle_degrees) >= 45.0:
        raise ValueError("outer rotation must be finite and strictly inside +/-45 degrees")
    angle = np.deg2rad(angle_degrees)
    target_axes = target_frame_for_axis(args.target_axis)
    target_axis = np.asarray(target_axes, dtype=float)[:, 2]
    rotation = trimesh.transformations.rotation_matrix(float(angle), [0.0, 1.0, 0.0])
    outer = rotated_box([0, 0, 0], [2.0, 1.6, 1.2], angle)
    if args.design_layout == "paired-boxes":
        design_parts = [
            rotated_box([0.36, 0.0, 0.24], [0.42, 0.80, 0.28], angle),
            rotated_box([-0.36, 0.0, -0.24], [0.42, 0.80, 0.28], angle),
        ]
    else:
        # A centered, symmetric interior core is fixed before the solve.  It
        # preserves the yo-yo center-of-mass constraint while giving the
        # optimizer a connected domain with enough volume to change Q_xz.
        if args.design_layout == "interior-core":
            design_parts = [rotated_box([0.0, 0.0, 0.0], [1.60, 1.20, 0.80], angle)]
        else:
            # Two equal local boxes occupy opposite x*z quadrants.  Their
            # first moments cancel, while their signed xz moment is large
            # enough to test a non-principal yo-yo target at five degrees.
            local_centers = ([0.56, 0.0, -0.37], [-0.56, 0.0, 0.37])
            extents = [0.75, 1.20, 0.42]
            design_parts = [
                rotated_box((rotation @ np.r_[center, 1.0])[:3], extents, angle)
                for center in local_centers
            ]
    design = trimesh.util.concatenate(design_parts)
    if not outer.contains(design.vertices).all():
        raise RuntimeError("frozen design boxes are not contained by the outer cuboid")

    shell_vertices = np.vstack((outer.vertices, design.vertices))
    shell_faces = np.vstack((outer.faces, design.faces[:, ::-1] + len(outer.vertices)))
    for name, mesh in (("outer", outer), ("design", design)):
        mesh.export(output / f"{name}.obj")
    trimesh.Trimesh(shell_vertices, shell_faces, process=False).export(output / "fixed-shell.obj")

    geometry = bridge.SourceGeometrySettings(
        root_origin=(-1.6, -1.6, -1.6), root_side=3.2, maximum_depth=5
    )
    solver = bridge.SourceSolverSettings(
        objective="YOYO",
        initial_level=int(args.initial_level),
        maximum_optimization_level=int(args.maximum_level),
        maximum_adaptive_passes=10,
        maximum_leaf_count=6000,
        initial_fill=0.5,
        initial_phi=0.0,
        laplacian_weight=0.001,
        target_axes=target_axes,
        slqp=bridge.AuditedSlqpSettings(maximum_iterations=240),
    )
    boundary = bridge.SourceBoundarySettings(maximum_input_cuboids=6000)
    config = {
        "asset": "RECONSTRUCTED_ROTATED_CUBOID_YOYO_V1_" + args.design_layout.upper(),
        "author_original": False,
        "rotation_about_y_degrees": angle_degrees,
        "target_axis": target_axis.tolist(),
        "target_origin": [0.0, 0.0, 0.0],
        "target_axes": [list(row) for row in target_axes],
        "design_layout": args.design_layout,
        "design_parts_local_centers": ([
            [0.36, 0.0, 0.24], [-0.36, 0.0, -0.24]
        ] if args.design_layout == "paired-boxes" else ([[0.0, 0.0, 0.0]] if args.design_layout == "interior-core" else [
            [0.56, 0.0, -0.37], [-0.56, 0.0, 0.37]
        ])),
        "design_parts_local_extents": ([
            [0.42, 0.80, 0.28], [0.42, 0.80, 0.28]
        ] if args.design_layout == "paired-boxes" else ([[1.60, 1.20, 0.80]] if args.design_layout == "interior-core" else [
            [0.75, 1.20, 0.42], [0.75, 1.20, 0.42]
        ])),
        "initial_level": int(args.initial_level),
        "maximum_level": int(args.maximum_level),
        "source": asdict(geometry),
        "solver": asdict(solver),
        "boundary": asdict(boundary),
        "minimum_outer_wall": 0.015,
        "solver_density": 1.0,
        "meters_per_model_unit": 0.025,
        "physical_density_kg_per_m3": 1200.0,
        "protocol": "NON_PRINCIPAL_YOYO_HOLLOWING_ALGORITHM_VALIDATION",
    }
    write(output / "config.json", config)
    root = Path(__file__).resolve().parents[2]
    identities = {
        str(path.resolve()): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in (Path(__file__), Path(native().__file__), root / "python/make_it_stand/verification/source_spin_material.py")
    }
    write(output / "source.json", {
        "head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip(),
        "files": identities,
        "scope": "LOADED_BINARY_AND_DECLARED_FILES_NOT_FULL_BUILD_ATTESTATION",
    })

    problem = SpinArtifactProblem(
        original=(outer.vertices, outer.faces),
        fixed_regions=[("shell", (shell_vertices, shell_faces), 1.0)],
        design_regions=[("design", (design.vertices, design.faces), 1.0)],
        length_unit="25mm",
        mass_unit="18.75g",
        target={"mode": "YO_YO", "axis": target_axis.tolist(), "origin": [0, 0, 0], "support_points": []},
        minimum_outer_wall=0.015,
    )
    write(output / "input.json", problem.identity)
    print("Running frozen non-principal cuboid yo-yo optimization", flush=True)
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
