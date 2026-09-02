#!/usr/bin/env python3
"""Reconstruct a strict outer surface from an official occupancy cache."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from skimage import measure
import trimesh

from make_it_stand.official_beta import ReleasedBetaArchive


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--beta-root", required=True, type=Path)
    parser.add_argument("--case", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--hausdorff-samples", type=int, default=100_000)
    arguments = parser.parse_args()
    root = arguments.beta_root.resolve()
    archive = ReleasedBetaArchive.load(
        root / "data" / "configs" / f"{arguments.case}.mis",
        beta_root=root,
    )
    resolution = archive.voxel_resolution
    occupied = archive.cell_ids_z_fast.reshape((resolution,) * 3) >= 0
    field = np.pad(occupied.astype(np.float32), 1)
    vertices, faces, _, _ = measure.marching_cubes(field, level=0.5)
    vertices = (vertices.astype(np.float64) - 0.5) / float(resolution)
    # marching_cubes returns float32 coordinates. Snap its half-cell lattice
    # before OBJ serialization so boundary points remain inside the exact
    # official deformation domain instead of drifting ~1e-8 across a face.
    vertices = np.rint(vertices * (2 * resolution)) / (2 * resolution)
    reconstructed = trimesh.Trimesh(vertices=vertices, faces=faces, process=True)
    components = reconstructed.split(only_watertight=False)
    reconstructed = max(components, key=lambda mesh: abs(float(mesh.volume)))
    if reconstructed.volume < 0.0:
        reconstructed.invert()
    if not reconstructed.is_watertight or not reconstructed.is_winding_consistent:
        raise ValueError("occupancy reconstruction is not a closed oriented surface")

    source = trimesh.Trimesh(
        vertices=archive.normalized_vertices,
        faces=archive.faces,
        process=False,
    )
    rng = np.random.default_rng(2013)
    source_points, _ = trimesh.sample.sample_surface(
        source,
        arguments.hausdorff_samples,
        seed=rng,
    )
    reconstructed_points, _ = trimesh.sample.sample_surface(
        reconstructed,
        arguments.hausdorff_samples,
        seed=rng,
    )
    _, source_to_reconstructed, _ = trimesh.proximity.closest_point(
        reconstructed,
        source_points,
    )
    _, reconstructed_to_source, _ = trimesh.proximity.closest_point(
        source,
        reconstructed_points,
    )

    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    obj_text = trimesh.exchange.obj.export_obj(reconstructed, digits=17)
    arguments.output.write_text(obj_text, encoding="utf-8")
    report = {
        "schema": "make-it-stand-occupancy-surface-reconstruction-v1",
        "schema_version": 1,
        "case": arguments.case,
        "classification": "reconstructed_from_official_occupancy_not_original_surface",
        "voxel_resolution": resolution,
        "occupied_cells": archive.occupied_cell_count,
        "source_vertices": int(len(source.vertices)),
        "source_faces": int(len(source.faces)),
        "reconstructed_vertices": int(len(reconstructed.vertices)),
        "reconstructed_faces": int(len(reconstructed.faces)),
        "reconstructed_watertight": bool(reconstructed.is_watertight),
        "reconstructed_winding_consistent": bool(
            reconstructed.is_winding_consistent
        ),
        "reconstructed_signed_volume": float(reconstructed.volume),
        "sample_count_each_direction": arguments.hausdorff_samples,
        "sampled_source_to_reconstructed_max": float(
            np.max(source_to_reconstructed)
        ),
        "sampled_reconstructed_to_source_max": float(
            np.max(reconstructed_to_source)
        ),
        "sampled_symmetric_hausdorff": float(max(
            np.max(source_to_reconstructed),
            np.max(reconstructed_to_source),
        )),
        "output_sha256": _sha256(arguments.output),
        "claim_boundary": (
            "The distance is a deterministic sampled Hausdorff estimate, not an "
            "exact continuous Hausdorff proof. The reconstructed mesh is for "
            "simulation/compatibility evidence and is not the untouched paper input."
        ),
    }
    report_path = arguments.output.with_suffix(".reconstruction.json")
    report_path.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
