#!/usr/bin/env python3
"""Generate deterministic voxel counterparts for every Spin-It 2014 example."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Callable


SIDE = 16
Predicate = Callable[[float, float, float], bool]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def encode_rle(values: list[int]) -> list[list[int]]:
    if not values:
        return []
    runs: list[list[int]] = []
    value, count = values[0], 1
    for current in values[1:]:
        if current == value:
            count += 1
        else:
            runs.append([value, count])
            value, count = current, 1
    runs.append([value, count])
    return runs


def decode_rle(runs: list[list[int]]) -> list[int]:
    result: list[int] = []
    for value, count in runs:
        if value not in (0, 1) or count <= 0:
            raise ValueError("invalid binary run-length encoding")
        result.extend([value] * count)
    return result


def sphere(x: float, y: float, z: float, cx: float, cy: float, cz: float, radius: float) -> bool:
    return (x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2 <= radius ** 2


def ellipsoid(x: float, y: float, z: float, rx: float, ry: float, rz: float) -> bool:
    return (x / rx) ** 2 + (y / ry) ** 2 + (z / rz) ** 2 <= 1.0


def occupancy(
    predicate: Predicate,
    *,
    protect_diagonal_contacts: bool = False,
) -> tuple[list[int], list[int]]:
    material: list[int] = []
    for z_index in range(SIDE):
        z = -1.0 + 2.0 * (z_index + 0.5) / SIDE
        for y_index in range(SIDE):
            y = -1.0 + 2.0 * (y_index + 0.5) / SIDE
            for x_index in range(SIDE):
                x = -1.0 + 2.0 * (x_index + 0.5) / SIDE
                material.append(int(predicate(x, y, z)))
    carveable = [0] * len(material)
    for z in range(1, SIDE - 1):
        for y in range(1, SIDE - 1):
            for x in range(1, SIDE - 1):
                index = x + SIDE * (y + SIDE * z)
                if protect_diagonal_contacts:
                    neighbors = tuple(
                        (x + dx) + SIDE * ((y + dy) + SIDE * (z + dz))
                        for dz in (-1, 0, 1)
                        for dy in (-1, 0, 1)
                        for dx in (-1, 0, 1)
                        if (dx, dy, dz) != (0, 0, 0)
                    )
                else:
                    neighbors = (
                        index - 1, index + 1, index - SIDE, index + SIDE,
                        index - SIDE * SIDE, index + SIDE * SIDE,
                    )
                carveable[index] = int(
                    material[index] and all(material[neighbor] for neighbor in neighbors)
                )
    if not any(carveable):
        raise RuntimeError("generated model has no one-voxel protected shell")
    return material, carveable


def predicate(kind: str) -> Predicate:
    if kind == "ellipsoid":
        return lambda x, y, z: ellipsoid(x, y, z, 0.82, 0.7, 0.5)
    if kind == "heart":
        return lambda x, y, z: (
            ((x * x + 1.15 * y * y + (z / 0.62) ** 2 - 0.68) ** 3 -
             0.42 * x * x * (z / 0.62) ** 3) <= 0.0
        )
    if kind == "armadillo-hand":
        return lambda x, y, z: (
            ellipsoid(x, y, z + 0.03, 0.72, 0.62, 0.52) or
            sphere(x, y, z, 0.28, 0.0, 0.28, 0.22) or
            sphere(x, y, z, -0.28, 0.0, 0.28, 0.22)
        )
    if kind == "armadillo-shell":
        return lambda x, y, z: (
            ellipsoid(x, y, z, 0.74, 0.62, 0.5) or
            sphere(x, y, z, -0.32, 0.0, 0.25, 0.22)
        )
    if kind == "teapot":
        return lambda x, y, z: (
            ellipsoid(x, y, z, 0.7, 0.62, 0.42) or
            ellipsoid(x - 0.53, y, z - 0.04, 0.34, 0.28, 0.27) or
            ellipsoid(x + 0.53, y, z - 0.04, 0.34, 0.28, 0.27)
        )
    if kind == "cuboid":
        return lambda x, y, z: abs(x) <= 0.75 and abs(y) <= 0.68 and abs(z) <= 0.38
    if kind == "woven-ring":
        return lambda x, y, z: ellipsoid(x, y, z, 0.82, 0.82, 0.34)
    if kind == "couple":
        return lambda x, y, z: ellipsoid(x, y, z, 0.82, 0.68, 0.5)
    if kind == "sumone":
        return lambda x, y, z: (
            ellipsoid(x, y, z + 0.12, 0.68, 0.58, 0.42) or
            ellipsoid(x, y, z - 0.32, 0.32, 0.3, 0.22)
        )
    if kind == "trex":
        return lambda x, y, z: (
            ellipsoid(x + 0.12, y, z + 0.1, 0.72, 0.34, 0.48) or
            ellipsoid(x + 0.3, y, z - 0.2, 0.43, 0.25, 0.3) or
            sphere(x, y, z, 0.45, 0.0, 0.25, 0.34)
        )
    if kind == "elephant":
        return lambda x, y, z: (
            ellipsoid(x, y, z, 0.72, 0.62, 0.42) or
            ellipsoid(x - 0.42, y, z - 0.12, 0.3, 0.28, 0.25) or
            ellipsoid(x + 0.42, y, z - 0.12, 0.3, 0.28, 0.25)
        )
    raise ValueError(kind)


CASES = (
    ("ellipsoid-top", "ellipsoid", "SPINNING_TOP", 1, 1, "Figure 6"),
    ("heart-top", "heart", "SPINNING_TOP", 1, 1, "Figure 7"),
    ("armadillo-hand", "armadillo-hand", "SPINNING_TOP", 1, 3, "Figure 8"),
    ("armadillo-shell", "armadillo-shell", "SPINNING_TOP", 5, 2, "Figure 8"),
    ("teapot-inertia", "teapot", "SPINNING_TOP", 2, 1, "Figure 9 left"),
    ("teapot-lowering", "teapot", "SPINNING_TOP", 2, 1, "Figure 9 middle"),
    ("teapot-mass-reduction", "teapot", "SPINNING_TOP", 2, 1, "Figure 9 right"),
    ("cuboid-yoyo", "cuboid", "YO_YO", 1, 1, "Figure 10 top"),
    ("woven-ring-yoyo", "woven-ring", "YO_YO", 1, 4, "Figure 10 bottom"),
    ("dancing-couple-hollowing", "couple", "SPINNING_TOP", 1, 5, "Figure 11"),
    ("dancing-couple-deformation", "couple", "SPINNING_TOP", 1, 5, "Figure 11"),
    ("sumone-dual-density", "sumone", "SPINNING_TOP", 1, 1, "Figure 5"),
    ("trex-static", "trex", "STATIC_BALANCE", 1, 1, "Figure 12"),
    ("elephant-accessibility", "elephant", "SPINNING_TOP", 3, 2, "Table 1"),
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    output = args.output.resolve()
    assets = output / "paper-assets"
    assets.mkdir(parents=True, exist_ok=True)
    rows = []
    for case_id, kind, mode, voids, cuts, role in CASES:
        material, carveable = occupancy(
            predicate(kind),
            protect_diagonal_contacts=kind in {"couple", "trex"},
        )
        if kind == "trex":
            # Keep the reconstructed T-Rex outer shell intentionally
            # off-center, while exposing one connected negative-x cavity.
            # The cavity's complete removal exactly cancels the shell's
            # lateral first moment on the 16^3 source lattice.  Because every
            # carveable cell has negative x, the binary balanced solution is
            # unique and remains valid under dyadic level-9 refinement.
            carveable = [
                value if index % SIDE <= 5 else 0
                for index, value in enumerate(carveable)
            ]
        asset_name = f"{case_id}.vox.json"
        # Top/static coordinates use the support/contact point as z=0, as in
        # the paper's body frame.  Yo-yos use the axle center as the origin.
        minimum_z, maximum_z = (-1.0, 1.0) if mode == "YO_YO" else (0.0, 2.0)
        asset = {
            "schema": "spin-it-2014-voxel-asset-v1",
            "dimensions": [SIDE, SIDE, SIDE],
            "minimum": [-1.0, -1.0, minimum_z],
            "maximum": [1.0, 1.0, maximum_z],
            "material_rle": encode_rle(material),
            "carveable_rle": encode_rle(carveable),
        }
        path = assets / asset_name
        path.write_text(json.dumps(asset, sort_keys=True, separators=(",", ":")) + "\n", encoding="utf-8")
        row = {
            "id": case_id,
            "paper_role": role,
            "mode": mode,
            "asset_class": "RECONSTRUCTED_ASSET",
            "author_exact": False,
            "physical_validation": False,
            "target_axis": [0.0, 0.0, 1.0],
            "center_of_mass_weight": 1.0,
            "initial_fill": 0.5,
            "laplacian_weight": 0.001,
            "paper_level_laplacian_weight": 0.00001,
            "laplacian_weight_evidence": "RECONSTRUCTED_PARAMETER",
            "maximum_octree_level": 4,
            "paper_maximum_octree_level": 9,
            "binary_constraint_tolerance": 0.005,
            "paper_binary_constraint_tolerance": 1e-6,
            "expected_void_components": voids,
            "manual_planar_cuts": cuts,
            "accessibility": {
                "paper_reported_manual_cuts": cuts,
                "automatic_cut_optimization_in_paper": False,
            },
            "occupancy": {"path": f"paper-assets/{asset_name}", "sha256": sha256(path)},
            "source_attribution": "project-generated structural counterpart; CC0-1.0",
        }
        if case_id == "teapot-inertia":
            row["center_of_mass_weight"] = 0.0
        elif case_id == "teapot-lowering":
            row["center_of_mass_weight"] = 2.0
        elif case_id == "teapot-mass-reduction":
            row["center_of_mass_weight"] = 1.0
            row["mass_reduction_extension"] = True
        elif case_id == "dancing-couple-deformation":
            row["cage"] = {
                "type": "AUTOMATIC_BBOX",
                "deformation_enabled": True,
                "maximum_handle_count": 2,
                "projected_gradient_tolerance": 1e-6,
            }
        elif case_id == "sumone-dual-density":
            row["densities_g_cm3"] = {"shell": 1.17, "insert": 8.1}
        elif case_id == "trex-static":
            # The exact binary solution lies at the upper box bound.  A
            # feasible warm start avoids asking the equality-constrained SLQP
            # line search to cross that bound from an infeasible midpoint.
            row["initial_fill"] = 1.0
            row["initial_fill_evidence"] = "RECONSTRUCTED_PARAMETER"
        rows.append(row)
    manifest = {
        "schema": "spin-it-2014-paper-cases-v1",
        "schema_version": 1,
        "paper_doi": "10.1145/2601097.2601157",
        "asset_class": "RECONSTRUCTED_ASSET",
        "author_exact": False,
        "physical_validation": False,
        "paper_resolution_reproduced": False,
        "cases": rows,
    }
    (output / "paper_cases.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
