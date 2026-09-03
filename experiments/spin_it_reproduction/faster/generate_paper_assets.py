#!/usr/bin/env python3
"""Generate deterministic structural counterparts for Faster paper examples."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Callable

import numpy as np


Mesh = tuple[list[tuple[float, float, float]], list[tuple[int, int, int]]]
BOX_FACES = (
    (0, 2, 1), (0, 3, 2),
    (4, 5, 6), (4, 6, 7),
    (0, 1, 5), (0, 5, 4),
    (3, 7, 6), (3, 6, 2),
    (0, 4, 7), (0, 7, 3),
    (1, 2, 6), (1, 6, 5),
)


def subdivided_cube(levels: int = 2) -> Mesh:
    vertices = [
        (-1.0, -1.0, -1.0), (1.0, -1.0, -1.0),
        (1.0, 1.0, -1.0), (-1.0, 1.0, -1.0),
        (-1.0, -1.0, 1.0), (1.0, -1.0, 1.0),
        (1.0, 1.0, 1.0), (-1.0, 1.0, 1.0),
    ]
    faces = list(BOX_FACES)
    for _ in range(levels):
        midpoints: dict[tuple[int, int], int] = {}

        def midpoint(first: int, second: int) -> int:
            edge = tuple(sorted((first, second)))
            if edge not in midpoints:
                a, b = vertices[first], vertices[second]
                midpoints[edge] = len(vertices)
                vertices.append(tuple((a[i] + b[i]) * 0.5 for i in range(3)))
            return midpoints[edge]

        refined: list[tuple[int, int, int]] = []
        for a, b, c in faces:
            ab, bc, ca = midpoint(a, b), midpoint(b, c), midpoint(c, a)
            refined.extend(((a, ab, ca), (ab, b, bc), (ca, bc, c), (ab, bc, ca)))
        faces = refined
    return vertices, faces


def warped_cube(scale: float, mapping: Callable[[float, float, float], tuple[float, float, float]]) -> Mesh:
    vertices, faces = subdivided_cube()
    return [mapping(scale * x, scale * y, scale * z) for x, y, z in vertices], faces


def affine_cube(
    scales: tuple[float, float, float],
    *,
    levels: int,
    angle: float = 0.0,
    center: tuple[float, float, float] = (0.0, 0.0, 0.0),
) -> Mesh:
    vertices, faces = subdivided_cube(levels)
    cosine, sine = math.cos(angle), math.sin(angle)
    mapped = []
    for x, y, z in vertices:
        x, y, z = scales[0] * x, scales[1] * y, scales[2] * z
        mapped.append((
            cosine * x - sine * y + center[0],
            sine * x + cosine * y + center[1],
            z + center[2],
        ))
    return mapped, faces


def torus(major: float, minor: float, nu: int = 48, nv: int = 20) -> Mesh:
    vertices: list[tuple[float, float, float]] = []
    for iu in range(nu):
        u = 2.0 * math.pi * iu / nu
        for iv in range(nv):
            v = 2.0 * math.pi * iv / nv
            radius = major + minor * math.cos(v)
            vertices.append((radius * math.cos(u), radius * math.sin(u), minor * math.sin(v)))
    faces: list[tuple[int, int, int]] = []
    for iu in range(nu):
        for iv in range(nv):
            a = iu * nv + iv
            b = ((iu + 1) % nu) * nv + iv
            c = ((iu + 1) % nu) * nv + (iv + 1) % nv
            d = iu * nv + (iv + 1) % nv
            faces.extend(((a, b, c), (a, c, d)))
    return vertices, faces


def radial_extrusion(scale: float, lobes: int, samples: int = 96) -> Mesh:
    bottom = 0
    top = 1
    vertices: list[tuple[float, float, float]] = [(0.0, 0.0, -0.22 * scale), (0.0, 0.0, 0.22 * scale)]
    for z in (-0.22 * scale, 0.22 * scale):
        for index in range(samples):
            angle = 2.0 * math.pi * index / samples
            radius = scale * (1.0 + 0.22 * math.cos(lobes * angle))
            vertices.append((radius * math.cos(angle), radius * math.sin(angle), z))
    faces: list[tuple[int, int, int]] = []
    bottom_offset = 2
    top_offset = 2 + samples
    for index in range(samples):
        following = (index + 1) % samples
        bi, bj = bottom_offset + index, bottom_offset + following
        ti, tj = top_offset + index, top_offset + following
        faces.extend(((bottom, bj, bi), (top, ti, tj)))
        faces.extend(((bi, bj, tj), (bi, tj, ti)))
    return vertices, faces


def write_obj(path: Path, mesh: Mesh) -> str:
    vertices, faces = mesh
    text = "".join(
        [f"v {x:.17g} {y:.17g} {z:.17g}\n" for x, y, z in vertices]
        + [f"f {a + 1} {b + 1} {c + 1}\n" for a, b, c in faces]
    )
    path.write_text(text, encoding="utf-8")
    return hashlib.sha256(text.encode()).hexdigest()


def raw_mass_properties(mesh: Mesh) -> tuple[float, np.ndarray, np.ndarray]:
    vertices, faces = mesh
    mass = 0.0
    first = np.zeros(3)
    second = np.zeros((3, 3))
    for ia, ib, ic in faces:
        a, b, c = (np.asarray(vertices[index], dtype=float) for index in (ia, ib, ic))
        volume = float(np.dot(a, np.cross(b, c)) / 6.0)
        mass += volume
        first += volume * (a + b + c) / 4.0
        vectors = (a, b, c)
        second += volume * sum(np.outer(value, value) for value in vectors) / 10.0
        second += volume * sum(
            np.outer(vectors[i], vectors[j]) + np.outer(vectors[j], vectors[i])
            for i in range(3) for j in range(i + 1, 3)
        ) / 20.0
    return mass, first, second


def principal_target_frame(mesh: Mesh, *, tip_origin: bool) -> tuple[list[float], list[list[float]]]:
    mass, first, second = raw_mass_properties(mesh)
    center = first / mass
    inertia = np.trace(second) * np.eye(3) - second
    inertia -= mass * ((center @ center) * np.eye(3) - np.outer(center, center))
    _eigenvalues, eigenvectors = np.linalg.eigh(inertia)
    third_index = int(np.argmax(np.abs(eigenvectors[2, :])))
    e3 = eigenvectors[:, third_index]
    if e3[2] < 0.0:
        e3 = -e3
    first_index = next(index for index in range(3) if index != third_index)
    e1 = eigenvectors[:, first_index]
    e1 -= e3 * float(e1 @ e3)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(e3, e1)
    axes = np.column_stack((e1, e2, e3))
    origin = center.copy()
    if tip_origin:
        projections = (np.asarray(mesh[0], dtype=float) - center) @ e3
        origin += float(np.min(projections)) * e3
    return origin.tolist(), axes.tolist()


def pair(kind: str) -> tuple[Mesh, Mesh, str]:
    if kind == "vase":
        def shape(x: float, y: float, z: float) -> tuple[float, float, float]:
            taper = 1.0 + 0.22 * z
            return (taper * x + 0.28 * z, taper * y, 1.45 * z)
        return warped_cube(1.0, shape), warped_cube(0.68, shape), "sheared tapered vessel"
    if kind == "tea-infuser":
        return torus(1.25, 0.48), torus(1.25, 0.30), "toroidal body with nonzero genus"
    if kind == "island-cat":
        def shape(x: float, y: float, z: float) -> tuple[float, float, float]:
            return (1.55 * x + 0.20 * z, 1.05 * y, 0.32 * z + 0.10 * x)
        return warped_cube(1.0, shape), warped_cube(0.67, shape), "asymmetric floating platform A"
    if kind == "island-dog-human":
        def shape(x: float, y: float, z: float) -> tuple[float, float, float]:
            return (1.75 * x - 0.18 * y, 0.92 * y + 0.12 * z, 0.38 * z - 0.08 * x)
        return warped_cube(1.0, shape), warped_cube(0.66, shape), "attachable asymmetric platform B"
    if kind == "heart":
        def shape(x: float, y: float, z: float) -> tuple[float, float, float]:
            bulge = 1.0 + 0.18 * (1.0 - z * z)
            return (bulge * x + 0.10 * z, 0.82 * bulge * y, 1.30 * z - 0.12 * abs(x))
        return warped_cube(1.0, shape), warped_cube(0.67, shape), "asymmetric heart-like top body"
    if kind == "noodles":
        def shape(x: float, y: float, z: float) -> tuple[float, float, float]:
            angle = 0.58 * z
            c, s = math.cos(angle), math.sin(angle)
            return (c * x - s * y, s * x + c * y, 1.35 * z)
        return warped_cube(1.0, shape), warped_cube(0.64, shape), "twisted high-complexity proxy"
    if kind == "jellyfish":
        # A stretched, slightly rotated faceted body preserves the paper
        # problem's off-axis starting frame while retaining an exact, simple
        # interior design domain for robust Boolean differentiation.
        outer = affine_cube(
            (1.30, 1.29983286, 1.83696459),
            levels=3,
            angle=0.04255563,
            center=(0.0, 0.0, -0.059454092034775544),
        )
        design = affine_cube((1.0, 1.0, 1.0), levels=0)
        return outer, design, "faceted elongated jellyfish-like yo-yo counterpart"
    if kind == "flower":
        outer = radial_extrusion(1.25, 7)
        design = affine_cube((0.55, 0.55, 0.15), levels=0)
        return outer, design, "seven-lobed planar spinner with a protected cuboid core"
    raise ValueError(kind)


CASE_DEFINITIONS = (
    ("static-tilted-vase", "vase", "STATIC_BALANCE"),
    ("floating-tea-infuser", "tea-infuser", "BUOYANT_STABILITY"),
    ("floating-island-cat", "island-cat", "BUOYANT_STABILITY"),
    ("floating-island-dog-human", "island-dog-human", "BUOYANT_STABILITY"),
    ("heart-spinning-top", "heart", "SPINNING_TOP"),
    ("noodles-spinning-top", "noodles", "SPINNING_TOP"),
    ("jellyfish-yoyo", "jellyfish", "YO_YO"),
    ("flower-fidget-spinner", "flower", "FIDGET_SPINNER"),
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    output = args.output.resolve()
    assets = output / "paper-assets"
    assets.mkdir(parents=True, exist_ok=True)
    cases = []
    for case_id, kind, mode in CASE_DEFINITIONS:
        outer, design, note = pair(kind)
        outer_name = f"{case_id}-outer.obj"
        design_name = f"{case_id}-design.obj"
        outer_hash = write_obj(assets / outer_name, outer)
        design_hash = write_obj(assets / design_name, design)
        case = {
            "id": case_id,
            "paper_section": "6.2",
            "mode": mode,
            "asset_class": "RECONSTRUCTED_ASSET",
            "author_exact": False,
            "physical_validation": False,
            "units": "normalized_unit",
            "target_origin": [0.0, 0.0, 0.0],
            "target_axes": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            "construction": note,
            "license": "project-generated deterministic geometry; CC0-1.0",
            "outer": {"path": f"paper-assets/{outer_name}", "sha256": outer_hash},
            "design_domain": {"path": f"paper-assets/{design_name}", "sha256": design_hash},
        }
        if mode == "STATIC_BALANCE":
            case["static_objective"] = "MAXIMUM_MASS"
        elif mode == "SPINNING_TOP":
            case["center_of_mass_weight"] = 0.25
            case["target_origin"], case["target_axes"] = principal_target_frame(
                outer,
                tip_origin=True,
            )
        elif mode == "YO_YO":
            case["target_origin"], case["target_axes"] = principal_target_frame(
                outer,
                tip_origin=False,
            )
        elif mode == "BUOYANT_STABILITY":
            # Targets are the independently integrated total moments of the
            # listed feasible plane.  This guarantees a nonempty equilibrium
            # manifold without embedding a solved optimum.
            feasible = {
                "floating-tea-infuser": (4.59625251660662, 0.1849352704595321, 0.03698633309075461),
                "floating-island-cat": (3.5121711334769143, 0.08558581166706312, 0.005404338280521669),
                "floating-island-dog-human": (4.325410392245944, 0.09129467930121796, 0.0002290752576555812),
            }[case_id]
            case["known_feasible_quadric"] = [0.2, 1.0, 0.2, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            case["buoyancy"] = {
                "liquid_density": 1.0,
                "displaced_liquid_volume": feasible[0],
                "center_of_buoyancy_target": [feasible[1], feasible[2], 0.0],
                "waterline_area": 4.0,
                "waterline_central_second_moment": [[1.0, 0.0], [0.0, 1.0]],
            }
        cases.append(case)
    payload = {
        "schema": "spin-it-faster-paper-cases-v1",
        "schema_version": 1,
        "paper_doi": "10.1145/3658194",
        "asset_class": "RECONSTRUCTED_ASSET",
        "author_exact": False,
        "physical_validation": False,
        "cases": cases,
    }
    (output / "paper_cases.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
