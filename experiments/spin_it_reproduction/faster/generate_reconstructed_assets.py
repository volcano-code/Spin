#!/usr/bin/env python3
"""Generate deterministic, license-free OBJ fixtures for the Faster matrix."""

from __future__ import annotations

import math
from pathlib import Path


FACES = (
    (0, 2, 1), (0, 3, 2),
    (4, 5, 6), (4, 6, 7),
    (0, 1, 5), (0, 5, 4),
    (3, 7, 6), (3, 6, 2),
    (0, 4, 7), (0, 7, 3),
    (1, 2, 6), (1, 6, 5),
)


def box(half_x: float, half_y: float, half_z: float, angle: float = 0.0,
        center: tuple[float, float, float] = (0.0, 0.0, 0.0)):
    cosine, sine = math.cos(angle), math.sin(angle)
    vertices = []
    for x, y, z in (
        (-half_x, -half_y, -half_z), (half_x, -half_y, -half_z),
        (half_x, half_y, -half_z), (-half_x, half_y, -half_z),
        (-half_x, -half_y, half_z), (half_x, -half_y, half_z),
        (half_x, half_y, half_z), (-half_x, half_y, half_z),
    ):
        vertices.append((
            cosine * x - sine * y + center[0],
            sine * x + cosine * y + center[1],
            z + center[2],
        ))
    return vertices, list(FACES)


def subdivide(vertices, faces, levels: int):
    vertices = list(vertices)
    faces = list(faces)
    for _ in range(levels):
        midpoints: dict[tuple[int, int], int] = {}

        def midpoint(first: int, second: int) -> int:
            edge = tuple(sorted((first, second)))
            if edge not in midpoints:
                a, b = vertices[first], vertices[second]
                midpoints[edge] = len(vertices)
                vertices.append(tuple((a[index] + b[index]) * 0.5 for index in range(3)))
            return midpoints[edge]

        refined = []
        for a, b, c in faces:
            ab, bc, ca = midpoint(a, b), midpoint(b, c), midpoint(c, a)
            refined.extend(((a, ab, ca), (ab, b, bc), (ca, bc, c), (ab, bc, ca)))
        faces = refined
    return vertices, faces


def write_obj(path: Path, vertices, faces) -> None:
    lines = [f"v {x:.17g} {y:.17g} {z:.17g}\n" for x, y, z in vertices]
    lines.extend(f"f {a + 1} {b + 1} {c + 1}\n" for a, b, c in faces)
    path.write_text("".join(lines), encoding="utf-8")


def main() -> int:
    output = Path(__file__).resolve().parent / "assets"
    output.mkdir(parents=True, exist_ok=True)
    design = box(1.0, 1.0, 1.0)
    common = subdivide(*box(1.2, 1.2, 1.2), levels=3)
    yoyo = subdivide(
        *box(
            1.3,
            1.29983286,
            1.83696459,
            angle=0.04255563,
            center=(-2.6266152708628194e-17,
                    -1.669577259936739e-17,
                    -0.059454092034775544),
        ),
        levels=3,
    )
    write_obj(output / "design-domain.obj", *design)
    write_obj(output / "faceted-rotor-shell.obj", *common)
    write_obj(output / "faceted-yoyo-shell.obj", *yoyo)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
