"""Lossless extraction of the triangle surface from a 3-D Medit mesh."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class MeditSurface:
    vertices: np.ndarray
    faces: np.ndarray
    medit_vertex_count: int
    surface_vertex_count: int
    triangle_count: int
    excluded_unreferenced_vertex_count: int


def _clean_lines(path: Path) -> list[str]:
    try:
        content = path.read_text(encoding="ascii")
    except (OSError, UnicodeError) as error:
        raise ValueError(f"cannot read ASCII Medit mesh: {path}") from error
    return [line.split("#", 1)[0].strip() for line in content.splitlines()
            if line.split("#", 1)[0].strip()]


def _section(lines: list[str], name: str, width: int) -> list[list[str]]:
    locations = [index for index, line in enumerate(lines) if line == name]
    if len(locations) != 1:
        raise ValueError("Medit mesh must have exactly one Vertices and Triangles section")
    start = locations[0]
    if start + 1 >= len(lines):
        raise ValueError(f"Medit {name} section is missing its record count")
    try:
        count = int(lines[start + 1])
    except ValueError as error:
        raise ValueError(f"Medit {name} section has an invalid record count") from error
    if count <= 0:
        raise ValueError(f"Medit {name} section must be nonempty")
    stop = start + 2 + count
    if stop > len(lines) or any(line in {"Vertices", "Triangles", "Tetrahedra", "End"}
                                 for line in lines[start + 2:stop]):
        raise ValueError(f"Medit {name} section is truncated")
    records = [line.split() for line in lines[start + 2:stop]]
    if any(len(record) != width for record in records):
        raise ValueError(f"Medit {name} records must have {width} fields")
    return records


def read_medit_surface(path: str | Path) -> MeditSurface:
    """Read vertices and triangles without repairing, reordering, or perturbing.

    Medit indices are converted from one-based to zero-based. Volume-only
    vertices are omitted in ascending original-index order, and every face
    keeps its original winding and corner order.
    """
    path = Path(path)
    lines = _clean_lines(path)
    dimensions = [line.split() for line in lines if line.split()[0] == "Dimension"]
    if len(dimensions) != 1 or len(dimensions[0]) != 2:
        raise ValueError("Medit mesh is missing Dimension")
    if dimensions[0][1] != "3":
        raise ValueError("only three-dimensional Medit meshes are supported")
    vertex_records = _section(lines, "Vertices", 4)
    triangle_records = _section(lines, "Triangles", 4)
    try:
        vertices_all = np.asarray([[float(value) for value in row[:3]]
                                   for row in vertex_records], dtype=np.float64)
        faces_all = np.asarray([[int(value) - 1 for value in row[:3]]
                                for row in triangle_records], dtype=np.int64)
    except (ValueError, OverflowError) as error:
        raise ValueError("Medit surface contains an invalid numeric record") from error
    if not np.isfinite(vertices_all).all():
        raise ValueError("Medit vertices must be finite")
    vertex_count = len(vertices_all)
    if (faces_all < 0).any() or (faces_all >= vertex_count).any():
        raise ValueError("triangle index outside vertex range")
    if (faces_all[:, 0] == faces_all[:, 1]).any() or \
            (faces_all[:, 1] == faces_all[:, 2]).any() or \
            (faces_all[:, 2] == faces_all[:, 0]).any():
        raise ValueError("Medit surface contains a triangle with repeated vertex indices")

    used = np.unique(faces_all.reshape(-1))
    remap = np.full(vertex_count, -1, dtype=np.int64)
    remap[used] = np.arange(len(used), dtype=np.int64)
    vertices = np.ascontiguousarray(vertices_all[used], dtype=np.float64)
    faces = np.ascontiguousarray(remap[faces_all], dtype=np.int32)
    vertices.flags.writeable = False
    faces.flags.writeable = False
    return MeditSurface(vertices, faces, vertex_count, len(used), len(faces),
                        vertex_count - len(used))
