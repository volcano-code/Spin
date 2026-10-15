#!/usr/bin/env python3
"""Diagnose exact-set mismatches after restoring Cherchi rational points.

This is an audit utility, not a mesh repair or a production Boolean backend.
It compares the native CGAL result with the Cherchi output whose retained
exact rational arrangement points have been rounded once to binary64.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from fractions import Fraction
import json
from pathlib import Path
import sys

import numpy as np


CASES = {
    "one_sheet_hyperboloid": "subtraction",
    "two_sheet_hyperboloid": "intersection",
    "paraboloid": "subtraction",
    "plane": "intersection",
}


def read_obj(path: Path) -> tuple[np.ndarray, np.ndarray]:
    vertices: list[tuple[float, float, float]] = []
    faces: list[tuple[int, int, int]] = []
    for line in path.read_text(encoding="ascii").splitlines():
        fields = line.split()
        if not fields:
            continue
        if fields[0] == "v" and len(fields) == 4:
            vertices.append(tuple(float(value) for value in fields[1:4]))
        elif fields[0] == "f" and len(fields) == 4:
            faces.append(tuple(int(value.split("/")[0]) - 1 for value in fields[1:4]))
        else:
            raise ValueError(f"expected triangular OBJ record in {path}")
    v = np.asarray(vertices, dtype=np.float64)
    f = np.asarray(faces, dtype=np.int32)
    if v.ndim != 2 or v.shape[1:] != (3,) or f.ndim != 2 or f.shape[1:] != (3,):
        raise ValueError(f"invalid triangular OBJ shape in {path}")
    if not np.isfinite(v).all() or (f < 0).any() or (f >= len(v)).any():
        raise ValueError(f"invalid coordinates or indices in {path}")
    return v, f


def _binary_integer(token: str) -> int:
    return (-1 if token.startswith("-") else 1) * int(token.lstrip("-"), 0)


def _rational(token: str) -> Fraction:
    numerator, denominator = token.split("/")
    return Fraction(_binary_integer(numerator), _binary_integer(denominator))


def exact_rounded_vertices(sidecar: Path, reference_vertices: np.ndarray) -> np.ndarray:
    """Map exact-point rows to OBJ vertices by their in-memory coordinate key."""
    by_memory: dict[tuple[float, float, float], tuple[float, float, float]] = {}
    with sidecar.open("r", encoding="ascii", newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            if row["status"] != "MATCHED":
                continue
            memory = tuple(float(row[f"memory_{axis}"]) for axis in "xyz")
            exact = tuple(float(_rational(row[f"exact_{axis}"])) for axis in "xyz")
            if memory in by_memory:
                raise ValueError("sidecar has duplicate in-memory coordinate keys")
            by_memory[memory] = exact
    if set(by_memory) != {tuple(map(float, row)) for row in reference_vertices}:
        raise ValueError("exact-point sidecar does not bijectively cover the Cherchi OBJ")
    return np.asarray([by_memory[tuple(map(float, row))] for row in reference_vertices])


def _face_key(vertices: np.ndarray, face: np.ndarray) -> tuple[tuple[float, float, float], ...]:
    return tuple(sorted(tuple(map(float, vertices[int(index)])) for index in face))


def _source_planes(paths: list[Path]) -> list[tuple[str, np.ndarray, np.ndarray]]:
    planes: list[tuple[str, np.ndarray, np.ndarray]] = []
    for path in paths:
        vertices, faces = read_obj(path)
        for index, face in enumerate(faces):
            a, b, c = vertices[face]
            normal = np.cross(b - a, c - a)
            length = float(np.linalg.norm(normal))
            if length > 0.0:
                planes.append((f"{path.name}:{index}", normal / length, a))
    return planes


def _unmatched_face_plane_audit(
    face_counts: Counter,
    source_planes: list[tuple[str, np.ndarray, np.ndarray]],
) -> dict[str, object]:
    total = 0
    near = 0
    max_best_distance = 0.0
    for face, count in face_counts.items():
        points = np.asarray(face, dtype=np.float64)
        best = min(float(np.max(np.abs((points - anchor) @ normal)))
                   for _, normal, anchor in source_planes)
        total += count
        near += count if best < 1e-14 else 0
        max_best_distance = max(max_best_distance, best)
    return {
        "unmatched_triangle_count": total,
        "unmatched_triangles_within_1e-14_of_an_input_triangle_plane": near,
        "maximum_of_best_input_plane_distances": max_best_distance,
        "meaning": "proximity diagnostic only; not an exact coplanarity or solid-equality proof",
    }


def analyze(reference_dir: Path, sidecar_dir: Path, native_module: object) -> dict[str, object]:
    cases: dict[str, object] = {}
    for name, _operation in CASES.items():
        cgal_v, cgal_f = read_obj(reference_dir / f"cgal-material-{name}.obj")
        cherchi_v, cherchi_f = read_obj(reference_dir / f"cherchi-{name}-run-1.obj")
        rational_v = exact_rounded_vertices(sidecar_dir / f"{name}.obj.exact.tsv", cherchi_v)
        cgal_keys = {tuple(map(float, row)) for row in cgal_v}
        rational_keys = {tuple(map(float, row)) for row in rational_v}

        cgal_faces = Counter(_face_key(cgal_v, face) for face in cgal_f)
        rational_faces = Counter(_face_key(rational_v, face) for face in cherchi_f)
        cgal_only = cgal_faces - rational_faces
        rational_only = rational_faces - cgal_faces
        source_planes = _source_planes([
            reference_dir / "design.obj",
            reference_dir / f"quadric-{name}.obj",
        ])

        cgal_region = native_module.make_exact_material_region(cgal_v.tolist(), cgal_f.tolist())
        rational_region = native_module.make_exact_material_region(rational_v.tolist(), cherchi_f.tolist())
        relation = native_module.relate_exact_material_regions(cgal_region, rational_region)
        volume = relation.get("symmetric_difference_volume_exact")
        cases[name] = {
            "cgal_vertex_count": len(cgal_v),
            "exact_rational_rounded_cherchi_vertex_count": len(rational_v),
            "vertex_coordinate_sets_equal": cgal_keys == rational_keys,
            "shared_vertex_coordinate_count": len(cgal_keys & rational_keys),
            "cgal_face_count": len(cgal_f),
            "cherchi_face_count": len(cherchi_f),
            "shared_unoriented_triangle_count": sum((cgal_faces & rational_faces).values()),
            "cgal_only_unoriented_triangles": _unmatched_face_plane_audit(cgal_only, source_planes),
            "cherchi_only_unoriented_triangles": _unmatched_face_plane_audit(rational_only, source_planes),
            "exact_polyhedral_set_equal": relation.get("equal"),
            "exact_symmetric_difference_volume": volume,
            "exact_symmetric_difference_volume_decimal": (
                float(Fraction(volume)) if isinstance(volume, str) else None
            ),
        }
    return {
        "schema": "cherchi-cgal-triangulation-discrepancy-diagnostic-v1",
        "classification": "SYNTHETIC_REGRESSION",
        "sidecar_point_policy": "correctly-rounded exact rational arrangement points; diagnostic only",
        "face_comparison": "unoriented triangles keyed by exact binary64 coordinate tuples",
        "cases": cases,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-dir", type=Path, required=True)
    parser.add_argument("--sidecar-dir", type=Path, required=True)
    parser.add_argument("--build-directory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    native_dir = args.build_directory.resolve(strict=True) / "python"
    sys.path.insert(0, str(native_dir))
    try:
        import _make_it_stand as native  # type: ignore[import-not-found]
    finally:
        sys.path.pop(0)
    report = analyze(args.reference_dir, args.sidecar_dir, native)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        name: {
            "vertex_coordinate_sets_equal": facts["vertex_coordinate_sets_equal"],
            "shared_triangles": facts["shared_unoriented_triangle_count"],
            "cgal_only": facts["cgal_only_unoriented_triangles"]["unmatched_triangle_count"],
            "cherchi_only": facts["cherchi_only_unoriented_triangles"]["unmatched_triangle_count"],
            "exact_set_equal": facts["exact_polyhedral_set_equal"],
            "symdiff_volume": facts["exact_symmetric_difference_volume_decimal"],
        }
        for name, facts in report["cases"].items()
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
