#!/usr/bin/env python3
"""Evidence-bound same-operand CGAL/Cherchi Boolean cross-check for Faster.

This runner deliberately has a narrow scope.  It generates a deterministic,
clipped ellipsoid from the native Faster evaluator, writes the *actual*
discrete design and quadric meshes as standards-compliant OBJ, then supplies
those identical operands to the pinned Cherchi reference component.  It is a
synthetic method fixture, not an author asset, a replacement Boolean backend,
or a paper-case result.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib
from itertools import product
import json
from fractions import Fraction
from pathlib import Path
import subprocess
import sys
import time
import traceback

import numpy as np

from make_it_stand import spin
from make_it_stand.verification.run_evidence import (
    EXPLICIT_EFFECTIVE_SOURCE_POLICY, begin_run, check_freshness, finish_run,
)
from make_it_stand.verification.cherchi_face_provenance import recover_face_sources
from make_it_stand.verification.cherchi_motion_differential import (
    CherchiDifferentialUnavailable,
    differentiate_polyhedral_raw_moments,
    recover_construction_vertex_motion,
)
from make_it_stand.verification.mass_moments import (
    TriangleMesh,
    compare_exact_rational_and_exported_moments,
    integrate_geometric_moments,
)
from make_it_stand.verification.spin_common import canonical, parse


ROOT = Path(__file__).resolve().parents[3]
_SCHEMA = "faster-reference-boolean-crosscheck-v2"
_TOLERANCE = 1e-12
_SCALE_VOLUME_THRESHOLD = Fraction(1, 2**46)  # 64 * binary64 machine epsilon.
_FIXTURES = (
    ("ellipsoid", (-0.31, 0.07, -0.09, 0.11, 0.83, 1.17, 0.72, 0.06, -0.04, 0.03)),
    ("one_sheet_hyperboloid", (-0.31, 0.07, -0.09, 0.11, 0.83, 1.17, -0.72, 0.06, -0.04, 0.03)),
    ("two_sheet_hyperboloid", (0.31, 0.07, -0.09, 0.11, 0.83, 1.17, -0.72, 0.06, -0.04, 0.03)),
    ("plane", (0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)),
    ("paraboloid", (0.0, 0.0, 0.0, -1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0)),
)


def _write(path: Path, value: object) -> None:
    with path.open("x", encoding="ascii") as stream:
        stream.write(canonical(value) + "\n")


def _box(lo: tuple[float, float, float], hi: tuple[float, float, float]) -> tuple[np.ndarray, np.ndarray]:
    x0, y0, z0 = lo
    x1, y1, z1 = hi
    vertices = np.asarray(((x0, y0, z0), (x1, y0, z0), (x1, y1, z0), (x0, y1, z0),
                           (x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1)), dtype=float)
    faces = np.asarray(((0, 2, 1), (0, 3, 2), (4, 5, 6), (4, 6, 7),
                        (0, 1, 5), (0, 5, 4), (1, 2, 6), (1, 6, 5),
                        (2, 3, 7), (2, 7, 6), (3, 0, 4), (3, 4, 7)), dtype=int)
    return vertices, faces


def _write_obj(path: Path, vertices: tuple[tuple[float, ...], ...] | np.ndarray,
               faces: tuple[tuple[int, ...], ...] | np.ndarray) -> None:
    """Write a portable OBJ; normalize signed zero without changing geometry."""
    with path.open("x", encoding="ascii", newline="\n") as stream:
        for row in vertices:
            # ``-0`` caused an earlier malformed diagnostic to obscure the
            # real issue.  Canonical +0 is numerically identical and leaves
            # connectivity, sampling and all nonzero coordinates untouched.
            values = (0.0 if float(value) == 0.0 else float(value) for value in row)
            stream.write("v " + " ".join(f"{value:.17g}" for value in values) + "\n")
        for row in faces:
            stream.write("f " + " ".join(str(int(value) + 1) for value in row) + "\n")


def _read_obj(path: Path) -> tuple[np.ndarray, np.ndarray]:
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
            raise ValueError("reference Boolean output is not a triangular OBJ")
    if not vertices or not faces:
        raise ValueError("reference Boolean produced an empty OBJ")
    return np.asarray(vertices, dtype=float), np.asarray(faces, dtype=np.int32)


def _triangle_boundary_signature(vertices, faces) -> str:
    """Hash the oriented triangle boundary independent of vertex/face order.

    Coordinates use the same 17-significant-digit, round-trip-safe serialization as the pinned
    Cherchi adapter. The signature compares actual triangles and winding, not
    only moments, bounds, or mesh counts.
    """
    coordinates = []
    for row in vertices:
        values = []
        for value in row:
            scalar = float(value)
            if not np.isfinite(scalar):
                raise ValueError("nonfinite vertex in Boolean boundary signature")
            values.append(format(0.0 if scalar == 0.0 else scalar, ".17g"))
        coordinates.append(tuple(values))
    triangles = []
    for face in faces:
        if len(face) != 3 or len(set(int(index) for index in face)) != 3 or any(
                int(index) < 0 or int(index) >= len(coordinates) for index in face):
            raise ValueError("invalid triangle in Boolean boundary signature")
        triangle = [coordinates[int(index)] for index in face]
        pivot = min(range(3), key=triangle.__getitem__)
        triangles.append(tuple(triangle[pivot:] + triangle[:pivot]))
    payload = json.dumps(sorted(triangles), separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def _exact_solid_relation_status(facts: dict[str, object]) -> str:
    """Classify exact represented-polyhedral-set comparison facts."""
    if facts.get("completed") is not True:
        return "UNKNOWN"
    if facts.get("equal") is True and facts.get("a_subset_b") is True and facts.get("b_subset_a") is True:
        return "SAME_EXACT_POLYHEDRAL_SOLID"
    if facts.get("equal") is False:
        return "DIFFERENT_EXACT_POLYHEDRAL_SOLIDS"
    return "UNKNOWN"


def _symmetric_difference_volume_comparison(
    facts: dict[str, object], vertices_a: object, vertices_b: object,
) -> dict[str, object]:
    """Report exact symmetric-difference volume against a fixed scale policy.

    The threshold is a diagnostic tolerance (2^-46 times the cube of the
    combined meshes' largest axis span), not a Hausdorff or topology bound and
    never upgrades exact-set inequality to equality.
    """
    raw = facts.get("symmetric_difference_volume_exact")
    if facts.get("completed") is not True or not isinstance(raw, str):
        return {"status": "UNKNOWN", "exact_set_equal": facts.get("equal"),
                "symmetric_difference_volume_exact": raw}
    a = np.asarray(vertices_a, dtype=float)
    b = np.asarray(vertices_b, dtype=float)
    if a.ndim != 2 or b.ndim != 2 or a.shape[1:] != (3,) or b.shape[1:] != (3,) or \
            len(a) == 0 or len(b) == 0 or not np.isfinite(a).all() or not np.isfinite(b).all():
        return {"status": "UNKNOWN", "exact_set_equal": facts.get("equal"),
                "symmetric_difference_volume_exact": raw}
    exact = Fraction(raw)
    coords = [[Fraction(float(value)) for value in row] for row in np.vstack((a, b))]
    spans = [max(point[axis] for point in coords) - min(point[axis] for point in coords)
             for axis in range(3)]
    scale = max(spans)
    if scale <= 0:
        return {"status": "UNKNOWN", "exact_set_equal": facts.get("equal"),
                "symmetric_difference_volume_exact": raw}
    normalized = exact / (scale ** 3)
    tolerance = _SCALE_VOLUME_THRESHOLD
    if facts.get("equal") is True and exact == 0:
        status = "EXACT_SET_MATCH_ZERO_VOLUME_DIFFERENCE"
    elif exact <= tolerance * (scale ** 3):
        status = "WITHIN_INPUT_SCALE_VOLUME_THRESHOLD_NOT_EXACT_MATCH"
    else:
        status = "OUTSIDE_INPUT_SCALE_VOLUME_THRESHOLD"
    return {
        "status": status,
        "exact_set_equal": facts.get("equal"),
        "symmetric_difference_volume_exact": raw,
        "symmetric_difference_volume_decimal": float(exact),
        "normalizing_length_exact": f"{scale.numerator}/{scale.denominator}",
        "normalized_symmetric_difference_volume_exact": f"{normalized.numerator}/{normalized.denominator}",
        "normalized_symmetric_difference_volume_decimal": float(normalized),
        "threshold_policy": "2^-46 * (combined-bounds-largest-axis-span)^3; diagnostic only",
        "threshold_normalized_exact": f"{tolerance.numerator}/{tolerance.denominator}",
        "within_threshold": normalized <= tolerance,
    }


def _parse_exact_binary_rational(token: str) -> Fraction:
    if token == "0/1":
        return Fraction(0)
    negative = token.startswith("-")
    value = token[1:] if negative else token
    numerator, denominator = value.split("/")

    def integer(part: str) -> int:
        return int(part[2:], 2) if part.startswith("0b") else int(part)

    result = Fraction(integer(numerator), integer(denominator))
    return -result if negative else result


def _fixed_box_vertex_audit(exact_rows: list[list[str]], design_vertices: np.ndarray) -> dict[str, object]:
    """Compare exact construction and exported vertices against this fixture's box bounds.

    Vertex bounds are a diagnostic for a fixed convex box, not an exact
    arbitrary-domain material containment or Boolean-equivalence oracle.
    """
    points = [tuple(Fraction(float(value)) for value in row) for row in design_vertices]
    bounds = tuple((min(point[axis] for point in points), max(point[axis] for point in points))
                   for axis in range(3))
    corners = set(product(*(bound for bound in bounds)))
    if len(points) != 8 or set(points) != corners:
        return {"status": "UNKNOWN", "scope": "AXIS_ALIGNED_FIXED_DESIGN_BOX_VERTICES_ONLY",
                "reason": "design vertices are not the eight corners of a box"}
    exact_outside, exported_outside = [], []
    for index, row in enumerate(exact_rows):
        if len(row) != 11 or int(row[0]) != index:
            return {"status": "UNKNOWN", "scope": "AXIS_ALIGNED_FIXED_DESIGN_BOX_VERTICES_ONLY",
                    "reason": "exact-coordinate sidecar is incomplete or reordered"}
        exact = tuple(_parse_exact_binary_rational(row[2 + axis]) for axis in range(3))
        exported = tuple(Fraction(float(row[8 + axis])) for axis in range(3))
        if any(exact[axis] < bounds[axis][0] or exact[axis] > bounds[axis][1]
               for axis in range(3)):
            exact_outside.append(index)
        if any(exported[axis] < bounds[axis][0] or exported[axis] > bounds[axis][1]
               for axis in range(3)):
            exported_outside.append(index)
    return {"status": "DIAGNOSTIC_ONLY", "scope": "AXIS_ALIGNED_FIXED_DESIGN_BOX_VERTICES_ONLY",
            "vertex_count": len(exact_rows),
            "exact_rational_outside_count": len(exact_outside),
            "exported_binary64_outside_count": len(exported_outside),
            "exact_rational_outside_indices": exact_outside,
            "exported_binary64_outside_indices": exported_outside}


def _exact_construction_relation(extension: object, exact_rows: list[list[str]] | None,
                                 exported_vertices: np.ndarray, exported_faces: np.ndarray,
                                 cgal_region: object) -> tuple[dict[str, object], object | None]:
    """Compare Cherchi's construction rationals, not its rounded OBJ vertices.

    This only relates two represented polyhedral sets. It does not prove a
    continuous quadric Boolean, CAD identity, or a usable motion derivative.
    """
    result: dict[str, object] = {
        "scope": "cherchi-exact-construction-vs-cgal-exported-polyhedron-only",
        "status": "UNKNOWN", "construction_region_valid": False,
        "completed": False, "equal": None,
    }
    if exact_rows is None or len(exact_rows) != len(exported_vertices):
        result["reason"] = "construction sidecar is absent or has a different vertex count"
        return result, None
    try:
        rational_vertices = []
        for index, row in enumerate(exact_rows):
            if len(row) != 11 or int(row[0]) != index:
                raise ValueError("construction sidecar is incomplete or reordered")
            if any(float(row[8 + axis]) != float(exported_vertices[index, axis])
                   for axis in range(3)):
                raise ValueError("construction sidecar does not bind the exported OBJ vertex")
            exact = [_parse_exact_binary_rational(row[2 + axis]) for axis in range(3)]
            rational_vertices.append([f"{value.numerator}/{value.denominator}" for value in exact])
        region = extension.make_exact_material_region_from_rationals(
            rational_vertices, exported_faces.tolist())
        result["construction_region_digest"] = region.facts["region_digest"]
        result["construction_region_valid"] = bool(region.facts["valid"])
        if not result["construction_region_valid"]:
            result["reason"] = region.facts["error"]
            return result, None
        relation = extension.relate_exact_material_regions(region, cgal_region)
        result["completed"] = bool(relation["completed"])
        if not result["completed"]:
            result["reason"] = relation["error"]
            return result, None
        result["equal"] = bool(relation["equal"])
        result["a_subset_b"] = bool(relation["a_subset_b"])
        result["b_subset_a"] = bool(relation["b_subset_a"])
        result["status"] = ("SAME_EXACT_POLYHEDRAL_SOLID" if result["equal"]
                            else "DIFFERENT_EXACT_POLYHEDRAL_SOLIDS")
    except (ValueError, RuntimeError, KeyError, IndexError) as error:
        result["reason"] = str(error)
        return result, None
    return result, region


def _summarize_vertex_construction_provenance(
    input_meshes: tuple[Path, ...], output_mesh: Path, sidecar_path: Path,
) -> dict[str, object]:
    """Map retained Cherchi construction operands to exact input vertex/face IDs.

    This is a fixed-arrangement provenance diagnostic, not a motion
    differential. Duplicate exact coordinates/facets remain explicitly
    ambiguous so later consumers cannot infer a unique moving source.
    """
    output_vertices, _ = _read_obj(output_mesh)
    lines = sidecar_path.read_text(encoding="ascii").splitlines()
    header = lines[0].split("\t") if lines else []
    expected_header = ["vertex_index", "point_type", "operand_count", "parameter_t"] + [
        f"operand{index}_{axis}" for index in range(9) for axis in "xyz"
    ]
    if header != expected_header or len(lines) - 1 != len(output_vertices):
        raise ValueError("Cherchi construction sidecar header/count does not match output mesh")

    source_vertices: dict[tuple[Fraction, Fraction, Fraction], list[dict[str, int]]] = {}
    source_facets: dict[tuple[tuple[Fraction, Fraction, Fraction], ...], list[dict[str, int]]] = {}
    for object_index, path in enumerate(input_meshes):
        vertices, faces = _read_obj(path)
        exact_vertices = [tuple(Fraction(float(value)) for value in row) for row in vertices]
        for vertex_index, point in enumerate(exact_vertices):
            source_vertices.setdefault(point, []).append({
                "input_object": object_index, "input_vertex": vertex_index,
            })
        for face_index, face in enumerate(faces):
            key = tuple(sorted(exact_vertices[int(index)] for index in face))
            source_facets.setdefault(key, []).append({
                "input_object": object_index, "input_triangle": face_index,
            })

    expected_operands = {"EXPLICIT3D": 1, "LPI": 5, "TPI": 9, "LNC": 2}
    point_type_counts: dict[str, int] = {}
    operand_source_counts = {"UNIQUE": 0, "AMBIGUOUS": 0, "UNMATCHED": 0}
    plane_source_counts = {"UNIQUE": 0, "AMBIGUOUS": 0, "UNMATCHED": 0}
    records: list[dict[str, object]] = []
    for expected_index, line in enumerate(lines[1:]):
        fields = line.split("\t")
        if len(fields) != len(expected_header) or int(fields[0]) != expected_index:
            raise ValueError("Cherchi construction sidecar row is malformed or reordered")
        point_type = fields[1]
        operand_count = int(fields[2])
        if point_type not in expected_operands or operand_count != expected_operands[point_type]:
            raise ValueError("Cherchi construction sidecar contains unsupported point construction")
        parameter = float(fields[3])
        if not np.isfinite(parameter):
            raise ValueError("Cherchi construction sidecar has a nonfinite parameter")
        if point_type == "LNC" and not 0.0 <= parameter <= 1.0:
            raise ValueError("Cherchi LNC parameter is outside the convex interpolation interval")
        operands = []
        for operand_index in range(9):
            coordinate_tokens = fields[4 + 3 * operand_index:7 + 3 * operand_index]
            if operand_index >= operand_count:
                if coordinate_tokens != ["-", "-", "-"]:
                    raise ValueError("unused Cherchi construction operand slot is not empty")
                continue
            point = tuple(_parse_exact_binary_rational(token) for token in coordinate_tokens)
            candidates = source_vertices.get(point, [])
            status = "UNMATCHED" if not candidates else "UNIQUE" if len(candidates) == 1 else "AMBIGUOUS"
            operand_source_counts[status] += 1
            operands.append({"source_vertices": candidates, "status": status})

        if point_type == "LPI":
            plane_groups = [(2, 3, 4)]
        elif point_type == "TPI":
            plane_groups = [(0, 1, 2), (3, 4, 5), (6, 7, 8)]
        else:
            plane_groups = []
        plane_sources = []
        for group in plane_groups:
            points = tuple(sorted(tuple(_parse_exact_binary_rational(fields[4 + 3 * index + axis])
                                        for axis in range(3)) for index in group))
            candidates = source_facets.get(points, [])
            status = "UNMATCHED" if not candidates else "UNIQUE" if len(candidates) == 1 else "AMBIGUOUS"
            plane_source_counts[status] += 1
            plane_sources.append({"source_triangles": candidates, "status": status})
        point_type_counts[point_type] = point_type_counts.get(point_type, 0) + 1
        operand_coordinates = []
        for operand_index in range(operand_count):
            coordinate_tokens = fields[4 + 3 * operand_index:7 + 3 * operand_index]
            exact_point = tuple(_parse_exact_binary_rational(token) for token in coordinate_tokens)
            operand_coordinates.append([float(value) for value in exact_point])
        records.append({
            "vertex_index": expected_index,
            "point_type": point_type,
            # Cherchi's LNC construction stores a fixed interpolation weight
            # in the pinned implicit point. Retain it so a fixed-stratum
            # differential can reproduce the exact construction semantics.
            "parameter_t": parameter,
            "operands": operands,
            "operand_coordinates": operand_coordinates,
            "surface_planes": plane_sources,
        })

    all_operands_are_exact_input_vertices = operand_source_counts["UNMATCHED"] == 0
    all_surface_planes_have_unique_sources = (
        plane_source_counts["AMBIGUOUS"] == 0 and plane_source_counts["UNMATCHED"] == 0
    )
    if operand_source_counts["UNMATCHED"] or plane_source_counts["UNMATCHED"]:
        source_identity_status = "UNMATCHED"
    elif operand_source_counts["AMBIGUOUS"] or plane_source_counts["AMBIGUOUS"]:
        source_identity_status = "AMBIGUOUS"
    else:
        source_identity_status = "UNIQUE"
    return {
        "schema": "cherchi-vertex-construction-source-diagnostic-v1",
        "classification": "SYNTHETIC_REGRESSION",
        "status": "DIAGNOSTIC_ONLY",
        "meaning": "exact construction operands mapped to candidate input vertex/face IDs; not a motion derivative",
        "vertex_count": len(records),
        "point_type_counts": dict(sorted(point_type_counts.items())),
        "operand_source_counts": operand_source_counts,
        "surface_plane_source_counts": plane_source_counts,
        "source_identity_status": source_identity_status,
        "all_operands_are_exact_input_vertices": all_operands_are_exact_input_vertices,
        "all_surface_planes_have_unique_sources": all_surface_planes_have_unique_sources,
        "records": records,
    }


def _independent_raw_moments(vertices: np.ndarray, faces: np.ndarray) -> np.ndarray:
    moments = integrate_geometric_moments(TriangleMesh(
        tuple(tuple(float(value) for value in row) for row in vertices),
        tuple(tuple(int(value) for value in face) for face in faces),
    ))
    second = moments.second_moment
    return np.asarray((
        moments.volume, *moments.first_moment,
        second[0][0], second[1][1], second[2][2],
        second[1][2], second[0][2], second[0][1],
    ), dtype=np.float64)


def _construction_topology_signature(summary: dict[str, object], faces: np.ndarray) -> str:
    records = summary.get("records")
    if not isinstance(records, list):
        raise ValueError("construction provenance lacks vertex records")
    descriptors: list[str] = []
    for index, record in enumerate(records):
        if record.get("vertex_index") != index:
            raise ValueError("construction provenance vertices are reordered")
        operand_ids = []
        for operand in record.get("operands", []):
            candidates = operand.get("source_vertices", [])
            if operand.get("status") != "UNIQUE" or len(candidates) != 1:
                raise ValueError("topology signature requires unique construction source vertices")
            candidate = candidates[0]
            operand_ids.append((candidate["input_object"], candidate["input_vertex"]))
        plane_ids = []
        for plane in record.get("surface_planes", []):
            candidates = plane.get("source_triangles", [])
            if plane.get("status") != "UNIQUE" or len(candidates) != 1:
                raise ValueError("topology signature requires unique source planes")
            candidate = candidates[0]
            plane_ids.append((candidate["input_object"], candidate["input_triangle"]))
        # LNC's fixed interpolation parameter is part of the construction
        # identity; identical source operands with different t are different
        # represented points and cannot share a motion-cache identity.
        parameter = record.get("parameter_t") if record.get("point_type") == "LNC" else None
        descriptors.append(json.dumps((record.get("point_type"), operand_ids, plane_ids, parameter),
                                      separators=(",", ":")))
    if len(set(descriptors)) != len(descriptors):
        raise ValueError("construction provenance does not uniquely identify output vertices")
    triangles = []
    for face in faces:
        indices = tuple(int(index) for index in face)
        if len(indices) != 3 or min(indices) < 0 or max(indices) >= len(descriptors):
            raise ValueError("Boolean output face is invalid")
        triangles.append(tuple(sorted(descriptors[index] for index in indices)))
    payload = json.dumps(sorted(triangles), separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def _cherchi_motion_finite_difference_audit(
    *, fixture_name: str, extension, probe: Path, operation: str, design_path: Path,
    quadric_path: Path, base_output_path: Path, coefficients: tuple[float, ...],
    request: spin.SpinItFasterRequest, base_evaluation: spin.SpinQuadricEvaluation,
    base_constructions: dict[str, object], output_directory: Path,
) -> dict[str, object]:
    """Compare recovered source-motion transport against the pinned Boolean map.

    This is intentionally a diagnostic artifact, not an optimizer result. Every
    perturbed Boolean must preserve its provenance-labeled topology; otherwise
    no Jacobian comparison is reported for that trial.
    """
    subdivisions, max_radius_doublings = 2, 4
    settings = extension.SpinTriangulationSettings()
    settings.icosphere_subdivisions = subdivisions
    settings.maximum_radius_doublings = max_radius_doublings
    base_source = extension.evaluate_spin_it_faster_source_motion(
        np.asarray(coefficients, dtype=np.float64),
        np.asarray(request.design_vertices, dtype=np.float64), settings)
    if not base_source["defined"]:
        return {"status": "SOURCE_MOTION_UNAVAILABLE", "reason": base_source["diagnostic"]}
    source_vertices = np.asarray(base_source["vertices"], dtype=np.float64)
    source_faces = np.asarray(base_source["faces"], dtype=np.int64)
    source_motion = np.asarray(base_source["vertex_velocity"], dtype=np.float64).reshape(-1, 3, 10)
    evaluation_vertices = np.asarray(base_evaluation.interface_vertices, dtype=np.float64)
    evaluation_faces = np.asarray(base_evaluation.interface_faces, dtype=np.int64)
    if source_faces.shape != evaluation_faces.shape or not np.array_equal(source_faces, evaluation_faces) or \
            source_vertices.shape != evaluation_vertices.shape or \
            not np.allclose(source_vertices, evaluation_vertices, rtol=0.0, atol=2e-13):
        return {"status": "SOURCE_MESH_MISMATCH", "source_triangulation_hash": base_source["triangulation_hash"]}

    base_vertices, base_faces = _read_obj(base_output_path)
    try:
        output_motion = recover_construction_vertex_motion(
            base_constructions["records"], request.design_vertices, source_vertices,
            source_motion, base_vertices)
        analytic_jacobian = differentiate_polyhedral_raw_moments(
            base_vertices, base_faces, output_motion, density=request.material_density)
        base_topology_signature = _construction_topology_signature(base_constructions, base_faces)
    except (CherchiDifferentialUnavailable, ValueError, KeyError, TypeError) as error:
        return {"status": "PROVENANCE_MOTION_UNAVAILABLE", "reason": str(error),
                "source_triangulation_hash": base_source["triangulation_hash"]}

    normalized = np.asarray(coefficients, dtype=np.float64)
    normalized /= np.linalg.norm(normalized)
    _, _, vh = np.linalg.svd(normalized.reshape(1, 10), full_matrices=True)
    directions = vh[1:]
    steps = (1e-4, 3e-5, 1e-5)
    trials: list[dict[str, object]] = []
    same_topology = True
    for direction_index, direction in enumerate(directions):
        exact_directional = analytic_jacobian @ direction
        for step in steps:
            perturbed_moments: dict[int, np.ndarray] = {}
            topology_hashes = {}
            for sign in (-1, 1):
                perturbed_coefficients = normalized + sign * step * direction
                perturbed_source = extension.evaluate_spin_it_faster_source_motion(
                    perturbed_coefficients, np.asarray(request.design_vertices, dtype=np.float64), settings)
                if not perturbed_source["defined"]:
                    return {"status": "PERTURBED_SOURCE_MOTION_UNAVAILABLE", "reason": perturbed_source["diagnostic"],
                            "direction_index": direction_index, "step": step,
                            "completed_trials": trials}
                vertices = np.asarray(perturbed_source["vertices"], dtype=np.float64)
                faces = np.asarray(perturbed_source["faces"], dtype=np.int64)
                if not np.array_equal(faces, source_faces) or \
                        perturbed_source["truncation_radius"] != base_source["truncation_radius"]:
                    return {"status": "QUADRIC_SOURCE_BRANCH_CHANGED", "direction_index": direction_index,
                            "step": step, "sign": sign, "completed_trials": trials}
                moved_source_path = output_directory / f"motion-{fixture_name}-q{direction_index}-h{step:g}-s{sign:+d}.obj"
                result_path = output_directory / f"motion-out-{fixture_name}-q{direction_index}-h{step:g}-s{sign:+d}.obj"
                _write_obj(moved_source_path, vertices, faces)
                completed = subprocess.run(
                    [str(probe), operation, str(design_path), str(moved_source_path),
                     str(result_path), "--dump-actual-provenance"],
                    cwd=output_directory, capture_output=True, text=True, check=False)
                if completed.returncode != 0:
                    return {"status": "CHERCHI_PERTURBATION_FAILED", "reason": completed.stderr.strip(),
                            "direction_index": direction_index, "step": step, "sign": sign,
                            "completed_trials": trials}
                output_vertices, output_faces = _read_obj(result_path)
                construction = _summarize_vertex_construction_provenance(
                    (design_path, moved_source_path), result_path,
                    Path(str(result_path) + ".construction.tsv"))
                try:
                    signature = _construction_topology_signature(construction, output_faces)
                except (ValueError, KeyError, TypeError) as error:
                    return {"status": "PERTURBED_TOPOLOGY_IDENTITY_UNAVAILABLE", "reason": str(error),
                            "direction_index": direction_index, "step": step, "sign": sign,
                            "completed_trials": trials}
                topology_hashes[sign] = signature
                if signature != base_topology_signature:
                    same_topology = False
                perturbed_moments[sign] = _independent_raw_moments(output_vertices, output_faces)
            if not same_topology or topology_hashes[-1] != topology_hashes[1]:
                return {"status": "TOPOLOGY_CHANGED", "direction_index": direction_index,
                        "step": step, "topology_hashes": {str(key): value for key, value in topology_hashes.items()},
                        "completed_trials": trials}
            finite_difference = (perturbed_moments[1] - perturbed_moments[-1]) / (2.0 * step)
            relative_error = float(np.linalg.norm(finite_difference - exact_directional) /
                                   max(np.linalg.norm(finite_difference), np.linalg.norm(exact_directional), 1e-14))
            trials.append({"direction_index": direction_index, "direction": direction.tolist(),
                           "step": step, "relative_error": relative_error,
                           "topology_sha256": topology_hashes[1]})
    maximum_relative_error = max(float(trial["relative_error"]) for trial in trials)
    return {
        "status": "PASS" if maximum_relative_error <= 1e-6 else "DIFFERENTIAL_MISMATCH",
        "classification": "SYNTHETIC_REGRESSION",
        "derivative_mapping": "Cherchi construction provenance -> native source motion -> Reynolds polyhedral transport",
        "native_source_motion_hash": base_source["scientific_hash"],
        "source_triangulation_hash": base_source["triangulation_hash"],
        "source_vertex_count": len(source_vertices), "boolean_vertex_count": len(base_vertices),
        "boolean_face_count": len(base_faces), "tested_direction_count": len(directions),
        "step_sizes": list(steps), "maximum_relative_error": maximum_relative_error,
        "threshold": 1e-6, "trials": trials,
    }


def _native_from_build(build: Path):
    directory = (build / "python").resolve(strict=True)
    candidates = sorted(directory.glob("_make_it_stand*.so"))
    if len(candidates) != 1:
        raise ValueError("declared build does not contain exactly one native extension")
    sys.path.insert(0, str(directory))
    try:
        extension = importlib.import_module("_make_it_stand")
    finally:
        sys.path.pop(0)
    if Path(getattr(extension, "__file__", "")).resolve() != candidates[0].resolve():
        raise ValueError("loaded extension differs from declared build directory")
    return extension


def _operand_necessary_conditions(extension, output_region, design_region,
                                  quadric_region, operation: str) -> dict[str, object]:
    """Check necessary exact-set conditions on a represented Boolean mesh.

    Passing these conditions is not sufficient to establish equality with the
    exact Boolean of the input operands.
    """
    if operation not in {"intersection", "subtraction"}:
        raise ValueError("unsupported reference Boolean operation")
    digests = {
        "output_digest": output_region.facts["region_digest"],
        "design_digest": design_region.facts["region_digest"],
        "quadric_digest": quadric_region.facts["region_digest"],
    }
    if not all(region.facts["valid"] for region in
               (output_region, design_region, quadric_region)):
        return {**digests, "operation": operation, "completed": False,
                "subset_of_design": None, "quadric_operation_condition": None,
                "status": "UNKNOWN", "reason": "INVALID_EXACT_REGION"}
    try:
        design_check = extension.check_exact_material_subset(output_region, design_region)
        if operation == "intersection":
            quadric_check = extension.check_exact_material_subset(output_region, quadric_region)
            quadric_condition = quadric_check["a_subset_b"] if quadric_check["completed"] else None
            condition_semantics = "output_subset_of_quadric_operand"
        else:
            quadric_check = extension.relate_exact_material_regions(output_region, quadric_region)
            overlap = quadric_check["positive_volume_overlap"] if quadric_check["completed"] else None
            quadric_condition = None if overlap is None else not overlap
            condition_semantics = "no_positive_volume_overlap_with_quadric_operand"
        subset = design_check["a_subset_b"] if design_check["completed"] else None
        completed = subset is not None and quadric_condition is not None
        status = ("PASS" if subset and quadric_condition else
                  "FAIL" if subset is False or quadric_condition is False else "UNKNOWN")
        return {**digests, "operation": operation, "completed": completed,
                "subset_of_design": subset,
                "quadric_operation_condition": quadric_condition,
                "quadric_condition_semantics": condition_semantics,
                "status": status, "design_subset_facts": dict(design_check),
                "quadric_operand_facts": dict(quadric_check)}
    except Exception as error:
        return {**digests, "operation": operation, "completed": False,
                "subset_of_design": None, "quadric_operation_condition": None,
                "status": "UNKNOWN", "reason": str(error)}


def _fixtures() -> tuple[spin.SpinItFasterRequest, tuple[tuple[str, tuple[float, ...], spin.SpinQuadricEvaluation], ...]]:
    outer_vertices, outer_faces = _box((-2.0, -2.0, -2.0), (2.0, 2.0, 2.0))
    design_vertices, design_faces = _box((-1.0, -1.0, -1.0), (1.0, 1.0, 1.0))
    request = spin.SpinItFasterRequest(
        outer_vertices.tolist(), outer_faces.tolist(),
        design_vertices.tolist(), design_faces.tolist(),
        mode="YO_YO",
        settings=spin.SpinSolverConfig(equation_profile="DUAL_MONOMIAL_RAW_V1"),
    )
    settings = spin.SpinEvaluatorConfig(
        derivative_method="DISCRETE_ELLIPSOID", icosphere_subdivisions=2,
        maximum_radius_doublings=4,
    )
    records = []
    for name, coefficients in _FIXTURES:
        evaluation = spin.evaluate_spin_it_faster_quadric(
            request, coefficients, evaluator_settings=settings,
        )
        if len(evaluation.interface_faces) == 0 or len(evaluation.material_faces) <= len(request.design_faces):
            raise RuntimeError(name + " fixture is not an actual clipped discrete-quadric Boolean")
        records.append((name, coefficients, evaluation))
    return request, tuple(records)


def run(output_directory: Path, *, build_directory: Path, cherchi_probe: Path,
        dependency_roots: dict[str, Path]) -> dict[str, object]:
    output = output_directory.absolute()
    output.mkdir(parents=False, exist_ok=False)
    build = build_directory.resolve()
    probe = cherchi_probe.resolve(strict=True)
    if not probe.is_file() or not probe.stat().st_mode & 0o111:
        raise ValueError("pinned Cherchi reference probe is unavailable or non-executable")
    extension = _native_from_build(build)
    request, fixtures = _fixtures()
    design_path = output / "design.obj"
    _write_obj(design_path, request.design_vertices, request.design_faces)
    design_vertices_on_disk, design_faces_on_disk = _read_obj(design_path)
    design_region = extension.make_exact_material_region(
        design_vertices_on_disk.tolist(), design_faces_on_disk.tolist())
    quadric_paths: dict[str, Path] = {}
    cgal_material_paths: dict[str, Path] = {}
    for name, _, evaluation in fixtures:
        quadric_paths[name] = output / f"quadric-{name}.obj"
        _write_obj(quadric_paths[name], evaluation.interface_vertices, evaluation.interface_faces)
        cgal_material_paths[name] = output / f"cgal-material-{name}.obj"
        _write_obj(cgal_material_paths[name], evaluation.material_vertices, evaluation.material_faces)
    started, tick, begin, events = datetime.now(timezone.utc).isoformat(), time.perf_counter(), None, []
    try:
        inputs = parse(canonical({
            "schema": _SCHEMA + "-inputs",
            "classification": "SYNTHETIC_REGRESSION",
            "fixtures": [{
                "name": name, "coefficients": list(coefficients),
                "native_geometry_hash": evaluation.geometry_hash,
                "native_topology_hash": evaluation.topology_hash,
                "native_operation": "intersection" if evaluation.quadric_solid_is_positive_region else "subtraction",
                "native_variable_raw_moments": list(evaluation.variable_raw_moments),
            } for name, coefficients, evaluation in fixtures],
            "cherchi_probe": {"path": str(probe), "sha256": hashlib.sha256(probe.read_bytes()).hexdigest()},
            "comparison_tolerance": _TOLERANCE,
        }))
        input_path = output / "inputs.json"
        _write(input_path, inputs)
        begin = begin_run(
            repository_root=ROOT,
            runner_files=[Path(__file__),
                          ROOT / "python/make_it_stand/verification/cherchi_face_provenance.py"],
            checker_files=[
                ROOT / "experiments/cherchi_reference_probe/probe.cpp",
                ROOT / "experiments/cherchi_reference_probe/CMakeLists.txt",
                ROOT / "src/spin/QuadricMomentEvaluator.cpp",
                ROOT / "src/verification/ExactBoundaryAudit.cpp",
                ROOT / "src/verification/ExactMaterialRegion.cpp",
                ROOT / "src/verification/ExactLibiglBoolean.cpp",
                ROOT / "include/verification/ExactLibiglBoolean.h",
                ROOT / "include/verification/ExactRationalCoordinates.h",
                ROOT / "bindings/pybind_module.cpp",
                ROOT / "libigl/include/igl/copyleft/cgal/mesh_boolean.cpp",
            ],
            asset_files=[input_path, design_path, *quadric_paths.values(),
                         *cgal_material_paths.values(), probe],
            build_configuration_files=[ROOT / "CMakeLists.txt", build / "CMakeCache.txt"],
            parameters=inputs,
            seed_protocol={"seed": 0, "protocol": "deterministic same-operand clipped quadric type matrix"},
            dependency_roots=dependency_roots,
            native_modules=[extension],
            pure_python=False,
            source_policy=EXPLICIT_EFFECTIVE_SOURCE_POLICY,
        )
        _write(output / "begin.json", begin.as_dict())
        scientific_fixtures = []
        for name, coefficients, evaluation in fixtures:
            operation = "intersection" if evaluation.quadric_solid_is_positive_region else "subtraction"
            runs = []
            expected = np.asarray(evaluation.variable_raw_moments, dtype=float)
            for repeat in (1, 2):
                result_path = output / f"cherchi-{name}-run-{repeat}.obj"
                completed = subprocess.run(
                    [str(probe), operation, str(design_path), str(quadric_paths[name]), str(result_path)],
                    cwd=output, capture_output=True, text=True, check=False,
                )
                if completed.returncode != 0:
                    raise RuntimeError(name + " Cherchi reference Boolean failed: " + completed.stderr.strip())
                vertices, faces = _read_obj(result_path)
                moments = np.asarray(spin.integrate_raw_mass_moments(vertices.tolist(), faces.tolist()), dtype=float)
                maximum_absolute_error = float(np.max(np.abs(moments - expected)))
                if maximum_absolute_error > _TOLERANCE:
                    raise RuntimeError(name + " independent Boolean raw moments exceed declared comparison tolerance")
                runs.append({
                    "path": str(result_path), "sha256": hashlib.sha256(result_path.read_bytes()).hexdigest(),
                    "stdout": completed.stdout.strip(), "vertex_count": int(vertices.shape[0]),
                    "face_count": int(faces.shape[0]), "raw_moments": [float(value) for value in moments],
                    "maximum_absolute_error": maximum_absolute_error,
                })
            if runs[0]["sha256"] != runs[1]["sha256"] or runs[0]["raw_moments"] != runs[1]["raw_moments"]:
                raise RuntimeError(name + " repeated Cherchi Boolean output is nondeterministic")
            provenance_path = output / f"cherchi-{name}-provenance.obj"
            provenance_run = subprocess.run(
                [str(probe), operation, str(design_path), str(quadric_paths[name]),
                 str(provenance_path), "--dump-actual-provenance"],
                cwd=output, capture_output=True, text=True, check=False,
            )
            if provenance_run.returncode != 0:
                raise RuntimeError(name + " Cherchi exact provenance diagnostic failed: " +
                                   provenance_run.stderr.strip())
            provenance_geometry_hash = hashlib.sha256(provenance_path.read_bytes()).hexdigest()
            provenance_summary: dict[str, object]
            construction_summary: dict[str, object]
            exact_rows: list[list[str]] | None = None
            if provenance_geometry_hash != runs[0]["sha256"]:
                provenance_summary = {
                    "status": "DIAGNOSTIC_OUTPUT_MISMATCH",
                    "geometry_sha256": provenance_geometry_hash,
                    "regular_geometry_sha256": runs[0]["sha256"],
                    "source_records": 0,
                    "status_counts": {},
                }
                construction_summary = {
                    "schema": "cherchi-vertex-construction-source-diagnostic-v1",
                    "classification": "SYNTHETIC_REGRESSION",
                    "status": "DIAGNOSTIC_OUTPUT_MISMATCH",
                    "source_identity_status": "OUTPUT_MISMATCH",
                    "vertex_count": 0,
                    "point_type_counts": {},
                    "operand_source_counts": {},
                    "surface_plane_source_counts": {},
                    "all_operands_are_exact_input_vertices": False,
                    "all_surface_planes_have_unique_sources": False,
                }
            else:
                exact_rows = [line.split("\t") for line in
                              Path(str(provenance_path) + ".exact.tsv").read_text(
                                  encoding="ascii").splitlines()[1:]]
                fixed_box_vertex_audit = _fixed_box_vertex_audit(
                    exact_rows, design_vertices_on_disk)
                face_rows = [line.split("\t") for line in
                             Path(str(provenance_path) + ".faces.tsv").read_text(
                                 encoding="ascii").splitlines()[1:]]
                approximation_disagreements = sum(
                    row[1] != "MATCHED" for row in exact_rows
                )
                source_records = []
                source_path = None
                if not approximation_disagreements:
                    source_records = recover_face_sources(
                        (design_path, quadric_paths[name]), provenance_path,
                        exact_rows, face_rows,
                    )
                    source_path = output / f"cherchi-{name}-face-sources.json"
                    _write(source_path, {
                        "schema": "cherchi-face-source-diagnostic-v1",
                        "classification": "SYNTHETIC_REGRESSION",
                        "meaning": "candidate source input triangle per output face; not a derivative map",
                        "records": source_records,
                    })
                status_counts: dict[str, int] = {}
                for record in source_records:
                    status = str(record["status"])
                    status_counts[status] = status_counts.get(status, 0) + 1
                construction_sidecar_path = Path(str(provenance_path) + ".construction.tsv")
                construction_summary = _summarize_vertex_construction_provenance(
                    (design_path, quadric_paths[name]), provenance_path,
                    construction_sidecar_path,
                )
                if name in {"ellipsoid", "one_sheet_hyperboloid"} and not approximation_disagreements:
                    motion_summary = _cherchi_motion_finite_difference_audit(
                        fixture_name=name, extension=extension, probe=probe,
                        operation="intersection" if evaluation.quadric_solid_is_positive_region else "subtraction",
                        design_path=design_path, quadric_path=quadric_paths[name],
                        base_output_path=provenance_path, coefficients=coefficients,
                        request=request, base_evaluation=evaluation,
                        base_constructions=construction_summary, output_directory=output,
                    )
                    motion_path = output / f"cherchi-{name}-motion-differential.json"
                    _write(motion_path, motion_summary)
                    construction_summary["motion_differential"] = {
                        key: value for key, value in motion_summary.items() if key != "trials"
                    }
                    construction_summary["motion_differential_record_path"] = str(motion_path)
                    construction_summary["motion_differential_record_sha256"] = hashlib.sha256(
                        motion_path.read_bytes()).hexdigest()
                construction_source_path = output / f"cherchi-{name}-vertex-constructions.json"
                _write(construction_source_path, construction_summary)
                construction_summary = {
                    key: value for key, value in construction_summary.items()
                    if key != "records"
                }
                construction_summary["source_record_path"] = str(construction_source_path)
                construction_summary["source_record_sha256"] = hashlib.sha256(
                    construction_source_path.read_bytes()).hexdigest()
                provenance_summary = {
                    "status": ("APPROXIMATE_OUTPUT_COORDINATES" if approximation_disagreements
                               else "DIAGNOSTIC_ONLY"),
                    "geometry_sha256": provenance_geometry_hash,
                    "regular_geometry_sha256": runs[0]["sha256"],
                    "source_record_path": str(source_path) if source_path is not None else None,
                    "source_record_sha256": (hashlib.sha256(source_path.read_bytes()).hexdigest()
                                             if source_path is not None else None),
                    "source_records": len(source_records),
                    "status_counts": status_counts,
                    "approximation_disagreement_count": approximation_disagreements,
                    "fixed_box_vertex_audit": fixed_box_vertex_audit,
                    "provenance_run_stdout": provenance_run.stdout.strip(),
                }
            reference_vertices, reference_faces = _read_obj(Path(runs[0]["path"]))
            cgal_signature = _triangle_boundary_signature(
                evaluation.material_vertices, evaluation.material_faces)
            cherchi_signature = _triangle_boundary_signature(reference_vertices, reference_faces)
            # Vertex/face signatures are diagnostic only: two valid
            # triangulations can describe the same represented polyhedral
            # solid.  Use the exact represented-set relation for the actual
            # geometry conclusion.
            cgal_vertices_on_disk, cgal_faces_on_disk = _read_obj(cgal_material_paths[name])
            cgal_region = extension.make_exact_material_region(
                cgal_vertices_on_disk.tolist(), cgal_faces_on_disk.tolist())
            rational_relation, construction_region = _exact_construction_relation(
                extension, exact_rows, reference_vertices, reference_faces, cgal_region)
            cherchi_region = extension.make_exact_material_region(
                reference_vertices.tolist(), reference_faces.tolist())
            quadric_vertices_on_disk, quadric_faces_on_disk = _read_obj(quadric_paths[name])
            quadric_region = extension.make_exact_material_region(
                quadric_vertices_on_disk.tolist(), quadric_faces_on_disk.tolist())
            input_boolean_operation = "intersection" if operation == "intersection" else "difference"
            cgal_exact = extension.libigl_exact_boolean_construction(
                design_vertices_on_disk, design_faces_on_disk,
                quadric_vertices_on_disk, quadric_faces_on_disk,
                input_boolean_operation)
            cgal_exact_path = output / f"cgal-{name}-exact-construction.json"
            _write(cgal_exact_path, {
                "schema": "pinned-libigl-exact-Boolean-construction-v1",
                "scope": "serialized-binary64-input-polyhedra-before-output-rounding",
                "operation": input_boolean_operation,
                "design_sha256": hashlib.sha256(design_path.read_bytes()).hexdigest(),
                "cutter_sha256": hashlib.sha256(quadric_paths[name].read_bytes()).hexdigest(),
                "succeeded": bool(cgal_exact["succeeded"]),
                "vertices": cgal_exact["vertices"],
                "faces": np.asarray(cgal_exact["faces"], dtype=int).tolist(),
                "birth_faces": np.asarray(cgal_exact["birth_faces"], dtype=int).tolist(),
            })
            if cgal_exact["succeeded"]:
                cgal_exact_region = extension.make_exact_material_region_from_rationals(
                    cgal_exact["vertices"], cgal_exact["faces"])
                if cgal_exact_region.facts["valid"]:
                    exact_relation = dict(extension.relate_exact_boolean_operation_to_material(
                        design_region, quadric_region, cgal_exact_region, input_boolean_operation))
                    cgal_exact_relation = {
                        "completed": exact_relation["completed"],
                        "equal": exact_relation["equal"],
                        "symmetric_difference_volume_exact": exact_relation["symmetric_difference_volume_exact"],
                        "error": exact_relation["error"],
                    }
                else:
                    cgal_exact_relation = {"completed": False, "equal": None,
                                           "symmetric_difference_volume_exact": None,
                                           "error": cgal_exact_region.facts["error"]}
            else:
                cgal_exact_relation = {"completed": False, "equal": None,
                                       "symmetric_difference_volume_exact": None,
                                       "error": "pinned libigl exact-output Boolean failed"}
            if cgal_exact["succeeded"] and cgal_exact_region.facts["valid"] and construction_region is not None:
                direct = extension.relate_exact_material_regions(
                    cgal_exact_region, construction_region)
                internal_construction_relation = {
                    "completed": direct["completed"],
                    "equal": direct["equal"],
                    "symmetric_difference_volume_exact": direct["symmetric_difference_volume_exact"],
                    "error": direct["error"],
                }
            else:
                internal_construction_relation = {"completed": False, "equal": None,
                    "symmetric_difference_volume_exact": None,
                    "error": "one or both exact Boolean constructions are unavailable"}
            if not cgal_exact["succeeded"] or cgal_exact_relation.get("equal") is not True:
                cgal_rounding_moments = {"status": "UNKNOWN", "reason":
                    "exact construction was not validated against the input Boolean"}
            else:
                input_length = max(
                    max(float(row[axis]) for row in design_vertices_on_disk) -
                    min(float(row[axis]) for row in design_vertices_on_disk)
                    for axis in range(3))
                cgal_rounding_moments = compare_exact_rational_and_exported_moments(
                    cgal_exact["vertices"], cgal_exact["faces"].tolist(),
                    cgal_vertices_on_disk.tolist(), cgal_faces_on_disk.tolist(), input_length)
            input_oracle = {
                "scope": "exact-Boolean-of-serialized-input-polyhedra-only",
                "cgal_construction": cgal_exact_relation,
                "cgal_cherchi_exact_construction_relation": internal_construction_relation,
                "cgal_output_rounding_moments": cgal_rounding_moments,
                "cgal_export": dict(extension.relate_exact_boolean_operation_to_material(
                    design_region, quadric_region, cgal_region, input_boolean_operation)),
                "cherchi_export": dict(extension.relate_exact_boolean_operation_to_material(
                    design_region, quadric_region, cherchi_region, input_boolean_operation)),
                "cherchi_construction": (
                    dict(extension.relate_exact_boolean_operation_to_material(
                        design_region, quadric_region, construction_region, input_boolean_operation))
                    if construction_region is not None else
                    {"completed": False, "equal": None,
                     "error": rational_relation.get("reason", "exact construction unavailable")}),
            }
            operand_checks = {
                "cherchi": _operand_necessary_conditions(
                    extension, cherchi_region, design_region, quadric_region, operation),
                "cgal": _operand_necessary_conditions(
                    extension, cgal_region, design_region, quadric_region, operation),
            }
            exact_relation = extension.relate_exact_material_regions(cgal_region, cherchi_region)
            geometry_status = _exact_solid_relation_status(exact_relation)
            volume_comparison = _symmetric_difference_volume_comparison(
                exact_relation, evaluation.material_vertices, reference_vertices)
            scientific_fixtures.append({
                "name": name, "coefficients": list(coefficients), "operation": operation,
                "quadric_type": evaluation.interface_type, "quadric_subtype": evaluation.interface_subtype,
                "derivative_defined": evaluation.derivative_defined,
                "quadric_vertices": len(evaluation.interface_vertices), "quadric_faces": len(evaluation.interface_faces),
                "cgal_material_vertices": len(evaluation.material_vertices), "cgal_material_faces": len(evaluation.material_faces),
                "cgal_raw_moments": list(evaluation.variable_raw_moments), "cherchi_runs": runs,
                "face_source_provenance": provenance_summary,
                "vertex_construction_provenance": construction_summary,
                "geometry_comparison": {
                    "status": geometry_status,
                    "cgal_mesh": str(cgal_material_paths[name]),
                    "cgal_mesh_sha256": hashlib.sha256(cgal_material_paths[name].read_bytes()).hexdigest(),
                    "cgal_signature_sha256": cgal_signature,
                    "cherchi_signature_sha256": cherchi_signature,
                    "same_oriented_triangle_set": cgal_signature == cherchi_signature,
                    "operand_necessary_conditions": operand_checks,
                    "exact_construction_relation": rational_relation,
                    "exact_input_boolean_oracle": input_oracle,
                    "exact_relation": dict(exact_relation),
                    "symmetric_difference_volume_comparison": volume_comparison,
                },
            })
        geometry_statuses = {row["name"]: row["geometry_comparison"]["status"]
                             for row in scientific_fixtures}
        if all(value == "SAME_EXACT_POLYHEDRAL_SOLID" for value in geometry_statuses.values()):
            overall_geometry_status = "EXACT_SOLID_MATCH"
            overall_status = "SYNTHETIC_REFERENCE_BOOLEAN_CROSSCHECK_VALIDATED"
        elif any(value == "DIFFERENT_EXACT_POLYHEDRAL_SOLIDS" for value in geometry_statuses.values()):
            overall_geometry_status = "EXACT_SOLID_DISCREPANCY"
            overall_status = "SYNTHETIC_REFERENCE_BOOLEAN_GEOMETRY_DISCREPANCY"
        else:
            overall_geometry_status = "UNKNOWN"
            overall_status = "SYNTHETIC_REFERENCE_BOOLEAN_GEOMETRY_UNKNOWN"
        scientific = {
            "schema": _SCHEMA,
            "classification": "SYNTHETIC_REGRESSION",
            "scope": "SAME_OPERAND_CLIPPED_QUADRIC_BOOLEAN_ONLY_NOT_AUTHOR_ASSET_OR_PAPER_CASE",
            "design_vertices": len(request.design_vertices), "design_faces": len(request.design_faces),
            "geometry_status": overall_geometry_status,
            "fixtures": scientific_fixtures,
        }
        _write(output / "scientific.json", scientific)
        events.append({"iteration": 0, "status": "ok", "data": {
            "kind": "FASTER_SAME_OPERAND_REFERENCE_BOOLEAN_TYPE_MATRIX",
            "fixture_hashes": {row["name"]: row["cherchi_runs"][0]["sha256"] for row in scientific_fixtures},
            "geometry_statuses": geometry_statuses,
            "overall_geometry_status": overall_geometry_status,
            "maximum_absolute_error": max(row["cherchi_runs"][0]["maximum_absolute_error"] for row in scientific_fixtures),
        }})
        evidence = finish_run(begin, events=events,
                              output_files=sorted(path for path in output.rglob("*")
                                                  if path.is_file() and path.name != "begin.json"),
                              measurements={"started_at": started, "elapsed_seconds": time.perf_counter() - tick})
        _write(output / "evidence.json", evidence.as_dict())
        freshness = check_freshness(evidence)
        summary = {
            "status": overall_status,
            "geometry_status": overall_geometry_status,
            "geometry_discrepancy_fixtures": sorted(
                name for name, value in geometry_statuses.items()
                if value == "DIFFERENT_EXACT_POLYHEDRAL_SOLIDS"),
            "classification": "SYNTHETIC_REGRESSION", "project_status": "PARTIAL_BLOCKED",
            "original_evidence_admission": "BLOCKED_MISSING_AUTHOR_ASSETS_AND_PARAMETERS",
            "paper_case_admission": "BLOCKED_MISSING_AUTHOR_ASSETS_AND_PARAMETERS",
            "freshness": {"status": freshness.status, "mismatches": list(freshness.mismatches)},
        }
        _write(output / "summary.json", summary)
        return summary
    except Exception as error:
        if begin is not None and not (output / "evidence.json").exists():
            events.append({"iteration": len(events), "status": "failed", "data": {"error": str(error)}})
            _write(output / "evidence.json", finish_run(
                begin, events=events,
                output_files=sorted(path for path in output.rglob("*")
                                    if path.is_file() and path.name != "begin.json"),
                measurements={"started_at": started, "elapsed_seconds": time.perf_counter() - tick},
            ).as_dict())
        _write(output / "failure.json", {"status": "FAILED", "error": str(error),
                                           "traceback": traceback.format_exc()})
        raise


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-directory", type=Path, required=True)
    parser.add_argument("--build-directory", type=Path, required=True)
    parser.add_argument("--cherchi-probe", type=Path, required=True)
    parser.add_argument("--dependency-root", action="append", required=True, metavar="NAME=PATH")
    args = parser.parse_args()
    roots: dict[str, Path] = {}
    for item in args.dependency_root:
        name, separator, value = item.partition("=")
        if not name or not separator or name in roots:
            parser.error("invalid or duplicate dependency root")
        roots[name] = Path(value).resolve()
    result = run(args.output_directory, build_directory=args.build_directory,
                 cherchi_probe=args.cherchi_probe, dependency_roots=roots)
    print(canonical(result))
    return 0 if result["freshness"]["status"] == "CURRENT" else 1


if __name__ == "__main__":
    raise SystemExit(main())
