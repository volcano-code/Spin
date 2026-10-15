#!/usr/bin/env python3
"""Fail-closed audit for the fixed Spin-It computational requirements.

Validates catalogue completeness and source evidence. Registered exact
checkers may verify only their explicitly bounded requirements. Asset
readiness and caller claims never grant verification or publication.
"""

from __future__ import annotations

import argparse
from contextlib import contextmanager
import hashlib
import importlib
import json
import math
import re
import sys
import numpy as np
from fractions import Fraction
from pathlib import Path
from typing import Any, NoReturn

from requirements_catalogue import (
    CATALOGUE_SCHEMA,
    CATALOGUE_VERSION,
    REQUIRED_REQUIREMENT_IDS,
    catalogue,
)


EVIDENCE_SCHEMA = "spin-it-strict-evidence-v1"
EVIDENCE_VERSION = 1
AUDIT_SCHEMA = "spin-it-strict-audit-v1"
AUDIT_VERSION = 1

_CATALOGUE_FIELDS = {
    "requirement_id",
    "paper_locator",
    "requirement_kind",
    "expected_outcome",
    "dependencies",
    "required_evidence",
}
_EVIDENCE_FIELDS = {
    "requirement_id",
    "status",
    "evidence_description",
    "details",
}
_ASSET_FIELDS = {
    "asset_id",
    "asset_class",
    "path",
    "sha256",
    "source_reference",
    "required_critical_parameters",
    "parameter_evidence",
}
_PARAMETER_FIELDS = {
    "parameter",
    "status",
    "evidence_description",
    "value",
}
_BLOCKED_REQUIREMENT_STATUSES = {"NOT_VERIFIED", "NOT_IMPLEMENTED", "BLOCKED"}
_ALLOWED_ASSET_CLASSES = {"AUTHOR_ORIGINAL", "PUBLISHED_EXACT_SPECIFICATION"}
_ALLOWED_PARAMETER_STATUSES = {"PUBLISHED", "AUTHOR_PROVIDED", "EXACTLY_DERIVED"}
CASE_ASSET_READINESS_SCHEMA = "spin-it-case-asset-readiness-v1"
CASE_ASSET_READINESS_VERSION = 1
_CASE_ASSET_READINESS_FIELDS = {
    "requirement_id", "asset_status", "missing_artifacts",
    "missing_critical_parameters", "source_summary",
}
_CASE_ASSET_READINESS_DOCUMENT_FIELDS = {
    "schema", "schema_version", "source_audit_sha256", "cases",
}
_READY_CASE_ASSET_STATUSES = {"AUTHOR_ORIGINAL", "PUBLISHED_EXACT_SPECIFICATION"}
_CASE_ASSET_STATUSES = {"MISSING"} | _READY_CASE_ASSET_STATUSES


def _raise_nonfinite(constant: str) -> NoReturn:
    raise ValueError(f"nonfinite JSON value is not allowed: {constant}")


def _reject_nonfinite(value: Any, location: str) -> None:
    if isinstance(value, float) and not math.isfinite(value):
        raise ValueError(f"nonfinite value at {location}")
    if isinstance(value, dict):
        for key, child in value.items():
            _reject_nonfinite(child, f"{location}.{key}")
    elif isinstance(value, list):
        for index, child in enumerate(value):
            _reject_nonfinite(child, f"{location}[{index}]")


def _validate_parameter_value(value: Any, location: str) -> None:
    if value is None:
        raise ValueError(f"invalid parameter value at {location}: null is not evidence")
    if isinstance(value, bool):
        raise ValueError(f"invalid parameter value at {location}: booleans are not evidence")
    if isinstance(value, (int, float)):
        if isinstance(value, float) and not math.isfinite(value):
            raise ValueError(f"nonfinite parameter value at {location}")
        return
    if isinstance(value, str):
        if not value.strip():
            raise ValueError(f"invalid parameter value at {location}: text is empty")
        return
    if isinstance(value, list):
        if not value:
            raise ValueError(f"invalid parameter value at {location}: array is empty")
        for index, child in enumerate(value):
            _validate_parameter_value(child, f"{location}[{index}]")
        return
    if isinstance(value, dict):
        if not value:
            raise ValueError(f"invalid parameter value at {location}: object is empty")
        for key, child in value.items():
            if not isinstance(key, str) or not key.strip():
                raise ValueError(
                    f"invalid parameter value at {location}: object key is empty"
                )
            _validate_parameter_value(child, f"{location}.{key}")
        return
    raise ValueError(
        f"invalid parameter value at {location}: unsupported type {type(value).__name__}"
    )


def _mapping(value: Any, location: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ValueError(f"malformed {location}: expected an object")
    return value


def _list(value: Any, location: str) -> list[Any]:
    if not isinstance(value, list):
        raise ValueError(f"malformed {location}: expected an array")
    return value


def _nonempty_text(value: Any, location: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"empty evidence or malformed text at {location}")
    return value


def _exact_fields(row: dict[str, Any], expected: set[str], location: str) -> None:
    missing = sorted(expected - row.keys())
    unknown = sorted(row.keys() - expected)
    if missing or unknown:
        parts = []
        if missing:
            parts.append(f"missing fields {missing}")
        if unknown:
            parts.append(f"unknown fields {unknown}")
        raise ValueError(f"malformed {location}: {'; '.join(parts)}")


def _planar_number(value: Any, location: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"malformed Faster planar component {location}: expected finite number")
    return float(value)


def _planar_vector(value: Any, size: int, location: str) -> list[float]:
    if not isinstance(value, list) or len(value) != size:
        raise ValueError(f"malformed Faster planar component {location}: expected {size} values")
    return [_planar_number(item, location) for item in value]


def _planar_true(row: dict[str, Any], name: str) -> None:
    if row.get(name) is not True:
        raise ValueError(f"Faster planar component witness is false: {name}")


def _validate_faster_planar_component_result(name: str, raw: Any) -> dict[str, Any]:
    """Validate one numerical method witness emitted by a protected native test."""
    row = _mapping(raw, f"Faster planar component {name}")
    common = {"schema", "classification", "component"}
    if row.get("schema") != "faster-planar-component-result-v1" or \
            row.get("classification") != "SYNTHETIC_REGRESSION" or row.get("component") != name:
        raise ValueError("Faster planar component identity or synthetic scope is invalid")

    if name == "dogleg":
        fields = common | {
            "status", "solved_center", "known_root_center", "center_error",
            "normalized_residual", "residual_tolerance", "iterations",
            "function_evaluations", "function_evaluation_budget", "exact_root_residual",
            "zero_radius_initial_converged", "budget_exhaustion_distinguished",
            "topology_transition_detected",
        }
        _exact_fields(row, fields, "Faster planar dogleg result")
        if row["status"] != "CONVERGED":
            raise ValueError("Faster planar dogleg status is not converged")
        solved = _planar_vector(row["solved_center"], 2, "solved_center")
        target = _planar_vector(row["known_root_center"], 2, "known_root_center")
        if target != [1.0, 0.0]:
            raise ValueError("Faster planar dogleg known root differs from its registered fixture")
        reported_error = _planar_number(row["center_error"], "center_error")
        measured_error = math.dist(solved, target)
        residual = _planar_number(row["normalized_residual"], "normalized_residual")
        tolerance = _planar_number(row["residual_tolerance"], "residual_tolerance")
        exact_residual = _planar_number(row["exact_root_residual"], "exact_root_residual")
        if reported_error < 0 or abs(reported_error - measured_error) > 1e-14 or reported_error > 2e-7:
            raise ValueError("Faster planar dogleg center differs from the balanced root")
        if tolerance != 1e-10 or residual < 0 or residual > tolerance:
            raise ValueError("Faster planar dogleg residual exceeds its registered tolerance")
        if exact_residual < 0 or exact_residual > 2e-16:
            raise ValueError("Faster planar dogleg exact-root residual is too large")
        iterations, evaluations, budget = (
            row["iterations"], row["function_evaluations"], row["function_evaluation_budget"]
        )
        if any(isinstance(value, bool) or not isinstance(value, int) for value in
               (iterations, evaluations, budget)) or iterations < 0 or \
                not 1 <= evaluations <= budget == 60:
            raise ValueError("Faster planar dogleg function-evaluation budget is invalid")
        for key in ("zero_radius_initial_converged", "budget_exhaustion_distinguished",
                    "topology_transition_detected"):
            _planar_true(row, key)
    elif name == "polygon_disk":
        fields = common | {
            "analytic_moment_cases", "maximum_scaled_absolute_moment_error",
            "scaled_absolute_error_tolerance", "covers_holes",
            "covers_empty_and_tangent_intersections", "tangent_topology_transition_detected",
        }
        _exact_fields(row, fields, "Faster polygon/disk result")
        cases = row["analytic_moment_cases"]
        error = _planar_number(row["maximum_scaled_absolute_moment_error"], "moment_error")
        tolerance = _planar_number(row["scaled_absolute_error_tolerance"], "moment_tolerance")
        # Includes the translated open-arc/chord analytic regression. Historical
        # eight-case evidence remains historical; it cannot cover this fix.
        if isinstance(cases, bool) or not isinstance(cases, int) or cases != 9 or \
                tolerance != 2e-12 or error < 0 or error > tolerance:
            raise ValueError("Faster polygon/disk analytic moment contract failed")
        for key in ("covers_holes", "covers_empty_and_tangent_intersections",
                    "tangent_topology_transition_detected"):
            _planar_true(row, key)
    elif name == "extrusion":
        fields = common | {
            "planar_boundary_loop_count", "surface_shell_count", "extrusion_height", "density",
            "planar_mass", "extruded_mass", "maximum_scaled_absolute_moment_error",
            "exact_arc_extrusion_mesh_convergence_error",
            "maximum_scaled_objective_constraint_error", "planar_objective", "spatial_objective",
            "planar_constraints", "spatial_constraints", "scaled_absolute_error_tolerance",
        }
        _exact_fields(row, fields, "Faster extrusion result")
        for key in ("planar_boundary_loop_count", "surface_shell_count"):
            value = row[key]
            if isinstance(value, bool) or not isinstance(value, int) or value != 2:
                raise ValueError("Faster extrusion boundary topology contract failed")
        height = _planar_number(row["extrusion_height"], "extrusion_height")
        density = _planar_number(row["density"], "density")
        planar_mass = _planar_number(row["planar_mass"], "planar_mass")
        extruded_mass = _planar_number(row["extruded_mass"], "extruded_mass")
        moment_error = _planar_number(row["maximum_scaled_absolute_moment_error"], "moment_error")
        arc_mesh_error = _planar_number(
            row["exact_arc_extrusion_mesh_convergence_error"], "arc_mesh_convergence_error"
        )
        objective_error = _planar_number(
            row["maximum_scaled_objective_constraint_error"], "objective_constraint_error"
        )
        tolerance = _planar_number(row["scaled_absolute_error_tolerance"], "tolerance")
        if height != 1.5 or density != 2.4 or planar_mass <= 0 or extruded_mass <= 0 or \
                tolerance != 2e-12 or min(moment_error, objective_error) < 0 or \
                max(moment_error, objective_error) > tolerance:
            raise ValueError("Faster planar-to-3D extrusion moment contract failed")
        if arc_mesh_error < 0 or arc_mesh_error > 1e-6:
            raise ValueError("Faster exact arc-extrusion mesh convergence exceeds its independent 1e-6 gate")
        mass_error = abs(planar_mass - extruded_mass) / max(1.0, abs(planar_mass))
        if moment_error + 1e-14 < mass_error:
            raise ValueError("Faster extrusion moment error understates the mass discrepancy")
        planar_objective = _planar_number(row["planar_objective"], "planar_objective")
        spatial_objective = _planar_number(row["spatial_objective"], "spatial_objective")
        planar_constraints = _planar_vector(row["planar_constraints"], 2, "planar_constraints")
        spatial_constraints = _planar_vector(row["spatial_constraints"], 2, "spatial_constraints")
        objective_differences = [
            abs(planar_objective - spatial_objective) / max(1.0, abs(planar_objective)),
            *(abs(a - b) / max(1.0, abs(a)) for a, b in
              zip(planar_constraints, spatial_constraints)),
        ]
        if objective_error + 1e-14 < max(objective_differences):
            raise ValueError("Faster extrusion objective/constraint error is understated")
    else:
        raise ValueError("unregistered Faster planar component")
    _reject_nonfinite(row, f"Faster planar component {name}")
    return row


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _oriented_triangle_boundary_sha256(path: Path) -> str:
    """Independently hash an OBJ's oriented triangle boundary, ignoring order."""
    vertices: list[tuple[str, str, str]] = []
    triangles: list[tuple[tuple[str, str, str], ...]] = []
    with path.open("r", encoding="ascii") as stream:
        for line in stream:
            fields = line.split()
            if not fields:
                continue
            if fields[0] == "v" and len(fields) == 4:
                values = []
                for token in fields[1:]:
                    value = float(token)
                    if not math.isfinite(value):
                        raise ValueError("nonfinite Boolean boundary vertex")
                    values.append(format(0.0 if value == 0.0 else value, ".17g"))
                vertices.append(tuple(values))
            elif fields[0] == "f" and len(fields) == 4:
                face = tuple(int(token.split("/")[0]) - 1 for token in fields[1:])
                if len(set(face)) != 3 or any(index < 0 or index >= len(vertices) for index in face):
                    raise ValueError("invalid Boolean boundary triangle index")
                triangle = [vertices[index] for index in face]
                pivot = min(range(3), key=triangle.__getitem__)
                triangles.append(tuple(triangle[pivot:] + triangle[:pivot]))
            else:
                raise ValueError("Boolean boundary OBJ is not triangular")
    if not vertices or not triangles:
        raise ValueError("Boolean boundary OBJ is empty")
    payload = json.dumps(sorted(triangles), separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def _obj_mesh_arrays(path: Path) -> tuple[list[list[float]], list[list[int]]]:
    vertices: list[list[float]] = []
    faces: list[list[int]] = []
    with path.open("r", encoding="ascii") as stream:
        for line in stream:
            fields = line.split()
            if not fields:
                continue
            if fields[0] == "v" and len(fields) == 4:
                row = [float(value) for value in fields[1:]]
                if not all(math.isfinite(value) for value in row):
                    raise ValueError("nonfinite Boolean mesh vertex")
                vertices.append(row)
            elif fields[0] == "f" and len(fields) == 4:
                face = [int(value.split("/")[0]) - 1 for value in fields[1:]]
                if len(set(face)) != 3 or any(index < 0 or index >= len(vertices) for index in face):
                    raise ValueError("invalid Boolean mesh triangle index")
                faces.append(face)
            else:
                raise ValueError("Boolean material OBJ is not triangular")
    if not vertices or not faces:
        raise ValueError("Boolean material OBJ is empty")
    return vertices, faces


def _reference_boolean_expected_operation(name: str) -> str:
    """Independently bind the frozen fixture's solid side to its Boolean operation."""
    if name in {"ellipsoid", "one_sheet_hyperboloid", "two_sheet_hyperboloid"}:
        return "subtraction"
    if name in {"plane", "paraboloid"}:
        return "intersection"
    raise ValueError("unknown reference Boolean fixture: " + name)


def _reference_boolean_operand_conditions(native, operation: str, design, quadric, output) -> dict[str, Any]:
    """Recompute necessary exact-set conditions independently of the runner.

    These predicates can reject an invalid result; even both passing cannot
    prove equality with the complete input Boolean operation.
    """
    if operation not in {"intersection", "subtraction"}:
        raise ValueError("unsupported reference Boolean operation")
    digests = {"output_digest": output.facts["region_digest"],
               "design_digest": design.facts["region_digest"],
               "quadric_digest": quadric.facts["region_digest"]}
    if not all(region.facts["valid"] for region in (design, quadric, output)):
        return {**digests, "operation": operation, "completed": False,
                "subset_of_design": None, "quadric_operation_condition": None,
                "status": "UNKNOWN", "reason": "INVALID_EXACT_REGION"}
    try:
        design_fact = dict(native.check_exact_material_subset(output, design))
        if operation == "intersection":
            quadric_fact = dict(native.check_exact_material_subset(output, quadric))
            quadric_ok = quadric_fact["a_subset_b"] if quadric_fact["completed"] else None
            semantics = "output_subset_of_quadric_operand"
        else:
            quadric_fact = dict(native.relate_exact_material_regions(output, quadric))
            overlap = quadric_fact["positive_volume_overlap"] if quadric_fact["completed"] else None
            quadric_ok = None if overlap is None else not overlap
            semantics = "no_positive_volume_overlap_with_quadric_operand"
        design_ok = design_fact["a_subset_b"] if design_fact["completed"] else None
        completed = design_ok is not None and quadric_ok is not None
        status = ("FAIL" if design_ok is False or quadric_ok is False else
                  "PASS" if completed else "UNKNOWN")
        return {**digests, "operation": operation, "completed": completed,
                "subset_of_design": design_ok, "quadric_operation_condition": quadric_ok,
                "quadric_condition_semantics": semantics, "status": status,
                "design_subset_facts": design_fact, "quadric_operand_facts": quadric_fact}
    except Exception as error:
        return {**digests, "operation": operation, "completed": False,
                "subset_of_design": None, "quadric_operation_condition": None,
                "status": "UNKNOWN", "reason": str(error)}


def _reference_boolean_exact_construction_relation(
    native, sidecar_path: Path, exported_vertices: list[list[float]],
    exported_faces: list[list[int]], cgal_region,
) -> dict[str, Any]:
    """Independently compare Cherchi's rational construction to CGAL output.

    This is a represented-polyhedron diagnostic, never a continuous-quadric
    correctness certificate or a replacement for the ordinary OBJ check.
    """
    result: dict[str, Any] = {
        "scope": "cherchi-exact-construction-vs-cgal-exported-polyhedron-only",
        "status": "UNKNOWN", "construction_region_valid": False,
        "completed": False, "equal": None,
    }
    try:
        lines = sidecar_path.read_text(encoding="ascii").splitlines()
        if not lines or lines[0].split("\t") != [
            "vertex_index", "status", "exact_x", "exact_y", "exact_z",
            "memory_x", "memory_y", "memory_z", "obj_x", "obj_y", "obj_z",
        ] or len(lines) - 1 != len(exported_vertices):
            raise ValueError("construction sidecar is absent or has a different vertex count")

        def coordinate(token: str) -> Fraction:
            numerator, denominator = token.split("/")
            def integer(value: str) -> int:
                sign = -1 if value.startswith("-") else 1
                unsigned = value[1:] if sign < 0 else value
                return sign * int(unsigned[2:], 2) if unsigned.startswith("0b") else sign * int(unsigned, 10)
            return Fraction(integer(numerator), integer(denominator))

        vertices = []
        for index, line in enumerate(lines[1:]):
            row = line.split("\t")
            if len(row) != 11 or int(row[0]) != index:
                raise ValueError("construction sidecar is incomplete or reordered")
            if any(float(row[8 + axis]) != float(exported_vertices[index][axis])
                   for axis in range(3)):
                raise ValueError("construction sidecar does not bind the exported OBJ vertex")
            exact = [coordinate(row[2 + axis]) for axis in range(3)]
            vertices.append([f"{value.numerator}/{value.denominator}" for value in exact])
        region = native.make_exact_material_region_from_rationals(vertices, exported_faces)
        result["construction_region_digest"] = region.facts["region_digest"]
        result["construction_region_valid"] = bool(region.facts["valid"])
        if not result["construction_region_valid"]:
            result["reason"] = region.facts["error"]
            return result
        relation = native.relate_exact_material_regions(region, cgal_region)
        result["completed"] = bool(relation["completed"])
        if not result["completed"]:
            result["reason"] = relation["error"]
            return result
        result["equal"] = bool(relation["equal"])
        result["a_subset_b"] = bool(relation["a_subset_b"])
        result["b_subset_a"] = bool(relation["b_subset_a"])
        result["status"] = ("SAME_EXACT_POLYHEDRAL_SOLID" if result["equal"]
                            else "DIFFERENT_EXACT_POLYHEDRAL_SOLIDS")
    except (ValueError, RuntimeError, KeyError, IndexError) as error:
        result["reason"] = str(error)
    return result


def _reference_boolean_input_oracle(native, root: Path, name: str,
                                    operation: str) -> dict[str, Any]:
    """Recompute the input Boolean from bound files, independent of runner facts."""
    python_root = str(Path(__file__).resolve().parents[3] / "python")
    sys.path.insert(0, python_root)
    try:
        from make_it_stand.verification.mass_moments import (
            compare_exact_rational_and_exported_moments,
        )
    finally:
        sys.path.pop(0)
    if operation not in {"intersection", "subtraction"}:
        raise ValueError("unsupported reference Boolean operation")
    def region(path: Path):
        vertices, faces = _obj_mesh_arrays(path)
        return native.make_exact_material_region(vertices, faces)

    design = region(root / "design.obj")
    cutter = region(root / f"quadric-{name}.obj")
    cgal = region(root / f"cgal-material-{name}.obj")
    cherchi = region(root / f"cherchi-{name}-run-1.obj")
    exported_vertices, exported_faces = _obj_mesh_arrays(
        root / f"cherchi-{name}-run-1.obj")
    sidecar = root / f"cherchi-{name}-provenance.obj.exact.tsv"
    lines = sidecar.read_text(encoding="ascii").splitlines()
    expected_header = ["vertex_index", "status", "exact_x", "exact_y", "exact_z",
                       "memory_x", "memory_y", "memory_z", "obj_x", "obj_y", "obj_z"]
    if not lines or lines[0].split("\t") != expected_header or \
            len(lines) - 1 != len(exported_vertices):
        raise ValueError("exact construction sidecar does not bind Cherchi output")

    def coordinate(token: str) -> Fraction:
        numerator, denominator = token.split("/")
        def integer(value: str) -> int:
            sign = -1 if value.startswith("-") else 1
            unsigned = value[1:] if sign < 0 else value
            return sign * int(unsigned[2:], 2) if unsigned.startswith("0b") else sign * int(unsigned, 10)
        return Fraction(integer(numerator), integer(denominator))

    rational_vertices = []
    for index, line in enumerate(lines[1:]):
        row = line.split("\t")
        if len(row) != 11 or int(row[0]) != index or any(
            float(row[8 + axis]) != float(exported_vertices[index][axis])
            for axis in range(3)
        ):
            raise ValueError("exact construction sidecar vertex does not bind Cherchi OBJ")
        point = [coordinate(row[2 + axis]) for axis in range(3)]
        rational_vertices.append([f"{value.numerator}/{value.denominator}" for value in point])
    constructed = native.make_exact_material_region_from_rationals(
        rational_vertices, exported_faces)
    if not all(item.facts["valid"] for item in (design, cutter, cgal, cherchi, constructed)):
        raise ValueError("exact input Boolean has an invalid operand or output region")
    selected = "difference" if operation == "subtraction" else "intersection"
    design_path = root / "design.obj"
    cutter_path = root / f"quadric-{name}.obj"
    design_vertices, design_faces = _obj_mesh_arrays(design_path)
    cutter_vertices, cutter_faces = _obj_mesh_arrays(cutter_path)
    cgal_exact = native.libigl_exact_boolean_construction(
        np.asarray(design_vertices, dtype=float), np.asarray(design_faces, dtype=int),
        np.asarray(cutter_vertices, dtype=float), np.asarray(cutter_faces, dtype=int),
        selected)
    cgal_exact_path = root / f"cgal-{name}-exact-construction.json"
    regenerated = {
        "schema": "pinned-libigl-exact-Boolean-construction-v1",
        "scope": "serialized-binary64-input-polyhedra-before-output-rounding",
        "operation": selected,
        "design_sha256": hashlib.sha256(design_path.read_bytes()).hexdigest(),
        "cutter_sha256": hashlib.sha256(cutter_path.read_bytes()).hexdigest(),
        "succeeded": bool(cgal_exact["succeeded"]),
        "vertices": cgal_exact["vertices"],
        "faces": cgal_exact["faces"].tolist(),
        "birth_faces": cgal_exact["birth_faces"].tolist(),
    }
    saved_construction = json.loads(cgal_exact_path.read_text(encoding="ascii"))
    for field in ("schema", "scope", "operation", "design_sha256", "cutter_sha256", "succeeded"):
        value = regenerated[field]
        if saved_construction.get(field) != value:
            raise ValueError("saved CGAL exact construction differs from fresh libigl output: " + field)
    saved_faces = saved_construction.get("faces")
    saved_birth_faces = saved_construction.get("birth_faces")
    if not isinstance(saved_faces, list) or not isinstance(saved_birth_faces, list) or \
            len(saved_birth_faces) != len(saved_faces) or any(
                not isinstance(index, int) or index < 0 or
                index >= len(design_faces) + len(cutter_faces)
                for index in saved_birth_faces):
        raise ValueError("saved CGAL exact construction has invalid face provenance")
    if cgal_exact["succeeded"]:
        saved_region = native.make_exact_material_region_from_rationals(
            saved_construction["vertices"], saved_faces)
        cgal_exact_region = native.make_exact_material_region_from_rationals(
            cgal_exact["vertices"], cgal_exact["faces"])
        if saved_region.facts["valid"] and cgal_exact_region.facts["valid"]:
            replay_relation = native.relate_exact_material_regions(
                saved_region, cgal_exact_region)
            if not replay_relation["completed"] or not replay_relation["equal"]:
                raise ValueError("saved CGAL exact construction differs as a represented solid")
            exact_relation = dict(native.relate_exact_boolean_operation_to_material(
                design, cutter, cgal_exact_region, selected))
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
    if cgal_exact["succeeded"] and saved_region.facts["valid"] and \
            cgal_exact_region.facts["valid"]:
        direct = native.relate_exact_material_regions(cgal_exact_region, constructed)
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
        input_length = max(max(float(row[axis]) for row in design_vertices) -
                           min(float(row[axis]) for row in design_vertices)
                           for axis in range(3))
        rounded_vertices, rounded_faces = _obj_mesh_arrays(
            root / f"cgal-material-{name}.obj")
        cgal_rounding_moments = compare_exact_rational_and_exported_moments(
            cgal_exact["vertices"], cgal_exact["faces"].tolist(),
            rounded_vertices, rounded_faces, input_length)
    return {
        "scope": "exact-Boolean-of-serialized-input-polyhedra-only",
        "cgal_construction": cgal_exact_relation,
        "cgal_cherchi_exact_construction_relation": internal_construction_relation,
        "cgal_output_rounding_moments": cgal_rounding_moments,
        "cgal_export": dict(native.relate_exact_boolean_operation_to_material(
            design, cutter, cgal, selected)),
        "cherchi_export": dict(native.relate_exact_boolean_operation_to_material(
            design, cutter, cherchi, selected)),
        "cherchi_construction": dict(native.relate_exact_boolean_operation_to_material(
            design, cutter, constructed, selected)),
    }


def _exact_symmetric_difference_volume_comparison(
    facts: dict[str, Any], vertices_a: list[list[float]], vertices_b: list[list[float]],
) -> dict[str, Any]:
    """Independently derive a scale-normalized exact-volume diagnostic.

    A within-threshold result is not exact set equality and does not satisfy
    the strict Boolean geometry checker. The threshold is 2^-46 * D^3, where D
    is the combined bounds' largest axis span; it is a fixed diagnostic policy,
    not a Hausdorff or topology guarantee.
    """
    raw = facts.get("symmetric_difference_volume_exact")
    if facts.get("completed") is not True or not isinstance(raw, str) or not vertices_a or not vertices_b:
        return {"status": "UNKNOWN", "exact_set_equal": facts.get("equal"),
                "symmetric_difference_volume_exact": raw}
    exact = Fraction(raw)
    points = [[Fraction(float(coordinate)) for coordinate in point]
              for point in [*vertices_a, *vertices_b]]
    spans = [max(point[axis] for point in points) - min(point[axis] for point in points)
             for axis in range(3)]
    scale = max(spans)
    if scale <= 0:
        return {"status": "UNKNOWN", "exact_set_equal": facts.get("equal"),
                "symmetric_difference_volume_exact": raw}
    normalized = exact / (scale ** 3)
    tolerance = Fraction(1, 2**46)
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


def _load_bound_native_module_and_freshness(evidence: Any, module_path: Path, check_freshness: Any):
    """Check run freshness against its observed extension without cache leakage."""
    current = sys.modules.get("_make_it_stand")
    if current is not None and Path(getattr(current, "__file__", "")).resolve() == module_path:
        return current, check_freshness(evidence)
    sys.modules.pop("_make_it_stand", None)
    sys.path.insert(0, str(module_path.parent))
    try:
        loaded = importlib.import_module("_make_it_stand")
        freshness = check_freshness(evidence)
        if Path(getattr(loaded, "__file__", "")).resolve() != module_path:
            raise ValueError("native module loaded from " + str(getattr(loaded, "__file__", "")) +
                             " but evidence binds " + str(module_path))
        return loaded, freshness
    finally:
        sys.path.pop(0)
        sys.modules.pop("_make_it_stand", None)
        if current is not None:
            sys.modules["_make_it_stand"] = current


def _complete_ctest_receipt(inventory_text: str, execution_text: str) -> tuple[int, list[int]]:
    """Independently require one passed receipt for every listed CTest number."""
    entries = [int(item) for item in re.findall(r"^\s*Test\s+#(\d+): .+$", inventory_text, flags=re.MULTILINE)]
    total_match = re.search(r"^Total Tests: (\d+)\s*$", inventory_text, flags=re.MULTILINE)
    if total_match is None:
        raise ValueError("CTest inventory lacks Total Tests")
    total = int(total_match.group(1))
    if entries != list(range(1, total + 1)):
        raise ValueError("CTest inventory is not a complete numbered list")
    if "Not Run" in execution_text:
        raise ValueError("CTest execution contains Not Run")
    receipts = re.findall(
        r"^\s*(\d+)/(\d+)\s+Test\s+#(\d+):.*\s+(Passed|Failed)\b.*$",
        execution_text,
        flags=re.MULTILINE,
    )
    if len(receipts) != total:
        raise ValueError("CTest receipt count differs from inventory")
    test_ids = [int(test_id) for _, _, test_id, _ in receipts]
    if test_ids != entries or any(int(denominator) != total for _, denominator, _, _ in receipts):
        raise ValueError("CTest execution does not cover the complete inventory")
    if any(status != "Passed" for _, _, _, status in receipts):
        raise ValueError("CTest execution has a non-passing receipt")
    summary = re.search(r"^100% tests passed, 0 tests failed out of (\d+)\s*$", execution_text,
                        flags=re.MULTILINE)
    if summary is None or int(summary.group(1)) != total:
        raise ValueError("CTest execution summary disagrees with inventory")
    return total, test_ids


def load_catalogue() -> dict[str, Any]:
    """Return and validate a fresh copy of the versioned fixed catalogue."""

    document = catalogue()
    if document.get("schema") != CATALOGUE_SCHEMA:
        raise ValueError("malformed fixed catalogue schema")
    if document.get("schema_version") != CATALOGUE_VERSION:
        raise ValueError("malformed fixed catalogue version")
    validate_catalogue(document.get("requirements"))
    return document


def validate_catalogue(rows: Any) -> None:
    """Reject any incomplete, duplicated, unknown, or altered catalogue."""

    entries = _list(rows, "catalogue requirements")
    ids: list[str] = []
    for index, raw_row in enumerate(entries):
        row = _mapping(raw_row, f"catalogue requirement {index}")
        _exact_fields(row, _CATALOGUE_FIELDS, f"catalogue requirement {index}")
        requirement_id = _nonempty_text(
            row["requirement_id"], f"catalogue requirement {index}.requirement_id"
        )
        ids.append(requirement_id)

    duplicates = sorted({item for item in ids if ids.count(item) > 1})
    if duplicates:
        raise ValueError(f"duplicate catalogue requirement IDs: {duplicates}")

    required_ids = set(REQUIRED_REQUIREMENT_IDS)
    actual_ids = set(ids)
    unknown = sorted(actual_ids - required_ids)
    missing = sorted(required_ids - actual_ids)
    if unknown:
        raise ValueError(f"unknown catalogue requirement IDs: {unknown}; missing: {missing}")
    if missing:
        raise ValueError(f"missing catalogue requirement IDs: {missing}")

    canonical = {
        row["requirement_id"]: row for row in catalogue()["requirements"]
    }
    for row in entries:
        requirement_id = row["requirement_id"]
        if row != canonical[requirement_id]:
            raise ValueError(
                f"catalogue requirement {requirement_id} differs from the fixed catalogue"
            )


def _validated_parameter_evidence(
    asset_id: str,
    declared_parameters: list[str],
    raw_evidence: Any,
) -> list[dict[str, Any]]:
    evidence_rows = _list(raw_evidence, f"asset {asset_id} parameter_evidence")
    evidence: list[dict[str, Any]] = []
    names: list[str] = []
    for index, raw_row in enumerate(evidence_rows):
        row = _mapping(raw_row, f"asset {asset_id} parameter evidence {index}")
        _exact_fields(
            row,
            _PARAMETER_FIELDS,
            f"asset {asset_id} parameter evidence {index}",
        )
        name = _nonempty_text(
            row["parameter"], f"asset {asset_id} parameter evidence {index}.parameter"
        )
        names.append(name)
        status = _nonempty_text(
            row["status"], f"asset {asset_id} parameter {name}.status"
        )
        if status == "UNKNOWN_CRITICAL":
            raise ValueError(f"asset {asset_id} parameter {name} is UNKNOWN_CRITICAL")
        if status not in _ALLOWED_PARAMETER_STATUSES:
            raise ValueError(
                f"asset {asset_id} parameter {name} has unsupported evidence status {status}"
            )
        _nonempty_text(
            row["evidence_description"],
            f"asset {asset_id} parameter {name}.evidence_description",
        )
        _validate_parameter_value(
            row["value"], f"asset {asset_id} parameter {name}.value"
        )
        evidence.append(dict(row))

    duplicates = sorted({name for name in names if names.count(name) > 1})
    if duplicates:
        raise ValueError(f"duplicate parameter evidence for asset {asset_id}: {duplicates}")
    declared = set(declared_parameters)
    actual = set(names)
    unknown = sorted(actual - declared)
    if unknown:
        raise ValueError(f"unknown parameter evidence for asset {asset_id}: {unknown}")
    missing = sorted(declared - actual)
    if missing:
        raise ValueError(f"missing critical parameter evidence for asset {asset_id}: {missing}")
    return evidence


def preflight_assets(assets: Any, asset_root: Path) -> dict[str, Any]:
    """Validate exact asset provenance, bytes, and critical parameters.

    ``READY`` means only that inputs are ready to be consumed by future
    independent checkers.  It is not a verification result.
    """

    rows = _list(assets, "assets")
    root = Path(asset_root).resolve()
    if not rows:
        return {
            "schema": "spin-it-strict-asset-preflight-v1",
            "schema_version": 1,
            "status": "BLOCKED",
            "ready": False,
            "source_authenticity_verified": False,
            "strict_requirement_verification": False,
            "assets": [],
            "blockers": ["no exact assets were declared"],
        }

    validated: list[dict[str, Any]] = []
    asset_ids: list[str] = []
    for index, raw_asset in enumerate(rows):
        asset = _mapping(raw_asset, f"asset {index}")
        _exact_fields(asset, _ASSET_FIELDS, f"asset {index}")
        asset_id = _nonempty_text(asset["asset_id"], f"asset {index}.asset_id")
        asset_ids.append(asset_id)

        asset_class = _nonempty_text(asset["asset_class"], f"asset {asset_id}.asset_class")
        if asset_class == "RECONSTRUCTED_ASSET":
            raise ValueError(f"reconstructed asset {asset_id} is not admissible")
        if asset_class not in _ALLOWED_ASSET_CLASSES:
            raise ValueError(f"unsupported asset class {asset_class} for asset {asset_id}")

        source_reference = _nonempty_text(
            asset["source_reference"], f"asset {asset_id}.source_reference"
        )
        declared_raw = _list(
            asset["required_critical_parameters"],
            f"asset {asset_id}.required_critical_parameters",
        )
        if not declared_raw:
            raise ValueError(f"asset {asset_id} declares no required critical parameters")
        declared = [
            _nonempty_text(value, f"asset {asset_id}.required_critical_parameters[{index}]")
            for index, value in enumerate(declared_raw)
        ]
        duplicates = sorted({name for name in declared if declared.count(name) > 1})
        if duplicates:
            raise ValueError(f"duplicate critical parameters for asset {asset_id}: {duplicates}")
        _validated_parameter_evidence(asset_id, declared, asset["parameter_evidence"])

        relative_path = Path(_nonempty_text(asset["path"], f"asset {asset_id}.path"))
        if relative_path.is_absolute():
            raise ValueError(f"asset {asset_id} path must be relative to the asset root")
        resolved = (root / relative_path).resolve()
        try:
            resolved.relative_to(root)
        except ValueError as error:
            raise ValueError(f"asset {asset_id} path escapes the asset root") from error
        if not resolved.is_file():
            raise ValueError(f"asset {asset_id} file does not exist: {relative_path}")

        declared_digest = _nonempty_text(asset["sha256"], f"asset {asset_id}.sha256")
        if len(declared_digest) != 64:
            raise ValueError(f"malformed SHA-256 digest for asset {asset_id}")
        try:
            int(declared_digest, 16)
        except ValueError as error:
            raise ValueError(f"malformed SHA-256 digest for asset {asset_id}") from error
        actual_digest = _sha256_file(resolved)
        if actual_digest != declared_digest.lower():
            raise ValueError(
                f"digest mismatch for asset {asset_id}: expected {declared_digest.lower()}, "
                f"got {actual_digest}"
            )

        validated.append(
            {
                "asset_id": asset_id,
                "asset_class": asset_class,
                "path": relative_path.as_posix(),
                "sha256": actual_digest,
                "source_reference": source_reference,
                "critical_parameters": declared,
                "status": "READY",
            }
        )

    duplicate_assets = sorted(
        {asset_id for asset_id in asset_ids if asset_ids.count(asset_id) > 1}
    )
    if duplicate_assets:
        raise ValueError(f"duplicate asset IDs: {duplicate_assets}")

    return {
        "schema": "spin-it-strict-asset-preflight-v1",
        "schema_version": 1,
        "status": "READY",
        "ready": True,
        "source_authenticity_verified": False,
        "strict_requirement_verification": False,
        "assets": validated,
        "blockers": [],
        "readiness_scope": (
            "syntax, declared provenance, and checked file-byte integrity only; "
            "not source authenticity or strict requirement verification"
        ),
    }


def _required_original_experiment_ids() -> set[str]:
    return {
        row["requirement_id"] for row in catalogue()["requirements"]
        if row["requirement_kind"] == "ORIGINAL_EXPERIMENT"
    }


def paper_scope_requirement_ids(paper: str) -> tuple[str, ...]:
    """Return an immutable scope derived from the complete fixed catalogue."""
    if paper not in {"faster", "spin2014", "all"}:
        raise ValueError(f"invalid paper scope {paper!r}; expected faster, spin2014, or all")
    ids = set(REQUIRED_REQUIREMENT_IDS)
    if paper == "all":
        selected = ids
    elif paper == "faster":
        selected = {
            requirement_id for requirement_id in ids
            if requirement_id.startswith("faster.")
            or requirement_id.startswith("shared.")
            or requirement_id == "physical.faster.measurements"
        }
    else:
        selected = {
            requirement_id for requirement_id in ids
            if requirement_id.startswith("spin2014.")
            or requirement_id.startswith("shared.")
            or requirement_id == "physical.2014.measurements"
        }
    return tuple(sorted(selected))


def derive_paper_status(paper: str, requirement_rows: list[dict[str, Any]],
                        case_readiness: dict[str, dict[str, Any]]) -> dict[str, Any]:
    """Derive a paper-level completion result only from audited rows and case readiness."""
    scoped_ids = set(paper_scope_requirement_ids(paper))
    row_by_id = {row["requirement_id"]: row for row in requirement_rows
                 if row["requirement_id"] in scoped_ids}
    blockers: list[str] = []
    for requirement_id in sorted(scoped_ids - row_by_id.keys()):
        blockers.append(f"{requirement_id}: missing audited requirement row")
    for requirement_id, row in sorted(row_by_id.items()):
        if row["status"] not in {"VERIFIED", "NOT_IN_SCOPE"}:
            blockers.append(
                f"{requirement_id}: {row['status']} — {row['evidence_description']}"
            )
    catalogue_rows = {row["requirement_id"]: row
                      for row in catalogue()["requirements"]}
    experiment_ids = {
        requirement_id for requirement_id in scoped_ids
        if catalogue_rows[requirement_id]["requirement_kind"] == "ORIGINAL_EXPERIMENT"
    }
    for requirement_id in sorted(experiment_ids):
        readiness = case_readiness.get(requirement_id)
        if readiness is None or readiness.get("asset_status") not in _READY_CASE_ASSET_STATUSES:
            blockers.append(
                f"{requirement_id}: original experiment lacks exact published/author asset readiness"
            )
    scoped_rows = [row for requirement_id, row in row_by_id.items()]
    non_experiment_ids = scoped_ids - experiment_ids - {
        requirement_id for requirement_id in scoped_ids
        if catalogue_rows[requirement_id]["requirement_kind"] == "PHYSICAL_MEASUREMENT"
    }
    formula_coverage = bool(non_experiment_ids) and all(
        row_by_id.get(requirement_id, {}).get("status") == "VERIFIED"
        for requirement_id in non_experiment_ids
    )
    experiment_coverage = bool(experiment_ids) and all(
        row_by_id.get(requirement_id, {}).get("status") == "VERIFIED"
        and case_readiness.get(requirement_id, {}).get("asset_status") in
        _READY_CASE_ASSET_STATUSES
        for requirement_id in experiment_ids
    )
    complete = not blockers and len(row_by_id) == len(scoped_ids)
    return {
        "paper": paper,
        "requirement_count": len(scoped_ids),
        "verified_count": sum(row.get("status") == "VERIFIED" for row in scoped_rows),
        "original_experiment_count": len(experiment_ids),
        "formula_coverage": formula_coverage,
        "experiment_coverage": experiment_coverage,
        "status": "VERIFIED_COMPUTATIONAL_REPRODUCTION" if complete else "PARTIAL_BLOCKED",
        "public_complete_reproduction": complete,
        "blockers": blockers,
    }


def _nonempty_text_list(value: Any, location: str) -> list[str]:
    rows = _list(value, location)
    values = [_nonempty_text(item, f"{location}[{index}]")
              for index, item in enumerate(rows)]
    if len(set(values)) != len(values):
        raise ValueError(f"duplicate values in {location}")
    return values


def load_case_asset_readiness(path: Path) -> dict[str, dict[str, Any]]:
    """Compatibility row view of the validated, versioned source record."""
    return load_case_asset_readiness_record(path)["cases"]


def load_case_asset_readiness_record(path: Path) -> dict[str, Any]:
    """Load an exhaustive, source-audit-bound eligibility inventory.

    The inventory expresses missing input facts only.  It is not an asset
    authenticity proof and cannot by itself verify any paper requirement.
    """

    document = _load_json(Path(path))
    _mapping(document, "case asset readiness")
    version = document.get("schema_version")
    if type(version) is not int or version not in (1, 2):
        raise ValueError("unsupported case asset readiness version")
    schema = CASE_ASSET_READINESS_SCHEMA if version == 1 else "spin-it-case-asset-readiness-v2"
    if document.get("schema") != schema:
        raise ValueError("unsupported case asset readiness schema/version pairing")
    fields = _CASE_ASSET_READINESS_DOCUMENT_FIELDS | ({"source_audit_path"} if version == 2 else set())
    _exact_fields(document, fields, "case asset readiness")
    if version == 1:
        source_audit = Path(__file__).resolve().parents[1] / "sources/strict_asset_audit.md"
        scope = "LEGACY_MUTABLE_SOURCE_AUDIT_NOT_ASSET_AUTHENTICITY"
    else:
        digest = document["source_audit_sha256"]
        if (not isinstance(digest, str) or len(digest) != 64
                or any(character not in "0123456789abcdef" for character in digest)):
            raise ValueError("invalid source audit SHA-256")
        expected_path = f"snapshots/{digest}.md"
        if document["source_audit_path"] != expected_path:
            raise ValueError("invalid case asset readiness snapshot path")
        root = Path(path).resolve(strict=True).parent
        source_audit = root / expected_path
        if source_audit.resolve(strict=True) != source_audit:
            raise ValueError("case asset readiness snapshot path must not use symlinks")
        scope = "FROZEN_SOURCE_SNAPSHOT_NOT_ASSET_AUTHENTICITY"
    expected_digest = _sha256_file(source_audit)
    if document["source_audit_sha256"] != expected_digest:
        raise ValueError("case asset readiness source audit SHA-256 mismatch")

    cases = _list(document["cases"], "case asset readiness cases")
    rows: dict[str, dict[str, Any]] = {}
    for index, raw_row in enumerate(cases):
        row = _mapping(raw_row, f"case asset readiness case {index}")
        _exact_fields(row, _CASE_ASSET_READINESS_FIELDS,
                      f"case asset readiness case {index}")
        requirement_id = _nonempty_text(
            row["requirement_id"], f"case asset readiness case {index}.requirement_id"
        )
        if requirement_id in rows:
            raise ValueError(f"duplicate case asset readiness ID: {requirement_id}")
        status = _nonempty_text(
            row["asset_status"], f"case asset readiness {requirement_id}.asset_status"
        )
        if status not in _CASE_ASSET_STATUSES:
            raise ValueError(
                f"unsupported case asset readiness status {status} for {requirement_id}"
            )
        artifacts = _nonempty_text_list(
            row["missing_artifacts"],
            f"case asset readiness {requirement_id}.missing_artifacts",
        )
        parameters = _nonempty_text_list(
            row["missing_critical_parameters"],
            f"case asset readiness {requirement_id}.missing_critical_parameters",
        )
        _nonempty_text(row["source_summary"],
                       f"case asset readiness {requirement_id}.source_summary")
        if status in _READY_CASE_ASSET_STATUSES and (artifacts or parameters):
            raise ValueError(
                f"ready original case {requirement_id} cannot list missing evidence"
            )
        if status == "MISSING" and (not artifacts or not parameters):
            raise ValueError(
                f"missing case {requirement_id} requires artifact and parameter gaps"
            )
        rows[requirement_id] = dict(row)

    expected = _required_original_experiment_ids()
    actual = set(rows)
    missing = sorted(expected - actual)
    unknown = sorted(actual - expected)
    if missing or unknown:
        raise ValueError(
            f"case asset readiness IDs differ from original experiments: missing {missing}; unknown {unknown}"
        )
    return {"cases": rows, "source_audit": {
        "path": str(source_audit.resolve(strict=True)), "sha256": expected_digest,
        "scope": scope, "binding_version": version}}


def _validate_requirement_evidence(raw_requirements: Any) -> list[dict[str, Any]]:
    rows = _list(raw_requirements, "requirement evidence")
    ids: list[str] = []
    validated: list[dict[str, Any]] = []
    requirement_kinds = {
        row["requirement_id"]: row["requirement_kind"]
        for row in catalogue()["requirements"]
    }
    for index, raw_row in enumerate(rows):
        row = _mapping(raw_row, f"requirement evidence {index}")
        _exact_fields(row, _EVIDENCE_FIELDS, f"requirement evidence {index}")
        requirement_id = _nonempty_text(
            row["requirement_id"], f"requirement evidence {index}.requirement_id"
        )
        ids.append(requirement_id)
        if requirement_id not in requirement_kinds:
            continue
        status = _nonempty_text(
            row["status"], f"requirement {requirement_id}.status"
        )
        if status == "VERIFIED":
            raise ValueError(
                f"VERIFIED is not admissible for {requirement_id}: caller status cannot replace independent checks"
            )
        if requirement_kinds[requirement_id] == "PHYSICAL_MEASUREMENT":
            if status != "NOT_IN_SCOPE":
                raise ValueError(
                    f"physical requirement {requirement_id} must remain NOT_IN_SCOPE"
                )
        elif status not in _BLOCKED_REQUIREMENT_STATUSES:
            raise ValueError(f"unsupported blocked status {status} for {requirement_id}")
        _nonempty_text(
            row["evidence_description"],
            f"requirement {requirement_id}.evidence_description",
        )
        details = _mapping(row["details"], f"requirement {requirement_id}.details")
        _reject_nonfinite(details, f"requirement {requirement_id}.details")
        validated.append(
            {
                "requirement_id": requirement_id,
                "status": status,
                "evidence_description": row["evidence_description"],
                "details": dict(details),
            }
        )

    duplicates = sorted({item for item in ids if ids.count(item) > 1})
    if duplicates:
        raise ValueError(f"duplicate requirement evidence IDs: {duplicates}")
    required = set(REQUIRED_REQUIREMENT_IDS)
    actual = set(ids)
    unknown = sorted(actual - required)
    missing = sorted(required - actual)
    if unknown:
        raise ValueError(f"unknown requirement evidence IDs: {unknown}; missing: {missing}")
    if missing:
        raise ValueError(f"missing requirement evidence IDs: {missing}")
    return validated


def _independent_checks(requests: Any) -> dict[str, dict[str, Any]]:
    requests = _mapping(requests, "independent evidence")
    registered = {
        "faster.section6_6",
        "faster.discrete_moments",
        "faster.reference_boolean",
        "faster.stratum_matrix",
        "faster.eq14_newton",
        "faster.moment_foundations",
        "faster.hydrostatic_geometry",
        "faster.planar_methods",
        "faster.static_application_contract",
        "faster.buoyant_application_contract",
        "faster.yoyo_application_contract",
        "faster.spinning_top_application_contract",
        "faster.cad_methods",
        "spin2014.givens_equations",
        "spin2014.source_cage",
        "spin2014.source_deformation",
        "spin2014.dual_material",
        "spin2014.source_method",
        "spin2014.source_adaptive_transitions",
        "spin2014.supplement_algorithms",
        "shared.registered_tests",
        "shared.independent_moments",
    }
    if set(requests) - registered:
        raise ValueError("unregistered independent checker")
    if not requests:
        return {}

    def repeated_paths(key: str, label: str) -> list[Path]:
        paths = requests[key]
        if not isinstance(paths, list) or len(paths) != 2 or not all(isinstance(p, str) for p in paths):
            raise ValueError(label + " requires two independent evidence paths")
        paths = [Path(p).resolve(strict=True) for p in paths]
        if paths[0].samefile(paths[1]):
            raise ValueError("repeated evidence is not two runs")
        return paths

    verified: dict[str, dict[str, Any]] = {}
    if "shared.registered_tests" in requests:
        paths = repeated_paths("shared.registered_tests", "registered CTest evidence")
        python_root = str(Path(__file__).resolve().parents[3] / "python")
        sys.path.insert(0, python_root)
        try:
            from make_it_stand.verification.run_evidence import PaperRunEvidenceV2, check_freshness
        finally:
            sys.path.pop(0)
        signatures = []
        for path in paths:
            evidence = PaperRunEvidenceV2.from_json(path.read_text(encoding="utf-8"))
            if check_freshness(evidence).status != "CURRENT":
                raise ValueError("STALE_RESULT: registered CTest evidence is not current")
            payload = evidence.as_dict()
            parameters = payload["begin"]["request"]["parameters"]
            if (parameters.get("schema") != "registered-ctest-inputs-v2" or
                    parameters.get("scope") != "REGISTERED_CTEST_COMPLETE_EXECUTION_GATE_ONLY_NOT_PAPER_CASE_EVIDENCE" or
                    parameters.get("mode") != "COMPLETE_EXECUTION"):
                raise ValueError("registered CTest evidence requires COMPLETE_EXECUTION")
            outputs = {Path(row["path"]).resolve() for row in payload["scientific"]["outputs"]}
            root = path.parent
            expected = {root / "inputs.json", root / "scientific.json", root / "ctest-inventory.txt",
                        root / "ctest-execution.txt"}
            if not {item.resolve() for item in expected} <= outputs:
                raise ValueError("registered CTest evidence does not bind complete execution outputs")
            scientific = _load_json(root / "scientific.json")
            if (set(scientific) != {"schema", "scope", "mode", "inventory", "execution"} or
                    scientific["schema"] != "registered-ctest-scientific-v2" or
                    scientific["scope"] != parameters["scope"] or
                    scientific["mode"] != "COMPLETE_EXECUTION"):
                raise ValueError("registered CTest scientific record has invalid scope")
            inventory = root / "ctest-inventory.txt"
            transcript = root / "ctest-execution.txt"
            total, test_ids = _complete_ctest_receipt(
                inventory.read_text(encoding="utf-8"), transcript.read_text(encoding="utf-8"),
            )
            execution = scientific["execution"]
            if (not isinstance(execution, dict) or execution.get("accepted") is not True or
                    execution.get("returncode") != 0 or execution.get("path") != str(transcript) or
                    execution.get("sha256") != _sha256_file(transcript) or
                    execution.get("summary") != {"inventory_total": total, "executed_test_ids": test_ids}):
                raise ValueError("registered CTest scientific execution record is inconsistent")
            inventory_record = scientific["inventory"]
            if (not isinstance(inventory_record, dict) or inventory_record.get("path") != str(inventory) or
                    inventory_record.get("sha256") != _sha256_file(inventory) or inventory_record.get("total") != total):
                raise ValueError("registered CTest scientific inventory record is inconsistent")
            events = payload["scientific"]["events"]
            if (len(events) != 1 or events[0]["status"] != "ok" or
                    events[0]["data"].get("kind") != "REGISTERED_CTEST_EXECUTION" or
                    events[0]["data"].get("mode") != "COMPLETE_EXECUTION" or
                    events[0]["data"].get("inventory_total") != total):
                raise ValueError("registered CTest event receipt is incomplete")
            ctest = Path(parameters["ctest"]["path"]).resolve(strict=True)
            if parameters["ctest"].get("sha256") != _sha256_file(ctest):
                raise ValueError("STALE_RESULT: registered CTest executable differs from observation")
            signatures.append(hashlib.sha256(json.dumps({
                "ctest_sha256": parameters["ctest"]["sha256"], "inventory_sha256": _sha256_file(inventory),
                "inventory_total": total, "executed_test_ids": test_ids,
            }, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")).hexdigest())
        if signatures[0] != signatures[1]:
            raise ValueError("repeated registered CTest scientific outputs differ")
        verified["shared.registered_tests"] = {
            "checker": {"path": str(Path(__file__).resolve()), "sha256": _sha256_file(Path(__file__)),
                        "size": Path(__file__).stat().st_size},
            "evidence": [{"path": str(path), "sha256": _sha256_file(path), "size": path.stat().st_size}
                         for path in paths],
            "scientific_data_sha256": signatures[0],
            "scope": "REGISTERED_CTEST_COMPLETE_EXECUTION_GATE_ONLY_NOT_PAPER_CASE_OR_AUTHOR_ASSET_EVIDENCE",
        }
    if "shared.independent_moments" in requests:
        paths = repeated_paths("shared.independent_moments", "independent mass moments")
        python_root = str(Path(__file__).resolve().parents[3] / "python")
        sys.path.insert(0, python_root)
        try:
            from make_it_stand.verification.run_evidence import PaperRunEvidenceV2, check_freshness
        finally:
            sys.path.pop(0)
        scope = (
            "SYNTHETIC_ANALYTIC_AND_COMPOSITE_INDEPENDENT_MASS_FIRST_SECOND_AND_"
            "CENTRAL_INERTIA_AUDIT_ONLY_NOT_AUTHOR_CASE"
        )
        signatures = []
        for path in paths:
            evidence = PaperRunEvidenceV2.from_json(path.read_text(encoding="utf-8"))
            if check_freshness(evidence).status != "CURRENT":
                raise ValueError("STALE_RESULT: independent mass-moment evidence is not current")
            root = path.parent
            scientific = _load_json(root / "scientific.json")
            summary = _load_json(root / "summary.json")
            outputs = {Path(row["path"]).resolve() for row in evidence.as_dict()["scientific"]["outputs"]}
            scientific_path = (root / "scientific.json").resolve()
            inputs_path = (root / "inputs.json").resolve()
            if scientific_path not in outputs or inputs_path not in outputs:
                raise ValueError("independent mass-moment evidence does not bind its scientific inputs")
            if scientific.get("schema") != "shared-independent-moments-scientific-v1" or \
                    scientific.get("classification") != "SYNTHETIC_REGRESSION" or \
                    scientific.get("scope") != scope or \
                    scientific.get("integrator") != "signed_tetrahedron_python_independent" or \
                    scientific.get("density_applied_once") is not True or \
                    scientific.get("non_diagonal_inertia_fixture") is not True or \
                    summary.get("status") != "SYNTHETIC_INDEPENDENT_MASS_MOMENTS_VALIDATED" or \
                    summary.get("original_evidence_admission") != "BLOCKED_MISSING_AUTHOR_ASSETS_AND_PARAMETERS" or \
                    summary.get("paper_case_admission") != "BLOCKED_MISSING_AUTHOR_ASSETS_AND_PARAMETERS":
                raise ValueError("independent mass-moment evidence has invalid scope or status")
            fixtures = scientific.get("fixtures")
            expected_names = {
                "unit_cube_density_2", "unit_tetrahedron_density_3",
                "two_translated_cubes_non_diagonal_composite",
            }
            if not isinstance(fixtures, list) or {row.get("name") for row in fixtures} != expected_names:
                raise ValueError("independent mass-moment fixture matrix is incomplete or duplicated")
            maximum_error = scientific.get("maximum_absolute_error")
            tolerance = scientific.get("tolerance")
            if not isinstance(maximum_error, (int, float)) or not math.isfinite(maximum_error) or \
                    not isinstance(tolerance, (int, float)) or tolerance != 5e-13 or \
                    maximum_error < 0 or maximum_error > tolerance:
                raise ValueError("independent mass-moment error exceeds its registered tolerance")
            for row in fixtures:
                if not isinstance(row, dict) or row.get("maximum_absolute_error") != maximum_error and \
                        not (isinstance(row.get("maximum_absolute_error"), (int, float)) and
                             0 <= row["maximum_absolute_error"] <= tolerance):
                    raise ValueError("independent mass-moment fixture error is malformed")
            signatures.append(hashlib.sha256(
                json.dumps(scientific, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")
            ).hexdigest())
        if signatures[0] != signatures[1]:
            raise ValueError("repeated independent mass-moment scientific outputs differ")
        record = {
            "checker": {"path": str(Path(__file__).resolve()), "sha256": _sha256_file(Path(__file__)),
                        "size": Path(__file__).stat().st_size},
            "evidence": [{"path": str(path), "sha256": _sha256_file(path), "size": path.stat().st_size}
                         for path in paths],
            "scientific_data_sha256": signatures[0],
            "scope": scope,
        }
        verified["shared.independent_moments"] = record
    if "faster.section6_6" in requests:
        paths = repeated_paths("faster.section6_6", "Section 6.6")
        python_root = str(Path(__file__).resolve().parents[3] / "python")
        sys.path.insert(0, python_root)
        try:
            from make_it_stand.verification import faster_section_6_6 as exact
        finally:
            sys.path.pop(0)
        observations = []
        for path in paths:
            assessment = exact.check_evidence_freshness(path)
            if assessment["status"] != "CURRENT":
                raise ValueError(assessment["status"] + ": " + "; ".join(assessment["mismatches"]))
            evidence = exact.load_spec(path)
            result = exact.load_spec(evidence["outputs"]["result"]["path"])
            observations.append((evidence, result))
        first, second = observations
        if Path(first[0]["outputs"]["result"]["path"]).samefile(second[0]["outputs"]["result"]["path"]):
            raise ValueError("two evidence files refer to one result")
        if exact.canonical_json(first[1]) != exact.canonical_json(second[1]):
            raise ValueError("repeated scientific results differ")
        for key in ("python_executable", "python_version"):
            if not first[0].get("runtime", {}).get(key) or first[0]["runtime"][key] != second[0].get("runtime", {}).get(key):
                raise ValueError("repeat runtime identity mismatch")
        record = {"checker": exact.file_record(Path(exact.__file__)),
                  "evidence": [exact.file_record(p) for p in paths],
                  "scientific_data_sha256": exact.scientific_sha256(first[1]["scientific_data"]),
                  "scope": "PUBLISHED_SECTION_6_6_EXACT_ANALYTIC_EXAMPLES_ONLY"}
        # Section 6.6 is the published, explicit p=0 branch evidence.  Its
        # scope remains limited to the two stated examples: this records their
        # degenerate preconditions/outcomes without inventing a universal
        # classification of all possible p=0 problems.
        verified.update({name: dict(record) for name in
                         ("faster.section6_6.example.1", "faster.section6_6.example.2",
                          "faster.theorem.degenerate_conditions")})

    if "faster.discrete_moments" in requests:
        paths = repeated_paths("faster.discrete_moments", "Faster discrete moment map")
        python_root = str(Path(__file__).resolve().parents[3] / "python")
        sys.path.insert(0, python_root)
        try:
            from make_it_stand.verification.run_evidence import PaperRunEvidenceV2, check_freshness
        finally:
            sys.path.pop(0)
        markers = (
            "historical-sphere supplemental_richardson_entry_relative=",
            "contained-positive supplemental_richardson_entry_relative=",
            "clipped-positive supplemental_richardson_entry_relative=",
            "contained-negative supplemental_richardson_entry_relative=",
            "clipped-negative supplemental_richardson_entry_relative=",
            "topology-transition below=",
            "topology-open-stratum label=below",
            "topology-open-stratum label=above",
            "hyperboloid cap-radius requested=",
            "hyperboloid discrete-map h=3.0000000000000001e-05 full_map_relative=",
            "hyperboloid discrete-map h=1.0000000000000001e-05 full_map_relative=",
            "hyperboloid discrete-map h=3.0000000000000001e-06 full_map_relative=",
            "[PASS] discrete quadric moment differential",
        )
        signatures = []
        for path in paths:
            evidence = PaperRunEvidenceV2.from_json(path.read_text(encoding="utf-8"))
            modules = evidence.as_dict()["begin"]["observation"]["runtime"]["native_modules"]
            module = modules.get("_make_it_stand") if isinstance(modules, dict) else None
            if not isinstance(module, dict):
                raise ValueError("Faster discrete-moment evidence lacks native module identity")
            module_path = Path(module.get("path", "")).resolve(strict=True)
            if module.get("sha256") != _sha256_file(module_path) or module.get("size") != module_path.stat().st_size:
                raise ValueError("STALE_RESULT: Faster discrete-moment native module differs from observation")
            if "_make_it_stand" not in sys.modules:
                sys.path.insert(0, str(module_path.parent))
                try:
                    importlib.import_module("_make_it_stand")
                finally:
                    sys.path.pop(0)
            if check_freshness(evidence).status != "CURRENT":
                raise ValueError("STALE_RESULT: Faster discrete-moment evidence is not current")
            root = path.parent
            scientific = _load_json(root / "scientific.json")
            summary = _load_json(root / "summary.json")
            if scientific.get("schema") != "faster-discrete-moment-scientific-v2" or \
                    scientific.get("classification") != "SYNTHETIC_REGRESSION" or \
                    summary.get("status") != "SYNTHETIC_DISCRETE_FULLRANK_QUADRIC_METHOD_VALIDATED" or \
                    summary.get("original_evidence_admission") != "BLOCKED_UNKNOWN" or \
                    summary.get("paper_case_admission") != "BLOCKED_MISSING_AUTHOR_ASSETS_AND_PARAMETERS":
                raise ValueError("Faster discrete-moment evidence has invalid scope")
            outputs = {Path(row["path"]).resolve() for row in evidence.as_dict()["scientific"]["outputs"]}
            transcript = root / "discrete-moments.stdout"
            required_outputs = {root / "scientific.json", transcript}
            if not {item.resolve() for item in required_outputs} <= outputs:
                raise ValueError("Faster discrete-moment evidence does not bind all scientific outputs")
            result = scientific.get("native_regression")
            if not isinstance(result, dict) or result.get("returncode") != 0 or \
                    result.get("required_markers") != list(markers):
                raise ValueError("Faster discrete-moment native result is malformed")
            recorded = Path(result.get("path", "")).resolve()
            if recorded != transcript.resolve() or recorded not in outputs or not recorded.is_file():
                raise ValueError("Faster discrete-moment transcript is not bound")
            data = recorded.read_bytes()
            if hashlib.sha256(data).hexdigest() != result.get("sha256") or \
                    any(marker not in data.decode("utf-8") for marker in markers):
                raise ValueError("Faster discrete-moment transcript does not establish the native result")
            signatures.append(result["sha256"])
        if signatures[0] != signatures[1]:
            raise ValueError("repeated Faster discrete-moment scientific outputs differ")
        record = {
            "checker": {"path": str(Path(__file__).resolve()), "sha256": _sha256_file(Path(__file__)),
                        "size": Path(__file__).stat().st_size},
            "evidence": [{"path": str(path), "sha256": _sha256_file(path), "size": path.stat().st_size}
                         for path in paths],
            "scientific_data_sha256": signatures[0],
            "scope": "SYNTHETIC_FIXED_TOPOLOGY_DISCRETE_MOMENT_AND_INVERSE_GAUSS_METHODS_ONLY_NOT_AUTHOR_CASE_OR_INDEPENDENT_BOOLEAN_BACKEND",
        }
        verified.update({name: dict(record) for name in
                         ("faster.discrete_moment_map.differentiable", "faster.geometry.inverse_gauss_map")})

    if "faster.reference_boolean" in requests:
        paths = repeated_paths("faster.reference_boolean", "Faster reference Boolean")
        python_root = str(Path(__file__).resolve().parents[3] / "python")
        sys.path.insert(0, python_root)
        try:
            from make_it_stand.verification.run_evidence import PaperRunEvidenceV2, check_freshness
        finally:
            sys.path.pop(0)
        signatures = []
        for path in paths:
            evidence = PaperRunEvidenceV2.from_json(path.read_text(encoding="utf-8"))
            modules = evidence.as_dict()["begin"]["observation"]["runtime"]["native_modules"]
            module = modules.get("_make_it_stand") if isinstance(modules, dict) else None
            if not isinstance(module, dict):
                raise ValueError("Faster reference-Boolean evidence lacks native module identity")
            module_path = Path(module.get("path", "")).resolve(strict=True)
            if module.get("sha256") != _sha256_file(module_path) or module.get("size") != module_path.stat().st_size:
                raise ValueError("STALE_RESULT: Faster reference-Boolean native module differs from observation")
            try:
                loaded, freshness = _load_bound_native_module_and_freshness(
                    evidence, module_path, check_freshness)
            except (ImportError, OSError, ValueError) as error:
                raise ValueError("Faster reference-Boolean cannot load the evidenced native module") from error
            if freshness.status != "CURRENT":
                raise ValueError("STALE_RESULT: Faster reference-Boolean evidence is not current")
            root = path.parent
            scientific = _load_json(root / "scientific.json")
            summary = _load_json(root / "summary.json")
            if scientific.get("schema") != "faster-reference-boolean-crosscheck-v2" or \
                    scientific.get("classification") != "SYNTHETIC_REGRESSION" or \
                    scientific.get("scope") != "SAME_OPERAND_CLIPPED_QUADRIC_BOOLEAN_ONLY_NOT_AUTHOR_ASSET_OR_PAPER_CASE" or \
                    summary.get("original_evidence_admission") != "BLOCKED_MISSING_AUTHOR_ASSETS_AND_PARAMETERS" or \
                    summary.get("paper_case_admission") != "BLOCKED_MISSING_AUTHOR_ASSETS_AND_PARAMETERS":
                raise ValueError("Faster reference-Boolean evidence has invalid scope")
            outputs = {Path(row["path"]).resolve() for row in evidence.as_dict()["scientific"]["outputs"]}
            expected_operations = {
                name: _reference_boolean_expected_operation(name)
                for name in ("ellipsoid", "one_sheet_hyperboloid", "two_sheet_hyperboloid",
                             "plane", "paraboloid")
            }
            required = {root / "scientific.json", root / "design.obj"}
            for name in expected_operations:
                required.update({root / f"quadric-{name}.obj", root / f"cherchi-{name}-run-1.obj",
                                 root / f"cherchi-{name}-run-2.obj", root / f"cgal-material-{name}.obj",
                                 root / f"cherchi-{name}-provenance.obj",
                                 root / f"cherchi-{name}-provenance.obj.exact.tsv"})
            if not {item.resolve() for item in required} <= outputs:
                raise ValueError("Faster reference-Boolean evidence does not bind all scientific outputs")
            fixtures = scientific.get("fixtures")
            if not isinstance(fixtures, list) or scientific.get("design_faces") != 12:
                raise ValueError("Faster reference-Boolean scientific record is malformed")
            by_name = {row.get("name"): row for row in fixtures if isinstance(row, dict)}
            if set(by_name) != set(expected_operations):
                raise ValueError("Faster reference-Boolean fixture matrix is incomplete or duplicated")
            fixture_signatures = []
            geometry_statuses = {}
            design_vertices, design_faces = _obj_mesh_arrays(root / "design.obj")
            design_region = loaded.make_exact_material_region(design_vertices, design_faces)
            for name, operation in expected_operations.items():
                fixture = by_name[name]
                runs = fixture.get("cherchi_runs")
                if fixture.get("operation") != operation or fixture.get("quadric_faces", 0) <= 0 or \
                        fixture.get("cgal_material_faces", 0) <= scientific.get("design_faces", 0) or \
                        not isinstance(fixture.get("cgal_raw_moments"), list) or len(fixture["cgal_raw_moments"]) != 10 or \
                        not isinstance(runs, list) or len(runs) != 2:
                    raise ValueError("Faster reference-Boolean fixture is malformed")
                geometry = fixture.get("geometry_comparison")
                cgal_mesh = root / f"cgal-material-{name}.obj"
                if not isinstance(geometry, dict) or \
                        geometry.get("cgal_mesh") != str(cgal_mesh) or \
                        cgal_mesh.resolve() not in outputs or not cgal_mesh.is_file() or \
                        _sha256_file(cgal_mesh) != geometry.get("cgal_mesh_sha256"):
                    raise ValueError("Faster reference-Boolean evidence lacks a bound CGAL material mesh")
                cgal_geometry_signature = _oriented_triangle_boundary_sha256(cgal_mesh)
                if cgal_geometry_signature != geometry.get("cgal_signature_sha256"):
                    raise ValueError("Faster reference-Boolean CGAL triangle boundary signature is invalid")
                cgal_vertices, cgal_faces = _obj_mesh_arrays(cgal_mesh)
                result_hashes = []
                for index, row in enumerate(runs, start=1):
                    if not isinstance(row, dict) or row.get("path") != str(root / f"cherchi-{name}-run-{index}.obj"):
                        raise ValueError("Faster reference-Boolean output path is malformed")
                    result_path = Path(row["path"]).resolve()
                    if result_path not in outputs or not result_path.is_file() or \
                            _sha256_file(result_path) != row.get("sha256") or \
                            row.get("face_count") <= 0 or row.get("vertex_count") <= 0 or \
                            not isinstance(row.get("raw_moments"), list) or len(row["raw_moments"]) != 10 or \
                            not isinstance(row.get("maximum_absolute_error"), (int, float)) or \
                            row["maximum_absolute_error"] > 1e-12:
                        raise ValueError("Faster reference-Boolean output is not a valid same-operand comparison")
                    cherchi_geometry_signature = _oriented_triangle_boundary_sha256(result_path)
                    if cherchi_geometry_signature != geometry.get("cherchi_signature_sha256"):
                        raise ValueError("Faster reference-Boolean Cherchi triangle boundary signature is invalid")
                    result_hashes.append(row["sha256"])
                reference_path = Path(runs[0]["path"]).resolve()
                reference_vertices, reference_faces = _obj_mesh_arrays(reference_path)
                cgal_region = loaded.make_exact_material_region(cgal_vertices, cgal_faces)
                cherchi_region = loaded.make_exact_material_region(reference_vertices, reference_faces)
                rational_relation = _reference_boolean_exact_construction_relation(
                    loaded, root / f"cherchi-{name}-provenance.obj.exact.tsv",
                    reference_vertices, reference_faces, cgal_region)
                if geometry.get("exact_construction_relation") != rational_relation:
                    raise ValueError("Faster reference-Boolean rational construction relation is not independently reproducible")
                input_oracle = _reference_boolean_input_oracle(
                    loaded, root, name, operation)
                if geometry.get("exact_input_boolean_oracle") != input_oracle:
                    raise ValueError("Faster reference-Boolean exact input operation is not independently reproducible")
                for construction_key in ("cgal_construction", "cherchi_construction",
                                         "cgal_cherchi_exact_construction_relation"):
                    construction = input_oracle.get(construction_key)
                    if not isinstance(construction, dict) or \
                            construction.get("completed") is not True or \
                            construction.get("equal") is not True or \
                            construction.get("symmetric_difference_volume_exact") != "0/1":
                        raise ValueError(
                            f"Faster reference-Boolean internal exact construction fails for {name}: "
                            f"{construction_key}")
                exact_relation = dict(loaded.relate_exact_material_regions(cgal_region, cherchi_region))
                if exact_relation != geometry.get("exact_relation"):
                    raise ValueError("Faster reference-Boolean exact geometry relation is not independently reproducible")
                volume_comparison = _exact_symmetric_difference_volume_comparison(
                    exact_relation, cgal_vertices, reference_vertices)
                if volume_comparison != geometry.get("symmetric_difference_volume_comparison"):
                    raise ValueError("Faster reference-Boolean symmetric-difference volume comparison is not reproducible")
                if exact_relation.get("completed") is not True:
                    geometry_status = "UNKNOWN"
                elif exact_relation.get("equal") is True and exact_relation.get("a_subset_b") is True and \
                        exact_relation.get("b_subset_a") is True:
                    geometry_status = "SAME_EXACT_POLYHEDRAL_SOLID"
                elif exact_relation.get("equal") is False:
                    geometry_status = "DIFFERENT_EXACT_POLYHEDRAL_SOLIDS"
                else:
                    geometry_status = "UNKNOWN"
                if geometry.get("status") != geometry_status or \
                        geometry.get("same_oriented_triangle_set") != (cgal_geometry_signature == cherchi_geometry_signature):
                    raise ValueError("Faster reference-Boolean geometry status disagrees with exact mesh checks")
                geometry_statuses[name] = geometry_status
                if geometry_status == "UNKNOWN":
                    raise ValueError(f"Faster reference-Boolean exported polyhedral comparison is unknown for {name}")
                quadric_vertices, quadric_faces = _obj_mesh_arrays(root / f"quadric-{name}.obj")
                quadric_region = loaded.make_exact_material_region(quadric_vertices, quadric_faces)
                independently_checked_operands = {
                    "cherchi": _reference_boolean_operand_conditions(
                        loaded, operation, design_region, quadric_region, cherchi_region),
                    "cgal": _reference_boolean_operand_conditions(
                        loaded, operation, design_region, quadric_region, cgal_region),
                }
                if geometry.get("operand_necessary_conditions") != independently_checked_operands:
                    raise ValueError("Faster reference-Boolean input-operand conditions are not independently reproducible")
                if result_hashes[0] != result_hashes[1] or runs[0]["raw_moments"] != runs[1]["raw_moments"]:
                    raise ValueError("Faster reference-Boolean repetitions differ")
                fixture_signatures.append((name, result_hashes[0]))
            all_exports_match = all(value == "SAME_EXACT_POLYHEDRAL_SOLID"
                                    for value in geometry_statuses.values())
            expected_geometry = "EXACT_SOLID_MATCH" if all_exports_match else "EXACT_SOLID_DISCREPANCY"
            expected_status = ("SYNTHETIC_REFERENCE_BOOLEAN_CROSSCHECK_VALIDATED" if all_exports_match
                               else "SYNTHETIC_REFERENCE_BOOLEAN_GEOMETRY_DISCREPANCY")
            expected_discrepancies = sorted(name for name, value in geometry_statuses.items()
                                            if value == "DIFFERENT_EXACT_POLYHEDRAL_SOLIDS")
            if scientific.get("geometry_status") != expected_geometry or \
                    summary.get("status") != expected_status or \
                    summary.get("geometry_status") != expected_geometry or \
                    sorted(summary.get("geometry_discrepancy_fixtures", [])) != expected_discrepancies:
                raise ValueError("Faster reference-Boolean exported geometry status is inconsistent")
            signatures.append(tuple(fixture_signatures))
        if signatures[0] != signatures[1]:
            raise ValueError("repeated Faster reference-Boolean scientific outputs differ")
        matrix_sha256 = hashlib.sha256(
            json.dumps(signatures[0], separators=(",", ":")).encode("ascii")
        ).hexdigest()
        record = {
            "checker": {"path": str(Path(__file__).resolve()), "sha256": _sha256_file(Path(__file__)),
                        "size": Path(__file__).stat().st_size},
            "evidence": [{"path": str(path), "sha256": _sha256_file(path), "size": path.stat().st_size}
                         for path in paths],
            "scientific_data_sha256": matrix_sha256,
            "scope": "SYNTHETIC_SAME_OPERAND_INTERNAL_EXACT_CONSTRUCTION_ONLY;BINARY64_EXPORT_SEPARATELY_BLOCKED;NOT_AUTHOR_CASE",
        }
        verified["faster.geometry.boolean_independent"] = record

    if "faster.moment_foundations" in requests:
        paths = repeated_paths("faster.moment_foundations", "Faster moment foundations")
        python_root = str(Path(__file__).resolve().parents[3] / "python")
        sys.path.insert(0, python_root)
        try:
            from make_it_stand.verification.run_evidence import PaperRunEvidenceV2, check_freshness
        finally:
            sys.path.pop(0)
        signatures = []
        expected = {
            "raw-moments": ((), (
                "[PASS] raw ten-moment volume and first-moment convention",
                "[PASS] scenario=quadric-moment-evaluator derivative_error=",
            )),
            "dual-monomial": (("--dual-coordinates",), (
                "[PASS] literal coefficient/raw-moment cross-term duality",
                "[PASS] literal coefficient-to-geometry differential",
                "[PASS] dual monomial coordinate Newton",
            )),
        }
        for path in paths:
            evidence = PaperRunEvidenceV2.from_json(path.read_text(encoding="utf-8"))
            if check_freshness(evidence).status != "CURRENT":
                raise ValueError("STALE_RESULT: Faster moment-foundation evidence is not current")
            root = path.parent
            scientific, summary = _load_json(root / "scientific.json"), _load_json(root / "summary.json")
            outputs = {Path(row["path"]).resolve() for row in evidence.as_dict()["scientific"]["outputs"]}
            commands = scientific.get("commands")
            if scientific.get("schema") != "faster-moment-foundations-scientific-v1" or \
                    scientific.get("classification") != "SYNTHETIC_REGRESSION" or \
                    scientific.get("scope") != "SYNTHETIC_RAW_TEN_MOMENT_AND_LITERAL_CROSS_TERM_COEFFICIENT_DUALITY_METHODS_ONLY_NOT_FASTER_THEOREMS_OR_AUTHOR_CASE" or \
                    summary.get("status") != "SYNTHETIC_MOMENT_FOUNDATIONS_VALIDATED" or \
                    summary.get("theorem_coverage") != "NOT_GRANTED_BY_THIS_RUNNER" or \
                    not isinstance(commands, list) or len(commands) != len(expected):
                raise ValueError("Faster moment-foundation evidence is malformed")
            observed = {}
            for command in commands:
                if not isinstance(command, dict) or command.get("name") not in expected:
                    raise ValueError("Faster moment-foundation command is malformed")
                name = command["name"]
                arguments, markers = expected[name]
                transcript = root / f"{name}.stdout"
                if name in observed or command.get("arguments") != list(arguments) or \
                        command.get("required_markers") != list(markers) or command.get("returncode") != 0 or \
                        Path(command.get("path", "")).resolve() != transcript.resolve() or \
                        transcript.resolve() not in outputs or not transcript.is_file() or \
                        _sha256_file(transcript) != command.get("sha256"):
                    raise ValueError("Faster moment-foundation transcript is not bound")
                text = transcript.read_text(encoding="utf-8")
                if any(marker not in text for marker in markers):
                    raise ValueError("Faster moment-foundation transcript lacks a required marker")
                observed[name] = command["sha256"]
            if set(observed) != set(expected):
                raise ValueError("Faster moment-foundation command set is incomplete")
            signatures.append(tuple(sorted(observed.items())))
        if signatures[0] != signatures[1]:
            raise ValueError("repeated Faster moment-foundation scientific outputs differ")
        record = {
            "checker": {"path": str(Path(__file__).resolve()), "sha256": _sha256_file(Path(__file__)),
                        "size": Path(__file__).stat().st_size},
            "evidence": [{"path": str(path), "sha256": _sha256_file(path), "size": path.stat().st_size}
                         for path in paths],
            "scientific_data_sha256": hashlib.sha256(
                json.dumps(signatures[0], separators=(",", ":")).encode("ascii")
            ).hexdigest(),
            "scope": "SYNTHETIC_RAW_TEN_MOMENT_AND_LITERAL_CROSS_TERM_COEFFICIENT_DUALITY_METHODS_ONLY_NOT_FASTER_THEOREMS_OR_AUTHOR_CASE",
        }
        verified.update({name: dict(record) for name in (
            "faster.moments.raw", "faster.moments.coefficient_duality",
        )})

    if "faster.hydrostatic_geometry" in requests:
        paths = repeated_paths("faster.hydrostatic_geometry", "Faster hydrostatic geometry")
        python_root = str(Path(__file__).resolve().parents[3] / "python")
        sys.path.insert(0, python_root)
        try:
            from make_it_stand.verification.run_evidence import PaperRunEvidenceV2, check_freshness
        finally:
            sys.path.pop(0)
        signatures = []
        for path in paths:
            evidence = PaperRunEvidenceV2.from_json(path.read_text(encoding="utf-8"))
            if check_freshness(evidence).status != "CURRENT":
                raise ValueError("STALE_RESULT: Faster hydrostatic geometry evidence is not current")
            root = path.parent
            scientific, summary = _load_json(root / "scientific.json"), _load_json(root / "summary.json")
            result = scientific.get("native_regression")
            outputs = {Path(row["path"]).resolve() for row in evidence.as_dict()["scientific"]["outputs"]}
            transcript = root / "hydrostatic.stdout"
            if scientific.get("schema") != "faster-hydrostatic-geometry-scientific-v2" or \
                    scientific.get("classification") != "SYNTHETIC_REGRESSION" or \
                    scientific.get("scope") != "SYNTHETIC_CLOSED_DISPLACEMENT_BOUNDARY_ACTUAL_WATERLINE_AND_VOLUME_EQUILIBRIUM_METHODS_ONLY_NOT_AUTHOR_FLOATING_CASE" or \
                    summary.get("status") != "SYNTHETIC_HYDROSTATIC_GEOMETRY_METHOD_VALIDATED" or \
                    not isinstance(result, dict) or result.get("returncode") != 0 or \
                    result.get("required_marker") != "[PASS] hydrostatic geometry actual mesh clipping and volume equilibrium" or \
                    Path(result.get("path", "")).resolve() != transcript.resolve() or \
                    transcript.resolve() not in outputs or not transcript.is_file() or \
                    _sha256_file(transcript) != result.get("sha256") or \
                    result["required_marker"] not in transcript.read_text(encoding="utf-8"):
                raise ValueError("Faster hydrostatic geometry evidence is malformed")
            signatures.append(result["sha256"])
        if signatures[0] != signatures[1]:
            raise ValueError("repeated Faster hydrostatic geometry outputs differ")
        record = {"checker": {"path": str(Path(__file__).resolve()), "sha256": _sha256_file(Path(__file__)),
                              "size": Path(__file__).stat().st_size},
                  "evidence": [{"path": str(path), "sha256": _sha256_file(path), "size": path.stat().st_size}
                               for path in paths], "scientific_data_sha256": signatures[0],
                  "scope": "SYNTHETIC_CLOSED_DISPLACEMENT_BOUNDARY_ACTUAL_WATERLINE_AND_VOLUME_EQUILIBRIUM_METHODS_ONLY_NOT_AUTHOR_FLOATING_CASE"}
        verified["faster.hydrostatics.actual"] = record

    if "faster.planar_methods" in requests:
        paths = repeated_paths("faster.planar_methods", "Faster planar methods")
        python_root = str(Path(__file__).resolve().parents[3] / "python")
        sys.path.insert(0, python_root)
        try:
            from make_it_stand.verification.run_evidence import PaperRunEvidenceV2, check_freshness
        finally:
            sys.path.pop(0)
        expected = {
            "dogleg": ("planar-fidget-dogleg-regression", (
                "FASTER_PLANAR_COMPONENT_JSON=",
            )),
            "polygon_disk": ("planar-disk-intersection-regression", (
                "FASTER_PLANAR_COMPONENT_JSON=",
            )),
            "moment_problem": ("planar-spin-moment-problem-regression", ()),
            "mass_moments": ("planar-mass-moments-regression", ()),
            "extrusion": ("planar-extrusion-regression", (
                "FASTER_PLANAR_COMPONENT_JSON=",
                "[PASS] planar ring moments equal independently integrated 3D extrusion",
                "[PASS] planar and 3D fidget objectives/constraints agree",
            )),
        }
        signatures = []
        for path in paths:
            evidence = PaperRunEvidenceV2.from_json(path.read_text(encoding="utf-8"))
            if check_freshness(evidence).status != "CURRENT":
                raise ValueError("STALE_RESULT: Faster planar-method evidence is not current")
            root = path.parent
            scientific, summary = _load_json(root / "scientific.json"), _load_json(root / "summary.json")
            outputs = {Path(row["path"]).resolve() for row in evidence.as_dict()["scientific"]["outputs"]}
            if scientific.get("schema") != "faster-planar-methods-scientific-v1" or \
                    scientific.get("classification") != "SYNTHETIC_REGRESSION" or \
                    scientific.get("scope") != "SYNTHETIC_PLANAR_DOGLEG_POLYGON_DISK_MOMENT_AND_EXTRUSION_METHODS_ONLY_NOT_AUTHOR_2D_CASES" or \
                    summary.get("status") != "SYNTHETIC_PLANAR_METHODS_VALIDATED":
                raise ValueError("Faster planar-method evidence has invalid scope")
            commands = scientific.get("commands")
            if not isinstance(commands, list) or {row.get("name") for row in commands} != set(expected):
                raise ValueError("Faster planar-method evidence misses a protected regression")
            normalized = []
            for row in commands:
                name = row["name"]
                binary_name, markers = expected[name]
                transcript = Path(row.get("path", ""))
                binary = Path(row.get("binary", ""))
                capture = row.get("capture")
                if binary.name != binary_name or row.get("returncode") != 0 or \
                        row.get("required_markers") != list(markers) or transcript.resolve() not in outputs or \
                        not transcript.is_file() or _sha256_file(transcript) != row.get("sha256") or \
                        not binary.is_file() or _sha256_file(binary) != row.get("binary_sha256"):
                    raise ValueError("Faster planar-method command is malformed")
                content = transcript.read_text(encoding="utf-8")
                if markers:
                    if capture != "VERBATIM_STDOUT_STDERR_V1" or any(marker not in content for marker in markers):
                        raise ValueError("Faster planar-method transcript lacks its required marker")
                elif capture != "EMPTY_STDOUT_STDERR_RECEIPT_V1" or _load_json(transcript) != {
                        "schema": "faster-planar-empty-output-receipt-v1",
                        "raw_stdout_stderr_bytes": 0,
                        "raw_stdout_stderr_sha256": hashlib.sha256(b"").hexdigest(),
                        "returncode": 0}:
                    raise ValueError("Faster planar-method transcript lacks its required marker")
                normalized.append((name, row["sha256"], row["binary_sha256"]))
            raw_components = scientific.get("component_results")
            if not isinstance(raw_components, dict) or set(raw_components) != {
                    "dogleg", "polygon_disk", "extrusion"}:
                raise ValueError("Faster planar evidence lacks per-component numerical results")
            components = {}
            for name in ("dogleg", "polygon_disk", "extrusion"):
                result = _validate_faster_planar_component_result(name, raw_components[name])
                transcript = root / f"{name}.stdout"
                marker_lines = [line for line in transcript.read_text(encoding="utf-8").splitlines()
                                if line.startswith("FASTER_PLANAR_COMPONENT_JSON=")]
                if len(marker_lines) != 1:
                    raise ValueError("Faster planar component transcript marker is absent or duplicated")
                try:
                    transcript_result = json.loads(marker_lines[0].split("=", 1)[1])
                except json.JSONDecodeError as error:
                    raise ValueError("Faster planar component transcript is malformed JSON") from error
                _reject_nonfinite(transcript_result, f"Faster planar transcript {name}")
                if transcript_result != result:
                    raise ValueError("Faster planar scientific result differs from its native transcript")
                components[name] = result
            normalized_signature = tuple(sorted(normalized))
            component_signature = json.dumps(components, sort_keys=True, separators=(",", ":"),
                                             allow_nan=False)
            signatures.append((normalized_signature, component_signature))
        if signatures[0] != signatures[1]:
            raise ValueError("repeated Faster planar-method scientific outputs differ")
        component_results = json.loads(signatures[0][1])
        record = {"checker": {"path": str(Path(__file__).resolve()), "sha256": _sha256_file(Path(__file__)),
                              "size": Path(__file__).stat().st_size},
                  "evidence": [{"path": str(path), "sha256": _sha256_file(path), "size": path.stat().st_size}
                               for path in paths], "scientific_data_sha256": hashlib.sha256(
                      json.dumps(signatures[0], separators=(",", ":")).encode("ascii")).hexdigest(),
                  "scope": "SYNTHETIC_PLANAR_DOGLEG_POLYGON_DISK_MOMENT_AND_EXTRUSION_METHODS_ONLY_NOT_AUTHOR_2D_CASES"}
        verified.update({name: dict(record) for name in
                         ("faster.2d.methods",
                          "faster.application.fidget_spinner")})
        component_scopes = {
            "dogleg": ("faster.2d.dogleg",
                       "SYNTHETIC_2D_DOGLEG_SOLVER_CONTRACT_ONLY_NOT_AUTHOR_FIDGET_CASE"),
            "polygon_disk": ("faster.2d.polygon_circle",
                             "SYNTHETIC_POLYGON_DISK_MOMENT_INTEGRATION_CONTRACT_ONLY_NOT_AUTHOR_FIDGET_CASE"),
            "extrusion": ("faster.2d.extrusion",
                          "SYNTHETIC_PLANAR_TO_3D_EXTRUSION_CROSSCHECK_ONLY_NOT_AUTHOR_FIDGET_CASE"),
        }
        for name, (requirement_id, scope) in component_scopes.items():
            component_record = dict(record)
            component_record["scope"] = scope
            component_record["component_result"] = component_results[name]
            component_record["scientific_data_sha256"] = hashlib.sha256(
                json.dumps(component_results[name], sort_keys=True, separators=(",", ":"),
                           allow_nan=False).encode("ascii")
            ).hexdigest()
            verified[requirement_id] = component_record

    if "faster.static_application_contract" in requests:
        paths = repeated_paths(
            "faster.static_application_contract", "Faster static-stability application contract"
        )
        python_root = str(Path(__file__).resolve().parents[3] / "python")
        sys.path.insert(0, python_root)
        try:
            from make_it_stand.verification.run_evidence import PaperRunEvidenceV2, check_freshness
        finally:
            sys.path.pop(0)
        scope = ("SYNTHETIC_STATIC_STABILITY_EQ7_EQ13_COEFFICIENT_CONDITIONS_ONLY_"
                 "NOT_AUTHOR_CASE")
        expected_markers = (
            "[PASS] static maximum-volume target frame, first-moment constraints, Eq. 7, and vertical-plane coefficients agree independently",
            "[PASS] static lowering target frame, first-moment constraints, Eq. 13, p!=0, and center-crossing nonvertical plane agree independently",
            "[PASS] static application scope is synthetic method validation only",
        )
        signatures = []
        for path in paths:
            evidence = PaperRunEvidenceV2.from_json(path.read_text(encoding="utf-8"))
            if check_freshness(evidence).status != "CURRENT":
                raise ValueError("STALE_RESULT: Faster static-stability application evidence is not current")
            root = path.parent
            scientific, summary = _load_json(root / "scientific.json"), _load_json(root / "summary.json")
            outputs = {Path(row["path"]).resolve() for row in evidence.as_dict()["scientific"]["outputs"]}
            command = scientific.get("command")
            transcript = root / "static-contract.stdout"
            if (
                scientific.get("schema") != "faster-static-stability-contract-scientific-v1" or
                scientific.get("classification") != "SYNTHETIC_REGRESSION" or
                scientific.get("scope") != scope or
                summary.get("status") != "SYNTHETIC_STATIC_STABILITY_APPLICATION_CONTRACT_VALIDATED" or
                summary.get("scope") != scope or
                not isinstance(command, dict) or
                command.get("name") != "static-contract" or
                Path(command.get("path", "")).resolve() != transcript.resolve() or
                transcript.resolve() not in outputs or not transcript.is_file() or
                command.get("returncode") != 0 or
                command.get("required_markers") != list(expected_markers) or
                command.get("sha256") != _sha256_file(transcript) or
                Path(command.get("binary", "")).name != "faster-static-stability-contract-regression" or
                not Path(command.get("binary", "")).is_file() or
                command.get("binary_sha256") != _sha256_file(Path(command["binary"]))
            ):
                raise ValueError("Faster static-stability application evidence is malformed")
            text = transcript.read_text(encoding="utf-8")
            if any(marker not in text for marker in expected_markers):
                raise ValueError("Faster static-stability application transcript lacks a required marker")
            signatures.append((command["sha256"], command["binary_sha256"]))
        if signatures[0] != signatures[1]:
            raise ValueError("repeated Faster static-stability application scientific outputs differ")
        record = {"checker": {"path": str(Path(__file__).resolve()), "sha256": _sha256_file(Path(__file__)),
                              "size": Path(__file__).stat().st_size},
                  "evidence": [{"path": str(path), "sha256": _sha256_file(path), "size": path.stat().st_size}
                               for path in paths], "scientific_data_sha256": hashlib.sha256(
                      json.dumps(signatures[0], separators=(",", ":")).encode("ascii")).hexdigest(),
                  "scope": scope}
        verified["faster.application.static_stability"] = record

    if "faster.buoyant_application_contract" in requests:
        if "faster.hydrostatics.actual" not in verified:
            raise ValueError("Faster buoyant-stability application requires verified actual hydrostatic geometry")
        paths = repeated_paths(
            "faster.buoyant_application_contract", "Faster buoyant-stability application contract"
        )
        python_root = str(Path(__file__).resolve().parents[3] / "python")
        sys.path.insert(0, python_root)
        try:
            from make_it_stand.verification.run_evidence import PaperRunEvidenceV2, check_freshness
        finally:
            sys.path.pop(0)
        scope = ("SYNTHETIC_BUOYANT_STABILITY_DENSITY_SHELL_WATERLINE_METACENTER_"
                 "EQ7_AND_EQ13_CONDITIONS_ONLY_NOT_AUTHOR_CASE")
        expected_markers = (
            "[PASS] prescribed-waterline plane solve and finite-precision Eq. 7 audit passed",
        )
        signatures = []
        for path in paths:
            evidence = PaperRunEvidenceV2.from_json(path.read_text(encoding="utf-8"))
            if check_freshness(evidence).status != "CURRENT":
                raise ValueError("STALE_RESULT: Faster buoyant-stability application evidence is not current")
            root = path.parent
            scientific, summary = _load_json(root / "scientific.json"), _load_json(root / "summary.json")
            outputs = {Path(row["path"]).resolve() for row in evidence.as_dict()["scientific"]["outputs"]}
            command = scientific.get("command")
            transcript = root / "buoyant-contract.stdout"
            if (
                scientific.get("schema") != "faster-buoyant-stability-contract-scientific-v1" or
                scientific.get("classification") != "SYNTHETIC_REGRESSION" or
                scientific.get("scope") != scope or
                summary.get("status") != "SYNTHETIC_BUOYANT_STABILITY_APPLICATION_CONTRACT_VALIDATED" or
                summary.get("scope") != scope or
                not isinstance(command, dict) or
                command.get("name") != "buoyant-plane-balance-contract" or
                Path(command.get("path", "")).resolve() != transcript.resolve() or
                transcript.resolve() not in outputs or not transcript.is_file() or
                command.get("returncode") != 0 or
                command.get("required_markers") != list(expected_markers) or
                command.get("sha256") != _sha256_file(transcript) or
                Path(command.get("binary", "")).name != "buoyant-plane-balance-regression" or
                not Path(command.get("binary", "")).is_file() or
                command.get("binary_sha256") != _sha256_file(Path(command["binary"]))
            ):
                raise ValueError("Faster buoyant-stability application evidence is malformed")
            text = transcript.read_text(encoding="utf-8")
            if any(marker not in text for marker in expected_markers):
                raise ValueError("Faster buoyant-stability application transcript lacks a required marker")
            signatures.append((command["sha256"], command["binary_sha256"]))
        if signatures[0] != signatures[1]:
            raise ValueError("repeated Faster buoyant-stability application scientific outputs differ")
        record = {"checker": {"path": str(Path(__file__).resolve()), "sha256": _sha256_file(Path(__file__)),
                              "size": Path(__file__).stat().st_size},
                  "evidence": [{"path": str(path), "sha256": _sha256_file(path), "size": path.stat().st_size}
                               for path in paths], "scientific_data_sha256": hashlib.sha256(
                      json.dumps(signatures[0], separators=(",", ":")).encode("ascii")).hexdigest(),
                  "scope": scope,
                  "depends_on": "faster.hydrostatics.actual"}
        verified["faster.application.buoyant_stability"] = record

    if "faster.yoyo_application_contract" in requests:
        paths = repeated_paths(
            "faster.yoyo_application_contract", "Faster yo-yo application contract"
        )
        python_root = str(Path(__file__).resolve().parents[3] / "python")
        sys.path.insert(0, python_root)
        try:
            from make_it_stand.verification.run_evidence import PaperRunEvidenceV2, check_freshness
        finally:
            sys.path.pop(0)
        scope = ("SYNTHETIC_YOYO_SPIN_FRAME_INERTIA_BLOCK_AND_EQ7_COEFFICIENT_"
                 "CONDITIONS_ONLY_NOT_AUTHOR_CASE")
        expected_markers = (
            "[PASS] yo-yo spin-frame objective and center/inertia-block constraints agree independently",
            "[PASS] yo-yo Eq. 7 coefficient condition a=0, A33<0, and p!=0 agree independently",
            "[PASS] yo-yo application scope is synthetic method validation only",
        )
        signatures = []
        for path in paths:
            evidence = PaperRunEvidenceV2.from_json(path.read_text(encoding="utf-8"))
            if check_freshness(evidence).status != "CURRENT":
                raise ValueError("STALE_RESULT: Faster yo-yo application evidence is not current")
            root = path.parent
            scientific, summary = _load_json(root / "scientific.json"), _load_json(root / "summary.json")
            outputs = {Path(row["path"]).resolve() for row in evidence.as_dict()["scientific"]["outputs"]}
            command = scientific.get("command")
            transcript = root / "yoyo-contract.stdout"
            if (
                scientific.get("schema") != "faster-yoyo-contract-scientific-v1" or
                scientific.get("classification") != "SYNTHETIC_REGRESSION" or
                scientific.get("scope") != scope or
                summary.get("status") != "SYNTHETIC_YOYO_APPLICATION_CONTRACT_VALIDATED" or
                summary.get("scope") != scope or
                not isinstance(command, dict) or
                command.get("name") != "yoyo-contract" or
                Path(command.get("path", "")).resolve() != transcript.resolve() or
                transcript.resolve() not in outputs or not transcript.is_file() or
                command.get("returncode") != 0 or
                command.get("required_markers") != list(expected_markers) or
                command.get("sha256") != _sha256_file(transcript) or
                Path(command.get("binary", "")).name != "spin-moment-problem-regression" or
                not Path(command.get("binary", "")).is_file() or
                command.get("binary_sha256") != _sha256_file(Path(command["binary"]))
            ):
                raise ValueError("Faster yo-yo application evidence is malformed")
            text = transcript.read_text(encoding="utf-8")
            if any(marker not in text for marker in expected_markers):
                raise ValueError("Faster yo-yo application transcript lacks a required marker")
            signatures.append((command["sha256"], command["binary_sha256"]))
        if signatures[0] != signatures[1]:
            raise ValueError("repeated Faster yo-yo application scientific outputs differ")
        record = {"checker": {"path": str(Path(__file__).resolve()), "sha256": _sha256_file(Path(__file__)),
                              "size": Path(__file__).stat().st_size},
                  "evidence": [{"path": str(path), "sha256": _sha256_file(path), "size": path.stat().st_size}
                               for path in paths], "scientific_data_sha256": hashlib.sha256(
                      json.dumps(signatures[0], separators=(",", ":")).encode("ascii")).hexdigest(),
                  "scope": scope}
        verified["faster.application.yoyo"] = record

    if "faster.spinning_top_application_contract" in requests:
        paths = repeated_paths("faster.spinning_top_application_contract", "Faster spinning-top application contract")
        python_root = str(Path(__file__).resolve().parents[3] / "python")
        sys.path.insert(0, python_root)
        try:
            from make_it_stand.verification.run_evidence import PaperRunEvidenceV2, check_freshness
        finally:
            sys.path.pop(0)
        scope = "SYNTHETIC_SPINNING_TOP_GRAVITY_AXIS_FREE_VERTICAL_CENTER_AND_EQ7_COEFFICIENT_CONDITIONS_ONLY_NOT_AUTHOR_CASE"
        expected_markers = (
            "[PASS] spinning-top gravity-axis objective, free vertical center, and inertia-block constraints agree independently",
            "[PASS] spinning-top Eq. 7 coefficient condition A33<0 and p!=0 agrees independently",
            "[PASS] spinning-top application scope is synthetic method validation only",
        )
        signatures = []
        for path in paths:
            evidence = PaperRunEvidenceV2.from_json(path.read_text(encoding="utf-8"))
            if check_freshness(evidence).status != "CURRENT":
                raise ValueError("STALE_RESULT: Faster spinning-top application evidence is not current")
            root = path.parent; scientific, summary = _load_json(root / "scientific.json"), _load_json(root / "summary.json")
            outputs = {Path(row["path"]).resolve() for row in evidence.as_dict()["scientific"]["outputs"]}
            command, transcript = scientific.get("command"), root / "spinning-top-contract.stdout"
            if (scientific.get("schema") != "faster-spinning-top-contract-scientific-v1" or scientific.get("classification") != "SYNTHETIC_REGRESSION" or scientific.get("scope") != scope or summary.get("status") != "SYNTHETIC_SPINNING_TOP_APPLICATION_CONTRACT_VALIDATED" or summary.get("scope") != scope or not isinstance(command, dict) or command.get("name") != "spinning-top-contract" or Path(command.get("path", "")).resolve() != transcript.resolve() or transcript.resolve() not in outputs or not transcript.is_file() or command.get("returncode") != 0 or command.get("required_markers") != list(expected_markers) or command.get("sha256") != _sha256_file(transcript) or Path(command.get("binary", "")).name != "spin-moment-problem-regression" or not Path(command.get("binary", "")).is_file() or command.get("binary_sha256") != _sha256_file(Path(command["binary"]))):
                raise ValueError("Faster spinning-top application evidence is malformed")
            if any(marker not in transcript.read_text(encoding="utf-8") for marker in expected_markers):
                raise ValueError("Faster spinning-top application transcript lacks a required marker")
            signatures.append((command["sha256"], command["binary_sha256"]))
        if signatures[0] != signatures[1]: raise ValueError("repeated Faster spinning-top application scientific outputs differ")
        record = {"checker": {"path": str(Path(__file__).resolve()), "sha256": _sha256_file(Path(__file__)), "size": Path(__file__).stat().st_size}, "evidence": [{"path": str(path), "sha256": _sha256_file(path), "size": path.stat().st_size} for path in paths], "scientific_data_sha256": hashlib.sha256(json.dumps(signatures[0], separators=(",", ":")).encode("ascii")).hexdigest(), "scope": scope}
        verified["faster.application.spinning_top"] = record

    if "faster.cad_methods" in requests:
        paths = repeated_paths("faster.cad_methods", "Faster pinned CAD methods")
        python_root = str(Path(__file__).resolve().parents[3] / "python")
        sys.path.insert(0, python_root)
        try:
            from make_it_stand.verification.run_evidence import PaperRunEvidenceV2, check_freshness
        finally:
            sys.path.pop(0)
        scope = "SYNTHETIC_PINNED_OCCT_7_9_3_BREP_STEP_ROUNDTRIP_AND_ANALYTIC_MASS_METHOD_ONLY_NOT_AUTHOR_CAD_CASE"
        expected_volume = 8.0 - math.pi / 6.0
        signatures = []
        for path in paths:
            evidence = PaperRunEvidenceV2.from_json(path.read_text(encoding="utf-8"))
            if check_freshness(evidence).status != "CURRENT":
                raise ValueError("STALE_RESULT: Faster pinned CAD-method evidence is not current")
            root = path.parent
            scientific, summary = _load_json(root / "scientific.json"), _load_json(root / "summary.json")
            outputs = {Path(row["path"]).resolve() for row in evidence.as_dict()["scientific"]["outputs"]}
            if scientific.get("schema") != "faster-cad-methods-scientific-v1" or \
                    scientific.get("classification") != "SYNTHETIC_REGRESSION" or \
                    scientific.get("scope") != scope or \
                    summary.get("status") != "SYNTHETIC_PINNED_OCCT_CAD_METHOD_VALIDATED":
                raise ValueError("Faster pinned CAD-method evidence has invalid scope")
            probe, artifacts, result = scientific.get("probe"), scientific.get("artifacts"), scientific.get("result")
            if not isinstance(probe, dict) or probe.get("returncode") != 0 or \
                    not isinstance(artifacts, dict) or not isinstance(result, dict) or \
                    result.get("classification") != "SYNTHETIC_REGRESSION" or \
                    result.get("status") != "SYNTHETIC_OCCT_STEP_BREP_ROUNDTRIP_VALIDATED" or \
                    not math.isclose(float(result.get("analytic_volume", math.nan)), expected_volume, abs_tol=1e-14):
                raise ValueError("Faster pinned CAD-method result is malformed")
            transcript = Path(probe.get("transcript", {}).get("path", ""))
            binary = Path(probe.get("path", ""))
            if not transcript.is_file() or transcript.resolve() not in outputs or \
                    _sha256_file(transcript) != probe.get("transcript", {}).get("sha256") or \
                    not binary.is_file() or _sha256_file(binary) != probe.get("sha256"):
                raise ValueError("Faster pinned CAD-method probe provenance is malformed")
            for name, filename in (("brep", "material.brep"), ("step", "step.step")):
                row = artifacts.get(name)
                artifact = Path(row.get("path", "")) if isinstance(row, dict) else Path()
                if artifact.name != filename or not artifact.is_file() or artifact.stat().st_size == 0 or \
                        artifact.resolve() not in outputs or _sha256_file(artifact) != row.get("sha256"):
                    raise ValueError("Faster pinned CAD-method artifact is malformed")
            moments = result.get("moments")
            if not isinstance(moments, dict) or set(moments) != {"original", "brep", "step"}:
                raise ValueError("Faster pinned CAD-method misses a material moment")
            normalized = []
            for name in ("original", "brep", "step"):
                row = moments[name]
                if not isinstance(row, dict) or row.get("center") != [1, 1, 1] or \
                        not math.isclose(float(row.get("volume", math.nan)), expected_volume, abs_tol=1e-10):
                    raise ValueError("Faster pinned CAD-method mass moment is malformed")
                normalized.append((name, row["volume"], tuple(row["center"])))
            signatures.append((probe["sha256"], result["analytic_volume"], tuple(normalized)))
        if signatures[0] != signatures[1]:
            raise ValueError("repeated Faster pinned CAD-method scientific outputs differ")
        record = {"checker": {"path": str(Path(__file__).resolve()), "sha256": _sha256_file(Path(__file__)),
                              "size": Path(__file__).stat().st_size},
                  "evidence": [{"path": str(path), "sha256": _sha256_file(path), "size": path.stat().st_size}
                               for path in paths], "scientific_data_sha256": hashlib.sha256(
                      json.dumps(signatures[0], separators=(",", ":")).encode("ascii")).hexdigest(),
                  "scope": scope}
        verified["faster.cad.methods"] = record

    if "faster.eq14_newton" in requests:
        paths = repeated_paths("faster.eq14_newton", "Faster Equation-14 Newton")
        python_root = str(Path(__file__).resolve().parents[3] / "python")
        sys.path.insert(0, python_root)
        try:
            from make_it_stand.verification.run_evidence import PaperRunEvidenceV2, check_freshness
        finally:
            sys.path.pop(0)
        expected = {
            "default": ((), ("[PASS] scenario=quadric-hollowing-solver iterations=", " residual=", " trace=")),
            "raw-stopping": (("--raw-stopping",), ()),
            "dual-pairing": (("--dual-pairing",), ()),
            "dual-coordinates": (("--dual-coordinates",), ()),
        }
        signatures = []
        for path in paths:
            evidence = PaperRunEvidenceV2.from_json(path.read_text(encoding="utf-8"))
            if check_freshness(evidence).status != "CURRENT":
                raise ValueError("STALE_RESULT: Faster Equation-14 Newton evidence is not current")
            root = path.parent
            scientific, summary = _load_json(root / "scientific.json"), _load_json(root / "summary.json")
            outputs = {Path(row["path"]).resolve() for row in evidence.as_dict()["scientific"]["outputs"]}
            if scientific.get("schema") != "faster-eq14-newton-scientific-v1" or \
                    scientific.get("classification") != "SYNTHETIC_REGRESSION" or \
                    summary.get("status") != "SYNTHETIC_EQ14_NEWTON_METHOD_VALIDATED":
                raise ValueError("Faster Equation-14 Newton evidence is malformed")
            commands = scientific.get("commands")
            if not isinstance(commands, list) or {row.get("name") for row in commands} != set(expected):
                raise ValueError("Faster Equation-14 Newton evidence misses a protected mode")
            normalized = []
            for row in commands:
                name = row["name"]
                arguments, markers = expected[name]
                transcript = Path(row.get("path", ""))
                if row.get("arguments") != list(arguments) or row.get("returncode") != 0 or \
                        row.get("required_markers") != list(markers) or transcript.resolve() not in outputs or \
                        not transcript.is_file() or _sha256_file(transcript) != row.get("sha256"):
                    raise ValueError("Faster Equation-14 Newton native command is malformed")
                content = transcript.read_text(encoding="utf-8")
                if any(marker not in content for marker in markers):
                    raise ValueError("Faster Equation-14 Newton transcript lacks its required marker")
                normalized.append((name, row["sha256"]))
            signatures.append(tuple(sorted(normalized)))
        if signatures[0] != signatures[1]:
            raise ValueError("repeated Faster Equation-14 Newton outputs differ")
        record = {
            "checker": {"path": str(Path(__file__).resolve()), "sha256": _sha256_file(Path(__file__)),
                        "size": Path(__file__).stat().st_size},
            "evidence": [{"path": str(path), "sha256": _sha256_file(path), "size": path.stat().st_size}
                         for path in paths],
            "scientific_data_sha256": hashlib.sha256(
                json.dumps(signatures[0], separators=(",", ":"), ensure_ascii=True).encode("ascii")
            ).hexdigest(),
            "scope": "SYNTHETIC_NATIVE_EQ14_NEWTON_GEODESIC_MULTIPLIER_AND_RESTART_METHODS_ONLY_NOT_AUTHOR_CASE_OR_GLOBAL_OPTIMALITY",
        }
        verified.update({name: dict(record) for name in
                         ("faster.eq14.newton", "faster.eq14.geodesic", "faster.eq14.multiplier",
                          "faster.eq14.restart")})

    if "faster.stratum_matrix" in requests:
        paths = repeated_paths("faster.stratum_matrix", "Faster stratum matrix")
        python_root = str(Path(__file__).resolve().parents[3] / "python")
        sys.path.insert(0, python_root)
        try:
            from make_it_stand.verification.run_evidence import PaperRunEvidenceV2, check_freshness
        finally:
            sys.path.pop(0)
        signatures=[]
        expected={"plane","elliptic-cylinder","elliptic-paraboloid","hyperbolic-paraboloid","parabolic-cylinder","cone","parallel-planes","intersecting-planes","double-plane","isolated-point","imaginary-ellipsoid","ellipsoid-control","one-sheet-hyperboloid-control","two-sheet-hyperboloid-control"}
        for path in paths:
            evidence=PaperRunEvidenceV2.from_json(path.read_text(encoding="utf-8"))
            if check_freshness(evidence).status!="CURRENT": raise ValueError("STALE_RESULT: Faster stratum evidence is not current")
            root=path.parent; scientific=_load_json(root/"scientific.json"); summary=_load_json(root/"summary.json"); matrix=_load_json(root/"stratum.json")
            outputs={Path(x["path"]).resolve() for x in evidence.as_dict()["scientific"]["outputs"]}
            if scientific.get("schema")!="faster-stratum-matrix-scientific-v1" or summary.get("status")!="SYNTHETIC_STRATUM_MATRIX_VALIDATED" or matrix.get("schema")!="m1-stratum-matrix-result-v1" or not matrix.get("overall_pass") or {x.get("name") for x in matrix.get("cases",[])}!=expected or any(x.get("violations") for x in matrix["cases"]) or (root/"stratum.json").resolve() not in outputs or _sha256_file(root/"stratum.json")!=scientific.get("matrix_sha256"):
                raise ValueError("Faster stratum evidence is malformed")
            signatures.append(scientific["matrix_sha256"])
        if signatures[0]!=signatures[1]: raise ValueError("repeated Faster stratum outputs differ")
        record={"checker":{"path":str(Path(__file__).resolve()),"sha256":_sha256_file(Path(__file__)),"size":Path(__file__).stat().st_size},"evidence":[{"path":str(p),"sha256":_sha256_file(p),"size":p.stat().st_size} for p in paths],"scientific_data_sha256":signatures[0],"scope":"SYNTHETIC_FROZEN_STRATUM_MATRIX_IDEAL_POINT_CAP_AND_DEGENERATE_TYPE_METHODS_ONLY"}
        verified.update({x:dict(record) for x in ("faster.geometry.ideal_points","faster.geometry.caps","faster.geometry.degenerate_types")})

    if "spin2014.source_cage" in requests:
        paths = repeated_paths("spin2014.source_cage", "Spin-It 2014 source cage")
        python_root = str(Path(__file__).resolve().parents[3] / "python")
        sys.path.insert(0, python_root)
        try:
            from make_it_stand.verification.run_evidence import PaperRunEvidenceV2, check_freshness
        finally:
            sys.path.pop(0)
        observations = []
        for path in paths:
            evidence = PaperRunEvidenceV2.from_json(path.read_text(encoding="utf-8"))
            modules = evidence.as_dict()["begin"]["observation"]["runtime"]["native_modules"]
            module = modules.get("_make_it_stand") if isinstance(modules, dict) else None
            if not isinstance(module, dict):
                raise ValueError("source cage evidence lacks native module identity")
            module_path = Path(module.get("path", "")).resolve(strict=True)
            if module.get("sha256") != _sha256_file(module_path) or module.get("size") != module_path.stat().st_size:
                raise ValueError("STALE_RESULT: source cage native module differs from observation")
            try:
                _, freshness = _load_bound_native_module_and_freshness(
                    evidence, module_path, check_freshness)
            except (ImportError, OSError, ValueError) as error:
                raise ValueError("source cage cannot load the evidenced native module: " + str(error)) from error
            if freshness.status != "CURRENT":
                raise ValueError("STALE_RESULT: source cage evidence is not current: " +
                                 "; ".join(freshness.mismatches))
            root = path.parent
            scientific, summary = _load_json(root / "scientific.json"), _load_json(root / "summary.json")
            outputs = {Path(row["path"]).resolve() for row in evidence.as_dict()["scientific"]["outputs"]}
            required = {root / "scientific.json", root / "surface.stdout", root / "bbw.stdout"}
            if not {item.resolve() for item in required} <= outputs:
                raise ValueError("source cage evidence does not bind all scientific outputs")
            if summary.get("status") != "SYNTHETIC_METHOD_VALIDATED" or \
                    summary.get("classification") != "SYNTHETIC_REGRESSION" or \
                    summary.get("original_evidence_admission") != "BLOCKED_UNKNOWN":
                raise ValueError("source cage summary has an invalid scope")
            if scientific.get("schema") != "source-spin2014-cage-scientific-v1" or \
                    scientific.get("classification") != "SYNTHETIC_REGRESSION":
                raise ValueError("source cage scientific record is malformed")
            records = []
            for field, markers in (("surface", ("level9 GEOMETRY:", "complete geometry only")),
                                   ("bbw", ("actual_702_variable_and_26_column_box_qps", "SYNTHETIC_REGRESSION: 0 failed groups"))):
                row = scientific.get(field)
                if not isinstance(row, dict) or row.get("returncode") != 0 or row.get("required_markers") != list(markers):
                    raise ValueError("source cage native result is malformed")
                transcript = Path(row.get("path", ""))
                if transcript.resolve() not in outputs or not transcript.is_file():
                    raise ValueError("source cage transcript is not bound")
                data = transcript.read_bytes()
                if hashlib.sha256(data).hexdigest() != row.get("sha256") or any(token not in data.decode("utf-8") for token in markers):
                    raise ValueError("source cage transcript does not establish its declared native result")
                records.append((row["sha256"], tuple(markers)))
            observations.append((path, tuple(records)))
        if observations[0][0].parent.samefile(observations[1][0].parent) or observations[0][1] != observations[1][1]:
            raise ValueError("repeated source cage scientific outputs differ")
        record = {"checker": {"path": str(Path(__file__).resolve()), "sha256": _sha256_file(Path(__file__)),
                              "size": Path(__file__).stat().st_size},
                  "evidence": [{"path": str(path), "sha256": _sha256_file(path), "size": path.stat().st_size}
                               for path in paths],
                  "scientific_data_sha256": observations[0][1][0][0] + observations[0][1][1][0],
                  "scope": "SYNTHETIC_COMPLETE_SOURCE_CAGE_AND_ACTUAL_BBW_METHODS_ONLY_NOT_AUTHOR_EXPERIMENT"}
        verified.update({name: dict(record) for name in
                         ("spin2014.cage.complete_octree", "spin2014.cage.bbw_interpolation")})

    if "spin2014.source_deformation" in requests:
        paths = repeated_paths("spin2014.source_deformation", "Spin-It 2014 joint deformation")
        python_root = str(Path(__file__).resolve().parents[3] / "python"); sys.path.insert(0, python_root)
        try:
            from make_it_stand.verification.run_evidence import PaperRunEvidenceV2, check_freshness
        finally: sys.path.pop(0)
        signatures=[]
        expected={"functional":"[PASS] scenario=spin-it-2014-cage-functional","solver":"[PASS] scenario=spin-it-2014-cage-solver","shape":"SYNTHETIC_REGRESSION source cage shape PASS","geometry":"actual-source geometry admission PASS","guarded":"actual-source guarded SLQP PASS"}
        for path in paths:
            evidence=PaperRunEvidenceV2.from_json(path.read_text(encoding="utf-8")); modules=evidence.as_dict()["begin"]["observation"]["runtime"]["native_modules"]; module=modules.get("_make_it_stand") if isinstance(modules,dict) else None
            if not isinstance(module,dict): raise ValueError("joint deformation lacks native module identity")
            mp=Path(module.get("path","")).resolve(strict=True)
            if module.get("sha256")!=_sha256_file(mp) or module.get("size")!=mp.stat().st_size: raise ValueError("STALE_RESULT: joint deformation native module changed")
            if "_make_it_stand" not in sys.modules:
                sys.path.insert(0,str(mp.parent))
                try: importlib.import_module("_make_it_stand")
                finally: sys.path.pop(0)
            if check_freshness(evidence).status!="CURRENT": raise ValueError("STALE_RESULT: joint deformation evidence is not current")
            root=path.parent; science=_load_json(root/"scientific.json"); summary=_load_json(root/"summary.json")
            if science.get("schema")!="source-spin2014-deformation-scientific-v1" or summary.get("status")!="SYNTHETIC_JOINT_DEFORMATION_VALIDATED": raise ValueError("joint deformation record malformed")
            outputs={Path(x["path"]).resolve() for x in evidence.as_dict()["scientific"]["outputs"]}; rows=science.get("checks")
            if not isinstance(rows,list) or len(rows)!=5: raise ValueError("joint deformation omitted a required check")
            sig=[]
            for row in rows:
                label=row.get("label"); marker=expected.get(label); transcript=Path(row.get("path","")).resolve()
                if marker is None or transcript not in outputs or row.get("returncode")!=0 or marker not in transcript.read_text(encoding="utf-8") or hashlib.sha256(transcript.read_bytes()).hexdigest()!=row.get("sha256"): raise ValueError("joint deformation check invalid")
                sig.append((label,row["sha256"]))
            signatures.append(tuple(sig))
        if signatures[0]!=signatures[1]: raise ValueError("repeated joint deformation outputs differ")
        record={"checker":{"path":str(Path(__file__).resolve()),"sha256":_sha256_file(Path(__file__)),"size":Path(__file__).stat().st_size},"evidence":[{"path":str(x),"sha256":_sha256_file(x),"size":x.stat().st_size} for x in paths],"scientific_data_sha256":"".join(x[1] for x in signatures[0]),"scope":"SYNTHETIC_JOINT_FILL_CAGE_DEFORMATION_AND_GEOMETRY_GUARDS_ONLY_NOT_AUTHOR_EXPERIMENT"}
        verified.update({name: dict(record) for name in
                         ("spin2014.deformation.regularizers",
                          "spin2014.deformation.joint")})

    if "spin2014.dual_material" in requests:
        paths = repeated_paths("spin2014.dual_material", "Spin-It 2014 dual material")
        python_root = str(Path(__file__).resolve().parents[3] / "python")
        sys.path.insert(0, python_root)
        try:
            from make_it_stand.verification.run_evidence import PaperRunEvidenceV2, check_freshness
        finally:
            sys.path.pop(0)
        signatures = []
        markers = (
            "[PASS] scenario=spin-it-2014-dual-material mass=",
            "paper_givens=validated",
            "independent_label_moments=validated",
        )
        for path in paths:
            evidence = PaperRunEvidenceV2.from_json(path.read_text(encoding="utf-8"))
            modules = evidence.as_dict()["begin"]["observation"]["runtime"]["native_modules"]
            module = modules.get("_make_it_stand") if isinstance(modules, dict) else None
            if not isinstance(module, dict):
                raise ValueError("dual-material evidence lacks native module identity")
            module_path = Path(module.get("path", "")).resolve(strict=True)
            if module.get("sha256") != _sha256_file(module_path) or module.get("size") != module_path.stat().st_size:
                raise ValueError("STALE_RESULT: dual-material native module changed")
            if "_make_it_stand" not in sys.modules:
                sys.path.insert(0, str(module_path.parent))
                try:
                    importlib.import_module("_make_it_stand")
                finally:
                    sys.path.pop(0)
            if check_freshness(evidence).status != "CURRENT":
                raise ValueError("STALE_RESULT: dual-material evidence is not current")
            root = path.parent
            scientific = _load_json(root / "scientific.json")
            summary = _load_json(root / "summary.json")
            if scientific.get("schema") != "spin2014-dual-material-scientific-v1" or \
                    scientific.get("classification") != "SYNTHETIC_REGRESSION" or \
                    summary.get("status") != "SYNTHETIC_DUAL_MATERIAL_VALIDATED" or \
                    summary.get("original_evidence_admission") != "BLOCKED_UNKNOWN" or \
                    summary.get("paper_case_admission") != "BLOCKED_MISSING_AUTHOR_SUMONE_ASSET_AND_PARAMETERS":
                raise ValueError("dual-material evidence has invalid scope")
            outputs = {Path(row["path"]).resolve() for row in evidence.as_dict()["scientific"]["outputs"]}
            transcript = root / "dual-material.stdout"
            required_outputs = {root / "scientific.json", transcript}
            if not {item.resolve() for item in required_outputs} <= outputs:
                raise ValueError("dual-material evidence does not bind all scientific outputs")
            result = scientific.get("native_regression")
            if not isinstance(result, dict) or result.get("returncode") != 0 or \
                    result.get("required_markers") != list(markers):
                raise ValueError("dual-material native result is malformed")
            recorded = Path(result.get("path", "")).resolve()
            if recorded != transcript.resolve() or recorded not in outputs or not recorded.is_file():
                raise ValueError("dual-material transcript is not bound")
            data = recorded.read_bytes()
            if hashlib.sha256(data).hexdigest() != result.get("sha256") or \
                    any(marker not in data.decode("utf-8") for marker in markers):
                raise ValueError("dual-material transcript does not establish the native result")
            signatures.append(result["sha256"])
        if signatures[0] != signatures[1]:
            raise ValueError("repeated dual-material scientific outputs differ")
        record = {
            "checker": {"path": str(Path(__file__).resolve()), "sha256": _sha256_file(Path(__file__)),
                        "size": Path(__file__).stat().st_size},
            "evidence": [{"path": str(path), "sha256": _sha256_file(path), "size": path.stat().st_size}
                         for path in paths],
            "scientific_data_sha256": signatures[0],
            "scope": "SYNTHETIC_DUAL_DENSITY_MATERIAL_LABEL_AND_MOMENT_SUM_METHOD_ONLY_NOT_SUMONE_OR_AUTHOR_EXPERIMENT",
        }
        verified["spin2014.material.dual_density"] = record

    if "spin2014.source_method" in requests:
        paths = repeated_paths("spin2014.source_method", "Spin-It 2014 source method")
        python_root = str(Path(__file__).resolve().parents[3] / "python")
        sys.path.insert(0, python_root)
        try:
            from make_it_stand.verification.run_evidence import (
                PaperRunEvidenceV2, check_freshness,
            )
        finally:
            sys.path.pop(0)

        observations = []
        bound_outputs: dict[Path, set[Path]] = {}
        for path in paths:
            observed = PaperRunEvidenceV2.from_json(path.read_text(encoding="utf-8"))
            runtime_modules = observed.as_dict()["begin"]["observation"]["runtime"]["native_modules"]
            module_record = runtime_modules.get("_make_it_stand") if isinstance(runtime_modules, dict) else None
            if not isinstance(module_record, dict):
                raise ValueError("source method evidence lacks the native module record")
            try:
                module_path = Path(module_record["path"]).resolve(strict=True)
            except (KeyError, OSError, TypeError) as error:
                raise ValueError("source method native module path is unavailable") from error
            if module_record.get("size") != module_path.stat().st_size or \
                    module_record.get("sha256") != _sha256_file(module_path):
                raise ValueError("STALE_RESULT: source method native module differs from its observation")
            try:
                _, freshness = _load_bound_native_module_and_freshness(
                    observed, module_path, check_freshness)
            except (ImportError, OSError, ValueError) as error:
                raise ValueError("source method cannot load the evidenced native module") from error
            if freshness.status != "CURRENT":
                raise ValueError("STALE_RESULT: " + "; ".join(freshness.mismatches))
            run_root = path.parent
            native_path, scientific_path, summary_path = (
                run_root / "native-run.json", run_root / "scientific.json",
                run_root / "summary.json",
            )
            try:
                native_run = _load_json(native_path)
                scientific = _load_json(scientific_path)
                summary = _load_json(summary_path)
            except (OSError, json.JSONDecodeError) as error:
                raise ValueError("source method outputs are unavailable") from error
            _reject_nonfinite(native_run, "source native run")
            _reject_nonfinite(scientific, "source scientific data")
            _reject_nonfinite(summary, "source summary")
            output_paths = [Path(row["path"]).resolve() for row in observed.as_dict()["scientific"]["outputs"]]
            if native_path.resolve() not in output_paths or scientific_path.resolve() not in output_paths:
                raise ValueError("source method evidence does not bind native/scientific output")
            bound_outputs[path] = set(output_paths)
            observations.append((path, native_run, scientific, summary))

        first, second = observations
        if first[0].parent.samefile(second[0].parent):
            raise ValueError("two source-method evidence files refer to one run directory")
        if json.dumps(first[2], sort_keys=True, separators=(",", ":"), allow_nan=False) != \
                json.dumps(second[2], sort_keys=True, separators=(",", ":"), allow_nan=False):
            raise ValueError("repeated source-method scientific outputs differ")
        native_run, scientific, summary = first[1], first[2], first[3]
        if not isinstance(native_run, dict) or not isinstance(scientific, dict) or not isinstance(summary, dict):
            raise ValueError("source method outputs are malformed")
        if scientific.get("classification") != "SYNTHETIC_REGRESSION" or \
                summary.get("classification") != "SYNTHETIC_REGRESSION" or \
                summary.get("status") != "BINARY_MOMENT_VALIDATED" or not summary.get("sealed"):
            raise ValueError("source method does not establish the bounded synthetic outcome")
        if summary.get("build_source_attestation") != "NOT_VERIFIED" or \
                summary.get("original_evidence_admission") != "BLOCKED_UNKNOWN":
            raise ValueError("source method scope is not explicitly bounded")
        solve = native_run.get("solve")
        settings = native_run.get("settings", {}).get("solver_settings")
        if not isinstance(solve, dict) or not isinstance(settings, dict) or \
                solve.get("status") != "BINARY_MOMENT_VALIDATED":
            raise ValueError("source method solve record is malformed")
        epsilon = settings.get("fractional_epsilon")
        tolerance = settings.get("binary_constraint_tolerance")
        initial_fill = settings.get("initial_fill")
        if not isinstance(epsilon, (int, float)) or not 0 < epsilon < .5 or \
                not isinstance(tolerance, (int, float)) or tolerance <= 0 or \
                not isinstance(initial_fill, (int, float)) or not 0 <= initial_fill <= 1:
            raise ValueError("source method fill controls are malformed")
        passes = solve.get("passes")
        relaxed, binary, rounding = (solve.get(key) for key in
                                     ("relaxed_fill", "binary_fill", "rounding_delta"))
        if not isinstance(passes, list) or not passes or not all(
                isinstance(value, (int, float)) and 0 <= value <= 1 for value in relaxed or ()) or \
                not relaxed or not isinstance(binary, list) or len(binary) != len(relaxed) or \
                not isinstance(rounding, list) or len(rounding) != len(relaxed):
            raise ValueError("source method fill vectors are malformed")
        if not any(0 < value < 1 for item in passes for value in item.get("input_beta", ())):
            raise ValueError("source method trace does not contain a continuous fill witness")
        for index, (continuous, discrete, delta) in enumerate(zip(relaxed, binary, rounding)):
            expected = 1.0 if continuous >= .5 else 0.0
            if discrete not in (0.0, 1.0) or discrete != expected or abs(delta - (discrete - continuous)) > 1e-14:
                raise ValueError("source method binary rounding rule differs at index " + str(index))
        binary_check = solve.get("binary_check")
        if not isinstance(binary_check, dict) or binary_check.get("passed") is not True or \
                binary_check.get("normalized_moment_error", float("inf")) > tolerance or \
                binary_check.get("normalized_constraint_residual", float("inf")) > tolerance:
            raise ValueError("source method independent binary check is malformed or exceeds tolerance")

        # The native numerical trace above says nothing about whether the
        # material object really preserves a fixed shell.  Reload each sealed
        # object through the independent exact-geometry checker, then inspect
        # its recomputed facts.  This grants only a synthetic *method* witness:
        # the bundle itself remains explicitly non-author and its unresolved
        # same-component thickness stays UNKNOWN.
        python_root = str(Path(__file__).resolve().parents[3] / "python")
        sys.path.insert(0, python_root)
        try:
            from make_it_stand.verification.spin_bundle import load_spin_artifact
        finally:
            sys.path.pop(0)
        seal_digests = []
        for evidence_path, _, repeated_scientific, repeated_summary in observations:
            root = evidence_path.parent
            manifest = root / "sealed" / "manifest.json"
            if manifest.resolve() not in bound_outputs[evidence_path]:
                raise ValueError("source method evidence does not bind sealed material manifest")
            identity = repeated_scientific.get("identity") if isinstance(repeated_scientific, dict) else None
            input_identity = identity.get("input_identity") if isinstance(identity, dict) else None
            expected_digest = identity.get("problem_digest") if isinstance(identity, dict) else None
            if not isinstance(expected_digest, str):
                raise ValueError("source method sealed material lacks locked problem digest")
            try:
                artifact = load_spin_artifact(root / "sealed", expected_problem_digest=expected_digest)
            except (OSError, ValueError, RuntimeError) as error:
                raise ValueError("source method sealed material cannot be independently reloaded") from error
            if artifact.scientific_digest != repeated_summary.get("artifact_scientific_digest") or \
                    artifact.scientific_digest != identity.get("artifact_scientific_digest"):
                raise ValueError("source method sealed material digest differs from observed solve")
            facts = artifact.facts
            if facts.get("geometry_status") != "MESH_MATERIAL_VALIDATED" or \
                    facts.get("original_evidence_admission") != "BLOCKED_UNKNOWN" or \
                    facts.get("manufacturing_admission") != "BLOCKED_UNKNOWN_WALL_THICKNESS":
                raise ValueError("source method sealed material has an invalid scope")
            if facts.get("outer") != input_identity.get("original") or artifact.outer != artifact.problem.original:
                raise ValueError("source method sealed material changed the outer surface")
            for partition_name in ("material_partition", "carved_partition"):
                partition = facts.get(partition_name)
                if not isinstance(partition, dict) or not all(partition.get(key) is value for key, value in (
                        ("completed", True), ("complete_partition", True),
                        ("union_equals_whole", True), ("no_positive_volume_overlap", True),
                        ("excess_union_content", False), ("missing_whole_content", False))):
                    raise ValueError("source method sealed material has an invalid " + partition_name)
            decompositions = facts.get("decompositions")
            if not isinstance(decompositions, list) or len(decompositions) != 1 or \
                    decompositions[0].get("id") != "design":
                raise ValueError("source method lacks an explicit design partition")
            design_partition = decompositions[0].get("partition")
            if not isinstance(design_partition, dict) or not all(design_partition.get(key) is value for key, value in (
                    ("completed", True), ("complete_partition", True),
                    ("union_equals_whole", True), ("no_positive_volume_overlap", True),
                    ("excess_union_content", False), ("missing_whole_content", False))):
                raise ValueError("source method sealed material has an invalid design partition")
            clearance = facts.get("clearance")
            if not isinstance(clearance, dict) or not all(clearance.get(key) is value for key, value in (
                    ("completed", True), ("applicable", True), ("at_least_minimum", True), ("contact", False))):
                raise ValueError("source method sealed material lacks verified outer-wall clearance")
            wall = facts.get("wall_thickness_audit")
            if wall != {"outer_to_newly_carved": "VERIFIED", "same_component_feature": "UNKNOWN",
                        "inter_cavity_web": "NOT_APPLICABLE", "overall": "UNKNOWN"}:
                raise ValueError("source method sealed material changes unresolved wall-thickness scope")
            regions = facts.get("regions")
            if not isinstance(regions, list) or [(row.get("id"), row.get("kind")) for row in regions] != \
                    [("shell", "fixed"), ("design", "retained")]:
                raise ValueError("source method sealed material does not preserve fixed and retained regions")
            comparison = facts.get("moment_comparison")
            if not isinstance(comparison, dict) or comparison.get("status") != "WITHIN_ENGINEERING_THRESHOLD" or \
                    not isinstance(comparison.get("absolute_errors"), list) or \
                    not all(error <= 1e-6 for error in comparison["absolute_errors"]):
                raise ValueError("source method sealed material has inconsistent independent mass moments")
            seal_digests.append(artifact.scientific_digest)
        if len(set(seal_digests)) != 1:
            raise ValueError("repeated source-method sealed material outputs differ")

        # This is deliberately a literal-loop audit of the *observed* native
        # trace, not a trust of its `final_kkt` or active-set labels.  The small
        # source fixture is sufficient to establish the three solver methods,
        # but remains insufficient for any paper case or geometric claim.
        slqp_controls = settings.get("slqp")
        if not isinstance(slqp_controls, dict):
            raise ValueError("source method SLQP controls are malformed")
        active_tolerance = slqp_controls.get("active_set_tolerance")
        maximum_iterations = slqp_controls.get("maximum_iterations")
        if slqp_controls.get("algorithm") != "AUDITED_KKT_RESTORATION_V1" or \
                not isinstance(active_tolerance, (int, float)) or not 0 <= active_tolerance < .5 or \
                not isinstance(maximum_iterations, int) or maximum_iterations < 1:
            raise ValueError("source method SLQP controls are malformed")

        def literal_kkt(evaluation: Any, state: Any, witness: Any, beta_count: Any) -> tuple[float, float, float, float, float]:
            if not isinstance(evaluation, dict) or not isinstance(witness, dict) or \
                    not isinstance(state, list) or not isinstance(beta_count, int) or \
                    beta_count < 0 or beta_count > len(state):
                raise ValueError("source method KKT record is malformed")
            gradient = evaluation.get("gradient")
            constraints = evaluation.get("constraints")
            jacobian = evaluation.get("constraint_jacobian")
            lambdas, alpha, eta = (witness.get(key) for key in ("lambda", "alpha", "eta"))
            dimension = len(state)
            if not all(isinstance(vector, list) for vector in (gradient, constraints, jacobian, lambdas, alpha, eta)) or \
                    len(gradient) != dimension or len(alpha) != dimension or len(eta) != dimension or \
                    len(lambdas) != len(constraints) or len(jacobian) != len(constraints) or \
                    any(not isinstance(row, list) or len(row) != dimension for row in jacobian):
                raise ValueError("source method KKT dimensions are malformed")
            equality = max((abs(value) for value in constraints), default=0.0)
            bound = max(0.0, *(max(-state[index], state[index] - 1.0) for index in range(beta_count)))
            stationarity = 0.0
            dual_sign = 0.0
            complementarity = 0.0
            for index in range(dimension):
                residual = gradient[index] - alpha[index] + eta[index] + sum(
                    jacobian[row][index] * lambdas[row] for row in range(len(constraints)))
                stationarity = max(stationarity, abs(residual))
                dual_sign = max(dual_sign, -alpha[index], -eta[index])
                if index < beta_count:
                    complementarity = max(complementarity, abs(alpha[index] * state[index]),
                                          abs(eta[index] * (1.0 - state[index])))
                elif alpha[index] != 0.0 or eta[index] != 0.0:
                    raise ValueError("source method applies box duals outside its beta prefix")
            return equality, bound, stationarity, dual_sign, complementarity

        saw_active_transition = False
        saw_lbfgs_update = False
        saw_accepted_step = False
        for solver_pass in passes:
            if not isinstance(solver_pass, dict) or solver_pass.get("algorithm") != "AUDITED_KKT_RESTORATION_V1" or \
                    solver_pass.get("formulation") != "PAPER_GIVENS" or solver_pass.get("status") != "CONVERGED":
                raise ValueError("source method does not expose a converged audited paper SLQP pass")
            beta = solver_pass.get("beta")
            trace = solver_pass.get("trace")
            final_kkt = solver_pass.get("final_kkt")
            final_state = solver_pass.get("final_state")
            evaluation = solver_pass.get("evaluation")
            if not isinstance(beta, list) or not isinstance(trace, list) or not trace or \
                    len(beta) != len(final_state or ()) or len(beta) != solver_pass.get("bounded_variable_count"):
                raise ValueError("source method SLQP pass does not preserve its bounded state")
            final_values = literal_kkt(evaluation, final_state, final_kkt, len(beta))
            recorded = (final_kkt.get("equality_primal"), final_kkt.get("bound_primal"),
                        final_kkt.get("stationarity"), final_kkt.get("dual_sign"),
                        final_kkt.get("complementarity"))
            if final_kkt.get("valid") is not True or any(abs(left - right) > 1e-10
                    for left, right in zip(final_values, recorded)) or max(final_values) > 1e-7:
                raise ValueError("source method final KKT is not independently reproducible")
            previous_membership = None
            for row_index, row in enumerate(trace):
                if not isinstance(row, dict):
                    raise ValueError("source method trace row is malformed")
                state, membership, witness = (row.get(key) for key in ("state", "active_membership", "kkt"))
                if not isinstance(state, list) or not isinstance(membership, list) or len(state) != len(beta) or \
                        len(membership) != len(beta):
                    raise ValueError("source method active-set trace dimensions are malformed")
                digest = hashlib.sha256(b"spin2014-active-membership-v1" + bytes(membership)).hexdigest()
                if row.get("membership_hash") != digest:
                    raise ValueError("source method active-set membership digest is invalid")
                for index, value in enumerate(membership):
                    if value not in (0, 1, 2) or (value == 1 and state[index] > active_tolerance) or \
                            (value == 2 and state[index] < 1.0 - active_tolerance):
                        raise ValueError("source method active-set membership is inconsistent with its bounded state")
                row_evaluation = {
                    "gradient": row.get("gradient"),
                    "constraints": row.get("constraints"),
                    "constraint_jacobian": row.get("constraint_jacobian"),
                }
                values = literal_kkt(row_evaluation, state, witness, len(beta))
                recorded = (witness.get("equality_primal"), witness.get("bound_primal"),
                            witness.get("stationarity"), witness.get("dual_sign"), witness.get("complementarity"))
                if witness.get("valid") is not True or any(abs(left - right) > 1e-10
                        for left, right in zip(values, recorded)):
                    raise ValueError("source method trace KKT is not independently reproducible")
                if previous_membership is not None and membership != previous_membership:
                    saw_active_transition = True
                    if "ACTIVE_MEMBERSHIP_CHANGED" not in str(row.get("history_reset_reason")) or \
                            (row.get("curvature_pair_accepted") is not True and row.get("history_after") != 0):
                        raise ValueError("source method does not reset L-BFGS history after an active-set transition")
                previous_membership = membership
                if row.get("accepted_step", 0.0) > 0:
                    saw_accepted_step = True
                    trials = row.get("objective_trials")
                    if not isinstance(trials, list) or not any(
                            trial.get("accepted") is True and abs(trial.get("step") - row["accepted_step"]) <= 1e-14
                            for trial in trials if isinstance(trial, dict)):
                        raise ValueError("source method accepted an unrecorded SLQP line-search step")
                if row.get("curvature_pair_accepted") is True:
                    if row_index + 1 >= len(trace) or row.get("history_after", 0) <= 0:
                        raise ValueError("source method accepted an unverifiable L-BFGS curvature pair")
                    following = trace[row_index + 1]
                    step = [following["state"][i] - state[i] for i in range(len(state))]
                    gradient_change = [following["gradient"][i] - row["gradient"][i]
                                       for i in range(len(state))]
                    if not all(math.isfinite(value) for value in step + gradient_change) or \
                            sum(left * right for left, right in zip(step, gradient_change)) <= 0.0:
                        raise ValueError("source method L-BFGS curvature update lacks positive observed curvature")
                    saw_lbfgs_update = True
        if not saw_accepted_step or not saw_active_transition or not saw_lbfgs_update:
            raise ValueError("source method trace lacks a required SLQP, active-set, or L-BFGS witness")
        record = {
            "checker": {"path": str(Path(__file__).resolve()),
                        "sha256": _sha256_file(Path(__file__)),
                        "size": Path(__file__).stat().st_size},
            "evidence": [
                {"path": str(path), "sha256": _sha256_file(path), "size": path.stat().st_size}
                for path in paths
            ],
            "scientific_data_sha256": scientific.get("scientific_digest"),
            "scope": "SYNTHETIC_NATIVE_FILL_AND_SLQP_METHODS_ONLY_NOT_AUTHOR_EXPERIMENT_OR_BUILD_ATTESTATION",
        }
        verified.update({name: dict(record) for name in
                         ("spin2014.fill.continuous", "spin2014.fill.binary",
                          "spin2014.solver.slqp", "spin2014.solver.active_set",
                          "spin2014.solver.lbfgs")})
        fixed_shell_record = dict(record)
        fixed_shell_record["sealed_material_scientific_digest"] = seal_digests[0]
        fixed_shell_record["scope"] = (
            "SYNTHETIC_EXACT_FIXED_SHELL_AND_DESIGN_PARTITION_METHOD_ONLY_"
            "NOT_AUTHOR_EXPERIMENT_NOT_COMPLETE_FEATURE_THICKNESS")
        verified["spin2014.hollowing.fixed_shell"] = fixed_shell_record

    if "spin2014.source_adaptive_transitions" in requests:
        paths = repeated_paths("spin2014.source_adaptive_transitions", "Spin-It 2014 adaptive transitions")
        if "spin2014.source_method" not in requests:
            raise ValueError("adaptive split evidence requires independently repeated merge evidence")
        python_root = str(Path(__file__).resolve().parents[3] / "python")
        sys.path.insert(0, python_root)
        try:
            from make_it_stand.verification.run_evidence import PaperRunEvidenceV2, check_freshness
        finally:
            sys.path.pop(0)
        observations = []
        for path in paths:
            observed = PaperRunEvidenceV2.from_json(path.read_text(encoding="utf-8"))
            modules = observed.as_dict()["begin"]["observation"]["runtime"]["native_modules"]
            module_record = modules.get("_make_it_stand") if isinstance(modules, dict) else None
            if not isinstance(module_record, dict):
                raise ValueError("adaptive transition evidence lacks native module identity")
            module_path = Path(module_record["path"]).resolve(strict=True)
            if module_record.get("sha256") != _sha256_file(module_path) or module_record.get("size") != module_path.stat().st_size:
                raise ValueError("STALE_RESULT: adaptive transition module differs from observation")
            try:
                _, freshness = _load_bound_native_module_and_freshness(
                    observed, module_path, check_freshness)
            except (ImportError, OSError, ValueError) as error:
                raise ValueError("adaptive transition cannot load the evidenced native module") from error
            if freshness.status != "CURRENT":
                raise ValueError("STALE_RESULT: adaptive transition evidence is not current")
            root = path.parent
            native_run, scientific, inspection = (_load_json(root / name) for name in
                                                   ("native-run.json", "scientific.json", "domain-inspection.json"))
            outputs = {Path(row["path"]).resolve() for row in observed.as_dict()["scientific"]["outputs"]}
            if not {root / "native-run.json", root / "scientific.json", root / "domain-inspection.json"} <= outputs:
                raise ValueError("adaptive transition evidence does not bind all scientific outputs")
            _reject_nonfinite(native_run, "adaptive native run")
            _reject_nonfinite(scientific, "adaptive scientific output")
            _reject_nonfinite(inspection, "adaptive domain inspection")
            observations.append((path, native_run, scientific, inspection))
        first, second = observations
        if first[0].parent.samefile(second[0].parent) or \
                json.dumps(first[2], sort_keys=True, separators=(",", ":"), allow_nan=False) != json.dumps(second[2], sort_keys=True, separators=(",", ":"), allow_nan=False) or \
                json.dumps(first[3], sort_keys=True, separators=(",", ":"), allow_nan=False) != json.dumps(second[3], sort_keys=True, separators=(",", ":"), allow_nan=False):
            raise ValueError("adaptive transition repeats are not independent and deterministic")
        native_run, inspection = first[1], first[3]
        solve, source = native_run.get("solve"), native_run.get("source")
        if not isinstance(solve, dict) or not isinstance(source, dict) or \
                solve.get("status") != "BINARY_VALIDATION_FAILED" or solve.get("termination_policy") != "PAPER_POST_TRANSITION_V1":
            raise ValueError("adaptive transition outcome is not the registered expected binary failure")
        if inspection.get("schema") != "source-spin2014-domain-inspection-v1" or inspection.get("source", {}).get("source_hash") != source.get("source_hash") or \
                inspection.get("source", {}).get("partition_hash") != source.get("partition_hash"):
            raise ValueError("source-tree inspection is not bound to the solved source surface")
        cells = inspection.get("source_cells")
        stats = inspection.get("source", {}).get("stats", {})
        if not isinstance(cells, list) or not cells or cells[0].get("index") != 0 or cells[0].get("level") != 0:
            raise ValueError("source-tree root record is malformed")
        boundaries = [cell for cell in cells if isinstance(cell, dict) and cell.get("kind") == "BOUNDARY"]
        if len(boundaries) != stats.get("boundary_nodes") or not boundaries or any(
                cell.get("level") != stats.get("maximum_visited_level") or cell.get("unresolved_volume", 0) <= 0 for cell in boundaries):
            raise ValueError("source-tree boundary subtrees are not auditable at the configured finest level")
        for index, cell in enumerate(cells):
            if not isinstance(cell, dict) or cell.get("index") != index:
                raise ValueError("source-tree cell indexing is noncanonical")
            children = cell.get("children")
            if cell.get("kind") == "INTERNAL":
                if not isinstance(children, list) or len(children) != 8 or children != list(range(children[0], children[0] + 8)):
                    raise ValueError("source-tree internal child block is malformed")
                for child in children:
                    if child >= len(cells) or cells[child].get("parent") != index or cells[child].get("level") != cell.get("level") + 1:
                        raise ValueError("source-tree parent/child lineage is malformed")
        leaves, edges, triplets = inspection.get("leaves"), inspection.get("edges"), inspection.get("laplacian_triplets")
        maximum = inspection.get("maximum_optimization_level")
        if not isinstance(leaves, list) or not isinstance(edges, list) or not isinstance(triplets, list) or not isinstance(maximum, int):
            raise ValueError("source adaptive graph inspection is malformed")
        keys = [(item.get("level"), tuple(item.get("coordinates", ())) if isinstance(item, dict) else None) for item in leaves]
        if len(keys) != inspection.get("leaf_count") or len(set(keys)) != len(keys):
            raise ValueError("source adaptive leaves are noncanonical")
        def face_adjacent(left, right):
            la, lb = (1 << (maximum - left[0])), (1 << (maximum - right[0]))
            a0, b0 = ([value * la for value in left[1]], [value * lb for value in right[1]])
            for axis in range(3):
                if a0[axis] + la == b0[axis] or b0[axis] + lb == a0[axis]:
                    others = [other for other in range(3) if other != axis]
                    if all(min(a0[o] + la, b0[o] + lb) > max(a0[o], b0[o]) for o in others): return True
            return False
        expected_edges = {(i, j) for i in range(len(keys)) for j in range(i + 1, len(keys)) if face_adjacent(keys[i], keys[j])}
        actual_edges = {tuple(edge) for edge in edges if isinstance(edge, list) and len(edge) == 2}
        if expected_edges != actual_edges:
            raise ValueError("source adaptive face adjacency differs from independent dyadic boxes")
        matrix = {(row, column): value for row, column, value in triplets}
        for i in range(len(keys)):
            if matrix.get((i, i)) != sum(i in edge for edge in actual_edges):
                raise ValueError("source Laplacian diagonal differs from independent edge degree")
        if any(matrix.get((i, j)) != (-1.0 if (min(i, j), max(i, j)) in actual_edges else 0.0)
               for i in range(len(keys)) for j in range(len(keys)) if i != j):
            raise ValueError("source Laplacian off-diagonal differs from independent edge graph")
        passes = solve.get("passes")
        if not isinstance(passes, list) or len(passes) != 3 or any(item.get("status") != "CONVERGED" for item in passes):
            raise ValueError("adaptive split run lacks three converged passes")
        saw_split = False
        for item in passes:
            before = {(row["level"], tuple(row["coordinates"])): beta for row, beta in zip(item.get("before", ()), item.get("beta", ())) }
            after = {(row["level"], tuple(row["coordinates"])): beta for row, beta in zip(item.get("after", ()), item.get("output_beta", ())) }
            changes = item.get("split", {}).get("changes", [])
            expected = set(before)
            for change in changes:
                parent = (change["parent"]["level"], tuple(change["parent"]["coordinates"]))
                children = {(child["level"], tuple(child["coordinates"])) for child in change["children"]}
                expected_children = {(parent[0] + 1, tuple(
                    2 * parent[1][axis] + ((position >> axis) & 1) for axis in range(3)))
                    for position in range(8)}
                if parent not in before or not children or not children <= expected_children:
                    raise ValueError("adaptive split lineage is malformed")
                if any(after.get(child) != before[parent] for child in children):
                    raise ValueError("adaptive split does not preserve the parent fill value")
                expected.remove(parent); expected.update(children); saw_split = True
            if expected != set(after) or item.get("merge", {}).get("changes"):
                raise ValueError("adaptive split pass has an unexpected structural transition")
        merge_paths = repeated_paths("spin2014.source_method", "Spin-It 2014 source method")
        saw_merge = any(_load_json(path.parent / "native-run.json").get("solve", {}).get("passes", [{}])[0].get("merge", {}).get("changes") for path in merge_paths)
        if not saw_split or not saw_merge or passes[-1].get("paper_resolution_complete") is not True:
            raise ValueError("adaptive evidence lacks registered split/merge or paper-resolution witnesses")
        record = {"checker": {"path": str(Path(__file__).resolve()), "sha256": _sha256_file(Path(__file__)), "size": Path(__file__).stat().st_size},
                  "evidence": [{"path": str(path), "sha256": _sha256_file(path), "size": path.stat().st_size} for path in paths],
                  "scope": "SYNTHETIC_SOURCE_TREE_SPLIT_MERGE_AND_UNIFORM_LAPLACIAN_ONLY_NOT_AUTHOR_EXPERIMENT"}
        verified.update({name: dict(record) for name in ("spin2014.octree.boundary_subtrees", "spin2014.octree.split_merge", "spin2014.regularizer.laplacian")})

    if "spin2014.supplement_algorithms" in requests:
        paths = repeated_paths("spin2014.supplement_algorithms", "Spin-It 2014 Supplement")
        strict_root = str(Path(__file__).resolve().parent)
        sys.path.insert(0, strict_root)
        try:
            from run_spin2014_supplement_algorithms import (  # local strict runner; no user data
                check_evidence_freshness, canonical_json, file_record,
            )
        finally:
            sys.path.pop(0)
        observations = []
        for path in paths:
            assessment = check_evidence_freshness(path)
            if assessment["status"] != "CURRENT":
                raise ValueError(assessment["status"] + ": " + "; ".join(assessment["mismatches"]))
            evidence = _load_json(path)
            result = _load_json(Path(evidence["outputs"]["result"]["path"]))
            observations.append((evidence, result))
        first, second = observations
        if Path(first[0]["outputs"]["result"]["path"]).samefile(second[0]["outputs"]["result"]["path"]):
            raise ValueError("two evidence files refer to one result")
        if canonical_json(first[1]) != canonical_json(second[1]):
            raise ValueError("repeated scientific results differ")
        for key in ("python_executable", "python_version"):
            if not first[0].get("runtime", {}).get(key) or first[0]["runtime"][key] != second[0].get("runtime", {}).get(key):
                raise ValueError("repeat runtime identity mismatch")
        scientific = first[1].get("scientific_data")
        if not isinstance(scientific, dict) or scientific.get("outcome") != "REFERENCE_AND_NATIVE_AGREE":
            raise ValueError("Supplement scientific outcome is malformed")
        equation4_checks = scientific.get("equation4_checks")
        expected_paper_order = ("M", "LX", "LY", "LZ", "QXY", "QYZ", "QXZ", "QXX", "QYY", "QZZ")
        if not isinstance(equation4_checks, list) or len(equation4_checks) < 2 or any(
                not isinstance(check, dict) or
                tuple(check.get("paper_term_order", ())) != expected_paper_order or
                check.get("reference_exact_max_absolute_error", float("inf")) > 2e-12 or
                check.get("native_exact_max_absolute_error", float("inf")) > 2e-12 or
                check.get("reference_native_max_absolute_error", float("inf")) > 2e-12
                for check in equation4_checks):
            raise ValueError("Equation (4) measurements are malformed or exceed tolerance")
        fixtures = scientific.get("fixtures")
        if not isinstance(fixtures, list) or len(fixtures) < 2 or any(
                not isinstance(fixture, dict) or fixture.get("moment_max_absolute_error", float("inf")) > 2e-12 or
                fixture.get("derivative_max_absolute_error", float("inf")) > 4e-12
                for fixture in fixtures):
            raise ValueError("Supplement fixture measurements are malformed or exceed tolerance")
        record = {"checker": file_record(Path(__file__).resolve().parent / "run_spin2014_supplement_algorithms.py"),
                  "evidence": [file_record(p) for p in paths],
                  "scientific_data_sha256": first[1]["scientific_data_sha256"],
                  "scope": "PUBLISHED_SPIN2014_EQ4_AND_SUPPLEMENT_ALGORITHMS_1_2_LOCAL_SURFACE_KERNEL_ONLY"}
        verified.update({name: dict(record) for name in
                         ("spin2014.eq.4", "spin2014.supplement.algorithm.1",
                          "spin2014.supplement.algorithm.2")})

    if "spin2014.givens_equations" in requests:
        paths = repeated_paths("spin2014.givens_equations", "Spin-It 2014 Givens equations")
        strict_root = str(Path(__file__).resolve().parent)
        sys.path.insert(0, strict_root)
        try:
            from run_spin2014_givens_equations import (  # local strict runner; no user data
                check_evidence_freshness, canonical_json, file_record, SCOPED_REQUIREMENTS,
            )
        finally:
            sys.path.pop(0)
        observations = []
        for path in paths:
            assessment = check_evidence_freshness(path)
            if assessment["status"] != "CURRENT":
                raise ValueError(assessment["status"] + ": " + "; ".join(assessment["mismatches"]))
            evidence = _load_json(path)
            result = _load_json(Path(evidence["outputs"]["result"]["path"]))
            observations.append((evidence, result))
        first, second = observations
        if Path(first[0]["outputs"]["result"]["path"]).samefile(second[0]["outputs"]["result"]["path"]):
            raise ValueError("two evidence files refer to one result")
        if canonical_json(first[1]) != canonical_json(second[1]):
            raise ValueError("repeated scientific results differ")
        for key in ("python_executable", "python_version"):
            if not first[0].get("runtime", {}).get(key) or first[0]["runtime"][key] != second[0].get("runtime", {}).get(key):
                raise ValueError("repeat runtime identity mismatch")
        scientific = first[1].get("scientific_data")
        if not isinstance(scientific, dict) or scientific.get("outcome") != "REFERENCE_AND_NATIVE_GIVENS_VALUES_AGREE" or \
                tuple(scientific.get("requirements", ())) != tuple(SCOPED_REQUIREMENTS):
            raise ValueError("Givens equation scientific outcome is malformed")
        stability_checks = scientific.get("stability_checks")
        if not isinstance(stability_checks, list) or len(stability_checks) < 2 or any(
                not isinstance(check, dict) or
                check.get("max_absolute_error", float("inf")) > 1e-13
                for check in stability_checks):
            raise ValueError("Eq. (1) stability measurements are malformed or exceed tolerance")
        checks = scientific.get("checks")
        expected_objectives = {
            "YOYO", "TOP_INERTIA_ONLY", "TOP_LOWERING_ONLY",
            "TOP_MASS_AND_LOWERING", "STATIC_LOWERING",
        }
        if not isinstance(checks, list) or len(checks) != 10 or any(
                not isinstance(check, dict) or check.get("objective_max_absolute_error", float("inf")) > 1e-13 or
                check.get("constraint_max_absolute_error", float("inf")) > 1e-13
                for check in checks):
            raise ValueError("Givens equation measurements are malformed or exceed tolerance")
        by_fixture: dict[str, set[str]] = {}
        for check in checks:
            fixture_id = check.get("fixture_id")
            objective = check.get("objective")
            if not isinstance(fixture_id, str) or objective not in expected_objectives:
                raise ValueError("Givens objective measurement lacks a known fixture/objective identity")
            by_fixture.setdefault(fixture_id, set()).add(objective)
        if len(by_fixture) != 2 or any(objectives != expected_objectives for objectives in by_fixture.values()):
            raise ValueError("Givens evidence does not cover all five objectives on both fixtures")
        differential = scientific.get("native_differential_measurement")
        if not isinstance(differential, dict) or set(differential) != expected_objectives or any(
                not isinstance(value, dict) or
                value.get("active_dimension") != (10 if name == "STATIC_LOWERING" else 11) or
                value.get("gradient_max_absolute_error", float("inf")) > 2e-7 or
                value.get("hessian_max_absolute_error", float("inf")) > 3e-6
                for name, value in differential.items()):
            raise ValueError("Givens native objective derivative measurements are incomplete or exceed tolerance")
        record = {"checker": file_record(Path(__file__).resolve().parent / "run_spin2014_givens_equations.py"),
                  "evidence": [file_record(p) for p in paths],
                  "scientific_data_sha256": first[1]["scientific_data_sha256"],
                  "scope": "PUBLISHED_SPIN2014_EQ1_TO_EQ3_EQ5_TO_EQ8_AND_STATIC_GIVENS_ONLY"}
        verified.update({name: dict(record) for name in SCOPED_REQUIREMENTS})
    return verified


def _saved_planar_material_support(raw_paths, *, static=False):
    if not isinstance(raw_paths, list) or len(raw_paths) != 2 or \
            not all(isinstance(path, str) for path in raw_paths):
        raise ValueError("saved planar material requires two independent checker evidence paths")
    paths = [Path(path).resolve(strict=True) for path in raw_paths]
    if paths[0].samefile(paths[1]):
        raise ValueError("repeated saved planar evidence is not two runs")
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]/"python"))
    try:
        from make_it_stand.verification.run_evidence import PaperRunEvidenceV2, check_freshness
        if static:
            from make_it_stand.verification.static_material import audit_saved_static_application as checker, SCOPE
        else:
            from make_it_stand.verification.planar_material import audit_saved_planar_application as checker, SCOPE
    finally:
        sys.path.pop(0)
    reports, inputs = [], []
    for path in paths:
        evidence = PaperRunEvidenceV2.from_json(path.read_text(encoding="utf-8"))
        if check_freshness(evidence).status != "CURRENT":
            raise ValueError("STALE_RESULT: saved planar checker evidence is not current")
        value = evidence.as_dict()
        request = value["begin"]["request"]
        operation = "READ_SAVED_STATIC_GEOMETRY" if static else "READ_SAVED_PLANAR_GEOMETRY"
        if request["parameters"] != {"scope": SCOPE, "operation": operation} or \
                request["pure_python"] is not True or len(request["asset_files"]) != 1:
            raise ValueError("saved planar checker has invalid operation or inputs")
        application = Path(request["asset_files"][0]).resolve(strict=True)
        inputs.append(application)
        actual = checker(application)
        scientific = path.parent/"scientific.json"
        outputs = {Path(row["path"]).resolve() for row in value["scientific"]["outputs"]}
        if scientific.resolve() not in outputs or _load_json(scientific) != actual:
            raise ValueError("saved planar report differs from independently recomputed geometry")
        # Do not accept a file containing only failures as successful supporting
        # geometry evidence. Failures remain readable in the underlying record.
        field = "terminal_numerics_qualified" if static else "balance_constraints_pass"
        if not any(row["moments_match"] and row[field]
                   for case in actual["cases"] for row in case["attempts"]):
            raise ValueError("saved planar evidence has no independently balanced material")
        reports.append(actual["cases"])
    if inputs[0].samefile(inputs[1]):
        raise ValueError("repeated checks need distinct saved application outputs")
    if reports[0] != reports[1]:
        raise ValueError("repeated independent saved planar results differ")
    return {"scope": SCOPE,
            "checker": {"path": str(Path(__file__).resolve()), "sha256": _sha256_file(Path(__file__))},
            "evidence": [{"path": str(path), "sha256": _sha256_file(path)} for path in paths],
            "recomputed_cases": reports[0],
            "verification_effect": "SUPPORTING_ONLY_NO_CATALOGUE_REQUIREMENT_VERIFIED"}


def _assert_distinct_cherchi_native_runs(applications,bundles):
    """Different wrappers/roots are insufficient: compare original run files."""
    seen_inodes=set();seen_records=set()
    for application,bundle in zip(applications,bundles,strict=True):
        root=Path(bundle).resolve(strict=True)
        for case in _load_json(Path(application))["cases"]:
            if case.get("solver_family")!="CHERCHI_RAW_EQ14_V1":continue
            for row in case["runs"]:
                if row["status"]=="INFEASIBLE_CERTIFIED":continue
                path=(root/row["relative_output"]/"evidence.json").resolve(strict=True)
                if not path.is_relative_to(root):raise ValueError("native evidence escapes registered root")
                stat=path.stat();inode=(stat.st_dev,stat.st_ino);digest=_sha256_file(path)
                if inode in seen_inodes or digest in seen_records:
                    raise ValueError("reused native evidence is the same native run, not an independent repeat")
                seen_inodes.add(inode);seen_records.add(digest)


@contextmanager
def _cherchi_bound_runtime(evidence,check_freshness):
    modules=evidence.as_dict()["begin"]["observation"]["runtime"]["native_modules"]
    if set(modules)!={"_make_it_stand"}:
        raise ValueError("saved Cherchi evidence lacks its supported native module identity")
    observed=modules["_make_it_stand"];path=Path(observed["path"]).resolve(strict=True)
    if observed["sha256"]!=_sha256_file(path) or observed["size"]!=path.stat().st_size:
        raise ValueError("STALE_RESULT: saved Cherchi native binary differs from evidence")
    previous=sys.modules.get("_make_it_stand")
    loaded,freshness=_load_bound_native_module_and_freshness(evidence,path,check_freshness)
    sys.modules["_make_it_stand"]=loaded
    try:
        if freshness.status!="CURRENT":
            raise ValueError("STALE_RESULT: saved Cherchi checker evidence is not current")
        yield
    finally:
        sys.modules.pop("_make_it_stand",None)
        if previous is not None:sys.modules["_make_it_stand"]=previous


def _saved_cherchi_material_support(raw_paths):
    if not isinstance(raw_paths, list) or len(raw_paths) != 2 or not all(isinstance(p,str) for p in raw_paths):
        raise ValueError("saved Cherchi material requires two independent checker evidence paths")
    paths=[Path(p).resolve(strict=True) for p in raw_paths]
    if paths[0].samefile(paths[1]):
        raise ValueError("repeated saved Cherchi evidence is not two runs")
    sys.path.insert(0,str(Path(__file__).resolve().parents[3]/"python"))
    try:
        from make_it_stand.verification.run_evidence import PaperRunEvidenceV2,check_freshness
        from make_it_stand.verification.cherchi_application import audit_saved_cherchi_application,SCOPE
    finally:
        sys.path.pop(0)
    reports=[]; applications=[]; bundles=[]
    for path in paths:
        evidence=PaperRunEvidenceV2.from_json(path.read_text())
        value=evidence.as_dict();request=value["begin"]["request"]
        params=request["parameters"]
        if set(params)!={"scope","operation","bundle_root"} or params["scope"]!=SCOPE or \
                params["operation"]!="READ_SAVED_CHERCHI_APPLICATION" or request["pure_python"] is not False or \
                len(request["asset_files"])!=1:
            raise ValueError("saved Cherchi checker has invalid operation or inputs")
        application=Path(request["asset_files"][0]).resolve(strict=True)
        bundle=Path(params["bundle_root"]).resolve(strict=True)
        with _cherchi_bound_runtime(evidence,check_freshness):
            actual=audit_saved_cherchi_application(application,bundle_root=bundle)
        outputs={Path(row["path"]).resolve() for row in value["scientific"]["outputs"]}
        scientific=path.parent/"scientific.json"
        if scientific.resolve() not in outputs or _load_json(scientific)!=actual:
            raise ValueError("saved Cherchi report differs from independently reloaded material")
        if not any(case["qualified_material_runs"]>0 for case in actual["cases"]):
            raise ValueError("saved Cherchi evidence has no independently qualified material")
        reports.append(actual);applications.append(application);bundles.append(bundle)
    if applications[0].samefile(applications[1]) or bundles[0].samefile(bundles[1]):
        raise ValueError("repeated checks require distinct saved applications and native runs")
    _assert_distinct_cherchi_native_runs(applications,bundles)
    if reports[0]!=reports[1]:
        raise ValueError("repeated independent saved Cherchi results differ")
    return {"scope":SCOPE,
        "checker":{"path":str(Path(__file__).resolve()),"sha256":_sha256_file(Path(__file__))},
        "evidence":[{"path":str(p),"sha256":_sha256_file(p)} for p in paths],
        "recomputed_cases":reports[0]["cases"],
        "verification_effect":"SUPPORTING_ONLY_NO_CATALOGUE_REQUIREMENT_VERIFIED"}


def _saved_buoyant_material_support(raw_paths):
    if not isinstance(raw_paths, list) or len(raw_paths) != 2 or not all(isinstance(p,str) for p in raw_paths):
        raise ValueError("saved buoyant material requires two independent checker evidence paths")
    paths = [Path(p).resolve(strict=True) for p in raw_paths]
    if paths[0].samefile(paths[1]):
        raise ValueError("repeated buoyant checker evidence is not two independent runs")
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]/"python"))
    try:
        from make_it_stand.verification.run_evidence import PaperRunEvidenceV2, check_freshness
        from make_it_stand.verification.buoyant_evidence import (
            audit_buoyant_execution, assert_independent_buoyant_runs, scientific_signature, SCOPE, OPERATION)
    finally:
        sys.path.pop(0)
    reports = []
    for path in paths:
        evidence = PaperRunEvidenceV2.from_json(path.read_text())
        value = evidence.as_dict()
        request = value['begin']['request']
        params = request['parameters']
        keys = {'application', 'bundle_root', 'execution_evidence', 'source_job_root'}
        if set(params) != keys | {'scope', 'operation'} or params['scope'] != SCOPE or \
                params['operation'] != OPERATION or request['pure_python'] is not False:
            raise ValueError('saved buoyant checker has invalid operation or inputs')
        # The old helper name is historical; it binds the same native module for
        # Cherchi and buoyancy rather than trusting the ambient Python import.
        with _cherchi_bound_runtime(evidence, check_freshness):
            report = audit_buoyant_execution(**{key:params[key] for key in keys})
        scientific = path.parent/'scientific.json'
        outputs = {Path(row['path']).resolve() for row in value['scientific']['outputs']}
        if scientific.resolve() not in outputs or _load_json(scientific) != report:
            raise ValueError('saved buoyant report differs from recomputed execution and material')
        if not any(case['qualified_runs'] > 0 for case in report['cases']):
            raise ValueError('saved buoyant evidence has no independently qualified material')
        reports.append(report)
    assert_independent_buoyant_runs([report['provenance'] for report in reports])
    if scientific_signature(reports[0]) != scientific_signature(reports[1]):
        raise ValueError('repeated independent buoyant scientific results differ')
    return {'scope':SCOPE, 'evidence':[{'path':str(p),'sha256':_sha256_file(p)} for p in paths],
            'recomputed_cases':reports[0]['cases'],
            'native_runs':[report['provenance'] for report in reports],
            'verification_effect':'SUPPORTING_ONLY_NO_CATALOGUE_REQUIREMENT_VERIFIED'}


def _supporting_checks(requests: Any) -> dict[str, dict[str, Any]]:
    """Validate evidence that is useful but insufficient to verify a requirement.

    Supporting evidence is deliberately a separate channel: records returned
    here never mutate a catalogue requirement to VERIFIED.  This prevents a
    finite numerical witness from being represented as a proof of the Faster
    quadric-reduction or optimality theorems.
    """
    requests = _mapping(requests, "supporting evidence")
    registered = {"faster.theorem_numerical_conditions", "faster.planar_saved_material", "faster.static_saved_material", "faster.cherchi_saved_material", "faster.buoyant_saved_material"}
    if set(requests) - registered:
        raise ValueError("unregistered supporting checker")
    if not requests:
        return {}

    supporting = {}
    if "faster.buoyant_saved_material" in requests:
        supporting["faster.buoyant_saved_material"] = _saved_buoyant_material_support(
            requests["faster.buoyant_saved_material"])
    if "faster.cherchi_saved_material" in requests:
        supporting["faster.cherchi_saved_material"] = _saved_cherchi_material_support(
            requests["faster.cherchi_saved_material"])
    if "faster.static_saved_material" in requests:
        supporting["faster.static_saved_material"] = _saved_planar_material_support(
            requests["faster.static_saved_material"], static=True)
    if "faster.planar_saved_material" in requests:
        supporting["faster.planar_saved_material"] = _saved_planar_material_support(
            requests["faster.planar_saved_material"])
    if "faster.theorem_numerical_conditions" not in requests:
        return supporting

    paths = requests["faster.theorem_numerical_conditions"]
    if not isinstance(paths, list) or len(paths) != 2 or not all(isinstance(p, str) for p in paths):
        raise ValueError("Faster theorem numerical conditions require two independent evidence paths")
    paths = [Path(path).resolve(strict=True) for path in paths]
    if paths[0].samefile(paths[1]):
        raise ValueError("repeated supporting evidence is not two runs")
    python_root = str(Path(__file__).resolve().parents[3] / "python")
    sys.path.insert(0, python_root)
    try:
        from make_it_stand.verification.run_evidence import PaperRunEvidenceV2, check_freshness
        from make_it_stand.verification.faster_theorem_exact import verify_certificate
    finally:
        sys.path.pop(0)
    scope = ("SYNTHETIC_EXACT_RATIONAL_THEOREM_INSTANCES_AND_CANONICAL_EQ7_"
             "MATERIAL_METHODS_ONLY_NOT_GENERAL_THEOREM_PROOF_OR_AUTHOR_CASE")
    expected = {
        "theorem-conditions": ((), (
            "[PASS] Eq. 7 numerical hypotheses accept exact fixture",
            "[PASS] Eq. 7 numerical hypotheses reject infeasible and nonpositive-mu fixtures",
            "[PASS] finite signed variation certificate proves a strict discrete increase",
            "[PASS] finite beta-zero cells never stand in for a measure-zero proof",
            "[PASS] theorem audit scope remains numerical only",
        )),
            "canonical-replay": (("replay",), (
                "same-q replay geometry=1 moment-diff=0 J-diff=0",
                "[PASS] replay",
            )),
    }
    signatures = []
    for path in paths:
        evidence = PaperRunEvidenceV2.from_json(path.read_text(encoding="utf-8"))
        if check_freshness(evidence).status != "CURRENT":
            raise ValueError("STALE_RESULT: Faster theorem numerical-condition evidence is not current")
        root = path.parent
        scientific, summary = _load_json(root / "scientific.json"), _load_json(root / "summary.json")
        exact_certificate_path, exact_result_path = root / "exact-certificate.json", root / "exact-instances.json"
        exact_certificate, exact_result = _load_json(exact_certificate_path), _load_json(exact_result_path)
        expected_exact_result = verify_certificate(exact_certificate)
        outputs = {Path(row["path"]).resolve() for row in evidence.as_dict()["scientific"]["outputs"]}
        commands = scientific.get("commands")
        if scientific.get("schema") != "faster-theorem-numerical-conditions-scientific-v2" or \
                scientific.get("classification") != "SYNTHETIC_REGRESSION" or \
                scientific.get("scope") != scope or \
                scientific.get("theorem_coverage") != "NOT_GRANTED_BY_THIS_RUNNER" or \
                scientific.get("exact_instance_coverage") != "EXACT_INSTANCE_CONDITIONS_VERIFIED" or \
                scientific.get("exact_instance_result_sha256") != _sha256_file(exact_result_path) or \
                exact_result != expected_exact_result or \
                exact_result.get("general_theorem_proof") != "NOT_CLAIMED" or \
                summary.get("status") != "SYNTHETIC_EXACT_THEOREM_INSTANCES_AND_NUMERICAL_CONDITIONS_VALIDATED" or \
                summary.get("scope") != scope or \
                summary.get("theorem_coverage") != "NOT_GRANTED_BY_THIS_RUNNER" or \
                summary.get("exact_instance_coverage") != "EXACT_INSTANCE_CONDITIONS_VERIFIED" or \
                exact_certificate_path.resolve() not in outputs or exact_result_path.resolve() not in outputs or \
                not isinstance(commands, list) or len(commands) != len(expected):
            raise ValueError("Faster theorem numerical-condition evidence is malformed")
        observed = {}
        for command in commands:
            if not isinstance(command, dict) or command.get("name") not in expected:
                raise ValueError("Faster theorem numerical-condition command is malformed")
            name = command["name"]
            arguments, markers = expected[name]
            transcript = root / f"{name}.stdout"
            if name in observed or command.get("arguments") != list(arguments) or \
                    command.get("required_markers") != list(markers) or command.get("returncode") != 0 or \
                    Path(command.get("path", "")).resolve() != transcript.resolve() or \
                    transcript.resolve() not in outputs or not transcript.is_file() or \
                    _sha256_file(transcript) != command.get("sha256"):
                raise ValueError("Faster theorem numerical-condition transcript is not bound")
            text = transcript.read_text(encoding="utf-8")
            if any(marker not in text for marker in markers):
                raise ValueError("Faster theorem numerical-condition transcript lacks a required marker")
            observed[name] = command["sha256"]
        if set(observed) != set(expected):
            raise ValueError("Faster theorem numerical-condition command set is incomplete")
        signatures.append((tuple(sorted(observed.items())), scientific["exact_instance_result_sha256"]))
    if signatures[0] != signatures[1]:
        raise ValueError("repeated Faster theorem numerical-condition scientific outputs differ")
    return {**supporting, "faster.theorem_numerical_conditions": {
        "checker": {"path": str(Path(__file__).resolve()), "sha256": _sha256_file(Path(__file__)),
                    "size": Path(__file__).stat().st_size},
        "evidence": [{"path": str(path), "sha256": _sha256_file(path), "size": path.stat().st_size}
                     for path in paths],
        "scientific_data_sha256": hashlib.sha256(
            json.dumps(signatures[0], separators=(",", ":"), ensure_ascii=True).encode("ascii")
        ).hexdigest(),
        "scope": scope,
        "verification_effect": "SUPPORTING_ONLY_NO_CATALOGUE_REQUIREMENT_VERIFIED",
    }}


def _requirement_status_index(
    requirements: list[dict[str, Any]],
    catalogue_rows: list[dict[str, Any]],
    case_asset_readiness: dict[str, dict[str, Any]],
) -> dict[str, dict[str, Any]]:
    """Return the four non-promotional status columns for every requirement.

    This is an inventory, rather than a second acceptance path.  In particular,
    a blocked caller row never becomes an implementation claim merely because
    its free-form details mention a test or an experiment.
    """

    catalogue_by_id = {
        row["requirement_id"]: row for row in catalogue_rows
    }
    index: dict[str, dict[str, Any]] = {}
    for requirement in requirements:
        requirement_id = requirement["requirement_id"]
        catalogue_row = catalogue_by_id[requirement_id]
        kind = catalogue_row["requirement_kind"]
        status = requirement["status"]
        dependencies = list(catalogue_row["dependencies"])

        if kind == "PHYSICAL_MEASUREMENT":
            index[requirement_id] = {
                "implementation": "NOT_APPLICABLE",
                "independent_verification": "NOT_IN_SCOPE",
                "original_experiment": "NOT_IN_SCOPE",
                "blocker_dependencies": dependencies,
            }
            continue

        independently_verified = status == "VERIFIED"
        if independently_verified:
            implementation = "IMPLEMENTED_AND_VERIFIED"
            verification = "VERIFIED"
        else:
            implementation = "UNASSESSED"
            verification = "NOT_VERIFIED"
            dependencies.append("independent_checker")

        if kind == "ORIGINAL_EXPERIMENT":
            readiness = case_asset_readiness.get(requirement_id)
            has_exact_input = (
                readiness is not None
                and readiness["asset_status"] in _READY_CASE_ASSET_STATUSES
            )
            if independently_verified and has_exact_input:
                original_experiment = "VERIFIED"
            elif has_exact_input:
                original_experiment = "ELIGIBLE_INPUT_AWAITING_VERIFICATION"
            else:
                original_experiment = "NOT_ELIGIBLE_MISSING_ORIGINAL_EVIDENCE"
                dependencies.extend(("asset_evidence", "critical_parameter_evidence"))
        else:
            original_experiment = "NOT_APPLICABLE"

        index[requirement_id] = {
            "implementation": implementation,
            "independent_verification": verification,
            "original_experiment": original_experiment,
            "blocker_dependencies": sorted(set(dependencies)),
        }
    return index


def audit(payload: Any, *, asset_root: Path, independent_evidence=None,
          supporting_evidence=None, case_asset_readiness: Path | None = None,
          paper: str = "all") -> dict[str, Any]:
    """Caller statuses remain untrusted; only registered checks may verify rows."""

    document = _mapping(payload, "evidence document")
    _reject_nonfinite(document, "evidence")
    _exact_fields(
        document,
        {"schema", "schema_version", "requirements", "assets"},
        "evidence document",
    )
    if document["schema"] != EVIDENCE_SCHEMA:
        raise ValueError(f"unsupported evidence schema: {document['schema']!r}")
    if document["schema_version"] != EVIDENCE_VERSION:
        raise ValueError(f"unsupported evidence schema version: {document['schema_version']!r}")

    catalogue_document = load_catalogue()
    requirements = _validate_requirement_evidence(document["requirements"])
    asset_preflight = preflight_assets(document["assets"], Path(asset_root))
    checked = _independent_checks({} if independent_evidence is None else independent_evidence)
    supporting = _supporting_checks({} if supporting_evidence is None else supporting_evidence)
    for row in requirements:
        if row["requirement_id"] in checked:
            row.update(status="VERIFIED", evidence_description="Exact published certificates recomputed; two run outputs agree.",
                       details=checked[row["requirement_id"]])

    readiness = ({"status": "NOT_SUPPLIED", "cases": {}}
                 if case_asset_readiness is None else {
                     "status": "CURRENT",
                     "source": {
                         "path": str(case_asset_readiness.resolve(strict=True)),
                         "sha256": _sha256_file(case_asset_readiness),
                         "size": case_asset_readiness.stat().st_size,
                     },
                     **load_case_asset_readiness_record(case_asset_readiness),
                 })
    requirement_status_index = _requirement_status_index(
        requirements, catalogue_document["requirements"], readiness["cases"]
    )

    paper_statuses = {
        name: derive_paper_status(name, requirements, readiness["cases"])
        for name in ("faster", "spin2014", "all")
    }
    selected_paper = paper_statuses[paper] if paper in paper_statuses else None
    if selected_paper is None:
        raise ValueError(f"invalid paper scope {paper!r}; expected faster, spin2014, or all")

    physical_count = sum(row["status"] == "NOT_IN_SCOPE" for row in requirements)
    blocked_requirements = [
        row for row in requirements if row["status"] not in {"NOT_IN_SCOPE", "VERIFIED"}
    ]
    all_blockers = [
        "independent requirement checkers are incomplete; input readiness cannot grant verification",
        "formula coverage remains blocked until each formula requirement has an independent checker result",
        "experiment coverage remains blocked until each original experiment has independently checked evidence",
    ]
    all_blockers.extend(
        f"{row['requirement_id']}: {row['status']} — {row['evidence_description']}"
        for row in blocked_requirements
    )
    all_blockers.extend(selected_paper["blockers"])
    all_blockers.extend(asset_preflight["blockers"])
    blockers = (all_blockers if paper == "all" else list(selected_paper["blockers"]))

    return {
        "schema": AUDIT_SCHEMA,
        "schema_version": AUDIT_VERSION,
        "catalogue_schema": catalogue_document["schema"],
        "catalogue_version": catalogue_document["schema_version"],
        "paper_scope": paper,
        "paper_statuses": paper_statuses,
        "status": selected_paper["status"],
        "public_complete_reproduction": selected_paper["public_complete_reproduction"],
        "formula_coverage": selected_paper["formula_coverage"],
        "experiment_coverage": selected_paper["experiment_coverage"],
        "asset_preflight": asset_preflight,
        "case_asset_readiness": readiness,
        "requirements": requirements,
        "requirement_status_index": requirement_status_index,
        "supporting_evidence": supporting,
        "requirement_summary": {
            "required": len(REQUIRED_REQUIREMENT_IDS),
            "blocked_or_not_implemented": len(blocked_requirements),
            "physical_measurements_out_of_scope": physical_count,
            "verified": len(checked),
        },
        "blockers": blockers,
        "limitations": [
            "The strict-audit foundation has partial checker integration, not full reproduction coverage.",
            "Only registered independent checkers can verify their explicitly scoped requirements.",
            "Physical fabricated-object measurements are retained explicitly and remain out of scope.",
            "READY asset preflight covers syntax, declared provenance, and checked file-byte integrity only; it does not prove source authenticity, strict verification, or publication readiness.",
        ],
    }


def _load_json(path: Path) -> Any:
    return json.loads(
        path.read_text(encoding="utf-8"),
        parse_constant=_raise_nonfinite,
    )


def _write_json(path: Path, document: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(document, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Validate strict Spin-It evidence without granting publication claims."
    )
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--asset-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--independent-evidence", type=Path,
                        help="JSON mapping registered checker IDs to evidence paths, never pass flags")
    parser.add_argument("--supporting-evidence", type=Path,
                        help="JSON mapping registered supporting checker IDs to evidence paths; never verifies a requirement")
    parser.add_argument("--case-asset-readiness", type=Path,
                        help="exhaustive original-case asset and parameter readiness inventory")
    parser.add_argument("--paper", choices=("faster", "spin2014", "all"), default="all",
                        help="report completion for one fixed paper scope; all remains the default")
    args = parser.parse_args(argv)

    try:
        result = audit(_load_json(args.evidence), asset_root=args.asset_root,
                       independent_evidence=_load_json(args.independent_evidence) if args.independent_evidence else None,
                       supporting_evidence=_load_json(args.supporting_evidence) if args.supporting_evidence else None,
                       case_asset_readiness=args.case_asset_readiness,
                       paper=args.paper)
    except (OSError, json.JSONDecodeError, ValueError) as error:
        result = {
            "schema": AUDIT_SCHEMA,
            "schema_version": AUDIT_VERSION,
            "status": "INVALID_BLOCKED",
            "public_complete_reproduction": False,
            "formula_coverage": False,
            "experiment_coverage": False,
            "blockers": [str(error)],
        }
    _write_json(args.output, result)
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0 if result.get("public_complete_reproduction") is True else 2


if __name__ == "__main__":
    raise SystemExit(main())
