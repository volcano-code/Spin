#!/usr/bin/env python3
"""Synthetic integration diagnostic: pinned Cherchi Boolean inside Eq.14 Newton.

This is not a paper case. Each Newton evaluation rebuilds the quadric source,
invokes the pinned reference Boolean, recovers construction motion, and
computes moments/Jacobian from that same Cherchi output. A converged run also
constructs the editable-domain complement through Cherchi and seals/reloads
the resulting locked material artifact; numerical convergence alone is not a
sealed geometry result.
"""
from __future__ import annotations

import argparse
import csv
from contextlib import contextmanager
from dataclasses import replace
from functools import cmp_to_key
from fractions import Fraction
import hashlib
import math
import json
import os
from pathlib import Path
import re
import signal
import struct
import subprocess
import sys
import time
from typing import Callable

import numpy as np

from make_it_stand import spin
from run_faster_reference_boolean import (
    _FIXTURES,
    _construction_topology_signature,
    _read_obj,
    _summarize_vertex_construction_provenance,
    _write_obj,
)
from make_it_stand.verification.cherchi_motion_differential import (
    CherchiDifferentialUnavailable,
    differentiate_polyhedral_raw_moments,
    integrate_rational_polyhedral_raw_moments,
    normalize_closed_component_orientation,
    recover_construction_vertex_motion,
)
from make_it_stand.verification.mass_moments import (
    TriangleMesh,
    compute_mass_properties,
)
from make_it_stand.verification.run_evidence import (
    EXPLICIT_EFFECTIVE_SOURCE_POLICY,
    RunInputsV2,
    begin_run,
    check_freshness,
    finish_run,
    native_build_configuration_files,
)
from make_it_stand.verification import (
    SpinArtifactProblem,
    load_spin_artifact,
    save_spin_artifact,
    solve_cherchi_geometry_and_seal_faster,
)
from make_it_stand.verification.faster_infeasibility import (
    verify_axis_aligned_yoyo_moment_capacity,
)
from make_it_stand.verification.faster_material import (
    _prepare, _raw_moments, _json_value, mesh_snapshot_from_exact_boolean_sidecar,
    read_exact_boolean_sidecar_coordinates, rotational_objective_definition,
    material_run_scientific_identity,
)
from make_it_stand.verification.spin_problem import MeshSnapshot
from make_it_stand.verification.rational_text import fraction_from_rational_text
from make_it_stand.verification.spin_common import digest

ROOT = Path(__file__).resolve().parents[3]


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _prepare_exact_motion_material(vertices, faces, sidecar):
    """Authenticate raw export, orient exact geometry, then validate its solid.

    Keep original OBJ coordinates and sidecar derivation bound together.
    Canonical display is only the carrier for rational geometry computations;
    it does not grant float-export or manufacturing approval. Area/orientation
    success never bypasses full native embedding and material validation.
    """
    rows, derivation = read_exact_boolean_sidecar_coordinates(vertices, sidecar)
    canonical = np.asarray([[float(fraction_from_rational_text(value)) for value in row]
                            for row in rows], dtype=np.float64)
    oriented, orientation = normalize_closed_component_orientation(
        canonical, faces, rational_vertices=rows)
    snapshot = MeshSnapshot.from_rationals(rows, oriented)
    derivation['rational_geometry_identity'] = snapshot.identity
    return snapshot, derivation, orientation


def _terminal_exact_material(vertices, faces, sidecar, density):
    """Authenticate the exact complement without granting its display validity.

    Terminal complement geometry has no Newton derivative. Preserve its
    original faces; the native rational checker must independently validate
    winding, closure and embedding before any moments are used. A collapsed
    float display stays explicitly rejected, not repaired or differentiated.
    """
    snapshot, derivation = mesh_snapshot_from_exact_boolean_sidecar(
        vertices, faces, sidecar, return_audit=True)
    try:
        _, float_orientation = normalize_closed_component_orientation(vertices, faces)
        float_check = {'status': 'ORIENTATION_ONLY_PASSED', 'audit': float_orientation}
    except CherchiDifferentialUnavailable as error:
        float_check = {'status': 'REJECTED', 'diagnostic': str(error)}
    orientation = {
        'schema': 'terminal-exact-orientation-audit-v1',
        'authority': 'EXACT_RATIONAL_BOUNDARY',
        'face_order_changed': False,
        'floating_orientation_check': float_check,
        'scope': 'EXACT_BOUNDARY_NOT_DISPLAY_OR_MANUFACTURING_ACCEPTANCE',
    }
    moments = np.asarray(_raw_moments(snapshot, density), dtype=np.float64)
    return snapshot, derivation, orientation, moments


def _terminal_audit_checker_paths() -> list[Path]:
    paths = [Path(__file__).resolve(), Path(spin._native_module().__file__).resolve(),
             Path(spin.__file__).resolve(),
             Path(__file__).with_name("run_faster_reference_boolean.py").resolve()]
    paths.extend((ROOT / "python/make_it_stand/verification").rglob("*.py"))
    return sorted({path.resolve() for path in paths})


def _save_terminal_audit_checkpoint(output: Path, payload: dict, *,
                                    dependencies: list[Path]) -> Path:
    """Persist audit inputs before the expensive exact terminal checks.

    This is an audit-only checkpoint, never a completed scientific run.
    All supplied inputs plus current checker sources and extension are bound.
    """
    checkpoint = output / "terminal-audit-checkpoint.json"
    if checkpoint.exists():
        raise ValueError("terminal audit checkpoint already exists")
    source_paths = _terminal_audit_checker_paths()
    paths = sorted({path.resolve() for path in [*dependencies, *source_paths]})
    record = {
        "schema": "faster-terminal-audit-checkpoint-v1",
        "status": "PENDING_TERMINAL_AUDIT",
        "complete_run_evidence": False,
        "checker_identity": [str(path) for path in source_paths],
        "payload": _json_safe(payload),
        "dependencies": [{"path": str(path), "sha256": _sha(path)} for path in paths],
    }
    record["checkpoint_digest"] = digest(record)
    temporary = checkpoint.with_suffix(".json.pending")
    with temporary.open("x", encoding="utf-8") as stream:
        stream.write(json.dumps(record, sort_keys=True, separators=(",", ":"),
                                allow_nan=False) + "\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, checkpoint)
    directory_fd = os.open(output, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(directory_fd)
    finally:
        os.close(directory_fd)
    return checkpoint


def _load_terminal_audit_checkpoint(checkpoint: Path) -> dict:
    record = json.loads(checkpoint.read_text(encoding="utf-8"))
    claimed_digest = record.pop("checkpoint_digest", None)
    if claimed_digest != digest(record):
        raise ValueError("terminal audit checkpoint digest mismatch")
    if (record.get("schema") != "faster-terminal-audit-checkpoint-v1" or
            record.get("status") != "PENDING_TERMINAL_AUDIT" or
            record.get("complete_run_evidence") is not False):
        raise ValueError("invalid terminal audit checkpoint state")
    dependencies = record.get("dependencies")
    if not isinstance(dependencies, list) or not dependencies:
        raise ValueError("terminal audit checkpoint lacks dependencies")
    if record.get("checker_identity") != [str(path) for path in _terminal_audit_checker_paths()]:
        raise ValueError("terminal audit checker identity mismatch")
    if not set(record["checker_identity"]).issubset(
            {dependency["path"] for dependency in dependencies}):
        raise ValueError("terminal audit checker identity lacks content dependencies")
    for dependency in dependencies:
        path = Path(dependency["path"])
        if not path.is_file() or _sha(path) != dependency["sha256"]:
            raise ValueError(f"terminal audit dependency changed: {path}")
    if not isinstance(record.get("payload"), dict):
        raise ValueError("terminal audit checkpoint lacks payload")
    return record["payload"]


def resume_terminal_audit(output: Path) -> dict:
    """Re-audit persisted geometry without rerunning Newton or issuing run evidence."""
    started = time.perf_counter()
    checkpoint = output / "terminal-audit-checkpoint.json"
    payload = _load_terminal_audit_checkpoint(checkpoint)
    audit = _audit_terminal_geometry_binding(
        payload, quadric_hash=spin._native_module().spin_it_faster_quadric_hash)
    # Fail closed if files changed while the audit was running as well.
    _load_terminal_audit_checkpoint(checkpoint)
    receipt = {
        "schema": "faster-resumed-terminal-audit-v1",
        "classification": "INDEPENDENT_TERMINAL_AUDIT_ONLY",
        "status": "TERMINAL_AUDIT_PASSED" if audit["status"] == "PASSED"
                  else "TERMINAL_AUDIT_FAILED",
        "complete_run_evidence": False,
        "checkpoint_sha256": _sha(checkpoint),
        "audit": audit,
        "wall_seconds": time.perf_counter() - started,
    }
    receipt_path = output / "resumed-terminal-audit.json"
    with receipt_path.open("x", encoding="utf-8") as stream:
        stream.write(json.dumps(_json_safe(receipt), sort_keys=True,
                                separators=(",", ":"), allow_nan=False) + "\n")
        stream.flush()
        os.fsync(stream.fileno())
    return receipt


def resume_run_finalization(output: Path) -> dict:
    """Finish a staged run only from its own checked context and persisted bundle."""
    started = time.perf_counter()
    output = output.absolute()
    if any((output / name).exists() for name in ("result.json", "evidence.json", "freshness.json")):
        raise ValueError("full-run finalization refuses existing completed outputs")
    checkpoint = output / "terminal-audit-checkpoint.json"
    payload = _load_terminal_audit_checkpoint(checkpoint)
    context = payload.get("run_finalization")
    if not isinstance(context, dict):
        raise ValueError("checkpoint lacks full-run finalization context")
    if context.get("schema") != "faster-run-finalization-context-v1":
        raise ValueError("unsupported full-run finalization context")
    if Path(context["output_directory"]).resolve() != output.resolve():
        raise ValueError("full-run finalization output identity mismatch")
    evidence_inputs = RunInputsV2(json.dumps(context["begin"], allow_nan=False))
    expected_digest = context["expected_artifact_digest"]
    artifact = None
    artifact_bundle = None
    if expected_digest is not None:
        artifact_bundle = output / "sealed-material-artifact"
        artifact = load_spin_artifact(
            artifact_bundle, expected_problem_digest=context["locked_problem_digest"])
        if artifact.scientific_digest != expected_digest:
            raise ValueError("full-run finalization sealed artifact digest mismatch")
    elif context["native_result"]["status"] == "CONVERGED":
        raise ValueError("converged full-run finalization lacks its sealed artifact")
    _load_terminal_audit_checkpoint(checkpoint)
    return _finalize_run(
        output, payload=payload, context=context, evidence_inputs=evidence_inputs,
        sealed_artifact=artifact, artifact_bundle=artifact_bundle,
        artifact_reload_digest=None if artifact is None else artifact.scientific_digest,
        started=started, prior_wall_seconds=context["solver_wall_seconds"],
        execution_protocol="CHECKPOINT_RESUMED_FINALIZATION_V1",
    )


def _evidence_events(result_record: dict[str, object]) -> list[dict[str, object]]:
    """Preserve every actual Boolean callback, including failed trial records."""
    evaluations = result_record.get("evaluations")
    if not isinstance(evaluations, list):
        raise ValueError("run evidence requires the ordered evaluation ledger")
    events: list[dict[str, object]] = []
    for fallback_index, raw in enumerate(evaluations):
        if not isinstance(raw, dict):
            raise ValueError("run evidence evaluation ledger contains a non-object row")
        index = raw.get("evaluation_index", fallback_index)
        if isinstance(index, bool) or not isinstance(index, int) or index < 0:
            raise ValueError("run evidence evaluation index must be a nonnegative integer")
        moment_map = raw.get("moment_map")
        successful = isinstance(moment_map, dict) and isinstance(moment_map.get("sha256"), str)
        data: dict[str, object] = {
            "evaluation_index": index,
            "quadric_coefficients": raw.get("quadric_coefficients"),
            "boolean_operation": raw.get("boolean_operation", raw.get("operation")),
            "geometry_hash": raw.get("geometry_hash"),
            "topology_hash": raw.get("topology_hash"),
        }
        if successful:
            data["moment_map_sha256"] = moment_map["sha256"]
            data["independent_moment_audit"] = raw.get("independent_moment_audit")
        else:
            data["stage"] = raw.get("stage", "evaluation_failed")
            data["diagnostic"] = raw.get("diagnostic", raw.get("error", "failed trial without diagnostic"))
        events.append({"iteration": index, "status": "ok" if successful else "failed", "data": data})
    if not events:
        events.append({"iteration": 0, "status": "failed", "data": {
            "stage": "solver",
            "diagnostic": str(result_record.get("diagnostic", "no Boolean evaluation was recorded")),
            "solver_status": str(result_record.get("solver_status", "UNKNOWN")),
        }})
    return events


def _compact_exact_fraction(value: Fraction, *, decimal_safe_bits: int = 1024) -> str:
    """Serialize exact rationals without triggering huge-int decimal limits."""
    numerator, denominator = value.numerator, value.denominator
    if max(abs(numerator).bit_length(), denominator.bit_length()) <= decimal_safe_bits:
        return f"{numerator}/{denominator}"
    digest = hashlib.sha256()
    for integer in (numerator, denominator):
        magnitude = abs(integer)
        encoded = magnitude.to_bytes(max(1, (magnitude.bit_length() + 7) // 8), "big")
        digest.update(b"-" if integer < 0 else b"+")
        digest.update(len(encoded).to_bytes(8, "big"))
        digest.update(encoded)
    return "sha256:" + digest.hexdigest()


def _audit_raw_boolean_topology(vertices: object, faces: object) -> dict[str, object]:
    """Summarize raw mesh incidence without repairing or validating the solid.

    Axis-aligned boundary-plane matches are floating-point diagnostics only;
    exact plane membership requires the retained exact-coordinate sidecar.
    """
    points = np.asarray(vertices, dtype=np.float64)
    triangles = np.asarray(faces)
    if points.ndim != 2 or points.shape[1:] != (3,) or not np.isfinite(points).all():
        raise ValueError("raw Boolean topology audit requires finite Nx3 vertices")
    if triangles.ndim != 2 or triangles.shape[1:] != (3,) or \
            triangles.dtype.kind not in "iu" or np.any(triangles < 0) or \
            np.any(triangles >= len(points)):
        raise ValueError("raw Boolean topology audit requires valid integer triangles")

    incidence: dict[tuple[int, int], list[tuple[int, int]]] = {}
    repeated_vertex_face_count = 0
    for face_index, face in enumerate(triangles):
        a, b, c = (int(value) for value in face)
        if len({a, b, c}) != 3:
            repeated_vertex_face_count += 1
            continue
        for first, second in ((a, b), (b, c), (c, a)):
            edge = (min(first, second), max(first, second))
            direction = 1 if (first, second) == edge else -1
            incidence.setdefault(edge, []).append((face_index, direction))

    boundary_edges = [edge for edge, adjacent in incidence.items() if len(adjacent) == 1]
    nonmanifold_edges = [edge for edge, adjacent in incidence.items() if len(adjacent) > 2]
    same_direction_edges = [
        edge for edge, adjacent in incidence.items()
        if len(adjacent) == 2 and adjacent[0][1] == adjacent[1][1]
    ]
    face_graph: list[list[tuple[int, int, tuple[int, int]]]] = [[] for _ in range(len(triangles))]
    for edge, adjacent in incidence.items():
        if len(adjacent) != 2:
            continue
        (first_face, first_direction), (second_face, second_direction) = adjacent
        required_difference = int(first_direction == second_direction)
        face_graph[first_face].append((second_face, required_difference, edge))
        face_graph[second_face].append((first_face, required_difference, edge))
    flip_state: list[int | None] = [None] * len(triangles)
    component_count = 0
    orientable_component_count = 0
    orientation_conflicts: set[tuple[int, int]] = set()
    for root in range(len(triangles)):
        if flip_state[root] is not None:
            continue
        component_count += 1
        flip_state[root] = 0
        stack = [root]
        component_conflicted = False
        while stack:
            face_index = stack.pop()
            assert flip_state[face_index] is not None
            for neighbor, difference, edge in face_graph[face_index]:
                expected = flip_state[face_index] ^ difference
                if flip_state[neighbor] is None:
                    flip_state[neighbor] = expected
                    stack.append(neighbor)
                elif flip_state[neighbor] != expected:
                    orientation_conflicts.add(edge)
                    component_conflicted = True
        if not component_conflicted:
            orientable_component_count += 1

    boundary_adjacency: dict[int, list[int]] = {}
    for first, second in boundary_edges:
        boundary_adjacency.setdefault(first, []).append(second)
        boundary_adjacency.setdefault(second, []).append(first)
    unseen = set(boundary_adjacency)
    component_sizes: list[int] = []
    component_degrees: list[list[int]] = []
    directed_boundary_cycle_count = 0
    boundary_out: dict[int, int] = {}
    boundary_in: dict[int, int] = {}
    for edge in boundary_edges:
        face_index, direction = incidence[edge][0]
        if flip_state[face_index] == 1:
            direction *= -1
        first, second = edge if direction == 1 else (edge[1], edge[0])
        boundary_out[first] = boundary_out.get(first, 0) + 1
        boundary_in[second] = boundary_in.get(second, 0) + 1
    while unseen:
        root = min(unseen)
        stack = [root]
        component: set[int] = set()
        unseen.remove(root)
        while stack:
            vertex = stack.pop()
            component.add(vertex)
            for neighbor in boundary_adjacency[vertex]:
                if neighbor in unseen:
                    unseen.remove(neighbor)
                    stack.append(neighbor)
        component_sizes.append(len(component))
        degrees = [len(boundary_adjacency[index]) for index in component]
        component_degrees.append(degrees)
        if all(boundary_out.get(index, 0) == 1 and boundary_in.get(index, 0) == 1
               for index in component):
            directed_boundary_cycle_count += 1

    boundary_vertices = sorted(boundary_adjacency)
    plane_candidates: list[dict[str, object]] = []
    if boundary_vertices:
        scale = max(1.0, float(np.max(np.abs(points[boundary_vertices]))))
        tolerance = 512.0 * np.finfo(np.float64).eps * scale
        for axis, name in enumerate(("x", "y", "z")):
            coordinates = points[boundary_vertices, axis]
            residual = float(np.max(coordinates) - np.min(coordinates))
            if residual <= tolerance:
                plane_candidates.append({
                    "axis": name,
                    "coordinate": float(np.mean(coordinates)),
                    "maximum_coordinate_spread": residual,
                    "tolerance": tolerance,
                    "meaning": "floating-point planar candidate only; not exact coplanarity proof",
                })

    if repeated_vertex_face_count:
        status = "REPEATED_VERTEX_FACE"
    elif nonmanifold_edges:
        status = "NONMANIFOLD_EDGE_INCIDENCE"
    elif boundary_edges:
        status = "OPEN_BOUNDARY"
    else:
        status = "CLOSED_EDGE_MANIFOLD"
    return {
        "schema": "raw-boolean-topology-diagnostic-v1",
        "status": status,
        "closed_and_manifold": status == "CLOSED_EDGE_MANIFOLD",
        "vertex_count": int(len(points)),
        "face_count": int(len(triangles)),
        "edge_count": int(len(incidence)),
        "boundary_edge_count": int(len(boundary_edges)),
        "boundary_component_sizes": sorted(component_sizes),
        "boundary_loop_sizes": sorted(
            size for size, degrees in zip(component_sizes, component_degrees)
            if all(degree == 2 for degree in degrees)
        ),
        "boundary_vertices_not_degree_two": int(sum(
            degree != 2 for degrees in component_degrees for degree in degrees
        )),
        "nonmanifold_edge_count": int(len(nonmanifold_edges)),
        "same_direction_two_face_edge_count": int(len(same_direction_edges)),
        "orientation_conflict_edge_count": int(len(orientation_conflicts)),
        "orientation_consistent_as_emitted": not same_direction_edges,
        "face_component_count": int(component_count),
        "orientable_face_component_count": int(orientable_component_count),
        "suggested_face_flip_count": int(sum(value == 1 for value in flip_state)),
        "directed_boundary_cycle_count_after_consistent_reorientation": int(directed_boundary_cycle_count),
        "repeated_vertex_face_count": repeated_vertex_face_count,
        "near_axis_aligned_boundary_planes": plane_candidates,
        "repair_performed": False,
    }


def _audit_exact_boolean_boundary(
    vertices: object, faces: object, exact_sidecar: Path, face_sidecar: Path,
) -> dict[str, object]:
    """Bind open boundary cycles to exact Cherchi coordinates and face labels."""
    points = np.asarray(vertices, dtype=np.float64)
    triangles = np.asarray(faces, dtype=np.int64)

    def exact_fraction(token: str) -> Fraction:
        if token == "0/1":
            return Fraction(0)
        negative = token.startswith("-")
        numerator, denominator = token[1:].split("/") if negative else token.split("/")

        def integer(value: str) -> int:
            return int(value[2:], 2) if value.startswith("0b") else int(value)

        result = Fraction(integer(numerator), integer(denominator))
        return -result if negative else result

    exact_by_memory: dict[tuple[float, float, float], tuple[Fraction, Fraction, Fraction]] = {}
    with exact_sidecar.open("r", encoding="ascii", newline="") as stream:
        rows = csv.DictReader(stream, delimiter="\t")
        required = {"status", "exact_x", "exact_y", "exact_z", "memory_x", "memory_y", "memory_z"}
        if rows.fieldnames is None or not required.issubset(rows.fieldnames):
            raise ValueError("exact-coordinate sidecar has an unsupported header")
        for row in rows:
            if row["status"] != "MATCHED":
                continue
            memory = tuple(float(row[f"memory_{axis}"]) for axis in "xyz")
            exact = tuple(exact_fraction(row[f"exact_{axis}"]) for axis in "xyz")
            if memory in exact_by_memory and exact_by_memory[memory] != exact:
                raise ValueError("one Cherchi memory coordinate maps to multiple exact points")
            exact_by_memory[memory] = exact
    exact_vertices: list[tuple[Fraction, Fraction, Fraction]] = []
    try:
        exact_vertices = [exact_by_memory[tuple(map(float, vertex))] for vertex in points]
    except KeyError as error:
        raise ValueError("raw Boolean vertex is absent from exact-coordinate sidecar") from error

    labels: list[str] = []
    with face_sidecar.open("r", encoding="ascii", newline="") as stream:
        rows = csv.DictReader(stream, delimiter="\t")
        if rows.fieldnames != ["face_index", "input_object_mask"]:
            raise ValueError("face-label sidecar has an unsupported header")
        for expected_index, row in enumerate(rows):
            if int(row["face_index"]) != expected_index:
                raise ValueError("face-label sidecar is not in output-face order")
            labels.append(row["input_object_mask"])
    if len(labels) != len(triangles):
        raise ValueError("face-label sidecar does not align with raw Boolean faces")

    edge_faces: dict[tuple[int, int], list[int]] = {}
    for face_index, face in enumerate(triangles):
        a, b, c = (int(value) for value in face)
        for first, second in ((a, b), (b, c), (c, a)):
            edge_faces.setdefault((min(first, second), max(first, second)), []).append(face_index)
    boundary_edges = [edge for edge, adjacent in edge_faces.items() if len(adjacent) == 1]
    adjacency: dict[int, list[int]] = {}
    for first, second in boundary_edges:
        adjacency.setdefault(first, []).append(second)
        adjacency.setdefault(second, []).append(first)

    unseen = set(adjacency)
    cycles: list[dict[str, object]] = []
    while unseen:
        root = min(unseen)
        stack = [root]
        component: set[int] = set()
        unseen.remove(root)
        while stack:
            vertex = stack.pop()
            component.add(vertex)
            for neighbor in adjacency[vertex]:
                if neighbor in unseen:
                    unseen.remove(neighbor)
                    stack.append(neighbor)
        if any(len(adjacency[vertex]) != 2 for vertex in component):
            cycles.append({"vertex_count": len(component), "status": "BOUNDARY_GRAPH_NOT_CYCLE"})
            continue
        start = min(component)
        first = min(adjacency[start])
        cycle = [start, first]
        previous, current = start, first
        while True:
            following = next(vertex for vertex in adjacency[current] if vertex != previous)
            if following == start:
                break
            if following in cycle:
                raise ValueError("boundary cycle traversal repeated a vertex")
            cycle.append(following)
            previous, current = current, following
        planes: list[dict[str, object]] = []
        for axis, name in enumerate("xyz"):
            values = {exact_vertices[vertex][axis] for vertex in cycle}
            if len(values) == 1:
                value = next(iter(values))
                planes.append({"axis": name, "coordinate_exact": f"{value.numerator}/{value.denominator}"})
        if len(planes) == 1:
            normal_axis = "xyz".index(str(planes[0]["axis"]))
            projection_axes = [axis for axis in range(3) if axis != normal_axis]
            exact_simple = _audit_exact_planar_cycle([
                (exact_vertices[vertex][projection_axes[0]],
                 exact_vertices[vertex][projection_axes[1]])
                for vertex in cycle
            ])
        else:
            exact_simple = {"status": "NOT_SINGLE_AXIS_PLANAR", "intersection_count": None}
        cycle_edges = {
            (min(cycle[index], cycle[(index + 1) % len(cycle)]),
             max(cycle[index], cycle[(index + 1) % len(cycle)]))
            for index in range(len(cycle))
        }
        adjacent_faces = [edge_faces[edge][0] for edge in cycle_edges]
        label_counts: dict[str, int] = {}
        for face_index in adjacent_faces:
            label_counts[labels[face_index]] = label_counts.get(labels[face_index], 0) + 1
        cycles.append({
            "vertex_count": len(cycle),
            "status": exact_simple["status"],
            "simple_planar_polygon_proven": exact_simple["status"] == "CERTIFIED_SIMPLE_POLYGON",
            "exact_polygon_simplicity": exact_simple,
            "exact_axis_aligned_planes": planes,
            "adjacent_face_input_object_mask_counts": label_counts,
        })
    return {
        "schema": "cherchi-exact-open-boundary-provenance-v1",
        "boundary_edge_count": len(boundary_edges),
        "components": cycles,
        "repair_performed": False,
    }


def _audit_exact_planar_cycle(
    points: list[tuple[Fraction, Fraction]], *, maximum_candidate_pairs: int = 1_000_000,
) -> dict[str, object]:
    """Certify simplicity of a rational polygon using exact segment predicates.

    An exact rational uniform-grid broad phase generates a superset of all
    intersecting segment pairs. Orientation, collinearity, and contact tests
    are then exact. This certifies only polygon simplicity, not which planar
    region belongs to a Boolean result.
    """
    count = len(points)
    if count < 3 or any(len(point) != 2 for point in points):
        raise ValueError("exact planar cycle requires at least three 2D rational points")
    if len(set(points)) != count:
        return {"status": "SELF_INTERSECTION_OR_TOUCH", "intersection_count": 1,
                "witness": "repeated_nonadjacent_vertex"}
    low = [min(point[axis] for point in points) for axis in range(2)]
    high = [max(point[axis] for point in points) for axis in range(2)]
    if any(low[axis] == high[axis] for axis in range(2)):
        return {"status": "DEGENERATE_PROJECTION", "intersection_count": None}
    nx = ny = min(128, max(16, int(np.ceil(np.sqrt(count)))))

    def cell(value: Fraction, axis: int) -> int:
        bins = (nx, ny)[axis]
        scaled = (value - low[axis]) * bins / (high[axis] - low[axis])
        return min(bins - 1, max(0, scaled.numerator // scaled.denominator))

    buckets: dict[tuple[int, int], list[int]] = {}
    for index, first in enumerate(points):
        second = points[(index + 1) % count]
        ix0, ix1 = cell(min(first[0], second[0]), 0), cell(max(first[0], second[0]), 0)
        iy0, iy1 = cell(min(first[1], second[1]), 1), cell(max(first[1], second[1]), 1)
        for ix in range(ix0, ix1 + 1):
            for iy in range(iy0, iy1 + 1):
                buckets.setdefault((ix, iy), []).append(index)
    candidates: set[tuple[int, int]] = set()
    for indices in buckets.values():
        ordered = sorted(set(indices))
        for offset, first in enumerate(ordered):
            for second in ordered[offset + 1:]:
                if second == first + 1 or (first == 0 and second == count - 1):
                    continue
                candidates.add((first, second))
                if len(candidates) > maximum_candidate_pairs:
                    return {"status": "PAIR_BUDGET_EXCEEDED", "candidate_pair_count": len(candidates),
                            "maximum_candidate_pairs": maximum_candidate_pairs,
                            "intersection_count": None}

    def orient(a, b, c) -> Fraction:
        return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])

    def on_segment(a, b, point) -> bool:
        return (min(a[0], b[0]) <= point[0] <= max(a[0], b[0]) and
                min(a[1], b[1]) <= point[1] <= max(a[1], b[1]))

    def intersects(first_index: int, second_index: int) -> bool:
        a, b = points[first_index], points[(first_index + 1) % count]
        c, d = points[second_index], points[(second_index + 1) % count]
        o1, o2, o3, o4 = orient(a, b, c), orient(a, b, d), orient(c, d, a), orient(c, d, b)
        if ((o1 > 0 > o2 or o2 > 0 > o1) and (o3 > 0 > o4 or o4 > 0 > o3)):
            return True
        return ((o1 == 0 and on_segment(a, b, c)) or
                (o2 == 0 and on_segment(a, b, d)) or
                (o3 == 0 and on_segment(c, d, a)) or
                (o4 == 0 and on_segment(c, d, b)))

    witnesses = [(first, second) for first, second in sorted(candidates)
                 if intersects(first, second)]
    if not witnesses:
        # Adjacent collinear edges that reverse direction overlap beyond their
        # shared endpoint and are not part of a simple polygon boundary.
        for index in range(count):
            previous, current, following = points[index - 1], points[index], points[(index + 1) % count]
            if orient(previous, current, following) == 0:
                left = (previous[0] - current[0], previous[1] - current[1])
                right = (following[0] - current[0], following[1] - current[1])
                if left[0] * right[0] + left[1] * right[1] > 0:
                    witnesses = [(index - 1, index)]
                    break
    return {
        "status": "SELF_INTERSECTION_OR_TOUCH" if witnesses else "CERTIFIED_SIMPLE_POLYGON",
        "intersection_count": len(witnesses),
        "candidate_pair_count": len(candidates),
        "witness_segment_pairs": [list(pair) for pair in witnesses[:8]],
        "predicate": "exact rational orientation and segment contact",
    }


def _exact_point_in_polygon(
    polygon: list[tuple[Fraction, Fraction]], query: tuple[Fraction, Fraction],
) -> str:
    """Classify a rational 2D point against a simple rational polygon."""
    if len(polygon) < 3 or len(query) != 2:
        raise ValueError("exact point-in-polygon requires a polygon and a 2D query")
    inside = False
    for index, first in enumerate(polygon):
        second = polygon[(index + 1) % len(polygon)]
        orientation = (second[0] - first[0]) * (query[1] - first[1]) - \
            (second[1] - first[1]) * (query[0] - first[0])
        if orientation == 0 and min(first[0], second[0]) <= query[0] <= max(first[0], second[0]) and \
                min(first[1], second[1]) <= query[1] <= max(first[1], second[1]):
            return "ON_BOUNDARY"
        if (first[1] > query[1]) != (second[1] > query[1]):
            x_at_query_y = first[0] + (query[1] - first[1]) * \
                (second[0] - first[0]) / (second[1] - first[1])
            if x_at_query_y > query[0]:
                inside = not inside
    return "INSIDE" if inside else "OUTSIDE"


def _exact_mesh_plane_section_cycles(
    vertices: object, faces: object, *, axis: str, plane: Fraction,
) -> dict[str, object]:
    """Intersect represented binary64 triangles with an exact rational plane.

    Returns all closed section cycles only for the generic case where each
    intersected triangle contributes one segment and every section vertex has
    degree two. Coplanar facets, tangencies, and non-cycles are reported as
    UNKNOWN rather than regularized.
    """
    if axis not in {"x", "y", "z"} or not isinstance(plane, Fraction):
        raise ValueError("exact plane section requires axis x/y/z and a Fraction plane")
    points = np.asarray(vertices, dtype=np.float64)
    triangles = np.asarray(faces)
    axis_index = {"x": 0, "y": 1, "z": 2}[axis]
    if points.ndim != 2 or points.shape[1:] != (3,) or not np.isfinite(points).all() or \
            triangles.ndim != 2 or triangles.shape[1:] != (3,) or \
            triangles.dtype.kind not in "iu" or np.any(triangles < 0) or \
            np.any(triangles >= len(points)):
        return {"status": "UNKNOWN_INVALID_MESH", "cycle_count": 0, "cycles": []}
    exact_points = [tuple(Fraction.from_float(float(value)) for value in point) for point in points]
    segment_set: set[tuple[tuple[Fraction, Fraction, Fraction],
                           tuple[Fraction, Fraction, Fraction]]] = set()
    tangent_points: set[tuple[Fraction, Fraction, Fraction]] = set()
    coplanar_face_count = 0
    for face in triangles:
        triangle = [exact_points[int(index)] for index in face]
        distances = [point[axis_index] - plane for point in triangle]
        if all(distance == 0 for distance in distances):
            coplanar_face_count += 1
            continue
        intersections: set[tuple[Fraction, Fraction, Fraction]] = set()
        for edge in range(3):
            first = triangle[edge]
            second = triangle[(edge + 1) % 3]
            d_first = distances[edge]
            d_second = distances[(edge + 1) % 3]
            if d_first == 0:
                intersections.add(first)
            if d_first * d_second < 0:
                parameter = d_first / (d_first - d_second)
                intersections.add(tuple(
                    first[coordinate] + parameter * (second[coordinate] - first[coordinate])
                    for coordinate in range(3)
                ))
        if len(intersections) == 2:
            first, second = sorted(intersections)
            segment_set.add((first, second))
        elif len(intersections) == 1:
            tangent_points.update(intersections)
        elif len(intersections) > 2:
            return {
                "status": "UNKNOWN_AMBIGUOUS_FACE_SECTION",
                "cycle_count": 0,
                "cycles": [],
                "source_face": int(len(segment_set)),
            }
    if coplanar_face_count:
        return {
            "status": "UNKNOWN_COPLANAR_SOURCE_FACES",
            "cycle_count": 0,
            "cycles": [],
            "coplanar_face_count": coplanar_face_count,
        }
    adjacency: dict[tuple[Fraction, Fraction, Fraction], set[tuple[Fraction, Fraction, Fraction]]] = {}
    for first, second in segment_set:
        adjacency.setdefault(first, set()).add(second)
        adjacency.setdefault(second, set()).add(first)
    if any(len(neighbors) != 2 for neighbors in adjacency.values()) or \
            any(point not in adjacency for point in tangent_points):
        return {
            "status": "UNKNOWN_NON_CYCLE_SECTION",
            "cycle_count": 0,
            "cycles": [],
            "section_vertex_count": len(adjacency),
            "segment_count": len(segment_set),
            "non_two_valent_vertex_count": sum(len(neighbors) != 2 for neighbors in adjacency.values()),
            "isolated_tangency_count": sum(point not in adjacency for point in tangent_points),
        }
    cycles: list[list[tuple[Fraction, Fraction, Fraction]]] = []
    unseen = set(adjacency)
    while unseen:
        start = min(unseen)
        cycle = [start]
        previous = None
        current = start
        while True:
            candidates = sorted(neighbor for neighbor in adjacency[current] if neighbor != previous)
            if not candidates:
                return {"status": "UNKNOWN_NON_CYCLE_SECTION", "cycle_count": 0, "cycles": []}
            following = candidates[0]
            if following == start:
                break
            if following in cycle:
                return {"status": "UNKNOWN_NON_CYCLE_SECTION", "cycle_count": 0, "cycles": []}
            cycle.append(following)
            previous, current = current, following
        cycles.append(cycle)
        unseen.difference_update(cycle)
    return {
        "status": "CLOSED_SECTION_CYCLES",
        "cycle_count": len(cycles),
        "cycles": cycles,
        "section_vertex_count": len(adjacency),
        "segment_count": len(segment_set),
        "predicate": "exact rational plane/triangle edge intersection on binary64 input coordinates",
    }


def _exact_cycle_matches_up_to_collinear_subdivision(
    candidate: list[tuple[Fraction, Fraction, Fraction]],
    reference: list[tuple[Fraction, Fraction, Fraction]], *, axis: str,
) -> dict[str, object]:
    """Prove two planar cycles have identical edges up to exact subdivisions."""
    if axis not in {"x", "y", "z"}:
        raise ValueError("cycle comparison axis must be x, y, or z")
    if len(candidate) < 3 or len(reference) < 3:
        return {"matches": False, "diagnostic": "cycle has fewer than three vertices"}
    projection = {"x": (1, 2), "y": (0, 2), "z": (0, 1)}[axis]
    candidate_2d = [(point[projection[0]], point[projection[1]]) for point in candidate]
    reference_2d = [(point[projection[0]], point[projection[1]]) for point in reference]

    def cross(a, b, p):
        return (b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0])

    def parameter(a, b, p):
        delta = b[0] - a[0]
        if delta != 0:
            return (p[0] - a[0]) / delta
        delta = b[1] - a[1]
        return (p[1] - a[1]) / delta

    output_edges = list(zip(candidate_2d, candidate_2d[1:] + candidate_2d[:1]))
    reference_edges = list(zip(reference_2d, reference_2d[1:] + reference_2d[:1]))
    intervals: list[list[tuple[Fraction, Fraction]]] = [
        [] for _ in reference_edges
    ]
    edge_to_reference: list[int] = []
    for first, second in output_edges:
        matches = []
        for index, (ref_first, ref_second) in enumerate(reference_edges):
            if cross(ref_first, ref_second, first) != 0 or cross(ref_first, ref_second, second) != 0:
                continue
            if all(min(ref_first[dimension], ref_second[dimension]) <= point[dimension] <=
                   max(ref_first[dimension], ref_second[dimension]) for point in (first, second)
                   for dimension in range(2)):
                matches.append(index)
        if len(matches) != 1:
            return {
                "matches": False,
                "diagnostic": "candidate edge is not contained in exactly one reference edge",
                "unmatched_or_ambiguous_candidate_edge": len(edge_to_reference),
                "matching_reference_edge_count": len(matches),
            }
        reference_index = matches[0]
        ref_first, ref_second = reference_edges[reference_index]
        parameters = sorted((parameter(ref_first, ref_second, first),
                             parameter(ref_first, ref_second, second)))
        intervals[reference_index].append((parameters[0], parameters[1]))
        edge_to_reference.append(reference_index)

    for index, covered in enumerate(intervals):
        covered.sort()
        if not covered or covered[0][0] != 0 or covered[-1][1] != 1 or \
                any(covered[j][1] != covered[j + 1][0] for j in range(len(covered) - 1)):
            return {
                "matches": False,
                "diagnostic": "reference edge is not exactly covered by candidate subedges",
                "uncovered_reference_edge": index,
                "candidate_subedge_count": len(covered),
            }
    reference_vertices = set(reference)
    extra_vertices = sum(point not in reference_vertices for point in candidate)
    return {
        "matches": True,
        "diagnostic": "all cycle edges match exactly after collinear subdivision",
        "reference_segment_count": len(reference_edges),
        "candidate_segment_count": len(output_edges),
        "extra_collinear_vertices": extra_vertices,
        "candidate_edge_reference_mapping_sha256": hashlib.sha256(
            ",".join(str(index) for index in edge_to_reference).encode("ascii")
        ).hexdigest(),
        "predicate": "exact rational collinearity, containment, and interval coverage",
    }


def _triangulate_planar_loop_from_kernel_point(
    loop: list[tuple[Fraction, Fraction, Fraction]], *, axis: str,
    plane: Fraction, desired_normal_sign: int,
) -> dict[str, object]:
    """Build a fan only when the rounded center is exactly in the polygon kernel."""
    if axis not in {"x", "y", "z"} or not isinstance(plane, Fraction) or \
            desired_normal_sign not in {-1, 1}:
        raise ValueError("fan requires an axis plane and desired normal sign +/-1")
    if len(loop) < 3 or any(point[{"x": 0, "y": 1, "z": 2}[axis]] != plane for point in loop):
        return {"status": "UNKNOWN_INVALID_PLANAR_LOOP", "triangles": []}
    projection = {"x": (1, 2), "y": (0, 2), "z": (0, 1)}[axis]
    factor = {"x": 1, "y": -1, "z": 1}[axis]
    polygon = [(point[projection[0]], point[projection[1]]) for point in loop]
    simplicity = _audit_exact_planar_cycle(polygon)
    if simplicity["status"] != "CERTIFIED_SIMPLE_POLYGON":
        return {"status": "UNKNOWN_NON_SIMPLE_LOOP", "triangles": [], "simplicity": simplicity}
    crosses = [
        first[0] * second[1] - second[0] * first[1]
        for first, second in zip(polygon, polygon[1:] + polygon[:1])
    ]
    twice_area = sum(crosses)
    if twice_area == 0:
        return {"status": "UNKNOWN_ZERO_AREA_LOOP", "triangles": []}
    center_2d = (
        sum((first[0] + second[0]) * cross
            for (first, second), cross in zip(zip(polygon, polygon[1:] + polygon[:1]), crosses)) /
        (3 * twice_area),
        sum((first[1] + second[1]) * cross
            for (first, second), cross in zip(zip(polygon, polygon[1:] + polygon[:1]), crosses)) /
        (3 * twice_area),
    )
    center_exact = [Fraction(0), Fraction(0), Fraction(0)]
    center_exact[{"x": 0, "y": 1, "z": 2}[axis]] = plane
    center_exact[projection[0]], center_exact[projection[1]] = center_2d
    center_binary64 = tuple(float(value) for value in center_exact)
    center_represented = tuple(Fraction.from_float(value) for value in center_binary64)
    center_represented_2d = (center_represented[projection[0]], center_represented[projection[1]])
    area_sign = 1 if twice_area > 0 else -1
    margins = [
        ((second[0] - first[0]) * (center_represented_2d[1] - first[1]) -
         (second[1] - first[1]) * (center_represented_2d[0] - first[0])) * area_sign
        for first, second in zip(polygon, polygon[1:] + polygon[:1])
    ]
    if any(margin <= 0 for margin in margins):
        return {
            "status": "UNKNOWN_ROUNDED_CENTER_OUTSIDE_KERNEL",
            "triangles": [],
            "kernel_halfplane_count": len(margins),
            "strict_kernel_halfplane_count": sum(margin > 0 for margin in margins),
        }
    if _exact_point_in_polygon(polygon, center_represented_2d) != "INSIDE":
        return {"status": "UNKNOWN_ROUNDED_CENTER_NOT_INTERIOR", "triangles": []}
    required_area_sign = desired_normal_sign * factor
    order = list(range(len(loop)))
    if area_sign != required_area_sign:
        order.reverse()
    triangles = [
        (len(loop), order[index], order[(index + 1) % len(loop)])
        for index in range(len(loop))
    ]
    exact_center_digest = hashlib.sha256()
    center_integer_bit_lengths = []
    for value in center_exact:
        for integer in (value.numerator, value.denominator):
            magnitude = abs(integer)
            encoded = magnitude.to_bytes(max(1, (magnitude.bit_length() + 7) // 8), "big")
            exact_center_digest.update(b"-" if integer < 0 else b"+")
            exact_center_digest.update(len(encoded).to_bytes(8, "big"))
            exact_center_digest.update(encoded)
            center_integer_bit_lengths.append(integer.bit_length())
    return {
        "status": "CERTIFIED_KERNEL_FAN",
        "triangles": triangles,
        "center_exact_sha256": exact_center_digest.hexdigest(),
        "center_exact_integer_bit_lengths": center_integer_bit_lengths,
        "center_binary64": list(center_binary64),
        "kernel_halfplane_count": len(margins),
        "strict_kernel_halfplane_count": sum(margin > 0 for margin in margins),
        "minimum_oriented_kernel_margin": float(min(margins)),
        "desired_normal_sign": desired_normal_sign,
        "projected_polygon_area_sign": area_sign,
        "projection_normal_factor": factor,
        "predicate": "exact rational kernel halfplanes after binary64 center rounding",
    }


def _triangulate_exact_star_shaped_polygon_difference(
    outer: list[tuple[Fraction, Fraction]],
    inner: list[tuple[Fraction, Fraction]], *,
    kernel_point: tuple[Fraction, Fraction],
) -> dict[str, object]:
    """Triangulate an exact planar annulus by a common-kernel radial sweep.

    Every ray through an input boundary vertex is an event.  If both simple
    polygons are strictly star-shaped about the supplied exact point and the
    inner loop is strictly inside the outer loop, each ray has one boundary
    hit per loop.  The intervening radial quadrilaterals partition the region;
    exact boundary-incidence and area checks certify the emitted triangles.
    """
    empty = {"status": "UNKNOWN_INVALID_POLYGON", "vertices": [], "triangles": []}
    if len(outer) < 3 or len(inner) < 3 or len(kernel_point) != 2 or \
            any(len(point) != 2 or not all(isinstance(value, Fraction) for value in point)
                for point in [*outer, *inner, kernel_point]):
        return empty

    def cross(a, b):
        return a[0] * b[1] - a[1] * b[0]

    def orient(a, b, c):
        return cross((b[0] - a[0], b[1] - a[1]),
                     (c[0] - a[0], c[1] - a[1]))

    def signed_area2(polygon):
        return sum(cross(first, second)
                   for first, second in zip(polygon, polygon[1:] + polygon[:1]))

    outer_audit = _audit_exact_planar_cycle(outer)
    inner_audit = _audit_exact_planar_cycle(inner)
    if outer_audit["status"] != "CERTIFIED_SIMPLE_POLYGON" or \
            inner_audit["status"] != "CERTIFIED_SIMPLE_POLYGON":
        return {**empty, "status": "UNKNOWN_NON_SIMPLE_INPUT",
                "outer_audit": outer_audit, "inner_audit": inner_audit}
    outer_sign = 1 if signed_area2(outer) > 0 else -1
    inner_sign = 1 if signed_area2(inner) > 0 else -1
    if signed_area2(outer) == 0 or signed_area2(inner) == 0:
        return {**empty, "status": "UNKNOWN_ZERO_AREA_INPUT"}
    if _exact_point_in_polygon(outer, kernel_point) != "INSIDE" or \
            _exact_point_in_polygon(inner, kernel_point) != "INSIDE":
        return {**empty, "status": "UNKNOWN_KERNEL_POINT_NOT_INTERIOR"}
    if any(_exact_point_in_polygon(outer, point) != "INSIDE" for point in inner):
        return {**empty, "status": "UNKNOWN_INNER_NOT_STRICTLY_INSIDE_OUTER"}
    for polygon, winding in ((outer, outer_sign), (inner, inner_sign)):
        for first, second in zip(polygon, polygon[1:] + polygon[:1]):
            if orient(first, second, kernel_point) * winding <= 0:
                return {**empty, "status": "UNKNOWN_KERNEL_POINT_OUTSIDE_POLYGON_KERNEL"}

    directions = []
    for point in (*outer, *inner):
        direction = (point[0] - kernel_point[0], point[1] - kernel_point[1])
        if direction == (0, 0):
            return {**empty, "status": "UNKNOWN_KERNEL_COINCIDES_WITH_VERTEX"}
        directions.append(direction)

    def half(direction):
        return 0 if direction[1] > 0 or (direction[1] == 0 and direction[0] >= 0) else 1

    def compare_directions(first, second):
        first_half, second_half = half(first), half(second)
        if first_half != second_half:
            return -1 if first_half < second_half else 1
        determinant = cross(first, second)
        return -1 if determinant > 0 else 1 if determinant < 0 else 0

    directions.sort(key=cmp_to_key(compare_directions))
    unique_directions = []
    for direction in directions:
        if unique_directions and cross(unique_directions[-1], direction) == 0 and \
                unique_directions[-1][0] * direction[0] + unique_directions[-1][1] * direction[1] > 0:
            continue
        unique_directions.append(direction)
    if len(unique_directions) < 3:
        return {**empty, "status": "UNKNOWN_TOO_FEW_RADIAL_EVENTS"}

    def ray_hit(polygon, direction):
        hits: dict[tuple[Fraction, Fraction], Fraction] = {}
        for first, second in zip(polygon, polygon[1:] + polygon[:1]):
            edge = (second[0] - first[0], second[1] - first[1])
            offset = (first[0] - kernel_point[0], first[1] - kernel_point[1])
            denominator = cross(direction, edge)
            if denominator == 0:
                if cross(offset, direction) == 0:
                    return None, "UNKNOWN_RAY_COLLINEAR_WITH_EDGE"
                continue
            ray_parameter = cross(offset, edge) / denominator
            edge_parameter = cross(offset, direction) / denominator
            if ray_parameter > 0 and 0 <= edge_parameter <= 1:
                point = (kernel_point[0] + ray_parameter * direction[0],
                         kernel_point[1] + ray_parameter * direction[1])
                hits[point] = ray_parameter
        if len(hits) != 1:
            return None, "UNKNOWN_RAY_HAS_NONUNIQUE_POLYGON_HIT"
        return next(iter(hits.items()))[0], "OK"

    outer_hits = []
    inner_hits = []
    for direction in unique_directions:
        outer_point, outer_status = ray_hit(outer, direction)
        inner_point, inner_status = ray_hit(inner, direction)
        if outer_status != "OK" or inner_status != "OK":
            return {**empty, "status": outer_status if outer_status != "OK" else inner_status,
                    "radial_event_count": len(unique_directions)}
        component = 0 if direction[0] != 0 else 1
        outer_parameter = (outer_point[component] - kernel_point[component]) / direction[component]
        inner_parameter = (inner_point[component] - kernel_point[component]) / direction[component]
        if inner_parameter <= 0 or outer_parameter <= inner_parameter:
            return {**empty, "status": "UNKNOWN_RADIAL_BOUNDARIES_NOT_NESTED"}
        outer_hits.append(outer_point)
        inner_hits.append(inner_point)

    vertex_list: list[tuple[Fraction, Fraction]] = []
    vertex_index: dict[tuple[Fraction, Fraction], int] = {}

    def intern(point):
        if point not in vertex_index:
            vertex_index[point] = len(vertex_list)
            vertex_list.append(point)
        return vertex_index[point]

    triangles: list[tuple[int, int, int]] = []
    triangle_area2_sum = Fraction(0)
    for index in range(len(unique_directions)):
        following = (index + 1) % len(unique_directions)
        oi, oj = intern(outer_hits[index]), intern(outer_hits[following])
        ii, ij = intern(inner_hits[index]), intern(inner_hits[following])
        wedge = ((oi, oj, ij), (oi, ij, ii))
        for triangle in wedge:
            a, b, c = (vertex_list[vertex] for vertex in triangle)
            area2 = orient(a, b, c)
            if area2 <= 0:
                return {**empty, "status": "UNKNOWN_RADIAL_WEDGE_NOT_POSITIVE"}
            centroid = ((a[0] + b[0] + c[0]) / 3, (a[1] + b[1] + c[1]) / 3)
            if _exact_point_in_polygon(outer, centroid) != "INSIDE" or \
                    _exact_point_in_polygon(inner, centroid) != "OUTSIDE":
                return {**empty, "status": "UNKNOWN_TRIANGLE_OUTSIDE_ANNULUS"}
            triangles.append(triangle)
            triangle_area2_sum += area2

    edge_counts: dict[tuple[int, int], int] = {}
    for a, b, c in triangles:
        for first, second in ((a, b), (b, c), (c, a)):
            edge = (min(first, second), max(first, second))
            edge_counts[edge] = edge_counts.get(edge, 0) + 1
    if any(count not in {1, 2} for count in edge_counts.values()):
        return {**empty, "status": "UNKNOWN_TRIANGLE_EDGE_INCIDENCE"}
    actual_boundary = {edge for edge, count in edge_counts.items() if count == 1}
    expected_boundary = set()
    for hits in (outer_hits, inner_hits):
        for index, first in enumerate(hits):
            second = hits[(index + 1) % len(hits)]
            expected_boundary.add(tuple(sorted((vertex_index[first], vertex_index[second]))))
    if actual_boundary != expected_boundary:
        return {**empty, "status": "UNKNOWN_BOUNDARY_EDGE_COVERAGE",
                "actual_boundary_edge_count": len(actual_boundary),
                "expected_boundary_edge_count": len(expected_boundary)}
    expected_area2 = abs(signed_area2(outer)) - abs(signed_area2(inner))
    if triangle_area2_sum != expected_area2:
        return {**empty, "status": "UNKNOWN_EXACT_AREA_MISMATCH",
                "triangle_area2": _compact_exact_fraction(triangle_area2_sum),
                "expected_area2": _compact_exact_fraction(expected_area2)}
    return {
        "status": "CERTIFIED_RADIAL_ANNULUS",
        "vertices": vertex_list,
        "triangles": triangles,
        "triangle_count": len(triangles),
        "boundary_edge_count": len(actual_boundary),
        "exact_area_encoding": _compact_exact_fraction(triangle_area2_sum / 2),
        "exact_area_float": float(triangle_area2_sum / 2),
        "radial_event_count": len(unique_directions),
        "proof_checks": {
            "both_loops_simple": True,
            "common_point_strictly_in_both_kernels": True,
            "inner_loop_strictly_inside_outer": True,
            "each_radial_boundary_hit_unique": True,
            "all_triangle_centroids_in_annulus": True,
            "exact_boundary_edge_set_matches": True,
            "exact_area_matches_polygon_difference": True,
        },
        "predicate": "exact rational radial sector partition with boundary and area certificates",
    }


def _build_diagnostic_cap_candidate(
    vertices: object, faces: object, *, patch_audit: dict[str, object],
    source_self_intersection_status: str,
    density: float, root_side: float,
    exact_boundary_coordinates: dict[int, tuple[Fraction, Fraction, Fraction]] | None = None,
) -> dict[str, object]:
    """Close one proven planar boundary for diagnostics, never as a Boolean proof.

    The cap is admitted only when the source is separately certified embedded,
    the patch audit identifies one exact interior disk, and the raw output has
    exactly one consistently orientable boundary cycle on that plane.  The
    completed candidate is checked by the exact binary64-boundary audit
    and by an independent signed-tetrahedron mass-moment audit using scales
    frozen from the input.  This still does not prove that the upstream
    Boolean selected the globally correct material region.
    """
    points = np.asarray(vertices, dtype=np.float64)
    triangles = np.asarray(faces, dtype=np.int64)
    if points.ndim != 2 or points.shape[1:] != (3,) or not np.isfinite(points).all() or \
            triangles.ndim != 2 or triangles.shape[1:] != (3,) or \
            np.any(triangles < 0) or np.any(triangles >= len(points)):
        raise ValueError("diagnostic cap requires finite vertices and valid triangles")
    if not np.isfinite(density) or density <= 0.0 or not np.isfinite(root_side) or root_side <= 0.0:
        raise ValueError("diagnostic cap requires positive input density and frozen root side")
    candidate_base_points = points.copy()
    if source_self_intersection_status != "INTERSECTION_FREE":
        return {"status": "UNKNOWN_SOURCE_SELF_INTERSECTION", "vertices": None, "faces": None}
    if patch_audit.get("selected_patch_candidate") == "FACE_COMPLEMENT":
        return {
            "status": "UNKNOWN_COMPLEMENT_PATCH_REQUIRES_HOLE_TRIANGULATION",
            "vertices": None,
            "faces": None,
        }
    if patch_audit.get("status") != "WITNESS_SUPPORTED_INTERIOR_DISK" or \
            patch_audit.get("selected_patch_candidate") != "INTERIOR_DISK" or \
            patch_audit.get("output_boundary_matches_unique_source_section") is not True:
        return {"status": "UNKNOWN_PATCH_NOT_PROVEN", "vertices": None, "faces": None}
    axis = patch_audit.get("axis")
    if axis not in {"x", "y", "z"}:
        return {"status": "UNKNOWN_PATCH_AXIS", "vertices": None, "faces": None}
    plane_token = str(patch_audit.get("plane_coordinate_exact", ""))
    try:
        plane = Fraction(plane_token)
    except (ValueError, ZeroDivisionError):
        return {"status": "UNKNOWN_PATCH_PLANE", "vertices": None, "faces": None}

    incidence: dict[tuple[int, int], list[tuple[int, int]]] = {}
    for face_index, face in enumerate(triangles):
        a, b, c = (int(value) for value in face)
        if len({a, b, c}) != 3:
            return {"status": "UNKNOWN_DEGENERATE_FACE", "vertices": None, "faces": None}
        for first, second in ((a, b), (b, c), (c, a)):
            edge = (min(first, second), max(first, second))
            incidence.setdefault(edge, []).append((face_index, 1 if (first, second) == edge else -1))
    if any(len(adjacent) > 2 for adjacent in incidence.values()):
        return {"status": "UNKNOWN_NONMANIFOLD_EDGE", "vertices": None, "faces": None}

    face_graph: list[list[tuple[int, int]]] = [[] for _ in range(len(triangles))]
    for adjacent in incidence.values():
        if len(adjacent) != 2:
            continue
        (first_face, first_direction), (second_face, second_direction) = adjacent
        flip_difference = int(first_direction == second_direction)
        face_graph[first_face].append((second_face, flip_difference))
        face_graph[second_face].append((first_face, flip_difference))
    flips: list[int | None] = [None] * len(triangles)
    component_count = 0
    for root in range(len(triangles)):
        if flips[root] is not None:
            continue
        component_count += 1
        flips[root] = 0
        stack = [root]
        while stack:
            current = stack.pop()
            assert flips[current] is not None
            for neighbor, difference in face_graph[current]:
                expected = flips[current] ^ difference
                if flips[neighbor] is None:
                    flips[neighbor] = expected
                    stack.append(neighbor)
                elif flips[neighbor] != expected:
                    return {"status": "UNKNOWN_NONORIENTABLE_RAW_MESH", "vertices": None, "faces": None}
    oriented = triangles.copy()
    for face_index, flip in enumerate(flips):
        if flip:
            oriented[face_index] = oriented[face_index, ::-1]

    directed_boundary: list[tuple[int, int]] = []
    for edge, adjacent in incidence.items():
        if len(adjacent) != 1:
            continue
        face_index, direction = adjacent[0]
        if flips[face_index]:
            direction *= -1
        directed_boundary.append(edge if direction == 1 else (edge[1], edge[0]))
    outgoing: dict[int, list[int]] = {}
    incoming: dict[int, int] = {}
    for first, second in directed_boundary:
        outgoing.setdefault(first, []).append(second)
        incoming[second] = incoming.get(second, 0) + 1
    boundary_vertices = set(outgoing) | set(incoming)
    if not boundary_vertices or any(len(outgoing.get(index, [])) != 1 or incoming.get(index) != 1
                                    for index in boundary_vertices):
        return {"status": "UNKNOWN_BOUNDARY_NOT_DIRECTED_CYCLES", "vertices": None, "faces": None}
    if len(directed_boundary) != len(boundary_vertices):
        return {"status": "UNKNOWN_BOUNDARY_MULTIPLE_CYCLES", "vertices": None, "faces": None}
    start = min(boundary_vertices)
    cycle = [start]
    current = start
    while True:
        current = outgoing[current][0]
        if current == start:
            break
        if current in cycle:
            return {"status": "UNKNOWN_BOUNDARY_CYCLE_REPEATS", "vertices": None, "faces": None}
        cycle.append(current)
    if len(cycle) != len(boundary_vertices):
        return {"status": "UNKNOWN_BOUNDARY_MULTIPLE_CYCLES", "vertices": None, "faces": None}

    axis_index = {"x": 0, "y": 1, "z": 2}[axis]
    replacement_count = 0
    maximum_displacement = 0.0
    if exact_boundary_coordinates is not None:
        if set(exact_boundary_coordinates) != boundary_vertices:
            return {"status": "UNKNOWN_EXACT_BOUNDARY_COORDINATE_COVERAGE", "vertices": None, "faces": None}
        for index, exact_point in exact_boundary_coordinates.items():
            if len(exact_point) != 3 or any(not isinstance(value, Fraction) for value in exact_point):
                return {"status": "UNKNOWN_INVALID_EXACT_BOUNDARY_COORDINATE", "vertices": None, "faces": None}
            replacement = np.asarray([float(value) for value in exact_point], dtype=np.float64)
            maximum_displacement = max(
                maximum_displacement, float(np.linalg.norm(replacement - candidate_base_points[index])),
            )
            if not np.array_equal(replacement, candidate_base_points[index]):
                replacement_count += 1
            candidate_base_points[index] = replacement
    exact_cycle = [
        exact_boundary_coordinates[index] if exact_boundary_coordinates is not None else
        tuple(Fraction.from_float(float(value)) for value in candidate_base_points[index])
        for index in cycle
    ]
    if any(point[axis_index] != plane for point in exact_cycle):
        return {"status": "UNKNOWN_BOUNDARY_NOT_ON_PROVEN_PLANE", "vertices": None, "faces": None}
    cap_cycle_indices = [cycle[0], *reversed(cycle[1:])]
    exact_by_vertex = dict(zip(cycle, exact_cycle))
    cap_cycle = [exact_by_vertex[index] for index in cap_cycle_indices]
    polygon_axes = {"x": (1, 2), "y": (0, 2), "z": (0, 1)}[axis]
    projected = [(point[polygon_axes[0]], point[polygon_axes[1]]) for point in cap_cycle]
    twice_area = sum(
        first[0] * second[1] - second[0] * first[1]
        for first, second in zip(projected, projected[1:] + projected[:1])
    )
    if twice_area == 0:
        return {"status": "UNKNOWN_ZERO_AREA_BOUNDARY", "vertices": None, "faces": None}
    projection_normal_factor = {"x": 1, "y": -1, "z": 1}[axis]
    fan = _triangulate_planar_loop_from_kernel_point(
        cap_cycle, axis=axis, plane=plane,
        desired_normal_sign=(1 if twice_area > 0 else -1) * projection_normal_factor,
    )
    if fan["status"] != "CERTIFIED_KERNEL_FAN":
        return {"status": fan["status"], "vertices": None, "faces": None,
                "fan_audit": fan}
    center_index = len(points)
    candidate_vertices = np.vstack((candidate_base_points, np.asarray(fan["center_binary64"], dtype=np.float64)))
    cap_faces = np.asarray([
        (center_index, cap_cycle_indices[first], cap_cycle_indices[second])
        for _, first, second in fan["triangles"]
    ], dtype=np.int64)
    candidate_faces = np.vstack((oriented, cap_faces))
    signed_volume = float(np.sum(
        np.einsum("ij,ij->i", candidate_vertices[candidate_faces[:, 0]],
                  np.cross(candidate_vertices[candidate_faces[:, 1]], candidate_vertices[candidate_faces[:, 2]]))
    ) / 6.0)
    if not math.isfinite(signed_volume) or abs(signed_volume) <= 1e-15:
        return {"status": "UNKNOWN_DEGENERATE_CANDIDATE_VOLUME", "vertices": None, "faces": None}
    if signed_volume < 0:
        candidate_faces = candidate_faces[:, ::-1]
        signed_volume = -signed_volume
    topology = _audit_raw_boolean_topology(candidate_vertices, candidate_faces)
    if topology["status"] != "CLOSED_EDGE_MANIFOLD" or topology["orientation_conflict_edge_count"] != 0 or \
            topology["face_component_count"] != 1:
        return {"status": "UNKNOWN_CANDIDATE_TOPOLOGY", "vertices": None, "faces": None,
                "topology_audit": topology}
    try:
        exact_audit = spin._native_module().audit_exact_triangle_boundary(
            candidate_vertices, candidate_faces, "OutwardPositive",
        )
    except Exception as error:
        return {
            "status": "UNKNOWN_CANDIDATE_EXACT_AUDIT_UNAVAILABLE",
            "vertices": candidate_vertices,
            "faces": candidate_faces,
            "topology_audit": topology,
            "exact_boundary_audit_diagnostic": f"{type(error).__name__}: {error}",
        }
    exact_audit = dict(exact_audit)
    if not (
        exact_audit.get("boundary_valid") is True
        and exact_audit.get("closed") is True
        and exact_audit.get("edge_manifold") is True
        and exact_audit.get("vertex_manifold") is True
        and exact_audit.get("consistently_oriented") is True
        and exact_audit.get("component_orientations_valid") is True
        and exact_audit.get("self_intersection_completed") is True
        and exact_audit.get("self_intersection_free") is True
    ):
        return {
            "status": "UNKNOWN_CANDIDATE_EXACT_AUDIT",
            "vertices": candidate_vertices,
            "faces": candidate_faces,
            "topology_audit": topology,
            "exact_boundary_audit": exact_audit,
            "source_self_intersection_status": source_self_intersection_status,
            "positive_signed_volume": signed_volume,
            "boolean_equivalence_proven": False,
            "candidate_self_intersection_status": (
                "INTERSECTIONS_FOUND"
                if exact_audit.get("self_intersection_completed") is True
                and exact_audit.get("self_intersection_free") is False
                else "UNKNOWN"
            ),
        }
    try:
        evaluator_moments = np.asarray(spin.integrate_raw_mass_moments(
            candidate_vertices.tolist(), candidate_faces.tolist(), density=float(density),
        ), dtype=np.float64)
        independent_moments = _independent_polyhedral_raw_moments(
            candidate_vertices, candidate_faces, float(density),
        )
        moment_audit = _audit_independent_moment_agreement(
            evaluator_moments, independent_moments,
            density=float(density), root_side=float(root_side),
        )
    except Exception as error:
        return {
            "status": "UNKNOWN_CANDIDATE_MOMENT_AUDIT",
            "vertices": candidate_vertices,
            "faces": candidate_faces,
            "topology_audit": topology,
            "exact_boundary_audit": exact_audit,
            "source_self_intersection_status": source_self_intersection_status,
            "positive_signed_volume": signed_volume,
            "boolean_equivalence_proven": False,
            "candidate_moment_audit_diagnostic": f"{type(error).__name__}: {error}",
        }
    return {
        "status": "CAP_CANDIDATE_GEOMETRY_AND_MOMENTS_PASSED",
        "vertices": candidate_vertices,
        "faces": candidate_faces,
        "topology_audit": topology,
        "exact_boundary_audit": exact_audit,
        "evaluator_raw_moments": evaluator_moments.tolist(),
        "independent_raw_moments": independent_moments.tolist(),
        "independent_moment_audit": moment_audit,
        "fan_audit": fan,
        "exact_boundary_coordinate_replacements": replacement_count,
        "maximum_boundary_coordinate_displacement": maximum_displacement,
        "source_self_intersection_status": source_self_intersection_status,
        "positive_signed_volume": signed_volume,
        "boolean_equivalence_proven": False,
        "candidate_exact_audit_face_count": int(exact_audit.get("face_count", len(candidate_faces))),
        "candidate_self_intersection_status": "INTERSECTION_FREE",
    }


def _classify_boolean_face_patch(
    operation: str, source_inside_loop_sample: str, source_inside_face_complement_sample: str,
) -> str:
    """Select a planar Boolean face region only from two strict side witnesses."""
    if operation not in {"intersection", "subtraction"}:
        return "UNKNOWN_BOOLEAN_OPERATION"
    inside_statuses = {"INSIDE_ODD_PARITY", "OUTSIDE_EVEN_PARITY"}
    if source_inside_loop_sample not in inside_statuses or \
            source_inside_face_complement_sample not in inside_statuses:
        return "UNKNOWN_FACE_REGION"
    retain_source = operation == "intersection"
    loop_kept = (source_inside_loop_sample == "INSIDE_ODD_PARITY") == retain_source
    complement_kept = (source_inside_face_complement_sample == "INSIDE_ODD_PARITY") == retain_source
    if loop_kept and not complement_kept:
        return "INTERIOR_DISK"
    if not loop_kept and complement_kept:
        return "FACE_COMPLEMENT"
    return "UNKNOWN_FACE_REGION"


def _audit_planar_boolean_patch_region(
    section: list[tuple[Fraction, Fraction, Fraction]], *, axis: str,
    design_vertices: object, design_faces: object,
    source_vertices: object, source_faces: object, operation: str,
) -> dict[str, object]:
    """Return witness-only evidence for a planar Boolean face patch.

    This identifies the candidate region bounded by one exact simple output
    loop. It does not prove global Boolean equivalence, source validity, or
    that an unobserved source intersection cannot split the face elsewhere.
    No cap or other geometry is generated.
    """
    if axis not in {"x", "y", "z"}:
        raise ValueError("planar patch axis must be x, y, or z")
    if operation not in {"intersection", "subtraction"}:
        return {"status": "UNKNOWN_BOOLEAN_OPERATION", "cap_generated": False}
    if len(section) < 3 or any(len(point) != 3 for point in section):
        return {"status": "UNKNOWN_INVALID_SECTION_LOOP", "cap_generated": False}
    axis_index = {"x": 0, "y": 1, "z": 2}[axis]
    projection = {"x": (1, 2), "y": (0, 2), "z": (0, 1)}[axis]
    plane = section[0][axis_index]
    if any(point[axis_index] != plane for point in section):
        return {"status": "UNKNOWN_NONPLANAR_SECTION", "cap_generated": False}
    section_2d = [(point[projection[0]], point[projection[1]]) for point in section]
    simple = _audit_exact_planar_cycle(section_2d)
    if simple["status"] != "CERTIFIED_SIMPLE_POLYGON":
        return {
            "status": "UNKNOWN_SECTION_NOT_SIMPLE",
            "section_cycle_audit": simple,
            "cap_generated": False,
        }

    domain_points = np.asarray(design_vertices, dtype=np.float64)
    domain_triangles = np.asarray(design_faces)
    if domain_points.ndim != 2 or domain_points.shape[1:] != (3,) or \
            domain_triangles.ndim != 2 or domain_triangles.shape[1:] != (3,) or \
            domain_triangles.dtype.kind not in "iu" or not np.isfinite(domain_points).all():
        return {"status": "UNKNOWN_INVALID_DESIGN_MESH", "cap_generated": False}
    exact_domain = [tuple(Fraction.from_float(float(value)) for value in point) for point in domain_points]
    coplanar_faces = [
        tuple(int(index) for index in face) for face in domain_triangles
        if all(exact_domain[int(index)][axis_index] == plane for index in face)
    ]
    if not coplanar_faces:
        return {"status": "UNKNOWN_NO_COPLANAR_DESIGN_FACE", "cap_generated": False}
    edge_counts: dict[tuple[int, int], int] = {}
    for face in coplanar_faces:
        for first, second in ((face[0], face[1]), (face[1], face[2]), (face[2], face[0])):
            edge = (min(first, second), max(first, second))
            edge_counts[edge] = edge_counts.get(edge, 0) + 1
    if any(count not in {1, 2} for count in edge_counts.values()):
        return {"status": "UNKNOWN_DESIGN_FACE_NONMANIFOLD", "cap_generated": False}
    boundary_edges = [edge for edge, count in edge_counts.items() if count == 1]
    adjacency: dict[int, list[int]] = {}
    for first, second in boundary_edges:
        adjacency.setdefault(first, []).append(second)
        adjacency.setdefault(second, []).append(first)
    if len(adjacency) < 3 or any(len(neighbors) != 2 for neighbors in adjacency.values()):
        return {"status": "UNKNOWN_DESIGN_FACE_NOT_ONE_CYCLE", "cap_generated": False}
    start = min(adjacency)
    cycle_indices = [start]
    previous = None
    current = start
    while True:
        candidates = [candidate for candidate in adjacency[current] if candidate != previous]
        if not candidates:
            return {"status": "UNKNOWN_DESIGN_FACE_NOT_ONE_CYCLE", "cap_generated": False}
        following = min(candidates) if previous is None else candidates[0]
        if following == start:
            break
        if following in cycle_indices:
            return {"status": "UNKNOWN_DESIGN_FACE_NOT_ONE_CYCLE", "cap_generated": False}
        cycle_indices.append(following)
        previous, current = current, following
    if len(cycle_indices) != len(adjacency):
        return {"status": "UNKNOWN_DESIGN_FACE_MULTIPLE_CYCLES", "cap_generated": False}
    face_polygon = [
        (exact_domain[index][projection[0]], exact_domain[index][projection[1]])
        for index in cycle_indices
    ]
    face_simple = _audit_exact_planar_cycle(face_polygon)
    if face_simple["status"] != "CERTIFIED_SIMPLE_POLYGON":
        return {
            "status": "UNKNOWN_DESIGN_FACE_NOT_SIMPLE",
            "design_face_cycle_audit": face_simple,
            "cap_generated": False,
        }
    turns = []
    for index, first in enumerate(face_polygon):
        second = face_polygon[(index + 1) % len(face_polygon)]
        third = face_polygon[(index + 2) % len(face_polygon)]
        turns.append((second[0] - first[0]) * (third[1] - second[1]) -
                     (second[1] - first[1]) * (third[0] - second[0]))
    nonzero_turns = [turn for turn in turns if turn]
    if not nonzero_turns or any((turn > 0) != (nonzero_turns[0] > 0) for turn in nonzero_turns):
        return {"status": "UNKNOWN_DESIGN_FACE_NOT_CONVEX", "cap_generated": False}
    loop_membership = [_exact_point_in_polygon(face_polygon, point) for point in section_2d]
    if any(status != "INSIDE" for status in loop_membership):
        return {
            "status": "UNKNOWN_SECTION_OUTSIDE_DESIGN_FACE",
            "section_vertices_in_design_face": loop_membership,
            "cap_generated": False,
        }

    twice_area = sum(
        first[0] * second[1] - second[0] * first[1]
        for first, second in zip(section_2d, section_2d[1:] + section_2d[:1])
    )
    if twice_area == 0:
        return {"status": "UNKNOWN_ZERO_AREA_SECTION", "cap_generated": False}
    centroid_x = sum(
        (first[0] + second[0]) * (first[0] * second[1] - second[0] * first[1])
        for first, second in zip(section_2d, section_2d[1:] + section_2d[:1])
    ) / (3 * twice_area)
    centroid_y = sum(
        (first[1] + second[1]) * (first[0] * second[1] - second[0] * first[1])
        for first, second in zip(section_2d, section_2d[1:] + section_2d[:1])
    ) / (3 * twice_area)
    if _exact_point_in_polygon(section_2d, (centroid_x, centroid_y)) != "INSIDE":
        return {"status": "UNKNOWN_SECTION_CENTROID_NOT_INTERIOR", "cap_generated": False}

    source_plane_section = _exact_mesh_plane_section_cycles(
        source_vertices, source_faces, axis=axis, plane=plane,
    )
    if source_plane_section["status"] != "CLOSED_SECTION_CYCLES":
        return {
            "status": "UNKNOWN_SOURCE_PLANE_SECTION",
            "source_plane_section_audit": source_plane_section,
            "cap_generated": False,
        }
    if source_plane_section["cycle_count"] != 1:
        return {
            "status": "UNKNOWN_NONUNIQUE_SOURCE_PLANE_SECTION",
            "source_plane_section_audit": source_plane_section,
            "cap_generated": False,
        }
    section_match = _exact_cycle_matches_up_to_collinear_subdivision(
        section, source_plane_section["cycles"][0], axis=axis,
    )
    if not section_match["matches"]:
        return {
            "status": "UNKNOWN_OUTPUT_BOUNDARY_SOURCE_SECTION_MISMATCH",
            "source_plane_section_audit": source_plane_section,
            "output_boundary_section_match": section_match,
            "cap_generated": False,
        }

    query_inside = [Fraction(0), Fraction(0), Fraction(0)]
    query_inside[axis_index] = plane
    query_inside[projection[0]] = centroid_x
    query_inside[projection[1]] = centroid_y
    inside_source = _classify_exact_mesh_point_parity(
        source_vertices, source_faces, tuple(query_inside),
    )
    complement_candidates = [point for point in face_polygon
                             if _exact_point_in_polygon(section_2d, point) == "OUTSIDE"]
    if not complement_candidates:
        return {
            "status": "UNKNOWN_NO_FACE_COMPLEMENT_WITNESS",
            "loop_interior_source_classification": inside_source,
            "cap_generated": False,
        }
    witness_2d = complement_candidates[0]
    query_complement = [Fraction(0), Fraction(0), Fraction(0)]
    query_complement[axis_index] = plane
    query_complement[projection[0]] = witness_2d[0]
    query_complement[projection[1]] = witness_2d[1]
    outside_source = _classify_exact_mesh_point_parity(
        source_vertices, source_faces, tuple(query_complement),
    )
    selected = _classify_boolean_face_patch(
        operation, str(inside_source["status"]), str(outside_source["status"]),
    )
    supported = selected != "UNKNOWN_FACE_REGION"
    return {
        "status": "WITNESS_SUPPORTED_" + selected if supported else "UNKNOWN_FACE_REGION",
        "selected_patch_candidate": selected,
        "operation": operation,
        "axis": axis,
        "plane_coordinate_exact": f"{plane.numerator}/{plane.denominator}",
        "section_cycle_audit": simple,
        "design_face_cycle_audit": face_simple,
        "source_plane_section_audit": {
            key: value for key, value in source_plane_section.items() if key != "cycles"
        },
        "output_boundary_matches_unique_source_section": bool(section_match["matches"]),
        "output_boundary_section_match": section_match,
        "loop_interior_source_classification": inside_source,
        "face_complement_source_classification": outside_source,
        "face_complement_witness_exact": [f"{value.numerator}/{value.denominator}" for value in witness_2d],
        "global_boolean_equivalence_proven": False,
        "source_self_intersection_audited": False,
        "cap_generated": False,
    }


def _classify_exact_mesh_point_parity(
    vertices: object, faces: object, query: tuple[Fraction, Fraction, Fraction],
) -> dict[str, object]:
    """Classify one point by exact ray/triangle parity on a closed mesh.

    A conservative padded floating AABB filter limits exact predicates to
    candidate triangles. Exact rational barycentric and ray-distance tests
    decide the result. This is odd/even parity, not a self-intersection proof
    or a signed-shell nesting interpretation.
    """
    points = np.asarray(vertices, dtype=np.float64)
    triangles = np.asarray(faces)
    if points.ndim != 2 or points.shape[1:] != (3,) or not np.isfinite(points).all():
        raise ValueError("exact point parity requires finite Nx3 mesh vertices")
    if triangles.ndim != 2 or triangles.shape[1:] != (3,) or \
            triangles.dtype.kind not in "iu" or np.any(triangles < 0) or \
            np.any(triangles >= len(points)):
        raise ValueError("exact point parity requires valid integer triangles")
    if len(query) != 3 or not all(isinstance(value, Fraction) for value in query):
        raise ValueError("exact point parity query must contain three Fraction coordinates")
    mesh_audit = _audit_raw_boolean_topology(points, triangles)
    if mesh_audit["status"] != "CLOSED_EDGE_MANIFOLD" or \
            mesh_audit["face_component_count"] != 1 or \
            mesh_audit["orientable_face_component_count"] != 1:
        return {
            "status": "UNKNOWN_SOURCE_NOT_CLOSED_EDGE_MANIFOLD",
            "mesh_topology_status": mesh_audit["status"],
            "mesh_component_count": mesh_audit["face_component_count"],
            "self_intersection_audited": False,
        }

    query_float = tuple(float(value) for value in query)
    projection_axes = (0, 2)  # cast along +y
    scale = max(1.0, float(np.max(np.abs(points))), *(abs(value) for value in query_float))
    padding_exact = max(
        abs(query[axis] - Fraction.from_float(query_float[axis])) for axis in projection_axes
    ) + 8 * Fraction.from_float(math.ulp(scale))
    padding = math.nextafter(float(padding_exact), math.inf)
    query_x, query_z = (query_float[axis] for axis in projection_axes)
    lower_x = math.nextafter(query_x - padding, -math.inf)
    upper_x = math.nextafter(query_x + padding, math.inf)
    lower_z = math.nextafter(query_z - padding, -math.inf)
    upper_z = math.nextafter(query_z + padding, math.inf)
    projected = points[triangles][:, :, projection_axes]
    candidate_mask = (
        (np.min(projected[:, :, 0], axis=1) <= upper_x) &
        (np.max(projected[:, :, 0], axis=1) >= lower_x) &
        (np.min(projected[:, :, 1], axis=1) <= upper_z) &
        (np.max(projected[:, :, 1], axis=1) >= lower_z)
    )
    candidates = np.flatnonzero(candidate_mask)
    point_2d = (query[0], query[2])

    def orient(a, b, c) -> Fraction:
        return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])

    def on_segment(a, b, p) -> bool:
        return min(a[0], b[0]) <= p[0] <= max(a[0], b[0]) and \
            min(a[1], b[1]) <= p[1] <= max(a[1], b[1])

    hits: list[Fraction] = []
    for face_index in candidates:
        face = triangles[int(face_index)]
        exact_triangle = [tuple(Fraction.from_float(float(value)) for value in points[int(index)])
                          for index in face]
        projected_triangle = [(point[0], point[2]) for point in exact_triangle]
        a, b, c = projected_triangle
        area = orient(a, b, c)
        if area == 0:
            if any(orient(first, second, point_2d) == 0 and on_segment(first, second, point_2d)
                   for first, second in ((a, b), (b, c), (c, a))):
                return {
                    "status": "UNKNOWN_RAY_COPLANAR_WITH_FACE",
                    "candidate_triangle_count": int(len(candidates)),
                    "witness_face": int(face_index),
                    "self_intersection_audited": False,
                }
            continue
        weights = (orient(b, c, point_2d), orient(c, a, point_2d), orient(a, b, point_2d))
        if area < 0:
            weights = tuple(-value for value in weights)
            denominator = -area
        else:
            denominator = area
        if all(value > 0 for value in weights):
            y_value = sum((weights[index] * exact_triangle[index][1] for index in range(3)), Fraction(0)) / denominator
            if y_value == query[1]:
                return {
                    "status": "UNKNOWN_QUERY_ON_MESH",
                    "candidate_triangle_count": int(len(candidates)),
                    "witness_face": int(face_index),
                    "self_intersection_audited": False,
                }
            if y_value > query[1]:
                hits.append(y_value)
        elif all(value >= 0 for value in weights) or all(value <= 0 for value in weights):
            return {
                "status": "UNKNOWN_RAY_EDGE_OR_VERTEX_HIT",
                "candidate_triangle_count": int(len(candidates)),
                "witness_face": int(face_index),
                "self_intersection_audited": False,
            }
    if len(set(hits)) != len(hits):
        return {
            "status": "UNKNOWN_COINCIDENT_RAY_CROSSINGS",
            "candidate_triangle_count": int(len(candidates)),
            "positive_ray_intersection_count": len(hits),
            "self_intersection_audited": False,
        }
    return {
        "status": "INSIDE_ODD_PARITY" if len(hits) % 2 else "OUTSIDE_EVEN_PARITY",
        "positive_ray_intersection_count": len(hits),
        "candidate_triangle_count": int(len(candidates)),
        "predicate": "exact rational ray/triangle barycentric and distance tests",
        "broad_phase": "binary64 AABBs padded by exact query-rounding error plus eight scale ulps",
        "self_intersection_audited": False,
    }


def _failed_evaluation_record(
    index: int, coefficients: object, error: Exception,
    *, artifact_directory: Path | None = None,
    operation: str | None = None,
) -> dict[str, object]:
    q = np.asarray(coefficients, dtype=np.float64)
    artifacts: dict[str, dict[str, str]] = {}
    names = {
        "source_geometry": f"source-{index:03d}.obj",
        "raw_boolean_geometry": f"material-{index:03d}.obj",
        "oriented_geometry": f"material-{index:03d}-oriented.obj",
        "construction_provenance": f"material-{index:03d}.obj.construction.tsv",
        "exact_construction_records": f"material-{index:03d}.obj.exact.tsv",
        "face_construction_records": f"material-{index:03d}.obj.faces.tsv",
    }
    if artifact_directory is not None:
        for key, name in names.items():
            path = artifact_directory / name
            if path.is_file():
                artifacts[key] = {"path": name, "sha256": _sha(path)}
    message = str(error).lower()
    selector_details = getattr(error, "selector_details", None)
    if isinstance(selector_details, dict):
        stage = "material_region_selection"
    elif "quadric source motion unavailable:" in message:
        stage = "source_motion_differential"
    elif (
        "quadric has exactly zero center level" in message
        or "unsupported quadratic signature" in message
        or "empty or unsupported quadratic signature" in message
    ):
        stage = "quadric_stratum_classification"
    else:
        stage = ""
    raw_topology: dict[str, object] | None = None
    exact_boundary: dict[str, object] | None = None
    raw_path = artifact_directory / f"material-{index:03d}.obj" if artifact_directory else None
    if raw_path is not None and raw_path.is_file():
        try:
            raw_vertices, raw_faces = _read_obj(raw_path)
            raw_topology = _audit_raw_boolean_topology(raw_vertices, raw_faces)
            exact_path = Path(str(raw_path) + ".exact.tsv")
            face_path = Path(str(raw_path) + ".faces.tsv")
            if exact_path.is_file() and face_path.is_file():
                exact_boundary = _audit_exact_boolean_boundary(
                    raw_vertices, raw_faces, exact_path, face_path,
                )
        except Exception as audit_error:
            raw_topology = {
                "schema": "raw-boolean-topology-diagnostic-v1",
                "status": "AUDIT_FAILED",
                "diagnostic": f"{type(audit_error).__name__}: {audit_error}",
                "repair_performed": False,
            }
    if not stage:
        if "source_geometry" not in artifacts:
            stage = "source_quadric_triangulation"
        elif "raw_boolean_geometry" not in artifacts:
            stage = "cherchi_boolean"
        elif "oriented_geometry" not in artifacts:
            stage = "orientation_normalization"
        elif "construction_provenance" not in artifacts:
            stage = "construction_provenance"
        else:
            stage = "moment_differential_or_integration"
    record = {
        "evaluation_index": int(index),
        "solver_evaluation_index": int(index) + 1,
        "quadric_coefficients": q.tolist(),
        "outcome": "failed",
        "failure_stage": stage,
        "failure_type": type(error).__name__,
        "failure_message": str(error),
        "geometry_hash": None,
        "topology_hash": None,
        "available_artifacts": artifacts,
        "raw_boolean_topology_diagnostic": raw_topology,
        "exact_boolean_boundary_provenance": exact_boundary,
    }
    if operation is not None:
        if operation not in {"intersection", "subtraction"}:
            raise ValueError(f"unsupported Boolean operation in failure evidence: {operation}")
        record["boolean_operation"] = operation
    if isinstance(selector_details, dict):
        record["closed_source_selector"] = selector_details
    return record


def _json_safe(value):
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, np.ndarray):
        return _json_safe(value.tolist())
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, (float, np.floating)) and not np.isfinite(value):
        return "+inf" if value > 0 else "-inf" if value < 0 else "nan"
    if isinstance(value, np.generic):
        return value.item()
    return value


_TRACEBACK_FRAME_LINE_NUMBER = re.compile(r"(?m)^(\s+[^\n]*\.py\()\d+(\): .*)$")


def _canonical_event_digest(events: object) -> str:
    """Hash solver events while ignoring only Python traceback frame line numbers.

    The original native event hash and event payload remain untouched. This
    additional digest lets independent runs compare event meaning when source
    insertions shift traceback locations without changing any diagnostic.
    """
    if not isinstance(events, list) or any(not isinstance(event, dict) for event in events):
        raise ValueError("canonical solver event digest requires a list of event objects")

    def canonicalize(value):
        if isinstance(value, dict):
            return {str(key): canonicalize(item) for key, item in value.items()}
        if isinstance(value, (list, tuple)):
            return [canonicalize(item) for item in value]
        if isinstance(value, str) and "\n\nAt:\n" in value:
            return _TRACEBACK_FRAME_LINE_NUMBER.sub(r"\1<line>\2", value)
        return _json_safe(value)

    payload = json.dumps(
        canonicalize(events), sort_keys=True, separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _scientific_event_digest_v2(events: object) -> str:
    """Hash complete solver events after removing only measured durations.

    The v1 event hashes and original event payload stay available for timing
    and historical audits. This digest gives repeated runs a stable event
    identity without discarding moment values, Newton systems, or failures.
    """
    if not isinstance(events, list) or any(not isinstance(event, dict) for event in events):
        raise ValueError("scientific solver event digest requires a list of event objects")
    scientific = _json_safe(events)
    for event in scientific:
        payload = event.get("evaluation_payload")
        if isinstance(payload, dict):
            payload.pop("evaluator_wall_seconds", None)
            payload.pop("residual_evaluation_seconds", None)
            payload.pop("derivative_evaluation_seconds", None)
    return _canonical_event_digest(scientific)


_EVALUATION_JOURNAL_SCHEMA = "faster-cherchi-evaluation-journal-v1"


class DurableJournalError(RuntimeError):
    """The run cannot safely continue after an evaluation-ledger write fails."""


class RunInterrupted(RuntimeError):
    def __init__(self, signal_number: int):
        super().__init__(f"run interrupted by {signal.Signals(signal_number).name}")
        self.signal_number = int(signal_number)


_PENDING_INTERRUPT_SIGNAL: int | None = None


def _raise_if_interrupt_pending() -> None:
    if _PENDING_INTERRUPT_SIGNAL is not None:
        raise RunInterrupted(_PENDING_INTERRUPT_SIGNAL)


def _durable_json_line(stream, value: object) -> None:
    try:
        stream.write(json.dumps(
            _json_safe(value), sort_keys=True, separators=(",", ":"),
            allow_nan=False,
        ) + "\n")
        stream.flush()
        os.fsync(stream.fileno())
    except Exception as error:
        raise DurableJournalError(
            f"durable evaluation journal write failed: {type(error).__name__}: {error}"
        ) from error


def _initialize_evaluation_journal(path: Path) -> None:
    try:
        with path.open("x", encoding="utf-8") as stream:
            _durable_json_line(stream, {
                "schema": _EVALUATION_JOURNAL_SCHEMA,
                "classification": "METHOD_RUN_EVALUATION_LEDGER",
            })
        directory_fd = os.open(path.parent, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
    except DurableJournalError:
        raise
    except Exception as error:
        raise DurableJournalError(
            f"evaluation journal initialization failed: {type(error).__name__}: {error}"
        ) from error


def _append_evaluation_journal(path: Path, record: dict[str, object]) -> None:
    index = record.get("evaluation_index")
    if isinstance(index, bool) or not isinstance(index, int) or index < 0:
        raise ValueError("evaluation journal requires a nonnegative integer evaluation_index")
    try:
        with path.open("a", encoding="utf-8") as stream:
            _durable_json_line(stream, {"record_type": "evaluation", "record": record})
    except DurableJournalError:
        raise
    except Exception as error:
        raise DurableJournalError(
            f"evaluation journal append failed: {type(error).__name__}: {error}"
        ) from error


_STAGE_TIMING_SCHEMA = "faster-stage-timings-v1"


def _initialize_stage_timing_journal(path: Path) -> None:
    """Create an explicitly non-scientific, durable stage timing ledger."""
    try:
        with path.open("x", encoding="utf-8") as stream:
            _durable_json_line(stream, {
                "schema": _STAGE_TIMING_SCHEMA,
                "classification": "NON_SCIENTIFIC_DIAGNOSTICS",
                "clock": "time.perf_counter_ns",
            })
        directory_fd = os.open(path.parent, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
    except DurableJournalError:
        raise
    except Exception as error:
        raise DurableJournalError(
            f"stage timing journal initialization failed: {type(error).__name__}: {error}"
        ) from error


def _append_stage_timing_record(path: Path, record: dict[str, object]) -> None:
    if record.get("event") not in {"stage_started", "stage_finished"}:
        raise ValueError("stage timing journal requires a start or finish event")
    if not isinstance(record.get("stage_id"), str) or not record["stage_id"]:
        raise ValueError("stage timing journal requires a stable stage_id")
    try:
        with path.open("a", encoding="utf-8") as stream:
            _durable_json_line(stream, record)
    except DurableJournalError:
        raise
    except Exception as error:
        raise DurableJournalError(
            f"stage timing journal append failed: {type(error).__name__}: {error}"
        ) from error


def _read_stage_timing_journal(path: Path) -> list[dict[str, object]]:
    rows = path.read_text(encoding="utf-8").splitlines()
    if not rows:
        raise ValueError("stage timing journal is empty")
    header = json.loads(rows[0])
    if not isinstance(header, dict) or header.get("schema") != _STAGE_TIMING_SCHEMA:
        raise ValueError("unsupported stage timing journal schema")
    records = []
    for line_number, line in enumerate(rows[1:], start=2):
        try:
            record = json.loads(line)
        except Exception as error:
            raise ValueError(f"invalid stage timing row {line_number}: {error}") from error
        if not isinstance(record, dict) or record.get("event") not in {
            "stage_started", "stage_finished",
        }:
            raise ValueError(f"invalid stage timing row {line_number}: malformed record")
        records.append(record)
    return records


def _open_stage_ids(path: Path) -> list[str]:
    if not path.is_file():
        return []
    opened: set[str] = set()
    for line in path.read_text(encoding="utf-8").splitlines()[1:]:
        try:
            record = json.loads(line)
        except (json.JSONDecodeError, UnicodeDecodeError):
            continue
        if not isinstance(record, dict) or not isinstance(record.get("stage_id"), str):
            continue
        if record.get("event") == "stage_started":
            opened.add(record["stage_id"])
        elif record.get("event") == "stage_finished":
            opened.discard(record["stage_id"])
    return sorted(opened)


@contextmanager
def _diagnostic_stage(path: Path, evaluation_index: int, stage: str):
    """Durably bracket one non-scientific stage and preserve interrupted work."""
    stage_id = f"{evaluation_index}:{stage}"
    _append_stage_timing_record(path, {
        "event": "stage_started", "stage_id": stage_id,
        "evaluation_index": evaluation_index, "stage": stage,
    })
    started_ns = time.perf_counter_ns()
    try:
        yield
    except BaseException as error:
        elapsed_ns = time.perf_counter_ns() - started_ns
        _append_stage_timing_record(path, {
            "event": "stage_finished", "stage_id": stage_id,
            "evaluation_index": evaluation_index, "stage": stage,
            "elapsed_ns": elapsed_ns,
            "outcome": "interrupted" if isinstance(error, (RunInterrupted, KeyboardInterrupt)) else "failed",
            "diagnostic": f"{type(error).__name__}: {error}",
        })
        raise
    else:
        elapsed_ns = time.perf_counter_ns() - started_ns
        _append_stage_timing_record(path, {
            "event": "stage_finished", "stage_id": stage_id,
            "evaluation_index": evaluation_index, "stage": stage,
            "elapsed_ns": elapsed_ns, "outcome": "succeeded",
        })


def _run_diagnostic_stage(path: Path, evaluation_index: int, stage: str, operation):
    """Run and time a callback as one durable non-scientific stage."""
    with _diagnostic_stage(path, evaluation_index, stage):
        return operation()


def _initialize_solver_event_journal(path: Path, *, inputs_sha256: str) -> None:
    with path.open("x", encoding="utf-8") as stream:
        _durable_json_line(stream, {
            "schema": "faster-streamed-solver-events-v2",
            "classification": "DETERMINISTIC_SOLVER_AND_CALLBACK_EVENTS",
            "inputs_sha256": inputs_sha256,
        })
    directory_fd = os.open(path.parent, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
    try:
        os.fsync(directory_fd)
    finally:
        os.close(directory_fd)


def _append_solver_journal_event(path: Path, record_type: str, record: dict[str, object]) -> None:
    if record_type not in {
        "solver_event", "geometry_evaluation_started", "geometry_evaluation_finished",
    }:
        raise ValueError("unsupported streamed solver journal record type")
    with path.open("a", encoding="utf-8") as stream:
        # Match the material bridge's typed, lossless diagnostic encoding.
        # In particular, failed Newton systems can legitimately report +inf;
        # keep that distinct from a textual diagnostic containing "+inf".
        _durable_json_line(stream, {"record_type": record_type, "record": _json_value(record)})


def _open_solver_evaluation_ids(path: Path) -> list[int]:
    if not path.is_file():
        return []
    open_ids: set[int] = set()
    for line in path.read_text(encoding="utf-8").splitlines()[1:]:
        try:
            row = json.loads(line)
        except (json.JSONDecodeError, UnicodeDecodeError):
            continue
        record = row.get("record") if isinstance(row, dict) else None
        if not isinstance(row, dict) or not isinstance(record, dict):
            continue
        index = record.get("evaluation_index")
        if isinstance(index, bool) or not isinstance(index, int) or index < 0:
            continue
        if row.get("record_type") == "geometry_evaluation_started":
            open_ids.add(index)
        elif row.get("record_type") == "geometry_evaluation_finished":
            open_ids.discard(index)
    return sorted(open_ids)


def _read_streamed_solver_events(path: Path) -> list[dict[str, object]]:
    rows = path.read_text(encoding="utf-8").splitlines()
    if not rows:
        raise ValueError("streamed solver event journal is empty")
    header = json.loads(rows[0])
    if not isinstance(header, dict) or header.get("schema") not in {
        "faster-streamed-solver-events-v1", "faster-streamed-solver-events-v2",
    }:
        raise ValueError("unsupported streamed solver event journal schema")
    events = []
    for line_number, line in enumerate(rows[1:], start=2):
        try:
            row = json.loads(line)
        except Exception as error:
            raise ValueError(f"invalid streamed solver journal row {line_number}: {error}") from error
        if not isinstance(row, dict) or not isinstance(row.get("record"), dict):
            raise ValueError(f"invalid streamed solver journal row {line_number}: malformed record")
        if row.get("record_type") == "solver_event":
            events.append(row["record"])
        elif row.get("record_type") not in {
            "geometry_evaluation_started", "geometry_evaluation_finished",
        }:
            raise ValueError(f"invalid streamed solver journal row {line_number}: unknown record type")
    return events


def _audit_streamed_solver_events(
    journal_path: Path, finalized_events: list[dict[str, object]],
) -> dict[str, object]:
    streamed = _read_streamed_solver_events(journal_path)
    streamed_digest = _canonical_event_digest(streamed)
    with journal_path.open(encoding="utf-8") as stream:
        schema = json.loads(stream.readline())["schema"]
    if schema == "faster-streamed-solver-events-v2":
        finalized_events = _json_value(finalized_events)
    finalized_digest = _canonical_event_digest(finalized_events)
    if streamed_digest != finalized_digest:
        raise ValueError("streamed native events differ from finalized native events")
    return {
        "status": "PASSED",
        "event_count": len(streamed),
        "canonical_event_sha256": streamed_digest,
        "finalized_native_event_sha256": finalized_digest,
        "sidecar_sha256": _sha(journal_path),
    }


def _read_evaluation_journal(path: Path) -> list[dict[str, object]]:
    rows = path.read_text(encoding="utf-8").splitlines()
    if not rows:
        raise ValueError("evaluation journal is empty")
    try:
        header = json.loads(rows[0])
    except Exception as error:
        raise ValueError(f"invalid evaluation journal header: {error}") from error
    if not isinstance(header, dict) or header.get("schema") != _EVALUATION_JOURNAL_SCHEMA:
        raise ValueError("unsupported evaluation journal schema")
    records: list[dict[str, object]] = []
    for row_number, line in enumerate(rows[1:], start=2):
        try:
            row = json.loads(line)
        except Exception as error:
            raise ValueError(f"invalid evaluation journal row {row_number}: {error}") from error
        record = row.get("record") if isinstance(row, dict) else None
        if not isinstance(row, dict) or row.get("record_type") != "evaluation" or not isinstance(record, dict):
            raise ValueError(f"invalid evaluation journal row {row_number}: malformed row")
        index = record.get("evaluation_index")
        if isinstance(index, bool) or not isinstance(index, int) or index != len(records):
            raise ValueError(f"invalid evaluation journal row {row_number}: noncontiguous evaluation index")
        records.append(record)
    return records


def _write_interruption_receipt(output: Path, signal_number: int) -> dict[str, object]:
    """Persist an explicitly incomplete receipt without inventing solver state."""
    output.mkdir(parents=True, exist_ok=True)
    inputs_path = output / "inputs.json"
    begin_path = output / "begin.json"
    journal_path = output / "evaluation-events.jsonl"
    stage_timing_path = output / "stage-timings.jsonl"
    solver_event_path = output / "solver-events.jsonl"
    journal_bytes = journal_path.read_bytes() if journal_path.is_file() else b""
    valid_evaluations: list[int] = []
    malformed_lines = 0
    for line in journal_bytes.splitlines():
        try:
            row = json.loads(line)
        except (json.JSONDecodeError, UnicodeDecodeError):
            malformed_lines += 1
            continue
        if not isinstance(row, dict) or row.get("record_type") != "evaluation":
            continue
        record = row.get("record")
        if not isinstance(record, dict):
            malformed_lines += 1
            continue
        index = record.get("evaluation_index")
        if isinstance(index, int) and not isinstance(index, bool) and index >= 0:
            valid_evaluations.append(index)
        else:
            malformed_lines += 1
    try:
        signal_name = signal.Signals(signal_number).name
    except ValueError:
        signal_name = f"SIGNAL_{signal_number}"
    receipt: dict[str, object] = {
        "schema": "faster-cherchi-interrupted-run-v1",
        "status": "INTERRUPTED",
        "termination_signal": signal_name,
        "signal_number": int(signal_number),
        "inputs_sha256": _sha(inputs_path) if inputs_path.is_file() else None,
        "begin_sha256": _sha(begin_path) if begin_path.is_file() else None,
        "evaluation_journal_sha256": (
            hashlib.sha256(journal_bytes).hexdigest() if journal_path.is_file() else None
        ),
        "stage_timing_journal_sha256": (
            _sha(stage_timing_path) if stage_timing_path.is_file() else None
        ),
        "open_stages": _open_stage_ids(stage_timing_path),
        "solver_event_journal_sha256": (
            _sha(solver_event_path) if solver_event_path.is_file() else None
        ),
        "open_solver_evaluation_indices": _open_solver_evaluation_ids(solver_event_path),
        "partial_evaluation_count": len(valid_evaluations),
        "last_evaluation_index": max(valid_evaluations) if valid_evaluations else None,
        "malformed_journal_line_count": malformed_lines,
        "solver_terminal_status": "NOT_AVAILABLE",
        "result_json_present": (output / "result.json").is_file(),
        "freshness": "NOT_GENERATED_INCOMPLETE_RUN",
        "sealed_artifact_status": "NOT_VERIFIED",
    }
    receipt_path = output / "interrupted-run.json"
    if receipt_path.exists():
        existing = json.loads(receipt_path.read_text(encoding="utf-8"))
        if existing != receipt:
            raise DurableJournalError("interruption receipt already exists with different contents")
        return existing
    with receipt_path.open("x", encoding="utf-8") as stream:
        stream.write(json.dumps(receipt, sort_keys=True, separators=(",", ":")) + "\n")
        stream.flush()
        os.fsync(stream.fileno())
    directory_fd = os.open(output, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
    try:
        os.fsync(directory_fd)
    finally:
        os.close(directory_fd)
    return receipt


def _supported_quadric_operation(coefficients: np.ndarray) -> str:
    """Pre-screen full-rank quadric strata and their algebraic-center sign.

    For definite forms, the wrong center-side sign describes an empty zero
    set and is rejected. Every nondegenerate indefinite form has a real
    hyperboloid. This return value is only a coarse preflight; the actual
    Boolean side must be selected from the integrated centroid of the current
    finite, capped source geometry by ``_operation_from_closed_source``.
    Rank-deficient and center-degenerate forms require their own charts.
    """
    q = np.asarray(coefficients, dtype=np.float64)
    if q.shape != (10,) or not np.isfinite(q).all():
        raise ValueError("quadric coefficients must be ten finite values")
    matrix = np.asarray(((q[4], q[9], q[8]),
                         (q[9], q[5], q[7]),
                         (q[8], q[7], q[6])))
    eigenvalues = np.linalg.eigvalsh(matrix)
    threshold = 1e-12 * max(1.0, float(np.max(np.abs(eigenvalues))))
    positive = int(np.count_nonzero(eigenvalues > threshold))
    negative = int(np.count_nonzero(eigenvalues < -threshold))
    if positive + negative != 3:
        raise ValueError("unsupported quadratic signature: rank-deficient or ambiguous")
    center_sign = _exact_quadric_center_level_sign(q)
    if center_sign == 0:
        raise ValueError("quadric has exactly zero center level (singular cone)")
    if positive == 3 and center_sign < 0:
        return "subtraction"
    if negative == 3 and center_sign > 0:
        return "intersection"
    if positive > 0 and negative > 0:
        return "intersection" if center_sign > 0 else "subtraction"
    raise ValueError("quadric has an empty or unsupported quadratic signature/center-side combination")


def _exact_quadric_center_level_sign(coefficients: np.ndarray) -> int:
    """Return the exact sign of the center level for binary64 coefficients.

    The center level is ``c - b.T A^-1 b / 4`` for
    ``x.T A x + b.T x + c``. Its sign determines the full-rank quadric
    stratum, including hyperboloids arbitrarily close to a cone. Use exact
    dyadic arithmetic for that topological decision; a floating tolerance
    would incorrectly collapse small but nonzero levels to a singular cone.
    """
    q = np.asarray(coefficients, dtype=np.float64)
    if q.shape != (10,) or not np.isfinite(q).all():
        raise ValueError("quadric coefficients must be ten finite values")
    a = [[Fraction.from_float(float(value)) for value in row] for row in (
        (q[4], q[9], q[8]),
        (q[9], q[5], q[7]),
        (q[8], q[7], q[6]),
    )]
    b = [Fraction.from_float(float(value)) for value in q[1:4]]
    c = Fraction.from_float(float(q[0]))
    det = (
        a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
    )
    if det == 0:
        raise ValueError("quadric quadratic form is exactly rank deficient")
    adjugate = (
        (a[1][1] * a[2][2] - a[1][2] * a[2][1],
         a[0][2] * a[2][1] - a[0][1] * a[2][2],
         a[0][1] * a[1][2] - a[0][2] * a[1][1]),
        (a[1][2] * a[2][0] - a[1][0] * a[2][2],
         a[0][0] * a[2][2] - a[0][2] * a[2][0],
         a[0][2] * a[1][0] - a[0][0] * a[1][2]),
        (a[1][0] * a[2][1] - a[1][1] * a[2][0],
         a[0][1] * a[2][0] - a[0][0] * a[2][1],
         a[0][0] * a[1][1] - a[0][1] * a[1][0]),
    )
    quadratic = sum(
        (b[row] * adjugate[row][column] * b[column]
         for row in range(3) for column in range(3)),
        Fraction(0),
    )
    numerator = 4 * c * det - quadratic
    if numerator == 0:
        return 0
    return 1 if (numerator > 0) == (det > 0) else -1


def _operation_from_closed_source(
    coefficients: np.ndarray,
    source_vertices: np.ndarray,
    source_faces: np.ndarray,
    face_is_cap: np.ndarray,
) -> str:
    """Mirror the native selector using oriented quadric-interface faces."""
    return str(_closed_source_selector_details(
        coefficients, source_vertices, source_faces, face_is_cap,
    )["operation"])


class _ClosedSourceSelectorError(ValueError):
    def __init__(self, message: str, details: dict[str, object]):
        super().__init__(message)
        self.selector_details = details


def _source_side_selection_details(source: dict[str, object]) -> dict[str, object]:
    """Use the native, winding-checked source-side classification for a trial."""
    if not bool(source.get("defined")):
        raise ValueError("quadric source motion unavailable: " + str(source.get("diagnostic")))
    if not bool(source.get("side_selection_defined")):
        details = {
            "status": "UNKNOWN_SIDE",
            "diagnostic": str(source.get("side_selection_diagnostic", "missing side witness")),
        }
        raise _ClosedSourceSelectorError(
            "native quadric-side witness is unavailable", details,
        )
    positive = bool(source["interior_is_positive"])
    operation = "intersection" if positive else "subtraction"
    return {
        "status": "SELECTED",
        "operation": operation,
        "interior_is_positive": positive,
        "method": str(source.get("side_selection_method", "unknown")),
        "diagnostic": str(source.get("side_selection_diagnostic", "")),
    }


def _closed_source_selector_details(
    coefficients: np.ndarray,
    source_vertices: np.ndarray,
    source_faces: np.ndarray,
    face_is_cap: np.ndarray,
) -> dict[str, object]:
    """Classify the source solid's interior sign from its outward interface.

    A volume centroid can lie outside a non-convex capped source. For each
    non-cap triangle, the outward area normal is compared with grad(q). A
    consistent positive alignment means the solid interior is q<0; negative
    alignment means q>0. Disagreement is rejected rather than guessed.
    """
    q = np.asarray(coefficients, dtype=np.float64)
    vertices = np.asarray(source_vertices, dtype=np.float64)
    faces = np.asarray(source_faces, dtype=np.int64)
    caps = np.asarray(face_is_cap, dtype=np.int8)
    if q.shape != (10,) or not np.isfinite(q).all() or vertices.ndim != 2 or \
            vertices.shape[1:] != (3,) or not np.isfinite(vertices).all() or \
            faces.ndim != 2 or faces.shape[1:] != (3,) or \
            caps.shape != (len(faces),) or np.any((caps != 0) & (caps != 1)) or \
            np.any(faces < 0) or np.any(faces >= len(vertices)):
        raise ValueError("closed quadric source selector requires finite geometry and face-cap provenance")
    quadratic = np.asarray(((q[4], q[9], q[8]),
                            (q[9], q[5], q[7]),
                            (q[8], q[7], q[6])))
    positive_exterior = 0
    negative_exterior = 0
    minimum_alignment = math.inf
    for face_index, face in enumerate(faces):
        if caps[face_index]:
            continue
        a, b, c = vertices[face]
        area_vector = np.cross(b - a, c - a)
        center = (a + b + c) / 3.0
        gradient = q[1:4] + 2.0 * (quadratic @ center)
        denominator = float(np.linalg.norm(area_vector) * np.linalg.norm(gradient))
        if not math.isfinite(denominator) or denominator <= 0.0:
            details = {"status": "DEGENERATE_INTERFACE", "face_index": int(face_index)}
            raise _ClosedSourceSelectorError("quadric-side witness is degenerate", details)
        alignment = float(np.dot(area_vector, gradient) / denominator)
        minimum_alignment = min(minimum_alignment, abs(alignment))
        if not math.isfinite(alignment) or abs(alignment) <= 1e-10:
            details = {
                "status": "NEAR_TRANSITION", "face_index": int(face_index),
                "normalized_alignment": alignment,
            }
            raise _ClosedSourceSelectorError("quadric-side witness is ambiguous", details)
        if alignment > 0.0:
            positive_exterior += 1
        else:
            negative_exterior += 1
    if positive_exterior + negative_exterior == 0:
        raise _ClosedSourceSelectorError(
            "quadric-side witness has no interface faces", {"status": "NO_INTERFACE"},
        )
    if positive_exterior and negative_exterior:
        raise _ClosedSourceSelectorError(
            "outward interface faces disagree on the quadric interior sign",
            {"status": "INCONSISTENT_INTERFACE", "positive_votes": positive_exterior,
             "negative_votes": negative_exterior},
        )
    positive_region = negative_exterior > 0
    operation = "intersection" if positive_region else "subtraction"
    normalized_alignment = minimum_alignment
    threshold = 1e-10
    margin = normalized_alignment - threshold
    details = {
        "status": "SELECTED",
        "operation": operation,
        "interior_is_positive": positive_region,
        "positive_exterior_face_count": positive_exterior,
        "negative_exterior_face_count": negative_exterior,
        "minimum_normalized_alignment": normalized_alignment,
        "threshold": threshold,
        "signed_margin": margin,
        "normalized_signed_margin": margin,
        "clearance": margin,
        "normalized_clearance": margin,
    }
    if not math.isfinite(normalized_alignment) or margin <= 0.0:
        details["status"] = "NEAR_TRANSITION"
        raise _ClosedSourceSelectorError(
            "quadric-side witness has insufficient interface alignment", details,
        )
    return details


def _paper_source_feature_digest(payload: object) -> str:
    """Authenticate saved source indexing, not coordinates or solid legality."""
    keys = {'schema', 'surface_type', 'vertex_count', 'normal_vertices',
            'normal_edges', 'projective_rays', 'faces'}
    if not isinstance(payload, dict) or set(payload) != keys or \
            payload.get('schema') != 'paper-projective-source-features-v1' or \
            payload.get('surface_type') not in {'ELLIPSOID', 'HYPERBOLOID'}:
        raise ValueError('invalid paper source feature payload')
    count = payload['vertex_count']
    if type(count) is not int or count < 3:
        raise ValueError('invalid paper source feature vertex count')
    for name, width in [('normal_vertices', 1), ('normal_edges', 2),
                        ('projective_rays', 4), ('faces', 3)]:
        rows = payload[name]
        if not isinstance(rows, list) or (name == 'faces' and not rows):
            raise ValueError('invalid paper source feature rows')
        if name != 'faces' and len(rows) != count and not (
                payload['surface_type'] == 'ELLIPSOID' and not rows):
            raise ValueError('incomplete paper source feature rows')
        for row in rows:
            values = [row] if width == 1 else row
            if not isinstance(values, list) or len(values) != width or any(
                    type(v) is not int or v < (-1 if name != 'faces' else 0)
                    for v in values):
                raise ValueError('invalid paper source feature indices')
            if name == 'faces' and any(v >= count for v in values):
                raise ValueError('out-of-range paper source feature face')
    return digest(payload)


def _pure_height_initial_quadrics(seed, count):
    """PCG64 normals in dual support, then convert to native coefficients.

    No fixture, saved solution or preflight point is an input to this protocol.
    NumPy/runtime identity is frozen with each run; this is a method protocol.
    """
    if type(seed) is not int or not 0 <= seed < 2**32 or type(count) is not int or count <= 0:
        raise ValueError('pure-height seed/count must be nonnegative uint32 and positive integer')
    rng=np.random.Generator(np.random.PCG64(seed))
    values=np.zeros((count,10)); values[:,[0,1,2,3,7,8]]=rng.standard_normal((count,6))
    values[:,7:]/=2.
    values/=np.linalg.norm(values,axis=1)[:,None]
    return tuple(tuple(float(x) for x in row) for row in values)


def _pure_height_source_motion(extension, coefficients, domain, settings):
    """Current native stratum motion chained to all five active directions.

    Does not supply unsupported ambient directions or change the current q.
    Planes/type transitions with no regular stratum derivative fail explicitly.
    """
    q=np.asarray(coefficients,dtype=float)
    u=np.eye(10)[:,[0,1,2,3,7,8]]
    if q.shape!=(10,) or not np.isfinite(q).all() or np.linalg.norm(q)==0 or \
            np.linalg.norm(q-u@(u.T@q))>1e-12:
        raise ValueError('pure-height source coefficients leave the supported coordinates')
    q=q/np.linalg.norm(q)
    source=dict(extension.evaluate_spin_it_faster_stratum_source_motion(q,domain,settings))
    if not source['defined']:
        raise ValueError('pure-height stratum source derivative unavailable: '+str(source['diagnostic']))
    if source['surface_type']!='PARABOLOID':
        raise ValueError('pure-height source currently requires a regular paraboloid stratum')
    h=u.T@q;h[0]+=1 if h[0]>=0 else -1;h/=np.linalg.norm(h)
    basis=u@(np.eye(6)-2*np.outer(h,h))[:,1:]
    native=np.asarray(source['coefficient_basis'],dtype=float)
    velocity=np.asarray(source['vertex_velocity'],dtype=float)
    if native.shape!=(10,8) or velocity.shape!=(len(source['vertices'])*3,8) or \
            not np.isfinite(native).all() or not np.isfinite(velocity).all() or \
            np.linalg.norm(native@(native.T@basis)-basis)>1e-12:
        raise ValueError('native stratum does not span all pure-height directions')
    source.update(native_tangent_basis=basis,vertex_velocity=velocity@(native.T@basis),
        directional_profile='PURE_HEIGHT_ACTIVE_STRATUM_V1',
        paraboloid_grid_resolution=int(settings.paraboloid_grid_resolution))
    return source


def _validate_directional_branch(branch):
    grid=branch.get('paraboloid_grid_resolution')
    if branch.get('directional_profile')!='PURE_HEIGHT_ACTIVE_STRATUM_V1' or \
            branch.get('surface_type')!='PARABOLOID' or type(grid) is not int or grid<1:
        raise ValueError('invalid directional source branch')
    if 'source_map' in branch:
        if 'finite_point_map' in branch or 'paper_source_features_sha256' in branch:
            raise ValueError('fixed-normal source cannot mix full-rank map identities')
        if branch['source_map']!='BOUNDED_FINITE_NORMAL_PARABOLOID_V1':
            raise ValueError('invalid directional source map')
        feature_hash=branch.get('finite_normal_features_sha256')
        if not isinstance(feature_hash,str) or len(feature_hash)!=64 or any(
                c not in '0123456789abcdef' for c in feature_hash):
            raise ValueError('invalid finite-normal source identity')


def _finite_normal_feature_digest(features):
    """Identity of the admitted bounded patch, not a complete ideal-point map."""
    if not isinstance(features,dict) or features.get('schema')!='finite-normal-paraboloid-features-v1' or \
            features.get('surface_type')!='PARABOLOID':
        raise ValueError('invalid finite-normal feature schema')
    normals=np.asarray(features.get('normals'),dtype=float)
    ids=np.asarray(features.get('normal_vertices'))
    faces=np.asarray(features.get('faces'))
    count=features.get('interface_vertex_count')
    if normals.ndim!=2 or normals.shape[1]!=3 or len(normals)==0 or not np.isfinite(normals).all() or \
            np.max(abs(np.linalg.norm(normals,axis=1)-1))>1e-12 or \
            ids.ndim!=1 or len(ids)==0 or not np.issubdtype(ids.dtype,np.integer) or \
            np.min(ids)<0 or np.max(ids)>=len(normals) or \
            faces.ndim!=2 or faces.shape[1]!=3 or len(faces)==0 or not np.issubdtype(faces.dtype,np.integer) or \
            np.min(faces)<0 or np.max(faces)>=len(ids) or type(count) is not int or count<=0 or count*2!=len(ids):
        raise ValueError('invalid finite-normal feature geometry')
    if not np.array_equal(ids[:count],ids[count:]):
        raise ValueError('finite-normal translated closure changed normal identities')
    return hashlib.sha256(json.dumps(features,sort_keys=True,separators=(',',':'),allow_nan=False).encode('ascii')).hexdigest()


def _source_derivative_branch(source: dict[str, object], *,
                              paper_projective_hyperboloid: bool = False,
                              algebraic_finite_map: bool = False):
    """Capture the actual source map; keep historical false-mode hashes intact."""
    branch = {
        'surface_type': str(source['surface_type']),
        'truncation_radius_hex': float(source['truncation_radius']).hex(),
        'initial_icosphere_subdivisions': int(source['initial_icosphere_subdivisions']),
        'effective_icosphere_subdivisions': int(source['effective_icosphere_subdivisions']),
    }
    if source.get('directional_profile') is not None:
        if source['directional_profile'] != 'PURE_HEIGHT_ACTIVE_STRATUM_V1' or paper_projective_hyperboloid or algebraic_finite_map:
            raise ValueError('unsupported or mixed directional source profile')
        branch['directional_profile'] = source['directional_profile']
        branch['paraboloid_grid_resolution'] = int(source['paraboloid_grid_resolution'])
    actual_fixed=source.get('fixed_normal_paraboloid',False)
    if type(actual_fixed) is not bool:
        raise ValueError('invalid fixed-normal source flag')
    if actual_fixed:
        if branch['surface_type']!='PARABOLOID' or paper_projective_hyperboloid or algebraic_finite_map or \
                source.get('paper_projective_hyperboloid',False) is not False or \
                source.get('algebraic_finite_map',False) is not False or \
                branch['initial_icosphere_subdivisions']!=branch['effective_icosphere_subdivisions']:
            raise ValueError('fixed-normal source requires an unmixed directional profile')
        if 'directional_profile' not in branch:
            # The existing 10x8 differential runner attaches and verifies its
            # own versioned fixed-normal feature payload. Preserve that schema;
            # it does not represent a supported-space Newton callback.
            return branch,None
        features={'schema':'finite-normal-paraboloid-features-v1','surface_type':branch['surface_type'],
            'normal_vertices':np.asarray(source['source_normal_vertices']).tolist(),
            'normals':np.asarray(source['fixed_normals']).tolist(),'faces':np.asarray(source['faces']).tolist(),
            'interface_vertex_count':int(source['interface_vertex_count'])}
        branch.update(source_map='BOUNDED_FINITE_NORMAL_PARABOLOID_V1',
            finite_normal_features_sha256=_finite_normal_feature_digest(features))
        return branch,features
    actual_paper = source.get('paper_projective_hyperboloid', False)
    actual_algebraic = source.get('algebraic_finite_map', False)
    if (type(algebraic_finite_map) is not bool or type(actual_algebraic) is not bool or
            actual_algebraic != algebraic_finite_map):
        raise ValueError('algebraic finite map differs from actual source parameterization')
    if algebraic_finite_map:
        branch['finite_point_map'] = 'ALGEBRAIC_V1'
    expected_paper = paper_projective_hyperboloid and branch['surface_type'] == 'HYPERBOLOID'
    if type(actual_paper) is not bool or actual_paper != expected_paper:
        raise ValueError('paper mode differs from actual source parameterization')
    if not paper_projective_hyperboloid:
        return branch, None
    if branch['initial_icosphere_subdivisions'] != branch['effective_icosphere_subdivisions']:
        raise ValueError('paper source requires fixed normal discretization')
    required = ('source_normal_vertices', 'source_normal_edges', 'source_projective_rays')
    if any(key not in source for key in required):
        raise ValueError('paper source feature provenance unavailable')
    features = {
        'schema': 'paper-projective-source-features-v1',
        'surface_type': branch['surface_type'], 'vertex_count': len(source['vertices']),
        'normal_vertices': np.asarray(source['source_normal_vertices']).tolist(),
        'normal_edges': np.asarray(source['source_normal_edges']).tolist(),
        'projective_rays': np.asarray(source['source_projective_rays']).tolist(),
        'faces': np.asarray(source['faces']).tolist(),
    }
    branch['paper_source_features_sha256'] = _paper_source_feature_digest(features)
    return branch, features


def _material_topology_identity(
    operation: str,
    construction: dict[str, object],
    oriented_faces: np.ndarray,
    *,
    derivative_branch: dict[str, object],
) -> str:
    """Bind Newton's derivative stratum to construction, CSG side and cap branch.

    Identical face connectivity/provenance under intersection and subtraction
    does not describe the same material map. A changed truncation cap or source
    discretization also invalidates the previous local Jacobian even when the
    output connectivity happens to remain identical.
    """
    construction_topology = _construction_topology_signature(
        construction, oriented_faces,
    )
    return _derivative_stratum_identity(
        operation, construction_topology, derivative_branch,
    )


def _derivative_stratum_identity(
    operation: str,
    construction_topology: str,
    derivative_branch: dict[str, object],
) -> str:
    """Recompute the Newton identity from fields retained in the run journal."""
    if operation not in {"intersection", "subtraction"}:
        raise ValueError("unsupported Cherchi Boolean operation for topology identity")
    if not isinstance(construction_topology, str) or len(construction_topology) != 64 or any(
        character not in "0123456789abcdef" for character in construction_topology
    ):
        raise ValueError("invalid Cherchi construction topology identity")
    required = {
        "surface_type", "truncation_radius_hex",
        "initial_icosphere_subdivisions", "effective_icosphere_subdivisions",
    }
    if 'directional_profile' in derivative_branch:
        required.update(('directional_profile', 'paraboloid_grid_resolution'))
        _validate_directional_branch(derivative_branch)
        if 'source_map' in derivative_branch:
            required.update(('source_map','finite_normal_features_sha256'))
    if 'finite_point_map' in derivative_branch:
        required.add('finite_point_map')
        if derivative_branch['finite_point_map'] != 'ALGEBRAIC_V1':
            raise ValueError('invalid finite point map identity')
    paper_features = derivative_branch.get('paper_source_features_sha256')
    if 'paper_source_features_sha256' in derivative_branch:
        required.add('paper_source_features_sha256')
        if not isinstance(paper_features, str) or len(paper_features) != 64 or any(
                c not in '0123456789abcdef' for c in paper_features):
            raise ValueError('invalid paper source feature identity')
    if set(derivative_branch) != required:
        raise ValueError("incomplete Cherchi derivative branch identity")
    if not isinstance(derivative_branch["surface_type"], str) or not derivative_branch["surface_type"]:
        raise ValueError("invalid Cherchi source surface type")
    radius = float.fromhex(str(derivative_branch["truncation_radius_hex"]))
    if not math.isfinite(radius) or radius <= 0.0:
        raise ValueError("invalid Cherchi truncation radius")
    initial = derivative_branch["initial_icosphere_subdivisions"]
    effective = derivative_branch["effective_icosphere_subdivisions"]
    if (type(initial) is not int or type(effective) is not int or
            ((initial,effective)!=(-1,-1) if 'directional_profile' in derivative_branch and 'source_map' not in derivative_branch
             else initial < 0 or effective < initial or effective>7)):
        raise ValueError("invalid Cherchi source discretization level")
    if 'source_map' in derivative_branch and initial!=effective:
        raise ValueError('fixed-normal source changed resolution')
    if paper_features is not None and initial != effective:
        raise ValueError('paper source requires fixed normal discretization')
    payload = json.dumps({
        "schema": "faster-cherchi-derivative-stratum-v2",
        "operation": operation,
        "construction_topology": construction_topology,
        "derivative_branch": derivative_branch,
    }, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("ascii")
    return hashlib.sha256(payload).hexdigest()


def _independent_polyhedral_raw_moments(
    vertices: object, faces: object, density: float,
) -> np.ndarray:
    """Integrate a closed Cherchi output independently via signed tetrahedra."""
    points = np.asarray(vertices, dtype=np.float64)
    triangles = np.asarray(faces)
    if points.ndim != 2 or points.shape[1] != 3 or \
            triangles.ndim != 2 or triangles.shape[1] != 3:
        raise ValueError("independent polyhedral integration requires Nx3 vertices and Mx3 faces")
    if triangles.dtype.kind not in "iu" or np.any(triangles < 0) or \
            np.any(triangles >= len(points)):
        raise ValueError("independent polyhedral integration requires valid integer faces")
    if not np.isfinite(points).all() or not np.isfinite(density) or density <= 0.0:
        raise ValueError("independent polyhedral integration requires finite geometry and positive density")
    mesh = TriangleMesh(
        tuple(tuple(float(coordinate) for coordinate in vertex) for vertex in points),
        tuple(tuple(int(index) for index in face) for face in triangles),
    )
    properties = compute_mass_properties(mesh, float(density))
    second = properties.second_moment
    return np.asarray((
        properties.mass,
        *properties.first_moment,
        second[0][0], second[1][1], second[2][2],
        second[1][2], second[0][2], second[0][1],
    ), dtype=np.float64)


def _audit_independent_moment_agreement(
    evaluator_moments: object,
    independent_moments: object,
    *, density: float,
    root_side: float,
    maximum_normalized_residual: float = 1e-6,
) -> dict[str, object]:
    """Compare both integrators using scales frozen from the input geometry."""
    evaluator = np.asarray(evaluator_moments, dtype=np.float64)
    independent = np.asarray(independent_moments, dtype=np.float64)
    if evaluator.shape != (10,) or independent.shape != (10,) or \
            not np.isfinite(evaluator).all() or not np.isfinite(independent).all():
        raise ValueError("independent raw-moment audit requires two finite ten-vectors")
    if not np.isfinite(density) or density <= 0.0 or \
            not np.isfinite(root_side) or root_side <= 0.0 or \
            not np.isfinite(maximum_normalized_residual) or maximum_normalized_residual <= 0.0:
        raise ValueError("independent raw-moment audit requires finite positive input scales")
    powers = np.asarray((3, 4, 4, 4, 5, 5, 5, 5, 5, 5), dtype=np.float64)
    scales = float(density) * np.power(float(root_side), powers)
    if not np.isfinite(scales).all() or np.any(scales <= 0.0):
        raise ValueError("independent raw-moment input scales are not representable")
    residuals = (independent - evaluator) / scales
    normalized_inf = float(np.max(np.abs(residuals)))
    facts = {
        "status": "PASSED" if normalized_inf <= maximum_normalized_residual else "FAILED",
        "normalization_source": "frozen input density and design-domain root side",
        "root_side": float(root_side),
        "density": float(density),
        "maximum_normalized_residual": float(maximum_normalized_residual),
        "normalized_residuals": residuals.tolist(),
        "normalized_inf_residual": normalized_inf,
        "evaluator_raw_moments": evaluator.tolist(),
        "independent_raw_moments": independent.tolist(),
    }
    if normalized_inf > maximum_normalized_residual:
        raise ValueError(
            f"independent raw-moment residual {normalized_inf:.17g} exceeds "
            f"{maximum_normalized_residual:.17g}"
        )
    return facts


def _moment_map_binding_sha256(
    geometry_hash: str,
    topology_hash: str,
    operation: str,
    raw_moments: object,
    moment_jacobian: object,
    *,
    quadric_coefficients: object,
    backend_identity: str,
    derivative_identity: str,
    source_identity: str,
    cap_identity: str,
    exact_geometry_identity: str = 'NO_EXACT_RATIONAL_SIDECAR',
    native_tangent_basis=None,
    support_basis=None,
) -> str:
    """Bind one map to coefficients, discrete geometry, backend and derivative stratum."""
    for label, value in (("geometry", geometry_hash), ("topology", topology_hash)):
        if not isinstance(value, str) or len(value) != 64 or any(
            character not in "0123456789abcdef" for character in value
        ):
            raise ValueError(f"{label} identity must be a lowercase SHA-256 digest")
    if operation not in {"intersection", "subtraction"}:
        raise ValueError("unsupported Boolean operation in moment-map binding")
    moments = np.asarray(raw_moments, dtype=np.float64)
    jacobian = np.asarray(moment_jacobian, dtype=np.float64)
    coefficients = np.asarray(quadric_coefficients, dtype=np.float64)
    directional = native_tangent_basis is not None or support_basis is not None
    if directional:
        b=np.asarray(native_tangent_basis,dtype=float);u=np.asarray(support_basis,dtype=float)
        if u.ndim!=2 or u.shape[0]!=10 or not 2<=u.shape[1]<=10 or \
                b.shape!=(10,u.shape[1]-1) or jacobian.shape!=b.shape or \
                not np.isfinite(u).all() or not np.isfinite(b).all() or \
                np.linalg.norm(u.T@u-np.eye(u.shape[1]))>1e-12 or \
                np.linalg.norm(b.T@b-np.eye(b.shape[1]))>1e-12:
            raise ValueError('invalid genuine directional basis/support dimensions or orthonormality')
    if moments.shape != (10,) or (not directional and jacobian.shape != (10, 10)) or \
            coefficients.shape != (10,) or not np.isfinite(moments).all() or \
            not np.isfinite(jacobian).all() or not np.isfinite(coefficients).all():
        raise ValueError("moment-map binding requires finite coefficient/moment 10-vectors and 10-by-10 Jacobian")
    identities = (backend_identity, derivative_identity, source_identity, cap_identity,
                  exact_geometry_identity)
    if any(not isinstance(value, str) or not value or ";" in value for value in identities):
        raise ValueError("moment-map backend, derivative, source, cap and exact geometry identities must be nonempty delimiter-safe strings")

    digest = hashlib.sha256()
    digest.update(b"faster-newton-map-v4-directional-exact-geometry\0" if directional else
                  b"faster-newton-map-v3-exact-geometry\0")
    for text in (geometry_hash, topology_hash, operation, *identities):
        encoded = text.encode("ascii")
        digest.update(struct.pack(">Q", len(encoded)))
        digest.update(encoded)
    digest.update(struct.pack(">10d", *coefficients.tolist()))
    digest.update(struct.pack(">Q", len(moments)))
    digest.update(struct.pack(">QQ", *jacobian.shape))
    digest.update(struct.pack(">10d", *moments.tolist()))
    digest.update(struct.pack(f">{jacobian.size}d", *jacobian.reshape(-1, order="C").tolist()))
    if directional:
        if np.linalg.norm(b.T@coefficients)>1e-12:
            raise ValueError('directional basis is not tangent at the bound coefficients')
        for matrix in (u,b):
            digest.update(struct.pack('>QQ',*matrix.shape))
            digest.update(struct.pack(f'>{matrix.size}d',*matrix.reshape(-1).tolist()))
    return digest.hexdigest()


def _audit_newton_evaluation_bindings(
    evaluations: list[dict[str, object]],
    events: list[dict[str, object]],
    *,
    quadric_hash: Callable[[object], str],
) -> dict[str, object]:
    """Rejoin evaluator responses to Newton events and recompute their native hashes."""
    successful: dict[int, dict[str, object]] = {}
    for evaluation in evaluations:
        if evaluation.get("outcome") == "failed" or evaluation.get("geometry_hash") is None:
            continue
        callback_index = int(evaluation["evaluation_index"])
        index = int(evaluation["solver_evaluation_index"])
        if index != callback_index + 1:
            raise ValueError("callback and Newton evaluation index convention mismatch")
        if index in successful:
            raise ValueError("duplicate successful evaluation index in moment-map records")
        mapping = evaluation.get("moment_map")
        if not isinstance(mapping, dict):
            raise ValueError("successful evaluation omitted its moment-map tuple")
        identity_fields = (
            "evaluation_map_identity", "boolean_backend_identity", "derivative_identity",
            "source_identity", "cap_identity",
        )
        if any(not isinstance(evaluation.get(field), str) or not evaluation[field]
               for field in identity_fields):
            raise ValueError("successful evaluation omitted complete map/backend/source/cap identities")
        if evaluation["evaluation_map_identity"] != mapping.get("sha256"):
            raise ValueError("evaluation map identity differs from its moment-map digest")
        recomputed = _moment_map_binding_sha256(
            str(evaluation["geometry_hash"]),
            str(evaluation["topology_hash"]),
            str(evaluation["boolean_operation"]),
            mapping.get("raw_moments"), mapping.get("moment_jacobian"),
            quadric_coefficients=evaluation.get("quadric_coefficients"),
            backend_identity=str(evaluation["boolean_backend_identity"]),
            derivative_identity=str(evaluation["derivative_identity"]),
            source_identity=str(evaluation["source_identity"]),
            cap_identity=str(evaluation["cap_identity"]),
            exact_geometry_identity=str(evaluation.get(
                "exact_geometry_identity", "NO_EXACT_RATIONAL_SIDECAR")),
            native_tangent_basis=mapping.get('native_tangent_basis'),
            support_basis=mapping.get('support_basis'),
        )
        if mapping.get("sha256") != recomputed:
            raise ValueError("moment-map binding hash mismatch")
        source_geometry_hash = evaluation.get("source_geometry_hash")
        source_triangulation_hash = evaluation.get("source_triangulation_hash")
        if not isinstance(source_geometry_hash, str) or len(source_geometry_hash) != 64 or \
                any(character not in "0123456789abcdef" for character in source_geometry_hash) or \
                not isinstance(source_triangulation_hash, str) or not source_triangulation_hash:
            raise ValueError("successful evaluation omitted source geometry/triangulation identity")
        expected_source_identity = hashlib.sha256(
            (source_geometry_hash + "\0" + source_triangulation_hash).encode("ascii")
        ).hexdigest()
        if evaluation["source_identity"] != expected_source_identity:
            raise ValueError("source identity does not match source geometry and triangulation")
        cap_payload = evaluation.get("cap_identity_payload")
        if not isinstance(cap_payload, dict) or \
                cap_payload.get("triangulation_hash") != source_triangulation_hash:
            raise ValueError("successful evaluation omitted a cap payload bound to its triangulation")
        required_cap_fields = {
            "surface_type", "truncation_radius_hex", "initial_icosphere_subdivisions",
            "effective_icosphere_subdivisions", "maximum_radius_doublings", "triangulation_hash",
        }
        if 'directional_profile' in cap_payload:
            required_cap_fields.update(('directional_profile','paraboloid_grid_resolution'))
            _validate_directional_branch(cap_payload)
            if 'source_map' in cap_payload:
                required_cap_fields.update(('source_map','finite_normal_features_sha256'))
        if 'finite_point_map' in cap_payload:
            required_cap_fields.add('finite_point_map')
            if cap_payload['finite_point_map'] != 'ALGEBRAIC_V1':
                raise ValueError('invalid cap finite point map identity')
        if 'paper_source_features_sha256' in cap_payload:
            required_cap_fields.add('paper_source_features_sha256')
            features = evaluation.get('source_feature_payload')
            if _paper_source_feature_digest(features) != cap_payload['paper_source_features_sha256'] or \
                    features['surface_type'] != cap_payload['surface_type']:
                raise ValueError('paper source feature identity mismatch')
        elif 'source_map' in cap_payload:
            if _finite_normal_feature_digest(evaluation.get('source_feature_payload'))!=cap_payload['finite_normal_features_sha256']:
                raise ValueError('finite-normal source feature identity mismatch')
        elif evaluation.get('source_feature_payload') is not None:
            raise ValueError('unbound paper source feature payload')
        if set(cap_payload) != required_cap_fields or \
                cap_payload.get("surface_type") not in (
                    {'PARABOLOID'} if 'directional_profile' in cap_payload else {'ELLIPSOID', 'HYPERBOLOID'}) or \
                not isinstance(cap_payload.get("truncation_radius_hex"), str) or \
                not isinstance(cap_payload.get("triangulation_hash"), str) or \
                len(str(cap_payload["triangulation_hash"])) != 64 or \
                any(character not in "0123456789abcdef" for character in cap_payload["triangulation_hash"]):
            raise ValueError("cap identity payload schema is invalid")
        try:
            cap_radius = float.fromhex(str(cap_payload["truncation_radius_hex"]))
        except ValueError as error:
            raise ValueError("cap identity payload schema is invalid") from error
        subdivisions = (
            cap_payload["initial_icosphere_subdivisions"],
            cap_payload["effective_icosphere_subdivisions"],
        )
        radius_doublings = cap_payload["maximum_radius_doublings"]
        if not math.isfinite(cap_radius) or cap_radius <= 0.0 or any(
            isinstance(value, bool) or not isinstance(value, int) or
            (value!=-1 if 'directional_profile' in cap_payload and 'source_map' not in cap_payload else not 0 <= value <= 7)
            for value in subdivisions
        ) or subdivisions[0] > subdivisions[1] or isinstance(radius_doublings, bool) or \
                not isinstance(radius_doublings, int) or not 0 <= radius_doublings <= 32:
            raise ValueError("cap identity payload schema is invalid")
        expected_cap_identity = hashlib.sha256(json.dumps(
            cap_payload, sort_keys=True, separators=(",", ":"), allow_nan=False,
        ).encode("ascii")).hexdigest()
        if evaluation["cap_identity"] != expected_cap_identity:
            raise ValueError("cap identity does not match cap payload")
        derivative_branch = {
            key: cap_payload[key] for key in (
                "surface_type", "truncation_radius_hex",
                "initial_icosphere_subdivisions", "effective_icosphere_subdivisions",
            )
        }
        if 'paper_source_features_sha256' in cap_payload:
            derivative_branch['paper_source_features_sha256'] = cap_payload['paper_source_features_sha256']
        if 'finite_point_map' in cap_payload:
            derivative_branch['finite_point_map'] = cap_payload['finite_point_map']
        if 'directional_profile' in cap_payload:
            derivative_branch.update({key:cap_payload[key] for key in
                ('directional_profile','paraboloid_grid_resolution')})
            if 'source_map' in cap_payload:
                derivative_branch.update({key:cap_payload[key] for key in ('source_map','finite_normal_features_sha256')})
        if evaluation.get("derivative_branch") != derivative_branch or \
                evaluation.get("topology_hash") != _derivative_stratum_identity(
                    str(evaluation.get("boolean_operation")),
                    evaluation.get("construction_topology_hash"),
                    derivative_branch,
                ):
            raise ValueError("Cherchi derivative stratum identity mismatch")
        successful[index] = evaluation

    if not successful:
        raise ValueError("no successful evaluation maps to audit")

    event_bindings: dict[int, dict[str, object]] = {}
    events_by_index: dict[int, dict[str, object]] = {}
    for event in events:
        if event.get("type") == "evaluation":
            index = int(event["evaluation_index"])
            if index in events_by_index:
                raise ValueError("duplicate Newton evaluation event index")
            events_by_index[index] = event
        geometry_hash = event.get("material_geometry_hash")
        topology_hash = event.get("material_topology_hash")
        if geometry_hash is None and topology_hash is None:
            continue
        if geometry_hash is None or topology_hash is None:
            raise ValueError("Newton event has incomplete material geometry identity")
        index = int(event["evaluation_index"])
        if index in event_bindings:
            raise ValueError("duplicate Newton event material identity for evaluation")
        if index not in successful:
            raise ValueError("Newton event material identity has no successful evaluator record")
        evaluation = successful[index]
        if evaluation["geometry_hash"] != geometry_hash or evaluation["topology_hash"] != topology_hash:
            raise ValueError("Newton event geometry identity mismatch")
        event_quadric_hash = event.get("quadric_hash")
        if not isinstance(event_quadric_hash, str) or len(event_quadric_hash) != 64 or any(
            character not in "0123456789abcdef" for character in event_quadric_hash
        ):
            raise ValueError("Newton event omitted quadric hash")
        evaluation_quadric_hash = evaluation.get("quadric_hash")
        if evaluation_quadric_hash != event_quadric_hash:
            raise ValueError("Newton event quadric hash mismatch")
        event_coefficients = np.asarray(event.get("quadric_coefficients"), dtype=np.float64)
        evaluation_coefficients = np.asarray(
            evaluation.get("quadric_coefficients"), dtype=np.float64,
        )
        if event_coefficients.shape != (10,) or evaluation_coefficients.shape != (10,) or \
                not np.isfinite(event_coefficients).all() or \
                not np.isfinite(evaluation_coefficients).all() or \
                not np.array_equal(event_coefficients, evaluation_coefficients):
            raise ValueError("Newton event quadric coefficients mismatch")
        recomputed_quadric_hash = quadric_hash(event_coefficients)
        if not isinstance(recomputed_quadric_hash, str) or \
                recomputed_quadric_hash != event_quadric_hash:
            raise ValueError("Newton event quadric hash does not match native coefficient digest")
        detail = event.get("detail")
        if not isinstance(detail, str):
            raise ValueError("Newton evaluation event omitted map identity detail")
        def event_identity(label: str) -> str:
            marker = ";" + label + "="
            if detail.count(marker) != 1:
                raise ValueError("Newton evaluation event omitted or duplicated " + label + " identity")
            return detail.split(marker, 1)[1].split(";", 1)[0]
        if event_identity("evaluation_map") != evaluation["evaluation_map_identity"]:
            raise ValueError("Newton event moment-map identity mismatch")
        if event_identity("backend") != evaluation["boolean_backend_identity"]:
            raise ValueError("Newton event backend identity mismatch")
        if event_identity("derivative") != evaluation["derivative_identity"]:
            raise ValueError("Newton event derivative identity mismatch")
        event_bindings[index] = event
    unbound = set(successful) - set(event_bindings)
    rejected_geometry_evaluations = []
    for index in sorted(unbound):
        event = events_by_index.get(index)
        if event is None or event.get("outcome") != "failed" or \
                event.get("reason") != "evaluation_exception" or \
                event.get("material_geometry_hash") is not None or \
                event.get("material_topology_hash") is not None or \
                not isinstance(event.get("detail"), str) or not event["detail"]:
            raise ValueError("successful evaluator records and Newton events do not match")
        rejected_geometry_evaluations.append(index)
    if set(event_bindings) | set(rejected_geometry_evaluations) != set(successful):
        raise ValueError("successful evaluator records and Newton events do not match")
    return {
        "status": "CONSISTENT",
        "bound_evaluation_count": len(event_bindings),
        "callback_map_count": len(successful),
        "solver_rejected_geometry_evaluations": rejected_geometry_evaluations,
        "evaluation_indices": sorted(event_bindings),
        "moment_map_hashes": [successful[index]["moment_map"]["sha256"]
                              for index in sorted(event_bindings)],
        "map_identities": [{
            "evaluation_index": index,
            "evaluation_map_identity": successful[index]["evaluation_map_identity"],
            "boolean_backend_identity": successful[index]["boolean_backend_identity"],
            "derivative_identity": successful[index]["derivative_identity"],
            "source_identity": successful[index]["source_identity"],
            "cap_identity": successful[index]["cap_identity"],
        } for index in sorted(event_bindings)],
    }


def _audit_material_domain_subset(
    material_vertices: np.ndarray,
    material_faces: np.ndarray,
    design_vertices: np.ndarray,
    design_faces: np.ndarray,
    *, material_snapshot=None,
) -> dict[str, object]:
    """Prove exact represented-set inclusion with the pinned exact Nef checker."""
    native = spin._native_module()
    make_region = getattr(native, "make_exact_material_region", None)
    check_subset = getattr(native, "check_exact_material_subset", None)
    if make_region is None or check_subset is None:
        raise ValueError("exact material-domain relation API is unavailable")
    if material_snapshot is None:
        material_region = make_region(material_vertices, material_faces)
    else:
        snapshot_vertices, snapshot_faces = material_snapshot.arrays()
        if not np.array_equal(snapshot_vertices, material_vertices) or \
                not np.array_equal(snapshot_faces, material_faces):
            raise ValueError("exact subset snapshot differs from the bound material mesh")
        material_region = material_snapshot.exact_region()
    design_region = make_region(design_vertices, design_faces)
    material_facts = dict(material_region.facts)
    design_facts = dict(design_region.facts)
    if not material_facts.get("completed") or not material_facts.get("valid"):
        raise ValueError("terminal material is not a valid exact represented set")
    if not design_facts.get("completed") or not design_facts.get("valid"):
        raise ValueError("request design domain is not a valid exact represented set")
    relation = dict(check_subset(material_region, design_region))
    included = relation.get("a_subset_b")
    if not relation.get("completed") or type(included) is not bool:
        raise ValueError(
            "exact material-domain subset query is incomplete: " +
            str(relation.get("error", "unknown exact Boolean failure"))
        )
    return {
        "status": "VERIFIED_SUBSET" if included else "OUTSIDE_DESIGN_DOMAIN",
        "completed": True,
        "material_subset_of_design_domain": included,
        "method": "EXACT_NEF_INCLUSIVE_SUBSET",
        "scope": ("exact-rational-material-vs-binary64-design-polyhedral-sets-only"
                  if material_snapshot is not None else
                  "exact-binary64-represented-polyhedral-sets-only"),
        "aabb_diagnostic_only": {
            "material_min": np.min(material_vertices, axis=0).tolist(),
            "material_max": np.max(material_vertices, axis=0).tolist(),
            "design_min": np.min(design_vertices, axis=0).tolist(),
            "design_max": np.max(design_vertices, axis=0).tolist(),
            "vertices_strictly_outside_design_aabb": int(np.count_nonzero(np.any(
                (material_vertices < np.min(design_vertices, axis=0)) |
                (material_vertices > np.max(design_vertices, axis=0)), axis=1,
            ))),
            "not_a_containment_proof": True,
        },
        "material_region_digest": material_facts.get("region_digest"),
        "design_domain_region_digest": design_facts.get("region_digest"),
        "relation": relation,
    }


def _audit_terminal_geometry_binding(
    result_record: dict[str, object],
    *,
    quadric_hash: Callable[[object], str],
) -> dict[str, object]:
    """Independently audit the exact accepted mesh at a converged Newton endpoint.

    This deliberately does not rewrite the native solver's
    ``final_geometry_verified`` claim. It proves that the accepted trace
    endpoint is bound to the same Cherchi material. Authenticated rational
    coordinates, when present, remain authoritative for boundary validity,
    inclusion and independent integration. The floating derived mesh is
    separately checked and may stay rejected. Neither route grants printing
    or physical validation.
    """
    if result_record.get("solver_status") != "CONVERGED":
        return {
            "status": "NOT_APPLICABLE",
            "geometry_verified": False,
            "reason": "solver did not converge",
            "material_domain_status": "NOT_VERIFIED",
        }
    if result_record.get("evaluation_map_binding_audit", {}).get("status") != "CONSISTENT":
        raise ValueError("terminal geometry audit requires a consistent Newton evaluation-map audit")
    trace = result_record.get("trace")
    evaluations = result_record.get("evaluations")
    if not isinstance(trace, list) or not trace or not isinstance(evaluations, list):
        raise ValueError("converged terminal geometry audit requires trace and evaluation records")
    terminal = trace[-1]
    if not isinstance(terminal, dict):
        raise ValueError("terminal Newton trace row is malformed")
    # ``function_evaluations`` is counted per restart in the native trace,
    # while external evaluation records use a run-global index. Bind the
    # terminal material by its scientific identities instead of comparing
    # those different counters.
    terminal_index = terminal.get("function_evaluations")
    if isinstance(terminal_index, bool) or not isinstance(terminal_index, int) or terminal_index <= 0:
        raise ValueError("terminal Newton trace omitted its per-restart evaluation count")
    matches = [
        row for row in evaluations
        if isinstance(row, dict)
        and row.get("quadric_hash") == terminal.get("quadric_hash")
        and row.get("geometry_hash") == terminal.get("material_geometry_hash")
        and row.get("topology_hash") == terminal.get("material_topology_hash")
    ]
    if len(matches) != 1:
        raise ValueError("terminal Newton trace has no unique material identity match")
    evaluation = matches[0]
    global_evaluation_index = evaluation.get("solver_evaluation_index")
    if isinstance(global_evaluation_index, bool) or not isinstance(global_evaluation_index, int) or \
            global_evaluation_index <= 0:
        raise ValueError("terminal material evaluation omitted its global solver evaluation index")
    coefficients = np.asarray(evaluation.get("quadric_coefficients"), dtype=np.float64)
    if coefficients.shape != (10,) or not np.isfinite(coefficients).all() or \
            quadric_hash(coefficients) != terminal.get("quadric_hash") or \
            evaluation.get("quadric_hash") != terminal.get("quadric_hash"):
        raise ValueError("terminal trace and material quadric identity differ")
    if evaluation.get("geometry_hash") != terminal.get("material_geometry_hash") or \
            evaluation.get("topology_hash") != terminal.get("material_topology_hash"):
        raise ValueError("terminal trace and material geometry identity differ")
    mesh_path = Path(str(evaluation.get("oriented_output_path", ""))).resolve(strict=True)
    geometry_digest = _sha(mesh_path)
    if geometry_digest != evaluation.get("geometry_hash"):
        raise ValueError("terminal oriented material bytes differ from the bound geometry hash")
    vertices, faces = _read_obj(mesh_path)
    if vertices.ndim != 2 or vertices.shape[1] != 3 or faces.ndim != 2 or faces.shape[1] != 3:
        raise ValueError("terminal oriented material must be a triangular 3D boundary")
    floating_audit = dict(spin._native_module().audit_exact_triangle_boundary(
        vertices, faces, "OutwardPositive",
    ))
    required = (
        "boundary_valid", "closed", "edge_manifold", "vertex_manifold",
        "consistently_oriented", "component_orientations_valid",
        "self_intersection_completed", "self_intersection_free",
    )
    float_valid = all(floating_audit.get(field) is True for field in required)
    exact_snapshot = None
    authority = "BINARY64_REPRESENTED_BOUNDARY"
    exact_audit = floating_audit
    raw_path = evaluation.get("output_path")
    sidecar = None if raw_path is None else Path(str(raw_path) + ".exact.tsv")
    has_exact_record = evaluation.get("exact_geometry_identity") is not None or \
        evaluation.get("exact_float_derivation_audit") is not None
    if has_exact_record or (sidecar is not None and sidecar.exists()):
        bound_derivation = evaluation.get("exact_float_derivation_audit")
        if sidecar is None or not sidecar.is_file() or not isinstance(bound_derivation, dict) or \
                _sha(sidecar) != bound_derivation.get("sidecar_sha256"):
            raise ValueError("terminal exact sidecar identity differs from the bound evaluation")
        exact_snapshot, derivation = mesh_snapshot_from_exact_boolean_sidecar(
            vertices, faces, sidecar, return_audit=True)
        if derivation != bound_derivation or digest({
                "rational_geometry": exact_snapshot.identity,
                "float_derivation_audit": derivation,
        }) != evaluation.get("exact_geometry_identity"):
            raise ValueError("terminal exact material identity differs from the bound evaluation")
        authority = "EXACT_RATIONAL_BOUNDARY"
        exact_audit = {
            "method": "NATIVE_RATIONAL_BOUNDARY_AND_SIGNED_NEF_VALIDATION",
            "region_facts": dict(exact_snapshot.exact_region().facts),
            "scope": "AUTHORITATIVE_EXACT_MATERIAL_NOT_FLOAT_OR_MANUFACTURING_APPROVAL",
        }
    elif not float_valid:
        raise ValueError("terminal material failed exact represented-boundary audit")
    design_vertices = np.asarray(result_record.get("design_domain_vertices"), dtype=np.float64)
    design_faces = np.asarray(result_record.get("design_domain_faces"), dtype=np.int64)
    if design_vertices.ndim != 2 or design_vertices.shape[1] != 3 or \
            design_faces.ndim != 2 or design_faces.shape[1] != 3:
        raise ValueError("terminal geometry audit requires the request's triangular design domain")
    domain_relation = _audit_material_domain_subset(
        # The sidecar already authenticated the original OBJ's bounded
        # conversion. Set inclusion concerns the authoritative rational
        # snapshot, not equality of its canonical display with that OBJ.
        exact_snapshot.arrays()[0] if exact_snapshot is not None else vertices,
        faces, design_vertices, design_faces,
        material_snapshot=exact_snapshot,
    )
    if domain_relation["material_subset_of_design_domain"] is not True:
        return {
            "status": "FAILED",
            "geometry_verified": False,
            "scope": "represented-boundary-and-exact-design-domain-subset-audit",
            "terminal_solver_evaluation_index": global_evaluation_index,
            "terminal_restart_function_evaluations": terminal_index,
            "terminal_quadric_hash": terminal["quadric_hash"],
            "terminal_geometry_hash": geometry_digest,
            "terminal_material_path": str(mesh_path),
            "terminal_material_sha256": geometry_digest,
            "exact_boundary_audit": exact_audit,
            "material_authority": authority,
            "floating_boundary_audit": floating_audit,
            "float_validation_status": "PASSED" if float_valid else "REJECTED",
            "material_domain_status": domain_relation["status"],
            "material_domain_relation": domain_relation,
            "diagnostic": "terminal material is not contained in the represented design domain",
            "native_solver_claim_preserved": result_record.get("final_geometry_verified"),
        }
    moment_map = evaluation.get("moment_map")
    independent_audit = evaluation.get("independent_moment_audit")
    if not isinstance(moment_map, dict) or not isinstance(independent_audit, dict) or \
            independent_audit.get("status") != "PASSED":
        raise ValueError("terminal material lacks a passed independent moment record")
    density = float(independent_audit.get("density", math.nan))
    root_side = float(independent_audit.get("root_side", math.nan))
    bound_moments = np.asarray(moment_map.get("raw_moments"), dtype=np.float64)
    if bound_moments.shape != (10,) or not np.isfinite(bound_moments).all():
        raise ValueError("terminal material evaluation has malformed moment map")
    if exact_snapshot is None:
        evaluator_moments = np.asarray(spin.integrate_raw_mass_moments(
            vertices.tolist(), faces.tolist(), density,
        ), dtype=np.float64)
        independent_moments = _independent_polyhedral_raw_moments(vertices, faces, density)
    else:
        # Recompute all ten terms from the authenticated rational boundary,
        # independently of the bound Newton map. The invalid floating display
        # must neither supply authoritative moments nor receive approval.
        evaluator_moments = bound_moments
        independent_moments = np.asarray(_raw_moments(exact_snapshot, density), dtype=np.float64)
    recomputed_moment_audit = _audit_independent_moment_agreement(
        evaluator_moments, independent_moments,
        density=density, root_side=root_side,
    )
    scales = density * np.power(root_side, np.asarray((3, 4, 4, 4, 5, 5, 5, 5, 5, 5)))
    map_residual = float(np.max(np.abs(independent_moments - bound_moments) / scales))
    if not math.isfinite(map_residual) or map_residual > 1e-6:
        raise ValueError("terminal boundary moments differ from the converged Newton map")
    return {
        "status": "PASSED",
        "geometry_verified": True,
        "scope": ("authoritative-exact-material-not-float-export-or-manufacturing"
                  if exact_snapshot is not None else
                  "represented-boundary-manifold-self-intersection-and-moment-audit"),
        "terminal_solver_evaluation_index": global_evaluation_index,
        "terminal_restart_function_evaluations": terminal_index,
        "terminal_quadric_hash": terminal["quadric_hash"],
        "terminal_geometry_hash": geometry_digest,
        "terminal_topology_hash": terminal["material_topology_hash"],
        "terminal_material_path": str(mesh_path),
        "terminal_material_sha256": geometry_digest,
        "vertices": int(len(vertices)),
        "faces": int(len(faces)),
        "exact_boundary_audit": exact_audit,
        "material_authority": authority,
        "floating_boundary_audit": floating_audit,
        "float_validation_status": "PASSED" if float_valid else "REJECTED",
        "independent_moment_audit": recomputed_moment_audit,
        "newton_map_moment_residual": map_residual,
        "material_domain_status": "VERIFIED_SUBSET",
        "material_domain_relation": domain_relation,
        "native_solver_claim_preserved": result_record.get("final_geometry_verified"),
    }


def _make_triangulation_configuration(
    extension, *, icosphere_subdivisions: int,
    maximum_icosphere_subdivisions: int = -1,
    maximum_radius_doublings: int = 6,
    paper_projective_hyperboloid: bool = False,
    algebraic_finite_map: bool = False,
):
    if type(algebraic_finite_map) is not bool:
        raise ValueError('algebraic_finite_map must be a boolean')
    _validate_paper_source_configuration(paper_projective_hyperboloid,
        icosphere_subdivisions, maximum_icosphere_subdivisions)
    if isinstance(icosphere_subdivisions, bool) or not isinstance(
        icosphere_subdivisions, int
    ) or not 0 <= icosphere_subdivisions <= 7:
        raise ValueError("icosphere subdivisions must be an integer in [0, 7]")
    if (isinstance(maximum_icosphere_subdivisions, bool) or
            not isinstance(maximum_icosphere_subdivisions, int) or
            maximum_icosphere_subdivisions not in range(-1, 8) or
            (maximum_icosphere_subdivisions != -1 and
             maximum_icosphere_subdivisions < icosphere_subdivisions)):
        raise ValueError(
            "maximum icosphere subdivisions must be -1 or an integer in "
            "[icosphere_subdivisions, 7]"
        )
    if isinstance(maximum_radius_doublings, bool) or not isinstance(
        maximum_radius_doublings, int
    ) or not 0 <= maximum_radius_doublings <= 16:
        raise ValueError("maximum radius doublings must be an integer in [0, 16]")
    source_settings = extension.SpinTriangulationSettings()
    source_settings.icosphere_subdivisions = icosphere_subdivisions
    source_settings.maximum_icosphere_subdivisions = maximum_icosphere_subdivisions
    source_settings.maximum_radius_doublings = maximum_radius_doublings
    if hasattr(source_settings, 'paper_projective_hyperboloid'):
        source_settings.paper_projective_hyperboloid = paper_projective_hyperboloid
    elif paper_projective_hyperboloid:
        raise RuntimeError('native paper-projective source support is unavailable')
    if hasattr(source_settings, 'algebraic_finite_map'):
        source_settings.algebraic_finite_map = algebraic_finite_map
    elif algebraic_finite_map:
        raise RuntimeError('native algebraic finite-map support is unavailable')
    evaluator_settings = spin.SpinEvaluatorConfig(
        derivative_method=('DISCRETE_FULL_RANK_QUADRIC' if paper_projective_hyperboloid or algebraic_finite_map
                           else 'DISCRETE_ELLIPSOID'),
        icosphere_subdivisions=icosphere_subdivisions,
        maximum_icosphere_subdivisions=maximum_icosphere_subdivisions,
        maximum_radius_doublings=maximum_radius_doublings,
        paper_projective_hyperboloid=paper_projective_hyperboloid,
        algebraic_finite_map=algebraic_finite_map,
    )
    return source_settings, evaluator_settings


def _validate_paper_source_configuration(enabled, initial, maximum):
    if type(enabled) is not bool:
        raise ValueError('paper_projective_hyperboloid must be a boolean')
    if enabled and maximum not in (-1, initial):
        raise ValueError('paper source requires fixed normal discretization')


def _solver_initial_quadrics(
    mode: str,
    fixture_coefficients,
    *,
    replay_coefficients=None,
) -> tuple[tuple[float, ...], ...]:
    """Return explicit starts, or an empty sequence to use the native seeded protocol."""
    coefficients = np.asarray(fixture_coefficients, dtype=np.float64)
    if coefficients.shape != (10,) or not np.all(np.isfinite(coefficients)):
        raise ValueError("fixture initialization must be a finite ten-coefficient quadric")
    if mode == "fixed_fixture":
        return (tuple(float(value) for value in coefficients),)
    if mode == "deterministic_mode_seed":
        # The C++ solver generates starts from mode + SpinSolverConfig.deterministic_seed.
        return ()
    if mode in {"exact_line_search_replay", "exact_newton_state_replay"}:
        if replay_coefficients is None:
            raise ValueError("diagnostic replay requires explicit replay coefficients")
        replay = np.asarray(replay_coefficients, dtype=np.float64)
        if replay.shape != (10,) or not np.all(np.isfinite(replay)):
            raise ValueError("replay coefficients must be a finite ten-coefficient quadric")
        norm = float(np.linalg.norm(replay))
        if not np.isfinite(norm) or abs(norm - 1.0) > 1e-12:
            raise ValueError("replay coefficients must have unit norm")
        return (tuple(float(value) for value in replay),)
    raise ValueError(f"unsupported initialization mode: {mode}")


def _quadric_coefficients_after_coordinate_scale(
    coefficients: object, coordinate_scale: float,
) -> np.ndarray:
    """Pull a normalized quadric through ``x_new = scale * x_old``.

    The coefficient convention is constant, three linear terms, then six
    quadratic terms. The returned vector is normalized because quadric scale
    does not change its zero set.
    """
    values = np.asarray(coefficients, dtype=np.float64)
    scale = float(coordinate_scale)
    if values.shape != (10,) or not np.isfinite(values).all():
        raise ValueError("quadric coefficients must be a finite ten-vector")
    if not math.isfinite(scale) or scale <= 0.0:
        raise ValueError("coordinate scale must be finite and positive")
    transformed = values * np.asarray(
        (1.0, *(1.0 / scale for _ in range(3)),
         *((1.0 / (scale * scale)) for _ in range(6))),
        dtype=np.float64,
    )
    norm = float(np.linalg.norm(transformed))
    if not math.isfinite(norm) or norm <= 0.0:
        raise ValueError("scaled quadric coefficients must have finite nonzero norm")
    return transformed / norm


def _opposite_boolean_operation(operation: str) -> str:
    """Return the other side of a closed quadric in the same design domain."""
    opposite = {"intersection": "subtraction", "subtraction": "intersection"}
    if operation not in opposite:
        raise ValueError("unsupported Cherchi Boolean side: " + str(operation))
    return opposite[operation]


def _synthetic_locked_problem(outer_vertices, outer_faces, design_vertices, design_faces):
    """Create a fixed, nonzero-density shell plus an editable cube domain.

    The shell is the exact partition remainder between the outer and inner
    cubes.  It is intentionally part of the locked input, not an inferred
    moment offset, so the Cherchi bridge must integrate and bind real fixed
    geometry.
    """
    outer_vertices = np.asarray(outer_vertices, dtype=np.float64)
    outer_faces = np.asarray(outer_faces, dtype=np.int32)
    design_vertices = np.asarray(design_vertices, dtype=np.float64)
    design_faces = np.asarray(design_faces, dtype=np.int32)
    fixed_vertices = np.vstack((outer_vertices, design_vertices))
    fixed_faces = np.vstack((
        outer_faces,
        design_faces[:, ::-1] + len(outer_vertices),
    ))
    return SpinArtifactProblem(
        original=(outer_vertices, outer_faces),
        fixed_regions=(("fixed-shell", (fixed_vertices, fixed_faces), 0.75),),
        design_regions=(("editable-domain", (design_vertices, design_faces), 1.0),),
        length_unit="m",
        mass_unit="kg",
        target={"mode": "YO_YO", "axis": (0., 0., 1.), "origin": (0., 0., 0.),
                "support_points": ()},
        minimum_outer_wall=0.1,
    )


def _mesh_vertices_in_target_frame(vertices, request):
    """Transform normalized model-frame vertices into the request target frame."""
    points = np.asarray(vertices, dtype=np.float64)
    origin = np.asarray(request.target_origin, dtype=np.float64)
    axes = np.asarray(request.target_axes, dtype=np.float64)
    if points.ndim != 2 or points.shape[1] != 3 or \
            origin.shape != (3,) or axes.shape != (3, 3) or \
            not np.isfinite(points).all() or not np.isfinite(origin).all() or \
            not np.isfinite(axes).all():
        raise ValueError("target-frame mesh transform received invalid coordinates")
    return (points - origin) @ axes


def _mesh_vertices_from_target_frame(vertices, request):
    """Invert the locked target-frame transform for normalized model points."""
    points = np.asarray(vertices, dtype=np.float64)
    origin = np.asarray(request.target_origin, dtype=np.float64)
    axes = np.asarray(request.target_axes, dtype=np.float64)
    if points.ndim != 2 or points.shape[1] != 3 or \
            origin.shape != (3,) or axes.shape != (3, 3) or \
            not np.isfinite(points).all() or not np.isfinite(origin).all() or \
            not np.isfinite(axes).all():
        raise ValueError("inverse target-frame mesh transform received invalid coordinates")
    if not np.allclose(axes.T @ axes, np.eye(3), rtol=0.0, atol=1e-10) or \
            abs(float(np.linalg.det(axes)) - 1.0) > 1e-10:
        raise ValueError("inverse target-frame transform requires a right-handed orthonormal frame")
    return points @ axes.T + origin


def _normalized_mesh_vertices_in_target_frame(vertices, normalization, request):
    """Map world/model geometry through locked normalization and target frame."""
    points = np.asarray(vertices, dtype=np.float64)
    origin = np.asarray(normalization.get("origin"), dtype=np.float64)
    length = float(normalization.get("L", 0.0))
    if points.ndim != 2 or points.shape[1] != 3 or origin.shape != (3,) or \
            not np.isfinite(points).all() or not np.isfinite(origin).all() or \
            not math.isfinite(length) or length <= 0.0:
        raise ValueError("normalized target-frame transform received invalid geometry or scale")
    return _mesh_vertices_in_target_frame((points - origin) / length, request)


def _load_problem_manifest(path: Path):
    """Load a caller-supplied, provenance-pinned locked geometry problem.

    Manifest paths are resolved relative to the manifest file. Geometry is
    parsed without repair; ``SpinArtifactProblem`` performs the exact
    partition and material validation before a solver request is built.
    """
    manifest_path = Path(path).resolve(strict=True)
    if not manifest_path.is_file():
        raise ValueError("problem manifest must be a regular file")
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except Exception as error:
        raise ValueError(f"invalid problem manifest JSON: {error}") from error
    required = {
        "schema", "outer_mesh", "fixed_regions", "design_regions",
        "length_unit", "mass_unit", "target", "minimum_outer_wall",
        "preflight_quadric_target_normalized",
    }
    if not isinstance(manifest, dict) or set(manifest) != required:
        missing = sorted(required - set(manifest)) if isinstance(manifest, dict) else sorted(required)
        extra = sorted(set(manifest) - required) if isinstance(manifest, dict) else []
        raise ValueError(f"problem manifest fields mismatch; missing={missing}, extra={extra}")
    if manifest["schema"] != "faster-cherchi-problem-manifest-v1":
        raise ValueError("unsupported problem manifest schema")
    preflight = np.asarray(
        manifest["preflight_quadric_target_normalized"], dtype=np.float64,
    )
    if preflight.shape != (10,) or not np.isfinite(preflight).all() or \
            abs(float(np.linalg.norm(preflight)) - 1.0) > 1e-12:
        raise ValueError(
            "preflight_quadric_target_normalized must be a finite unit ten-vector"
        )
    asset_paths = {manifest_path}

    def load_mesh(raw_path, label):
        if not isinstance(raw_path, str) or not raw_path.strip():
            raise ValueError(f"{label} mesh path must be a nonempty string")
        mesh_path = Path(raw_path)
        if not mesh_path.is_absolute():
            mesh_path = manifest_path.parent / mesh_path
        mesh_path = mesh_path.resolve(strict=True)
        if not mesh_path.is_file():
            raise ValueError(f"{label} mesh path must be a regular file")
        asset_paths.add(mesh_path)
        vertices, faces = _read_obj(mesh_path)
        return vertices, faces

    def load_regions(raw_regions, label):
        if not isinstance(raw_regions, list):
            raise ValueError(f"{label} must be an explicit list")
        result = []
        for index, region in enumerate(raw_regions):
            if not isinstance(region, dict) or set(region) != {"id", "mesh", "density"}:
                raise ValueError(f"{label}[{index}] must contain exactly id, mesh, density")
            result.append((
                region["id"], load_mesh(region["mesh"], f"{label}[{index}]"),
                region["density"],
            ))
        return result

    if not isinstance(manifest["target"], dict):
        raise ValueError("target must be an object")
    target_keys = {"mode", "axis", "origin", "support_points"}
    if set(manifest["target"]) != target_keys:
        raise ValueError("target must contain exactly mode, axis, origin, support_points")
    problem = SpinArtifactProblem(
        original=load_mesh(manifest["outer_mesh"], "outer"),
        fixed_regions=load_regions(manifest["fixed_regions"], "fixed_regions"),
        design_regions=load_regions(manifest["design_regions"], "design_regions"),
        length_unit=manifest["length_unit"],
        mass_unit=manifest["mass_unit"],
        target=manifest["target"],
        minimum_outer_wall=manifest["minimum_outer_wall"],
    )
    return problem, tuple(sorted(asset_paths)), tuple(float(x) for x in preflight)


def _preflight_and_solver_initial_quadrics(
    mode: str,
    fixture_coefficients,
    *,
    replay_coefficients=None,
) -> tuple[tuple[float, ...], tuple[tuple[float, ...], ...]]:
    """Keep the frozen backend-parity probe independent from the solver start."""
    fixture = np.asarray(fixture_coefficients, dtype=np.float64)
    if fixture.shape != (10,) or not np.all(np.isfinite(fixture)):
        raise ValueError("preflight fixture must be a finite ten-coefficient quadric")
    preflight = tuple(float(value) for value in fixture)
    starts = _solver_initial_quadrics(
        mode, fixture, replay_coefficients=replay_coefficients,
    )
    return preflight, starts


def _argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--build", type=Path)
    parser.add_argument("--build-receipt", type=Path,
                        help="observed native build receipt; defaults to build/build-receipt.json when present")
    parser.add_argument("--cherchi-probe", type=Path)
    parser.add_argument("--resume-terminal-audit", action="store_true",
                        help="audit a persisted checkpoint only; does not complete an interrupted run")
    parser.add_argument("--resume-finalization", action="store_true",
                        help="finalize checked run context and saved bundle without rerunning Newton")
    parser.add_argument("--defer-finalization", action="store_true",
                        help="persist solver/bundle state, leaving final audits and evidence for a separate job")
    parser.add_argument("--problem-manifest", type=Path,
                        help="locked arbitrary geometry manifest; omitting uses the synthetic regression fixture")
    parser.add_argument("--maximum-function-evaluations", type=int, default=60)
    parser.add_argument("--maximum-restarts", type=int, default=1)
    parser.add_argument("--icosphere-subdivisions", type=int, default=5)
    parser.add_argument("--maximum-icosphere-subdivisions", type=int, default=-1)
    parser.add_argument('--paper-projective-hyperboloid', action='store_true',
                        help='opt in to fixed-normal paper projective hyperboloid triangles; legacy default unchanged')
    parser.add_argument('--algebraic-finite-map', action='store_true',
                        help='opt in to the algebraic full-rank finite-point map; legacy default unchanged')
    parser.add_argument(
        "--maximum-radius-doublings", type=int, choices=range(17), default=6,
        help="maximum certified hyperboloid-cap radius doublings (strict profile default: 6)",
    )
    parser.add_argument("--cherchi-timeout-seconds", type=float, default=180.0)
    parser.add_argument("--center-of-mass-weight", type=float, default=0.0,
                        help="normalized spinning-top objective weight; other modes require zero")
    parser.add_argument("--spinning-top-objective", default="INERTIA_AND_LOWERING",
                        choices=("INERTIA_AND_LOWERING", "LOWER_CENTER_OF_MASS"))
    parser.add_argument('--directional-profile',default=None,choices=('PURE_HEIGHT_ACTIVE_STRATUM_V1',))
    parser.add_argument('--paraboloid-grid-resolution',type=int,default=8)
    parser.add_argument('--fixed-normal-paraboloid',action='store_true',
        help='explicit bounded finite-normal reduced source; not complete ideal-point support')
    parser.add_argument('--truncation-radius',type=float,default=4.)
    parser.add_argument(
        "--initialization-mode",
        choices=("fixed_fixture", "deterministic_mode_seed",
                 "exact_line_search_replay", "exact_newton_state_replay", "supported_mode_seed"),
        default="fixed_fixture",
    )
    parser.add_argument(
        "--initial-quadric-coefficients", type=float, nargs=10, default=None,
        metavar=("C0", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9"),
    )
    parser.add_argument("--deterministic-seed", type=int, default=0)
    return parser


def run(output: Path, *, build: Path, probe: Path,
        build_receipt: Path | None = None,
        problem_manifest: Path | None = None,
        maximum_function_evaluations: int = 60,
        maximum_restarts: int = 1,
        icosphere_subdivisions: int = 5,
        maximum_icosphere_subdivisions: int = -1,
        maximum_radius_doublings: int = 6,
        paper_projective_hyperboloid: bool = False,
        algebraic_finite_map: bool = False,
        initialization_mode: str = "fixed_fixture",
        initial_quadric_coefficients=None,
        deterministic_seed: int = 0,
        center_of_mass_weight: float = 0.0,
        spinning_top_objective: str = "INERTIA_AND_LOWERING",
        directional_profile: str | None = None,
        paraboloid_grid_resolution: int = 8,
        fixed_normal_paraboloid: bool = False,
        truncation_radius: float = 4.0,
        cherchi_timeout_seconds: float = 180.0,
        defer_finalization: bool = False) -> dict[str, object]:
    if type(algebraic_finite_map) is not bool:
        raise ValueError('algebraic_finite_map must be a boolean')
    if type(fixed_normal_paraboloid) is not bool or (fixed_normal_paraboloid and (
            directional_profile!='PURE_HEIGHT_ACTIVE_STRATUM_V1' or
            maximum_icosphere_subdivisions not in (-1,icosphere_subdivisions))):
        raise ValueError('fixed-normal paraboloid requires an explicit fixed-resolution directional profile')
    _validate_paper_source_configuration(paper_projective_hyperboloid,
        icosphere_subdivisions, maximum_icosphere_subdivisions)
    if isinstance(center_of_mass_weight, bool) or not isinstance(center_of_mass_weight, (int, float)) or \
            not math.isfinite(center_of_mass_weight) or center_of_mass_weight < 0.0:
        raise ValueError("center_of_mass_weight must be finite and nonnegative")
    center_of_mass_weight = float(center_of_mass_weight)
    mode = (json.loads(Path(problem_manifest).read_text())["target"]["mode"]
            if problem_manifest is not None else "YO_YO")
    support=None
    if directional_profile is not None:
        if directional_profile!='PURE_HEIGHT_ACTIVE_STRATUM_V1' or mode!='SPINNING_TOP' or \
                spinning_top_objective!='LOWER_CENTER_OF_MASS' or center_of_mass_weight!=0. or \
                initialization_mode!='supported_mode_seed' or paper_projective_hyperboloid or algebraic_finite_map:
            raise ValueError('pure-height directional profile requires top endpoint, supported seeds and unmixed stratum source')
        if type(paraboloid_grid_resolution) is not int or paraboloid_grid_resolution<1 or \
                type(truncation_radius) not in (int,float) or not math.isfinite(truncation_radius) or truncation_radius<=0:
            raise ValueError('pure-height grid/radius must be fixed positive finite settings')
        support=np.eye(10)[:,[0,1,2,3,7,8]]
    elif initialization_mode=='supported_mode_seed':
        raise ValueError('supported seeds require an explicit pure-height directional profile')
    objective_definition = rotational_objective_definition(
        mode, center_of_mass_weight, spinning_top_objective)
    if type(defer_finalization) is not bool:
        raise ValueError("defer finalization must be a boolean execution choice")
    if not 1 <= maximum_function_evaluations <= 60:
        raise ValueError("strict Eq.14 diagnostic budget must be in [1, 60]")
    if isinstance(maximum_restarts, bool) or not isinstance(maximum_restarts, int) or \
            not 1 <= maximum_restarts <= np.iinfo(np.int32).max:
        raise ValueError("maximum restarts must be a positive signed 32-bit integer")
    if not np.isfinite(cherchi_timeout_seconds) or not 1.0 <= cherchi_timeout_seconds <= 300.0:
        raise ValueError("Cherchi timeout must be finite and in [1, 300] seconds")
    if isinstance(deterministic_seed, bool) or not isinstance(deterministic_seed, int):
        raise ValueError("deterministic seed must be an integer")
    if initialization_mode not in (
        "fixed_fixture", "deterministic_mode_seed", "exact_line_search_replay",
        "exact_newton_state_replay",
        "supported_mode_seed",
    ):
        raise ValueError(f"unsupported initialization mode: {initialization_mode}")
    if initialization_mode in {"exact_line_search_replay", "exact_newton_state_replay"}:
        if initial_quadric_coefficients is None:
            raise ValueError("diagnostic replay requires explicit replay coefficients")
    elif initial_quadric_coefficients is not None:
        raise ValueError(
            "initial quadric coefficients are only valid for exact line-search replay"
        )
    output = output.absolute()
    if problem_manifest is not None:
        feasibility = verify_axis_aligned_yoyo_moment_capacity(problem_manifest)
        if feasibility.status == "INFEASIBLE_CERTIFIED":
            output.mkdir(parents=False, exist_ok=False)
            report_path = output / "preflight-feasibility-report.json"
            report_path.write_text(
                json.dumps(feasibility.to_dict(), sort_keys=True,
                           separators=(",", ":"), allow_nan=False) + "\n",
                encoding="utf-8",
            )
            result_record = {
                "schema": "faster-cherchi-newton-preflight-result-v1",
                "classification": "RECONSTRUCTED_METHOD_VALIDATION",
                "status": "INFEASIBLE_CERTIFIED",
                "solver_status": "NOT_RUN_PREFLIGHT_CERTIFIED",
                "diagnostic": "continuous moment-capacity certificate proves the requested center-of-mass constraint is infeasible",
                "function_evaluations": 0,
                "restart_count": 0,
                "wall_seconds": 0.0,
                "preflight_report": str(report_path),
                "input_hash": feasibility.input_hash,
                "state": "CONTINUOUS_GEOMETRIC_INFEASIBILITY_CERTIFIED_BEFORE_NATIVE_SOLVER",
            }
            (output / "result.json").write_text(
                json.dumps(result_record, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n",
                encoding="utf-8",
            )
            return result_record
    build = build.resolve(strict=True)
    output.mkdir(parents=False, exist_ok=False)
    extension_dir = (build / "python").resolve(strict=True)
    candidates = sorted(extension_dir.glob("_make_it_stand*.so"))
    if len(candidates) != 1:
        raise ValueError("build must contain exactly one native extension")
    sys.path.insert(0, str(extension_dir))
    try:
        extension = __import__("_make_it_stand")
    finally:
        sys.path.pop(0)
    if Path(extension.__file__).resolve() != candidates[0].resolve():
        raise ValueError("loaded native extension differs from declared build")
    build_configuration = native_build_configuration_files(
        build, extension.__file__, receipt=build_receipt)
    probe = probe.resolve(strict=True)
    if not probe.is_file() or not probe.stat().st_mode & 0o111:
        raise ValueError("pinned Cherchi probe is unavailable or non-executable")
    differential_settings, evaluator_settings = _make_triangulation_configuration(
        extension, icosphere_subdivisions=icosphere_subdivisions,
        maximum_icosphere_subdivisions=maximum_icosphere_subdivisions,
        maximum_radius_doublings=maximum_radius_doublings,
        paper_projective_hyperboloid=paper_projective_hyperboloid,
        algebraic_finite_map=algebraic_finite_map,
    )
    if support is not None:
        differential_settings.paraboloid_grid_resolution=paraboloid_grid_resolution
        differential_settings.truncation_radius=truncation_radius
        if fixed_normal_paraboloid:
            if not hasattr(differential_settings,'fixed_normal_paraboloid'):
                raise RuntimeError('native fixed-normal paraboloid support is unavailable')
            differential_settings.fixed_normal_paraboloid=True
        evaluator_settings=replace(evaluator_settings,derivative_method='DISCRETE_REDUCED_STRATUM',
            paraboloid_grid_resolution=paraboloid_grid_resolution,truncation_radius=truncation_radius,
            fixed_normal_paraboloid=fixed_normal_paraboloid)

    manifest_asset_paths: tuple[Path, ...] = ()
    manifest_preflight = None
    if problem_manifest is None:
        outer_vertices, outer_faces = (
            np.asarray(value) for value in (
                ((-2., -2., -2.), (2., -2., -2.), (2., 2., -2.), (-2., 2., -2.),
                 (-2., -2., 2.), (2., -2., 2.), (2., 2., 2.), (-2., 2., 2.)),
                ((0, 2, 1), (0, 3, 2), (4, 5, 6), (4, 6, 7), (0, 1, 5), (0, 5, 4),
                 (1, 2, 6), (1, 6, 5), (2, 3, 7), (2, 7, 6), (3, 0, 4), (3, 4, 7)),
            )
        )
        design_vertices = outer_vertices / 2.0
        problem = _synthetic_locked_problem(
            outer_vertices, outer_faces, design_vertices, outer_faces,
        )
    else:
        problem, manifest_asset_paths, manifest_preflight = _load_problem_manifest(
            problem_manifest,
        )
    solver_settings = spin.SpinSolverConfig(
        maximum_iterations=maximum_function_evaluations,
        maximum_restarts=maximum_restarts,
        maximum_function_evaluations_per_restart=maximum_function_evaluations,
        equation_profile="DUAL_MONOMIAL_RAW_V1",
        deterministic_seed=deterministic_seed,
        capture_evaluation_payload=True,
    )
    problem, request, problem_identity = _prepare(
        problem, evaluator_settings, solver_settings, allow_reduced_stratum=support is not None,
    )
    if request.mode != 'SPINNING_TOP' and center_of_mass_weight != 0.0:
        raise ValueError('center_of_mass_weight is supported only for spinning tops')
    request = replace(request, center_of_mass_weight=center_of_mass_weight,
                      spinning_top_objective=spinning_top_objective)
    coordinate_scale = 1.0 / float(problem.normalization["L"])
    fixture_coefficients = (
        np.asarray(manifest_preflight, dtype=np.float64)
        if manifest_preflight is not None else
        _quadric_coefficients_after_coordinate_scale(
            np.asarray(_FIXTURES[0][1], dtype=np.float64), coordinate_scale,
        )
    )
    replay_coefficients = initial_quadric_coefficients
    if replay_coefficients is not None:
        replay_coefficients = (
            np.asarray(replay_coefficients, dtype=np.float64)
            if manifest_preflight is not None else
            _quadric_coefficients_after_coordinate_scale(
                replay_coefficients, coordinate_scale,
            )
        )
    if support is None:
        preflight_coefficients, initial_quadrics = _preflight_and_solver_initial_quadrics(
            initialization_mode, fixture_coefficients,replay_coefficients=replay_coefficients)
    else:
        preflight_coefficients=tuple(float(x) for x in fixture_coefficients)
        initial_quadrics=_pure_height_initial_quadrics(deterministic_seed,maximum_restarts)
    classification = (
        "DIAGNOSTIC_REPLAY"
        if initialization_mode in {"exact_line_search_replay", "exact_newton_state_replay"}
        else "RECONSTRUCTED_METHOD_VALIDATION"
        if problem_manifest is not None else "SYNTHETIC_REGRESSION"
    )
    # The preflight point is fixed for backend parity; it is not necessarily the
    # solver start under native deterministic initialization.
    request_record = {
        "schema": "faster-cherchi-newton-integration-input-v7",
        "classification": classification,
        "execution_protocol": "STAGED_FINALIZATION_V1" if defer_finalization else "SINGLE_JOB",
        "mode": request.mode,
        "center_of_mass_weight": center_of_mass_weight,
        "spinning_top_objective": spinning_top_objective,
        "objective_definition": objective_definition,
        "initialization_mode": initialization_mode,
        "deterministic_seed": deterministic_seed,
        "initial_coefficients": list(initial_quadrics[0])
        if initial_quadrics else None,
        "preflight_coefficients": list(preflight_coefficients),
        "seed_protocol": (
            ("explicit preflight quadric frozen in the geometry manifest; no target solution read"
             if problem_manifest is not None else
             "explicit frozen regression fixture; no target solution read")
            if initialization_mode == "fixed_fixture"
            else (
                "captured frozen line-search trial supplied explicitly; diagnostic replay only; no target solution read"
            if initialization_mode == "exact_line_search_replay"
                else "frozen accepted Newton state supplied explicitly; diagnostic replay only; no target solution read"
                if initialization_mode == "exact_newton_state_replay"
                else "native C++ deterministicInitialQuadrics(mode, settings) using deterministic_seed; no target solution read"
            )
        ),
        "maximum_iterations": request.settings.maximum_iterations,
        "maximum_restarts": request.settings.maximum_restarts,
        "maximum_function_evaluations_per_restart": request.settings.maximum_function_evaluations_per_restart,
        "equation_profile": request.settings.equation_profile,
        "capture_evaluation_payload": bool(request.settings.capture_evaluation_payload),
        "quadric_icosphere_subdivisions": icosphere_subdivisions,
        "maximum_icosphere_subdivisions": maximum_icosphere_subdivisions,
        "maximum_radius_doublings": int(maximum_radius_doublings),
        "paper_projective_hyperboloid": paper_projective_hyperboloid,
        "algebraic_finite_map": algebraic_finite_map,
        "locked_problem_digest": problem.input_digest,
        "locked_problem_identity": problem.identity,
        "problem_manifest": (
            None if problem_manifest is None else {
                "path": str(Path(problem_manifest).resolve(strict=True)),
                "sha256": _sha(Path(problem_manifest).resolve(strict=True)),
                "preflight_quadric_target_normalized": list(manifest_preflight),
            }
        ),
        "normalization": dict(problem.normalization),
        "fixed_regions": [region.identity for region in problem.fixed_regions],
        "design_regions": [region.identity for region in problem.design_regions],
        "fixed_shell_moments_input_normalized": list(request.fixed_shell_moments),
        "quadric_coordinate_scale": (
            1.0 if manifest_preflight is not None else coordinate_scale
        ),
        "quadric_coordinate_frame": (
            "MANIFEST_TARGET_NORMALIZED" if manifest_preflight is not None
            else "LEGACY_FIXTURE_INPUT_NORMALIZED"
        ),
        "cherchi_boolean_timeout_seconds": cherchi_timeout_seconds,
        "quadric_branch_policy": "real full-rank ellipsoid/hyperboloid; use centroid sign only when winding verifies the centroid is inside, otherwise require unanimous outward interface-normal/quadric-gradient side witness; ambiguous strata rejected",
        "moment_map_binding": "faster-newton-map-v3-exact-geometry: geometry/topology/operation/quadric/backend/derivative/source/cap/exact rational geometry identity/float derivation audit/raw ten moments/full 10x10 native-coefficient Jacobian",
        "cherchi_probe_sha256": _sha(probe),
        "native_extension_sha256": _sha(candidates[0]),
        "runner_sha256": _sha(Path(__file__).resolve()),
        "orientation_normalizer_sha256": _sha(
            Path(normalize_closed_component_orientation.__code__.co_filename).resolve()
        ),
        "spin_python_api_sha256": _sha(Path(spin.__file__).resolve()),
    }
    if support is not None:
        request_record.update(
            directional_profile=directional_profile,support_basis=support.tolist(),
            initial_quadrics=[list(row) for row in initial_quadrics],
            paraboloid_grid_resolution=paraboloid_grid_resolution,truncation_radius=truncation_radius,
            seed_protocol='PCG64_STANDARD_NORMAL_DUAL_SUPPORT_0_1_2_3_7_8_V1; native cross entries halved then unit normalized; no fixture/answer input',
            quadric_branch_policy='regular fixed-grid paraboloid stratum only; complete active directions required; unsupported trials rejected',
            moment_map_binding='faster-newton-map-v4-directional-exact-geometry; genuine10x5 Jacobian plus native basis and dual support',
            paper_source_status='METHOD_FIXED_GRID_NOT_PAPER_FIXED_NORMAL')
        if fixed_normal_paraboloid:
            request_record.update(fixed_normal_paraboloid=True,
                source_map='BOUNDED_FINITE_NORMAL_PARABOLOID_V1',
                quadric_branch_policy='regular bounded finite-normal paraboloid patch; domain coverage and remote closure required; unsupported trials rejected',
                paper_source_status='METHOD_BOUNDED_FINITE_NORMAL_NOT_COMPLETE_IDEAL_MAP')
    (output / "inputs.json").write_text(
        json.dumps(request_record, sort_keys=True, separators=(",", ":")) + "\n",
        encoding="utf-8",
    )
    design_path = output / "design.obj"
    target_design_vertices = _mesh_vertices_in_target_frame(
        request.design_vertices, request,
    )
    _write_obj(design_path, target_design_vertices, request.design_faces)
    dependency_root = ROOT / "external/reference-sources-20260917/sources/InteractiveAndRobustMeshBooleans-7bd6c2697695f555f9a15d7e154e084d841316f0"
    if not dependency_root.is_dir():
        raise ValueError("pinned Cherchi dependency source tree is unavailable")
    evidence_inputs = begin_run(
        repository_root=ROOT,
        runner_files=[Path(__file__), Path(__file__).with_name("run_faster_reference_boolean.py")],
        checker_files=[
            ROOT / "python/make_it_stand/verification/run_evidence.py",
            ROOT / "python/make_it_stand/verification/cherchi_motion_differential.py",
            ROOT / "python/make_it_stand/verification/mass_moments.py",
            ROOT / "python/make_it_stand/verification/rational_text.py",
            ROOT / "python/make_it_stand/verification/faster_material.py",
            ROOT / "python/make_it_stand/verification/spin_problem.py",
            ROOT / "python/make_it_stand/verification/spin_artifact.py",
            ROOT / "python/make_it_stand/verification/spin_bundle.py",
            ROOT / "experiments/cherchi_reference_probe/probe.cpp",
            ROOT / "experiments/cherchi_reference_probe/CMakeLists.txt",
            ROOT / "src/spin/QuadricTriangulator.cpp",
            ROOT / "src/spin/QuadricDiscreteMomentDifferential.cpp",
            ROOT / "src/spin/QuadricMomentEvaluator.cpp",
            ROOT / "bindings/pybind_module.cpp",
        ],
        asset_files=[output / "inputs.json", design_path, probe,
                     *manifest_asset_paths],
        build_configuration_files=[ROOT / "CMakeLists.txt", *build_configuration],
        parameters=request_record,
        seed_protocol={
            "seed": deterministic_seed,
            "protocol": str(request_record["seed_protocol"]),
        },
        dependency_roots={"cherchi_reference": dependency_root},
        native_modules=[extension],
        pure_python=False,
        source_policy=EXPLICIT_EFFECTIVE_SOURCE_POLICY,
    )
    (output / "begin.json").write_text(
        json.dumps(evidence_inputs.as_dict(), sort_keys=True, separators=(",", ":")) + "\n",
        encoding="utf-8",
    )
    evaluation_journal_path = output / "evaluation-events.jsonl"
    _initialize_evaluation_journal(evaluation_journal_path)
    stage_timing_path = output / "stage-timings.jsonl"
    _initialize_stage_timing_journal(stage_timing_path)
    solver_event_path = output / "solver-events.jsonl"
    _initialize_solver_event_journal(
        solver_event_path, inputs_sha256=_sha(output / "inputs.json"),
    )
    evaluations: list[dict[str, object]] = []
    evaluator_call_count = 0
    evaluation_operations: dict[int, str] = {}

    # One CGAL preflight checks that the frozen initial source/operation matches
    # the existing evaluator. All in-loop moments and derivatives below are
    # then computed from Cherchi only; CGAL is not mixed into the callback map.
    with _diagnostic_stage(stage_timing_path, -1, "cgal_preflight_geometry"):
        initial_native = spin.evaluate_spin_it_faster_quadric(
            request, preflight_coefficients,
            evaluator_settings=evaluator_settings,
        )
    with _diagnostic_stage(stage_timing_path, -1, "preflight_source_motion"):
        initial_source = (_pure_height_source_motion(extension,preflight_coefficients,
            target_design_vertices,differential_settings) if support is not None else
            extension.evaluate_spin_it_faster_source_motion(
                preflight_coefficients, target_design_vertices, differential_settings))
    if not initial_source["defined"] or \
            not np.array_equal(np.asarray(initial_source["faces"]), initial_native.interface_faces) or \
        not np.allclose(initial_source["vertices"], initial_native.interface_vertices,
                        rtol=0.0, atol=2e-13):
        raise ValueError("initial source-motion mesh differs from native preflight quadric interface")
    _source_derivative_branch(initial_source,
        paper_projective_hyperboloid=paper_projective_hyperboloid,
        algebraic_finite_map=algebraic_finite_map)
    if support is None:
        _supported_quadric_operation(preflight_coefficients)
    initial_operation = str(_source_side_selection_details(initial_source)["operation"])
    if initial_operation != (
        "intersection" if initial_native.quadric_solid_is_positive_region else "subtraction"
    ):
        raise ValueError("closed-source selector disagrees with initial native Boolean operation")

    def evaluate(q_values: tuple[float, ...]) -> dict[str, object]:
        nonlocal evaluator_call_count
        _raise_if_interrupt_pending()
        index = evaluator_call_count
        evaluator_call_count += 1
        q = np.asarray(q_values, dtype=np.float64)
        _append_solver_journal_event(solver_event_path, "geometry_evaluation_started", {
            "evaluation_index": index, "quadric_coefficients": q.tolist(),
        })
        try:
            evaluated = evaluate_success(q, index)
        except DurableJournalError:
            # Never turn a failed persistence attempt into a scientific
            # evaluation failure: the current run must stop without a possibly
            # duplicated or silently missing ledger row.
            raise
        except RunInterrupted:
            raise
        except Exception as error:
            # Keep failed trials alongside successful geometries. The solver's
            # terminal diagnostic alone loses the trial index and coefficients.
            record = _failed_evaluation_record(
                index, q, error, artifact_directory=output,
                operation=evaluation_operations.get(index),
            )
            _append_evaluation_journal(evaluation_journal_path, record)
            evaluations.append(record)
            _append_solver_journal_event(solver_event_path, "geometry_evaluation_finished", {
                "evaluation_index": index, "outcome": "failed",
                "diagnostic": f"{type(error).__name__}: {error}",
            })
            raise
        _append_solver_journal_event(solver_event_path, "geometry_evaluation_finished", {
            "evaluation_index": index, "outcome": "returned",
        })
        return evaluated

    def persist_solver_event(event: dict[str, object]) -> None:
        _raise_if_interrupt_pending()
        _append_solver_journal_event(solver_event_path, "solver_event", event)

    def evaluate_success(q: np.ndarray, index: int) -> dict[str, object]:
        with _diagnostic_stage(stage_timing_path, index, "quadric_preflight"):
            if support is None:
                _supported_quadric_operation(q)
        with _diagnostic_stage(stage_timing_path, index, "source_motion"):
            source = (_pure_height_source_motion(extension,q,target_design_vertices,differential_settings)
                if support is not None else extension.evaluate_spin_it_faster_source_motion(
                    q,target_design_vertices,differential_settings))
        if not source["defined"]:
            raise ValueError("quadric source motion unavailable: " + str(source["diagnostic"]))
        derivative_branch, source_features = _source_derivative_branch(source,
            paper_projective_hyperboloid=paper_projective_hyperboloid,
            algebraic_finite_map=algebraic_finite_map)
        source_vertices = np.asarray(source["vertices"], dtype=np.float64)
        source_faces = np.asarray(source["faces"], dtype=np.int64)
        selector = _source_side_selection_details(source)
        operation = str(selector["operation"])
        evaluation_operations[index] = operation
        dimension=10 if support is None else support.shape[1]-1
        source_motion = np.asarray(source["vertex_velocity"], dtype=np.float64).reshape(-1, 3, dimension)
        directional_binding=({} if support is None else {
            'native_tangent_basis':source['native_tangent_basis'], 'support_basis':support})
        source_path = output / f"source-{index:03d}.obj"
        material_path = output / f"material-{index:03d}.obj"
        with _diagnostic_stage(stage_timing_path, index, "source_obj_write"):
            _write_obj(source_path, source_vertices, source_faces)
        with _diagnostic_stage(stage_timing_path, index, "cherchi_boolean"):
            completed = subprocess.run(
                [str(probe), operation, str(design_path), str(source_path),
                 str(material_path), "--dump-exact-rational"],
                cwd=output, capture_output=True, text=True, check=False,
                timeout=cherchi_timeout_seconds,
            )
        if completed.returncode != 0:
            raise ValueError("Cherchi Boolean failed: " + completed.stderr.strip())
        with _diagnostic_stage(stage_timing_path, index, "output_mesh_read_and_orientation"):
            material_vertices, material_faces = _read_obj(material_path)
        # Cherchi may emit a geometrically closed, orientable manifold whose
        # face windings are locally inconsistent. Preserve its original OBJ,
        # normalize only triangle winding, and retain the exact flip witness.
        # The normalizer fails closed for open/nonmanifold/nonorientable input.
        with _diagnostic_stage(stage_timing_path, index, "orientation_normalization_and_write"):
            exact_snapshot, exact_float_derivation_audit, orientation_audit = \
                _prepare_exact_motion_material(material_vertices, material_faces,
                    Path(str(material_path) + ".exact.tsv"))
            exact_vertices, oriented_faces = exact_snapshot.arrays()
            # Preserve authenticated raw coordinates in the derived OBJ.
            # Exact values and derivatives below use the rational snapshot,
            # not these possibly collapsed display coordinates.
            oriented_path = output / f"material-{index:03d}-oriented.obj"
            _write_obj(oriented_path, material_vertices, oriented_faces)
        with _diagnostic_stage(stage_timing_path, index, "exact_boolean_rational_snapshot"):
            exact_geometry_identity = digest({
                "rational_geometry": exact_snapshot.identity,
                "float_derivation_audit": exact_float_derivation_audit,
            })
        with _diagnostic_stage(stage_timing_path, index, "provenance_reconstruction"):
            construction = _summarize_vertex_construction_provenance(
                (design_path, source_path), material_path,
                Path(str(material_path) + ".construction.tsv"),
            )
        with _diagnostic_stage(stage_timing_path, index, "motion_recovery"):
            output_motion = recover_construction_vertex_motion(
                construction["records"], target_design_vertices, source_vertices,
                source_motion, exact_vertices,
            )
        with _diagnostic_stage(stage_timing_path, index, "discrete_moment_jacobian"):
            jacobian = differentiate_polyhedral_raw_moments(
                exact_vertices, oriented_faces, output_motion,
                density=request.material_density,
                rational_vertices=exact_snapshot.rational_vertices,
            )
        with _diagnostic_stage(stage_timing_path, index, "primary_moment_integration"):
            raw_moments = integrate_rational_polyhedral_raw_moments(
                exact_vertices, oriented_faces, density=request.material_density,
                rational_vertices=exact_snapshot.rational_vertices)
        design_domain = np.asarray(request.design_vertices, dtype=np.float64)
        design_root_side = float(np.max(np.ptp(design_domain, axis=0)))
        with _diagnostic_stage(stage_timing_path, index, "independent_moment_integration_and_audit"):
            # Signed tetrahedra are independent of the primary divergence /
            # surface integral and differentiate the same rational material.
            independent_raw_moments = np.asarray(
                _raw_moments(exact_snapshot, request.material_density), dtype=np.float64)
            independent_moment_audit = _audit_independent_moment_agreement(
                raw_moments, independent_raw_moments,
                density=request.material_density, root_side=design_root_side,
            )
        construction_topology_hash = _construction_topology_signature(
            construction, oriented_faces,
        )
        topology_hash = _material_topology_identity(
            operation, construction, oriented_faces,
            derivative_branch=derivative_branch,
        )
        boolean_backend_identity = (
            "cherchi-reference-7bd6c2697695f555f9a15d7e154e084d841316f0:"
            + _sha(probe)
        )
        derivative_identity = (
            "cherchi-provenance-rational-reynolds-surface-raw-moments-v2"
        )
        source_identity = hashlib.sha256(
            (_sha(source_path) + "\0" + str(source["triangulation_hash"])).encode("ascii")
        ).hexdigest()
        cap_payload = {
            **derivative_branch,
            "maximum_radius_doublings": int(
                differential_settings.maximum_radius_doublings
            ),
            "triangulation_hash": str(source["triangulation_hash"]),
        }
        cap_identity = hashlib.sha256(json.dumps(
            cap_payload, sort_keys=True, separators=(",", ":"), allow_nan=False,
        ).encode("ascii")).hexdigest()
        quadric_hash = str(extension.spin_it_faster_quadric_hash(q))
        moment_map_sha256 = _moment_map_binding_sha256(
            _sha(oriented_path), topology_hash, operation, raw_moments, jacobian,
            quadric_coefficients=q,
            backend_identity=boolean_backend_identity,
            derivative_identity=derivative_identity,
            source_identity=source_identity,
            cap_identity=cap_identity,
            exact_geometry_identity=exact_geometry_identity,
            **directional_binding,
        )
        cherchi_geometry_hash = _sha(material_path)
        geometry_hash = _sha(oriented_path)
        record = {
            "evaluation_index": index,
            "solver_evaluation_index": index + 1,
            "quadric_coefficients": q.tolist(),
            "operation": operation,
            "boolean_operation": operation,
            "closed_source_selector": selector,
            "geometry_hash": geometry_hash,
            "cherchi_geometry_hash": cherchi_geometry_hash,
            "construction_topology_hash": construction_topology_hash,
            "topology_hash": topology_hash,
            "derivative_branch": derivative_branch,
            **({'source_feature_payload': source_features} if source_features is not None else {}),
            "source_geometry_hash": _sha(source_path),
            "source_triangulation_hash": str(source["triangulation_hash"]),
            "source_surface_type": str(source["surface_type"]),
            "boolean_backend_identity": boolean_backend_identity,
            "derivative_identity": derivative_identity,
            "source_identity": source_identity,
            "cap_identity": cap_identity,
            "exact_geometry_identity": exact_geometry_identity,
            "exact_float_derivation_audit": exact_float_derivation_audit,
            "evaluation_map_identity": moment_map_sha256,
            "quadric_hash": quadric_hash,
            "cap_identity_payload": cap_payload,
            "source_truncation_radius": float(source["truncation_radius"]),
            "source_initial_icosphere_subdivisions": int(
                source["initial_icosphere_subdivisions"]
            ),
            "source_effective_icosphere_subdivisions": int(
                source["effective_icosphere_subdivisions"]
            ),
            "material_vertex_count": int(len(material_vertices)),
            "material_face_count": int(len(material_faces)),
            "oriented_output_path": str(oriented_path),
            "orientation_normalization": orientation_audit,
            "jacobian_shape": list(jacobian.shape),
            "moment_map": {
                "raw_moments": raw_moments.tolist(),
                "moment_jacobian": jacobian.tolist(),
                **{key:np.asarray(value).tolist() for key,value in directional_binding.items()},
                "sha256": moment_map_sha256,
            },
            "independent_moment_audit": independent_moment_audit,
            "output_path": str(material_path),
            "stdout": completed.stdout.strip(),
            "stderr": completed.stderr.strip(),
        }
        _append_evaluation_journal(evaluation_journal_path, record)
        evaluations.append(record)
        return {
            "variable_vertices": material_vertices,
            "variable_vertices_model": _mesh_vertices_from_target_frame(
                material_vertices, request),
            "variable_rational_vertices": (exact_snapshot.rational_vertices
                                           if exact_snapshot is not None else ()),
            "variable_exact_snapshot": exact_snapshot,
            "exact_geometry_identity": exact_geometry_identity,
            "exact_float_derivation_audit": exact_float_derivation_audit,
            "variable_faces": oriented_faces,
            "variable_raw_moments": raw_moments,
            "quadric_coefficients": q,
            "boolean_operation": operation,
            "quadric_source_vertices": source_vertices,
            "quadric_source_faces": source_faces,
            "quadric_source_identity": source_identity,
            "raw_moments": raw_moments,
            **({'moment_jacobian':jacobian} if support is None else {
                'native_tangent_basis':source['native_tangent_basis'],
                'native_directional_jacobian':jacobian}),
            "derivative_defined": True,
            "derivative_diagnostic": "Cherchi provenance + explicit closed-component orientation normalization + source motion + Reynolds transport",
            "geometry_hash": geometry_hash,
            "topology_hash": topology_hash,
            "moment_map_sha256": moment_map_sha256,
            "evaluation_map_identity": moment_map_sha256,
            "quadric_hash": quadric_hash,
            "boolean_backend_identity": boolean_backend_identity,
            "derivative_identity": derivative_identity,
            "source_identity": source_identity,
            "cap_identity": cap_identity,
        }

    def complement_evaluator(retained_raw: dict[str, object], context: dict[str, object]):
        design = output / "design.obj"
        source_path = output / "terminal-quadric-source.obj"
        carved_path = output / "terminal-complement.obj"
        source_vertices = np.asarray(retained_raw["quadric_source_vertices"], dtype=np.float64)
        source_faces = np.asarray(retained_raw["quadric_source_faces"], dtype=np.int32)
        source_identity = str(retained_raw["quadric_source_identity"])
        if not source_identity or source_vertices.ndim != 2 or source_vertices.shape[1:] != (3,) or \
                source_faces.ndim != 2 or source_faces.shape[1:] != (3,):
            raise ValueError("terminal Newton material omitted its bound quadric source geometry")
        _write_obj(source_path, source_vertices, source_faces)
        source_file_sha256 = _sha(source_path)
        opposite_operation = _opposite_boolean_operation(
            str(retained_raw["boolean_operation"]),
        )
        with _diagnostic_stage(stage_timing_path, evaluator_call_count, "terminal_complement_boolean"):
            completed = subprocess.run(
                [str(probe), opposite_operation, str(design), str(source_path),
                 str(carved_path), "--dump-exact-rational"],
                cwd=output, capture_output=True, text=True, check=False,
                timeout=cherchi_timeout_seconds,
            )
        if completed.returncode != 0:
            raise ValueError("Cherchi terminal complement failed: " + completed.stderr.strip())
        vertices, faces = _read_obj(carved_path)
        if len(vertices):
            exact_snapshot, exact_float_derivation_audit, orientation, moments = _terminal_exact_material(
                vertices, faces, Path(str(carved_path) + ".exact.tsv"), request.material_density,
            )
        else:
            exact_snapshot = None
            orientation = {'authority': 'EXACT_EMPTY_MATERIAL', 'face_order_changed': False}
            moments = np.zeros(10, dtype=np.float64)
            exact_float_derivation_audit = {
                "schema": "exact-boolean-float-derivation-audit-v1",
                "status": "EXPLICIT_EMPTY_MATERIAL",
                "vertex_count": 0,
                "approximation_status_vertex_count": 0,
                "coordinate_scale": 1.0,
                "maximum_absolute_error": 0.0,
                "absolute_error_bound": 64.0 * np.finfo(np.float64).eps,
                "sidecar_sha256": _sha(Path(str(carved_path) + ".exact.tsv")),
            }
        oriented_path = output / "terminal-complement-oriented.obj"
        _write_obj(oriented_path, vertices, faces)
        exact_geometry_identity = digest({
            "rational_geometry": exact_snapshot.identity if exact_snapshot is not None
            else {"explicit_empty_material": True},
            "float_derivation_audit": exact_float_derivation_audit,
        })
        topology = hashlib.sha256(
            ("cherchi-terminal-complement-topology-v1\0" +
             opposite_operation + "\0").encode("ascii") +
            np.asarray(faces, dtype="<i4").tobytes(order="C")
        ).hexdigest()
        backend = (
            "cherchi-reference-7bd6c2697695f555f9a15d7e154e084d841316f0:"
            + _sha(probe)
        )
        geometry_hash = _sha(oriented_path)
        derivative_identity = "NOT_USED_TERMINAL_COMPLEMENT_NO_NEWTON_DERIVATIVE"
        map_identity = hashlib.sha256(json.dumps({
            "geometry": geometry_hash,
            "topology": topology,
            "backend": backend,
            "operation": opposite_operation,
            "quadric_source_identity": source_identity,
            "moments": moments.tolist(),
            "exact_geometry_identity": exact_geometry_identity,
        }, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("ascii")).hexdigest()
        complement_record = {
            "operation": opposite_operation,
            "construction": "same-triangulated-quadric-opposite-side; exact partition independently checked by SealedSpinArtifact",
            "retained_operation": retained_raw["boolean_operation"],
            "same_quadric_source_identity": source_identity,
            "quadric_source_file_sha256": source_file_sha256,
            "backend": backend,
            "geometry_hash": geometry_hash,
            "topology_hash": topology,
            "moment_map_identity": map_identity,
            "exact_geometry_identity": exact_geometry_identity,
            "exact_float_derivation_audit": exact_float_derivation_audit,
            "moments": moments.tolist(),
            "vertices": int(len(vertices)),
            "faces": int(len(faces)),
            "orientation_normalization": orientation,
            "stdout": completed.stdout.strip(),
            "stderr": completed.stderr.strip(),
        }
        (output / "terminal-complement.json").write_text(
            json.dumps(_json_safe(complement_record), sort_keys=True,
                       separators=(",", ":"), allow_nan=False) + "\n",
            encoding="utf-8",
        )
        return {
            "variable_vertices": vertices,
            "variable_vertices_model": _mesh_vertices_from_target_frame(
                vertices, request),
            "variable_rational_vertices": (exact_snapshot.rational_vertices
                                           if exact_snapshot is not None else ()),
            "variable_exact_snapshot": exact_snapshot,
            "exact_geometry_identity": exact_geometry_identity,
            "exact_float_derivation_audit": exact_float_derivation_audit,
            "variable_faces": faces,
            "variable_raw_moments": moments,
            "boolean_operation": opposite_operation,
            "quadric_source_identity": source_identity,
            "quadric_source_vertices": source_vertices,
            "quadric_source_faces": source_faces,
            **({'moment_jacobian':np.zeros((10,10),dtype=np.float64)} if support is None else {}),
            "derivative_defined": False,
            "derivative_diagnostic": derivative_identity,
            "geometry_hash": geometry_hash,
            "topology_hash": topology,
            "exact_geometry_identity": exact_geometry_identity,
            "evaluation_map_identity": map_identity,
            "boolean_backend_identity": backend,
            "derivative_identity": derivative_identity,
        }

    started = time.perf_counter()
    bridge_result = _run_diagnostic_stage(
        stage_timing_path, -1, "newton_bridge_and_geometry_callbacks",
        lambda: solve_cherchi_geometry_and_seal_faster(
            problem,
            geometry_evaluator=lambda coefficients, _context: evaluate(coefficients),
            complement_evaluator=complement_evaluator,
            initial_quadrics=initial_quadrics,
            evaluator_settings=evaluator_settings,
            solver_settings=solver_settings,
            center_of_mass_weight=center_of_mass_weight,
            spinning_top_objective=spinning_top_objective,
            **({} if support is None else {'support_basis':support}),
            event_sink=persist_solver_event,
            stage_observer=lambda stage, operation: _run_diagnostic_stage(
                stage_timing_path, evaluator_call_count, stage, operation,
            ),
        ),
    )
    bridge_identity = bridge_result.identity
    result = bridge_identity.get("native_run")
    if not isinstance(result, dict):
        raise DurableJournalError("locked Cherchi bridge omitted the native Newton record")
    solver_event_audit = _audit_streamed_solver_events(
        solver_event_path, list(result["events"]),
    )
    sealed_artifact = bridge_result.artifact
    artifact_bundle = None
    artifact_reload_digest = None
    journal_records = _read_evaluation_journal(evaluation_journal_path)
    if len(journal_records) != len(evaluations) or [
        record["evaluation_index"] for record in journal_records
    ] != [record["evaluation_index"] for record in evaluations]:
        raise DurableJournalError(
            "completed solver ledger differs from durable per-evaluation journal"
        )
    try:
        map_binding_audit = _audit_newton_evaluation_bindings(
            evaluations, list(result["events"]),
            quadric_hash=extension.spin_it_faster_quadric_hash,
        )
    except Exception as error:
        map_binding_audit = {
            "status": "FAILED",
            "diagnostic": f"{type(error).__name__}: {error}",
        }
    terminal_audit_payload = {
        "solver_status": result["status"],
        "final_geometry_verified": result["final_geometry_verified"],
        "evaluation_map_binding_audit": map_binding_audit,
        "trace": result["trace"],
        "evaluations": evaluations,
        "design_domain_vertices": target_design_vertices,
        "design_domain_faces": request.design_faces,
    }
    finalization_context = {
        "schema": "faster-run-finalization-context-v1",
        "output_directory": str(output),
        "begin": evidence_inputs.as_dict(),
        "native_result": result,
        "classification": classification,
        "initialization_mode": initialization_mode,
        "deterministic_seed": deterministic_seed,
        "solver_initial_quadric_count": len(initial_quadrics),
        "locked_problem_digest": problem.input_digest,
        "bridge_status": bridge_result.status,
        "bridge_scientific_digest": bridge_result.scientific_digest,
        "bridge_identity": bridge_identity,
        "bridge_seal_error": bridge_result.seal_error,
        "expected_artifact_digest": None if sealed_artifact is None else sealed_artifact.scientific_digest,
        "solver_wall_seconds": time.perf_counter() - started,
    }
    terminal_audit_payload["run_finalization"] = finalization_context
    _save_terminal_audit_checkpoint(
        output, terminal_audit_payload,
        dependencies=[path for path in output.rglob("*") if path.is_file()
                      and path.name != "stage-timings.jsonl"],
    )
    if sealed_artifact is not None:
        artifact_bundle = output / "sealed-material-artifact"
        reloaded_artifact = save_spin_artifact(sealed_artifact, artifact_bundle)
        artifact_reload_digest = reloaded_artifact.scientific_digest
    if defer_finalization:
        return {
            "classification": classification,
            "status": "PENDING_FINALIZATION",
            "solver_status": result["status"],
            "diagnostic": "solver output persisted; final geometry and run evidence not yet complete",
            "function_evaluations": result["function_evaluations"],
            "restart_count": result["restart_count"],
            "wall_seconds": time.perf_counter() - started,
            "state": "AWAITING_INDEPENDENT_RUN_FINALIZATION",
        }
    return _finalize_run(
        output, payload=terminal_audit_payload, context=finalization_context,
        evidence_inputs=evidence_inputs, sealed_artifact=sealed_artifact,
        artifact_bundle=artifact_bundle, artifact_reload_digest=artifact_reload_digest,
        started=started, prior_wall_seconds=0.0, execution_protocol="SINGLE_JOB",
    )


def _bridge_identity_attachment(context: dict) -> dict:
    """Preserve the actual bridge identity; never reconstruct a legacy one.

    The digest alone cannot rehydrate a FasterMaterialResult for CAD. Old
    checkpoints remain readable but explicitly lack this evidence attachment.
    This binding validates persistence, not convergence or CAD suitability.
    """
    if "bridge_identity" not in context:
        return {"schema": "cherchi-bridge-identity-attachment-v1",
                "status": "LEGACY_IDENTITY_NOT_SAVED", "identity": None,
                "digest": context["bridge_scientific_digest"]}
    identity = context["bridge_identity"]
    required = {"operation", "native_run", "status", "seal_error",
                "problem_digest", "artifact_scientific_digest"}
    if not isinstance(identity, dict) or not required <= identity.keys():
        raise ValueError("bridge identity must be the original complete mapping")
    if (digest(material_run_scientific_identity(identity)) != context["bridge_scientific_digest"] or
            identity.get("operation") != "SOLVE_CHERCHI_EXTERNAL_GEOMETRY" or
            identity.get("native_run") != context["native_result"] or
            identity.get("status") != context["bridge_status"] or
            identity.get("seal_error") != context["bridge_seal_error"] or
            identity.get("problem_digest") != context["locked_problem_digest"] or
            identity.get("artifact_scientific_digest") != context["expected_artifact_digest"]):
        raise ValueError("bridge identity differs from original solver/problem/material binding")
    return {"schema": "cherchi-bridge-identity-attachment-v1",
            "status": "BOUND_IDENTITY", "identity": identity,
            "digest": context["bridge_scientific_digest"]}


def _finalize_run(output: Path, *, payload: dict, context: dict,
                  evidence_inputs: RunInputsV2, sealed_artifact,
                  artifact_bundle: Path | None, artifact_reload_digest: str | None,
                  started: float, prior_wall_seconds: float,
                  execution_protocol: str) -> dict:
    bridge_attachment = _bridge_identity_attachment(context)
    result = context["native_result"]
    evaluations = payload["evaluations"]
    extension = spin._native_module()
    evaluation_journal_path = output / "evaluation-events.jsonl"
    journal_records = _read_evaluation_journal(evaluation_journal_path)
    if _json_safe(journal_records) != _json_safe(evaluations):
        raise DurableJournalError("finalization evaluation ledger differs from its persisted journal")
    solver_event_audit = _audit_streamed_solver_events(
        output / "solver-events.jsonl", list(result["events"]))
    try:
        map_binding_audit = _audit_newton_evaluation_bindings(
            evaluations, list(result["events"]),
            quadric_hash=extension.spin_it_faster_quadric_hash)
    except Exception as error:
        map_binding_audit = {
            "status": "FAILED", "diagnostic": f"{type(error).__name__}: {error}"}
    terminal_audit_payload = {**payload, "evaluation_map_binding_audit": map_binding_audit}
    try:
        terminal_geometry_audit = _audit_terminal_geometry_binding(
            terminal_audit_payload,
            quadric_hash=extension.spin_it_faster_quadric_hash,
        )
    except Exception as error:
        terminal_geometry_audit = {
            "status": "FAILED",
            "geometry_verified": False,
            "diagnostic": f"{type(error).__name__}: {error}",
            "material_domain_status": "NOT_VERIFIED",
        }
    run_status = result["status"] if map_binding_audit["status"] == "CONSISTENT" else \
        "EVALUATION_MAP_BINDING_FAILED"
    if run_status == "CONVERGED" and context["bridge_status"] != "CONVERGED":
        run_status = context["bridge_status"]
    if run_status == "CONVERGED" and terminal_geometry_audit["status"] != "PASSED":
        run_status = "FINAL_GEOMETRY_AUDIT_FAILED"
    if run_status == "CONVERGED" and (
            sealed_artifact is None or artifact_reload_digest != sealed_artifact.scientific_digest):
        run_status = "SEALED_ARTIFACT_RELOAD_FAILED"
    result_record = {
        "schema": "faster-cherchi-newton-integration-result-v7",
        "classification": context["classification"],
        "initialization_mode": context["initialization_mode"],
        "deterministic_seed": context["deterministic_seed"],
        "solver_initial_quadric_count": context["solver_initial_quadric_count"],
        "execution_protocol": execution_protocol,
        "status": run_status,
        "solver_status": result["status"],
        "locked_problem_digest": context["locked_problem_digest"],
        "bridge_status": context["bridge_status"],
        "bridge_scientific_digest": context["bridge_scientific_digest"],
        "bridge_seal_error": context["bridge_seal_error"],
        "bridge_identity_attachment": bridge_attachment,
        "sealed_artifact": {
            "bundle_path": None if artifact_bundle is None else str(artifact_bundle),
            "scientific_digest": artifact_reload_digest,
            "problem_digest": (
                None if sealed_artifact is None else sealed_artifact.problem.input_digest
            ),
            "fixed_region_count": (
                None if sealed_artifact is None else len(sealed_artifact.problem.fixed_regions)
            ),
            "facts": None if sealed_artifact is None else sealed_artifact.facts,
        },
        "diagnostic": result["diagnostic"],
        "scope": result["scope"],
        "final_geometry_verified": result["final_geometry_verified"],
        "function_evaluations": result["function_evaluations"],
        "restart_count": result["restart_count"],
        "restart_function_evaluations": result["restart_function_evaluations"],
        "event_stream_version": result.get("event_stream_version"),
        "equation_profile": result["equation_profile"],
        "normalized_equation_residual": result["normalized_equation_residual"],
        "event_hash": result["event_hash"],
        "canonical_event_hash": _canonical_event_digest(result["events"]),
        "canonical_event_hash_algorithm": (
            "sha256-json-sorted-python-traceback-frame-line-normalized-v1"
        ),
        "scientific_event_hash_v2": _scientific_event_digest_v2(result["events"]),
        "scientific_event_hash_algorithm_v2": (
            "sha256-json-sorted-python-traceback-lines-normalized-"
            "evaluation-payload-durations-omitted-v2"
        ),
        "trace_hash": result["trace_hash"],
        "trace": result["trace"],
        "events": result["events"],
        "evaluations": evaluations,
        "evaluation_journal": {
            "path": str(evaluation_journal_path),
            "record_count": len(evaluations),
            "sha256": _sha(evaluation_journal_path),
        },
        "streamed_solver_event_audit": solver_event_audit,
        "evaluation_map_binding_audit": map_binding_audit,
        "external_terminal_geometry_audit": terminal_geometry_audit,
        "wall_seconds": prior_wall_seconds + time.perf_counter() - started,
        "wall_seconds_scope": (
            "current_job" if execution_protocol == "SINGLE_JOB" else
            "completed_solver_checkpoint_plus_current_finalization; prior post-checkpoint work excluded"
        ),
        "state": (
            "LOCKED_PROBLEM_SEALED_AND_INDEPENDENT_AUDITS_PASSED"
            if run_status == "CONVERGED" else
            "FINAL_GEOMETRY_AUDIT_FAILED"
            if result["status"] == "CONVERGED" else
            "NUMERICAL_RESULT_REQUIRES_INDEPENDENT_FINAL_GEOMETRY_AUDIT"
        ),
    }
    # Revalidate all checkpoint dependencies after potentially expensive checks.
    _load_terminal_audit_checkpoint(output / "terminal-audit-checkpoint.json")
    (output / "result.json").write_text(
        json.dumps(_json_safe(result_record), sort_keys=True,
                   separators=(",", ":"), allow_nan=False) + "\n",
        encoding="utf-8",
    )
    output_files = sorted(
        path for path in output.rglob("*")
        if path.is_file() and path.name not in {
            "evidence.json", "freshness.json", "stage-timings.jsonl",
        }
    )
    evidence = finish_run(
        evidence_inputs,
        events=_evidence_events(result_record),
        output_files=output_files,
        measurements={
            "wall_seconds": float(result_record["wall_seconds"]),
            "function_evaluations": int(result_record["function_evaluations"]),
            "solver_status": str(result_record["solver_status"]),
            "execution_protocol": execution_protocol,
        },
    )
    (output / "evidence.json").write_text(evidence.to_json() + "\n", encoding="utf-8")
    freshness = check_freshness(evidence)
    (output / "freshness.json").write_text(
        json.dumps({"status": freshness.status, "mismatches": list(freshness.mismatches),
                    "contextual_changes": list(freshness.contextual_changes)},
                   sort_keys=True, separators=(",", ":")) + "\n",
        encoding="utf-8",
    )
    if freshness.status != "CURRENT":
        raise ValueError("full-run finalization evidence is " + freshness.status)
    return result_record


def main() -> None:
    args = _argument_parser().parse_args()
    if args.resume_terminal_audit and args.resume_finalization:
        _argument_parser().error("choose only one resume operation")
    if args.resume_finalization:
        result = resume_run_finalization(args.output)
        print(json.dumps({key: result[key] for key in (
            "classification", "status", "execution_protocol", "wall_seconds", "state")},
            sort_keys=True))
        return
    if args.resume_terminal_audit:
        result = resume_terminal_audit(args.output)
        print(json.dumps({key: result[key] for key in (
            "classification", "status", "complete_run_evidence", "wall_seconds")},
            sort_keys=True))
        return
    if args.build is None or args.cherchi_probe is None:
        _argument_parser().error("--build and --cherchi-probe are required for a new solve")
    previous_sigterm_handler = signal.getsignal(signal.SIGTERM)

    def terminate_with_receipt(signal_number, _frame):
        global _PENDING_INTERRUPT_SIGNAL
        _PENDING_INTERRUPT_SIGNAL = int(signal_number)
        raise RunInterrupted(signal_number)

    global _PENDING_INTERRUPT_SIGNAL
    _PENDING_INTERRUPT_SIGNAL = None
    signal.signal(signal.SIGTERM, terminate_with_receipt)
    try:
        result = run(args.output, build=args.build, probe=args.cherchi_probe,
                     build_receipt=args.build_receipt,
                     problem_manifest=args.problem_manifest,
                     maximum_function_evaluations=args.maximum_function_evaluations,
                     maximum_restarts=args.maximum_restarts,
                     icosphere_subdivisions=args.icosphere_subdivisions,
                     maximum_icosphere_subdivisions=args.maximum_icosphere_subdivisions,
                     maximum_radius_doublings=args.maximum_radius_doublings,
                     paper_projective_hyperboloid=args.paper_projective_hyperboloid,
                     algebraic_finite_map=args.algebraic_finite_map,
                     initialization_mode=args.initialization_mode,
                     initial_quadric_coefficients=args.initial_quadric_coefficients,
                     deterministic_seed=args.deterministic_seed,
                     center_of_mass_weight=args.center_of_mass_weight,
                     spinning_top_objective=args.spinning_top_objective,
                     directional_profile=args.directional_profile,
                     paraboloid_grid_resolution=args.paraboloid_grid_resolution,
                     fixed_normal_paraboloid=args.fixed_normal_paraboloid,
                     truncation_radius=args.truncation_radius,
                     cherchi_timeout_seconds=args.cherchi_timeout_seconds,
                     defer_finalization=args.defer_finalization)
    except RunInterrupted as interruption:
        receipt = _write_interruption_receipt(args.output, interruption.signal_number)
        print(json.dumps({
            "status": receipt["status"],
            "termination_signal": receipt["termination_signal"],
            "partial_evaluation_count": receipt["partial_evaluation_count"],
            "solver_terminal_status": receipt["solver_terminal_status"],
        }, sort_keys=True), file=sys.stderr)
        raise SystemExit(128 + interruption.signal_number)
    except KeyboardInterrupt:
        receipt = _write_interruption_receipt(args.output, signal.SIGINT)
        print(json.dumps({
            "status": receipt["status"],
            "termination_signal": receipt["termination_signal"],
            "partial_evaluation_count": receipt["partial_evaluation_count"],
            "solver_terminal_status": receipt["solver_terminal_status"],
        }, sort_keys=True), file=sys.stderr)
        raise SystemExit(128 + signal.SIGINT)
    finally:
        signal.signal(signal.SIGTERM, previous_sigterm_handler)
        _PENDING_INTERRUPT_SIGNAL = None
    print(json.dumps({key: result[key] for key in (
        "classification", "status", "diagnostic", "function_evaluations",
        "restart_count", "wall_seconds", "state",
    )}, sort_keys=True))


if __name__ == "__main__":
    main()
