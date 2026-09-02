"""Post-flatten geometry and mass-moment seal consumed by simulation tools."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math
import os
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from .io import write_json_atomic, write_text_atomic
from .paper_compat import SupportFlattening, _faces, _vertices
from .paper_evidence import EvidenceClass


def _canonical_hash(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(
            value,
            allow_nan=False,
            separators=(",", ":"),
            sort_keys=True,
        ).encode("utf-8")
    ).hexdigest()


def _array_hash(*arrays: np.ndarray) -> str:
    digest = hashlib.sha256()
    for array in arrays:
        contiguous = np.ascontiguousarray(array)
        digest.update(str(contiguous.dtype).encode("ascii"))
        digest.update(json.dumps(list(contiguous.shape)).encode("ascii"))
        digest.update(contiguous.tobytes(order="C"))
    return digest.hexdigest()


def _obj_text(vertices: np.ndarray, faces: np.ndarray) -> str:
    lines = [
        f"v {point[0]:.17g} {point[1]:.17g} {point[2]:.17g}"
        for point in vertices
    ]
    lines.extend(
        f"f {face[0] + 1} {face[1] + 1} {face[2] + 1}"
        for face in faces
    )
    return "\n".join(lines) + "\n"


def _surface_moments(
    vertices: np.ndarray,
    faces: np.ndarray,
) -> tuple[float, np.ndarray, np.ndarray]:
    if faces.shape[0] == 0:
        return 0.0, np.zeros(3), np.zeros((3, 3))
    triangles = vertices[faces]
    reference = vertices.mean(axis=0)
    first = triangles[:, 0] - reference
    second = triangles[:, 1] - reference
    third = triangles[:, 2] - reference
    signed_mass = np.einsum(
        "ij,ij->i",
        first,
        np.cross(second, third),
    ) / 6.0
    references = np.broadcast_to(reference, (triangles.shape[0], 1, 3))
    tetrahedra = np.concatenate((references, triangles), axis=1)
    sums = tetrahedra.sum(axis=1)
    self_outer = np.einsum("mvi,mvj->mij", tetrahedra, tetrahedra)
    second_moments = signed_mass[:, None, None] * (
        np.einsum("mi,mj->mij", sums, sums) + self_outer
    ) / 20.0
    return (
        float(signed_mass.sum()),
        np.einsum("m,mi->i", signed_mass, sums) / 4.0,
        second_moments.sum(axis=0),
    )


def _compact_referenced_vertices(
    vertices: np.ndarray,
    faces: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Discard grid corners that are not part of the extracted boundary."""

    if faces.shape[0] == 0:
        return np.empty((0, 3), dtype=np.float64), np.empty((0, 3), dtype=np.int32)
    referenced = np.unique(faces.reshape(-1))
    remap = np.full(vertices.shape[0], -1, dtype=np.int64)
    remap[referenced] = np.arange(referenced.size, dtype=np.int64)
    compact_faces = remap[faces]
    if np.any(compact_faces < 0):
        raise ValueError("mesh compaction produced an invalid face index")
    return (
        np.ascontiguousarray(vertices[referenced], dtype=np.float64),
        np.ascontiguousarray(compact_faces, dtype=np.int32),
    )


def _triangulate_simple_polygon(
    vertices: np.ndarray,
    loop: Sequence[int],
    normal: Sequence[float],
) -> np.ndarray:
    """Triangulate one planar simple loop with deterministic ear clipping."""

    indices = [int(value) for value in loop]
    if len(indices) < 3 or len(set(indices)) != len(indices):
        raise ValueError("planar patch loop must contain at least three unique vertices")
    if min(indices) < 0 or max(indices) >= vertices.shape[0]:
        raise ValueError("planar patch loop index is outside the vertex array")
    direction = np.asarray(normal, dtype=np.float64)
    if direction.shape != (3,) or not np.isfinite(direction).all():
        raise ValueError("planar patch normal must contain three finite values")
    length = float(np.linalg.norm(direction))
    if length <= 1.0e-15:
        raise ValueError("planar patch normal must be nonzero")
    direction /= length
    reference = np.asarray((1.0, 0.0, 0.0), dtype=np.float64)
    if abs(float(direction @ reference)) > 0.9:
        reference = np.asarray((0.0, 1.0, 0.0), dtype=np.float64)
    first_axis = np.cross(direction, reference)
    first_axis /= np.linalg.norm(first_axis)
    second_axis = np.cross(direction, first_axis)
    points = vertices[indices]
    projected = np.column_stack((points @ first_axis, points @ second_axis))
    extent = float(np.linalg.norm(np.ptp(projected, axis=0)))
    tolerance = max(1.0, extent * extent) * 256.0 * np.finfo(np.float64).eps

    def cross2(first: np.ndarray, second: np.ndarray, third: np.ndarray) -> float:
        return float(
            (second[0] - first[0]) * (third[1] - first[1])
            - (second[1] - first[1]) * (third[0] - first[0])
        )

    signed_twice_area = sum(
        cross2(
            np.zeros(2),
            projected[index],
            projected[(index + 1) % len(indices)],
        )
        for index in range(len(indices))
    )
    if abs(signed_twice_area) <= tolerance:
        raise ValueError("planar patch loop has zero projected area")
    if signed_twice_area < 0.0:
        indices.reverse()
        projected = projected[::-1].copy()

    remaining = list(range(len(indices)))
    triangles: list[tuple[int, int, int]] = []
    while len(remaining) > 3:
        accepted = False
        for ordinal, current in enumerate(remaining):
            previous = remaining[ordinal - 1]
            following = remaining[(ordinal + 1) % len(remaining)]
            a, b, c = projected[[previous, current, following]]
            if cross2(a, b, c) <= tolerance:
                continue
            contains_vertex = False
            for candidate in remaining:
                if candidate in (previous, current, following):
                    continue
                point = projected[candidate]
                if (
                    cross2(a, b, point) >= -tolerance
                    and cross2(b, c, point) >= -tolerance
                    and cross2(c, a, point) >= -tolerance
                ):
                    contains_vertex = True
                    break
            if contains_vertex:
                continue
            triangles.append(
                (indices[previous], indices[current], indices[following])
            )
            del remaining[ordinal]
            accepted = True
            break
        if not accepted:
            raise ValueError(
                "planar support patch is not a triangulable simple polygon"
            )
    triangles.append(tuple(indices[index] for index in remaining))
    return np.ascontiguousarray(triangles, dtype=np.int32)


def _flattened_patch_loops(
    faces: np.ndarray,
    patch_face_indices: np.ndarray,
) -> list[tuple[int, ...]]:
    edge_uses: dict[tuple[int, int], int] = {}
    for face in faces[patch_face_indices]:
        values = tuple(int(value) for value in face)
        for first, second in (
            (values[0], values[1]),
            (values[1], values[2]),
            (values[2], values[0]),
        ):
            edge = tuple(sorted((first, second)))
            edge_uses[edge] = edge_uses.get(edge, 0) + 1
    if any(count not in (1, 2) for count in edge_uses.values()):
        raise ValueError("flattened support patch is not edge-manifold")
    boundary = {edge for edge, count in edge_uses.items() if count == 1}
    adjacency: dict[int, list[int]] = {}
    for first, second in boundary:
        adjacency.setdefault(first, []).append(second)
        adjacency.setdefault(second, []).append(first)
    if not boundary or any(len(neighbors) != 2 for neighbors in adjacency.values()):
        raise ValueError("flattened support patch boundary is not a union of loops")
    for neighbors in adjacency.values():
        neighbors.sort()

    unused = set(boundary)
    loops: list[tuple[int, ...]] = []
    while unused:
        start, current = min(unused)
        loop = [start]
        previous = start
        unused.remove((start, current))
        while current != start:
            if current in loop:
                raise ValueError("flattened support patch boundary self-repeats")
            loop.append(current)
            candidates = [
                neighbor
                for neighbor in adjacency[current]
                if tuple(sorted((current, neighbor))) in unused
            ]
            if not candidates:
                raise ValueError("flattened support patch boundary is open")
            following = min(candidates)
            edge = tuple(sorted((current, following)))
            unused.remove(edge)
            previous, current = current, following
        if len(loop) < 3:
            raise ValueError("flattened support patch boundary is degenerate")
        loops.append(tuple(loop))
    return sorted(loops)


def _remesh_flattened_support_patches(
    vertices: np.ndarray,
    faces: np.ndarray,
    flattenings: Sequence[SupportFlattening],
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    """Replace collapsed support triangles by one triangulation per boundary loop."""

    current_faces = np.ascontiguousarray(faces, dtype=np.int32)
    patches: list[dict[str, Any]] = []
    for flattening_index, flattening in enumerate(flattenings):
        flat_vertices = set(int(value) for value in flattening.vertex_indices)
        if not flattening.applied or len(flat_vertices) < 3:
            continue
        patch_mask = np.asarray(
            [all(int(value) in flat_vertices for value in face) for face in current_faces],
            dtype=bool,
        )
        patch_face_indices = np.flatnonzero(patch_mask)
        if patch_face_indices.size == 0:
            continue
        loops = _flattened_patch_loops(current_faces, patch_face_indices)
        normal = np.asarray(flattening.gravity, dtype=np.float64)
        replacements = [
            _triangulate_simple_polygon(vertices, loop, normal) for loop in loops
        ]
        replacement_faces = np.concatenate(replacements, axis=0)
        retained_faces = current_faces[~patch_mask]
        current_faces = np.ascontiguousarray(
            np.concatenate((retained_faces, replacement_faces), axis=0),
            dtype=np.int32,
        )
        patches.append(
            {
                "flattening_index": flattening_index,
                "removed_face_count": int(patch_face_indices.size),
                "replacement_face_count": int(replacement_faces.shape[0]),
                "boundary_loop_lengths": [len(loop) for loop in loops],
            }
        )
    if not patches:
        raise ValueError("no flattened support patch could be remeshed")
    compact_vertices, compact_faces = _compact_referenced_vertices(
        vertices,
        current_faces,
    )
    return compact_vertices, compact_faces, {
        "algorithm": "deterministic_planar_patch_ear_clipping_v1",
        "status": "APPLIED",
        "patches": patches,
    }


def _structured_dimensions(
    dimensions: Sequence[int],
) -> tuple[int, int, int]:
    if len(dimensions) != 3:
        raise ValueError("structured grid dimensions must contain three values")
    result = tuple(int(value) for value in dimensions)
    if any(value <= 0 for value in result):
        raise ValueError("structured grid dimensions must be positive")
    if any(isinstance(value, bool) for value in dimensions):
        raise ValueError("structured grid dimensions must be integers")
    return result


def _structured_corner_index(
    dimensions: tuple[int, int, int],
    corner: tuple[int, int, int],
) -> int:
    nx, ny, _nz = (value + 1 for value in dimensions)
    x, y, z = corner
    return int(x + nx * (y + ny * z))


def _structured_cavity_boundary(
    dimensions: Sequence[int],
    empty_cells: set[tuple[int, int, int]],
) -> np.ndarray:
    """Extract inward material-boundary triangles from structured void cells."""

    dims = _structured_dimensions(dimensions)
    for cell in empty_cells:
        if (
            len(cell) != 3
            or any(
                coordinate < 0 or coordinate >= dims[axis]
                for axis, coordinate in enumerate(cell)
            )
        ):
            raise ValueError("structured cavity cell is outside the grid")
    directions = (
        (-1, 0, 0), (1, 0, 0),
        (0, -1, 0), (0, 1, 0),
        (0, 0, -1), (0, 0, 1),
    )
    faces: list[tuple[int, int, int]] = []
    for x, y, z in sorted(empty_cells):
        for dx, dy, dz in directions:
            neighbor = (x + dx, y + dy, z + dz)
            if neighbor in empty_cells:
                continue
            if dx:
                plane = x + (1 if dx > 0 else 0)
                quad = (
                    (plane, y, z),
                    (plane, y + 1, z),
                    (plane, y + 1, z + 1),
                    (plane, y, z + 1),
                )
                positive_normal = dx < 0
            elif dy:
                plane = y + (1 if dy > 0 else 0)
                quad = (
                    (x, plane, z),
                    (x, plane, z + 1),
                    (x + 1, plane, z + 1),
                    (x + 1, plane, z),
                )
                positive_normal = dy < 0
            else:
                plane = z + (1 if dz > 0 else 0)
                quad = (
                    (x, y, plane),
                    (x + 1, y, plane),
                    (x + 1, y + 1, plane),
                    (x, y + 1, plane),
                )
                positive_normal = dz < 0
            indices = tuple(
                _structured_corner_index(dims, corner) for corner in quad
            )
            pair = (
                (indices[0], indices[1], indices[2]),
                (indices[0], indices[2], indices[3]),
            )
            if not positive_normal:
                pair = tuple((a, c, b) for a, b, c in pair)
            faces.extend(pair)
    return np.ascontiguousarray(faces, dtype=np.int32).reshape((-1, 3))


def _structured_occupancy_hash(
    dimensions: Sequence[int],
    empty_cells: set[tuple[int, int, int]],
) -> str:
    """Match the native FNV alpha hash after deterministic cavity repair."""

    dims = _structured_dimensions(dimensions)
    mask = (1 << 64) - 1
    value = 1469598103934665603
    ordinal = 0
    for x in range(dims[0]):
        for y in range(dims[1]):
            for z in range(dims[2]):
                value ^= ordinal
                value = (value * 1099511628211) & mask
                ordinal += 1
                value ^= 0 if (x, y, z) in empty_cells else 1
                value = (value * 1099511628211) & mask
    return str(value)


def _structured_corner_coordinates(
    dimensions: tuple[int, int, int],
    index: int,
) -> tuple[int, int, int]:
    nx = dimensions[0] + 1
    ny = dimensions[1] + 1
    plane = nx * ny
    return index % nx, (index // nx) % ny, index // plane


def _infer_structured_cavity(
    dimensions: Sequence[int],
    faces: np.ndarray,
) -> set[tuple[int, int, int]]:
    """Recover the exact void mask by x-ray parity over the grid topology."""

    dims = _structured_dimensions(dimensions)
    expected_corners = int(np.prod(np.asarray(dims, dtype=np.int64) + 1))
    if faces.size and (int(faces.min()) < 0 or int(faces.max()) >= expected_corners):
        raise ValueError("structured cavity face index is outside the grid")
    crossings: dict[tuple[int, int], dict[int, int]] = {}
    for face in faces:
        corners = [
            _structured_corner_coordinates(dims, int(index)) for index in face
        ]
        constant_axes = [
            axis for axis in range(3)
            if len({corner[axis] for corner in corners}) == 1
        ]
        if len(constant_axes) != 1:
            raise ValueError("structured cavity contains a non-grid-aligned face")
        if constant_axes[0] != 0:
            continue
        key = (
            min(corner[1] for corner in corners),
            min(corner[2] for corner in corners),
        )
        plane = corners[0][0]
        row = crossings.setdefault(key, {})
        row[plane] = row.get(plane, 0) + 1

    empty: set[tuple[int, int, int]] = set()
    for (y, z), counts in crossings.items():
        if any(count != 2 for count in counts.values()):
            raise ValueError("structured cavity x-ray contains an incomplete quad")
        planes = sorted(counts)
        if len(planes) % 2:
            raise ValueError("structured cavity x-ray has odd parity")
        for first, second in zip(planes[0::2], planes[1::2]):
            if first >= second:
                raise ValueError("structured cavity has a degenerate x interval")
            empty.update((x, y, z) for x in range(first, second))
    return empty


def _structured_boundary_issues(
    faces: np.ndarray,
) -> tuple[list[tuple[int, int]], list[int]]:
    edge_faces: dict[tuple[int, int], list[int]] = {}
    vertex_faces: dict[int, list[int]] = {}
    vertex_edges: dict[int, list[tuple[int, int]]] = {}
    for face_index, face in enumerate(faces):
        values = tuple(int(value) for value in face)
        for vertex in values:
            vertex_faces.setdefault(vertex, []).append(face_index)
        for first, second in (
            (values[0], values[1]),
            (values[1], values[2]),
            (values[2], values[0]),
        ):
            edge = tuple(sorted((first, second)))
            uses = edge_faces.setdefault(edge, [])
            if not uses:
                vertex_edges.setdefault(edge[0], []).append(edge)
                vertex_edges.setdefault(edge[1], []).append(edge)
            uses.append(face_index)
    bad_edges = sorted(edge for edge, incident in edge_faces.items() if len(incident) != 2)
    bad_vertices: list[int] = []
    for vertex, incident in vertex_faces.items():
        adjacency = {face: set() for face in incident}
        for edge in vertex_edges.get(vertex, []):
            uses = edge_faces[edge]
            for first in uses:
                adjacency[first].update(
                    second for second in uses if second != first
                )
        remaining = set(incident)
        components = 0
        while remaining:
            components += 1
            pending = [remaining.pop()]
            while pending:
                face = pending.pop()
                discovered = remaining.intersection(adjacency[face])
                remaining.difference_update(discovered)
                pending.extend(discovered)
        if components != 1:
            bad_vertices.append(vertex)
    return bad_edges, sorted(bad_vertices)


def _incident_empty_cells(
    dimensions: tuple[int, int, int],
    feature: tuple[int, ...],
    empty: set[tuple[int, int, int]],
) -> list[tuple[int, int, int]]:
    corners = [_structured_corner_coordinates(dimensions, value) for value in feature]
    candidates: list[tuple[int, int, int]] = []
    for dx in (-1, 0):
        for dy in (-1, 0):
            for dz in (-1, 0):
                cell = (corners[0][0] + dx, corners[0][1] + dy, corners[0][2] + dz)
                if (
                    cell in empty
                    and all(
                        cell[axis] <= corner[axis] <= cell[axis] + 1
                        for corner in corners
                        for axis in range(3)
                    )
                ):
                    candidates.append(cell)
    return sorted(set(candidates))


def _repair_structured_cavity(
    dimensions: Sequence[int],
    empty_cells: set[tuple[int, int, int]],
) -> tuple[set[tuple[int, int, int]], dict[str, Any]]:
    dims = _structured_dimensions(dimensions)
    source = set(empty_cells)
    repaired = set(source)
    refilled: list[tuple[int, int, int]] = []
    topology_transitions: list[dict[str, Any]] = []
    while True:
        faces = _structured_cavity_boundary(dims, repaired)
        bad_edges, bad_vertices = _structured_boundary_issues(faces)
        if not bad_edges and not bad_vertices:
            break
        feature = bad_edges[0] if bad_edges else (bad_vertices[0],)
        candidates = _incident_empty_cells(dims, feature, repaired)
        if not candidates:
            raise ValueError("non-manifold cavity has no refill candidate")
        current_score = len(bad_edges) + len(bad_vertices)
        choices: list[tuple[int, int, tuple[int, int, int]]] = []
        for cell in candidates:
            trial = set(repaired)
            trial.remove(cell)
            trial_edges, trial_vertices = _structured_boundary_issues(
                _structured_cavity_boundary(dims, trial)
            )
            choices.append((len(trial_edges), len(trial_vertices), cell))
        best = min(
            choices,
            key=lambda item: (item[0] + item[1], item[0], item[1], item[2]),
        )
        if best[0] + best[1] > current_score:
            raise ValueError(
                "deterministic cavity refill cannot avoid increasing topology defects"
            )
        repaired.remove(best[2])
        refilled.append(best[2])
        topology_transitions.append(
            {
                "refilled_cell_xyz": list(best[2]),
                "before": {
                    "nonmanifold_edges": len(bad_edges),
                    "nonmanifold_vertices": len(bad_vertices),
                },
                "after": {
                    "nonmanifold_edges": best[0],
                    "nonmanifold_vertices": best[1],
                },
            }
        )
        if len(refilled) > len(source):
            raise RuntimeError("structured cavity repair did not terminate")
    return repaired, {
        "algorithm": "deterministic_nonincreasing_refill_v2",
        "source_cavity_cell_count": len(source),
        "repaired_cavity_cell_count": len(repaired),
        "refilled_cells_xyz": [list(cell) for cell in refilled],
        "topology_transitions": topology_transitions,
        "source_occupancy_hash": _structured_occupancy_hash(dims, source),
        "repaired_occupancy_hash": _structured_occupancy_hash(dims, repaired),
    }


@dataclass(frozen=True)
class SealedBalanceArtifact:
    """Immutable geometry/moment snapshot after beta support projection."""

    _outer_vertices: np.ndarray
    _outer_faces: np.ndarray
    _inner_vertices: np.ndarray
    _inner_faces: np.ndarray
    _support_polygon: np.ndarray
    source_occupancy_hash: str
    occupancy_hash: str
    occupancy_repair_audit: Mapping[str, Any]
    outer_flattening_repair_audit: Mapping[str, Any]
    source_geometry_hash: str
    geometry_sha256: str
    material_boundary_sha256: str
    support_polygon_sha256: str
    gravity: tuple[float, float, float]
    evidence_class: EvidenceClass
    material_mass: float
    center_of_mass: tuple[float, float, float]
    inertia_tensor_about_center: tuple[tuple[float, float, float], ...]
    mass_properties_sha256: str
    safe_region_gamma_radians: float | None
    safe_region_state: str
    safe_region_residual: float | None
    safe_region_inset: float | None
    bounding_box_diagonal: float
    source_outer_boundary_audit: Mapping[str, Any]
    outer_boundary_audit: Mapping[str, Any]
    inner_boundary_audit: Mapping[str, Any]
    wall_thickness_audit: Mapping[str, Any]
    containment_audit: Mapping[str, Any]
    failures: tuple[str, ...]

    @classmethod
    def create(
        cls,
        *,
        outer_vertices: Any,
        outer_faces: Any,
        inner_vertices: Any,
        inner_faces: Any,
        flattenings: Sequence[SupportFlattening],
        occupancy_hash: str,
        source_geometry_hash: str,
        gravity: Sequence[float],
        support_polygon: Any,
        evidence_class: EvidenceClass | str,
        native_module: Any,
        wall_thickness_audit: Mapping[str, Any] | None = None,
        minimum_wall_thickness: float = 0.0,
        structured_grid_dimensions: Sequence[int] | None = None,
        repair_nonmanifold_cavity: bool = False,
        safe_region_gamma_radians: float | None = None,
    ) -> "SealedBalanceArtifact":
        outer = _vertices(outer_vertices)
        faces = _faces(outer_faces)
        inner = _vertices(inner_vertices, "inner_vertices", allow_empty=True)
        raw_inner_faces = np.asarray(inner_faces)
        if raw_inner_faces.size == 0:
            inner_faces_array = np.empty((0, 3), dtype=np.int32)
        else:
            inner_faces_array = _faces(raw_inner_faces)
        if faces.size and int(faces.max()) >= outer.shape[0]:
            raise ValueError("outer face index is outside outer vertices")
        if inner_faces_array.size and int(inner_faces_array.max()) >= inner.shape[0]:
            raise ValueError("inner face index is outside inner vertices")
        source_occupancy_hash = str(occupancy_hash)
        occupancy_repair_audit: dict[str, Any] = {
            "algorithm": "not_requested",
            "source_occupancy_hash": source_occupancy_hash,
            "repaired_occupancy_hash": source_occupancy_hash,
            "refilled_cells_xyz": [],
        }
        if repair_nonmanifold_cavity:
            if structured_grid_dimensions is None:
                raise ValueError(
                    "structured grid dimensions are required for cavity repair"
                )
            dimensions = _structured_dimensions(structured_grid_dimensions)
            expected_corners = int(
                np.prod(np.asarray(dimensions, dtype=np.int64) + 1)
            )
            if inner.shape[0] != expected_corners:
                raise ValueError(
                    "structured cavity repair requires the full grid-corner array"
                )
            source_empty = _infer_structured_cavity(
                dimensions,
                inner_faces_array,
            )
            inferred_hash = _structured_occupancy_hash(dimensions, source_empty)
            if inferred_hash != source_occupancy_hash:
                raise ValueError(
                    "structured cavity does not match the source occupancy hash"
                )
            repaired_empty, occupancy_repair_audit = _repair_structured_cavity(
                dimensions,
                source_empty,
            )
            inner_faces_array = _structured_cavity_boundary(
                dimensions,
                repaired_empty,
            )
            occupancy_hash = occupancy_repair_audit["repaired_occupancy_hash"]
        inner, inner_faces_array = _compact_referenced_vertices(
            inner,
            inner_faces_array,
        )
        source_outer_audit = dict(
            native_module.audit_triangle_boundary(outer, faces)
        )
        flattened = outer.copy()
        for flattening in flattenings:
            flattened = flattening.project(flattened)

        outer_audit = dict(native_module.audit_triangle_boundary(flattened, faces))
        outer_flattening_repair_audit: dict[str, Any] = {
            "algorithm": "deterministic_planar_patch_ear_clipping_v1",
            "status": "NOT_REQUIRED",
            "patches": [],
        }
        repairable_flattening_failure = (
            bool(source_outer_audit.get("strict_boundary_valid"))
            and
            not bool(outer_audit.get("strict_boundary_valid"))
            and outer_audit.get("self_intersection_completed") is True
            and outer_audit.get("self_intersection_free") is False
            and outer_audit.get("closed") is True
            and outer_audit.get("edge_manifold") is True
            and outer_audit.get("vertex_manifold") is True
            and outer_audit.get("consistently_oriented") is True
            and int(outer_audit.get("degenerate_face_count", 0)) == 0
        )
        if repairable_flattening_failure:
            source_audit = outer_audit
            try:
                flattened, faces, outer_flattening_repair_audit = (
                    _remesh_flattened_support_patches(
                        flattened,
                        faces,
                        flattenings,
                    )
                )
                outer_audit = dict(
                    native_module.audit_triangle_boundary(flattened, faces)
                )
                outer_flattening_repair_audit["source_audit"] = source_audit
                outer_flattening_repair_audit["result_audit"] = outer_audit
                if not bool(outer_audit.get("strict_boundary_valid")):
                    outer_flattening_repair_audit["status"] = "FAILED"
            except (RuntimeError, ValueError) as error:
                outer_flattening_repair_audit = {
                    "algorithm": "deterministic_planar_patch_ear_clipping_v1",
                    "status": "FAILED",
                    "reason": str(error),
                    "patches": [],
                    "source_audit": source_audit,
                }

        gravity_array = np.asarray(gravity, dtype=np.float64)
        if (
            gravity_array.shape != (3,)
            or not np.isfinite(gravity_array).all()
            or float(np.linalg.norm(gravity_array)) <= 1.0e-15
        ):
            raise ValueError("gravity must be a finite nonzero 3-vector")
        gravity_array /= np.linalg.norm(gravity_array)
        polygon = _vertices(support_polygon, "support_polygon")
        for flattening in flattenings:
            flattening_gravity = np.asarray(flattening.gravity, dtype=np.float64)
            flattening_gravity /= np.linalg.norm(flattening_gravity)
            if np.allclose(
                flattening_gravity,
                gravity_array,
                rtol=0.0,
                atol=64.0 * np.finfo(np.float64).eps,
            ):
                heights = -(polygon @ flattening_gravity)
                polygon += (
                    (flattening.support_height - heights)[:, None]
                    * (-flattening_gravity)[None, :]
                )

        if inner_faces_array.shape[0]:
            # Cavity faces are inward for the material boundary. Reverse a
            # copy solely for the void's outward-oriented topology audit.
            outward_void_faces = inner_faces_array[:, [0, 2, 1]]
            inner_audit = dict(
                native_module.audit_triangle_boundary(inner, outward_void_faces)
            )
        else:
            inner_audit = {
                "strict_boundary_valid": True,
                "empty_cavity": True,
                "face_count": 0,
            }

        outer_mass, outer_first, outer_second = _surface_moments(flattened, faces)
        inner_mass, inner_first, inner_second = _surface_moments(
            inner,
            inner_faces_array,
        )
        mass = outer_mass + inner_mass
        first_moment = outer_first + inner_first
        second_moment = outer_second + inner_second
        failures: list[str] = []
        if not bool(source_outer_audit.get("strict_boundary_valid")):
            failures.append("source_outer_boundary_invalid")
        if not bool(outer_audit.get("strict_boundary_valid")):
            failures.append("outer_boundary_invalid")
        if not bool(inner_audit.get("strict_boundary_valid")):
            failures.append("inner_boundary_invalid")
        if not math.isfinite(mass) or mass <= 0.0:
            failures.append("nonpositive_material_mass")
            center = np.zeros(3)
            inertia = np.zeros((3, 3))
        else:
            center = first_moment / mass
            inertia_origin = np.trace(second_moment) * np.eye(3) - second_moment
            inertia = inertia_origin - mass * (
                float(center @ center) * np.eye(3) - np.outer(center, center)
            )
            inertia = 0.5 * (inertia + inertia.T)
            if (
                not np.isfinite(center).all()
                or not np.isfinite(inertia).all()
                or float(np.linalg.eigvalsh(inertia).min()) <= 0.0
            ):
                failures.append("invalid_material_moments")

        bounding_box_diagonal = float(np.linalg.norm(np.ptp(flattened, axis=0)))
        safe_state = "NOT_REQUESTED"
        safe_residual: float | None = None
        safe_inset: float | None = None
        if safe_region_gamma_radians is not None:
            gamma = float(safe_region_gamma_radians)
            if not math.isfinite(gamma) or gamma < 0.0 or gamma >= 0.5 * math.pi:
                raise ValueError(
                    "safe-region gamma must be finite in [0, pi/2)"
                )
            if failures or not (bounding_box_diagonal > 0.0):
                safe_state = "INVALID_GEOMETRY"
            else:
                from .stability_oracle import evaluate_quasi_static_tilt

                nominal = evaluate_quasi_static_tilt(
                    polygon,
                    center,
                    gravity_array,
                    tilt_degrees=0.0,
                    azimuth_degrees=0.0,
                )
                safe_inset = nominal.center_height * math.tan(gamma)
                safe_residual = max(
                    0.0,
                    safe_inset - nominal.nominal_signed_margin,
                )
                tolerance = bounding_box_diagonal * 1.0e-6
                safe_state = "SAFE" if safe_residual <= tolerance else "OUTSIDE"

        requested_wall_thickness = float(minimum_wall_thickness)
        if (
            not math.isfinite(requested_wall_thickness)
            or requested_wall_thickness < 0.0
        ):
            raise ValueError(
                "minimum wall thickness must be finite and non-negative"
            )
        caller_wall_context = dict(wall_thickness_audit or {})
        if not inner_faces_array.shape[0]:
            wall_audit = {
                "status": "PASS",
                "reason": "no_cavity",
                "requested_thickness": requested_wall_thickness,
                "caller_context": caller_wall_context,
            }
            containment_audit = {
                "status": "PASS",
                "reason": "no_cavity",
            }
        else:
            material_shell_audit = dict(
                native_module.audit_material_shells(
                    flattened,
                    faces,
                    inner,
                    inner_faces_array,
                    requested_wall_thickness,
                )
            )
            wall_audit = dict(material_shell_audit["wall_clearance"])
            wall_pass = bool(
                wall_audit.get("completed")
                and wall_audit.get("intersection_free")
                and wall_audit.get("meets_requested_thickness")
            )
            wall_audit["status"] = "PASS" if wall_pass else "FAIL"
            wall_audit["caller_context"] = caller_wall_context
            containment_audit = dict(material_shell_audit["containment"])
            containment_pass = bool(
                containment_audit.get("completed")
                and containment_audit.get("all_cavities_inside_outer")
            )
            containment_audit["status"] = (
                "PASS" if containment_pass else "FAIL"
            )
            if not wall_pass:
                failures.append("wall_clearance_invalid")
            if not containment_pass:
                failures.append("cavity_containment_invalid")

        geometry_hash = _array_hash(flattened, faces)
        material_hash = _array_hash(
            flattened,
            faces,
            inner,
            inner_faces_array,
        )
        outer_text = _obj_text(flattened, faces)
        collision_hash = hashlib.sha256(outer_text.encode("utf-8")).hexdigest()
        entry = {
            "occupancy_hash": str(occupancy_hash),
            "geometry_hash": geometry_hash,
            "collision_mesh_sha256": collision_hash,
            "material_mass_native": mass,
            "center_of_mass_native": [float(item) for item in center],
            "inertia_tensor_about_center_native": [
                [float(item) for item in row] for row in inertia
            ],
            "gravity_direction_native": [float(item) for item in gravity_array],
        }
        from .pybullet_runner import _mass_properties_fingerprint

        return cls(
            _outer_vertices=np.ascontiguousarray(flattened),
            _outer_faces=np.ascontiguousarray(faces),
            _inner_vertices=np.ascontiguousarray(inner),
            _inner_faces=np.ascontiguousarray(inner_faces_array),
            _support_polygon=np.ascontiguousarray(polygon),
            source_occupancy_hash=source_occupancy_hash,
            occupancy_hash=str(occupancy_hash),
            occupancy_repair_audit=occupancy_repair_audit,
            outer_flattening_repair_audit=outer_flattening_repair_audit,
            source_geometry_hash=str(source_geometry_hash),
            geometry_sha256=geometry_hash,
            material_boundary_sha256=material_hash,
            support_polygon_sha256=_canonical_hash(polygon.tolist()),
            gravity=tuple(float(item) for item in gravity_array),
            evidence_class=EvidenceClass(evidence_class),
            material_mass=float(mass),
            center_of_mass=tuple(float(item) for item in center),
            inertia_tensor_about_center=tuple(
                tuple(float(item) for item in row) for row in inertia
            ),
            mass_properties_sha256=_mass_properties_fingerprint(entry),
            safe_region_gamma_radians=(
                None
                if safe_region_gamma_radians is None
                else float(safe_region_gamma_radians)
            ),
            safe_region_state=safe_state,
            safe_region_residual=safe_residual,
            safe_region_inset=safe_inset,
            bounding_box_diagonal=bounding_box_diagonal,
            source_outer_boundary_audit=source_outer_audit,
            outer_boundary_audit=outer_audit,
            inner_boundary_audit=inner_audit,
            wall_thickness_audit=wall_audit,
            containment_audit=containment_audit,
            failures=tuple(failures),
        )

    @property
    def ready(self) -> bool:
        return not self.failures

    @property
    def manufacturing_ready(self) -> bool:
        return self.ready and self.wall_thickness_audit.get("status") == "PASS"

    @property
    def outer_vertices(self) -> np.ndarray:
        return self._outer_vertices.copy()

    @property
    def outer_faces(self) -> np.ndarray:
        return self._outer_faces.copy()

    @property
    def inner_vertices(self) -> np.ndarray:
        return self._inner_vertices.copy()

    @property
    def inner_faces(self) -> np.ndarray:
        return self._inner_faces.copy()

    @property
    def support_polygon(self) -> np.ndarray:
        return self._support_polygon.copy()

    def pybullet_entry(
        self,
        *,
        name: str = "sealed-case",
        collision_mesh: str = "outer-flattened.obj",
    ) -> dict[str, Any]:
        if not self.ready:
            raise RuntimeError(
                "artifact is not sealed: " + ", ".join(self.failures)
            )
        outer_text = _obj_text(self._outer_vertices, self._outer_faces)
        entry: dict[str, Any] = {
            "name": name,
            "collision_mesh": collision_mesh,
            "visual_mesh": collision_mesh,
            "collision_mesh_sha256": hashlib.sha256(
                outer_text.encode("utf-8")
            ).hexdigest(),
            "occupancy_hash": self.occupancy_hash,
            "source_occupancy_hash": self.source_occupancy_hash,
            "occupancy_repair": dict(self.occupancy_repair_audit),
            "outer_flattening_repair": dict(
                self.outer_flattening_repair_audit
            ),
            "geometry_hash": self.geometry_sha256,
            "material_mass_native": self.material_mass,
            "center_of_mass_native": list(self.center_of_mass),
            "inertia_tensor_about_center_native": [
                list(row) for row in self.inertia_tensor_about_center
            ],
            "gravity_direction_native": list(self.gravity),
            "expected_stable": True,
            "provenance": "sealed-post-flatten-material-boundary-v1",
            "evidence_class": self.evidence_class.value,
            "sealed_artifact_status": "SEALED",
            "manufacturing_ready": self.manufacturing_ready,
            "wall_thickness_audit": dict(self.wall_thickness_audit),
            "containment_audit": dict(self.containment_audit),
            "safe_region_state": self.safe_region_state,
            "safe_region_residual": self.safe_region_residual,
            "safe_region_inset": self.safe_region_inset,
            "bounding_box_diagonal": self.bounding_box_diagonal,
        }
        if self._inner_faces.shape[0]:
            inner_text = _obj_text(self._inner_vertices, self._inner_faces)
            entry["inner_void_mesh"] = "inner-cavity.obj"
            entry["inner_void_mesh_sha256"] = hashlib.sha256(
                inner_text.encode("utf-8")
            ).hexdigest()
        entry["mass_properties_sha256"] = self.mass_properties_sha256
        return entry

    def manifest(self) -> dict[str, Any]:
        entry = self.pybullet_entry() if self.ready else None
        return {
            "schema": "make-it-stand-sealed-balance-artifact-v1",
            "schema_version": 1,
            "status": "SEALED" if self.ready else "REJECTED",
            "simulation_ready": self.ready,
            "manufacturing_ready": self.manufacturing_ready,
            "failures": list(self.failures),
            "evidence_class": self.evidence_class.value,
            "occupancy_hash": self.occupancy_hash,
            "source_occupancy_hash": self.source_occupancy_hash,
            "occupancy_repair": dict(self.occupancy_repair_audit),
            "outer_flattening_repair": dict(
                self.outer_flattening_repair_audit
            ),
            "source_geometry_hash": self.source_geometry_hash,
            "geometry_sha256": self.geometry_sha256,
            "material_boundary_sha256": self.material_boundary_sha256,
            "support_polygon_sha256": self.support_polygon_sha256,
            "support_polygon": self._support_polygon.tolist(),
            "gravity_direction_native": list(self.gravity),
            "mass_properties_sha256": self.mass_properties_sha256,
            "safe_region": {
                "gamma_radians": self.safe_region_gamma_radians,
                "state": self.safe_region_state,
                "residual": self.safe_region_residual,
                "inset": self.safe_region_inset,
                "bounding_box_diagonal": self.bounding_box_diagonal,
                "residual_limit": self.bounding_box_diagonal * 1.0e-6,
            },
            "source_outer_boundary_audit": dict(
                self.source_outer_boundary_audit
            ),
            "outer_boundary_audit": dict(self.outer_boundary_audit),
            "inner_boundary_audit": dict(self.inner_boundary_audit),
            "wall_thickness_audit": dict(self.wall_thickness_audit),
            "containment_audit": dict(self.containment_audit),
            "pybullet_entry": entry,
        }

    def write_directory(self, directory: os.PathLike[str] | str) -> Path:
        root = Path(directory)
        root.mkdir(parents=True, exist_ok=True)
        write_text_atomic(
            root / "outer-flattened.obj",
            _obj_text(self._outer_vertices, self._outer_faces),
        )
        if self._inner_faces.shape[0]:
            write_text_atomic(
                root / "inner-cavity.obj",
                _obj_text(self._inner_vertices, self._inner_faces),
            )
        manifest_path = root / "sealed-artifact.json"
        write_json_atomic(manifest_path, self.manifest())
        return manifest_path


def load_sealed_pybullet_entry(
    manifest_path: os.PathLike[str] | str,
) -> dict[str, Any]:
    """Load one sealed entry and revalidate every referenced byte hash."""

    path = Path(manifest_path).expanduser().resolve()
    payload = json.loads(path.read_text(encoding="utf-8"))
    if (
        not isinstance(payload, Mapping)
        or payload.get("schema") != "make-it-stand-sealed-balance-artifact-v1"
        or payload.get("schema_version") != 1
    ):
        raise ValueError("unsupported sealed artifact manifest")
    if payload.get("status") != "SEALED" or payload.get("simulation_ready") is not True:
        raise ValueError("sealed artifact is not simulation-ready")
    if payload.get("manufacturing_ready") is not True:
        raise ValueError("sealed artifact has not passed material-shell audits")
    if payload.get("wall_thickness_audit", {}).get("status") != "PASS":
        raise ValueError("sealed artifact wall-thickness audit did not pass")
    if payload.get("containment_audit", {}).get("status") != "PASS":
        raise ValueError("sealed artifact containment audit did not pass")
    entry = payload.get("pybullet_entry")
    if not isinstance(entry, Mapping):
        raise ValueError("sealed artifact has no PyBullet entry")
    result = dict(entry)
    root = path.parent
    for field, hash_field in (
        ("collision_mesh", "collision_mesh_sha256"),
        ("visual_mesh", "collision_mesh_sha256"),
        ("inner_void_mesh", "inner_void_mesh_sha256"),
    ):
        raw = result.get(field)
        if raw is None:
            continue
        candidate = Path(str(raw))
        if candidate.is_absolute() or ".." in candidate.parts:
            raise ValueError(f"{field} must remain inside the sealed directory")
        resolved = (root / candidate).resolve()
        try:
            resolved.relative_to(root)
        except ValueError as error:
            raise ValueError(f"{field} escapes the sealed directory") from error
        if not resolved.is_file():
            raise ValueError(f"sealed {field} is missing")
        expected = result.get(hash_field)
        actual = hashlib.sha256(resolved.read_bytes()).hexdigest()
        if actual != expected:
            label = "collision mesh" if field != "inner_void_mesh" else "inner mesh"
            raise ValueError(f"sealed {label} hash mismatch")
        result[field] = str(resolved)
    if result.get("geometry_hash") != payload.get("geometry_sha256"):
        raise ValueError("sealed geometry hash mismatch")
    if result.get("mass_properties_sha256") != payload.get(
        "mass_properties_sha256"
    ):
        raise ValueError("sealed mass-properties hash mismatch")
    from .pybullet_runner import _mass_properties_fingerprint
    if result["mass_properties_sha256"] != _mass_properties_fingerprint(result):
        raise ValueError("sealed mass-properties payload mismatch")
    return result


__all__ = ["SealedBalanceArtifact", "load_sealed_pybullet_entry"]
