"""Independent signed-tetrahedron mass and moment integration.

Only Python's standard library is used.  In particular, this module does not
call the native solver's mass-moment integrator.  Validation covers finite
coordinates, triangle validity, and oriented edge and vertex manifoldness.
It intentionally does not certify self-intersection freedom, nesting or
containment of negative components, wall thickness, or overlap between
density regions.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterable, TypeAlias


Vector3: TypeAlias = tuple[float, float, float]
Matrix3: TypeAlias = tuple[Vector3, Vector3, Vector3]
Face: TypeAlias = tuple[int, int, int]


class GeometryValidationError(ValueError):
    """Raised when an input cannot be accepted as oriented material geometry."""


@dataclass(frozen=True)
class TriangleMesh:
    """Immutable triangle surface input; semantic checks occur at integration."""

    vertices: tuple[Vector3, ...]
    faces: tuple[Face, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "vertices", tuple(tuple(vertex) for vertex in self.vertices))
        object.__setattr__(self, "faces", tuple(tuple(face) for face in self.faces))


@dataclass(frozen=True)
class DensityRegion:
    """An oriented material region whose constant density is applied once."""

    mesh: TriangleMesh
    density: float


@dataclass(frozen=True)
class GeometricMoments:
    """Density-free volume, first moment, and symmetric raw second moment."""

    volume: float
    first_moment: Vector3
    second_moment: Matrix3


@dataclass(frozen=True)
class MassProperties:
    """Density-weighted raw moments and inertia about the center of mass."""

    mass: float
    first_moment: Vector3
    second_moment: Matrix3
    center_of_mass: Vector3
    central_inertia: Matrix3


@dataclass(frozen=True)
class _PreparedMesh:
    vertices: tuple[Vector3, ...]
    faces: tuple[Face, ...]
    reference: Vector3
    scale: float


def _cross(left: Vector3, right: Vector3) -> Vector3:
    return (
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    )


def _dot(left: Vector3, right: Vector3) -> float:
    return left[0] * right[0] + left[1] * right[1] + left[2] * right[2]


def _prepare(mesh: TriangleMesh) -> _PreparedMesh:
    if not isinstance(mesh, TriangleMesh):
        raise TypeError("mesh must be a TriangleMesh")
    if not mesh.vertices:
        raise GeometryValidationError("mesh has no vertices")
    if not mesh.faces:
        raise GeometryValidationError("mesh has no faces")

    vertices: list[Vector3] = []
    for vertex_index, vertex in enumerate(mesh.vertices):
        if len(vertex) != 3:
            raise GeometryValidationError(f"vertex {vertex_index} does not have three coordinates")
        try:
            converted = (float(vertex[0]), float(vertex[1]), float(vertex[2]))
        except (TypeError, ValueError, OverflowError) as error:
            raise GeometryValidationError(f"vertex {vertex_index} is not numeric") from error
        if not all(math.isfinite(value) for value in converted):
            raise GeometryValidationError(f"vertex {vertex_index} is not finite")
        vertices.append(converted)

    minima = tuple(min(vertex[axis] for vertex in vertices) for axis in range(3))
    maxima = tuple(max(vertex[axis] for vertex in vertices) for axis in range(3))
    spans = tuple(maxima[axis] - minima[axis] for axis in range(3))
    scale = max(spans)
    if not math.isfinite(scale) or scale <= 0.0:
        raise GeometryValidationError("mesh has zero or unrepresentable spatial extent")
    reference = tuple(minima[axis] + spans[axis] * 0.5 for axis in range(3))
    normalized = tuple(
        tuple((vertex[axis] - reference[axis]) / scale for axis in range(3))
        for vertex in vertices
    )

    faces: list[Face] = []
    seen_faces: set[tuple[int, int, int]] = set()
    edge_uses: dict[tuple[int, int], list[tuple[int, int, int]]] = {}
    incident_faces: list[set[int]] = [set() for _ in vertices]
    for face_index, face in enumerate(mesh.faces):
        if len(face) != 3:
            raise GeometryValidationError(f"face {face_index} does not have three indices")
        if any(not isinstance(index, int) or isinstance(index, bool) for index in face):
            raise GeometryValidationError(f"face {face_index} has a non-integer index")
        converted_face = (face[0], face[1], face[2])
        if any(index < 0 or index >= len(vertices) for index in converted_face):
            raise GeometryValidationError(f"face {face_index} has an out-of-range index")
        if len(set(converted_face)) != 3:
            raise GeometryValidationError(f"face {face_index} repeats a vertex")
        face_key = tuple(sorted(converted_face))
        if face_key in seen_faces:
            raise GeometryValidationError(f"face {face_index} duplicates another triangle")
        seen_faces.add(face_key)

        a, b, c = (normalized[index] for index in converted_face)
        ab = tuple(b[axis] - a[axis] for axis in range(3))
        ac = tuple(c[axis] - a[axis] for axis in range(3))
        normal = _cross(ab, ac)
        if _dot(normal, normal) == 0.0:
            raise GeometryValidationError(f"face {face_index} is geometrically degenerate")

        faces.append(converted_face)
        for vertex_index in converted_face:
            incident_faces[vertex_index].add(face_index)
        for start, end in (
            (converted_face[0], converted_face[1]),
            (converted_face[1], converted_face[2]),
            (converted_face[2], converted_face[0]),
        ):
            edge_uses.setdefault(tuple(sorted((start, end))), []).append(
                (start, end, face_index)
            )

    for edge, uses in edge_uses.items():
        if len(uses) != 2:
            raise GeometryValidationError(
                f"edge {edge} is not manifold: expected two incident faces, got {len(uses)}"
            )
        first, second = uses
        if first[:2] != (second[1], second[0]):
            raise GeometryValidationError(f"edge {edge} has inconsistent face orientation")

    # An edge-manifold surface can still have two surface fans pinched at one
    # vertex.  Require each vertex link to be a single connected cycle.
    vertex_neighbors: list[dict[int, set[int]]] = [dict() for _ in vertices]
    for edge, uses in edge_uses.items():
        first_face, second_face = uses[0][2], uses[1][2]
        for vertex_index in edge:
            graph = vertex_neighbors[vertex_index]
            graph.setdefault(first_face, set()).add(second_face)
            graph.setdefault(second_face, set()).add(first_face)
    for vertex_index, incident in enumerate(incident_faces):
        if not incident:
            raise GeometryValidationError(f"vertex {vertex_index} is unreferenced")
        pending = [next(iter(incident))]
        reached: set[int] = set()
        while pending:
            face_index = pending.pop()
            if face_index in reached:
                continue
            reached.add(face_index)
            pending.extend(vertex_neighbors[vertex_index].get(face_index, ()))
        if reached != incident:
            raise GeometryValidationError(
                f"vertex {vertex_index} is non-manifold: incident faces form multiple fans"
            )

    return _PreparedMesh(normalized, tuple(faces), reference, scale)


def _normalized_moments(prepared: _PreparedMesh) -> tuple[float, Vector3, Matrix3]:
    volumes: list[float] = []
    first_terms = [[], [], []]
    second_terms = [[[] for _ in range(3)] for _ in range(3)]
    for face in prepared.faces:
        a, b, c = (prepared.vertices[index] for index in face)
        volume = _dot(a, _cross(b, c)) / 6.0
        volumes.append(volume)
        total = tuple(a[axis] + b[axis] + c[axis] for axis in range(3))
        for row in range(3):
            first_terms[row].append(volume * total[row] / 4.0)
            for column in range(3):
                diagonal_sum = (
                    a[row] * a[column]
                    + b[row] * b[column]
                    + c[row] * c[column]
                )
                second_terms[row][column].append(
                    volume * (diagonal_sum + total[row] * total[column]) / 20.0
                )
    return (
        math.fsum(volumes),
        tuple(math.fsum(axis_terms) for axis_terms in first_terms),
        tuple(
            tuple(math.fsum(second_terms[row][column]) for column in range(3))
            for row in range(3)
        ),
    )


def _scaled_local_moments(prepared: _PreparedMesh) -> tuple[float, Vector3, Matrix3]:
    normalized_volume, normalized_first, normalized_second = _normalized_moments(prepared)
    if normalized_volume <= 0.0:
        raise GeometryValidationError("mesh has zero or negative total oriented volume")
    try:
        scale3 = prepared.scale**3
        scale4 = prepared.scale**4
        scale5 = prepared.scale**5
    except OverflowError as error:
        raise GeometryValidationError("moments are not representable as finite floats") from error
    volume = normalized_volume * scale3
    first = tuple(value * scale4 for value in normalized_first)
    second = tuple(tuple(value * scale5 for value in row) for row in normalized_second)
    values = (volume, *first, *(value for row in second for value in row))
    if not all(math.isfinite(value) for value in values) or volume <= 0.0:
        raise GeometryValidationError("moments are not representable as finite positive-volume floats")
    return volume, first, second


def integrate_geometric_moments(mesh: TriangleMesh) -> GeometricMoments:
    """Integrate density-free signed moments of an oriented closed surface.

    A positive outer component and consistently reversed inner components
    produce subtraction by orientation.  Positive net volume is required,
    but component containment is outside this checker's scope.
    """

    prepared = _prepare(mesh)
    local_volume, local_first, local_second = _scaled_local_moments(prepared)
    reference = prepared.reference
    first = tuple(
        local_first[axis] + local_volume * reference[axis] for axis in range(3)
    )
    second = tuple(
        tuple(
            local_second[row][column]
            + reference[row] * local_first[column]
            + local_first[row] * reference[column]
            + local_volume * reference[row] * reference[column]
            for column in range(3)
        )
        for row in range(3)
    )
    values = (local_volume, *first, *(value for row in second for value in row))
    if not all(math.isfinite(value) for value in values):
        raise GeometryValidationError("moments are not representable as finite floats")
    return GeometricMoments(local_volume, first, second)


def _mass_properties_from_anchor(
    mass: float,
    first_at_anchor: Vector3,
    second_at_anchor: Matrix3,
    anchor: Vector3,
) -> MassProperties:
    if not math.isfinite(mass) or mass <= 0.0:
        raise GeometryValidationError("total material mass must be finite and positive")
    center_at_anchor = tuple(value / mass for value in first_at_anchor)
    center = tuple(anchor[axis] + center_at_anchor[axis] for axis in range(3))
    centered_second = tuple(
        tuple(
            second_at_anchor[row][column]
            - mass * center_at_anchor[row] * center_at_anchor[column]
            for column in range(3)
        )
        for row in range(3)
    )
    trace = sum(centered_second[axis][axis] for axis in range(3))
    inertia = tuple(
        tuple(
            (trace if row == column else 0.0) - centered_second[row][column]
            for column in range(3)
        )
        for row in range(3)
    )
    first = tuple(
        first_at_anchor[axis] + mass * anchor[axis] for axis in range(3)
    )
    second = tuple(
        tuple(
            second_at_anchor[row][column]
            + anchor[row] * first_at_anchor[column]
            + first_at_anchor[row] * anchor[column]
            + mass * anchor[row] * anchor[column]
            for column in range(3)
        )
        for row in range(3)
    )
    values = (
        *first,
        *(value for row in second for value in row),
        *center,
        *(value for row in inertia for value in row),
    )
    if not all(math.isfinite(value) for value in values):
        raise GeometryValidationError("mass properties are not representable as finite floats")
    return MassProperties(mass, first, second, center, inertia)


def aggregate_mass_properties(regions: Iterable[DensityRegion]) -> MassProperties:
    """Aggregate constant-density regions without certifying their overlap."""

    region_list = tuple(regions)
    if not region_list:
        raise GeometryValidationError("at least one density region is required")
    prepared_regions: list[tuple[float, _PreparedMesh, float, Vector3, Matrix3]] = []
    for region_index, region in enumerate(region_list):
        if not isinstance(region, DensityRegion):
            raise TypeError(f"region {region_index} must be a DensityRegion")
        try:
            density = float(region.density)
        except (TypeError, ValueError, OverflowError) as error:
            raise GeometryValidationError(f"region {region_index} density is not numeric") from error
        if not math.isfinite(density) or density <= 0.0:
            raise GeometryValidationError(f"region {region_index} density must be finite and positive")
        prepared = _prepare(region.mesh)
        volume, local_first, local_second = _scaled_local_moments(prepared)
        prepared_regions.append((density, prepared, volume, local_first, local_second))

    anchor = prepared_regions[0][1].reference
    masses: list[float] = []
    first_terms = [[], [], []]
    second_terms = [[[] for _ in range(3)] for _ in range(3)]
    for density, prepared, volume, local_first, local_second in prepared_regions:
        offset = tuple(prepared.reference[axis] - anchor[axis] for axis in range(3))
        first_at_anchor = tuple(
            local_first[axis] + volume * offset[axis] for axis in range(3)
        )
        second_at_anchor = tuple(
            tuple(
                local_second[row][column]
                + offset[row] * local_first[column]
                + local_first[row] * offset[column]
                + volume * offset[row] * offset[column]
                for column in range(3)
            )
            for row in range(3)
        )
        masses.append(density * volume)
        for row in range(3):
            first_terms[row].append(density * first_at_anchor[row])
            for column in range(3):
                second_terms[row][column].append(density * second_at_anchor[row][column])
    return _mass_properties_from_anchor(
        math.fsum(masses),
        tuple(math.fsum(terms) for terms in first_terms),
        tuple(
            tuple(math.fsum(second_terms[row][column]) for column in range(3))
            for row in range(3)
        ),
        anchor,
    )


def compute_mass_properties(mesh: TriangleMesh, density: float = 1.0) -> MassProperties:
    """Compute constant-density properties for one oriented surface mesh."""

    return aggregate_mass_properties((DensityRegion(mesh, density),))
