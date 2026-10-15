"""Independent signed-tetrahedron mass and moment integration.

The default uses Python's standard library; an explicitly selected, pinned
GMP backend accelerates exact rational integration. Neither calls the native
solver's mass-moment integrator. Validation covers finite
coordinates, triangle validity, and oriented edge and vertex manifoldness.
It intentionally does not certify self-intersection freedom, nesting or
containment of negative components, wall thickness, or overlap between
density regions.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from .rational_text import fraction_from_rational_text, fraction_to_rational_text
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


@dataclass(frozen=True)
class _PreparedRegion:
    mesh: _PreparedMesh
    mass: float
    center_local: Vector3
    centered_second: Matrix3


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
        # Internal affine transforms may provide exact rational coordinates.
        # Do not round these before either predicates or integration.
        vertices.append(tuple(value if isinstance(value, Fraction) else number
                              for value, number in zip(vertex, converted)))

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

    def exact_normalized():
        return tuple(tuple((Fraction(value) - Fraction(reference[axis])) / Fraction(scale)
                           for axis, value in enumerate(vertex)) for vertex in vertices)

    exact = any(isinstance(value, Fraction) for vertex in vertices for value in vertex)
    if exact:
        normalized = exact_normalized()

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
            # Binary64 centering can merge distinct source vertices. Escalate
            # the WHOLE integral, not just the predicate: otherwise a valid
            # sliver would be accepted but its material silently discarded.
            if not exact:
                normalized = exact_normalized()
                exact = True
            a, b, c = (normalized[index] for index in converted_face)
            normal = _cross(tuple(b[i]-a[i] for i in range(3)),
                            tuple(c[i]-a[i] for i in range(3)))
            if not any(normal):
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

    if exact:
        # Recheck earlier faces as well after escalation to exact arithmetic.
        for face_index, face in enumerate(faces):
            a, b, c = (normalized[index] for index in face)
            if not any(_cross(tuple(b[i]-a[i] for i in range(3)),
                              tuple(c[i]-a[i] for i in range(3)))):
                raise GeometryValidationError(f"face {face_index} is geometrically degenerate")
    return _PreparedMesh(normalized, tuple(faces), tuple(map(float, reference)), float(scale))


def _normalized_moments(prepared: _PreparedMesh) -> tuple[float, Vector3, Matrix3]:
    exact = isinstance(prepared.vertices[0][0], Fraction)
    accumulate = (lambda values: float(sum(values))) if exact else math.fsum
    volumes: list[float] = []
    first_terms = [[], [], []]
    second_terms = [[[] for _ in range(3)] for _ in range(3)]
    for face in prepared.faces:
        a, b, c = (prepared.vertices[index] for index in face)
        volume = _dot(a, _cross(b, c)) / 6
        volumes.append(volume)
        total = tuple(a[axis] + b[axis] + c[axis] for axis in range(3))
        for row in range(3):
            first_terms[row].append(volume * total[row] / 4)
            for column in range(3):
                diagonal_sum = (
                    a[row] * a[column]
                    + b[row] * b[column]
                    + c[row] * c[column]
                )
                second_terms[row][column].append(
                    volume * (diagonal_sum + total[row] * total[column]) / 20
                )
    return (
        accumulate(volumes),
        tuple(accumulate(axis_terms) for axis_terms in first_terms),
        tuple(
            tuple(accumulate(second_terms[row][column]) for column in range(3))
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


def integrate_exact_rational_raw_moments(vertices, faces, *,
        arithmetic_backend='FRACTION_V1') -> tuple[Fraction, ...]:
    """Integrate ten density-free moments without rounding rational vertices.

    This signed-tetrahedron checker is independent of the native integrator.
    It does *not* validate closure, nesting, self-intersection, or whether the
    represented solid equals a Boolean input operation; callers must check
    those separately. Coordinates may be finite binary64 numbers or exact
    ``numerator/denominator`` strings (including CGAL binary integers).
    The order is M,LX,LY,LZ,QXX,QYY,QZZ,QYZ,QXZ,QXY.
    The explicit GMP backend keeps the same exact formula and public Fraction
    return contract. It never falls back if its pinned dependency is missing.
    """
    if arithmetic_backend == 'FRACTION_V1':
        convert=lambda value:value
        zero=Fraction(0)
    elif arithmetic_backend == 'GMPY2_2_2_1_EXACT_V1':
        import gmpy2
        if gmpy2.version()!='2.2.1' or gmpy2.mp_version()!='GMP 6.3.0':
            raise ValueError('rational integration GMP dependency version mismatch')
        convert=lambda value:gmpy2.mpq(value.numerator,value.denominator)
        zero=gmpy2.mpq(0)
    else:
        raise ValueError('unsupported rational integration backend')
    def integer(token: str) -> int:
        sign = -1 if token.startswith('-') else 1
        digits = token[1:] if sign < 0 else token
        if digits.startswith('0b'):
            return sign * int(digits[2:], 2)
        return sign * int(digits, 10)

    def rational(value) -> Fraction:
        if isinstance(value, bool):
            raise GeometryValidationError('Boolean is not a rational coordinate')
        if isinstance(value, str) and '/' in value:
            numerator, denominator = value.split('/')
            return (fraction_from_rational_text(value) if not numerator.lstrip('-').startswith('0b')
                    else Fraction(integer(numerator), integer(denominator)))
        if isinstance(value, float) and not math.isfinite(value):
            raise GeometryValidationError('nonfinite rational coordinate')
        try:
            return Fraction(value)
        except (TypeError, ValueError, ZeroDivisionError) as error:
            raise GeometryValidationError('invalid rational coordinate') from error

    points = []
    for point in vertices:
        if len(point) != 3:
            raise GeometryValidationError('rational vertex needs three coordinates')
        points.append(tuple(convert(rational(value)) for value in point))
    if not points:
        raise GeometryValidationError('rational material has no vertices')
    # One ten-moment block per occupied power-of-two face count. Streaming
    # balanced reduction avoids repeatedly growing a global denominator, while
    # keeping the signed-tetrahedron formula and exact rational result intact.
    # This checker does not call the solver's surface integral or accumulator.
    partials = []
    face_count = 0
    for face in faces:
        if len(face) != 3 or any(isinstance(index, bool) or not isinstance(index, int)
                                 or index < 0 or index >= len(points) for index in face) or \
                len(set(face)) != 3:
            raise GeometryValidationError('invalid rational material face')
        a, b, c = (points[index] for index in face)
        cross = (b[1]*c[2]-b[2]*c[1],
                 b[2]*c[0]-b[0]*c[2],
                 b[0]*c[1]-b[1]*c[0])
        volume = sum(a[axis]*cross[axis] for axis in range(3)) / 6
        total = tuple(a[axis]+b[axis]+c[axis] for axis in range(3))
        contribution = [volume]
        for axis in range(3):
            contribution.append(volume*total[axis]/4)
        for row, column in ((0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)):
            diagonal = sum(point[row]*point[column] for point in (a, b, c))
            contribution.append(volume*(diagonal+total[row]*total[column])/20)
        level = 0
        while level < len(partials) and partials[level] is not None:
            contribution = [left+right for left, right in zip(partials[level], contribution)]
            partials[level] = None
            level += 1
        if level == len(partials):
            partials.append(contribution)
        else:
            partials[level] = contribution
        face_count += 1
    moments = [sum((block[j] for block in partials if block is not None), zero)
               for j in range(10)]
    if not face_count or moments[0] <= 0:
        raise GeometryValidationError('rational material has nonpositive volume')
    if arithmetic_backend == 'GMPY2_2_2_1_EXACT_V1':
        return tuple(Fraction(int(value.numerator),int(value.denominator)) for value in moments)
    return tuple(moments)


def compare_exact_rational_and_exported_moments(
    exact_vertices, exact_faces, exported_vertices, exported_faces, normalization_length,
) -> dict[str, object]:
    """Quantify Boolean output rounding using exact represented-mesh moments.

    A small moment difference is only a *numeric* observation; it cannot
    certify domain containment, exact set equality, or solver convergence.
    ``normalization_length`` is fixed from the input design domain, never
    selected from the exported result.
    """
    if isinstance(normalization_length, bool):
        raise GeometryValidationError('invalid input normalization length')
    try:
        length = Fraction(normalization_length)
    except (TypeError, ValueError, OverflowError) as error:
        raise GeometryValidationError('invalid input normalization length') from error
    if length <= 0:
        raise GeometryValidationError('input normalization length must be positive')
    exact = integrate_exact_rational_raw_moments(exact_vertices, exact_faces)
    exported = integrate_exact_rational_raw_moments(exported_vertices, exported_faces)
    scales = (length**3,) + (length**4,)*3 + (length**5,)*6
    errors = tuple(abs(a-b)/scale for a, b, scale in zip(exact, exported, scales))
    maximum = max(errors)
    return {
        'scope': 'ten-geometric-moments-only;not-exact-solid-equivalence',
        'order': ['M', 'LX', 'LY', 'LZ', 'QXX', 'QYY', 'QZZ', 'QYZ', 'QXZ', 'QXY'],
        'normalization_length_exact': fraction_to_rational_text(length),
        'maximum_absolute_normalized_difference_exact':
            fraction_to_rational_text(maximum),
        'maximum_absolute_normalized_difference': float(maximum),
        'within_1e_6_numeric_gate': maximum <= Fraction(1, 1000000),
    }


def _mass_properties_from_centered(
    mass: float,
    center_at_anchor: Vector3,
    centered_second: Matrix3,
    anchor: Vector3,
) -> MassProperties:
    if not math.isfinite(mass) or mass <= 0.0:
        raise GeometryValidationError("total material mass must be finite and positive")
    center = tuple(anchor[axis] + center_at_anchor[axis] for axis in range(3))
    trace = sum(centered_second[axis][axis] for axis in range(3))
    inertia = tuple(
        tuple(
            (trace if row == column else 0.0) - centered_second[row][column]
            for column in range(3)
        )
        for row in range(3)
    )
    first_at_anchor = tuple(mass * value for value in center_at_anchor)
    first = tuple(first_at_anchor[axis] + mass * anchor[axis] for axis in range(3))
    second_at_anchor = tuple(
        tuple(
            centered_second[row][column]
            + mass * center_at_anchor[row] * center_at_anchor[column]
            for column in range(3)
        )
        for row in range(3)
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
    prepared_regions: list[_PreparedRegion] = []
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
        mass = density * volume
        if not math.isfinite(mass) or mass <= 0.0:
            raise GeometryValidationError(
                f"region {region_index} material mass must be finite and positive"
            )
        center_local = tuple(value / volume for value in local_first)
        centered_second = tuple(
            tuple(
                density
                * (
                    local_second[row][column]
                    - volume * center_local[row] * center_local[column]
                )
                for column in range(3)
            )
            for row in range(3)
        )
        prepared_regions.append(
            _PreparedRegion(prepared, mass, center_local, centered_second)
        )

    dominant = max(prepared_regions, key=lambda prepared_region: prepared_region.mass)
    anchor = dominant.mesh.reference
    total_mass = math.fsum(region.mass for region in prepared_regions)
    centers_at_anchor = tuple(
        tuple(
            region.mesh.reference[axis]
            - anchor[axis]
            + region.center_local[axis]
            for axis in range(3)
        )
        for region in prepared_regions
    )
    center_at_anchor = tuple(
        math.fsum(
            region.mass * center[axis]
            for region, center in zip(prepared_regions, centers_at_anchor)
        )
        / total_mass
        for axis in range(3)
    )
    centered_terms = [[[] for _ in range(3)] for _ in range(3)]
    for region, region_center in zip(prepared_regions, centers_at_anchor):
        delta = tuple(
            region_center[axis] - center_at_anchor[axis] for axis in range(3)
        )
        for row in range(3):
            for column in range(3):
                centered_terms[row][column].append(
                    region.centered_second[row][column]
                    + region.mass * delta[row] * delta[column]
                )
    centered_second = tuple(
        tuple(
            math.fsum(centered_terms[row][column]) for column in range(3)
        )
        for row in range(3)
    )
    return _mass_properties_from_centered(
        total_mass,
        center_at_anchor,
        centered_second,
        anchor,
    )


def compute_mass_properties(mesh: TriangleMesh, density: float = 1.0) -> MassProperties:
    """Compute constant-density properties for one oriented surface mesh."""

    return aggregate_mass_properties((DensityRegion(mesh, density),))
