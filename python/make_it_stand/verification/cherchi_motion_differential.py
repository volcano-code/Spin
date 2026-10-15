"""Local motion differentials for supported Cherchi point constructions.

These routines transport a fixed Cherchi arrangement through smooth input
vertex motion. They do not certify Boolean solid equivalence or differentiable
behavior across a combinatorial/topology event. Unsupported constructions and
ill-conditioned intersections fail closed.
"""

from __future__ import annotations

from fractions import Fraction
from math import fsum

import numpy as np

from make_it_stand.verification.rational_text import fraction_from_rational_text, fraction_to_rational_text


class CherchiDifferentialUnavailable(ValueError):
    """The retained construction does not define a safe local differential."""


def _validated_rational_vertices(points, rational_vertices):
    if rational_vertices is None:
        return None
    try:
        rows = tuple(tuple(value if isinstance(value, Fraction) else
                           fraction_from_rational_text(value) for value in row)
                     for row in rational_vertices)
        if len(rows) != len(points) or any(len(row) != 3 for row in rows):
            raise ValueError('rational coordinate dimensions differ')
        display = np.asarray([[float(value) for value in row] for row in rows])
        if not np.isfinite(display).all() or not np.array_equal(display, points):
            raise ValueError('rational coordinates do not match their canonical display')
    except (TypeError, ValueError, OverflowError) as error:
        raise CherchiDifferentialUnavailable(str(error)) from error
    return rows


def _rational_area_vector(a, b, c):
    u = tuple(b[i] - a[i] for i in range(3))
    v = tuple(c[i] - a[i] for i in range(3))
    return (u[1]*v[2] - u[2]*v[1], u[2]*v[0] - u[0]*v[2], u[0]*v[1] - u[1]*v[0])


def _rational_component_volume(points, faces, indices):
    origin = points[int(faces[indices[0]][0])]
    total = Fraction(0)
    for index in indices:
        a, b, c = (tuple(points[int(vertex)][i] - origin[i] for i in range(3))
                   for vertex in faces[index])
        cross = _rational_area_vector((Fraction(0),)*3, b, c)
        total += sum(a[i] * cross[i] for i in range(3)) / 6
    return total


def _exact_binary64_triangle_area_nonzero(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> bool:
    """Test exact nonzero area of the represented binary64 coordinates.

    This slow path is used only when the floating cross product is below the
    ordinary scale-relative threshold. Fraction(float) captures each input
    binary64 value exactly, so this distinguishes true collinearity from a
    nonzero sliver without silently deleting or perturbing the triangle.
    """
    aa = tuple(Fraction(float(value)) for value in a)
    u = tuple(Fraction(float(b[i])) - aa[i] for i in range(3))
    v = tuple(Fraction(float(c[i])) - aa[i] for i in range(3))
    cross = (
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    )
    return any(component != 0 for component in cross)


def normalize_closed_component_orientation(
    vertices: object,
    faces: object,
    *, rational_vertices=None,
) -> tuple[np.ndarray, dict[str, object]]:
    """Coherently orient a closed, orientable triangle mesh without changing it.

    The input face order is preserved. For each connected component, edge
    adjacency determines relative face winding; the emitted component-wide
    signed-volume sign is retained so cavity shells remain inward-oriented.
    If local reorientation would reverse that sign, its global orientation is
    ambiguous and the input is rejected. Open/nonmanifold/nonorientable or
    numerically zero-volume components are rejected. Persist the returned flip
    list beside the original Boolean output; this is not an unrecorded repair.
    """
    points = np.asarray(vertices, dtype=np.float64)
    triangles = np.asarray(faces)
    if points.ndim != 2 or points.shape[1] != 3 or points.shape[0] < 4 or not np.isfinite(points).all():
        raise CherchiDifferentialUnavailable("orientation input vertices are invalid")
    if triangles.ndim != 2 or triangles.shape[1] != 3 or triangles.shape[0] < 4 or \
            triangles.dtype.kind not in "iu" or np.any(triangles < 0) or np.any(triangles >= len(points)):
        raise CherchiDifferentialUnavailable("orientation input faces are invalid")
    rational = _validated_rational_vertices(points, rational_vertices)

    edge_incidence: dict[tuple[int, int], list[tuple[int, int]]] = {}
    scale = max(1.0, float(np.linalg.norm(points, ord=np.inf)))
    area_floor = 64.0 * np.finfo(float).eps * scale * scale
    near_degenerate_faces: list[int] = []
    exact_areas = [] if rational is not None else None
    for face_index, face in enumerate(triangles):
        a, b, c = (int(index) for index in face)
        if len({a, b, c}) != 3:
            raise CherchiDifferentialUnavailable("orientation input contains a repeated face vertex")
        cross_norm = float(np.linalg.norm(np.cross(points[b] - points[a], points[c] - points[a])))
        if rational is not None:
            area = _rational_area_vector(rational[a], rational[b], rational[c])
            if not any(area):
                raise CherchiDifferentialUnavailable("orientation input contains a zero-area face")
            exact_areas.append(area)
            if cross_norm <= area_floor:
                near_degenerate_faces.append(face_index)
        elif cross_norm <= area_floor:
            if not _exact_binary64_triangle_area_nonzero(points[a], points[b], points[c]):
                raise CherchiDifferentialUnavailable("orientation input contains a zero-area face")
            near_degenerate_faces.append(face_index)
        for first, second in ((a, b), (b, c), (c, a)):
            key = (min(first, second), max(first, second))
            edge_incidence.setdefault(key, []).append((face_index, 1 if first < second else -1))
    boundary_edges = [edge for edge, incident in edge_incidence.items() if len(incident) == 1]
    nonmanifold_edges = [edge for edge, incident in edge_incidence.items() if len(incident) > 2]
    if boundary_edges or nonmanifold_edges:
        samples = boundary_edges[:8] or nonmanifold_edges[:8]
        raise CherchiDifferentialUnavailable(
            "orientation input is open or nonmanifold: "
            f"boundary_edges={len(boundary_edges)}; "
            f"nonmanifold_edges={len(nonmanifold_edges)}; "
            f"sample_edges={samples}"
        )

    adjacency: list[list[tuple[int, int]]] = [[] for _ in range(len(triangles))]
    for incident in edge_incidence.values():
        (first, first_direction), (second, second_direction) = incident
        # If emitted edge directions agree, exactly one adjacent face must flip.
        parity = int(first_direction == second_direction)
        adjacency[first].append((second, parity))
        adjacency[second].append((first, parity))

    assigned: list[int | None] = [None] * len(triangles)
    components: list[list[int]] = []
    for seed in range(len(triangles)):
        if assigned[seed] is not None:
            continue
        assigned[seed] = 0
        stack = [seed]
        component: list[int] = []
        while stack:
            current = stack.pop()
            component.append(current)
            for neighbor, parity in adjacency[current]:
                expected = int(assigned[current]) ^ parity
                if assigned[neighbor] is None:
                    assigned[neighbor] = expected
                    stack.append(neighbor)
                elif assigned[neighbor] != expected:
                    raise CherchiDifferentialUnavailable("orientation input is nonorientable")
        components.append(component)

    oriented = np.asarray(triangles, dtype=np.int64).copy()
    for face_index, should_flip in enumerate(assigned):
        if should_flip:
            oriented[face_index, [1, 2]] = oriented[face_index, [2, 1]]

    def signed_volume(face_indices: list[int], faces_array: np.ndarray) -> float:
        # Translating each closed component to a local origin preserves its
        # signed volume while avoiding cancellation from a distant global
        # origin.  fsum keeps the final accumulation stable for mixed signs.
        component_vertices = np.unique(faces_array[face_indices].reshape(-1))
        origin = points[int(component_vertices[0])]
        terms = []
        for index in face_indices:
            a, b, c = (points[int(vertex)] - origin for vertex in faces_array[index])
            terms.append(float(np.dot(a, np.cross(b, c)) / 6.0))
        return fsum(terms)

    component_audits: list[dict[str, object]] = []
    for component in components:
        component_vertices = np.unique(triangles[component].reshape(-1))
        component_points = points[component_vertices]
        component_extent = (np.asarray([float(max(rational[int(v)][axis] for v in component_vertices) -
            min(rational[int(v)][axis] for v in component_vertices)) for axis in range(3)])
            if rational is not None else np.ptp(component_points, axis=0))
        component_scale = float(np.max(component_extent))
        if not np.isfinite(component_scale) or component_scale <= 0.0:
            raise CherchiDifferentialUnavailable(
                "orientation component has zero or nonfinite spatial extent"
            )
        volume_floor = (
            256.0 * np.finfo(float).eps * component_scale**3 * max(1, len(component))
        )
        exact_emitted = exact_coherent = None
        if rational is not None:
            # det(a-o,b-o,c-o) = (a-o).cross(b-a,c-a). A local
            # winding flip negates this same exact term; do not recompute area.
            reference = rational[int(triangles[component[0]][0])]
            emitted_sum = _ExactFractionSum()
            has_local_flips = any(assigned[index] for index in component)
            coherent_sum = _ExactFractionSum() if has_local_flips else None
            for index in component:
                first = rational[int(triangles[index][0])]
                term = sum((first[axis]-reference[axis])*exact_areas[index][axis]
                           for axis in range(3))/6
                emitted_sum.add(term)
                if coherent_sum is not None:
                    coherent_sum.add(-term if assigned[index] else term)
            exact_emitted = emitted_sum.total()
            exact_coherent = exact_emitted if coherent_sum is None else coherent_sum.total()
            if exact_emitted == 0 or exact_coherent == 0:
                raise CherchiDifferentialUnavailable("orientation component has zero exact signed volume")
            emitted_volume, coherent_volume = float(exact_emitted), float(exact_coherent)
            if not np.isfinite(emitted_volume) or not np.isfinite(coherent_volume) or emitted_volume == 0 or coherent_volume == 0:
                raise CherchiDifferentialUnavailable("exact signed volume cannot be represented finitely")
        else:
            emitted_volume = signed_volume(component, triangles)
            coherent_volume = signed_volume(component, oriented)
            if abs(emitted_volume) <= volume_floor or abs(coherent_volume) <= volume_floor:
                raise CherchiDifferentialUnavailable("orientation component has numerically ambiguous signed volume")
        if (emitted_volume < 0.0) != (coherent_volume < 0.0):
            raise CherchiDifferentialUnavailable(
                "orientation component global sign is ambiguous after relative reorientation"
            )
        final_volume = coherent_volume
        local_flips = sum(int(assigned[index] or 0) for index in component)
        component_audits.append({
            "face_count": len(component),
            "component_scale": component_scale,
            "signed_volume_floor": volume_floor,
            "emitted_signed_volume": emitted_volume,
            "normalized_signed_volume": final_volume,
            "orientation_sign": 1 if final_volume > 0.0 else -1,
            "relative_flip_count": local_flips,
            "component_orientation_preserved": True,
            "relative_orientation_correction": abs(final_volume - emitted_volume) / abs(final_volume),
            **({'exact_signed_volume': fraction_to_rational_text(exact_coherent),
                'exact_emitted_signed_volume': fraction_to_rational_text(exact_emitted)}
               if rational is not None else {}),
        })

    flipped = np.any(oriented != np.asarray(triangles, dtype=np.int64), axis=1)
    # Recheck the returned contract explicitly; callers can keep using the
    # stricter differential validator without weakening its preconditions.
    for edge_key, incident in edge_incidence.items():
        directions = []
        for face_index, _ in incident:
            face = oriented[face_index]
            for first, second in ((int(face[0]), int(face[1])),
                                  (int(face[1]), int(face[2])),
                                  (int(face[2]), int(face[0]))):
                if (min(first, second), max(first, second)) == edge_key:
                    directions.append(1 if first < second else -1)
        if len(directions) != 2 or sum(directions) != 0:
            raise CherchiDifferentialUnavailable("orientation normalization failed its directed-edge audit")

    return oriented, {
        "schema": "closed-component-orientation-normalization-v1",
        "input_faces_preserved": True,
        "vertex_positions_changed": False,
        "component_count": len(components),
        "flipped_face_count": int(np.count_nonzero(flipped)),
        "flipped_face_indices": np.flatnonzero(flipped).astype(int).tolist(),
        "near_degenerate_face_indices": near_degenerate_faces,
        "near_degenerate_predicate": ("exact_rational_area_nonzero_v1" if rational is not None else
                                      "exact_binary64_dyadic_area_nonzero_v1"),
        "coordinate_authority": "RATIONAL_BOUNDARY" if rational is not None else "BINARY64_BOUNDARY",
        "floating_export_approved": False,
        "components": component_audits,
    }


def _validated_motion(positions: object, velocities: object, expected: int):
    points = np.asarray(positions, dtype=np.float64)
    motion = np.asarray(velocities, dtype=np.float64)
    if points.shape != (expected, 3) or not np.isfinite(points).all():
        raise CherchiDifferentialUnavailable("construction point positions have invalid shape or values")
    if motion.ndim != 3 or motion.shape[:2] != (expected, 3) or motion.shape[2] <= 0:
        raise CherchiDifferentialUnavailable("construction point velocities have invalid dimensions")
    if not np.isfinite(motion).all():
        raise CherchiDifferentialUnavailable("construction point velocities are nonfinite")
    return points, motion


class _ExactFractionSum:
    """Streaming balanced sum: exact arithmetic, O(log n) partial sums.

    Each occupied slot contains a power-of-two-sized block. Balancing avoids
    repeatedly adding small contributions to an ever-growing denominator;
    unlike a floating reduction, it cannot change the final rational value.
    """

    def __init__(self):
        self._partials = []

    def add(self, value):
        level = 0
        while level < len(self._partials) and self._partials[level] is not None:
            value = self._partials[level] + value
            self._partials[level] = None
            level += 1
        if level == len(self._partials):
            self._partials.append(value)
        else:
            self._partials[level] = value

    def total(self):
        return sum((value for value in self._partials if value is not None), Fraction(0))


def integrate_rational_polyhedral_raw_moments(vertices, faces, density=1.0, *, rational_vertices):
    """Primary exact surface-integral map; not the independent tetrahedron checker.

    Integrate div(r*f)/(3+degree(f)) for the ten mass monomials. This
    preserves rational geometry through integration, rounding only the final
    moments. Local orientation/topology checks are not solid embedding,
    containment, Boolean identity or manufacturing certification.
    """
    points = np.asarray(vertices, dtype=float)
    triangles = np.asarray(faces)
    if points.ndim != 2 or points.shape[1:] != (3,) or not np.isfinite(points).all():
        raise CherchiDifferentialUnavailable('surface moment vertices are invalid')
    if not np.isfinite(density) or density <= 0:
        raise CherchiDifferentialUnavailable('surface moment density must be finite and positive')
    rational = _validated_rational_vertices(points, rational_vertices)
    if rational is None:
        raise CherchiDifferentialUnavailable('surface moment map requires rational authority')
    oriented, _ = normalize_closed_component_orientation(points, triangles, rational_vertices=rational)
    if not np.array_equal(oriented, triangles):
        raise CherchiDifferentialUnavailable('surface moment boundary is not consistently oriented')
    moments = [_ExactFractionSum() for _ in range(10)]
    pairs = ((0,0),(1,1),(2,2),(1,2),(0,2),(0,1))
    for face in triangles:
        vertices = tuple(rational[int(index)] for index in face)
        area = _rational_area_vector(*vertices)
        sums = tuple(sum(point[c] for point in vertices) for c in range(3))
        second = {(p,q): sum(point[p]*point[q] for point in vertices)
                  for p in range(3) for q in range(p,3)}
        third = {(p,q,c): sum(point[p]*point[q]*point[c] for point in vertices)
                 for p in range(3) for q in range(p,3) for c in range(q,3)}
        def d(p,q):
            return second[tuple(sorted((p,q)))]
        moments[0].add(sum(area[c]*sums[c] for c in range(3)) / 18)
        for p in range(3):
            moments[p+1].add(sum(area[c]*(sums[p]*sums[c]+d(p,c)) for c in range(3)) / 96)
        for m,(p,q) in enumerate(pairs,4):
            moments[m].add(sum(area[c]*(sums[p]*sums[q]*sums[c]+sums[p]*d(q,c)+
                sums[q]*d(p,c)+sums[c]*d(p,q)+2*third[tuple(sorted((p,q,c)))])
                for c in range(3)) / 600)
    rho = Fraction(float(density))
    result = np.asarray([float(value.total()*rho) for value in moments])
    if not np.isfinite(result).all() or result[0] <= 0:
        raise CherchiDifferentialUnavailable('surface moment map has nonpositive or nonfinite mass')
    return result


def differentiate_cherchi_point(
    point_type: str,
    operand_positions: object,
    operand_velocities: object,
    *,
    parameter_t: float = 0.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Return the represented point and its d-dimensional local velocity.

    Supported source relations are EXPLICIT3D, LPI (line-plane intersection),
    TPI (three-plane intersection), and LNC (fixed-parameter linear
    combination). LNC's parameter is treated as fixed within one arrangement
    stratum, matching the pinned Cherchi construction object; topology changes
    still require a newly evaluated construction and differential.
    """
    kind = str(point_type).upper()
    counts = {"EXPLICIT3D": 1, "LPI": 5, "TPI": 9, "LNC": 2}
    if kind not in counts:
        raise CherchiDifferentialUnavailable(f"unsupported Cherchi point construction: {kind}")
    points, motion = _validated_motion(operand_positions, operand_velocities, counts[kind])
    dimensions = motion.shape[2]
    scale = max(1.0, float(np.linalg.norm(points, ord=np.inf)))
    if kind == "EXPLICIT3D":
        return points[0].copy(), motion[0].copy()

    if kind == "LNC":
        if isinstance(parameter_t, (bool, np.bool_)) or not isinstance(
            parameter_t, (int, float, np.integer, np.floating),
        ) or not np.isfinite(parameter_t) or not 0.0 <= parameter_t <= 1.0:
            raise CherchiDifferentialUnavailable("LNC parameter must be finite and in [0, 1]")
        parameter = float(parameter_t)
        point = (1.0 - parameter) * points[0] + parameter * points[1]
        velocity = (1.0 - parameter) * motion[0] + parameter * motion[1]
        if not np.isfinite(point).all() or not np.isfinite(velocity).all():
            raise CherchiDifferentialUnavailable("LNC motion produced nonfinite output")
        return point, velocity

    if kind == "LPI":
        line_a, line_b, plane_a, plane_b, plane_c = points
        va, vb, vpa, vpb, vpc = motion
        direction = line_b - line_a
        normal = np.cross(plane_b - plane_a, plane_c - plane_a)
        denominator = float(normal @ direction)
        denominator_scale = float(np.linalg.norm(normal) * np.linalg.norm(direction))
        if denominator_scale == 0.0 or abs(denominator) <= 128.0 * np.finfo(float).eps * denominator_scale:
            raise CherchiDifferentialUnavailable("LPI line is parallel or ill-conditioned against its source plane")
        numerator = float(normal @ (plane_a - line_a))
        parameter = numerator / denominator
        point = line_a + parameter * direction
        velocity = np.empty((3, dimensions), dtype=np.float64)
        for axis in range(dimensions):
            dnormal = np.cross(vpb[:, axis] - vpa[:, axis], plane_c - plane_a) + \
                np.cross(plane_b - plane_a, vpc[:, axis] - vpa[:, axis])
            ddirection = vb[:, axis] - va[:, axis]
            dnumerator = float(dnormal @ (plane_a - line_a) + normal @ (vpa[:, axis] - va[:, axis]))
            ddenominator = float(dnormal @ direction + normal @ ddirection)
            dparameter = (dnumerator * denominator - numerator * ddenominator) / (denominator * denominator)
            velocity[:, axis] = va[:, axis] + dparameter * direction + parameter * ddirection
        return point, velocity

    # TPI: each consecutive triple defines one moving plane n_i.x = b_i.
    matrix = np.empty((3, 3), dtype=np.float64)
    rhs = np.empty(3, dtype=np.float64)
    normals: list[np.ndarray] = []
    offsets: list[np.ndarray] = []
    for plane_index in range(3):
        start = plane_index * 3
        a, b, c = points[start:start + 3]
        normal = np.cross(b - a, c - a)
        normals.append(normal)
        offsets.append(normal @ a)
        matrix[plane_index] = normal
        rhs[plane_index] = offsets[-1]
    singular_values = np.linalg.svd(matrix, compute_uv=False)
    if singular_values[0] == 0.0 or singular_values[-1] <= 1e-12 * singular_values[0]:
        raise CherchiDifferentialUnavailable("TPI source planes are singular or ill-conditioned")
    point = np.linalg.solve(matrix, rhs)
    velocity = np.empty((3, dimensions), dtype=np.float64)
    for axis in range(dimensions):
        dmatrix = np.empty((3, 3), dtype=np.float64)
        drhs = np.empty(3, dtype=np.float64)
        for plane_index in range(3):
            start = plane_index * 3
            a, b, c = points[start:start + 3]
            va, vb, vc = motion[start:start + 3, :, axis]
            dnormal = np.cross(vb - va, c - a) + np.cross(b - a, vc - va)
            dmatrix[plane_index] = dnormal
            drhs[plane_index] = dnormal @ a + normals[plane_index] @ va
        velocity[:, axis] = np.linalg.solve(matrix, drhs - dmatrix @ point)
    if not np.isfinite(point).all() or not np.isfinite(velocity).all():
        raise CherchiDifferentialUnavailable("TPI motion solve produced nonfinite output")
    return point, velocity


def _raw_monomials(point: np.ndarray) -> np.ndarray:
    x, y, z = point
    return np.asarray((1.0, x, y, z, x*x, y*y, z*z, y*z, x*z, x*y))


def differentiate_polyhedral_raw_moments(
    vertices: object,
    faces: object,
    vertex_velocities: object,
    density: float = 1.0,
    *, rational_vertices=None,
) -> np.ndarray:
    """Differentiate ten raw moments of a closed oriented polyhedral boundary.

    The Reynolds transport integrand is integrated with a degree-three exact
    triangle rule. This is a local fixed-topology boundary differential; it
    does not audit self-intersection or prove a Boolean's set semantics.
    """
    points = np.asarray(vertices, dtype=np.float64)
    triangles = np.asarray(faces)
    motion = np.asarray(vertex_velocities, dtype=np.float64)
    if points.ndim != 2 or points.shape[1] != 3 or points.shape[0] < 4 or not np.isfinite(points).all():
        raise CherchiDifferentialUnavailable("polyhedral vertices are invalid")
    if triangles.ndim != 2 or triangles.shape[1] != 3 or triangles.shape[0] < 4 or \
            triangles.dtype.kind not in "iu" or np.any(triangles < 0) or np.any(triangles >= len(points)):
        raise CherchiDifferentialUnavailable("polyhedral faces are invalid")
    if motion.ndim != 3 or motion.shape[:2] != points.shape or motion.shape[2] <= 0 or \
            not np.isfinite(motion).all():
        raise CherchiDifferentialUnavailable("polyhedral vertex motion is invalid")
    if not np.isfinite(density) or density <= 0.0:
        raise CherchiDifferentialUnavailable("polyhedral density must be finite and positive")
    rational = _validated_rational_vertices(points, rational_vertices)

    # A differential of a material volume is valid only for a closed oriented
    # two-manifold. Check every directed edge pair before integrating.
    edge_directions: dict[tuple[int, int], list[int]] = {}
    scale = max(1.0, float(np.linalg.norm(points, ord=np.inf)))
    area_floor = 64.0 * np.finfo(float).eps * scale * scale
    for face in triangles:
        if len(set(int(index) for index in face)) != 3:
            raise CherchiDifferentialUnavailable("polyhedral boundary contains a repeated face vertex")
        a, b, c = (int(index) for index in face)
        cross_norm = float(np.linalg.norm(np.cross(points[b] - points[a], points[c] - points[a])))
        if rational is not None:
            if not any(_rational_area_vector(rational[a], rational[b], rational[c])):
                raise CherchiDifferentialUnavailable("polyhedral boundary contains a zero-area face")
        elif cross_norm <= area_floor and not _exact_binary64_triangle_area_nonzero(
                points[a], points[b], points[c]):
            raise CherchiDifferentialUnavailable("polyhedral boundary contains a zero-area face")
        for first, second in ((a, b), (b, c), (c, a)):
            key = (min(first, second), max(first, second))
            edge_directions.setdefault(key, []).append(1 if first < second else -1)
    boundary_edges = [edge for edge, directions in edge_directions.items()
                      if len(directions) == 1]
    nonmanifold_edges = [edge for edge, directions in edge_directions.items()
                         if len(directions) > 2]
    inconsistent_edges = [edge for edge, directions in edge_directions.items()
                          if len(directions) == 2 and sum(directions) != 0]
    if boundary_edges or nonmanifold_edges or inconsistent_edges:
        sample_edges = (boundary_edges + nonmanifold_edges + inconsistent_edges)[:8]
        raise CherchiDifferentialUnavailable(
            "polyhedral boundary is not closed and consistently oriented: "
            f"boundary_edges={len(boundary_edges)}; "
            f"nonmanifold_edges={len(nonmanifold_edges)}; "
            f"inconsistent_edges={len(inconsistent_edges)}; "
            f"sample_edges={sample_edges}"
        )

    if rational is not None:
        # Exact geometry coefficients of Reynolds transport. For normalized
        # triangle barycentrics, E[lambda_j]=1/3,
        # E[lambda_j*x]=(sum(x)+x_j)/12, and
        # E[lambda_j*x*y]=(sum(x)sum(y)+sum(x*y)+x_j sum(y)
        #                  +y_j sum(x)+2*x_j*y_j)/60.
        # Geometry is rational; inherited source velocities and returned
        # Newton Jacobian remain binary64. Round coefficients only after the
        # exact face integral, never replace the authority with its display.
        contributions = [[[] for _ in range(motion.shape[2])] for _ in range(10)]
        pairs = ((0,0),(1,1),(2,2),(1,2),(0,2),(0,1))
        rho = Fraction(float(density))
        for face in triangles:
            face_points = tuple(rational[int(vertex)] for vertex in face)
            area = _rational_area_vector(*face_points)
            sums = tuple(sum(point[c] for point in face_points) for c in range(3))
            pair_sums = tuple(sum(point[p]*point[q] for point in face_points) for p,q in pairs)
            coefficients = np.empty((3,10,3), dtype=float)
            for j, point in enumerate(face_points):
                expectations = (Fraction(1,3), *( (sums[c]+point[c])/12 for c in range(3)),
                    *((sums[p]*sums[q]+pair_sum+point[p]*sums[q]+point[q]*sums[p]+2*point[p]*point[q])/60
                      for (p,q),pair_sum in zip(pairs,pair_sums)))
                for m, expected in enumerate(expectations):
                    for c in range(3):
                        value = rho * area[c] * expected / 2
                        rounded = float(value)
                        if not np.isfinite(rounded) or (value and rounded == 0):
                            raise CherchiDifferentialUnavailable('rational transport coefficient cannot be represented finitely')
                        coefficients[j,m,c] = rounded
            term = np.einsum('jmc,jcd->md', coefficients, motion[face])
            for m in range(10):
                for d in range(motion.shape[2]):
                    contributions[m][d].append(float(term[m,d]))
        result = np.asarray([[fsum(values) for values in row] for row in contributions])
        if not np.isfinite(result).all():
            raise CherchiDifferentialUnavailable('rational polyhedral moment differential is nonfinite')
        return result

    barycentric = np.asarray(((1/3, 1/3, 1/3), (0.6, 0.2, 0.2),
                              (0.2, 0.6, 0.2), (0.2, 0.2, 0.6)))
    quadrature_weights = np.asarray((-27/48, 25/48, 25/48, 25/48))
    jacobian = np.zeros((10, motion.shape[2]), dtype=np.float64)
    for face in triangles:
        a, b, c = (int(index) for index in face)
        face_points = points[[a, b, c]]
        face_motion = motion[[a, b, c]]
        area_vector = np.cross(face_points[1] - face_points[0], face_points[2] - face_points[0])
        for bary, weight in zip(barycentric, quadrature_weights):
            point = bary @ face_points
            point_motion = np.einsum("i,ijd->jd", bary, face_motion)
            flux = area_vector @ point_motion
            jacobian += (0.5 * density * weight * flux) * _raw_monomials(point)[:, None]
    if not np.isfinite(jacobian).all():
        raise CherchiDifferentialUnavailable("polyhedral moment differential is nonfinite")
    return jacobian


def recover_construction_vertex_motion(
    records: object,
    design_vertices: object,
    source_vertices: object,
    source_vertex_velocity: object,
    output_vertices: object,
    *,
    coordinate_tolerance: float = 1e-10,
) -> np.ndarray:
    """Recover output-vertex motion only from unique retained Cherchi operands.

    Input object 0 is the fixed design mesh and object 1 is the moving source
    surface. LNC is differentiated only with its retained fixed interpolation
    parameter; the parameter is part of the local construction identity.
    """
    design = np.asarray(design_vertices, dtype=np.float64)
    source = np.asarray(source_vertices, dtype=np.float64)
    source_motion = np.asarray(source_vertex_velocity, dtype=np.float64)
    output = np.asarray(output_vertices, dtype=np.float64)
    if design.ndim != 2 or design.shape[1] != 3 or not np.isfinite(design).all():
        raise CherchiDifferentialUnavailable("design source vertices are invalid")
    if source.ndim != 2 or source.shape[1] != 3 or not np.isfinite(source).all():
        raise CherchiDifferentialUnavailable("moving source vertices are invalid")
    if source_motion.ndim != 3 or source_motion.shape[:2] != source.shape or \
            source_motion.shape[2] <= 0 or not np.isfinite(source_motion).all():
        raise CherchiDifferentialUnavailable("moving source vertex Jacobian is invalid")
    if output.ndim != 2 or output.shape[1] != 3 or not np.isfinite(output).all():
        raise CherchiDifferentialUnavailable("Boolean output vertices are invalid")
    if not isinstance(records, (list, tuple)) or len(records) != len(output):
        raise CherchiDifferentialUnavailable("construction records do not match Boolean output")
    if not np.isfinite(coordinate_tolerance) or coordinate_tolerance <= 0:
        raise CherchiDifferentialUnavailable("coordinate tolerance must be finite and positive")

    scale = max(1.0, float(np.max(np.abs(np.vstack((design, source, output))))))
    result = np.empty((len(output), 3, source_motion.shape[2]), dtype=np.float64)
    for index, record in enumerate(records):
        if not isinstance(record, dict) or record.get("vertex_index") != index:
            raise CherchiDifferentialUnavailable("construction records are malformed or reordered")
        point_type = record.get("point_type")
        operands = record.get("operands")
        coordinates = record.get("operand_coordinates")
        expected_count = {"EXPLICIT3D": 1, "LPI": 5, "TPI": 9, "LNC": 2}.get(point_type)
        if expected_count is None:
            raise CherchiDifferentialUnavailable(f"unsupported Cherchi construction type: {point_type}")
        if not isinstance(operands, (list, tuple)) or len(operands) != expected_count or \
                not isinstance(coordinates, (list, tuple)) or len(coordinates) != expected_count:
            raise CherchiDifferentialUnavailable("construction operands or coordinates are incomplete")
        points = np.empty((expected_count, 3), dtype=np.float64)
        motions = np.empty((expected_count, 3, source_motion.shape[2]), dtype=np.float64)
        for operand_index, operand in enumerate(operands):
            if not isinstance(operand, dict) or operand.get("status") != "UNIQUE":
                raise CherchiDifferentialUnavailable("Cherchi operand source is ambiguous or unmatched")
            candidates = operand.get("source_vertices")
            if not isinstance(candidates, (list, tuple)) or len(candidates) != 1:
                raise CherchiDifferentialUnavailable("Cherchi operand lacks a unique source vertex")
            candidate = candidates[0]
            if not isinstance(candidate, dict):
                raise CherchiDifferentialUnavailable("Cherchi source identity is malformed")
            object_index, vertex_index = candidate.get("input_object"), candidate.get("input_vertex")
            if not isinstance(object_index, int) or not isinstance(vertex_index, int):
                raise CherchiDifferentialUnavailable("Cherchi source identity is not integral")
            point = np.asarray(coordinates[operand_index], dtype=np.float64)
            if point.shape != (3,) or not np.isfinite(point).all():
                raise CherchiDifferentialUnavailable("Cherchi operand coordinate is invalid")
            if object_index == 0 and 0 <= vertex_index < len(design):
                points[operand_index] = design[vertex_index]
                motions[operand_index] = 0.0
            elif object_index == 1 and 0 <= vertex_index < len(source):
                points[operand_index] = source[vertex_index]
                motions[operand_index] = source_motion[vertex_index]
            else:
                raise CherchiDifferentialUnavailable("Cherchi source references an unknown input mesh vertex")
            if np.linalg.norm(points[operand_index] - point, ord=np.inf) > coordinate_tolerance * scale:
                raise CherchiDifferentialUnavailable("operand coordinate disagrees with its mapped source vertex")

        if point_type == "EXPLICIT3D":
            reconstructed, velocity = points[0], motions[0]
        elif point_type == "LPI":
            reconstructed, velocity = differentiate_cherchi_point("LPI", points, motions)
        elif point_type == "LNC":
            parameter = record.get("parameter_t")
            if isinstance(parameter, bool) or not isinstance(parameter, (int, float)):
                raise CherchiDifferentialUnavailable("LNC construction lacks a finite retained parameter")
            reconstructed, velocity = differentiate_cherchi_point(
                "LNC", points, motions, parameter_t=float(parameter),
            )
        else:
            reconstructed, velocity = differentiate_cherchi_point("TPI", points, motions)
        if np.linalg.norm(reconstructed - output[index], ord=np.inf) > coordinate_tolerance * scale:
            raise CherchiDifferentialUnavailable("reconstructed Cherchi point disagrees with output vertex")
        result[index] = velocity
    return result
