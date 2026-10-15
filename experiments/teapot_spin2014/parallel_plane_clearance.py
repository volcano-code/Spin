"""Exact conservative bound for distinct parallel triangle supporting planes.

This is NOT a wall-thickness definition or a global manufacturing certificate.
No projection-overlap, material-side, edge/vertex or nonparallel feature claim is
made. A small bound is inconclusive, not a thin-material witness. Coplanar pairs
cannot bound a positive-width parallel sheet and are counted separately.
"""
from fractions import Fraction
import math
from numbers import Integral, Real


def _number(value):
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError('finite real coordinates required')
    if not math.isfinite(value):
        raise ValueError('nonfinite coordinate or threshold')
    return Fraction(int(value)) if isinstance(value, Integral) else Fraction(float(value))


def audit_parallel_planes(vertices, faces, threshold):
    """Use represented binary64 values exactly; reject degenerate triangles.

    Plane normals are divided by their first nonzero component. Sorting offsets
    in each parallel family then checking adjacent distinct offsets gives the
    minimum positive plane separation in that family without quadratic storage.
    Triangle distance is at least supporting-plane distance. Thus a lower bound
    above threshold excludes subthreshold positive-width parallel face sheets,
    including pairs belonging to the same connected boundary component.
    """
    minimum = _number(threshold)
    if minimum <= 0:
        raise ValueError('positive threshold required')
    points = [tuple(_number(x) for x in row) for row in vertices]
    if not points or any(len(row) != 3 for row in points):
        raise ValueError('nonempty three-dimensional vertices required')
    planes = {}
    count = 0
    for ordinal, face in enumerate(faces):
        if len(face) != 3 or any(isinstance(i, bool) or not isinstance(i, Integral)
                                 or i < 0 or i >= len(points) for i in face):
            raise ValueError('triangle index invalid')
        a, b, c = [points[int(i)] for i in face]
        u, v = [b[i]-a[i] for i in range(3)], [c[i]-a[i] for i in range(3)]
        normal = (u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0])
        pivot = next((x for x in normal if x), None)
        if pivot is None:
            raise ValueError('degenerate triangle')
        direction = tuple(x/pivot for x in normal)
        offset = sum(x*y for x, y in zip(direction, a))
        offsets = planes.setdefault(direction, {})
        number, first = offsets.get(offset, (0, ordinal))
        offsets[offset] = (number+1, first)
        count += 1
    if not count:
        raise ValueError('nonempty triangles required')
    best = None
    witness = None
    distinct_pairs = coplanar_pairs = 0
    for normal, offsets in planes.items():
        population = sum(n for n, _ in offsets.values())
        coplanar = sum(n*(n-1)//2 for n, _ in offsets.values())
        coplanar_pairs += coplanar
        distinct_pairs += population*(population-1)//2 - coplanar
        ordered = sorted(offsets)
        norm2 = sum(x*x for x in normal)
        for left, right in zip(ordered, ordered[1:]):
            squared = (right-left)**2/norm2
            pair = sorted((offsets[left][1], offsets[right][1]))
            if best is None or squared < best or (squared == best and pair < witness):
                best, witness = squared, pair
    status = ('NO_DISTINCT_PARALLEL_PLANES' if best is None else
              'PARALLEL_PLANE_LOWER_BOUND_AT_LEAST_THRESHOLD' if best >= minimum**2 else
              'LOWER_BOUND_INCONCLUSIVE')
    return {'schema': 'parallel-support-plane-bound-v1', 'status': status,
            'scope': 'DISTINCT_PARALLEL_TRIANGLE_SUPPORTING_PLANES_ONLY',
            'complete_wall_thickness': 'UNKNOWN', 'geometry_admission': 'NOT_GRANTED',
            'thin_material_witness': False,
            'threshold_exact': str(minimum), 'triangle_count': count,
            'parallel_direction_count': len(planes), 'distinct_plane_triangle_pairs': distinct_pairs,
            'coplanar_triangle_pairs_excluded': coplanar_pairs,
            'minimum_squared_plane_separation_exact': None if best is None else str(best),
            'minimum_plane_separation_display': None if best is None else math.sqrt(float(best)),
            'supporting_plane_face_indices': witness,
            'unchecked': ['nonparallel faces', 'edge/vertex thin features',
                          'projection overlap', 'material-side classification',
                          'global topology and wall definition']}
