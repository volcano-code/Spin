"""Exact local certificate for two triangulations of a convex quadrilateral.

This checks one represented patch at ONE parameter value. It does not assert
that materials at different parameters are equal or certify differentiability,
source provenance, an entire solid, solver convergence or manufacturing.
"""
from collections import Counter
from fractions import Fraction

from .rational_text import fraction_to_rational_text


def _sub(a, b):
    return tuple(x-y for x, y in zip(a, b))


def _cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


def _dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def certify_convex_quad_retriangulation(points, first_faces, second_faces):
    """Return a rational certificate, or reject an unsupported/unequal patch.

    Both inputs must be two triangles covering the same strictly convex,
    coplanar four-corner patch with equal orientation. Integer/Fraction
    coordinates are accepted; floats are deliberately not exact authority.
    The common boundary and positive triangle areas prove complete coverage;
    the internal edge must occur once in each direction. No epsilon is used.
    """
    if not points or any(len(p) != 3 or any(
            isinstance(x, bool) or not isinstance(x, (int, Fraction)) for x in p)
            for p in points):
        raise ValueError('exact rational patch coordinates required')
    p = [tuple(Fraction(x) for x in row) for row in points]

    def triangles(rows):
        if len(rows) != 2 or any(len(f) != 3 for f in rows):
            raise ValueError('two triangles per quadrilateral patch required')
        if any(type(i) is not int or not 0 <= i < len(p) for f in rows for i in f):
            raise ValueError('invalid patch index')
        if any(len(set(f)) != 3 for f in rows):
            raise ValueError('degenerate patch triangle')
        return [tuple(f) for f in rows]

    first, second = triangles(first_faces), triangles(second_faces)
    corners = set(i for f in first for i in f)
    if len(corners) != 4 or set(i for f in second for i in f) != corners:
        raise ValueError('same four patch vertices required')
    a, b, c = (p[i] for i in first[0])
    n = _cross(_sub(b, a), _sub(c, a))
    if not any(n):
        raise ValueError('degenerate patch triangle')
    if any(_dot(n, _sub(p[i], a)) != 0 for i in corners):
        raise ValueError('patch vertices are not exactly coplanar')

    def boundary(rows):
        edges = Counter((f[k], f[(k+1) % 3]) for f in rows for k in range(3))
        if any(count != 1 for count in edges.values()):
            raise ValueError('duplicate/overlapping patch triangle')
        internal = [(i, j) for i, j in edges if (j, i) in edges]
        outer = [(i, j) for i, j in edges if (j, i) not in edges]
        if len(internal) != 2 or len(outer) != 4:
            raise ValueError('invalid oriented quadrilateral boundary')
        outgoing = dict(outer)
        if len(outgoing) != 4 or len(set(outgoing.values())) != 4:
            raise ValueError('nonmanifold patch boundary')
        cycle = [min(corners)]
        for _ in range(3):
            cycle.append(outgoing[cycle[-1]])
        if set(cycle) != corners or outgoing[cycle[-1]] != cycle[0]:
            raise ValueError('patch boundary is not a four-corner cycle')
        return cycle

    cycle = boundary(first)
    if boundary(second) != cycle:
        raise ValueError('patch orientation or boundary differs')
    for k in range(4):
        x, y, z = (p[cycle[(k+j) % 4]] for j in range(3))
        if _dot(n, _cross(_sub(y, x), _sub(z, y))) <= 0:
            raise ValueError('patch must be strictly convex')

    def area_vector(rows):
        total = [Fraction(0)]*3
        for f in rows:
            x, y, z = (p[i] for i in f)
            normal = _cross(_sub(y, x), _sub(z, x))
            if _dot(n, normal) <= 0:
                raise ValueError('degenerate or reversed patch triangle')
            total = [v+w for v, w in zip(total, normal)]
        return total

    area = area_vector(first)
    if area_vector(second) != area:
        raise ValueError('patch area/coverage differs')
    pivot = next(x for x in n if x)
    plane = [x/pivot for x in (*n, -_dot(n, a))]
    return dict(schema='exact-convex-quad-retriangulation-v1',
                boundary_cycle=cycle,
                plane=[fraction_to_rational_text(x) for x in plane],
                twice_area_vector=[fraction_to_rational_text(x) for x in area])
