"""Exact dyadic shared-edge refinement for the registered planar benchmark."""
from fractions import Fraction
import math


def refine_mesh(vertices, faces, *, rounds=0, split_longest=False):
    if type(rounds) is not int or rounds < 0:
        raise ValueError('INVALID_REFINEMENT_ROUNDS')
    points = [[float(x) for x in p] for p in vertices]
    triangles = [[int(i) for i in f] for f in faces]
    if any(len(p) != 3 or not all(math.isfinite(x) for x in p) for p in points):
        raise ValueError('INVALID_VERTICES')
    if any(len(f) != 3 or len(set(f)) != 3 or any(i < 0 or i >= len(points) for i in f) for f in triangles):
        raise ValueError('INVALID_TRIANGLES')

    def midpoint_factory():
        cache = {}

        def midpoint(u, v):
            key = tuple(sorted((u, v)))
            if key not in cache:
                exact = [(Fraction(a) + Fraction(b)) / 2 for a, b in zip(points[u], points[v])]
                floating = [float(x) for x in exact]
                if any(Fraction(y) != x for x, y in zip(exact, floating)):
                    raise ValueError('MIDPOINT_NOT_EXACT_BINARY64')
                cache[key] = len(points)
                points.append(floating)
            return cache[key]
        return midpoint

    for _ in range(rounds):
        midpoint = midpoint_factory()
        refined = []
        for a, b, c in triangles:
            ab, bc, ca = midpoint(a, b), midpoint(b, c), midpoint(c, a)
            refined.extend(([a, ab, ca], [ab, b, bc], [ca, bc, c], [ab, bc, ca]))
        triangles = refined

    if split_longest:
        exact = [tuple(Fraction(x) for x in p) for p in points]
        selected = set()
        for face in triangles:
            edges = [tuple(sorted((face[i], face[(i + 1) % 3]))) for i in range(3)]
            edge = max(edges, key=lambda uv: (sum((a - b) ** 2 for a, b in zip(exact[uv[0]], exact[uv[1]])), uv))
            selected.add(edge)
        midpoint = midpoint_factory()
        refined = []
        for face in triangles:
            marked = [i for i in range(3) if tuple(sorted((face[i], face[(i + 1) % 3]))) in selected]
            # Registered cube triangles pair along the hypotenuse. Other meshes
            # needing two- or three-edge conformity are outside this fixture helper.
            if len(marked) != 1:
                raise ValueError('LONGEST_EDGE_PAIRING_UNSUPPORTED')
            i = marked[0]
            a, b, c = face[i], face[(i + 1) % 3], face[(i + 2) % 3]
            m = midpoint(a, b)
            refined.extend(([a, m, c], [m, b, c]))
        triangles = refined
    return points, triangles
