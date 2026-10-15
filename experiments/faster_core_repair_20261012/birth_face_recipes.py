"""Exact incidence recipes for same-call CGAL tapes, not solver admission.

Birth faces restrict candidates. Every defining edge/triangle must contain
the exact material point. No proximity matching or floating predicates.
"""
from exact_planar_patches import subtract, cross, dot, point_on_triangle


def replay_vertex_recipe(recipe, operands):
    if recipe[0] == 'EXPLICIT':
        return tuple(operands[recipe[1]][recipe[2]])
    if recipe[0] != 'EDGE_PLANE':
        raise ValueError('UNSUPPORTED_VERTEX_RECIPE')
    _, edge_object, edge, plane_object, plane = recipe
    a, b = [operands[edge_object][i] for i in edge]
    p, q, r = [operands[plane_object][i] for i in plane]
    n = cross(subtract(q, p), subtract(r, p))
    direction = subtract(b, a)
    denominator = dot(n, direction)
    if denominator == 0:
        raise ValueError('SINGULAR_EDGE_PLANE_RECIPE')
    t = dot(n, subtract(p, a)) / denominator
    if not 0 <= t <= 1:
        raise ValueError('RECIPE_LEFT_SOURCE_EDGE')
    result = tuple(x + t * d for x, d in zip(a, direction))
    if not point_on_triangle(result, (p, q, r)):
        raise ValueError('RECIPE_LEFT_SOURCE_TRIANGLE')
    return result


def recover_vertex_recipe(point, birth_faces, operands, operand_faces):
    point = tuple(point)
    facets = []
    for birth in sorted(set(birth_faces)):
        obj = int(birth >= len(operand_faces[0]))
        index = birth - (len(operand_faces[0]) if obj else 0)
        if not 0 <= index < len(operand_faces[obj]):
            raise ValueError('BIRTH_FACE_OUT_OF_RANGE')
        face = tuple(operand_faces[obj][index])
        triangle = [operands[obj][v] for v in face]
        if not point_on_triangle(point, triangle):
            raise ValueError('BIRTH_FACE_INCIDENCE_MISMATCH')
        facets.append((obj, face))
    explicit = sorted({('EXPLICIT', obj, vertex) for obj, face in facets
                       for vertex in face if tuple(operands[obj][vertex]) == point})
    if explicit:
        # A unique source vertex; coincident independent vertices can separate.
        if len(explicit) != 1:
            raise ValueError('AMBIGUOUS_EXPLICIT_SOURCE_VERTEX')
        return explicit[0]
    candidates = set()
    for obj, face in facets:
        for i in range(3):
            edge = tuple(sorted((face[i], face[(i + 1) % 3])))
            a, b = [operands[obj][v] for v in edge]
            if a == b or cross(subtract(point, a), subtract(b, a)) != (0, 0, 0) or dot(subtract(point, a), subtract(point, b)) > 0:
                continue
            for other, plane in facets:
                if obj == other:
                    continue
                candidate = ('EDGE_PLANE', obj, edge, other, tuple(sorted(plane)))
                try:
                    if replay_vertex_recipe(candidate, operands) == point:
                        candidates.add(candidate)
                except ValueError:
                    continue
    if not candidates:
        raise ValueError('UNRESOLVED_EXACT_VERTEX_RECIPE')
    # This is a base-point witness. It must be replayed/checked at every sample;
    # choosing it does not itself prove persistent support or differentiability.
    return min(candidates)
