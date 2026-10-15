"""Finite exact supporting-facet connection witnesses, not neighborhood proofs.

Inputs must separately pass the production embedded-boundary audit. This
module only explains coplanar triangulation and collinear storage differences.
"""
from collections import Counter, defaultdict
from fractions import Fraction
import hashlib
from itertools import combinations
import json


def _digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True,
        separators=(',', ':')).encode()).hexdigest()


def _points(rows):
    result = [tuple(Fraction(x) for x in row) for row in rows]
    if any(len(row) != 3 for row in result):
        raise ValueError('INVALID_EXACT_POINT')
    return result


def _cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


def _plane(points, face):
    a, b, c = [points[int(i)] for i in face]
    normal = _cross(tuple(b[i]-a[i] for i in range(3)), tuple(c[i]-a[i] for i in range(3)))
    if not any(normal):
        raise ValueError('DEGENERATE_FACE')
    values = (*normal, -sum(normal[i]*a[i] for i in range(3)))
    first = next(value for value in values if value)
    return tuple(value/first for value in values)


def _mesh_identity(points, faces):
    return dict(vertices=[[str(x) for x in row] for row in points],
                faces=[list(map(int, face)) for face in faces])


def facet_connection_witness(vertices, faces, input_meshes):
    points = _points(vertices)
    planes = defaultdict(list)
    explicit = defaultdict(list)
    input_identity = []
    for obj, (input_vertices, input_faces) in enumerate(input_meshes):
        source = _points(input_vertices)
        input_identity.append(_mesh_identity(source, input_faces))
        local = defaultdict(list)
        for index, face in enumerate(input_faces):
            local[_plane(source, face)].append(index)
        for plane, indices in local.items():
            planes[plane].append((obj, tuple(indices)))
        for index, point in enumerate(source):
            explicit[point].append((obj, index))
    grouped = defaultdict(list)
    facet_normals = {}
    edge_counts = Counter()
    support = defaultdict(set)
    for face in faces:
        face = tuple(map(int, face))
        if len(face) != 3 or len(set(face)) != 3 or any(i < 0 or i >= len(points) for i in face):
            raise ValueError('INVALID_FACE')
        plane = _plane(points, face)
        if plane not in planes:
            raise ValueError('UNEXPLAINED_MATERIAL_PLANE')
        facet = tuple(sorted(planes[plane]))
        facet_normals[facet] = plane[:3]
        grouped[facet].append(face)
        for vertex in face:
            support[vertex].add(facet)
        for a, b in zip(face, face[1:]+face[:1]):
            edge_counts[(a, b)] += 1
    for (a, b), count in edge_counts.items():
        if count != 1 or edge_counts[(b, a)] != 1:
            raise ValueError('OPEN_OR_INCONSISTENT_BOUNDARY')
    keys = {}
    used_keys = {}
    def vertex_key(index):
        if index not in keys:
            if points[index] in explicit:
                key = ('INPUT_VERTEX', tuple(sorted(explicit[points[index]])))
            else:
                facets = tuple(sorted(support[index]))
                if len(facets) < 3:
                    raise ValueError('NONUNIQUE_CORNER_SUPPORT')
                normals = [facet_normals[facet] for facet in facets]
                if not any(sum(a[i]*_cross(b,c)[i] for i in range(3))
                           for a,b,c in combinations(normals,3)):
                    raise ValueError('NONUNIQUE_CORNER_SUPPORT')
                key = ('SUPPORT_INTERSECTION', facets)
            if key in used_keys and used_keys[key] != points[index]:
                raise ValueError('NONUNIQUE_CORNER_SUPPORT')
            used_keys[key] = points[index]
            keys[index] = json.dumps(key, separators=(',', ':'))
        return keys[index]
    connections = []
    for facet, triangles in sorted(grouped.items()):
        edges = Counter((a,b) for face in triangles for a,b in zip(face,face[1:]+face[:1]))
        boundary = {(a,b) for a,b in edges if not edges[(b,a)]}
        owners = {(a,b):index for index,face in enumerate(triangles)
                  for a,b in zip(face,face[1:]+face[:1])}
        successors = {}
        for edge in boundary:
            a,b = edge
            index = owners[edge]
            visited = set()
            while True:
                if index in visited:
                    raise ValueError('CLOSED_FAN_WITH_BOUNDARY_EDGE')
                visited.add(index)
                face = triangles[index]
                position = face.index(b)
                if face[(position-1)%3] != a:
                    raise ValueError('INCONSISTENT_FACET_FAN')
                c = face[(position+1)%3]
                if (b,c) in boundary:
                    successors[edge] = (b,c)
                    break
                if (c,b) not in owners:
                    raise ValueError('OPEN_FACET_BOUNDARY')
                a,index = c,owners[(c,b)]
        if set(successors.values()) != boundary:
            raise ValueError('AMBIGUOUS_FACET_BOUNDARY')
        cycles = []
        while boundary:
            start = min(boundary)
            cycle, current = [], start
            while True:
                cycle.append(current[0])
                if current not in boundary:
                    raise ValueError('INVALID_FACET_CYCLE')
                boundary.remove(current)
                current = successors[current]
                if current == start:
                    break
            changed = True
            while changed and len(cycle) > 3:
                changed = False
                for i in range(len(cycle)):
                    a,b,c = [points[cycle[j % len(cycle)]] for j in (i-1,i,i+1)]
                    u,v = tuple(b[k]-a[k] for k in range(3)), tuple(c[k]-b[k] for k in range(3))
                    if not any(_cross(u,v)) and sum(u[k]*v[k] for k in range(3)) > 0:
                        cycle.pop(i); changed = True; break
            labels = tuple(vertex_key(i) for i in cycle)
            cycles.append(min(labels[i:]+labels[:i] for i in range(len(labels))))
        connections.append(dict(facet=facet, oriented_cycles=sorted(cycles)))
    return dict(schema='exact-supporting-facet-connections-v1', connections=connections,
        material_boundary_sha256=_digest(_mesh_identity(points, faces)),
        input_boundary_sha256=_digest(input_identity), neighborhood_certification=False,
        scope='Finite oriented facet connections; separately audited inputs required')
