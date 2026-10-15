"""Exact local witness for Boolean triangulation changes, not a derivative gate.

At each parameter point transplant the base triangulation by unique source
feature labels. Prove that it and the actual triangulation have the same
oriented boundary on each parent plane, and consistently oriented nonzero
triangles. Compact planar 2-chains with equal boundary have equal winding
multiplicity. This certifies retriangulation at THAT point, not equality of
entities at different parameters, nor differentiability across events.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

import gmpy2


def _coordinates(rows):
    if any(not isinstance(x, str) for row in rows for x in row):
        raise ValueError('AUTHORITY_COORDINATES_MUST_BE_RATIONAL_TEXT')
    return [tuple(gmpy2.mpq(x) for x in row) for row in rows]


def _sub(a, b):
    return tuple(x-y for x,y in zip(a,b))


def _cross(a,b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


def _dot(a,b):
    return sum(x*y for x,y in zip(a,b))


def _patches(evaluation, faces=None, births=None):
    vertices = _coordinates(evaluation['material_vertices_target_authority'])
    design = _coordinates(evaluation['design_vertices_target_authority'])
    source = _coordinates(evaluation['source_vertices_target_authority'])
    labels = evaluation['material_vertex_provenance']
    if len(labels) != len(vertices) or len(set(labels)) != len(labels):
        raise ValueError('AMBIGUOUS_VERTEX_PROVENANCE')
    faces = evaluation['material_faces'] if faces is None else faces
    births = evaluation['material_birth_faces'] if births is None else births
    if not faces or len(faces) != len(births):
        raise ValueError('MISSING_FACE_OR_PARENT')
    patches = {}
    for face, birth in zip(faces,births):
        if len(face) != 3 or len(set(face)) != 3 or birth < 0:
            raise ValueError('INVALID_TRIANGLE_OR_PARENT')
        if birth < len(evaluation['design_faces']):
            parent = [design[j] for j in evaluation['design_faces'][birth]]
        else:
            parent = [source[j] for j in evaluation['source_faces'][birth-len(evaluation['design_faces'])]]
        normal = _cross(_sub(parent[1],parent[0]), _sub(parent[2],parent[0]))
        if not any(normal):
            raise ValueError('DEGENERATE_PARENT')
        pts = [vertices[j] for j in face]
        if any(_dot(_sub(v,parent[0]),normal) != 0 for v in pts):
            raise ValueError('EXACT_NONCOPLANAR_TRIANGLE')
        orientation = gmpy2.sign(_dot(_cross(_sub(pts[1],pts[0]),_sub(pts[2],pts[0])),normal))
        if orientation == 0:
            raise ValueError('DEGENERATE_TRIANGLE')
        patch = patches.setdefault(birth, dict(orientation=orientation, boundary=Counter()))
        if patch['orientation'] != orientation:
            raise ValueError('INCONSISTENT_PATCH_ORIENTATION')
        for a,b in zip(face,face[1:]+face[:1]):
            la,lb = labels[a],labels[b]
            edge = tuple(sorted((la,lb)))
            patch['boundary'][edge] += 1 if la < lb else -1
    for patch in patches.values():
        patch['boundary'] = {k:v for k,v in patch['boundary'].items() if v}
    return patches


def compare_retriangulation(base, current):
    result = dict(schema='exact-parent-plane-retriangulation-witness-v1', proved=False,
        scope='BASE_CONNECTIVITY_TRANSPLANTED_AT_CURRENT_PARAMETER_ONLY',
        derivative_gate_granted=False, geometric_tolerance_used=False)
    try:
        for key in ('source_faces','design_faces'):
            if base[key] != current[key]:
                raise ValueError('PARENT_CONNECTIVITY_CHANGED')
        base_labels = base['material_vertex_provenance']
        labels = current['material_vertex_provenance']
        if Counter(base_labels) != Counter(labels):
            raise ValueError('SOURCE_FEATURE_SET_CHANGED')
        original = _patches(base)
        actual = _patches(current)
        index = {label:i for i,label in enumerate(labels)}
        transplanted_faces = [[index[base_labels[j]] for j in face] for face in base['material_faces']]
        transplanted = _patches(current, transplanted_faces, base['material_birth_faces'])
        if original != actual or transplanted != actual:
            raise ValueError('ORIENTED_PARENT_PATCH_BOUNDARY_CHANGED')
        result.update(proved=True, parent_patch_count=len(actual),
            triangle_count=len(current['material_faces']), vertex_count=len(labels),
            theorem='compact coplanar oriented 2-chains with equal boundary have equal winding multiplicity')
    except (ValueError, KeyError, IndexError, TypeError) as exc:
        result['reason'] = str(exc)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--points-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    paths = [args.points_root/f'split-production-point-04-{i:03}/output/point.json' for i in range(81)]
    base = json.loads(paths[0].read_text())['evaluation']
    records = []
    for i,path in enumerate(paths):
        point = json.loads(path.read_text())
        if point['index'] != i:
            raise ValueError('POINT_INDEX_MISMATCH')
        record = compare_retriangulation(base,point['evaluation'])
        record.update(index=i, point_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
        records.append(record)
    report = dict(schema='81-point-retriangulation-audit-v1',
        witness_implementation_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        records=records, all_local_witnesses_proved=all(r['proved'] for r in records),
        production_derivative_gate_granted=False)
    with args.output.open('x') as out:
        json.dump(report,out,indent=2,allow_nan=False)
    print(json.dumps({k:v for k,v in report.items() if k!='records'}))


if __name__ == '__main__':
    main()
