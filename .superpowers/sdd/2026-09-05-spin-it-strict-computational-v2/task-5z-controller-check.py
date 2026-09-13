"""Independent fixed-protocol checks: rational tetrahedra, oriented chains and
exact triangle distances on actual binary64 coordinates. No native imports,
no arbitrary-shape topology certification, no seal/paper/physics grant.
"""
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path
import struct
import sys

path = Path(__file__).with_name('task-5y-controller-check.py')
spec = importlib.util.spec_from_file_location('binary_oracle', path)
binary = importlib.util.module_from_spec(spec)
spec.loader.exec_module(binary)
geo = binary.geo


def mesh_digest(mesh):
    h = hashlib.sha256(b'exact-boundary-input-v1\0' + bytes([2]))
    for n in (len(mesh['vertices']), len(mesh['faces'])): h.update(struct.pack('>Q', n))
    for p in mesh['vertices']:
        for v in p: h.update(struct.pack('>d', v))
    for f in mesh['faces']:
        for v in f: h.update(struct.pack('>Q', v))
    return h.hexdigest()


def digest(fields):
    return hashlib.sha256(b''.join(x.encode()+b'\0' for x in fields)).hexdigest()


def selected_mesh_digest(mesh):
    label = b'actual-selected-oriented-primitive-mesh-v1'
    h = hashlib.sha256(struct.pack('>Q',len(label))+label)
    for key in ('vertices','faces'):
        h.update(struct.pack('>QQ',len(mesh[key]),3))
        for row in mesh[key]:
            for value in row: h.update(struct.pack('>d',float(value)))
    return h.hexdigest()


def check_region(r, mesh, material_components, void_components):
    assert r['completed'] is True and r['valid'] is True
    empty = not mesh['faces']
    assert r['empty'] is empty
    assert r['material_components'] == material_components
    assert r['void_components'] == void_components
    if empty:
        assert mesh['vertices'] == [] and r['boundary_input'] is None
        assert r['input'] == '' and not r['parents'] and not r['signs']
        identity = ['explicit-empty']
    else:
        expected = mesh_digest(mesh)
        assert r['input'] == r['boundary_input'] == expected
        identity = ['indexed-boundary', expected]
        assert len(r['parents']) == len(r['depths']) == len(r['signs']) == material_components+void_components
        assert sorted(r['depths']) == ([0, 1] if void_components else [0])
        assert sorted(r['signs']) == ([-1, 1] if void_components else [1])
    assert r['digest'] == digest([r['algorithm'], r['backend'], r['semantics'], *identity])


def check_partition(p, whole, parts):
    assert p['completed'] is True
    assert p['whole'] == whole['digest'] and p['whole_empty'] is whole['empty']
    assert p['parts'] == [r['digest'] for r in parts]
    assert p['parts_empty'] == [r['empty'] for r in parts]
    assert p['complete_partition'] is True and p['equal_union'] is True and p['no_overlap'] is True
    assert p['missing'] is False and p['excess'] is False
    assert p['subsets'] == [True]*len(parts) and len(parts) > 0


def compact_primitive(s, p):
    ids = sorted({i for f in p['faces'] for i in f})
    index = {old: new for new, old in enumerate(ids)}
    return {'vertices': [s['points'][i] for i in ids],
            'faces': [[index[i] for i in f] for f in p['faces']]}


def check_regions(c, r, index):
    assert r['full_thickness'] == 'NOT_ESTABLISHED'
    assert len(r['scientific_hash']) == 64
    assert r['minimum'] == (.13 if index == 1 else .05)
    assert r['nef_limit'] == (1 if index == 2 else 1000000)
    assert r['binary_input_hash'] == c['scientific_hash']
    if c['id'] == 'mixed_rounding_failure':
        assert r['outcome'] == 'NOT_ADMITTED' and r['published'] is False
        assert r['admitted'] is False and r['same_binary'] is False
        return
    assert r['admitted'] is True and r['identities'] is True
    assert r['partition_hash'] == c['source']['partition_hash']
    assert r['source_hash'] == c['source']['source_hash']
    assert r['evaluation_hash'] == c['binary_value_hash'] or index == 2
    if index == 2:
        assert r['outcome'] == 'RESOURCE_LIMIT' and r['published'] is False and r['same_binary'] is False
        assert r['nef_exceeded'] is True
        return
    empty = c['empty_cavity']
    check_region(r['outer'], c['outer'], 1, 0)
    check_region(r['carved'], c['cavity'], 0 if empty else 1, 0)
    check_region(r['material'], c['material'], 1, 0 if empty else 1)
    check_partition(r['material_partition'], r['outer'], [r['material'], r['carved']])
    s = c['source']
    selected = [(i,p) for i,p in enumerate(s['primitives']) if s['state'][p['owner']] == 1]
    assert r['selected_count'] == len(r['selected']) == len(selected)
    for (i,p), actual in zip(selected, r['selected']):
        assert actual['primitive'] == i and actual['owner'] == p['owner']
        assert actual['beta'] == 1 and actual['level'] == p['level'] and actual['key'] == p['key']
        mesh = compact_primitive(s, p)
        assert actual['mesh_hash'] == selected_mesh_digest(mesh)
        check_region(actual['region'], mesh, 1, 0)
        assert actual['subset_completed'] is True and actual['subset'] is True
        assert actual['subset_a'] == actual['region']['digest'] and actual['subset_b'] == r['outer']['digest']
        # Fixed protocol only: outer is convex, so all selected vertices
        # satisfy every outward face halfspace; triangles are linear.
        points = [geo.vec(v) for v in c['outer']['vertices']]
        for face in c['outer']['faces']:
            a,b,d = [points[j] for j in face]
            n = geo.cross(geo.sub(b,a),geo.sub(d,a))
            assert all(geo.dot(n,geo.sub(geo.vec(v),a)) <= 0 for v in mesh['vertices'])
    parts = [p['region'] for p in r['selected']] or [r['carved']]
    check_partition(r['selected_partition'], r['carved'], parts)
    assert r['design_member'] is True
    q = r['clearance']
    assert q['completed'] is True
    assert q['a'] == r['outer']['digest'] and q['b'] == r['carved']['digest']
    assert q['a_input'] == r['outer']['input'] and q['b_input'] == r['carved']['input']
    assert F(q['minimum_exact']) == F(r['minimum'])
    if empty:
        assert q['state'] == 'NOT_APPLICABLE' and q['applicable'] is False
        assert all(q[k] is None for k in ('contact','at_least','squared_distance','witness'))
        passed = True
    else:
        a, b = c['outer'], c['cavity']
        points = [geo.vec(v) for v in a['vertices']+b['vertices']]
        fb = [[i+len(a['vertices']) for i in f] for f in b['faces']]
        distance = geo.boundary_distance(points, a['faces'], fb)
        assert distance == F(q['squared_distance']) and distance > 0
        fa, fc = q['witness']
        witness = geo.triangle_distance([geo.vec(a['vertices'][i]) for i in a['faces'][fa]],
                                       [geo.vec(b['vertices'][i]) for i in b['faces'][fc]])
        assert witness == distance
        passed = distance >= F(r['minimum'])**2
        assert q['applicable'] is True and q['contact'] is False and q['at_least'] is passed
    assert r['published'] is passed and r['same_binary'] is passed
    assert r['outcome'] == ('EXACT_MATERIAL_REGIONS_VALIDATED_PRE_SEAL' if passed else 'EXACT_MATERIAL_REGIONS_REJECTED')
    assert r['nef_exceeded'] is False and 0 < r['nef_peak'] <= r['nef_limit']


def check_file(path):
    rows, error = binary.check_file(path)
    for c in rows:
        assert len(c['regions']) == (3 if c['id'] == 'full_affine' else 1)
        for i, r in enumerate(c['regions']): check_regions(c,r,i)
    return rows, error


if __name__ == '__main__':
    _, error = check_file(sys.argv[1])
    print(f'PASS: 4 actual binary inputs / 6 native region-gate calls; independent mesh hashes, tetrahedra, selected chains, convex subset, rational distance and witness; max moment error {error:.17g}; no seal or original-paper grant')
