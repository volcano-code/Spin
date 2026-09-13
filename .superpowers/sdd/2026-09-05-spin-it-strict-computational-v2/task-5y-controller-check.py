"""Independent checking of actual-deformed binary exports, never native imports.

Uses exact rational signed tetrahedra, exported full BBW interpolation and
oriented triangle cancellation. Limited to this fixed synthetic protocol;
does not grant final material partition, full thickness, a seal or physics.
"""
from collections import Counter
from fractions import Fraction as F
import importlib.util
import json
import math
from pathlib import Path
import sys
import numpy as np

spec = importlib.util.spec_from_file_location('source_oracle', Path(__file__).with_name('task-5w-controller-check.py'))
source = importlib.util.module_from_spec(spec)
spec.loader.exec_module(source)
geo = source.geo


def oriented_triangle(points, face):
    a = tuple(tuple(points[i]) for i in face)
    return min(a, a[1:]+a[:1], a[2:]+a[:2])


def boundary_chain(points, faces):
    return Counter(oriented_triangle(points, f) for f in faces)


def selected_boundary(s):
    chain = Counter()
    for p in s['primitives']:
        if s['state'][p['owner']] == 1:
            chain.update(oriented_triangle(s['points'], f) for f in p['faces'])
    result = Counter()
    for tri, count in chain.items():
        reverse = oriented_triangle(tri, (0, 2, 1))
        net = count-chain[reverse]
        assert net in (-1, 0, 1)
        if net > 0: result[tri] = net
    return result


def integrate(mesh, rho):
    v, f = mesh['vertices'], mesh['faces']
    if not f:
        assert not v
        return np.zeros(10)
    assert sorted(set(i for t in f for i in t)) == list(range(len(v)))
    exact = geo.moments([geo.vec(p) for p in v], f)
    return np.array([float(F(rho)*x) for x in exact])


def transform(values, origin, rotation):
    m = values[0]
    l = values[1:4]
    q = np.array([[values[4], values[9], values[8]],
                  [values[9], values[5], values[7]],
                  [values[8], values[7], values[6]]])
    o, r = np.array(origin), np.array(rotation)
    assert np.max(abs(r @ r.T-np.eye(3))) < 1e-12 and np.linalg.det(r) > 0
    shifted = q - np.outer(o, l)-np.outer(l, o)+m*np.outer(o, o)
    q = r @ shifted @ r.T
    l = r @ (l-m*o)
    return np.array([m, *l, *q.diagonal(), q[1,2], q[0,2], q[0,1]])


def constraints(s, angle, static):
    if static: return np.array(s[1:3])
    c, t = math.cos(angle), math.sin(angle)
    return np.array([s[1],s[2],s[8],s[7],c*t*(s[4]-s[5])+(c*c-t*t)*s[9]])


def verify(c):
    assert c['scope'] == 'KNOWN_SOLUTION_VALIDATION'
    assert c['sealed'] is False and c['paper_reproduction'] is False
    assert c['id'] in {'empty_affine','full_affine','tie_affine','mixed_rounding_failure'}
    assert c['adaptive_outcome'] == 'RELAXED_ADAPTIVE_COMPLETE'
    s, relaxed = c['source'], np.array(c['relaxed_state'])
    error = source.check_source(s)
    beta = s['betas']
    rounded = np.array(s['state'])
    assert len(relaxed) == len(rounded)
    assert np.array_equal(rounded[:beta], (relaxed[:beta] >= .5).astype(float))
    assert np.array_equal(rounded[beta:], relaxed[beta:])
    assert np.array_equal(np.array(c['delta']), rounded-relaxed)
    assert s['density'] == 1.17 and s['handles'] == 26 and beta == 8
    scale, translation = 1.02, np.array([.011,-.007,.019])
    for h, center in enumerate(s['centers']):
        offset = s['handle_offset']+4*h
        assert rounded[offset] == scale
        assert np.max(abs(rounded[offset+1:offset+4]-((scale-1)*np.array(center)+translation))) < 1e-15
    assert np.max(abs(np.array(s['points'])-(scale*np.array(s['rest'])+translation))) < 1e-12
    expected_beta = {'empty_affine':.25,'full_affine':.75,'tie_affine':.5}.get(c['id'])
    if expected_beta is not None: assert np.all(relaxed[:beta] == expected_beta)
    else: assert set(relaxed[:beta]) == {.4,.6}
    outer, cavity, material = c['outer'], c['cavity'], c['material']
    n = len(outer['vertices'])
    assert outer['vertices'] == s['points'][:n]
    assert outer['faces'] == s['outer_faces']
    expected_cavity = selected_boundary(s)
    assert boundary_chain(cavity['vertices'], cavity['faces']) == expected_cavity
    assert c['empty_cavity'] == (not cavity['faces'])
    expected_material = boundary_chain(outer['vertices'],outer['faces'])
    expected_material.update(oriented_triangle(t,(0,2,1)) for t in expected_cavity.elements())
    assert boundary_chain(material['vertices'],material['faces']) == expected_material
    rho = s['density']
    world = integrate(material,rho)
    subtract = integrate(outer,rho)-integrate(cavity,rho)
    assert world[0] > 0
    for actual in (subtract,c['world_integrated'],s['moments'],c['binary_functional_moments']):
        err = float(np.max(abs(np.array(actual)-world))/rho)
        assert math.isfinite(err) and err <= 1e-10
        error = max(error,err)
    spin = transform(world,c['origin'],c['rotation'])
    assert np.max(abs(spin-np.array(c['spin_integrated'])))/rho <= 1e-10
    phi = rounded[beta] if s['angle'] else 0.
    raw = constraints(spin,phi,not s['angle'])
    assert np.max(abs(raw-np.array(c['constraints'])))/rho <= 1e-10
    assert np.max(abs(raw/rho-np.array(c['normalized_constraints']))) <= 1e-10
    residual = float(np.max(abs(raw))/rho)  # fixed input root side=1
    assert abs(residual-c['constraint_residual']) < 1e-10
    assert c['constraint_tolerance'] == 1e-6
    q = np.array([[world[4],world[9],world[8]], [world[9],world[5],world[7]], [world[8],world[7],world[6]]])
    central_q = q-np.outer(world[1:4],world[1:4])/world[0]
    inertia = np.trace(central_q)*np.eye(3)-central_q
    assert np.max(abs(inertia-np.array(c['central_inertia'])))/rho <= 1e-10
    assert np.linalg.eigvalsh(inertia).min() > 0
    assert np.max(abs(np.linalg.eigvalsh(inertia)-np.array(c['inertia_eigenvalues'])))/rho < 1e-10
    assert c['inertia_psd'] is True
    assert c['geometry_passed'] is True and c['material_boundary_valid'] is True
    failure = c['id'] == 'mixed_rounding_failure'
    assert (residual > 1e-6) == failure
    assert c['validated'] is (not failure)
    assert c['outcome'] == ('BINARY_VALIDATION_FAILED' if failure else 'BINARY_MATERIAL_VALIDATED_PRE_SEAL')
    assert len(c['scientific_hash']) == 64
    return error


def analytic_self_test():
    lo, hi, rho = np.array([-.2,.3,.4]), np.array([.7,.9,1.2]), 1.17
    vertices = [[hi[k] if i & (1<<k) else lo[k] for k in range(3)] for i in range(8)]
    faces = [[0,2,1],[1,2,3],[4,5,6],[5,7,6],[0,1,4],[1,5,4],
             [2,6,3],[3,6,7],[0,4,2],[2,4,6],[1,3,5],[3,7,5]]
    raw = integrate({'vertices':vertices,'faces':faces},rho)
    mass = rho*np.prod(hi-lo)
    mean = (lo+hi)/2
    q = mass*(lo*lo+lo*hi+hi*hi)/3
    expected = np.array([mass,*(mass*mean),*q,mass*mean[1]*mean[2],mass*mean[0]*mean[2],mass*mean[0]*mean[1]])
    assert np.max(abs(raw-expected)) < 1e-14
    rotation = np.array([[0.,-1.,0.],[1.,0.,0.],[0.,0.,1.]])
    origin = np.array([.17,-.29,.08])
    moved = [(rotation@(np.array(v)-origin)).tolist() for v in vertices]
    reference = integrate({'vertices':moved,'faces':faces},rho)
    assert np.max(abs(transform(raw,origin,rotation)-reference)) < 1e-14


def check_file(path):
    geo.self_test()
    analytic_self_test()
    rows = [json.loads(line) for line in Path(path).read_text().splitlines()]
    assert len(rows) == 4 and len({r['id'] for r in rows}) == 4
    result = max(verify(r) for r in rows)
    return rows, result


if __name__ == '__main__':
    _, error = check_file(sys.argv[1])
    print(f'PASS: 4 fixed actual-deformed binary cases, independent BBW/triangle cancellation/rational tetrahedra/constraints/inertia; max error {error:.17g}; no seal or paper grant')
