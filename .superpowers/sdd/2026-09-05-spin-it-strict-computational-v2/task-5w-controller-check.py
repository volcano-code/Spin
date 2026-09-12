"""Independent raw geometry/BBW/ownership check of Task5w synthetic exports.

The native transaction is not imported. Exact signed-tetrahedral integration
uses the separately implemented controller oracle from Task5u.
"""
from fractions import Fraction as F
import importlib.util
import json
import math
from pathlib import Path
import sys

spec = importlib.util.spec_from_file_location(
    "independent_geometry", Path(__file__).with_name("task-5u-controller-check.py"))
geo = importlib.util.module_from_spec(spec)
spec.loader.exec_module(geo)


def ancestor(a, b):
    return (a[0] <= b[0] and all(x == y >> (b[0]-a[0]) for x, y in zip(a[1:], b[1:])))


def key(leaf):
    return (leaf['level'], *leaf['key'])


def close(a, b, tolerance=1e-10):
    assert len(a) == len(b) and len(a) > 0
    assert all(math.isfinite(float(x)) for x in (*a, *b))
    scale = max(max(abs(float(x)) for x in b), 1e-30)
    error = max(abs(float(x)-float(y)) for x, y in zip(a, b)) / scale
    assert error <= tolerance, (error, tolerance)
    return error


def blend_point(p, source):
    """Independent closed-cell lookup and 8-node trilinear field interpolation."""
    nodes, cells, field = source['grid_nodes'], source['grid_cells'], source['weights']
    hits = []
    for cell in cells:
        lo, hi = nodes[cell[0]], nodes[cell[7]]
        if all(lo[k] <= p[k] <= hi[k] for k in range(3)):
            hits.append((tuple(lo), cell))
    assert hits
    _, cell = min(hits)
    lo, hi = nodes[cell[0]], nodes[cell[7]]
    t = [(p[k]-lo[k])/(hi[k]-lo[k]) for k in range(3)]
    weights = [0.] * source['handles']
    for bits, node in enumerate(cell):
        coefficient = math.prod(t[k] if bits & (1 << k) else 1-t[k] for k in range(3))
        for h in range(len(weights)):
            weights[h] += coefficient*field[node][h]
    assert abs(sum(weights)-1) < 1e-8
    params = source['state'][source['handle_offset']:]
    return [sum(weights[h] * (center[k] + params[4*h]*(p[k]-center[k]) + params[4*h+k+1])
                for h, center in enumerate(source['centers'])) for k in range(3)]


def check_source(s):
    beta = s['state'][:s['betas']]
    assert len(beta) == len(s['leaves']) > 0
    assert all(0 <= v <= 1 and math.isfinite(v) for v in beta)
    assert s['handles'] == 26
    assert s['handle_offset'] == s['betas'] + int(s['angle'])
    assert len(s['state']) == s['handle_offset'] + 4*s['handles']
    assert len(s['points']) == len(s['rest']) > 8
    for rest, point in zip(s['rest'], s['points']):
        close(point, blend_point(rest, s), 2e-13)
    points = [geo.vec(p) for p in s['points']]
    result = geo.moments(points, s['outer_faces'])
    primitive_moments = []
    jac = [[F(0) for _ in beta] for _ in range(10)]
    rest_per_owner = [[F(0)]*10 for _ in beta]
    rest = [geo.vec(p) for p in s['rest']]
    keys = [key(leaf) for leaf in s['leaves']]
    assert len(keys) == len(set(keys))
    assert not any(ancestor(a, b) for i, a in enumerate(keys) for j, b in enumerate(keys) if i != j)
    for primitive in s['primitives']:
        owner = primitive['owner']
        assert 0 <= owner < len(beta)
        assert ancestor(keys[owner], key(primitive))
        raw = geo.moments(points, primitive['faces'])
        assert raw[0] > 0
        primitive_moments.append(raw)
        raw_rest = geo.moments(rest, primitive['faces'])
        for j in range(10):
            result[j] -= F(beta[owner])*raw[j]
            jac[j][owner] -= F(s['density'])*raw[j]
            rest_per_owner[owner][j] += raw_rest[j]
    result = [x*F(s['density']) for x in result]
    error = close(s['moments'], result)
    for j in range(10):
        close(s['beta_jacobian'][j], jac[j])
    for i, leaf in enumerate(s['leaves']):
        close(leaf['moments'], rest_per_owner[i])
    assert result[0] > 0
    return error


def verify(c):
    assert c['scope'] == 'SYNTHETIC_REGRESSION'
    assert c['sealed'] is False and c['converged'] is False
    assert c['checkpoint'] is True
    old, new = c['before'], c['after']
    error = max(check_source(old), check_source(new))
    assert old['bbw_hash'] == new['bbw_hash'] and old['source_hash'] == new['source_hash']
    assert old['centers'] == new['centers'] and old['weights'] == new['weights']
    assert old['density'] == new['density'] == 1.17
    assert old['state'][old['betas']:] == new['state'][new['betas']:]
    for i, leaf in enumerate(new['leaves']):
        covering = [j for j, source in enumerate(old['leaves']) if ancestor(key(source), key(leaf))]
        if covering:
            assert len(covering) == 1
            assert new['state'][i] == old['state'][covering[0]]
        else:
            contained = [j for j, source in enumerate(old['leaves']) if ancestor(key(leaf), key(source))]
            assert contained and all(old['state'][j] == new['state'][i] for j in contained)
    if c['id'].startswith('split'):
        assert old['betas'] == 1 and new['betas'] == 8
        assert len(old['primitives']) == 1 and len(new['primitives']) == 8
        assert len(new['points']) > len(old['points'])
        assert old['partition_hash'] != new['partition_hash']
        assert old['upstream_hash'] != new['upstream_hash']
    elif c['id'].startswith('merge'):
        assert old['betas'] == 8 and new['betas'] == 1
    elif c['id'] == 'unequal_no_merge':
        assert old['betas'] == new['betas'] == 8
        assert old['partition_hash'] == new['partition_hash']
    else:
        raise AssertionError('unregistered controller case')
    if c['id'] in ('split_rest', 'split_affine', 'merge_equal'):
        close(new['moments'], old['moments'])
    # General non-affine moment changes are recorded, not forced to zero.
    if 'delta' in c:
        expected = [y-x for x, y in zip(old['moments'], new['moments'])]
        assert all(abs(x-y) <= 1e-15 for x, y in zip(c['delta'], expected))
    if c['id'] == 'split_static':
        consumer = c['consumer']
        assert consumer['initial'] == new['state']
        assert consumer['trace0']['state'] == new['state']
        assert consumer['trace0']['history_before'] == 0
        assert consumer['trace0']['history_reset_reason'] == 'NEW_SOLVE_OR_PARTITION'
        assert consumer['guard_calls'] >= 2 and consumer['final_rechecks'] == 1
        assert consumer['final_status'] == 1
        assert consumer['status'] in ('MAXIMUM_ITERATIONS', 'RESTORATION_STALLED', 'CONVERGED', 'LINE_SEARCH_FAILED', 'ITERATION_LIMIT')
    return error


if __name__ == '__main__':
    geo.self_test()
    rows = [json.loads(line) for line in Path(sys.argv[1]).read_text().splitlines()]
    assert len(rows) == 6 and len({c['id'] for c in rows}) == 6
    assert {c['id'] for c in rows} == {'split_rest', 'split_affine', 'split_nonaffine',
                                     'split_static', 'merge_equal', 'unequal_no_merge'}
    errors = [verify(c) for c in rows]
    print(f'PASS: {len(rows)} actual rebind cases; independent BBW interpolation, ownership, signed tetrahedral moments and beta derivatives; max error {max(errors):.17g}')
