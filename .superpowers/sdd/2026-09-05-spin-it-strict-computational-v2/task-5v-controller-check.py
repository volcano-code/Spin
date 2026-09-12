"""Independent arithmetic/linkage checks of actual guarded source exports.

No native admission/solver imports. This does not validate arbitrary meshes or
grant a sealed, feasible, original-paper or physical result. Null trial residuals
are explicit unavailable observations, never numerical success.
"""
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import struct
import sys
import numpy as np


def load(name):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(name+'.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


geometry = load('task-5u-controller-check')
arithmetic = load('task-5t-controller-check')
EXPECTED = {**{f'mode{i}': i for i in range(5)},
            'initial_reject': 2, 'guard_cap': 2, 'source_cap': 2,
            'restoration': 0, 'exhaustion': 3}
ABSENT = 2**64-1
ORDER = [(s, ABSENT, ABSENT) for s in (0, 1, 2)]
ORDER += [(s, i, ABSENT) for i in range(8) for s in (3, 4, 5)]
ORDER += [(6, i, j) for i in range(8) for j in range(i+1, 8)]


def sha(value):
    assert isinstance(value, str) and len(value) == 64
    assert set(value) <= set('0123456789abcdef')


def finite(value):
    if isinstance(value, dict):
        for v in value.values(): finite(v)
    elif isinstance(value, list):
        for v in value: finite(v)
    elif isinstance(value, float):
        assert math.isfinite(value)


def state_hash(value):
    return hashlib.sha256(struct.pack('>QQ', len(value), 1)
                          + struct.pack('>'+str(len(value))+'d', *value)).hexdigest()


def verify(c, quiet=False):
    finite(c)
    assert c['id'] in EXPECTED and c['mode'] == EXPECTED[c['id']]
    assert c['scope'] == 'SYNTHETIC_REGRESSION' and c['sealed'] is False
    assert c['feasibility'] == c['legacy_geometry'] == 0
    sha(c['input_hash']); sha(c['result_hash'])
    n = 112 + int(c['mode'] != 4)
    initial = arithmetic.array(c['initial'], (n,))
    calls, guards, stages, trace = (c[k] for k in ('calls', 'guards', 'stages', 'trace'))
    assert calls and len(calls) == len(guards)
    assert calls[0]['state'] == c['initial']
    assert c['callbacks'] == sum(r['purpose'] == 0 for r in calls)
    assert c['source_evaluations'] == sum(r['source_invoked'] for r in calls)
    assert c['checker_calls'] == sum(g['invoked'] for g in guards)
    assert c['checker_calls'] <= c['maximum_guard_calls']
    assert c['geometry_values_complete'] is True
    assert c['geometry_values'] == sum(g['values_invoked'] for g in guards)
    assert c['final_rechecks'] == sum(g['invoked'] and r['purpose'] == 1 for r,g in zip(calls, guards))
    assert [r['purpose'] for r in calls] in ([0]*len(calls), [0]*(len(calls)-1)+[1])
    cursor, latched = 0, False
    passing_states = set()
    for i, (r, g) in enumerate(zip(calls, guards)):
        assert r['sequence'] == i and g['callback'] == i
        arithmetic.array(r['state'], (n,))
        assert r['state_hash'] == state_hash(r['state'])
        assert r['status'] in (1, 2, 3)
        assert g['first_stage'] == cursor
        local = stages[cursor:cursor+g['stage_count']]
        assert len(local) == g['stage_count']
        cursor += len(local)
        assert g['outcome'] in (0, 1, 2, 3)
        assert not g['values_invoked'] or g['invoked']
        if latched:
            assert not g['invoked'] and not r['source_invoked'] and not g['values_invoked']
            assert r['status'] == 3 and g['outcome'] == 2
        if g['completed']:
            assert g['invoked'] and local
            sha(g['input_hash']); sha(g['result_hash'])
            assert [(s['stage'], s['first'], s['second']) for s in local] == ORDER[:len(local)]
            assert g['values_invoked'] == any(s['stage'] == 1 for s in local)
            assert all(s['outcome'] == 1 or (s['stage'] == 6 and s['outcome'] == 5) for s in local[:-1])
            if g['outcome'] == 0:
                assert len(local) == len(ORDER)
                assert local[-1]['outcome'] in (1, 5)
                sha(g['value_hash'])
            else:
                assert local[-1]['outcome'] == {1: 2, 2: 3, 3: 4}[g['outcome']]
                assert local[-1]['reason']
            source_stages = [s for s in local if s['stage'] == 1]
            if source_stages and source_stages[0]['outcome'] == 1:
                assert source_stages[0]['evidence_hash'] == g['value_hash']
        else:
            assert not local and g['outcome'] != 0
        if r['source_invoked']:
            assert g['completed'] and g['outcome'] == 0
        if r['status'] == 1:
            assert r['source_invoked'] and not r['reason']
            sha(r['output_hash'])
            passing_states.add(tuple(r['state']))
        else:
            assert r['reason'] and not r['output_hash']
            sha(r['reason_hash'])
        if r['status'] == 3 or g['outcome'] == 2:
            latched = True
    assert cursor == len(stages)
    assert (c['source_outcome'] == 1) == latched
    assert c['numerical_available'] == bool(trace)
    trial_cursor, accepted = 1, 0
    geometry_rejected, shrunk = 0, 0
    for i, t in enumerate(trace):
        assert t['iteration'] == i and tuple(t['state']) in passing_states
        x = arithmetic.array(t['state'], (n,))
        residual = arithmetic.kkt(t, 8)
        assert (x[:8] >= 0).all() and (x[:8] <= 1).all()
        assert t['final_state'] == (i == len(trace)-1)
        membership = t['active_membership']
        assert len(membership) == n and set(membership) <= {0, 1, 2}
        assert membership[8:] == [0]*(n-8)
        assert t['membership_hash'] == hashlib.sha256(b'spin2014-active-membership-v1'+bytes(membership)).hexdigest()
        assert 0 <= t['history_after'] <= c['controls']['lbfgs_memory']
        assert not (t['curvature_pair_accepted'] and t['curvature_pair_rejected'])
        merit = t['objective'] + c['controls']['merit_penalty']*float(np.abs(t['constraints']).sum())
        assert math.isclose(t['merit'], merit, rel_tol=1e-12, abs_tol=1e-12)
        all_trials = t['objective_trials']+t['restoration']['trials']
        assert len(all_trials) == t['line_search_evaluations']
        accepted_calls = []
        for group in (t['objective_trials'], t['restoration']['trials']):
            previous = None
            for trial in group:
                assert trial_cursor < len(calls)
                r, g = calls[trial_cursor], guards[trial_cursor]
                trial_cursor += 1
                assert r['purpose'] == 0 and trial['valid'] == (r['status'] == 1)
                assert 0 < trial['step'] <= 1 and trial['reason']
                if trial['valid']:
                    assert trial['constraint_l1'] is not None and trial['constraint_l1'] >= 0
                else:
                    assert trial['constraint_l1'] is None
                assert not trial['accepted'] or trial['valid']
                if g['outcome'] == 1: geometry_rejected += 1
                delta = arithmetic.array(r['state'])-x
                if group is t['restoration']['trials']:
                    assert np.allclose(delta, trial['step']*arithmetic.array(t['restoration']['direction']), atol=2e-13, rtol=2e-12)
                if previous and trial['step'] < previous[0]['step']:
                    ratio = trial['step']/previous[0]['step']
                    assert ratio == c['controls']['line_search_reduction']
                    # Unbounded handle/angle coordinates must follow the same direction.
                    assert np.allclose(delta[8:], ratio*previous[1][8:], atol=2e-13, rtol=2e-12)
                    if previous[2]['outcome'] == 1: shrunk += 1
                previous = (trial, delta, g)
                if trial['accepted']:
                    accepted_calls.append(r)
                    assert i+1 < len(trace) and trace[i+1]['state'] == r['state']
                    assert trial['step'] == t['accepted_step']
                    assert math.isclose(trial['constraint_l1'], sum(abs(v) for v in trace[i+1]['constraints']), rel_tol=1e-10, abs_tol=1e-12)
        assert len(accepted_calls) == int(t['accepted_step'] > 0)
        accepted += len(accepted_calls)
        if not accepted_calls:
            assert not t['curvature_pair_accepted'] and not t['curvature_pair_rejected']
        if t['restoration']['accepted']:
            assert accepted_calls and t['restoration']['trials'][-1]['accepted']
            old_l1 = sum(abs(v) for v in t['constraints'])
            predicted = t['restoration']['predicted_residual']
            new_l1 = t['restoration']['trials'][-1]['constraint_l1']
            assert old_l1 > predicted and new_l1 < old_l1
            assert new_l1 <= old_l1-c['controls']['armijo']*t['accepted_step']*(old_l1-predicted)+1e-12
    if trace:
        assert c['numerical_status'] in ('CONVERGED', 'MAXIMUM_ITERATIONS', 'LINE_SEARCH_EXHAUSTED', 'RESTORATION_STALLED', 'NUMERICAL_UNKNOWN')
        assert trace[0]['state'] == c['initial'] and trace[-1]['state'] == c['final']
        assert trial_cursor == len(calls)-1 and calls[-1]['purpose'] == 1
        assert calls[-1]['state'] == c['final']
        if c['numerical_status'] == 'CONVERGED':
            assert residual['equality_primal'] <= c['controls']['equality_tolerance']
            assert max(residual[k] for k in ('bound_primal','stationarity','dual_sign','complementarity')) <= c['controls']['projected_gradient_tolerance']
    else:
        assert len(calls) == 1 and c['final'] == [] and c['final_rechecks'] == 0
        assert c['source_outcome'] == 2 and c['final_status'] == 0
        assert c['numerical_status'] == 'INVALID_PROBLEM'  # Unavailable default, not a termination claim.
    if c['final_status'] == 1:
        assert c['final_evaluation_available'] and 'material' in c and c['final_rechecks'] == 1
        assert calls[-1]['status'] == 1 and guards[-1]['completed'] and guards[-1]['outcome'] == 0
        m = c['material']; points = [geometry.vec(p) for p in m['points']]
        total = geometry.moments(points, m['outer_faces'])
        assert total[0] > 0 and len(m['primitives']) == 8 and m['density'] > 0
        local = stages[guards[-1]['first_stage']:]
        for i, p in enumerate(m['primitives']):
            assert p['owner'] == i
            part = geometry.moments(points, p['faces']); assert part[0] > 0
            d = geometry.boundary_distance(points, m['outer_faces'], p['faces'])
            assert d == geometry.F(local[3+3*i+2]['distance'])
            assert d > 0 and d >= geometry.F(float(c['minimum_clearance']))**2
            total = [a-geometry.F(float(c['final'][i]))*b for a,b in zip(total,part)]
        expected = [float(v*geometry.F(float(m['density']))) for v in total]
        measured = arithmetic.array(m['moments'], (10,))
        error = max(abs(measured-expected))/max(abs(v) for v in expected)
        assert expected[0] > 0 and error <= 1e-10
    else:
        assert not c['final_evaluation_available'] and 'material' not in c
        error = None
    if c['id'].startswith('mode'):
        assert c['source_outcome'] == 0 and c['final_status'] == 1 and accepted > 0
    if c['id'] in ('guard_cap', 'source_cap'):
        assert latched and c['final_status'] == 4 and accepted == 0
        assert c['initial'] == c['final']
    if c['id'] == 'initial_reject':
        assert guards[0]['outcome'] == 1 and c['source_evaluations'] == 0
    if c['id'] == 'restoration':
        assert accepted == 1 and any(t['restoration']['accepted'] for t in trace)
        assert geometry_rejected > 0 and shrunk > 0 and c['final_status'] == 1
    if c['id'] == 'exhaustion':
        assert accepted == 0 and c['initial'] == c['final']
        assert c['numerical_status'] == 'RESTORATION_STALLED' and c['final_status'] == 1
        assert trace[0]['restoration']['attempted'] and not trace[0]['restoration']['accepted']
        assert trace[0]['history_before'] == trace[0]['history_after'] == 0
        assert geometry_rejected > 0 and shrunk > 0
    if not quiet:
        print(f"PASS {c['id']} calls={len(calls)} accepted={accepted} geometry_reject={geometry_rejected} shrink_after_geometry_reject={shrunk} final={c['final_status']} moments_error={error}")
    return geometry_rejected, shrunk


def verify_all(cases, quiet=False):
    ids = [c['id'] for c in cases]
    assert len(ids) == len(EXPECTED) and set(ids) == set(EXPECTED)
    counts = [verify(c, quiet) for c in cases]
    assert sum(c[0] for c in counts) > 0 and sum(c[1] for c in counts) > 0


if __name__ == '__main__':
    geometry.self_test()
    verify_all([json.loads(line) for line in Path(sys.argv[1]).read_text().splitlines()])
    print('PASS 10 actual-source audits; independent linkage/KKT/moments/distances; no seal or reproduction grant')
