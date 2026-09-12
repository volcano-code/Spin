"""Independent arithmetic audit of exported actual-source SLQP trajectories.

No native optimizer, witness LP, hashAuditedSlqp or geometry integrator import.
Source endpoint integration is checked separately by the Task5s QR/tetra oracle.
These synthetic trajectories cannot grant geometry or paper reproduction.
"""
import copy
import hashlib
import json
import math
import struct
import sys
from pathlib import Path
import numpy as np


def array(value, shape=None):
    result = np.asarray(value, dtype=float)
    assert np.isfinite(result).all()
    if shape is not None:
        assert result.shape == shape
    return result


def inf(value):
    return float(np.max(np.abs(value))) if len(value) else 0.


def kkt(entry, bounded):
    x, g, c = (array(entry[k]) for k in ('state', 'gradient', 'constraints'))
    n, m = len(x), len(c)
    assert x.shape == g.shape == (n,) and n > bounded > 0
    j = array(entry['constraint_jacobian'], (m, n))
    witness = entry['kkt']
    assert witness['valid'] is True
    lam, alpha, eta = (array(witness[k], (count,)) for k, count in
                       [('lambda', m), ('alpha', n), ('eta', n)])
    assert np.array_equal(alpha[bounded:], np.zeros(n-bounded))
    assert np.array_equal(eta[bounded:], np.zeros(n-bounded))
    stationarity = g + j.T @ lam - alpha + eta
    expected = {
        'equality_primal': inf(c),
        'bound_primal': max(0., float(np.max(-x[:bounded])), float(np.max(x[:bounded]-1))),
        'stationarity': inf(stationarity),
        'dual_sign': max(0., float(np.max(-alpha)), float(np.max(-eta))),
        'complementarity': max(inf(alpha[:bounded]*x[:bounded]),
                               inf(eta[:bounded]*(1-x[:bounded]))),
    }
    expected['full'] = max(expected.values())
    for key, want in expected.items():
        got = float(witness[key])
        assert math.isfinite(got) and abs(got-want) <= 1e-10*(1+want), (key, got, want)
    return expected


def verify(case):
    fixture, mode = case['fixture'], case['mode']
    assert fixture in (0, 1) and mode in range(5)
    nb, nh = case['beta_count'], case['handle_count']
    assert (nb, nh) == (8, 26)
    angle = mode != 4
    n = nb + int(angle) + 4*nh
    assert case['geometry_status'] == 'NOT_EVALUATED'
    assert case['feasibility_status'] == 'UNKNOWN'
    assert case['algorithm'] == 'AUDITED_KKT_RESTORATION_V1'
    initial, final = (array(case[k], (n,)) for k in ('initial', 'final'))
    controls = case['controls']
    eqtol, ktol = controls['equality_tolerance'], controls['projected_gradient_tolerance']
    assert eqtol == 1e-8 and ktol == 1e-7
    for name in ('source_hash', 'config_hash', 'result_hash', 'trace_hash'):
        assert len(case[name]) == 64 and all(c in '0123456789abcdef' for c in case[name])
    trace = case['trace']
    assert trace and len(trace) <= controls['maximum_iterations']+1
    assert np.array_equal(array(trace[0]['state']), initial)
    assert np.array_equal(array(trace[-1]['state']), final)
    assert trace[-1]['final_state'] is True
    assert trace[0]['history_before'] == 0
    assert 'NEW_SOLVE_OR_PARTITION' in trace[0]['history_reset_reason']
    calls = case['calls']
    assert calls and calls[-1]['purpose'] == 1 and all(c['purpose'] == 0 for c in calls[:-1])
    assert case['callback_count'] == len(calls)-1
    assert case['source_evaluation_count'] == sum(c['source_invoked'] for c in calls)
    assert np.array_equal(array(calls[-1]['state']), final)
    valid_states = set()
    for index, call in enumerate(calls):
        assert call['sequence'] == index and call['status'] in (1, 2)
        state = array(call['state'], (n,))
        digest = hashlib.sha256(struct.pack('>QQ', n, 1)+struct.pack('>'+str(n)+'d', *state)).hexdigest()
        assert call['state_hash'] == digest and call['source_invoked'] is True
        if call['status'] == 1:
            assert len(call['output_hash']) == 64 and not call['reason']
            valid_states.add(tuple(state))
        else:
            assert call['reason'] and len(call['reason_hash']) == 64 and not call['output_hash']
    assert calls[-1]['status'] == 1
    accepted = 0
    for i, t in enumerate(trace):
        x = array(t['state'], (n,))
        assert tuple(x) in valid_states
        residuals = kkt(t, nb)
        assert (x[:nb] >= 0).all() and (x[:nb] <= 1).all()
        assert t['phi'] == (x[nb] if angle else 0.)
        assert t['iteration'] == i
        assert abs(t['kkt_infinity']-residuals['full']) <= 1e-10*(1+residuals['full'])
        objective = float(t['objective'])
        assert math.isfinite(objective)
        merit = objective + controls['merit_penalty']*float(np.abs(array(t['constraints'])).sum())
        assert math.isclose(t['merit'], merit, rel_tol=1e-12, abs_tol=1e-12)
        membership = t['active_membership']
        assert len(membership) == n and set(membership) <= {0, 1, 2}
        assert membership[nb:] == [0]*(n-nb)
        assert t['membership_hash'] == hashlib.sha256(
            b'spin2014-active-membership-v1'+bytes(membership)).hexdigest()
        for index, active in enumerate(membership[:nb]):
            if active == 1: assert x[index] <= controls['active_set_tolerance']
            if active == 2: assert x[index] >= 1-controls['active_set_tolerance']
        if i and membership != trace[i-1]['active_membership']:
            assert 'ACTIVE_MEMBERSHIP_CHANGED' in t['history_reset_reason']
        assert not (t['curvature_pair_accepted'] and t['curvature_pair_rejected'])
        assert 0 <= t['history_after'] <= controls['lbfgs_memory']
        objective_trials = t['objective_trials']
        restoration = t['restoration']
        all_trials = objective_trials+restoration['trials']
        assert t['line_search_evaluations'] == len(all_trials)
        for trial in all_trials:
            assert 0 < trial['step'] <= 1 and trial['reason']
            assert not trial['accepted'] or trial['valid']
            if trial['valid']: assert math.isfinite(trial['constraint_l1'])
        if t['accepted_step'] > 0:
            assert i+1 < len(trace)
            accepted += 1
            nxt = trace[i+1]
            accepted_trials = [v for v in all_trials if v['accepted']]
            assert len(accepted_trials) == 1
            assert accepted_trials[0]['step'] == t['accepted_step']
            next_l1 = float(np.abs(array(nxt['constraints'])).sum())
            assert math.isclose(accepted_trials[0]['constraint_l1'], next_l1, rel_tol=1e-10, abs_tol=1e-12)
            if restoration['accepted']:
                direction = array(restoration['direction'], (n,))
                assert np.allclose(array(nxt['state']), x+t['accepted_step']*direction, rtol=1e-12, atol=1e-12)
                old_l1 = float(np.abs(array(t['constraints'])).sum())
                predicted = restoration['predicted_residual']
                assert old_l1 > predicted and next_l1 < old_l1
                assert next_l1 <= old_l1-controls['armijo']*t['accepted_step']*(old_l1-predicted)+1e-12
            else:
                assert nxt['merit'] < t['merit']
        else:
            assert not any(v['accepted'] for v in all_trials)
    final_kkt = kkt(trace[-1], nb)
    assert case['initial_phi'] == (initial[nb] if angle else 0.)
    assert case['final_phi'] == (final[nb] if angle else 0.)
    if case['status'] == 'CONVERGED':
        assert final_kkt['equality_primal'] <= eqtol
        assert all(final_kkt[k] <= ktol for k in ('bound_primal', 'stationarity', 'dual_sign', 'complementarity'))
    else:
        assert case['status'] in ('MAXIMUM_ITERATIONS', 'LINE_SEARCH_EXHAUSTED',
                                  'RESTORATION_STALLED', 'NUMERICAL_UNKNOWN')
    assert accepted > 0, 'fixture must exercise real progress, not only a known stationary start'
    assert inf(final[nb+int(angle):]-initial[nb+int(angle):]) > 0
    assert case['numerical_only'] is True
    print(f'PASS fixture={fixture} mode={mode} states={len(trace)} accepted={accepted} '
          f'status={case["status"]} final_kkt={final_kkt["full"]:.6g}')
    return fixture, mode


def main(path):
    cases = [json.loads(line) for line in Path(path).read_text().splitlines()]
    seen = [verify(c) for c in cases]
    assert len(seen) == 10 and set(seen) == {(f, m) for f in (0, 1) for m in range(5)}
    print('PASS 10 actual-source joint traces; numerical-only, no geometry or reproduction grant.')


if __name__ == '__main__':
    main(sys.argv[1])
