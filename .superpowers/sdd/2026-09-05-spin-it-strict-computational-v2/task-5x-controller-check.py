"""Independent arithmetic checks for actual adaptive loop exports.

This imports controller-only tetrahedral/BBW and KKT arithmetic, never the native
adaptive loop or its acceptance predicates. Known stationary control-flow cases
do not certify nontrivial optimization, a binary object, or paper reproduction.
"""
import hashlib
import importlib.util
import math
import json
from pathlib import Path
import struct
import sys


def load(name):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(name+'.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


material = load('task-5w-controller-check')
numeric = load('task-5t-controller-check')


def finite(value):
    if isinstance(value, dict):
        for item in value.values():
            finite(item)
    elif isinstance(value, list):
        for item in value:
            finite(item)
    elif isinstance(value, float):
        assert math.isfinite(value), 'nonfinite scientific record'


def sha(value):
    assert isinstance(value, str) and len(value) == 64
    assert set(value) <= set('0123456789abcdef')


def state_hash(value):
    return hashlib.sha256(struct.pack('>QQ', len(value), 1) +
                          struct.pack('>'+str(len(value))+'d', *value)).hexdigest()


def adaptive_state_hash(value):
    tag = b'source-adaptive-full-state-v1'
    return hashlib.sha256(struct.pack('>Q', len(tag)) + tag +
                          struct.pack('>Q', len(value)) +
                          struct.pack('>'+str(len(value))+'d', *value)).hexdigest()


def verify_guarded(g, source):
    finite(g)
    sha(g['hash'])
    assert g['initial'] == source['state']
    calls, guards, stages, trace = (g[k] for k in ('calls', 'guards', 'stages', 'trace'))
    assert len(calls) == len(guards)
    assert g['callback_count'] == sum(c['purpose'] == 0 for c in calls)
    assert g['source_calls'] == sum(c['invoked'] for c in calls)
    assert g['checker_calls'] == sum(c['invoked'] for c in guards)
    assert g['final_rechecks'] == sum(c['purpose'] == 1 and a['invoked'] for c, a in zip(calls, guards))
    n, cursor, valid_states = len(source['state']), 0, set()
    absent = 2**64-1
    pcount = len(source['primitives'])
    order = [(s, absent, absent) for s in (0, 1, 2)]
    order += [(s, p, absent) for p in range(pcount) for s in (3, 4, 5)]
    order += [(6, p, q) for p in range(pcount) for q in range(p+1, pcount)]
    for index, (c, a) in enumerate(zip(calls, guards)):
        assert c['sequence'] == a['callback'] == index
        assert len(c['state']) == n and c['state_hash'] == state_hash(c['state'])
        assert a['first_stage'] == cursor
        local = stages[cursor:cursor+a['stage_count']]
        cursor += len(local)
        assert len(local) == a['stage_count']
        if a['completed']:
            assert a['invoked'] and local
            assert [(s['stage'], s['first'], s['second']) for s in local] == order[:len(local)]
            if a['outcome'] == 0:
                assert len(local) == len(order)
                assert all(s['outcome'] == 1 or (s['stage'] == 6 and s['outcome'] == 5) for s in local)
            else:
                assert local[-1]['outcome'] == {1: 2, 2: 3, 3: 4}[a['outcome']]
        else:
            assert not local
        if c['status'] == 1:
            assert c['invoked'] and a['completed'] and a['outcome'] == 0
            sha(c['output_hash'])
            valid_states.add(tuple(c['state']))
        else:
            assert c['status'] in (2, 3) and c['reason'] and not c['output_hash']
    assert cursor == len(stages)
    assert g['numerical_available'] == bool(trace)
    if not trace:
        assert g['status'] == 'UNAVAILABLE' and not g['final_evaluation_available']
        return 0.
    assert trace[0]['state'] == g['initial']
    assert trace[0]['history_before'] == 0
    assert trace[0]['history_reset_reason'] == 'NEW_SOLVE_OR_PARTITION'
    assert trace[-1]['state'] == g['final']
    for index, t in enumerate(trace):
        assert t['iteration'] == index and tuple(t['state']) in valid_states
        assert t['final_state'] == (index == len(trace)-1)
        numeric.kkt(t, source['betas'])
        trials = t['objective_trials'] + t['restoration']['trials']
        assert len(trials) == t['line_search_evaluations']
        assert sum(trial['accepted'] for trial in trials) == int(t['accepted_step'] > 0)
        for trial in trials:
            assert 0 < trial['step'] <= 1
            assert not trial['accepted'] or trial['valid']
    kkt = numeric.kkt(trace[-1], source['betas'])
    if g['status'] == 'CONVERGED':
        assert kkt['equality_primal'] <= g['controls']['equality_tolerance']
        assert all(kkt[k] <= g['controls']['projected_gradient_tolerance'] for k in
                   ('bound_primal', 'stationarity', 'dual_sign', 'complementarity'))
    if g['final_evaluation_available']:
        assert calls[-1]['purpose'] == 1 and calls[-1]['status'] == 1
        assert calls[-1]['state'] == g['final'] == g['material']['state']
        assert calls[-1]['output_hash'] == g['material']['evaluation_hash']
        return material.check_source(g['material'])
    assert 'material' not in g
    return 0.


def verify_rebound(old, new):
    assert old['source_hash'] == new['source_hash']
    assert old['bbw_hash'] == new['bbw_hash']
    assert old['weights'] == new['weights'] and old['centers'] == new['centers']
    assert old['density'] == new['density']
    assert old['state'][old['betas']:] == new['state'][new['betas']:]
    for i, leaf in enumerate(new['leaves']):
        ancestors = [j for j, item in enumerate(old['leaves']) if material.ancestor(material.key(item), material.key(leaf))]
        if ancestors:
            assert len(ancestors) == 1 and new['state'][i] == old['state'][ancestors[0]]
        else:
            children = [j for j, item in enumerate(old['leaves']) if material.ancestor(material.key(leaf), material.key(item))]
            assert children and all(old['state'][j] == new['state'][i] for j in children)
    return material.check_source(new)


EXPECTED = {
    'stationary_multi': 'RELAXED_ADAPTIVE_COMPLETE',
    'stationary_pass_limit': 'MAXIMUM_PASSES',
    'nontrivial_stop': 'SOLVE_NOT_CONVERGED',
    'initial_geometry_reject': 'SOLVE_GEOMETRY_REJECTED',
    'rebind_geometry_reject': 'REBIND_GEOMETRY_REJECTED',
    'aggregate_limit': 'RESOURCE_LIMIT',
}


def verify(c):
    finite(c)
    assert c['id'] in EXPECTED and c['outcome'] == EXPECTED[c['id']]
    assert c['scope'] == ('SYNTHETIC_REGRESSION' if c['id'] == 'nontrivial_stop' else 'KNOWN_SOLUTION_VALIDATION')
    assert c['sealed'] is False and c['paper_reproduction'] is False
    sha(c['hash'])
    assert c['reason'] and 0 < len(c['passes']) <= c['maximum_passes']
    assert 0 <= c['guarded_reserved'] + c['rebind_reserved'] <= c['known_buffer_bytes'] <= c['known_buffer_cap']
    assert 0 < c['retained_pairs'] <= c['admitted_pairs'] <= c['maximum_passes']+1
    previous, solved, pending, cumulative, error = None, None, None, 0, 0.
    for index, p in enumerate(c['passes']):
        assert p['number'] == index+1
        assert cumulative <= p['known_buffer_bytes'] <= c['known_buffer_bytes']
        cumulative = p['known_buffer_bytes']
        source = p['source']
        assert source['density'] == 1.17
        error = max(error, material.check_source(source))
        assert p['old_partition'] == source['partition_hash']
        sha(p['old_functional']); sha(p['old_state_hash'])
        assert p['old_state_hash'] == adaptive_state_hash(source['state'])
        if previous is not None:
            assert previous['rebind']['source'] == source
            assert previous['old_partition'] != previous['new_partition']
            assert previous['new_functional'] == p['old_functional']
            assert previous['new_state_hash'] == p['old_state_hash']
        if 'guarded' not in p:
            assert index == len(c['passes'])-1 and c['outcome'] == 'RESOURCE_LIMIT'
            assert 'rebind' not in p and not p['fresh_history']
            break
        g = p['guarded']
        error = max(error, verify_guarded(g, source))
        assert p['fresh_history'] == g['numerical_available']
        authorized = g['status'] == 'CONVERGED' and g['final_status'] == 1 and g['final_evaluation_available']
        if authorized:
            assert g['final_rechecks'] == 1 and g['outcome'] == 0
            solved, pending = g['material'], None
            sha(p['solved_state_hash'])
            assert p['solved_state_hash'] == adaptive_state_hash(solved['state'])
        if 'rebind' in p:
            assert authorized
            r = p['rebind'];sha(r['hash'])
            if 'source' in r:
                assert r['outcome'] == 0
                error = max(error, verify_rebound(g['material'], r['source']))
                assert r['source']['partition_hash'] == p['new_partition']
                sha(p['new_functional']);sha(p['new_state_hash'])
                assert p['new_state_hash'] == adaptive_state_hash(r['source']['state'])
                if p['new_partition'] != p['old_partition']:
                    pending = r['source']
                else:
                    assert r['source'] == g['material']
                    assert p['old_functional'] == p['new_functional']
                    assert p['solved_state_hash'] == p['new_state_hash']
            else:
                assert r['outcome'] in (1, 2, 3)
                assert index == len(c['passes'])-1
        else:
            assert index == len(c['passes'])-1
        previous = p
    assert c.get('last_solved') == solved
    assert c.get('pending') == pending
    if c['outcome'] == 'RELAXED_ADAPTIVE_COMPLETE':
        assert solved is not None and pending is None and c['counts_available']
        frac = [(leaf, beta) for leaf, beta in zip(solved['leaves'], solved['state'][:solved['betas']])
                if c['epsilon'] <= beta <= 1-c['epsilon']]
        refinable = [leaf for leaf, _ in frac if leaf['level'] < c['maximum_level']]
        assert c['fractional'] == len(frac) and c['refinable'] == len(refinable) == 0
        # Two genuine partition changes with identical state dimensions are the
        # specific independent regression; a beta-count-only loop must fail it.
        assert len(c['passes']) == 3
        assert len({p['source']['betas'] for p in c['passes']}) == 1
        assert len({p['old_partition'] for p in c['passes']}) == 3
    if c['outcome'] == 'MAXIMUM_PASSES':
        assert solved is not None and pending is not None
        assert solved['partition_hash'] != pending['partition_hash']
        assert len(c['passes']) == c['maximum_passes'] == 1
    if c['outcome'] == 'SOLVE_NOT_CONVERGED':
        assert solved is None and pending is None
        assert len(c['passes']) == 1 and 'rebind' not in c['passes'][0]
        assert c['passes'][0]['guarded']['status'] != 'CONVERGED'
    if c['outcome'] == 'SOLVE_GEOMETRY_REJECTED':
        assert len(c['passes']) == 1 and 'rebind' not in c['passes'][0]
        g = c['passes'][0]['guarded']
        assert not g['numerical_available'] and g['guards'][0]['outcome'] == 1
        assert solved is None and pending is None
    if c['outcome'] == 'REBIND_GEOMETRY_REJECTED':
        assert len(c['passes']) == 1
        r = c['passes'][0]['rebind']
        assert r['outcome'] == 1 and r['last_stage'] == 1 and 'source' not in r
        assert solved is not None and pending is None
    if c['outcome'] == 'RESOURCE_LIMIT':
        # Fixed 320MiB protocol stops at first guarded reservation, not after
        # one pass. Keep this input; do not tune it to obtain a desired failure.
        assert solved is None and pending is None
        assert len(c['passes']) == 1 and 'guarded' not in c['passes'][0]
        assert c['guarded_reserved'] == c['rebind_reserved'] == 0
    return error


def verify_rows(rows):
    assert len(rows) == len(EXPECTED) and {c['id'] for c in rows} == set(EXPECTED)
    assert len({c['id'] for c in rows}) == len(rows)
    return max(verify(c) for c in rows)


if __name__ == '__main__':
    material.geo.self_test()
    rows = [json.loads(line) for line in Path(sys.argv[1]).read_text().splitlines()]
    error = verify_rows(rows)
    print(f'PASS: {len(rows)} actual adaptive cases; fresh history/KKT, raw guard callbacks, partition/remap chain, independent BBW/tetrahedral moments; max error {error:.17g}')
