"""Exercise two Newton steps through the loaded production native entry.

This is an algebraic solver regression, not a material or Heart optimization.
The second step must match a fresh solve on the new affine branch and differ
from the same branch supplied with the previous branch's stale Jacobian.
"""
import json
from dataclasses import replace
from pathlib import Path

import numpy as np

from make_it_stand import spin

CORE = Path(__file__).resolve().parents[2]/'experiments/faster_core_repair_20261012'


def test_branch_transition_next_step_consumes_new_jacobian():
    config = json.loads((CORE/'core-benchmark-inputs-job-01/output/small.json').read_text())
    request_data = dict(config['request'])
    settings = spin.SpinSolverConfig(maximum_iterations=2, maximum_restarts=1,
        maximum_function_evaluations_per_restart=60, equation_profile='DUAL_MONOMIAL_RAW_V1',
        capture_evaluation_payload=True)
    request_data.update(mode='STATIC_BALANCE', static_objective='MAXIMUM_MASS', settings=settings)
    request = spin.SpinItFasterRequest(**request_data)
    initial = np.array([1., .12, -.08, .09, -.07, .05, -.04, .03, -.02, .01])
    calls = []
    anchor = None
    def material(q, branch, raw, jacobian):
        return dict(raw_moments=raw, moment_jacobian=jacobian, derivative_defined=True,
            geometry_hash=spin._native_module().spin_it_faster_quadric_hash(q), topology_hash=branch,
            evaluation_map_identity='a'*64, boolean_backend_identity='algebraic-branch-regression',
            derivative_identity='affine-branch-jacobian')
    def crossing(q):
        nonlocal anchor
        calls.append(q)
        point = np.asarray(q)
        if len(calls) == 1:
            return material(q, 'branch-A', point, np.eye(10))
        if anchor is None:
            anchor = point.copy()
        return material(q, 'branch-B', 2*point-anchor, 2*np.eye(10))
    result = spin.solve_spin_it_faster_with_external_moment_evaluator(request, crossing,
        initial_quadrics=(initial,))
    accepted = [row for row in result['events'] if row['type'] == 'evaluation' and
                row['outcome'] == 'accepted' and (row['step'] or 0) > 0]
    assert len(accepted) == 2
    assert accepted[0]['reason'] == 'topology_transition_relinearized'
    np.testing.assert_array_equal(accepted[0]['quadric_coefficients'], anchor)
    one_step = replace(request, settings=replace(settings, maximum_iterations=1))
    def reference(jacobian):
        run = spin.solve_spin_it_faster_with_external_moment_evaluator(one_step,
            lambda q: material(q, 'branch-B', 2*np.asarray(q)-anchor, jacobian),
            initial_quadrics=(anchor,))
        updates = [row for row in run['events'] if row['type'] == 'evaluation' and
                   row['outcome'] == 'accepted' and (row['step'] or 0) > 0]
        assert len(updates) == 1
        return np.asarray(updates[0]['quadric_coefficients'])
    actual = np.asarray(accepted[1]['quadric_coefficients'])
    fresh = reference(2*np.eye(10)); stale = reference(np.eye(10))
    np.testing.assert_allclose(actual, fresh, rtol=0, atol=1e-13)
    assert np.max(np.abs(actual-stale)) > 1e-6, 'fixture cannot distinguish a stale Jacobian'
