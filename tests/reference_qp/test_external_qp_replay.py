"""External QP comparisons do not grant nonlinear or paper acceptance."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'experiments/teapot_spin2014'))


def fixture():
    # Strictly convex analytic test, first coordinate boxed, second free.
    return {'history': [], 'terminal': {
        'state': [.5, 0.], 'gradient': [-2., 0.], 'constraints': [0.],
        'constraint_jacobian': [[1., 1.]],
        'finite_box': {'checked': True, 'status': 'BOX_QP_FIRST_ORDER_CHECKED',
                       'direction': [.5, -.5], 'multipliers': [.5, 1.], 'working_set': [1]}}}


def test_upper_box_equality_and_free_variable():
    from external_qp_replay import solve_replay
    result = solve_replay(fixture(), bounded=1)
    assert result['verification']['status'] == 'FIXED_QP_FIRST_ORDER_VERIFIED'
    np.testing.assert_allclose(result['direction'], [.5, -.5], atol=1e-10)
    assert result['nonlinear_convergence'] == 'NOT_CHECKED'
    assert result['scope'] == 'EXTERNAL_FIXED_QP_COMPARISON_ONLY'


def test_nonidentity_bfgs_lower_bound():
    from external_qp_replay import solve_replay
    value = fixture()
    value['history'] = [{'s': [1., 0.], 'y': [2., 0.]}]
    value['terminal']['gradient'] = [3., 0.]
    result = solve_replay(value, bounded=1)
    np.testing.assert_allclose(result['direction'], [-.5, .5], atol=1e-10)
    assert result['box_dual'][0] < 0


def test_asymmetric_multisecant_hessian_matches_direct_action():
    from external_qp_replay import dense_hessian
    from direction_replay import hessian_action
    pairs = [(np.array([1., 2., -.4]), np.array([2., 1., -.3])),
             (np.array([-.2, 1., 1.]), np.array([.1, 3., 2.]))]
    h = dense_hessian(3, pairs)
    for vector in (np.array([.3, -.7, 1.]), np.array([-1., 2., -.8])):
        np.testing.assert_allclose(h @ vector, hessian_action(vector, pairs), atol=1e-13)
    assert np.linalg.eigvalsh(h).min() > 0


@pytest.mark.parametrize('mutation', ['nan', 'bad_history', 'bad_shape', 'bad_box', 'bounded_bool'])
def test_malformed_input_rejected_before_external_solver(mutation):
    from external_qp_replay import solve_replay
    value = fixture()
    bounded = 1
    if mutation == 'nan': value['terminal']['gradient'][0] = float('nan')
    elif mutation == 'bad_history': value['history'] = [{'s': [1., 0.], 'y': [-1., 0.]}]
    elif mutation == 'bad_shape': value['terminal']['constraint_jacobian'] = [[1.]]
    elif mutation == 'bad_box': value['terminal']['state'][0] = 1.1
    else: bounded = True
    with pytest.raises(ValueError): solve_replay(value, bounded=bounded)


@pytest.mark.parametrize('mutation', ['primal', 'dual', 'equality', 'free_dual', 'nan'])
def test_external_status_cannot_override_independent_kkt(mutation):
    from external_qp_replay import check_external_solution
    value = fixture()
    d, y, z = [.5, -.5], [.5], [1., 0.]
    if mutation == 'primal': d = [.6, -.6]
    elif mutation == 'dual': z = [-1., 0.]
    elif mutation == 'equality': y = [.6]
    elif mutation == 'free_dual': z[1] = .1
    else: d[0] = float('nan')
    with pytest.raises(ValueError):
        check_external_solution(value, 1, d, y, z)


def test_does_not_read_saved_native_solution_as_initial_guess():
    from external_qp_replay import solve_replay
    value = fixture()
    value['terminal']['finite_box'] = {'direction': [float('nan')], 'checked': False}
    result = solve_replay(value, bounded=1)
    assert result['verification']['status'] == 'FIXED_QP_FIRST_ORDER_VERIFIED'


def test_dense_allocation_budget_is_enforced_before_solving():
    from external_qp_replay import solve_replay
    with pytest.raises(ValueError, match='DENSE_DIMENSION_LIMIT'):
        solve_replay(fixture(), bounded=1, maximum_dimension=1)


def test_result_revalidation_recomputes_kkt_not_recorded_boolean():
    from external_qp_replay import solve_replay, verify_scientific_result
    result = solve_replay(fixture(), bounded=1)
    result['direction'][0] += .1
    result['verification'] = {'status': 'FIXED_QP_FIRST_ORDER_VERIFIED'}
    with pytest.raises(ValueError, match='INDEPENDENT_QP_KKT_FAILED'):
        verify_scientific_result(result, fixture())


def test_file_snapshot_rejects_changed_source(tmp_path):
    from external_qp_replay import file_hash, check_snapshot
    source = tmp_path / 'checker.py'
    source.write_text('original')
    snapshot = {str(source): file_hash(source)}
    check_snapshot(snapshot)
    source.write_text('changed')
    with pytest.raises(ValueError, match='STALE_RESULT'):
        check_snapshot(snapshot)


def test_binary_inventory_includes_transitive_non_python_libraries():
    from external_qp_replay import mapped_shared_objects
    paths = mapped_shared_objects()
    assert any(Path(path).name.startswith('libc.so') for path in paths)
    assert all(Path(path).is_file() for path in paths)
