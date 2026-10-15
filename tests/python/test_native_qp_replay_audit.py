"""Independent fixed-QP checks cannot grant nonlinear convergence."""
import copy
import sys
from pathlib import Path
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'experiments/teapot_spin2014'))


def fixture():
    # Min -2*d0 + .5*(d0^2+d1^2), d0+d1=0, -.5<=d0<=.5.
    # The exact minimizer is (.5,-.5), equality multiplier .5, upper multiplier 1.
    return {'history': [], 'terminal': {
        'state': [.5, 0.], 'gradient': [-2., 0.], 'constraints': [0.],
        'constraint_jacobian': [[1., 1.]],
        'finite_box': {'checked': True, 'status': 'BOX_QP_FIRST_ORDER_CHECKED',
                       'direction': [.5, -.5], 'multipliers': [.5, 1.], 'working_set': [1]}}}


def test_replay_direction_has_independent_full_box_kkt_certificate():
    from qp_replay_audit import check_replay
    result = check_replay(fixture(), bounded=1)
    assert result['status'] == 'FIXED_QP_FIRST_ORDER_VERIFIED'
    assert result['box_error'] == result['stationarity'] == result['equality_error'] == 0.
    assert result['nonlinear_convergence'] == 'NOT_CHECKED'


@pytest.mark.parametrize('mutation', ['wrong_direction', 'wrong_multiplier', 'duplicate_face',
                                     'bad_history', 'unknown', 'nan', 'invalid_state'])
def test_invalid_qp_certificate_is_rejected(mutation):
    from qp_replay_audit import check_replay
    value = copy.deepcopy(fixture()); t = value['terminal']; q = t['finite_box']
    if mutation == 'wrong_direction': q['direction'] = [.6, -.6]
    elif mutation == 'wrong_multiplier': q['multipliers'][1] = -1.
    elif mutation == 'duplicate_face': q['working_set'].append(1); q['multipliers'].append(0.)
    elif mutation == 'bad_history': value['history'] = [{'s': [1., 0.], 'y': [-1., 0.]}]
    elif mutation == 'unknown': q['checked'] = False
    elif mutation == 'nan': t['gradient'][0] = float('nan')
    else: t['state'][0] = 1.1
    with pytest.raises(ValueError): check_replay(value, bounded=1)


def test_valid_nonidentity_bfgs_metric_and_lower_face():
    from qp_replay_audit import check_replay
    # One secant pair y=2*s gives B=2*identity in this direction.
    p = fixture(); p['history'] = [{'s': [1., 0.], 'y': [2., 0.]}]
    p['terminal']['gradient'] = [3., 0.]
    p['terminal']['finite_box'].update(direction=[-.5, .5], multipliers=[-1., -1.], working_set=[0])
    assert check_replay(p, bounded=1)['stationarity'] == 0.
    p['terminal']['finite_box']['multipliers'][1] = 1.
    # Keep stationarity exact, leaving only the incorrect lower-face dual sign.
    p['terminal']['gradient'][0] = 1.
    with pytest.raises(ValueError, match='INDEPENDENT_QP_KKT_FAILED'):
        check_replay(p, bounded=1)


def test_unlisted_face_cannot_hide_box_violation():
    from qp_replay_audit import check_replay
    p = fixture(); t = p['terminal']
    t['gradient'] = [-2., 0.]
    t['finite_box'].update(direction=[1., -1.], multipliers=[1.], working_set=[])
    # Equality and stationarity both exact; only the unlisted box face fails.
    with pytest.raises(ValueError, match='INDEPENDENT_QP_KKT_FAILED'):
        check_replay(p, bounded=1)


def test_equality_residual_is_checked_without_changing_stationarity():
    from qp_replay_audit import check_replay
    p = fixture(); p['terminal']['constraints'] = [1e-8]
    with pytest.raises(ValueError, match='INDEPENDENT_QP_KKT_FAILED'):
        check_replay(p, bounded=1)
