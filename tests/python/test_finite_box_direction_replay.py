import copy
import numpy as np
import pytest
from test_slqp_direction_replay import api, replay_fixture
from make_it_stand.verification.source_trace_archive import write_archive


def fixture():
    events = replay_fixture()[:4]
    events[0]['trace_counts'] = [2]
    events[0]['payload']['settings']['solver_settings']['slqp']['direction_policy'] = 'FINITE_BOX_QP_V1'
    for i, event in enumerate(events[2:]):
        event['payload'].update(state=([.8,.2,0.] if i==0 else [.9,.1,0.]),
            gradient=([0.,2.,0.] if i==0 else [.1,1.9,0.]),
            constraints=[0.], constraint_jacobian=[[1.,1.,0.]], kkt={'lambda':[-1.]},
            active_membership=[0,0,0], history_before=i, lbfgs_history_size=i,
            final_state=bool(i), accepted_step=.5 if i==0 else 0.,
            history_reset_reason='NEW_SOLVE_OR_PARTITION' if i==0 else '')
    p=events[2]['payload']
    p['finite_box']=dict(attempted=True,checked=True,selected=True,accepted=True,
        status='BOX_QP_FIRST_ORDER_CHECKED',direction=[.2,-.2,0.],multipliers=[-1.8,1.6],
        working_set=[1],equality_residual=0.,stationarity=0.,trial_count=1,fallback_attempted=False,
        trace=[dict(working_set=[1],action='CHECKED',changed_index=-1,equality_residual=0.,stationarity=0.)])
    p['objective_trials']=[dict(step=.5,accepted=True,valid=True,reason='ACCEPTED',constraint_l1=0.)]
    return events


def test_replay_uses_actual_full_box_direction_not_old_clipped_qn(tmp_path):
    path=tmp_path/'box.jsonl';write_archive(path,fixture())
    result=api().replay(path)['passes'][0]
    assert result['accepted_step_comparisons']==1
    assert result['finite_box_accepted_step_comparisons']==1
    assert result['max_step_replay_error']<1e-14


@pytest.mark.parametrize('field,value',[
    ('direction',[.3,-.3,0.]),('multipliers',[-1.8,-1.6]),('working_set',[0]),
    ('working_set',[1,1]),('checked',False),('status','BOX_QP_UNKNOWN'),
    ('stationarity',float('nan')),('trial_count',0),('trace',[])])
def test_full_box_false_certificate_rejected(field,value):
    p=fixture()[2]['payload'];p['finite_box'][field]=value
    assert hasattr(api(),'checked_finite_box_direction'), 'missing independent finite-box verifier'
    with pytest.raises(ValueError):
        api().checked_finite_box_direction(p,np.array(p['state']),np.array(p['gradient']),
            np.array(p['constraint_jacobian']),np.array(p['constraints']),2,[],1e-10,1e-9)


def test_box_event_cannot_be_relabelled_as_default_profile(tmp_path):
    events=fixture();events[0]['payload']['settings']['solver_settings']['slqp'].pop('direction_policy')
    path=tmp_path/'wrong-profile.jsonl';write_archive(path,events)
    with pytest.raises(ValueError):api().replay(path)


def test_full_policy_cannot_accept_legacy_steps_by_clearing_both_qp_flags(tmp_path):
    events=replay_fixture()
    events[0]['payload']['settings']['solver_settings']['slqp']['direction_policy']='FINITE_BOX_QP_V1'
    for event in events[2:]:
        if event['payload']['accepted_step']>0:
            event['payload']['finite_box']=dict(attempted=False,checked=False,selected=False,
                accepted=False,status='BOX_QP_UNKNOWN')
    path=tmp_path/'cleared-flags.jsonl';write_archive(path,events)
    with pytest.raises(ValueError):api().replay(path)


@pytest.mark.parametrize('active,residual',[(1e-10,float('nan')),(float('nan'),1e-9),(-1.,1e-9),(1e-10,0.)])
def test_checker_rejects_invalid_tolerances(active,residual):
    p=fixture()[2]['payload']
    with pytest.raises(ValueError):
        api().checked_finite_box_direction(p,np.array(p['state']),np.array(p['gradient']),
            np.array(p['constraint_jacobian']),np.array(p['constraints']),2,[],active,residual)


@pytest.mark.parametrize('d0',[.5-5e-11,np.nextafter(.5,0.)])
def test_small_face_slack_with_large_multiplier_cannot_hide_complementarity(d0):
    p=fixture()[2]['payload']
    p['finite_box'].update(direction=[d0,0.,0.],multipliers=[0.,1e10])
    with pytest.raises(ValueError):
        api().checked_finite_box_direction(p,np.array([.5,.5,0.]),np.array([-1e10-d0,0.,0.]),
            np.array([[0.,0.,1.]]),np.array([0.]),2,[],1e-10,1e-9)
