import importlib.util
import copy
import sys
from pathlib import Path
import numpy as np
import pytest
from make_it_stand.verification.source_trace_archive import write_archive

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'experiments/teapot_spin2014'))
def api():
    assert importlib.util.find_spec('direction_replay'), 'missing audited direction replay'
    import direction_replay
    return direction_replay

def test_spd_metric_can_point_out_of_box_despite_inward_gradient():
    m=api()
    # One exact BFGS inverse update from I: H=[[1,2],[2,9]].
    pairs=[(np.array([1.,2.]),np.array([1.,0.]))]
    np.testing.assert_allclose(m.inverse_action(np.eye(2),pairs),[[1,2],[2,9]],rtol=0,atol=1e-14)
    np.testing.assert_allclose(m.hessian_action(np.eye(2),pairs),[[1.8,-.4],[-.4,.2]],rtol=0,atol=1e-14)
    d=m.quasi_newton_direction(np.array([-1.,1.]),np.empty((0,2)),np.empty(0),[0,0],pairs)
    np.testing.assert_allclose(d,[-1,-7],rtol=0,atol=1e-14)
    assert np.dot([-1.,1.],d)<0
    assert m.bound_step(np.array([0.,.5]),d,2)==0
    assert m.bound_step(np.array([0.,.5]),np.array([1.,-1.]),2)==.5

def test_linearized_equality_and_active_coordinate_preserved():
    m=api()
    d=m.quasi_newton_direction(np.array([3.,2.,4.]),np.array([[1.,1.,1.]]),np.array([.2]),[1,0,0],[])
    np.testing.assert_allclose(d,[0,.9,-1.1],rtol=0,atol=1e-14)
    assert m.bound_step(np.array([0.,.5,8.]),np.array([0.,.25,-20.]),2)==1

def test_nonfinite_or_nonpositive_curvature_cannot_be_replayed():
    m=api()
    for y in (np.array([-1.,0.]),np.array([np.nan,0.])):
        with pytest.raises(ValueError):m.inverse_action(np.ones(2),[(np.array([1.,0.]),y)])

def test_history_keeps_old_multipliers_in_secant_difference():
    m=api()
    a=dict(state=[.2,.4],gradient=[1.,2.],constraint_jacobian=[[1.,0.]],kkt={'lambda':[3.]})
    b=dict(state=[.3,.2],gradient=[2.,4.],constraint_jacobian=[[2.,1.]],kkt={'lambda':[999.]})
    s,y=m.secant_pair(a,b)
    np.testing.assert_allclose(s,[.1,-.2],atol=1e-15)
    np.testing.assert_allclose(y,[4.,5.],atol=1e-15)

@pytest.mark.parametrize('upper',[False,True])
def test_tangent_qp_handles_bound_in_the_same_spd_metric(upper):
    m=api()
    assert hasattr(m,'tangent_qp_direction'), 'missing metric-consistent box working set'
    sign=-1. if upper else 1.
    pairs=[(np.array([1.,sign*2.]),np.array([1.,0.]))]
    x=np.array([float(upper),.5]);g=np.array([-sign,1.])
    result=m.tangent_qp_direction(x,g,np.empty((0,2)),np.empty(0),2,pairs)
    assert result['status']=='QP_FIRST_ORDER_CHECKED'
    np.testing.assert_allclose(result['direction'],[0.,-5.],atol=1e-12)
    assert m.bound_step(x,result['direction'],2)==pytest.approx(.1)

def test_tangent_qp_releases_wrong_identity_metric_active_face():
    m=api();assert hasattr(m,'tangent_qp_direction')
    # Identity gradient would freeze x0 (g0>0); the SPD metric optimum is inward.
    r=m.tangent_qp_direction(np.array([0.,.5]),np.array([1.,-1.]),
        np.empty((0,2)),np.empty(0),2,[(np.array([1.,2.]),np.array([1.,0.]))])
    assert r['status']=='QP_FIRST_ORDER_CHECKED'
    np.testing.assert_allclose(r['direction'],[1.,7.],atol=1e-12)

def test_tangent_qp_preserves_equality_and_does_not_certify_infeasibility():
    m=api();assert hasattr(m,'tangent_qp_direction')
    r=m.tangent_qp_direction(np.array([0.,.5]),np.array([-1.,1.]),
        np.array([[1.,1.]]),np.array([0.]),2,[(np.array([1.,2.]),np.array([1.,0.]))])
    assert r['status']=='QP_FIRST_ORDER_CHECKED'
    np.testing.assert_allclose(r['direction'],[5/7,-5/7],atol=1e-12)
    r=m.tangent_qp_direction(np.array([0.]),np.array([1.]),np.array([[1.]]),np.array([1.]),1,[])
    assert r['status']=='QP_UNKNOWN'

def test_fixed_direction_is_exactly_zero_not_roundoff_outward():
    m=api()
    s=np.array([.1257302210933933,-.1321048632913019,.6404226504432821])
    y=np.array([.1257302210933933,-.2642097265826038,1.9212679513298463])
    g=np.array([.10490011715303971,-.535669373161111,.36159505490948474])
    x=np.array([0.,.5,.5]);a=np.array([[0.,1.,1.]])
    r=m.tangent_qp_direction(x,g,a,np.array([0.]),3,[(s,y)])
    assert r['status']=='QP_FIRST_ORDER_CHECKED'
    assert r['direction'][0]==0
    assert np.max(np.abs(a@r['direction']))<1e-12
    assert m.bound_step(x,r['direction'],3)>0

def test_tangent_qp_rejects_invalid_problem_instead_of_issuing_kkt_status():
    m=api()
    for x,g in ((np.array([np.nan]),np.array([1.])),(np.array([-.1]),np.array([1.])),
                (np.array([0.]),np.array([1.,2.]))):
        with pytest.raises(ValueError):m.tangent_qp_direction(x,g,np.empty((0,1)),np.empty(0),1,[])

def replay_fixture():
    events=[dict(kind='header',trace_counts=[3],payload={'settings':{'solver_settings':{'slqp':{'lbfgs_memory':8}}}}),
        dict(kind='pass',pass_index=0,payload={'index':0,'status':'CONVERGED','bounded_variable_count':2})]
    # Literal quadratic f=x0^2+(x1-1)^2, equality x0+x1=1.
    for i,(x,g,step,hist,applied,reason,active) in enumerate([
            ([.5,.5],[1.,-1.],.25,0,0,'NEW_SOLVE_OR_PARTITION',[0,0]),
            ([.25,.75],[.5,-.5],1.,1,1,'',[0,0]),
            ([0.,1.],[0.,0.],0.,2,0,'ACTIVE_MEMBERSHIP_CHANGED',[1,2])]):
        p=dict(state=x,gradient=g,constraint_jacobian=[[1.,1.]],constraints=[0.],
            kkt={'lambda':[0.]},history_before=hist,lbfgs_history_size=applied,
            history_reset_reason=reason,final_state=i==2,active_membership=active,
            used_first_order_fallback=False,restoration={'accepted':False},
            curvature_pair_accepted=i!=2,accepted_step=step,iteration=i)
        events.append(dict(kind='iteration',pass_index=0,trace_index=i,payload=p))
    return events

def test_complete_replay_and_history_or_state_mismatch_rejection(tmp_path):
    m=api();events=replay_fixture();path=tmp_path/'trace.jsonl';write_archive(path,events)
    result=m.replay(path)
    assert result['verified_records']==5
    assert result['passes'][0]['accepted_step_comparisons']==2
    assert result['passes'][0]['max_step_replay_error']<1e-14
    assert result['passes'][0]['history_resets']==2
    wrong=copy.deepcopy(events);wrong[3]['payload']['state']=[.26,.74]
    bad=tmp_path/'state.jsonl';write_archive(bad,wrong)
    with pytest.raises(ValueError,match='REPLAY_DISAGREES'):m.replay(bad)
    wrong=copy.deepcopy(events);wrong[3]['payload']['history_before']=7
    bad=tmp_path/'history.jsonl';write_archive(bad,wrong)
    with pytest.raises(ValueError,match='HISTORY_BEFORE'):m.replay(bad)

def test_replay_cannot_accept_missing_archive_footer(tmp_path):
    m=api();path=tmp_path/'complete.jsonl';write_archive(path,replay_fixture())
    bad=tmp_path/'truncated.jsonl';bad.write_text('\n'.join(path.read_text().splitlines()[:-1])+'\n')
    with pytest.raises(ValueError):m.replay(bad)

def test_replay_rejects_broadcastable_short_final_state(tmp_path):
    m=api();events=replay_fixture()[:4];events[0]['trace_counts']=[2]
    p=events[2]['payload'];p.update(gradient=[1.,1.],constraint_jacobian=[[1.,-1.]],
        accepted_step=.5,curvature_pair_accepted=False)
    p=events[3]['payload'];p.update(state=[0.],gradient=[0.,0.],constraint_jacobian=[[1.,-1.]],
        final_state=True,history_before=0,lbfgs_history_size=0,accepted_step=0.,curvature_pair_accepted=False)
    path=tmp_path/'short-final-state.jsonl';write_archive(path,events)
    with pytest.raises(ValueError,match='DIMENSION'):m.replay(path)

def metric_rescue_fixture():
    # SPD quadratic Hessian [[3.6,-.8],[-.8,.4]] plus a decoupled z=0
    # equality coordinate; only x0 bounded. No rank-deficient dummy row.
    # First secant s=(-1,-2), y=(-2,0) gives H=[[.5,1],[1,4.5]].
    # At x=(0,.5), g=(-1,2), ordinary d=(-1.5,-8) leaves the box;
    # the checked tangent QP instead gives d=(0,-5), multiplier=-3.
    events=replay_fixture()
    events[1]['payload']['bounded_variable_count']=1
    for i,(x,g) in enumerate([([1.,2.5],[1.,2.]),([0.,.5],[-1.,2.]),([0.,0.],[-.6,1.8])]):
        p=events[i+2]['payload']
        p.update(state=x+[0.],gradient=g+[0.],constraint_jacobian=[[0.,0.,1.]],constraints=[0.],kkt={'lambda':[0.]},
            active_membership=[0,0,0],accepted_step=[1.,.1,0.][i],
            history_before=i,lbfgs_history_size=i,history_reset_reason='NEW_SOLVE_OR_PARTITION' if i==0 else '')
    events[3]['payload']['metric_tangent']=dict(attempted=True,checked=True,
        selected=True,accepted=True,status='QP_FIRST_ORDER_CHECKED',direction=[0.,-5.,0.],
        working_set=[0],multipliers=[0.,-3.],fallback_attempted=False,trial_count=1,
        equality_residual=0.,stationarity=0.,
        trace=[dict(working_set=[0],action='CHECKED',changed_index=-1,equality_residual=0.,stationarity=0.)])
    events[3]['payload']['objective_trials']=[dict(step=.1,accepted=True,valid=True,
        constraint_l1=0.,reason='ACCEPTED')]
    return events

def test_replay_uses_checked_metric_rescue_not_unselected_qn_direction(tmp_path):
    path=tmp_path/'rescue.jsonl';write_archive(path,metric_rescue_fixture())
    row=api().replay(path)['passes'][0]
    assert row['accepted_step_comparisons']==2
    assert row['metric_accepted_step_comparisons']==1
    assert row['max_step_replay_error']<1e-12

@pytest.mark.parametrize('field,value',[
    ('checked',False),('selected',False),('direction',[0.]),
    ('direction',[0.,None,0.]),('multipliers',[0.,3.]),
    ('working_set',[0,0]),('working_set',[1]),('direction',[0.,-4.,0.]),
    ('attempted',False),('accepted',False),('status','QP_UNKNOWN'),
    ('fallback_attempted',True),('trial_count',0),('trial_count',2),
    ('equality_residual',1.),('stationarity',1.),('trace',[]),
])
def test_replay_rejects_false_or_invalid_metric_certificate(tmp_path,field,value):
    events=metric_rescue_fixture();events[3]['payload']['metric_tangent'][field]=value
    path=tmp_path/'invalid-rescue.jsonl';write_archive(path,events)
    with pytest.raises(ValueError,match='METRIC'):api().replay(path)

def test_replay_rejects_conflicting_rescue_and_fallback_acceptance(tmp_path):
    events=metric_rescue_fixture();events[3]['payload']['used_first_order_fallback']=True
    path=tmp_path/'conflict.jsonl';write_archive(path,events)
    with pytest.raises(ValueError,match='METRIC'):api().replay(path)

def test_replay_rejects_rank_deficient_native_checked_claim(tmp_path):
    events=metric_rescue_fixture()
    for event in events[2:]:event['payload']['constraint_jacobian']=[[0.,0.,0.]]
    path=tmp_path/'rank-deficient.jsonl';write_archive(path,events)
    with pytest.raises(ValueError,match='METRIC.*RANK'):api().replay(path)

@pytest.mark.parametrize('change',['step','accepted_trial','terminal_action'])
def test_replay_rejects_rescue_trial_trace_or_step_mismatch(tmp_path,change):
    events=metric_rescue_fixture();p=events[3]['payload']
    if change=='step':p['accepted_step']=2.;p['objective_trials'][-1]['step']=2.
    if change=='accepted_trial':p['objective_trials'][-1]['accepted']=False
    if change=='terminal_action':p['metric_tangent']['trace'][-1]['action']='RANK_DEFICIENT'
    path=tmp_path/'bad-trace.jsonl';write_archive(path,events)
    with pytest.raises(ValueError,match='METRIC'):api().replay(path)

def test_metric_report_includes_tolerated_face_dual_residual():
    m=api();p=metric_rescue_fixture()[3]['payload']
    p['metric_tangent']['multipliers']=[0.,5e-11]
    pairs=[(np.array([-1.,-2.,0.]),np.array([-2.,0.,0.]))]
    _,residual=m.checked_metric_direction(p,np.array([0.,.5,0.]),
        np.array([-4.-5e-11,2.,0.]),np.array([[0.,0.,1.]]),np.array([0.]),
        1,pairs,1e-10,1e-9)
    assert residual==pytest.approx(5e-11,rel=1e-12,abs=1e-20)
