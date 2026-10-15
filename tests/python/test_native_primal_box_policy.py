import pytest
from test_native_tangent_policy import arguments,solve
from test_native_finite_box_policy import policy as box_policy
from test_finite_box_direction_replay import fixture
from test_slqp_direction_replay import api
from make_it_stand.verification.source_trace_archive import write_archive

def policy():
    return dict(box_policy(),direction_policy='FINITE_BOX_PRIMAL_V2')

def test_native_primal_policy_is_explicit_and_reproducible():
    p,kw=arguments();kw['solver_settings']['slqp'].update(policy())
    try:r=solve(p,kw)
    except ValueError as exc:pytest.fail('missing primal-box policy: '+str(exc))
    assert r['settings']['solver_settings']['slqp']==kw['solver_settings']['slqp']
    assert r==solve(p,kw)
    assert any('finite_box' in t for p in r['solve']['passes'] for t in p['trace'])
    kw['solver_settings']['slqp']['direction_policy']='FINITE_BOX_QP_V1'
    assert r['solve']['input_hash']!=solve(p,kw)['solve']['input_hash']

def test_independent_replay_accepts_explicit_primal_profile(tmp_path):
    events=fixture();events[0]['payload']['settings']['solver_settings']['slqp']['direction_policy']='FINITE_BOX_PRIMAL_V2'
    path=tmp_path/'primal.jsonl';write_archive(path,events)
    assert api().replay(path)['passes'][0]['finite_box_accepted_step_comparisons']==1

def test_primal_disk_roundtrip(tmp_path):
    from test_source_online_trace import run,reconstruct
    p,kw=arguments();controls=kw['solver_settings']['slqp'];controls.update(policy())
    disk=run(tmp_path/'primal',slqp=controls)
    assert disk.diagnostics==run(slqp=controls).diagnostics==reconstruct(disk.iter_diagnostics())

def test_primal_rank_audit_uses_scaled_system_but_checks_raw_constraints(tmp_path):
    events=fixture();events[0]['payload']['settings']['solver_settings']['slqp']['direction_policy']='FINITE_BOX_PRIMAL_V2'
    events[1]['payload']['bounded_variable_count']=3
    for i,e in enumerate(events[2:]):
        e['payload'].update(state=[.5,.5,.5+.05*i],gradient=[0.,0.,-.1+.05*i],
            constraints=[0.,0.],constraint_jacobian=[[1e-8,0.,0.],[0.,1.,0.]],kkt={'lambda':[0.,0.]})
    q=events[2]['payload']['finite_box'];q.update(direction=[0.,0.,.1],multipliers=[0.,0.],working_set=[])
    q['trace'][-1]['working_set']=[]
    path=tmp_path/'scaled.jsonl';write_archive(path,events)
    assert api().replay(path)['passes'][0]['finite_box_accepted_step_comparisons']==1
    # Scaling may not hide a nonzero raw equality residual.
    events[2]['payload']['constraints'][1]=1e-4
    path=tmp_path/'bad-raw.jsonl';write_archive(path,events)
    with pytest.raises(ValueError):api().replay(path)
