from dataclasses import asdict
import copy
import pytest
import _make_it_stand as native
from test_source_spin_material_bridge import problem,settings
from make_it_stand.verification import source_spin_material as bridge

def arguments():
    p=problem();cfg=settings(bridge)
    kw={k:asdict(v) for k,v in cfg.items()}
    kw.update(density=1.,minimum_outer_wall=.5)
    return p,kw

def policy():
    return dict(direction_policy='METRIC_TANGENT_RESCUE_V1',metric_tangent_maximum_faces=256,
        metric_tangent_maximum_iterations=1028,metric_tangent_residual_tolerance=1e-9)

def solve(p,kw):
    return native.solve_source_spin2014_material(*p.original.arrays(),*p.design_regions[0].mesh.arrays(),**kw).diagnostics

def test_native_opt_in_is_roundtrippable_and_hash_bound():
    p,kw=arguments();kw['solver_settings']['slqp'].update(policy())
    try:r=solve(p,kw)
    except ValueError as e:pytest.fail('missing explicit native tangent policy: '+str(e))
    assert r['settings']['solver_settings']['slqp']==kw['solver_settings']['slqp']
    assert r==solve(p,kw)
    old=copy.deepcopy(kw)
    for k in policy():old['solver_settings']['slqp'].pop(k)
    baseline=solve(p,old)
    assert baseline['solve']['input_hash']!=r['solve']['input_hash']
    assert 'direction_policy' not in baseline['settings']['solver_settings']['slqp']
    for a in baseline['solve']['passes']:
        assert all('metric_tangent' not in t for t in a['trace'])
    changed=copy.deepcopy(kw);changed['solver_settings']['slqp']['metric_tangent_maximum_faces']=128
    assert solve(p,changed)['solve']['input_hash']!=r['solve']['input_hash']

@pytest.mark.parametrize('change',[{'direction_policy':'FORGED'},
    {'metric_tangent_maximum_faces':0},{'metric_tangent_residual_tolerance':float('nan')}])
def test_invalid_policy_controls_rejected(change):
    p,kw=arguments();kw['solver_settings']['slqp'].update(policy());kw['solver_settings']['slqp'].update(change)
    with pytest.raises(ValueError):solve(p,kw)
