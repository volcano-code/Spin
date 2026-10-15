import copy
import pytest
from test_native_tangent_policy import arguments,solve
from test_native_primal_box_policy import policy

def test_problem_sized_capacity_roundtrip_hash_and_disk(tmp_path):
    from test_source_online_trace import run,reconstruct
    p,kw=arguments();controls=kw['solver_settings']['slqp'];controls.update(policy())
    controls['finite_box_maximum_faces']='BOUNDED_VARIABLES'
    result=solve(p,kw)
    assert result['settings']['solver_settings']['slqp']==controls
    assert result==solve(p,kw)
    legacy=copy.deepcopy(kw);legacy['solver_settings']['slqp']['finite_box_maximum_faces']=256
    assert result['solve']['input_hash']!=solve(p,legacy)['solve']['input_hash']
    disk=run(tmp_path/'automatic',slqp=controls)
    assert disk.diagnostics==run(slqp=controls).diagnostics==reconstruct(disk.iter_diagnostics())

@pytest.mark.parametrize('value',['AUTOMATIC',-1,0,True,1.5])
def test_invalid_external_capacity_rejected(value):
    p,kw=arguments();kw['solver_settings']['slqp'].update(policy())
    kw['solver_settings']['slqp']['finite_box_maximum_faces']=value
    with pytest.raises((ValueError,TypeError)):solve(p,kw)

def test_automatic_mode_is_not_implicitly_enabled_for_legacy_solver():
    p,kw=arguments();kw['solver_settings']['slqp'].update(policy())
    kw['solver_settings']['slqp'].update(direction_policy='FINITE_BOX_QP_V1',finite_box_maximum_faces='BOUNDED_VARIABLES')
    with pytest.raises((ValueError,TypeError)):solve(p,kw)

def test_auto_capacity_survives_stage_resume_and_rejects_policy_swap(tmp_path):
    import _make_it_stand as native
    from test_source_spin_material_bridge import box,problem
    p,kw=arguments();kw['solver_settings']['slqp'].update(policy())
    # Analytically stationary symmetric UNIT fixture for checkpoint transport,
    # not an optimization-success or paper-case benchmark. The generic asymmetric
    # fixture can stop at a numerical UNKNOWN before a checkpoint exists.
    p=problem(outer=box((-1,-1,-1),(1,1,1)),design=box((-.5,-.5,-.5),(.5,.5,.5)))
    kw['source_settings']['root_origin']=(-1.,-1.,-1.)
    kw['solver_settings']['slqp']['finite_box_maximum_faces']='BOUNDED_VARIABLES'
    arrays=(*p.original.arrays(),*p.design_regions[0].mesh.arrays())
    baseline=native.solve_source_spin2014_material(*arrays,**kw)
    staged=native.solve_source_spin2014_material(*arrays,**kw,maximum_new_passes=1,
        iteration_store_directory=str(tmp_path/'first'))
    assert staged.checkpoint is not None and not staged.has_geometry
    swapped=copy.deepcopy(kw);swapped['solver_settings']['slqp']['finite_box_maximum_faces']=256
    with pytest.raises((ValueError,RuntimeError),match='checkpoint'):
        native.solve_source_spin2014_material(*arrays,**swapped,maximum_new_passes=1,
            iteration_store_directory=str(tmp_path/'rejected'),resume_checkpoint=staged.checkpoint)
    assert not (tmp_path/'rejected').exists()
    for index in range(8):
        if staged.checkpoint is None:break
        staged=native.solve_source_spin2014_material(*arrays,**kw,maximum_new_passes=1,
            iteration_store_directory=str(tmp_path/f'resume-{index}'),resume_checkpoint=staged.checkpoint)
    assert staged.checkpoint is None
    assert staged.diagnostics==baseline.diagnostics
