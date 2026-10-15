"""Stage resume must reach the same independently sealed material, not trust JSON."""
from dataclasses import replace
from pathlib import Path
import numpy as np
import pytest
from test_source_spin_material_bridge import problem, settings
from make_it_stand.verification import source_spin_material as bridge


def test_resumed_native_stages_keep_geometry_trace_and_seal(tmp_path):
    p=problem();cfg=settings(bridge)
    baseline=bridge.solve_and_seal_source_spin2014(p,**cfg)
    first=bridge.solve_and_seal_source_spin2014(p,**cfg,maximum_new_passes=1,
        iteration_store_directory=str(tmp_path/'first'))
    assert first.status=='STAGE_BOUNDARY_PAUSED_NOT_VALIDATED'
    assert first.artifact is None and first.material is None
    assert isinstance(first.checkpoint,bytes) and first.checkpoint
    checkpoint=first.checkpoint
    for index in range(8):
        result=bridge.solve_and_seal_source_spin2014(p,**cfg,maximum_new_passes=1,
            iteration_store_directory=str(tmp_path/f'resume-{index}'),resume_checkpoint=checkpoint)
        if result.checkpoint is None:break
        assert result.artifact is None
        checkpoint=result.checkpoint
    assert result.status=='BINARY_MOMENT_VALIDATED'
    assert result.artifact is not None
    assert result.scientific_digest==baseline.scientific_digest
    assert result.artifact.scientific_digest==baseline.artifact.scientific_digest
    for name in ['outer','material','cavity']:
        for a,b in zip(getattr(result,name),getattr(baseline,name)):
            assert a.tobytes()==b.tobytes()
    assert result.identity['native_run']['solve']==baseline.identity['native_run']['solve']


def test_changed_request_cannot_restore_or_export(tmp_path):
    p=problem();cfg=settings(bridge)
    first=bridge.solve_and_seal_source_spin2014(p,**cfg,maximum_new_passes=1,
        iteration_store_directory=str(tmp_path/'first'))
    cfg['solver_settings']=replace(cfg['solver_settings'],center_weight=2.)
    with pytest.raises((ValueError,RuntimeError)):
        bridge.solve_and_seal_source_spin2014(p,**cfg,maximum_new_passes=1,
            iteration_store_directory=str(tmp_path/'rejected'),resume_checkpoint=first.checkpoint)
    assert not (tmp_path/'rejected').exists()


def test_invalid_checkpoint_cannot_restore_or_export(tmp_path):
    with pytest.raises((ValueError,RuntimeError)):
        bridge.solve_and_seal_source_spin2014(problem(),**settings(bridge),maximum_new_passes=1,
            iteration_store_directory=str(tmp_path/'rejected'),resume_checkpoint=b'forged')
    assert not (tmp_path/'rejected').exists()


@pytest.mark.parametrize('changed',['wall','boundary_budget'])
def test_material_checkpoint_binds_entire_export_request(tmp_path,changed):
    p=problem();cfg=settings(bridge)
    first=bridge.solve_and_seal_source_spin2014(p,**cfg,maximum_new_passes=1,
        iteration_store_directory=str(tmp_path/'first'))
    if changed=='wall':p=problem(wall=.4)
    else:cfg['boundary_settings']=replace(cfg['boundary_settings'],maximum_output_vertices=900000)
    with pytest.raises((ValueError,RuntimeError),match='material checkpoint'):
        bridge.solve_and_seal_source_spin2014(p,**cfg,maximum_new_passes=1,
            iteration_store_directory=str(tmp_path/'rejected'),resume_checkpoint=first.checkpoint)
    assert not (tmp_path/'rejected').exists()


@pytest.mark.parametrize('extra',[
    {'maximum_new_passes':True}, {'maximum_new_passes':0},
    {'maximum_new_passes':1}, {'resume_checkpoint':b'fake'},
])
def test_invalid_stage_controls_fail_before_native_output(extra):
    with pytest.raises(ValueError):
        bridge.solve_and_seal_source_spin2014(problem(),**settings(bridge),**extra)


def test_changed_trace_budget_is_not_a_fresh_resume_budget(tmp_path):
    p=problem();cfg=settings(bridge)
    first=bridge.solve_and_seal_source_spin2014(p,**cfg,maximum_new_passes=1,
        iteration_store_directory=str(tmp_path/'first'))
    with pytest.raises((ValueError,RuntimeError)):
        bridge.solve_and_seal_source_spin2014(p,**cfg,maximum_new_passes=1,
            iteration_store_directory=str(tmp_path/'rejected'),resume_checkpoint=first.checkpoint,
            iteration_store_max_bytes=2*1073741824)
    assert not (tmp_path/'rejected').exists()


def test_missing_historical_trace_refuses_material_resume(tmp_path):
    p=problem();cfg=settings(bridge)
    first=bridge.solve_and_seal_source_spin2014(p,**cfg,maximum_new_passes=1,
        iteration_store_directory=str(tmp_path/'first'))
    trace=next((tmp_path/'first').glob('*.slqptrace'))
    trace.rename(trace.with_suffix('.preserved'))
    with pytest.raises((ValueError,RuntimeError)):
        bridge.solve_and_seal_source_spin2014(p,**cfg,maximum_new_passes=1,
            iteration_store_directory=str(tmp_path/'rejected'),resume_checkpoint=first.checkpoint)
    assert not (tmp_path/'rejected').exists()


def test_spliced_material_envelope_is_detected(tmp_path):
    p=problem();cfg=settings(bridge)
    first=bridge.solve_and_seal_source_spin2014(p,**cfg,maximum_new_passes=1,
        iteration_store_directory=str(tmp_path/'first'))
    cfg['boundary_settings']=replace(cfg['boundary_settings'],maximum_output_vertices=900000)
    second=bridge.solve_and_seal_source_spin2014(p,**cfg,maximum_new_passes=1,
        iteration_store_directory=str(tmp_path/'second'))
    spliced=second.checkpoint[:72]+first.checkpoint[72:]
    with pytest.raises((ValueError,RuntimeError),match='material checkpoint'):
        bridge.solve_and_seal_source_spin2014(p,**cfg,maximum_new_passes=1,
            iteration_store_directory=str(tmp_path/'rejected'),resume_checkpoint=spliced)
    assert not (tmp_path/'rejected').exists()
