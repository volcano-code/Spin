"""Evidence failure paths must not disappear from the experiment ledger."""
import copy
import sys
from pathlib import Path
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'experiments/teapot_spin2014'))
import summarize_resolution as report


def test_no_pass_resource_result_has_unavailable_metrics():
    row=report.pass_metrics([])
    assert row['last_kkt'] is None
    assert row['maximum_solved_leaf_level'] is None
    assert row['largest_solved_partition'] is None


def test_empty_partition_unknown_is_preserved():
    row=report.pass_metrics([dict(before=[],status='UNKNOWN',trace=[])])
    assert row['maximum_solved_leaf_level'] is None
    assert row['largest_solved_partition']==0
    assert row['last_solve_status']=='UNKNOWN'


def test_missing_first_repeat_is_not_false_determinism():
    row=report.repeat_comparison({'status':'NOT_COMPLETED'},{'native_record_sha256':'abc'})
    assert row['repeat_native_record_bytes_equal'] is None
    assert row['repeat_comparison_status']=='UNAVAILABLE'
    assert report.repeat_comparison({'native_record_sha256':'a'},
        {'native_record_sha256':'b'})['repeat_native_record_bytes_equal'] is False


def test_required_output_cannot_be_omitted_from_manifest(tmp_path):
    complete={str((tmp_path/name).resolve()):'digest' for name in
              ('native-run.json','config.json','summary.json','source.json')}
    report.verify_required_outputs(tmp_path,complete)
    for path in complete:
        incomplete={k:v for k,v in complete.items() if k!=path}
        with pytest.raises(ValueError,match='MISSING_EVIDENCE_COVERAGE'):
            report.verify_required_outputs(tmp_path,incomplete)


def test_actual_settings_must_match_every_declared_numeric_control():
    cfg=dict(source_settings={'maximum_depth':6},solver_settings={
        'termination_policy':'PAPER_POST_TRANSITION_V1','binary_constraint_tolerance':1e-6,
        'slqp':{'maximum_iterations':960,'projected_gradient_tolerance':1e-7}},
        boundary_settings={'maximum_input_cuboids':12000},density=1.,minimum_outer_wall=.015)
    report.verify_native_settings(cfg,cfg)
    for key in cfg:
        bad=copy.deepcopy(cfg);bad[key]=None
        with pytest.raises(ValueError,match='NATIVE_SETTINGS_MISMATCH'):
            report.verify_native_settings(bad,cfg)
    bad=copy.deepcopy(cfg);bad['solver_settings']['slqp']['maximum_iterations']=240
    with pytest.raises(ValueError,match='NATIVE_SETTINGS_MISMATCH'):
        report.verify_native_settings(bad,cfg)
