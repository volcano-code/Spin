"""Timing accounts for nested stages without treating missing work as free."""
import importlib

import pytest


def event(kind, stage, start, elapsed=None):
    row = dict(kind=kind, stage=stage, monotonic_ns=start)
    if elapsed is not None:
        row['elapsed_ns'] = elapsed
    return row


def test_nested_stages_have_exclusive_costs_and_explicit_unclassified_time():
    api = importlib.import_module('make_it_stand.verification.faster_timing')
    rows = [event('STAGE_STARTED', 'native_newton_solver', 0),
            event('STAGE_STARTED', 'production_geometry_evaluation', 10),
            event('STAGE_STARTED', 'boolean', 20),
            event('STAGE_RETURNED', 'boolean', 40, 20),
            event('STAGE_STARTED', 'intersection_motion', 40),
            event('STAGE_RETURNED', 'intersection_motion', 50, 10),
            event('STAGE_RETURNED', 'production_geometry_evaluation', 60, 50),
            event('STAGE_STARTED', 'evaluation_checkpoint', 60),
            event('STAGE_RETURNED', 'evaluation_checkpoint', 80, 20),
            event('STAGE_RETURNED', 'native_newton_solver', 100, 100)]
    report = api.summarize_solver_timing(rows)
    assert report['solver_completed']
    assert report['solver_wall_ns'] == 100
    assert report['measured_algorithm_ns'] == 30
    assert report['recording_ns'] == 20
    assert report['unclassified_ns'] == 50
    assert sum(report['exclusive_stage_ns'].values()) == 100
    assert report['paper_comparable_seconds'] is None


def test_interrupted_solver_is_never_given_a_time_to_convergence():
    api = importlib.import_module('make_it_stand.verification.faster_timing')
    rows = [event('STAGE_STARTED', 'native_newton_solver', 0),
            event('STAGE_STARTED', 'boolean', 10),
            event('STAGE_RETURNED', 'boolean', 30, 20)]
    report = api.summarize_solver_timing(rows)
    assert not report['solver_completed']
    assert report['solver_wall_ns'] is None
    assert report['measured_algorithm_ns'] == 20
    assert report['unclassified_ns'] is None
    assert report['paper_comparable_seconds'] is None


def test_crossed_stage_records_are_rejected():
    api = importlib.import_module('make_it_stand.verification.faster_timing')
    with pytest.raises(ValueError, match='STAGE_NESTING'):
        api.summarize_solver_timing([
            event('STAGE_STARTED', 'native_newton_solver', 0),
            event('STAGE_STARTED', 'boolean', 10),
            event('STAGE_RETURNED', 'native_newton_solver', 30, 30)])


def test_offline_trial_profile_is_rejected_before_loading_large_input(tmp_path, monkeypatch):
    from pathlib import Path
    root = Path(__file__).resolve().parents[2]
    monkeypatch.syspath_prepend(str(root/'experiments/faster_core_repair_20261012'))
    worker = importlib.import_module('run_noodles_terminal_recovery')
    output = tmp_path/'uncreated'
    with pytest.raises(ValueError, match='OFFLINE_RECOVERY_REQUIRES_INDEPENDENT_BOUNDARY'):
        worker.main(output, tmp_path/'missing-source', classified_boundary=True)
    assert not output.exists()
