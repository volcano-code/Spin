"""Resource rollback must be auditable without becoming scientific success."""
import copy
import hashlib
import importlib
import json
import subprocess
import sys
from pathlib import Path

import pytest
from test_stage_audit import ROOT, fixture, failed_pass_fixture


def auditor():
    # Missing implementation is an explicit assertion failure, not collection failure.
    assert importlib.util.find_spec('stage_audit_v2') is not None, 'versioned rollback auditor missing'
    return importlib.import_module('stage_audit_v2')


def rollback_fixture():
    header, stage, terminal = failed_pass_fixture()
    stage['status'] = 'CONVERGED'
    header['solve'].update(status='RESOURCE_LIMIT_REACHED_UNKNOWN',
        diagnostic='RESOURCE_LIMIT: adaptive source leaves',
        terminal_relaxation_kkt_applies=True, binary_cuboids=None)
    header['settings']['solver_settings'].update(maximum_optimization_level=2,
        maximum_leaf_count=8, slqp={'projected_gradient_tolerance': 1e-7,
                                  'equality_tolerance': 1e-8})
    # Two leaves, one splittable: conservative pre-merge bound 2+7=9 exceeds 8.
    return header, stage, terminal


def test_resource_rollback_preserves_converged_kkt_but_never_admits_material():
    header, stage, terminal = rollback_fixture()
    original = copy.deepcopy((header, stage, terminal))
    result = auditor().audit_components(header, [(stage, terminal)])
    assert (header, stage, terminal) == original
    item = result['passes'][0]
    assert item['solver_status'] == 'CONVERGED'
    assert item['transition']['state'] == 'ROLLED_BACK_RESOURCE_LIMIT'
    assert item['transition']['resource_rollback']['one_split_upper_bound'] == 9
    assert item['transition']['resource_rollback']['actual_split_size'] is None
    assert item['kkt_evidence']['terminal_partition'] == 'APPLICABLE_UNCHANGED_TERMINAL_STATE'
    assert item['split_merge_moment_delta'] == [0.0]*10
    assert result['terminal_rounding']['available'] is False
    assert result['geometry_export_stage'] == 'NOT_AVAILABLE'
    assert result['scientific_status'] == 'RESOURCE_LIMIT_REACHED_UNKNOWN'
    assert result['paper_acceptance'] == result['sealed_artifact'] == 'NOT_GRANTED'


@pytest.mark.parametrize('mutation', [
    'wrong_status', 'wrong_cause', 'partition', 'fill', 'terminal_state', 'terminal_phi',
    'final_fill', 'final_phi', 'final_partition', 'partial_hash', 'changes', 'removed',
    'resolution', 'binary_fill', 'binary_moments', 'binary_check', 'binary_cuboids',
    'rounding_delta', 'geometry', 'invalid_kkt', 'excess_kkt', 'negative_kkt',
    'nonfinite_kkt', 'applicability', 'impossible_cap', 'invalid_leaf', 'duplicate_leaf',
    'bad_index', 'missing_binary_fill', 'invalid_tolerance'])
def test_rejects_forged_or_inconsistent_resource_rollback(mutation):
    api = auditor()
    header, stage, terminal = rollback_fixture()
    solve = header['solve']; cfg = header['settings']['solver_settings']
    if mutation == 'wrong_status': solve['status'] = 'BINARY_MOMENT_VALIDATED'
    elif mutation == 'wrong_cause': solve['diagnostic'] = 'RESOURCE_LIMIT: wall time'
    elif mutation == 'partition': stage['after'][0]['coordinates'][0] = 1
    elif mutation == 'fill': stage['output_beta'][0] += 1e-13
    elif mutation == 'terminal_state': terminal['state'][0] += 1e-13
    elif mutation == 'terminal_phi': terminal['phi'] += 1e-13
    elif mutation == 'final_fill': solve['relaxed_fill'] = [0.2500000000001, .75]
    elif mutation == 'final_phi': solve['final_phi'] += 1e-13
    elif mutation == 'final_partition': solve['final_partition'] = []
    elif mutation == 'partial_hash': stage['split']['before_hash'] = stage['before_hash']
    elif mutation == 'changes': stage['merge']['changes'] = [{}]
    elif mutation == 'removed': stage['split']['removed_after'][0] = 1e-13
    elif mutation == 'resolution': stage['paper_resolution_complete'] = True
    elif mutation in ('binary_fill', 'rounding_delta'): solve[mutation] = [0.0]
    elif mutation in ('binary_moments', 'binary_cuboids'): solve[mutation] = []
    elif mutation == 'binary_check': solve['binary_check'] = {'passed': True}
    elif mutation == 'geometry': header['geometry']['status'] = 'EXPORTED'
    elif mutation == 'invalid_kkt': stage['final_kkt']['valid'] = False
    elif mutation == 'excess_kkt': stage['final_kkt']['full'] = 1.01e-7
    elif mutation == 'negative_kkt': stage['final_kkt']['stationarity'] = -1e-8
    elif mutation == 'nonfinite_kkt': stage['final_kkt']['full'] = float('nan')
    elif mutation == 'applicability': solve['terminal_relaxation_kkt_applies'] = False
    elif mutation == 'impossible_cap': cfg['maximum_leaf_count'] = 9
    elif mutation == 'invalid_leaf': stage['before'][0]['coordinates'] = [2, 0, 0]
    elif mutation == 'duplicate_leaf':
        for cell in (stage['before'][1], stage['after'][1]):
            cell.update(level=1, coordinates=[0, 0, 0])
    elif mutation == 'bad_index': stage['index'] = False
    elif mutation == 'missing_binary_fill': del solve['binary_fill']
    elif mutation == 'invalid_tolerance': cfg['slqp']['projected_gradient_tolerance'] = float('inf')
    with pytest.raises(ValueError):
        api.audit_components(header, [(stage, terminal)])


def test_rollback_cannot_be_followed_by_another_stage():
    header, stage, terminal = rollback_fixture()
    following = copy.deepcopy(stage); following['index'] = 1
    with pytest.raises(ValueError, match='ROLLBACK_NOT_TERMINAL'):
        auditor().audit_components(header, [(stage, terminal), (following, terminal)])


def test_equality_tolerance_is_not_replaced_by_projected_gradient_tolerance():
    header, stage, terminal = rollback_fixture()
    header['settings']['solver_settings']['slqp']['equality_tolerance'] = 1e-6
    stage['final_kkt'].update(equality_primal=5e-7, full=5e-7)
    result = auditor().audit_components(header, [(stage, terminal)])
    assert result['passes'][0]['transition']['state'] == 'ROLLED_BACK_RESOURCE_LIMIT'
    assert result['paper_acceptance'] == 'NOT_GRANTED'


def test_equality_residual_above_its_own_tolerance_is_rejected():
    header, stage, terminal = rollback_fixture()
    stage['final_kkt']['equality_primal'] = 5e-8
    with pytest.raises(ValueError):
        auditor().audit_components(header, [(stage, terminal)])


def test_conflicting_geometry_state_cannot_hide_behind_not_attempted_status():
    header, stage, terminal = rollback_fixture()
    header['geometry']['state'] = 'ATTEMPTED'
    with pytest.raises(ValueError):
        auditor().audit_components(header, [(stage, terminal)])


@pytest.mark.parametrize('kind', ['executed', 'solver_failed'])
def test_existing_transition_algebra_is_preserved(kind):
    from stage_audit import audit_components
    header, stage, terminal = fixture() if kind == 'executed' else failed_pass_fixture()
    previous = audit_components(header, [(stage, terminal)])
    current = auditor().audit_components(header, [(stage, terminal)])
    for name in ('passes', 'terminal_rounding', 'source_digests', 'geometry_export_stage'):
        assert current[name] == previous[name]


def test_archive_cli_checks_footer_and_refuses_overwrite(tmp_path):
    from make_it_stand.verification.source_trace_archive import write_archive
    api = auditor()
    header, stage, terminal = rollback_fixture()
    events = [{'kind': 'header', 'trace_counts': [1], 'payload': header},
              {'kind': 'pass', 'pass_index': 0, 'payload': stage},
              {'kind': 'iteration', 'pass_index': 0, 'trace_index': 0, 'payload': terminal}]
    archive = tmp_path/'trace.jsonl'; output = tmp_path/'audit.json'
    write_archive(archive, events)
    original = archive.read_bytes()
    archive.write_bytes(original.rsplit(b'\n', 2)[0]+b'\n')
    with pytest.raises(ValueError): api.audit_archive(archive)
    archive.write_bytes(original)
    command = [sys.executable, str(Path(api.__file__)), str(archive), str(output)]
    run = subprocess.run(command, capture_output=True, text=True)
    assert run.returncode == 0, run.stderr
    result = json.loads(output.read_text())
    assert result['schema'] == 'spin2014-stage-audit-v2'
    assert result['input_archive_sha256'] == hashlib.sha256(original).hexdigest()
    assert result['archive_integrity'] == 'VERIFIED_COMPLETE_SOURCE_TRACE_ARCHIVE_V1'
    assert subprocess.run(command, capture_output=True).returncode != 0


def test_native_event_audit_matches_archive_audit_without_materializing_jsonl(tmp_path):
    from make_it_stand.verification.source_trace_archive import write_archive
    api = auditor()
    header, stage, terminal = rollback_fixture()
    events = [{'kind': 'header', 'trace_counts': [1], 'payload': header},
              {'kind': 'pass', 'pass_index': 0, 'payload': stage},
              {'kind': 'iteration', 'pass_index': 0, 'trace_index': 0, 'payload': terminal}]
    archive = tmp_path/'native-trace.jsonl'
    write_archive(archive, events)
    expected = api.audit_archive(archive)
    native_line = lambda event: (json.dumps(event, sort_keys=True,
        separators=(',', ':'), allow_nan=False) + '\n').encode()
    actual = api.audit_native_diagnostics(
        iter(native_line(event) for event in events), expected_pass_count=1,
        expected_record_count=1)
    assert actual['passes'] == expected['passes']
    assert actual['terminal_rounding'] == expected['terminal_rounding']
    assert actual['source_digests'] == expected['source_digests']
    assert actual['archive_integrity'] == 'VERIFIED_NATIVE_TRACE_SEQUENCE_AND_RECORD_CHECKSUMS'
    assert 'input_archive_sha256' not in actual


def test_compact_iteration_markers_preserve_full_stage_audit_result(tmp_path):
    from make_it_stand.verification.source_trace_archive import write_archive
    api = auditor(); header, stage, terminal = rollback_fixture()
    header_event = {'kind': 'header', 'trace_counts': [2], 'payload': header}
    full_events = [header_event, {'kind': 'pass', 'pass_index': 0, 'payload': stage},
                   {'kind': 'iteration', 'pass_index': 0, 'trace_index': 0,
                    'payload': {'iteration': 0, 'gradient': [1.0]}},
                   {'kind': 'iteration', 'pass_index': 0, 'trace_index': 1, 'payload': terminal}]
    archive = tmp_path/'full-trace.jsonl'; write_archive(archive, full_events)
    expected = api.audit_archive(archive)
    compact_events = [header_event, full_events[1],
        {'kind': 'iteration', 'pass_index': 0, 'trace_index': 0, 'payload': {}},
        full_events[3]]
    encoded = ((json.dumps(event, sort_keys=True, separators=(',', ':'), allow_nan=False) + '\n').encode()
               for event in compact_events)
    actual = api.audit_native_diagnostics(encoded, expected_pass_count=1, expected_record_count=2)
    assert actual['passes'] == expected['passes']
    assert actual['terminal_rounding'] == expected['terminal_rounding']
    assert actual['source_digests'] == expected['source_digests']


@pytest.mark.parametrize('corruption,reason', [
    ('iteration_sequence', 'ITERATION_SEQUENCE'),
    ('missing_terminal', 'INCOMPLETE_ARCHIVE'),
    ('truncated_json', 'TRUNCATED_NATIVE_EVENT'),
])
def test_native_event_audit_rejects_incomplete_or_misordered_stream(corruption, reason):
    api = auditor(); header, stage, terminal = rollback_fixture()
    events = [{'kind': 'header', 'trace_counts': [1], 'payload': header},
              {'kind': 'pass', 'pass_index': 0, 'payload': stage},
              {'kind': 'iteration', 'pass_index': 0, 'trace_index': 0, 'payload': terminal}]
    if corruption == 'iteration_sequence': events[-1]['trace_index'] = 1
    elif corruption == 'missing_terminal': events.pop()
    raw = [json.dumps(event, sort_keys=True, separators=(',', ':'), allow_nan=False).encode() + b'\n'
           for event in events]
    if corruption == 'truncated_json': raw[-1] = raw[-1][:-2]
    with pytest.raises(ValueError, match=reason):
        api.audit_native_diagnostics(iter(raw), expected_pass_count=1, expected_record_count=1)
