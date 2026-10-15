"""Fail-closed, full-record storage comparisons, not paper admission tests."""
import importlib.util
from pathlib import Path
import json
import shutil
import pytest
from make_it_stand.verification.source_trace_archive import write_archive


def api():
    path = Path(__file__).with_name('online_compare.py')
    assert path.is_file(), 'online storage comparison checker is missing'
    spec = importlib.util.spec_from_file_location('online_compare', path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def events(value=.5):
    return [
        {'kind': 'header', 'trace_counts': [1], 'payload': {'solve': {
            'status': 'BINARY_VALIDATION_FAILED', 'trace_hash': 'same-native-hash'}}},
        {'kind': 'pass', 'pass_index': 0, 'payload': {'index': 0}},
        {'kind': 'iteration', 'pass_index': 0, 'trace_index': 0,
         'payload': {'state': [value], 'final_state': True}},
    ]


def test_comparison_checks_records_even_when_native_hash_claim_matches(tmp_path):
    left, right = tmp_path/'left.jsonl', tmp_path/'right.jsonl'
    write_archive(left, events()); write_archive(right, events(.6))
    with pytest.raises(ValueError, match='RECORD_MISMATCH'):
        api().compare_archives(left, right)


def test_signed_zero_is_not_silently_normalized(tmp_path):
    left, right = tmp_path/'left.jsonl', tmp_path/'right.jsonl'
    write_archive(left, events(-0.0)); write_archive(right, events(0.0))
    with pytest.raises(ValueError, match='RECORD_MISMATCH'):
        api().compare_archives(left, right)


def test_equal_failed_solves_are_only_storage_equivalent(tmp_path):
    left, right = tmp_path/'left.jsonl', tmp_path/'right.jsonl'
    write_archive(left, events()); write_archive(right, events())
    result = api().compare_archives(left, right)
    assert result['records'] == 3
    assert result['archive_equivalence'] == 'VERIFIED'
    assert 'storage_equivalence' not in result
    assert result['numerical_status'] == 'BINARY_VALIDATION_FAILED'
    assert result['paper_reproduction'] == 'NOT_EVALUATED'


def test_matching_prefix_without_footer_is_rejected(tmp_path):
    left, right = tmp_path/'left.jsonl', tmp_path/'right.jsonl'
    write_archive(left, events()); write_archive(right, events())
    raw = right.read_bytes(); right.write_bytes(raw.rsplit(b'\n', 2)[0] + b'\n')
    with pytest.raises(ValueError, match='MISSING_COMPLETION_FOOTER'):
        api().compare_archives(left, right)


def test_input_change_cannot_be_hidden_by_matching_output(tmp_path):
    path = tmp_path/'input.obj'; path.write_bytes(b'old input')
    mod = api(); observed = {str(path): mod.sha256(path)}
    path.write_bytes(b'new input')
    with pytest.raises(ValueError, match='STALE_RESULT'):
        mod.verify_files(observed)


def test_empty_and_identical_path_cannot_count_as_comparison(tmp_path):
    left, right = tmp_path/'left.jsonl', tmp_path/'right.jsonl'
    left.touch(); right.touch()
    with pytest.raises(ValueError, match='MISSING_COMPLETION_FOOTER'):
        api().compare_archives(left, right)
    with pytest.raises(ValueError, match='DISTINCT_RUN_ARCHIVES_REQUIRED'):
        api().compare_archives(left, left)


def fixture_inputs(tmp_path):
    root = tmp_path/'source'; root.mkdir()
    for directory, filename in [('src', 'solver.cpp'), ('include', 'solver.h'),
                                ('bindings', 'module.cpp'), ('python', 'helper.py')]:
        folder = root/directory; folder.mkdir(); (folder/filename).write_text('source v1')
    (root/'CMakeLists.txt').write_text('project(test)')
    inputs = tmp_path/'input'; inputs.mkdir()
    for name in ('outer.obj', 'design.obj'): (inputs/name).write_text('geometry')
    (inputs/'manifest.json').write_text('{}')
    config = tmp_path/'config.json'; config.write_text(json.dumps({'input': str(inputs)}))
    build = tmp_path/'build'; (build/'python').mkdir(parents=True)
    (build/'CMakeCache.txt').write_text('build controls')
    module = build/'python'/'native.so'; module.write_bytes(b'native binary A')
    return root, config, module


def saved_run(mod, root, config, module, output, storage='memory', value=.5,
              storage_verification='COMPLETE_ARCHIVE_READBACK'):
    output.mkdir()
    p = mod.make_provenance(config, module, output, storage, 1000000, source_root=root)
    mod.write_new(output/'provenance.json', p)
    e = events(value)
    e[0]['payload']['solve'].update(input_hash='input', output_hash='output')
    facts = write_archive(output/'native-run.jsonl', e)
    mod.write_new(output/'header.json', e[0])
    mod.write_new(output/'summary.json', dict(
        run_id=p['run_id'], output_directory=str(output.resolve()), storage=storage,
        provenance_sha256=mod.sha256(output/'provenance.json'),
        storage_verification=storage_verification,
        numerical_status='BINARY_VALIDATION_FAILED', input_hash='input',
        trace_hash='same-native-hash', output_hash='output', trace_counts=[1], archive=facts,
        native_storage=({'mode': 'DISK_ITERATIONS_V1', 'resident_pass_count': 1, 'record_count': 1}
                        if storage == 'disk' else {'mode': 'MEMORY', 'resident_pass_count': 1, 'record_count': 1})))
    if storage == 'disk':
        (output/'iteration-store').mkdir()
        (output/'iteration-store'/'pass-000000.slqptrace').write_bytes(b'durable iteration file')
    mod.write_new(output/'outputs.json', mod.output_hashes(output))
    return output


@pytest.mark.parametrize('kind,modes,binary_change', [
    ('cross-mode', ('memory', 'disk'), False),
    ('deterministic-repeat', ('disk', 'disk'), False),
    ('legacy-binary', ('memory', 'memory'), True)])
def test_run_comparison_requires_declared_valid_kind(tmp_path, kind, modes, binary_change):
    mod = api(); root, config, module = fixture_inputs(tmp_path)
    left = saved_run(mod, root, config, module, tmp_path/'left', modes[0])
    second_module = module
    if binary_change:
        other = tmp_path/'other-build'; (other/'python').mkdir(parents=True)
        (other/'CMakeCache.txt').write_text('old build controls')
        second_module = other/'python'/'native.so'; second_module.write_bytes(b'native binary B')
    right = saved_run(mod, root, config, second_module, tmp_path/'right', modes[1])
    result = mod.compare_runs(left, right, kind=kind)
    assert result['comparison_kind'] == kind
    assert result['run_equivalence'] == 'VERIFIED'
    assert result['archive_equivalence'] == 'VERIFIED'
    assert result['numerical_status'] == 'BINARY_VALIDATION_FAILED'


def test_same_mode_and_copied_run_cannot_be_cross_mode(tmp_path):
    mod = api(); root, config, module = fixture_inputs(tmp_path)
    left = saved_run(mod, root, config, module, tmp_path/'left')
    right = saved_run(mod, root, config, module, tmp_path/'right')
    with pytest.raises(ValueError, match='MODE'): mod.compare_runs(left, right, kind='cross-mode')
    copied = tmp_path/'copied'; shutil.copytree(left, copied)
    with pytest.raises(ValueError, match='DIRECTORY|DISTINCT'): mod.compare_runs(left, copied, kind='legacy-binary')
    archive_only = tmp_path/'archive-only'; archive_only.mkdir()
    shutil.copyfile(left/'native-run.jsonl', archive_only/'native-run.jsonl')
    with pytest.raises(ValueError, match='MANIFEST'): mod.compare_runs(left, archive_only, kind='cross-mode')


@pytest.mark.parametrize('mutation', ['config', 'binary', 'source', 'untracked', 'archive',
                                    'durable_modified', 'durable_missing', 'build'])
def test_changed_run_evidence_is_rejected(tmp_path, mutation):
    mod = api(); root, config, module = fixture_inputs(tmp_path)
    left = saved_run(mod, root, config, module, tmp_path/'left')
    right = saved_run(mod, root, config, module, tmp_path/'right', 'disk')
    if mutation == 'config': config.write_text('{}')
    elif mutation == 'binary': module.write_bytes(b'changed native')
    elif mutation == 'source': (root/'src'/'solver.cpp').write_text('changed source')
    elif mutation == 'untracked': (root/'src'/'new.cpp').write_text('untracked effective source')
    elif mutation == 'archive': (right/'native-run.jsonl').write_bytes(b'changed archive')
    elif mutation == 'durable_modified': (right/'iteration-store'/'pass-000000.slqptrace').write_bytes(b'changed record')
    elif mutation == 'durable_missing': (right/'iteration-store'/'pass-000000.slqptrace').unlink()
    else: (module.parent.parent/'CMakeCache.txt').write_text('changed build')
    with pytest.raises(ValueError): mod.compare_runs(left, right, kind='cross-mode')


def test_run_identity_is_automatic_and_source_observation_includes_untracked(tmp_path):
    mod = api(); root, config, module = fixture_inputs(tmp_path)
    a = mod.make_provenance(config, module, tmp_path/'a', 'disk', 1000, source_root=root)
    b = mod.make_provenance(config, module, tmp_path/'b', 'disk', 1000, source_root=root)
    assert a['run_id'] != b['run_id']
    assert str(root/'src'/'solver.cpp') in a['source_observation']['files']
    assert a['source_observation']['dirty_source_diff_sha256']
    assert 'NOT_COMPILATION_ATTESTATION' in a['source_binding']


def test_distinct_valid_manifests_still_require_matching_config_and_binary(tmp_path):
    mod = api(); root, config, module = fixture_inputs(tmp_path)
    left = saved_run(mod, root, config, module, tmp_path/'left')
    other_config = tmp_path/'other.json'; other_config.write_text(json.dumps({'input': json.loads(config.read_text())['input'], 'changed': True}))
    right = saved_run(mod, root, other_config, module, tmp_path/'right', 'disk')
    with pytest.raises(ValueError, match='IDENTITY'): mod.compare_runs(left, right, kind='cross-mode')


def rewrite_json(path, value):
    path.write_text(json.dumps(value, sort_keys=True, indent=2))


def test_summary_archive_digest_must_match_actual_archive(tmp_path):
    mod = api(); root, config, module = fixture_inputs(tmp_path)
    left = saved_run(mod, root, config, module, tmp_path/'left')
    right = saved_run(mod, root, config, module, tmp_path/'right', 'disk')
    summary = json.loads((right/'summary.json').read_text())
    summary['archive']['archive_sha256'] = '0' * 64
    rewrite_json(right/'summary.json', summary)
    rewrite_json(right/'outputs.json', mod.output_hashes(right))
    with pytest.raises(ValueError, match='ARCHIVE'): mod.compare_runs(left, right, kind='cross-mode')


def test_native_writer_evidence_uses_bound_hash_and_footer_without_archive_readback(tmp_path, monkeypatch):
    mod = api(); root, config, module = fixture_inputs(tmp_path)
    run = saved_run(mod, root, config, module, tmp_path/'native', 'disk',
                    storage_verification='COMPLETE_NATIVE_STREAM_FOOTER_DIGEST')
    monkeypatch.setattr(mod.archive, 'read_archive',
                        lambda *_args, **_kwargs: pytest.fail('native writer path must not parse full archive'))

    _provenance, summary = mod.verified_run(run)

    assert summary['storage_verification'] == 'COMPLETE_NATIVE_STREAM_FOOTER_DIGEST'
    assert summary['archive']['records'] == 3


def test_native_writer_evidence_rejects_footer_digest_mismatch(tmp_path, monkeypatch):
    mod = api(); root, config, module = fixture_inputs(tmp_path)
    run = saved_run(mod, root, config, module, tmp_path/'native', 'disk',
                    storage_verification='COMPLETE_NATIVE_STREAM_FOOTER_DIGEST')
    archive_path = run/'native-run.jsonl'
    raw = archive_path.read_bytes()
    lines = raw.splitlines(keepends=True)
    footer = json.loads(lines[-1])
    footer['payload_sha256'] = '0' * 64
    lines[-1] = (json.dumps(footer, sort_keys=True, separators=(',', ':')) + '\n').encode()
    archive_path.write_bytes(b''.join(lines))
    summary = json.loads((run/'summary.json').read_text())
    summary['archive']['archive_sha256'] = mod.sha256(archive_path)
    rewrite_json(run/'summary.json', summary)
    rewrite_json(run/'outputs.json', mod.output_hashes(run))
    monkeypatch.setattr(mod.archive, 'read_archive',
                        lambda *_args, **_kwargs: pytest.fail('native writer path must not parse full archive'))

    with pytest.raises(ValueError, match='ARCHIVE_PAYLOAD_DIGEST_MISMATCH'):
        mod.verified_run(run)


def test_native_writer_evidence_rejects_boolean_trace_count(tmp_path):
    mod = api(); root, config, module = fixture_inputs(tmp_path)
    run = saved_run(mod, root, config, module, tmp_path/'native', 'disk',
                    storage_verification='COMPLETE_NATIVE_STREAM_FOOTER_DIGEST')
    archive_path = run/'native-run.jsonl'
    lines = archive_path.read_bytes().splitlines(keepends=True)
    header = json.loads(lines[0]); header['trace_counts'] = [True]
    lines[0] = (json.dumps(header, sort_keys=True, separators=(',', ':')) + '\n').encode()
    payload = b''.join(lines[:-1])
    footer = json.loads(lines[-1]); footer['payload_sha256'] = __import__('hashlib').sha256(payload).hexdigest()
    lines[-1] = (json.dumps(footer, sort_keys=True, separators=(',', ':')) + '\n').encode()
    archive_path.write_bytes(b''.join(lines))
    summary = json.loads((run/'summary.json').read_text())
    summary['trace_counts'] = [True]
    summary['archive']['payload_sha256'] = footer['payload_sha256']
    summary['archive']['archive_sha256'] = mod.sha256(archive_path)
    rewrite_json(run/'header.json', header)
    rewrite_json(run/'summary.json', summary)
    rewrite_json(run/'outputs.json', mod.output_hashes(run))

    with pytest.raises(ValueError, match='NATIVE_TRACE_COUNT_MISMATCH'):
        mod.verified_run(run)


def test_duplicate_invocation_and_distinct_binary_are_rejected(tmp_path):
    mod = api(); root, config, module = fixture_inputs(tmp_path)
    left = saved_run(mod, root, config, module, tmp_path/'left')
    right = saved_run(mod, root, config, module, tmp_path/'right', 'disk')
    p = json.loads((right/'provenance.json').read_text())
    p['run_id'] = json.loads((left/'provenance.json').read_text())['run_id']
    rewrite_json(right/'provenance.json', p)
    s = json.loads((right/'summary.json').read_text())
    s.update(run_id=p['run_id'], provenance_sha256=mod.sha256(right/'provenance.json'))
    rewrite_json(right/'summary.json', s)
    rewrite_json(right/'outputs.json', mod.output_hashes(right))
    with pytest.raises(ValueError, match='DISTINCT'): mod.compare_runs(left, right, kind='cross-mode')
    other = tmp_path/'other-build'; (other/'python').mkdir(parents=True)
    (other/'CMakeCache.txt').write_text('different build')
    native = other/'python'/'native.so'; native.write_bytes(b'different binary')
    third = saved_run(mod, root, config, native, tmp_path/'third', 'disk')
    with pytest.raises(ValueError, match='BINARY'): mod.compare_runs(left, third, kind='cross-mode')


def test_cli_compares_run_directories_without_changing_them(tmp_path):
    import subprocess
    import sys
    mod = api(); root, config, module = fixture_inputs(tmp_path)
    left = saved_run(mod, root, config, module, tmp_path/'left')
    right = saved_run(mod, root, config, module, tmp_path/'right', 'disk')
    before = (mod.output_hashes(left), mod.output_hashes(right))
    output = tmp_path/'comparison.json'
    completed = subprocess.run([sys.executable, mod.__file__, 'compare', str(left), str(right),
                                '--kind', 'cross-mode', '--output', str(output)], capture_output=True, text=True)
    assert completed.returncode == 0, completed.stderr
    assert json.loads(output.read_text())['comparison_kind'] == 'cross-mode'
    assert before == (mod.output_hashes(left), mod.output_hashes(right))
