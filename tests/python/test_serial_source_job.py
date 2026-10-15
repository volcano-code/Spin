"""The lock must guard the workload, not merely its launcher's lifetime."""
import json
import subprocess
import sys

import pytest

from make_it_stand.verification import serial_source_job as serial


def test_project_lock_excludes_another_process(tmp_path):
    with serial.project_job_lock(tmp_path):
        code = "from pathlib import Path; from make_it_stand.verification.serial_source_job import project_job_lock;\nwith project_job_lock(Path(__import__('sys').argv[1])): pass"
        other = subprocess.run([sys.executable, '-c', code, str(tmp_path)], capture_output=True, text=True)
        assert other.returncode != 0
        assert 'PROJECT_JOB_BUSY' in other.stderr
    with serial.project_job_lock(tmp_path):
        pass


def test_controller_exit_does_not_make_live_worker_safe(tmp_path):
    unit = dict(WorkingDirectory=str(tmp_path/'job'), ActiveState='active',
                SubState='running', MainPID='0', ControlGroup='/test/job')
    with pytest.raises(RuntimeError, match='PROJECT_WORKLOAD_ACTIVE'):
        serial.assert_project_idle(tmp_path, [('spin2014-job-test.service', unit)], lambda _: True)


def test_unrelated_project_is_not_blocked(tmp_path):
    unit = dict(WorkingDirectory=str(tmp_path.parent/'other'), ActiveState='active',
                SubState='running', MainPID='123', ControlGroup='/test/job')
    serial.assert_project_idle(tmp_path, [('spin2014-job-test.service', unit)], lambda _: True)


def test_missing_result_blocks_even_after_worker_exits(tmp_path):
    job = tmp_path/'job'; job.mkdir()
    (job/'config.json').write_text(json.dumps(dict(unit='spin2014-job-test.service')))
    unit = dict(WorkingDirectory=str(job), ActiveState='active', SubState='exited',
                MainPID='0', ControlGroup='/test/job')
    with pytest.raises(RuntimeError, match='PROJECT_JOB_UNFINALIZED'):
        serial.assert_project_idle(tmp_path, [('spin2014-job-test.service', unit)], lambda _: False)


@pytest.mark.parametrize('status', ['PROCESS_EXITED', 'RESOURCE_LIMIT_WALL', 'PROCESS_FAILED'])
def test_terminal_record_and_empty_group_allow_next_job(tmp_path, status):
    job = tmp_path/'job'; job.mkdir()
    for name in ('config', 'result'):
        (job/f'{name}.json').write_text(json.dumps(dict(unit='spin2014-job-test.service', status=status)))
    unit = dict(WorkingDirectory=str(job), ActiveState='inactive', SubState='dead',
                MainPID='0', ControlGroup='/test/job')
    serial.assert_project_idle(tmp_path, [('spin2014-job-test.service', unit)], lambda _: False)


def test_terminal_record_never_overrides_populated_group(tmp_path):
    job = tmp_path/'job'; job.mkdir()
    (job/'result.json').write_text(json.dumps(dict(unit='spin2014-job-test.service', status='PROCESS_EXITED')))
    unit = dict(WorkingDirectory=str(job), ActiveState='inactive', SubState='dead',
                MainPID='0', ControlGroup='/test/job')
    with pytest.raises(RuntimeError, match='PROJECT_WORKLOAD_ACTIVE'):
        serial.assert_project_idle(tmp_path, [('spin2014-job-test.service', unit)], lambda _: True)


def test_pending_controller_record_blocks_before_systemd_start(tmp_path):
    pending = tmp_path/'pending.json'
    pending.write_text(json.dumps(dict(directory=str(tmp_path/'not-started'))))
    with pytest.raises(RuntimeError, match='PROJECT_JOB_UNFINALIZED'):
        serial.assert_pending_finalized(tmp_path, pending, lambda _: False)


def test_unknown_group_observation_fails_closed(tmp_path):
    unit = dict(WorkingDirectory=str(tmp_path/'job'), ActiveState='inactive',
                SubState='dead', MainPID='0', ControlGroup='/test/job')
    def unavailable(_):
        raise OSError('observation unavailable')
    with pytest.raises(OSError):
        serial.assert_project_idle(tmp_path, [('spin2014-job-test.service', unit)], unavailable)


def test_wrapper_keeps_lock_until_resource_result_is_published(tmp_path, monkeypatch):
    from make_it_stand.verification import source_job
    monkeypatch.setattr(serial, 'observe_units', lambda: [])
    def resource_job(directory, argv, identities, **limits):
        with pytest.raises(RuntimeError, match='PROJECT_JOB_BUSY'):
            with serial.project_job_lock(tmp_path):
                pass
        directory.mkdir()
        for name in ('config', 'result'):
            (directory/f'{name}.json').write_text(json.dumps(dict(
                unit='spin2014-job-test.service', status='PROCESS_EXITED')))
        return 'published'
    monkeypatch.setattr(source_job, 'run_source_job', resource_job)
    assert serial.run_serial_source_job(tmp_path, tmp_path/'job', ['/bin/true'], {'asset': __file__}) == 'published'
    assert (tmp_path/'.faster-job-scheduler/pending.json').exists()


def test_abandoned_admission_requires_reconciliation_not_automatic_retry(tmp_path, monkeypatch):
    from make_it_stand.verification import source_job
    monkeypatch.setattr(serial, 'observe_units', lambda: [])
    def failure(*args, **kwargs):
        raise RuntimeError('simulated controller loss')
    monkeypatch.setattr(source_job, 'run_source_job', failure)
    with pytest.raises(RuntimeError, match='simulated controller loss'):
        serial.run_serial_source_job(tmp_path, tmp_path/'job', ['/bin/true'], {'asset': __file__})
    with pytest.raises(RuntimeError, match='PROJECT_JOB_UNFINALIZED'):
        serial.run_serial_source_job(tmp_path, tmp_path/'other-job', ['/bin/true'], {'asset': __file__})


def test_result_unit_mismatch_is_not_terminal_proof(tmp_path):
    job = tmp_path/'job'; job.mkdir()
    (job/'config.json').write_text(json.dumps(dict(unit='spin2014-job-real.service')))
    (job/'result.json').write_text(json.dumps(dict(unit='spin2014-job-other.service', status='PROCESS_EXITED')))
    pending = tmp_path/'pending.json'
    pending.write_text(json.dumps(dict(directory=str(job))))
    with pytest.raises(RuntimeError, match='PROJECT_JOB_UNFINALIZED'):
        serial.assert_pending_finalized(tmp_path, pending, lambda _: False)


def test_registered_legacy_missing_result_blocks_when_unit_has_disappeared(tmp_path, monkeypatch):
    job = tmp_path/'old-job'; job.mkdir()
    (job/'config.json').write_text(json.dumps(dict(
        schema='source-job-v1', unit='spin2014-job-old.service', boot_id='prior-boot')))
    register = getattr(serial, 'register_legacy_source_jobs', None)
    assert callable(register), 'legacy work must be registered before migrating launchers'
    register(tmp_path, [job])
    monkeypatch.setattr(serial, 'observe_units', lambda: [])
    with pytest.raises(RuntimeError, match='PROJECT_JOB_UNFINALIZED'):
        serial.run_serial_source_job(tmp_path, tmp_path/'new-job', ['/bin/true'], {'asset': __file__})
    assert not (tmp_path/'new-job').exists()


def test_reconciliation_preserves_missing_result_and_allows_only_terminated_work(tmp_path, monkeypatch):
    from make_it_stand.verification import source_job
    job = tmp_path/'old-job'; job.mkdir()
    boot = serial.Path('/proc/sys/kernel/random/boot_id').read_text().strip()
    (job/'config.json').write_text(json.dumps(dict(
        schema='source-job-v1', unit='spin2014-job-old.service', boot_id=boot)))
    reconcile = getattr(serial, 'reconcile_terminated_source_job', None)
    assert callable(reconcile), 'stopped legacy work needs a separate durable observation'
    monkeypatch.setattr(source_job, '_show', lambda _: dict(
        LoadState='loaded', ActiveState='failed', SubState='failed', MainPID='0',
        ControlGroup='/test/old', Result='signal', ExecMainStatus='9'))
    monkeypatch.setattr(serial, 'cgroup_populated', lambda _: False)
    record = reconcile(tmp_path, job)
    assert record['scientific_validation'] == 'NOT_GRANTED'
    assert not (job/'result.json').exists()
    pending = tmp_path/'pending.json'
    pending.write_text(json.dumps(dict(directory=str(job))))
    serial.assert_pending_finalized(tmp_path, pending, lambda _: False)
    with pytest.raises(RuntimeError, match='PROJECT_WORKLOAD_ACTIVE'):
        serial.assert_pending_finalized(tmp_path, pending, lambda _: True)


def test_reconciliation_refuses_live_group(tmp_path, monkeypatch):
    from make_it_stand.verification import source_job
    job = tmp_path/'job'; job.mkdir()
    (job/'config.json').write_text(json.dumps(dict(
        schema='source-job-v1', unit='spin2014-job-old.service',
        boot_id=serial.Path('/proc/sys/kernel/random/boot_id').read_text().strip())))
    reconcile = getattr(serial, 'reconcile_terminated_source_job', None)
    assert callable(reconcile), 'reconciliation must observe actual workload termination'
    monkeypatch.setattr(source_job, '_show', lambda _: dict(
        LoadState='loaded', ActiveState='failed', SubState='failed', MainPID='0',
        ControlGroup='/test/old', Result='signal'))
    monkeypatch.setattr(serial, 'cgroup_populated', lambda _: True)
    with pytest.raises(RuntimeError, match='PROJECT_WORKLOAD_ACTIVE'):
        reconcile(tmp_path, job)
