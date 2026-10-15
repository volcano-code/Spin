"""Project-wide admission for existing systemd resource jobs.

The advisory lock is not evidence of worker termination. A durable pending
record and fresh unit/cgroup observations survive loss of the controller.
Legacy launchers must be migrated; they cannot be forced to respect this lock.
"""
from contextlib import contextmanager
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import time
import uuid


_TERMINAL = frozenset(('PROCESS_EXITED', 'PROCESS_FAILED', 'RESOURCE_LIMIT_WALL',
                       'RESOURCE_LIMIT_MEMORY', 'RESOURCE_LIMIT_DISK',
                       'RESOURCE_LIMIT_FREE_SPACE'))


@contextmanager
def project_job_lock(project):
    project = Path(project).resolve(strict=True)
    path = project / '.faster-project-job.lock'
    fd = os.open(path, os.O_CREAT | os.O_RDWR | os.O_CLOEXEC | os.O_NOFOLLOW, 0o600)
    try:
        try:
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError('PROJECT_JOB_BUSY') from exc
        yield
    finally:
        os.close(fd)


def cgroup_populated(group):
    if not group:
        return False
    relative = Path(group)
    if not relative.is_absolute() or '..' in relative.parts:
        raise ValueError('INVALID_CONTROL_GROUP')
    path = Path('/sys/fs/cgroup') / str(relative).lstrip('/')
    try:
        values = dict(line.split() for line in (path/'cgroup.events').read_text().splitlines())
    except FileNotFoundError:
        if path.exists():
            raise RuntimeError('CONTROL_GROUP_OBSERVATION_INCOMPLETE')
        return False
    if values.get('populated') not in ('0', '1'):
        raise RuntimeError('CONTROL_GROUP_OBSERVATION_INCOMPLETE')
    return values['populated'] == '1'


def _boot_id():
    return Path('/proc/sys/kernel/random/boot_id').read_text().strip()


def _reconciliation_path(project, directory, config_bytes):
    key = hashlib.sha256(str(directory).encode() + config_bytes + _boot_id().encode()).hexdigest()
    return Path(project)/'.faster-job-scheduler'/('reconciliation-'+key+'.json')


def _terminal_record(directory, project=None):
    try:
        config_bytes = (directory/'config.json').read_bytes()
        config = json.loads(config_bytes)
        if not (directory/'result.json').exists() and project is not None:
            record = json.loads(_reconciliation_path(project, directory, config_bytes).read_text())
            if (record.get('schema') != 'source-job-termination-observation-v1' or
                    record.get('directory') != str(directory) or
                    record.get('unit') != config.get('unit') or
                    record.get('configuration_sha256') != hashlib.sha256(config_bytes).hexdigest() or
                    record.get('observation_boot_id') != _boot_id() or
                    record.get('status') != 'RECONCILED_STOPPED_WITHOUT_RESULT' or
                    record.get('scientific_validation') != 'NOT_GRANTED'):
                raise ValueError('invalid termination observation')
            return record
        result = json.loads((directory/'result.json').read_text())
    except (OSError, ValueError) as exc:
        raise RuntimeError('PROJECT_JOB_UNFINALIZED: '+str(directory)) from exc
    if (result.get('status') not in _TERMINAL or not config.get('unit') or
            config['unit'] != result.get('unit')):
        raise RuntimeError('PROJECT_JOB_UNFINALIZED: '+str(directory))
    return result


def assert_project_idle(project, units, populated):
    project = Path(project).resolve(strict=True)
    for name, state in units:
        working = state.get('WorkingDirectory', '')
        if not working or not Path(working).is_absolute():
            raise RuntimeError('PROJECT_UNIT_OBSERVATION_INCOMPLETE: '+name)
        directory = Path(working).resolve()
        if not directory.is_relative_to(project):
            continue
        group_live = populated(state.get('ControlGroup', ''))
        if (group_live or state.get('MainPID', '0') != '0' or
                state.get('ActiveState') in ('activating', 'deactivating', 'reloading') or
                state.get('SubState') not in ('dead', 'exited', 'failed')):
            raise RuntimeError('PROJECT_WORKLOAD_ACTIVE: '+name)
        result = _terminal_record(directory, project)
        if result['unit'] != name:
            raise RuntimeError('PROJECT_UNIT_IDENTITY_MISMATCH: '+name)


def assert_pending_finalized(project, pending, populated):
    pending = Path(pending)
    if not pending.exists():
        return
    try:
        record = json.loads(pending.read_text())
        directory = Path(record['directory']).resolve()
    except (OSError, ValueError, KeyError, TypeError) as exc:
        raise RuntimeError('PROJECT_PENDING_RECORD_INVALID') from exc
    if not directory.is_relative_to(Path(project).resolve(strict=True)):
        raise RuntimeError('PROJECT_PENDING_RECORD_INVALID')
    result = _terminal_record(directory, Path(project).resolve(strict=True))
    groups = result.get('observed_empty_control_groups', []) + [result.get('service', {}).get('ControlGroup', '')]
    if any(populated(group) for group in groups if group):
        raise RuntimeError('PROJECT_WORKLOAD_ACTIVE: '+result['unit'])


def register_legacy_source_jobs(project, directories):
    """Persist explicit migration inventory; absent systemd units do not erase it."""
    from .source_job import _publish
    project = Path(project).resolve(strict=True)
    with project_job_lock(project):
        state = project/'.faster-job-scheduler'
        state.mkdir(exist_ok=True, mode=0o700)
        for directory in directories:
            directory = Path(directory).resolve(strict=True)
            if directory == project or not directory.is_relative_to(project):
                raise ValueError('JOB_OUTSIDE_PROJECT')
            data = (directory/'config.json').read_bytes()
            config = json.loads(data)
            if config.get('schema') != 'source-job-v1' or not re.fullmatch(
                    r'spin2014-job-[a-zA-Z0-9_-]+\.service', str(config.get('unit'))):
                raise ValueError('INVALID_LEGACY_SOURCE_JOB')
            key = hashlib.sha256(str(directory).encode()).hexdigest()
            record = dict(directory=str(directory), unit=config['unit'],
                          configuration_sha256=hashlib.sha256(data).hexdigest())
            path = state/('legacy-'+key+'.json')
            if path.exists():
                if json.loads(path.read_text()) != record:
                    raise RuntimeError('LEGACY_CONFIGURATION_CHANGED')
            else:
                _publish(path, record)


def reconcile_terminated_source_job(project, directory):
    """Record fresh termination evidence separately; never invent a job result."""
    from .source_job import _show, _publish
    project = Path(project).resolve(strict=True)
    directory = Path(directory).resolve(strict=True)
    if directory == project or not directory.is_relative_to(project):
        raise ValueError('JOB_OUTSIDE_PROJECT')
    with project_job_lock(project):
        data = (directory/'config.json').read_bytes()
        config = json.loads(data)
        if (config.get('schema') != 'source-job-v1' or not config.get('boot_id') or
                not re.fullmatch(r'spin2014-job-[a-zA-Z0-9_-]+\.service', str(config.get('unit')))):
            raise ValueError('INVALID_LEGACY_SOURCE_JOB')
        if (directory/'result.json').exists():
            raise ValueError('JOB_ALREADY_HAS_RESULT')
        service = _show(config['unit'])
        same_boot = config['boot_id'] == _boot_id()
        groups = {service.get('ControlGroup', '')}
        kernel_path = directory/'kernel.json'
        if same_boot and kernel_path.exists():
            kernel = json.loads(kernel_path.read_text())
            if kernel.get('boot_id') != config['boot_id']:
                raise RuntimeError('JOB_KERNEL_BOOT_MISMATCH')
            groups.add(kernel.get('cgroup', ''))
        if (any(cgroup_populated(group) for group in groups if group) or
                service.get('MainPID') != '0' or
                service.get('ActiveState') not in ('inactive', 'failed') or
                service.get('SubState') not in ('dead', 'exited', 'failed')):
            raise RuntimeError('PROJECT_WORKLOAD_ACTIVE: '+config['unit'])
        if same_boot and not any(groups):
            raise RuntimeError('PROJECT_JOB_RECONCILIATION_UNPROVEN')
        record = dict(schema='source-job-termination-observation-v1', directory=str(directory),
            unit=config['unit'], configuration_sha256=hashlib.sha256(data).hexdigest(),
            observation_boot_id=_boot_id(), observed_boottime_ns=time.clock_gettime_ns(time.CLOCK_BOOTTIME),
            service=dict(service, ControlGroup=next((g for g in groups if g), '')),
            observed_empty_control_groups=sorted(g for g in groups if g),
            status='RECONCILED_STOPPED_WITHOUT_RESULT', scientific_validation='NOT_GRANTED')
        state = project/'.faster-job-scheduler'; state.mkdir(exist_ok=True, mode=0o700)
        path = _reconciliation_path(project, directory, data)
        if not path.exists():
            _publish(path, record)
        return json.loads(path.read_text())


def _assert_legacy_finalized(project, populated):
    for path in (project/'.faster-job-scheduler').glob('legacy-*.json'):
        record = json.loads(path.read_text())
        directory = Path(record['directory']).resolve()
        if directory == project or not directory.is_relative_to(project):
            raise RuntimeError('PROJECT_LEGACY_RECORD_INVALID')
        if hashlib.sha256((directory/'config.json').read_bytes()).hexdigest() != record['configuration_sha256']:
            raise RuntimeError('LEGACY_CONFIGURATION_CHANGED')
        result = _terminal_record(directory, project)
        if result['unit'] != record['unit']:
            raise RuntimeError('PROJECT_UNIT_IDENTITY_MISMATCH')
        groups = result.get('observed_empty_control_groups', []) + [result.get('service', {}).get('ControlGroup', '')]
        if any(populated(group) for group in groups if group):
            raise RuntimeError('PROJECT_WORKLOAD_ACTIVE: '+result['unit'])


def observe_units():
    listing = subprocess.run(['/usr/bin/systemctl', '--user', 'list-units',
        'spin2014-job-*.service', '--all', '--plain', '--no-legend', '--no-pager'],
        capture_output=True, text=True, timeout=5, check=True)
    units = []
    for line in listing.stdout.splitlines():
        name = line.split()[0]
        if not re.fullmatch(r'spin2014-job-[a-zA-Z0-9_-]+\.service', name):
            raise RuntimeError('PROJECT_UNIT_LIST_INVALID')
        shown = subprocess.run(['/usr/bin/systemctl', '--user', 'show', name,
            '--property=WorkingDirectory,ActiveState,SubState,MainPID,ControlGroup'],
            capture_output=True, text=True, timeout=5, check=True)
        state = dict(line.split('=', 1) for line in shown.stdout.splitlines() if '=' in line)
        if not all(key in state for key in ('WorkingDirectory', 'ActiveState', 'SubState', 'MainPID', 'ControlGroup')):
            raise RuntimeError('PROJECT_UNIT_OBSERVATION_INCOMPLETE')
        units.append((name, state))
    return units


def _register(pending, directory):
    # Retain every registration, including a crash before source_job creates its
    # config. Only the pointer is replaced; historical registrations remain.
    registration = pending.parent / ('registration-'+uuid.uuid4().hex+'.json')
    with registration.open('x') as out:
        json.dump(dict(directory=str(directory)), out)
        out.flush()
        os.fsync(out.fileno())
    pointer = pending.parent / ('pending-'+uuid.uuid4().hex)
    os.link(registration, pointer)
    os.replace(pointer, pending)
    fd = os.open(pending.parent, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def run_serial_source_job(project, directory, argv, identity_files, **limits):
    """Fail closed on contention; caller may retry admission, never the job root."""
    from .source_job import run_source_job
    project = Path(project).resolve(strict=True)
    directory = Path(directory).resolve()
    if directory == project or not directory.is_relative_to(project):
        raise ValueError('JOB_OUTSIDE_PROJECT')
    if directory.exists():
        raise FileExistsError(directory)
    with project_job_lock(project):
        state = project / '.faster-job-scheduler'
        state.mkdir(exist_ok=True, mode=0o700)
        pending = state/'pending.json'
        _assert_legacy_finalized(project, cgroup_populated)
        assert_pending_finalized(project, pending, cgroup_populated)
        assert_project_idle(project, observe_units(), cgroup_populated)
        identities = dict(identity_files)
        if 'project_serial_scheduler' in identities:
            raise ValueError('RESERVED_SCHEDULER_IDENTITY')
        identities['project_serial_scheduler'] = Path(__file__).resolve()
        _register(pending, directory)
        result = run_source_job(directory, argv, identities, **limits)
        assert_pending_finalized(project, pending, cgroup_populated)
        assert_project_idle(project, observe_units(), cgroup_populated)
        return result
