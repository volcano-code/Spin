"""Real transient user services; no mocked cgroup enforcement or paper claims."""
import importlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import pytest


@pytest.fixture
def api():
    try:return importlib.import_module('make_it_stand.verification.source_job')
    except ModuleNotFoundError as exc:
        if exc.name!='make_it_stand.verification.source_job':raise
        class Missing:
            def __getattr__(self,name):pytest.fail('whole-job cgroup envelope is not implemented')
        return Missing()


def run(api,tmp_path,code,**changes):
    source=tmp_path/'input.txt';source.touch(exist_ok=True)
    options=dict(maximum_wall_ns=5_000_000_000,maximum_memory_bytes=256*1024**2,
        maximum_disk_bytes=8*1024**2,minimum_start_free_bytes=0,minimum_free_bytes=0)
    options.update(changes)
    return api.run_source_job(tmp_path/'job',['/usr/bin/python3','-c',code],{'input':source},**options)


def test_kernel_limits_and_whole_workflow_are_recorded(api,tmp_path):
    r=run(api,tmp_path,'import time;print("solve");time.sleep(.1);print("seal");time.sleep(.1);print("done")')
    assert r.status=='PROCESS_EXITED' and r.returncode==0
    assert r.elapsed_ns>=200_000_000
    assert (r.directory/'stdout.log').read_text().splitlines()==['solve','seal','done']
    kernel=json.loads((r.directory/'kernel.json').read_text())
    assert kernel['memory_max']==256*1024**2 and kernel['swap_max']==0
    assert r.unit in kernel['cgroup'] and kernel['boot_id']
    assert json.loads((r.directory/'result.json').read_text())['scientific_validation']=='NOT_GRANTED'


def test_post_solve_sealing_is_inside_same_deadline(api,tmp_path):
    r=run(api,tmp_path,'import time;print("SOLVE_DONE",flush=True);time.sleep(10)',maximum_wall_ns=1_000_000_000)
    assert r.status=='RESOURCE_LIMIT_WALL'
    assert 'SOLVE_DONE' in (r.directory/'stdout.log').read_text()
    assert r.elapsed_ns<4_000_000_000


def test_memory_limit_is_shared_across_forked_children(api,tmp_path):
    r=run(api,tmp_path,'''import os,time
for i in range(2):
    if os.fork()==0:
        allocation=bytearray(192*1024**2)
        time.sleep(10)
        os._exit(0)
time.sleep(10)
''')
    assert r.status=='RESOURCE_LIMIT_MEMORY'
    assert r.service['Result']=='oom-kill'


def test_failed_process_is_not_scientific_success(api,tmp_path):
    r=run(api,tmp_path,'print("failed evidence");raise SystemExit(7)')
    assert r.status=='PROCESS_FAILED' and r.returncode==7
    assert 'failed evidence' in (r.directory/'stdout.log').read_text()


def test_existing_job_cannot_be_restarted_with_new_budget(api,tmp_path):
    r=run(api,tmp_path,'pass');before=(r.directory/'result.json').read_bytes()
    with pytest.raises((FileExistsError,ValueError)):run(api,tmp_path,'pass',maximum_wall_ns=8_000_000_000)
    assert (r.directory/'result.json').read_bytes()==before


def test_disk_limit_stops_job_and_retains_output(api,tmp_path):
    r=run(api,tmp_path,'import sys,time;sys.stdout.write("x"*1000000);sys.stdout.flush();time.sleep(10)',maximum_disk_bytes=200000)
    assert r.status=='RESOURCE_LIMIT_DISK'
    assert (r.directory/'stdout.log').stat().st_size>=1000000


def test_missing_free_space_blocks_payload(api,tmp_path):
    with pytest.raises(ValueError,match='FREE_SPACE'):
        run(api,tmp_path,'print("must not execute")',minimum_start_free_bytes=2**62)
    assert not (tmp_path/'job/stdout.log').exists()


def test_same_group_and_new_session_descendants_are_stopped(api,tmp_path):
    r=run(api,tmp_path,'''import os,time
if os.fork()==0:
    os.setsid()
    print(os.getpid(),flush=True)
    time.sleep(20)
else:time.sleep(.1)
''')
    assert r.status=='PROCESS_EXITED'
    pid=int((r.directory/'stdout.log').read_text())
    until=time.monotonic()+2
    while time.monotonic()<until:
        stat=Path(f'/proc/{pid}/stat')
        if not stat.exists() or stat.read_text().split(') ',1)[1][0]=='Z':break
        time.sleep(.01)
    else:pytest.fail('descendant escaped whole-job cleanup')


@pytest.mark.parametrize('activation_delay',[0,.5])
def test_service_deadline_survives_controller_sigkill(api,tmp_path,activation_delay):
    source=tmp_path/'input.txt';source.touch()
    root=tmp_path/'job'
    program='''from make_it_stand.verification.source_job import run_source_job
import sys,subprocess,time
original=subprocess.run
def delayed(command,*args,**kw):
    if command[0]=='/usr/bin/systemd-run':time.sleep(float(sys.argv[3]))
    return original(command,*args,**kw)
subprocess.run=delayed
run_source_job(sys.argv[1],['/usr/bin/python3','-c','import os,time;print(os.getpid(),flush=True);time.sleep(20)'],
 {'input':sys.argv[2]},maximum_wall_ns=1500000000,maximum_memory_bytes=268435456,
 maximum_disk_bytes=8388608,minimum_start_free_bytes=0,minimum_free_bytes=0)
'''
    child=subprocess.Popen([sys.executable,'-c',program,str(root),str(source),str(activation_delay)])
    unit=None
    try:
        until=time.monotonic()+5
        while time.monotonic()<until:
            assert child.poll() is None
            if (root/'stdout.log').exists() and (root/'stdout.log').stat().st_size:
                unit=json.loads((root/'config.json').read_text())['unit']
                pid=int((root/'stdout.log').read_text());break
            time.sleep(.005)
        else:pytest.fail('payload did not start')
        child.kill();child.wait(timeout=3)
        deadline=json.loads((root/'config.json').read_text())['deadline_boottime_ns']
        until=deadline+250_000_000
        while time.clock_gettime_ns(time.CLOCK_BOOTTIME)<until:
            p=Path(f'/proc/{pid}/stat')
            # The independently enforced deadline can reap the process between
            # stat/open/read. ESRCH/ENOENT here mean the expected exit, not a
            # failed containment test; permission and other I/O errors propagate.
            try: state=p.read_text().split(') ',1)[1][0]
            except (FileNotFoundError,ProcessLookupError):break
            if state=='Z':break
            time.sleep(.02)
        else:pytest.fail('payload survived independent job deadline')
        assert not (root/'complete.json').exists()
        with pytest.raises(FileExistsError):run(api,tmp_path,'pass')
    finally:
        if child.poll() is None:child.kill();child.wait(timeout=3)
        if unit:
            subprocess.run(['/usr/bin/systemctl','--user','stop',unit],capture_output=True,timeout=3)
            subprocess.run(['/usr/bin/systemctl','--user','reset-failed',unit],capture_output=True,timeout=3)


def test_pinned_worker_imports_ignore_payload_pythonpath(api,tmp_path):
    shadow=tmp_path/'shadow';package=shadow/'make_it_stand';package.mkdir(parents=True)
    marker=tmp_path/'shadow-executed'
    (package/'__init__.py').write_text(f'from pathlib import Path;Path({str(marker)!r}).touch();raise RuntimeError("shadow package")')
    r=run(api,tmp_path,'print("correct payload")',environment={'PYTHONPATH':str(shadow)})
    assert r.status=='PROCESS_EXITED' and not marker.exists()


def test_retargeted_executable_alias_cannot_run_unverified_bytes(api,tmp_path,monkeypatch):
    good=tmp_path/'good';bad=tmp_path/'bad';alias=tmp_path/'alias'
    good.write_text('#!/bin/sh\nprintf GOOD');good.chmod(0o700)
    bad.write_text('#!/bin/sh\nprintf UNVERIFIED');bad.chmod(0o700);alias.symlink_to(good)
    publish=api._publish
    def retarget(path,data):
        publish(path,data)
        if path.name=='config.json':alias.unlink();alias.symlink_to(bad)
    monkeypatch.setattr(api,'_publish',retarget)
    with pytest.raises(ValueError,match='IDENTITY'):
        api.run_source_job(tmp_path/'job',[str(alias)],{'input':good},
            maximum_wall_ns=5_000_000_000,maximum_memory_bytes=256*1024**2,
            maximum_disk_bytes=8*1024**2,minimum_start_free_bytes=0,minimum_free_bytes=0)
    assert b'UNVERIFIED' not in (tmp_path/'job/stdout.log').read_bytes()


@pytest.mark.parametrize('key,value',[
    ('maximum_wall_ns',True),('maximum_wall_ns',0),('maximum_memory_bytes',-1),
    ('maximum_disk_bytes',0),('minimum_start_free_bytes',-1)])
def test_invalid_limits_do_not_start_service(api,tmp_path,key,value):
    with pytest.raises(ValueError):run(api,tmp_path,'pass',**{key:value})
    assert not (tmp_path/'job/stdout.log').exists()
