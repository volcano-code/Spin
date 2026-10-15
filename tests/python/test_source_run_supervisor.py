"""Real subprocess budgets: synthetic fixtures, no paper success claims."""
import importlib
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

import pytest


@pytest.fixture
def api():
    try:
        return importlib.import_module('make_it_stand.verification.source_run_supervisor')
    except ModuleNotFoundError as exc:
        if exc.name=='make_it_stand.verification.source_run_supervisor':
            class Missing:
                def __getattr__(self,name):
                    pytest.fail('cumulative subprocess supervisor is not implemented')
            return Missing()
        raise


def supervisor(api,tmp_path,**kw):
    source=tmp_path/'input';source.touch(exist_ok=True)
    options=dict(maximum_wall_ns=3_000_000_000,maximum_disk_bytes=8_000_000,
                 maximum_memory_bytes=512*1024**2,minimum_free_bytes=0,poll_seconds=.005)
    options.update(kw)
    return api.StageRunSupervisor(tmp_path/'run',{'input':source},**options)


def command(code):
    return [sys.executable,'-c',code]


def test_disk_scan_tolerates_file_removed_after_enumeration(api,tmp_path,monkeypatch):
    temporary=tmp_path/'CMakeCache.txt.tmp'
    temporary.write_bytes(b'temporary')
    stable=tmp_path/'stable'
    stable.write_bytes(b'x'*8192)
    original=Path.lstat
    def unlink_before_stat(path,*args,**kwargs):
        if path==temporary and temporary.exists(): temporary.unlink()
        return original(path,*args,**kwargs)
    monkeypatch.setattr(Path,'lstat',unlink_before_stat)
    assert api._disk_bytes(tmp_path)>=stable.stat().st_size


def test_disk_scan_counts_internal_library_links_once(api,tmp_path):
    library=tmp_path/'lib.so.1.2';library.write_bytes(b'x'*8192)
    (tmp_path/'lib.so.1').symlink_to(library.name)
    (tmp_path/'lib.so').symlink_to('lib.so.1')
    paths=[tmp_path,*tmp_path.iterdir()]
    expected=sum(max(p.lstat().st_size,p.lstat().st_blocks*512) for p in paths)
    assert api._disk_bytes(tmp_path)==expected


@pytest.mark.parametrize('kind',['external','dangling','directory','loop'])
def test_disk_scan_rejects_unsafe_links(api,tmp_path,kind):
    root=tmp_path/'run';root.mkdir()
    if kind=='external':
        target=tmp_path/'outside';target.write_bytes(b'x')
    elif kind=='directory':
        target=root/'directory';target.mkdir()
    elif kind=='loop':target=root/'link'
    else:target=root/'missing'
    (root/'link').symlink_to(target)
    with pytest.raises(ValueError,match='UNSUPPORTED_RUN_ARTIFACT'):
        api._disk_bytes(root)


def test_disk_scan_does_not_hide_permission_error(api,tmp_path,monkeypatch):
    denied=tmp_path/'denied';denied.touch()
    original=Path.lstat
    def denied_stat(path,*args,**kwargs):
        if path==denied: raise PermissionError('denied')
        return original(path,*args,**kwargs)
    monkeypatch.setattr(Path,'lstat',denied_stat)
    with pytest.raises(PermissionError):api._disk_bytes(tmp_path)


def test_failed_attempt_time_and_outputs_survive_new_supervisor(api,tmp_path):
    s=supervisor(api,tmp_path)
    a=s.run(command('import time;time.sleep(.06);print("failed evidence");raise SystemExit(7)'),attempt_wall_ns=1_000_000_000)
    assert a.returncode==7 and a.status=='PROCESS_FAILED'
    assert a.charged_wall_ns>=60_000_000
    assert a.stdout.read_text().strip()=='failed evidence'
    fresh=supervisor(api,tmp_path)
    b=fresh.run(command('print("success")'),attempt_wall_ns=1_000_000_000)
    assert b.index==1 and b.status=='PROCESS_EXITED'
    assert b.cumulative_wall_ns>=a.charged_wall_ns+b.charged_wall_ns
    assert fresh.inspect()['attempts']==2
    assert a.stdout.read_text().strip()=='failed evidence'


def test_timeout_counts_and_restart_cannot_reset_global_budget(api,tmp_path):
    s=supervisor(api,tmp_path,maximum_wall_ns=200_000_000)
    a=s.run(command('import time;time.sleep(9)'),attempt_wall_ns=10_000_000_000)
    assert a.status=='RESOURCE_LIMIT_WALL' and a.returncode<0
    assert a.cumulative_wall_ns>=200_000_000
    with pytest.raises(ValueError,match='RESOURCE_LIMIT_WALL'):
        supervisor(api,tmp_path,maximum_wall_ns=200_000_000).run(command('print("must not run")'))
    assert not (s.directory/'attempt-00000001').exists()


def test_external_sigkill_is_charged_and_retry_uses_new_directory(api,tmp_path):
    s=supervisor(api,tmp_path)
    a=s.run(command('import os,signal,time;time.sleep(.04);os.kill(os.getpid(),signal.SIGKILL)'))
    assert a.returncode==-signal.SIGKILL and a.charged_wall_ns>=40_000_000
    b=supervisor(api,tmp_path).run(command('print("retry")'))
    assert b.index==1 and b.cumulative_wall_ns>a.cumulative_wall_ns


def test_disk_limit_includes_failed_output_and_blocks_restart(api,tmp_path):
    s=supervisor(api,tmp_path,maximum_disk_bytes=100_000)
    a=s.run(command('import sys,time;sys.stdout.write("x"*300000);sys.stdout.flush();time.sleep(5)'))
    assert a.status=='RESOURCE_LIMIT_DISK' and a.peak_disk_bytes>=300000
    assert a.stdout.stat().st_size>=300000
    with pytest.raises(ValueError,match='RESOURCE_LIMIT_DISK'):
        supervisor(api,tmp_path,maximum_disk_bytes=100_000).run(command('pass'))


def test_final_result_publication_is_included_in_disk_accounting(api,tmp_path):
    s=supervisor(api,tmp_path);a=s.run(command('pass'))
    # Independent filesystem observation, not the supervisor's accounting helper.
    paths=[s.directory,*s.directory.rglob('*')];inodes={}
    for path in paths:
        info=path.stat();inodes[(info.st_dev,info.st_ino)]=info.st_blocks*512
    actual=sum(inodes.values())
    assert a.peak_disk_bytes>=actual
    assert s.inspect()['peak_disk_bytes']>=actual


def test_post_result_check_failure_cannot_be_resumed_as_complete(api,tmp_path,monkeypatch):
    s=supervisor(api,tmp_path);original=api._disk_bytes
    def fail_after_result(root):
        if (root/'attempt-00000000/result.json').exists():
            raise OSError('injected post-result filesystem failure')
        return original(root)
    with monkeypatch.context() as patch:
        patch.setattr(api,'_disk_bytes',fail_after_result)
        with pytest.raises(OSError,match='post-result'):s.run(command('pass'))
    with pytest.raises(ValueError,match='UNFINISHED_ATTEMPT'):
        supervisor(api,tmp_path).run(command('pass'))
    assert not (s.directory/'attempt-00000001').exists()


def test_memory_rlimit_reaches_actual_child(api,tmp_path):
    a=supervisor(api,tmp_path).run(command('import resource;print(resource.getrlimit(resource.RLIMIT_AS)[0])'))
    assert a.returncode==0 and int(a.stdout.read_text())==512*1024**2


def test_allocation_over_cap_fails_without_losing_attempt(api,tmp_path):
    s=supervisor(api,tmp_path)
    a=s.run(command('x=bytearray(1024**3)'))
    assert a.returncode!=0 and b'MemoryError' in a.stderr.read_bytes()
    assert s.inspect()['attempts']==1


def test_identity_changed_during_execution_is_not_returned_as_success(api,tmp_path):
    s=supervisor(api,tmp_path)
    with pytest.raises(ValueError,match='STALE_SUPERVISOR_IDENTITY'):
        s.run(command(f'from pathlib import Path;Path({str(tmp_path/"input")!r}).write_text("changed during run")'))
    record=json.loads((s.directory/'attempt-00000000/result.json').read_text())['data']
    assert record['charged_wall_ns']>0 and record['status']=='SUPERVISOR_ERROR'


@pytest.mark.parametrize('change',['limits','identity','journal'])
def test_changed_identity_limits_or_corrupt_journal_block_launch(api,tmp_path,change):
    s=supervisor(api,tmp_path);s.run(command('pass'))
    if change=='identity':(tmp_path/'input').write_text('changed')
    if change=='journal':(s.directory/'attempt-00000000/result.json').write_text('{}')
    options={'maximum_wall_ns':4_000_000_000} if change=='limits' else {}
    with pytest.raises(ValueError):supervisor(api,tmp_path,**options).run(command('pass'))
    assert not (s.directory/'attempt-00000001').exists()


def test_minimum_free_space_prevents_any_launch(api,tmp_path):
    s=supervisor(api,tmp_path,minimum_free_bytes=2**62)
    with pytest.raises(ValueError,match='RESOURCE_LIMIT_FREE_SPACE'):s.run(command('pass'))
    assert not (s.directory/'attempt-00000000').exists()


def test_launch_error_still_produces_charged_failed_attempt(api,tmp_path):
    s=supervisor(api,tmp_path)
    a=s.run(['/does/not/exist'])
    assert a.returncode!=0 and a.charged_wall_ns>0
    assert s.inspect()['attempts']==1


@pytest.mark.parametrize('field,value',[
    ('maximum_wall_ns',True),('maximum_disk_bytes',0),('maximum_memory_bytes',-1),
    ('poll_seconds',float('nan')),('poll_seconds',0),('minimum_free_bytes',-1)])
def test_invalid_budgets_rejected(api,tmp_path,field,value):
    with pytest.raises(ValueError):supervisor(api,tmp_path,**{field:value})


def test_child_cannot_retain_inherited_supervisor_lock(api,tmp_path):
    a=supervisor(api,tmp_path).run(command('import os;print(len(os.listdir("/proc/self/fd")))'))
    assert a.returncode==0 and int(a.stdout.read_text())<=4


@pytest.mark.parametrize('change',['delete','modify','add','receipt-delete','receipt-modify'])
def test_completed_attempt_evidence_cannot_disappear_before_retry(api,tmp_path,change):
    s=supervisor(api,tmp_path);a=s.run(command('print("preserve me")'))
    if change=='delete':a.stdout.unlink()
    elif change=='modify':a.stdout.write_text('rewritten')
    elif change=='add':(a.stdout.parent/'unexpected').write_text('added after completion')
    elif change=='receipt-delete':(a.stdout.parent/'complete.json').unlink()
    else:(a.stdout.parent/'complete.json').write_text('{}')
    with pytest.raises((ValueError,OSError)):supervisor(api,tmp_path).run(command('pass'))
    assert not (s.directory/'attempt-00000001').exists()


def test_same_group_descendant_does_not_survive_leader_completion(api,tmp_path):
    a=supervisor(api,tmp_path).run(command('''import os,time
pid=os.fork()
if pid==0:time.sleep(20)
else:print(pid,flush=True)
'''))
    assert a.returncode==0
    pid=int(a.stdout.read_text())
    until=time.monotonic()+2
    while time.monotonic()<until:
        info=Path(f'/proc/{pid}/stat')
        if not info.exists() or info.read_text().split(') ',1)[1][0]=='Z':break
        time.sleep(.01)
    else:pytest.fail('same-group descendant still executing after supervisor returned')


def wait_file(path,process):
    until=time.monotonic()+5
    while time.monotonic()<until:
        if path.exists():return
        assert process.poll() is None,'supervisor exited before handshake'
        time.sleep(.005)
    pytest.fail('handshake timeout')


def test_supervisor_death_reserves_budget_and_fails_closed(api,tmp_path):
    # Losing the measurement process must NOT silently refund a reservation.
    source=tmp_path/'input';source.touch()
    program='''from make_it_stand.verification.source_run_supervisor import StageRunSupervisor
import sys
s=StageRunSupervisor(sys.argv[1],{'input':sys.argv[2]},maximum_wall_ns=3000000000,maximum_disk_bytes=8000000,maximum_memory_bytes=536870912,minimum_free_bytes=0,poll_seconds=.005)
s.run([sys.executable,'-c','import time;time.sleep(30)'],attempt_wall_ns=1000000000)
'''
    process=subprocess.Popen([sys.executable,'-c',program,str(tmp_path/'run'),str(source)])
    child=None
    try:
        path=tmp_path/'run/attempt-00000000/process.json';wait_file(path,process)
        child=json.loads(path.read_text())['data']['pid']
        with pytest.raises(ValueError,match='SUPERVISOR_BUSY'):supervisor(api,tmp_path).run(command('pass'))
        process.kill();process.wait(timeout=5)
        with pytest.raises(ValueError,match='UNFINISHED_ATTEMPT'):supervisor(api,tmp_path).run(command('pass'))
        assert not (tmp_path/'run/attempt-00000001').exists()
    finally:
        if process.poll() is None:process.kill();process.wait(timeout=5)
        if child:
            try:os.killpg(child,signal.SIGKILL)
            except ProcessLookupError:pass
