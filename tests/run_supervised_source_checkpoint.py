"""Native mid-stage SIGKILL + sealed-boundary recovery; synthetic cases only."""
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import time

from make_it_stand.verification.source_checkpoint_store import StageCheckpointStore
from make_it_stand.verification.source_run_supervisor import StageRunSupervisor


def setup(spec_path,directory):
    spec=json.loads(spec_path.read_text())
    repo=Path(__file__).resolve().parents[1]
    identity={'spec':spec_path,'binary':spec['binary'],'checker':Path(__file__).resolve(),
              'codec':repo/'src/spin2014/SourceCheckpointCodec.cpp',
              'solver':repo/'src/spin2014/SourceAdaptiveSpinSolver.cpp',
              'fixture':repo/'tests/spin_it_2014_checkpoint_codec.cpp',
              'store':repo/'python/make_it_stand/verification/source_checkpoint_store.py'}
    supervisor=StageRunSupervisor(directory,identity,maximum_wall_ns=30_000_000_000,
        maximum_disk_bytes=256*1024**2,maximum_memory_bytes=2*1024**3,
        minimum_free_bytes=20*1024**3,poll_seconds=.005)
    return spec,identity,supervisor


def attempt(spec_path,directory,mode,snapshot=None):
    spec,_,supervisor=setup(spec_path,directory)
    command=[spec['binary'],mode,str(spec['fixture']),'{attempt}/native']
    if snapshot is not None:command.append(str(snapshot))
    result=supervisor.run(command,attempt_wall_ns=10_000_000_000)
    print(json.dumps(dict(index=result.index,status=result.status,returncode=result.returncode,
        charged_wall_ns=result.charged_wall_ns,cumulative_wall_ns=result.cumulative_wall_ns,
        peak_disk_bytes=result.peak_disk_bytes,stdout=str(result.stdout),stderr=str(result.stderr))))


def check(binary,root):
    root.mkdir(exist_ok=False)
    env=dict(os.environ,PYTHONPATH=str(Path(__file__).resolve().parents[1]/'python'))
    outcomes=[]
    for fixture in [0,2]:
        for repeat in range(2):
            work=root/f'fixture-{fixture}-repeat-{repeat}';work.mkdir()
            spec_path=work/'spec.json'
            spec_path.write_text(json.dumps({'binary':str(binary),'fixture':fixture,'scope':'SYNTHETIC_REGRESSION'}))
            def launch(directory,mode,snapshot=None):
                command=[sys.executable,str(Path(__file__).resolve()),'--attempt',str(spec_path),str(directory),mode]
                if snapshot:command.append(str(snapshot))
                return subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE,env=env)
            def finish(child):
                out,err=child.communicate(timeout=35)
                assert child.returncode==0,(out,err)
                return json.loads(out)
            baseline=finish(launch(work/'baseline','baseline'))
            assert baseline['returncode']==0 and baseline['status']=='PROCESS_EXITED',baseline
            expected=Path(baseline['stdout']).read_bytes();assert expected[:1]==b'F'
            directory=work/'resumed'
            initial=finish(launch(directory,'start'))
            _,identity,supervisor=setup(spec_path,directory)
            store=StageCheckpointStore(directory/'checkpoints',identity,
                maximum_wall_ns=30_000_000_000,maximum_trace_bytes=256*1024**2,
                maximum_payload_bytes=64*1024**2)
            boundaries=0
            def seal(result):
                nonlocal boundaries
                data=Path(result['stdout']).read_bytes();assert data[:1]==b'P'
                # Failed partial traces are retained under the supervised root,
                # but never inserted into a scientifically valid checkpoint.
                traces={str(p.relative_to(directory)):p for p in directory.glob('attempt-*/native/*.slqptrace')
                        if p.read_bytes()[-24:-16]==b'SLQPEND1'}
                store.publish(boundaries,payload=data[1:],cumulative_wall_ns=result['cumulative_wall_ns'],
                    cumulative_trace_bytes=sum(p.stat().st_size for p in traces.values()),artifacts=traces)
                loaded=store.load_latest();assert loaded.stage==boundaries
                snapshot=directory/f'snapshot-{boundaries}.bin';snapshot.write_bytes(loaded.payload)
                boundaries+=1
                return snapshot
            snapshot=seal(initial)
            old_checkpoint=store.load_latest().sha256
            child=launch(directory,'resume',snapshot)
            process_record=directory/'attempt-00000001/process.json'
            trace=directory/'attempt-00000001/native/pass-000001.slqptrace'
            until=time.monotonic()+10;pid=None;killed=False
            try:
                while time.monotonic()<until:
                    assert child.poll() is None,'native stage finished before SIGKILL handshake'
                    if process_record.exists() and trace.exists() and trace.stat().st_size>=96:
                        pid=json.loads(process_record.read_text())['data']['pid']
                        os.kill(pid,signal.SIGSTOP)
                        prefix=trace.read_bytes()
                        assert prefix[-24:-16]!=b'SLQPEND1','SIGKILL must occur before stage seal'
                        os.kill(pid,signal.SIGKILL);killed=True;break
                    time.sleep(.0005)
                assert killed,'no active native SLQP trace observed'
                failed=finish(child)
            finally:
                if child.poll() is None:
                    if pid:
                        try:os.killpg(pid,signal.SIGKILL)
                        except ProcessLookupError:pass
                    child.kill();child.wait(timeout=5)
            assert failed['returncode']==-signal.SIGKILL and failed['status']=='PROCESS_FAILED',failed
            assert failed['charged_wall_ns']>0 and trace.read_bytes()==prefix
            assert store.load_latest().sha256==old_checkpoint
            result=finish(launch(directory,'resume',snapshot))
            while Path(result['stdout']).read_bytes()[:1]==b'P':
                snapshot=seal(result)
                result=finish(launch(directory,'resume',snapshot))
            assert result['returncode']==0 and result['status']=='PROCESS_EXITED',result
            actual=Path(result['stdout']).read_bytes()
            assert actual==expected,'recovered scientific output or complete iteration bytes changed'
            state=supervisor.inspect()
            assert state['attempts']==boundaries+2,'failed attempt disappeared or stage counter reset'
            assert result['cumulative_wall_ns']>failed['cumulative_wall_ns']>initial['cumulative_wall_ns']
            outcome=dict(fixture=fixture,repeat=repeat,boundaries=boundaries,
                attempts=state['attempts'],killed_attempt=failed,final_cumulative_wall_ns=state['cumulative_wall_ns'],
                scientific_sha256=hashlib.sha256(actual).hexdigest(),
                partial_trace_sha256=hashlib.sha256(prefix).hexdigest())
            (work/'result.json').write_text(json.dumps(outcome,indent=2)+'\n');outcomes.append(outcome)
    report=dict(scope='SYNTHETIC_NATIVE_MID_STAGE_SIGKILL_RECOVERY',runs=outcomes,
        original_paper_result='NOT_CLAIMED',supervisor_death='FAIL_CLOSED_NOT_AUTO_RECOVERY')
    (root/'result.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))


if __name__=='__main__':
    if sys.argv[1]=='--attempt':
        attempt(Path(sys.argv[2]),Path(sys.argv[3]),sys.argv[4],Path(sys.argv[5]) if len(sys.argv)>5 else None)
    else:
        binary=Path(sys.argv[1]).resolve(strict=True)
        if len(sys.argv)==3:check(binary,Path(sys.argv[2]).absolute())
        else:
            with tempfile.TemporaryDirectory(prefix='spin-supervised-native-') as temporary:
                check(binary,Path(temporary)/'evidence')
