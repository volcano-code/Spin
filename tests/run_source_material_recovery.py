"""Cross-process recovery of the locked reconstructed G6 teapot protocol.

This is an expected binary-validation failure, never a successful paper case.
The test asserts the SAME failed numerical outcome and absence of material.
Successful independently sealed material is tested separately on simple bodies.
"""
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

from make_it_stand.verification.source_checkpoint_store import StageCheckpointStore
from make_it_stand.verification.source_run_supervisor import StageRunSupervisor

ROOT=Path(__file__).resolve().parents[1]


def encode(value):
    return json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()


def worker(spec_path,output,checkpoint=None):
    import trimesh
    import _make_it_stand as native
    from make_it_stand.verification.source_checkpoint_store import _fingerprints
    spec=json.loads(spec_path.read_text())
    assert _fingerprints(spec['files'])==spec['identity'],'input/source/binary changed'
    assert Path(native.__file__).resolve()==Path(spec['files']['native']).resolve()
    cfg=json.loads(Path(spec['files']['config']).read_text())
    inputs=Path(cfg['input'])
    outer=trimesh.load(inputs/'outer.obj',force='mesh',process=False)
    design=trimesh.load(inputs/'design.obj',force='mesh',process=False)
    output.mkdir(exist_ok=False)
    controls={} if checkpoint=='BASELINE' else dict(maximum_new_passes=1,
        resume_checkpoint=None if checkpoint is None else Path(checkpoint).read_bytes())
    run=native.solve_source_spin2014_material(outer.vertices,outer.faces,design.vertices,design.faces,
        source_settings=cfg['geometry'],solver_settings=cfg['solver'],boundary_settings=cfg['boundary'],
        density=1.,minimum_outer_wall=cfg['minimum_wall'],iteration_store_directory=str(output/'trace'),
        iteration_store_max_bytes=256*1024**2,**controls)
    (output/'diagnostics.json').write_bytes(encode(run.diagnostics))
    if run.checkpoint is not None:(output/'checkpoint.bin').write_bytes(run.checkpoint)
    if run.has_geometry:
        raise AssertionError('locked negative regression unexpectedly has geometry; review before admission')
    assert _fingerprints(spec['files'])==spec['identity'],'changed during execution'
    print(json.dumps(dict(status=run.diagnostics['solve']['status'],has_geometry=run.has_geometry)))


def campaign(root,native):
    from make_it_stand.verification.source_checkpoint_store import _fingerprints
    import runpy
    # This is the exact file named in the preceding G6 run's provenance, not
    # the older metric-tangent protocol with a different Laplacian weight.
    cfg_path=ROOT/'build/spin2014-primal-box-20260915/low.json'
    cfg=json.loads(cfg_path.read_text())
    # Reuse manifest admission; no changes to model/pose/design/solver settings.
    adapter=runpy.run_path(str(ROOT/'experiments/teapot_spin2014/supervised.py'))
    admitted=adapter['prepare'](cfg['input'],'inertia',native,'/usr/lib/x86_64-linux-gnu/libstdc++.so.6')
    files={k:str(Path(v).resolve()) for k,v in admitted['identity_files'].items()}
    files.update(config=str(cfg_path),checker=str(Path(__file__).resolve()))
    root.mkdir(exist_ok=False)
    spec=root/'spec.json';spec.write_bytes(encode(dict(files=files,identity=_fingerprints(files))))
    files['spec']=str(spec)
    outcomes=[]
    for repeat in range(2):
        work=root/f'repeat-{repeat}';work.mkdir()
        attempts=work/'attempts'
        supervision=dict(maximum_wall_ns=300_000_000_000,maximum_disk_bytes=1024**3,
            maximum_memory_bytes=24*1024**3,minimum_free_bytes=20*1024**3,poll_seconds=.005)
        supervisor=StageRunSupervisor(attempts,files,**supervision)
        def command(checkpoint):
            args=['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
                sys.executable,str(Path(__file__).resolve()),'worker',str(spec),'{attempt}/payload']
            if checkpoint is not None:args.append(str(checkpoint))
            return args
        def run(checkpoint=None):
            result=supervisor.run(command(checkpoint),attempt_wall_ns=120_000_000_000)
            assert result.status=='PROCESS_EXITED' and result.returncode==0,(result,result.stderr.read_text())
            return result, result.stdout.parent/'payload'
        baseline,baseline_dir=run('BASELINE')
        expected=(baseline_dir/'diagnostics.json').read_bytes()
        assert json.loads(expected)['solve']['status']=='BINARY_VALIDATION_FAILED'
        initial,initial_dir=run()
        assert (initial_dir/'checkpoint.bin').is_file()
        store=StageCheckpointStore(work/'checkpoints',files,maximum_wall_ns=300_000_000_000,
            maximum_trace_bytes=256*1024**2,maximum_payload_bytes=64*1024**2)
        boundaries=0
        def publish(result,directory):
            nonlocal boundaries
            # Baseline traces and incomplete failed attempts do not enter state.
            traces={str(p.relative_to(attempts)):p for p in attempts.glob('attempt-*/payload/trace/*.slqptrace')
                if baseline_dir not in p.parents and p.read_bytes()[-24:-16]==b'SLQPEND1'}
            store.publish(boundaries,payload=(directory/'checkpoint.bin').read_bytes(),
                cumulative_wall_ns=result.cumulative_wall_ns,
                cumulative_trace_bytes=sum(p.stat().st_size for p in traces.values()),artifacts=traces)
            loaded=store.load_latest();assert loaded.stage==boundaries
            path=work/f'resume-{boundaries}.bin';path.write_bytes(loaded.payload);boundaries+=1
            return path
        checkpoint=publish(initial,initial_dir);old=store.load_latest().sha256
        # Separate supervisor process remains alive to account for killed worker.
        child=subprocess.Popen([sys.executable,str(Path(__file__).resolve()),'attempt',str(spec),str(attempts),str(checkpoint)],
            stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        process_record=attempts/'attempt-00000002/process.json'
        trace=attempts/'attempt-00000002/payload/trace/pass-000001.slqptrace'
        until=time.monotonic()+100;pid=None;partial=None
        try:
            while time.monotonic()<until:
                assert child.poll() is None,'stage exited before kill handshake'
                if process_record.exists() and trace.exists() and trace.stat().st_size>=96:
                    pid=json.loads(process_record.read_text())['data']['pid']
                    os.kill(pid,signal.SIGSTOP);partial=trace.read_bytes()
                    assert partial[-24:-16]!=b'SLQPEND1','stage already sealed'
                    os.kill(pid,signal.SIGKILL);break
                time.sleep(.0005)
            assert partial is not None,'no in-progress native trace'
            stdout,stderr=child.communicate(timeout=15)
            assert child.returncode==0,(stdout,stderr)
            failed=json.loads(stdout)
            assert failed['returncode']==-9 and failed['status']=='PROCESS_FAILED'
        finally:
            if child.poll() is None:
                if pid:
                    try:os.killpg(pid,signal.SIGKILL)
                    except ProcessLookupError:pass
                child.kill();child.wait(timeout=5)
        assert trace.read_bytes()==partial and store.load_latest().sha256==old
        for _ in range(cfg['solver']['maximum_adaptive_passes']):
            result,directory=run(checkpoint)
            if not (directory/'checkpoint.bin').exists():break
            checkpoint=publish(result,directory)
        actual=(directory/'diagnostics.json').read_bytes()
        assert actual==expected,'full resumed numerical/geometry diagnostics differ'
        assert supervisor.inspect()['attempts']==boundaries+3
        outcomes.append(dict(repeat=repeat,status=json.loads(actual)['solve']['status'],
            has_geometry=False,sealed_artifact=False,boundaries=boundaries,
            attempts=supervisor.inspect()['attempts'],killed_attempt=failed,
            scientific_sha256=hashlib.sha256(actual).hexdigest(),partial_sha256=hashlib.sha256(partial).hexdigest()))
    assert outcomes[0]['scientific_sha256']==outcomes[1]['scientific_sha256']
    report=dict(scope='RECONSTRUCTED_G6_CHECKPOINT_RECOVERY_NEGATIVE_REGRESSION',runs=outcomes,
        original_paper_result='NOT_CLAIMED',optimization_success=False)
    (root/'result.json').write_bytes(encode(report));print(json.dumps(report,indent=2))


if __name__=='__main__':
    if sys.argv[1]=='worker':worker(Path(sys.argv[2]),Path(sys.argv[3]),sys.argv[4] if len(sys.argv)>4 else None)
    elif sys.argv[1]=='attempt':
        spec=Path(sys.argv[2]);data=json.loads(spec.read_text());files=dict(data['files'],spec=str(spec))
        sup=StageRunSupervisor(Path(sys.argv[3]),files,maximum_wall_ns=300_000_000_000,
            maximum_disk_bytes=1024**3,maximum_memory_bytes=24*1024**3,minimum_free_bytes=20*1024**3,poll_seconds=.005)
        r=sup.run(['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',sys.executable,
            str(Path(__file__).resolve()),'worker',str(spec),'{attempt}/payload',sys.argv[4]],attempt_wall_ns=120_000_000_000)
        print(json.dumps(dict(returncode=r.returncode,status=r.status,charged_wall_ns=r.charged_wall_ns)))
    else:campaign(Path(sys.argv[1]).absolute(),Path(sys.argv[2]).resolve(strict=True))
