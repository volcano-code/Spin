"""Wait for the active rational audit, then build and measure the minimal material fix."""
import json
import os
from pathlib import Path
import sys
import time

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job

sequence=HERE/'core-material-followup-sequence-01';sequence.mkdir(exist_ok=False)
with (sequence/'events.jsonl').open('x') as events:
    def event(status,**values):
        events.write(json.dumps(dict(status=status,**values))+'\n');events.flush();os.fsync(events.fileno())
    previous=HERE/'noodles-rational-authority-diagnosis-job-01'
    deadline=time.monotonic()+3700
    event('WAIT_ACTIVE_RATIONAL_AUTHORITY_AUDIT')
    while not (previous/'complete.json').exists():
        if time.monotonic()>deadline:
            event('QUEUE_WAIT_LIMIT_NO_BUILD_STARTED');raise SystemExit(0)
        time.sleep(5)
    event('MINIMAL_EXACT_MATERIAL_BUILD_STARTED')
    launcher=HERE/'launch_exact_core_material_build.py'
    # This child uses the same capped source-job supervisor; no unbounded build.
    import subprocess
    subprocess.run(['/home/lin/anaconda3/bin/python3','-B',str(launcher)],check=True)
    native=HERE/'exact-core-material-build-job-01/output/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
    if not native.is_file():
        event('BUILD_FAILED_NO_MEASUREMENT_STARTED');raise SystemExit(0)
    heart=ROOT/'experiments/faster_heart_contact_20261011'
    request=heart/'rounded-clearance-v3-faster-frozen-v2/request.json'
    worker=HERE/'diagnose_heart_exact_material_gate.py'
    files=dict(launcher=Path(__file__),worker=worker,native=native,request=request,
        source_manifest=HERE/'repaired-heart-production-moments-job-02/output/manifest.json',
        settings=heart/'core_faster_heart.py',timer=HERE/'phase_timing.py',serializer=HERE/'run_small_inprocess_solver.py')
    files.update({f'input_{i}':Path(path) for i,path in enumerate(json.loads(request.read_text())['input_sources'])})
    event('EXACT_MATERIAL_FIXED_POINT_MEASUREMENT_STARTED')
    result=run_source_job(HERE/'heart-exact-material-gate-diagnosis-job-02',
        ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6','/home/lin/anaconda3/bin/python3',
         '-B',str(worker),'--output','{job}/output'],files,
        environment={'PYTHONPATH':str(ROOT/'python')+':'+str(native.parent)})
    event('EXACT_MATERIAL_FIXED_POINT_MEASUREMENT_TERMINAL',status_code=result.status,returncode=result.returncode)
    event('FOLLOWUP_SEQUENCE_TERMINAL_NO_OPTIMIZATION_STARTED')
