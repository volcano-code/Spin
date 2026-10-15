"""Wait for genuine moment/capacity completion, then one bounded input export."""
from pathlib import Path
import json
import os
import sys
import time
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
sequence=HERE/'noodles-balanced-shared-export-sequence-01';sequence.mkdir(exist_ok=False)
prior=HERE/'noodles-balanced-input-moments-job-01'
with (sequence/'events.jsonl').open('x') as stream:
    def event(status,**values):
        stream.write(json.dumps(dict(status=status,**values))+'\n');stream.flush();os.fsync(stream.fileno())
    event('WAIT_INDEPENDENT_INPUT_MOMENTS_CAPACITY')
    deadline=time.monotonic()+3700
    while not (prior/'complete.json').exists():
        if time.monotonic()>deadline:event('WAIT_LIMIT_NO_EXPORT_STARTED');sys.exit(1)
        time.sleep(5)
    report=prior/'output/report.json'
    if json.loads((prior/'result.json').read_text())['returncode']!=0 or not report.exists():
        event('MOMENTS_FAILED_NO_EXPORT_STARTED');sys.exit(1)
    data=json.loads(report.read_text())
    if data['capacity_status']=='EXCLUDED_BY_FIRST_MOMENT_BOUND':
        event('CAPACITY_EXCLUDED_NO_EXPORT_OR_OPTIMIZATION_STARTED');sys.exit(1)
    demo=ROOT/'experiments/quadric_physics_demo_20261011'
    heart=ROOT/'experiments/faster_heart_contact_20261011'
    files=dict(launcher=Path(__file__),worker=HERE/'export_balanced_shared_input.py',
        assembly=HERE/'shared_input_boundary.py',representation=demo/'core_boolean_mesh.py',
        integration=heart/'rounded_input.py',report=report,capacity=prior/'output/capacity.json',
        manifest=demo/'noodles-dense-v5/manifest.json',
        native=HERE/'exact-core-material-build-job-01/output/_make_it_stand.cpython-313-x86_64-linux-gnu.so')
    for name in ('design','outer','fixed'):
        files[name]=HERE/f'noodles-balanced-full-input-job-01/output/{name}-rational.json'
        files[name+'-moments']=prior/f'output/{name}-moments-exact.json'
    event('SHARED_INPUT_EXPORT_STARTED')
    result=run_source_job(HERE/'noodles-balanced-shared-export-job-01',
        ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
         '/home/lin/anaconda3/bin/python3','-B',str(files['worker']),'{job}/output'],files,
        environment={'PYTHONPATH':str(ROOT/'python')+':'+str(files['native'].parent)})
    event('SHARED_INPUT_EXPORT_TERMINAL',status_code=result.status,returncode=result.returncode,
        optimization_started=False)
