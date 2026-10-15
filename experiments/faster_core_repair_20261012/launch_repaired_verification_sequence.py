"""Serialize new production-point checks and the Noodles represented-boundary check."""
import json
import os
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job

sequence=HERE/'repaired-production-verification-sequence-01';sequence.mkdir(exist_ok=False)
native=HERE/'projective-belt-build-job-04/output/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
request=ROOT/'experiments/faster_heart_contact_20261011/rounded-clearance-v3-faster-frozen-v2/request.json'
with (sequence/'events.jsonl').open('x') as stream:
    def event(status,**values):
        stream.write(json.dumps(dict(status=status,**values))+'\n');stream.flush();os.fsync(stream.fileno())
    for ordinal,name in enumerate(('heart-first-cap-failure','heart-first-small-face-rejection'),1):
        worker=HERE/'probe_repaired_heart_moments.py'
        files=dict(launcher=Path(__file__),worker=worker,native=native,request=request,
            states=HERE/'frozen-heart-states.json',integrator=HERE/'independent_tetra_moments.py',
            weights=HERE/'probe_inprocess_derivatives.py',serializer=HERE/'run_small_inprocess_solver.py',
            core_settings=ROOT/'experiments/faster_heart_contact_20261011/core_faster_heart.py')
        files.update({f'input_{i}':Path(path) for i,path in enumerate(json.loads(request.read_text())['input_sources'])})
        job=HERE/f'repaired-heart-production-moments-job-{ordinal:02d}'
        event('PRODUCTION_DERIVATIVE_JOB_STARTED',state=name,job=str(job))
        result=run_source_job(job,
            ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6','/home/lin/anaconda3/bin/python3',
             '-B',str(worker),'--output','{job}/output','--state',name],files,
            environment={'PYTHONPATH':str(ROOT/'python')+':'+str(native.parent)})
        event('PRODUCTION_DERIVATIVE_JOB_TERMINAL',state=name,status_code=result.status,returncode=result.returncode)
        report=job/'output/report.json'
        if not report.exists() or not json.loads(report.read_text())['numerical_checks_passed']:
            event('PRODUCTION_DERIVATIVE_GATE_NOT_PASSED',state=name)
            break
    demo=ROOT/'experiments/quadric_physics_demo_20261011'
    worker=demo/'audit_noodles_core_design.py'
    source=demo/'noodles-dense-v5-core-rational-suffix-v1/domain-job-01/output'
    sequence_job=demo/'noodles-dense-v5-core-design-audit-v1';sequence_job.mkdir(exist_ok=False)
    event('NOODLES_DESIGN_AUDIT_STARTED')
    result=run_source_job(sequence_job/'audit-job-01',
        ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6','/home/lin/anaconda3/bin/python3',
         '-B',str(worker),'--input',str(source),'--output','{job}/output'],
        dict(launcher=Path(__file__),worker=worker,native=native,domain=source/'domain-081.json',
             manifest=demo/'noodles-dense-v5/manifest.json',integrator=HERE/'independent_tetra_moments.py'),
        environment={'PYTHONPATH':str(ROOT/'python')+':'+str(native.parent)})
    event('NOODLES_DESIGN_AUDIT_TERMINAL',status_code=result.status,returncode=result.returncode)
    event('VERIFICATION_SEQUENCE_TERMINAL_NO_OPTIMIZATION_STARTED')
