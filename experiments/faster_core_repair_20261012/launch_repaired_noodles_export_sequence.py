from pathlib import Path
import json
import os
import sys
import time
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
prior=HERE/'noodles-repaired-rational-audit-job-01'
sequence=HERE/'noodles-repaired-export-sequence-01';sequence.mkdir(exist_ok=False)
with (sequence/'events.jsonl').open('x') as stream:
    def event(status,**values):
        stream.write(json.dumps(dict(status=status,**values))+'\n');stream.flush();os.fsync(stream.fileno())
    event('WAIT_ACTIVE_RATIONAL_BOUNDARY_AUDIT');deadline=time.monotonic()+3700
    while not (prior/'complete.json').exists():
        if time.monotonic()>deadline:event('WAIT_LIMIT');sys.exit(1)
        time.sleep(5)
    audit=prior/'output/report.json'
    if not audit.exists() or not json.loads(audit.read_text())['boundary_valid']:
        event('RATIONAL_BOUNDARY_FAILED_NO_EXPORT_STARTED');sys.exit(1)
    demo=ROOT/'experiments/quadric_physics_demo_20261011';case=demo/'noodles-dense-v5'
    manifest=case/'manifest.json';data=json.loads(manifest.read_text())
    native=HERE/'exact-core-material-build-job-01/output/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
    worker=HERE/'export_repaired_noodles_design.py'
    files=dict(launcher=Path(__file__),worker=worker,native=native,audit=audit,manifest=manifest,
        authority=HERE/'noodles-repaired-rational-suffix-job-01/output/domain-081-rational.json',
        helper=demo/'core_boolean_mesh.py',integrator=HERE/'independent_tetra_moments.py')
    files.update({name:case/data['parts'][name]['model_obj']
        for name in ['outer-ellipsoid']+data['construction']['fixed_parts']})
    event('DERIVED_DESIGN_EXPORT_AUDIT_STARTED')
    result=run_source_job(HERE/'noodles-repaired-derived-design-job-01',
        ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
         '/home/lin/anaconda3/bin/python3','-B',str(worker),'{job}/output'],files,
        environment={'PYTHONPATH':str(ROOT/'python')+':'+str(native.parent)})
    event('DERIVED_DESIGN_JOB_TERMINAL',status_code=result.status,returncode=result.returncode)
    event('EXPORT_SEQUENCE_TERMINAL_NO_OPTIMIZATION_STARTED')
