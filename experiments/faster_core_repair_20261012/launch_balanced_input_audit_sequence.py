"""Serial boundary checks after a complete bounded construction; no optimizer."""
from pathlib import Path
import json
import os
import sys
import time
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
sequence=HERE/'noodles-balanced-input-verification-sequence-01';sequence.mkdir(exist_ok=False)
prior=HERE/'noodles-balanced-full-input-job-01'
with (sequence/'events.jsonl').open('x') as stream:
    def event(status,**values):
        stream.write(json.dumps(dict(status=status,**values))+'\n');stream.flush();os.fsync(stream.fileno())
    event('WAIT_FULL_BALANCED_CONSTRUCTION')
    deadline=time.monotonic()+3700
    while not (prior/'complete.json').exists():
        if time.monotonic()>deadline:
            event('WAIT_LIMIT_NO_AUDIT_STARTED');sys.exit(1)
        time.sleep(5)
    result=json.loads((prior/'result.json').read_text())
    prep=prior/'output/report.json'
    if result['returncode']!=0 or not prep.exists() or not json.loads(prep.read_text())['complete_81_part_construction']:
        event('CONSTRUCTION_FAILED_NO_AUDIT_STARTED');sys.exit(1)
    worker=HERE/'audit_repaired_noodles_rational.py'
    exe=HERE/'noodles-rational-authority-diagnosis-job-01/output/inspect-rational-boundary'
    for name in ['design','outer','fixed']:
        source=prior/'output'/f'{name}-rational.json'
        job=HERE/f'noodles-balanced-{name}-audit-job-01'
        event('EXACT_BOUNDARY_AUDIT_STARTED',material=name)
        result=run_source_job(job,
            ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
             '/home/lin/anaconda3/bin/python3','-B',str(worker),'{job}/output',str(source)],
            dict(launcher=Path(__file__),worker=worker,input=source,executable=exe,prep=prep))
        report=job/'output/report.json'
        if result.returncode!=0 or not report.exists() or not json.loads(report.read_text())['boundary_valid']:
            event('BOUNDARY_FAILED_NO_OPTIMIZATION_STARTED',material=name,result=result.status);sys.exit(1)
        event('BOUNDARY_PASSED',material=name)
    event('THREE_BOUNDARIES_PASSED_PARTITION_MOMENTS_CAPACITY_AND_OPTIMIZATION_STILL_PENDING')
