from pathlib import Path
import json
import os
import sys
import time
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
prior=HERE/'noodles-all-original-rational-job-02'
sequence=HERE/'noodles-all-original-input-sequence-01';sequence.mkdir(exist_ok=False)
with (sequence/'events.jsonl').open('x') as stream:
    def event(status,**values):
        stream.write(json.dumps(dict(status=status,**values))+'\n');stream.flush();os.fsync(stream.fileno())
    event('WAIT_ACTIVE_ALL_ORIGINAL_RATIONAL_CONSTRUCTION');deadline=time.monotonic()+3700
    while not (prior/'complete.json').exists():
        if time.monotonic()>deadline:event('WAIT_LIMIT_NO_AUDIT_STARTED');sys.exit(1)
        time.sleep(5)
    source=prior/'output/domain-081-rational.json'
    if not source.exists() or not (prior/'output/preparation-report.json').exists():
        event('CONSTRUCTION_FAILED_OR_INCOMPLETE_NO_AUDIT_STARTED');sys.exit(1)
    worker=HERE/'audit_repaired_noodles_rational.py'
    exe=HERE/'noodles-rational-authority-diagnosis-job-01/output/inspect-rational-boundary'
    event('ALL_ORIGINAL_RATIONAL_BOUNDARY_AUDIT_STARTED')
    result=run_source_job(HERE/'noodles-all-original-rational-audit-job-01',
        ['/home/lin/anaconda3/bin/python3','-B',str(worker),'{job}/output',str(source)],
        dict(launcher=Path(__file__),worker=worker,input=source,executable=exe,
            prep=prior/'output/preparation-report.json'))
    event('ALL_ORIGINAL_RATIONAL_AUDIT_TERMINAL',status_code=result.status,returncode=result.returncode)
    event('ALL_ORIGINAL_INPUT_SEQUENCE_TERMINAL_NO_OPTIMIZATION_STARTED')
