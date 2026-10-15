"""Wait for the one active suffix job; never overlap its native audit."""
from pathlib import Path
import json
import os
import sys
import time
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
sequence=HERE/'noodles-repaired-input-sequence-01';sequence.mkdir(exist_ok=False)
prior=HERE/'noodles-repaired-rational-suffix-job-01'
with (sequence/'events.jsonl').open('x') as stream:
    def event(status,**values):
        stream.write(json.dumps(dict(status=status,**values))+'\n');stream.flush();os.fsync(stream.fileno())
    event('WAIT_ACTIVE_RATIONAL_SUFFIX')
    deadline=time.monotonic()+3700
    while not (prior/'complete.json').exists():
        if time.monotonic()>deadline:
            event('SUFFIX_WAIT_LIMIT_NO_AUDIT_STARTED');sys.exit(1)
        time.sleep(5)
    source=prior/'output/domain-081-rational.json'
    if not source.exists() or not (prior/'output/preparation-report.json').exists():
        event('SUFFIX_FAILED_NO_AUDIT_STARTED');sys.exit(1)
    worker=HERE/'audit_repaired_noodles_rational.py'
    exe=HERE/'noodles-rational-authority-diagnosis-job-01/output/inspect-rational-boundary'
    event('REPAIRED_RATIONAL_BOUNDARY_AUDIT_STARTED')
    result=run_source_job(HERE/'noodles-repaired-rational-audit-job-01',
        ['/home/lin/anaconda3/bin/python3','-B',str(worker),'{job}/output'],
        dict(launcher=Path(__file__),worker=worker,input=source,executable=exe,
            prep=prior/'output/preparation-report.json'))
    event('REPAIRED_RATIONAL_AUDIT_TERMINAL',status_code=result.status,returncode=result.returncode)
    event('INPUT_SEQUENCE_TERMINAL_NO_OPTIMIZATION_STARTED')
