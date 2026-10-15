"""Wait for the existing tape job, then one small closure job; no case autorun."""
from pathlib import Path
import hashlib
import json
import os
import subprocess
import sys
import time
HERE=Path(__file__).resolve().parent
prior=HERE/'noodles-native-tape-certification-job-01'
out=HERE/'short-closeout-sequence-01';out.mkdir(exist_ok=False)
bound=[HERE/'launch_small_proof_closure.py',HERE/'run_small_proof_closure.py']
digests={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in bound}
with (out/'events.jsonl').open('x') as stream:
    def event(status,**values):
        stream.write(json.dumps(dict(status=status,**values))+'\n');stream.flush();os.fsync(stream.fileno())
    event('WAIT_EXISTING_TAPE_CAPTURE',job=str(prior),next_scripts=digests)
    deadline=time.monotonic()+3700
    while not (prior/'complete.json').exists():
        if time.monotonic()>deadline:
            event('WAIT_LIMIT_NO_JOB_STARTED');sys.exit(1)
        time.sleep(5)
    result=json.loads((prior/'result.json').read_text())
    report=prior/'output/report.json'
    if result.get('returncode')!=0 or not report.exists() or not json.loads(report.read_text()).get('complete_frozen_sample_certification'):
        event('STOP_TAPE_CERTIFICATION_NOT_PASSED');sys.exit(1)
    if any(hashlib.sha256(Path(p).read_bytes()).hexdigest()!=expected for p,expected in digests.items()):
        event('STOP_NEXT_WORKER_CHANGED');sys.exit(1)
    event('SMALL_CLOSURE_STARTED')
    completed=subprocess.run([sys.executable,'-B',str(bound[0])],check=False)
    event('SMALL_CLOSURE_TERMINAL',returncode=completed.returncode,
          next_action='Inspect actual closure report before any full case optimization')
