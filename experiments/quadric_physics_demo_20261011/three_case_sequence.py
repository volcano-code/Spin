"""Serial approved scientific jobs, with durable queue state and case gates.

Every scientific worker has its own existing source_job resource envelope.
The controller waits for the already-running Heart job; it never launches a
second large job beside it. Failed scientific gates remain failed.
"""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
PYTHON='/home/lin/anaconda3/bin/python3'
NATIVE=ROOT/'build/spin2014-objectives-quadric-20261007/python/_make_it_stand.cpython-313-x86_64-linux-gnu.so'


def main():
    queue=HERE/'three-case-sequence-20261012-01'
    queue.mkdir(exist_ok=False)
    stream=(queue/'events.jsonl').open('x')
    def event(stage,**values):
        row=dict(stage=stage,**values)
        stream.write(json.dumps(row,allow_nan=False)+'\n')
        stream.flush();os.fsync(stream.fileno())
        print(json.dumps(row),flush=True)
    heart=ROOT/'experiments/faster_heart_contact_20261011/rounded-clearance-v3-bounded-admission-job-03'
    deadline=time.monotonic()+3600
    event('WAIT_EXISTING_HEART_JOB',job=str(heart))
    while not (heart/'complete.json').exists():
        if time.monotonic()>deadline:
            raise RuntimeError('EXISTING_JOB_STILL_RUNNING_NO_CONCURRENT_LARGE_JOB_STARTED')
        time.sleep(2)
    event('HEART_JOB_EXITED',receipt=json.loads((heart/'result.json').read_text()))
    flower=HERE/'paper-flower-method-v2'
    def run(argv,label):
        event('JOB_START',operation=label,argv=argv)
        with (queue/(label+'.log')).open('xb') as log:
            result=subprocess.run(argv,stdout=log,stderr=subprocess.STDOUT)
        event('JOB_EXIT',operation=label,returncode=result.returncode)
        return result.returncode==0
    def flower_job(operation,job):
        argv=[PYTHON,'-B',str(HERE/'dual_track_jobs.py'),operation,
              '--input',str(flower),'--job',str(flower/job)]
        if operation=='saved-check':argv+=['--application',str(flower/'solve-job-01/output/application.json')]
        return run(argv,'flower-'+operation)
    gate_pass=False
    if flower_job('saved-check','saved-independent-job-01'):
        gate=json.loads((flower/'saved-independent-job-01/output/saved-independent-check.json').read_text())
        gate_pass=gate.get('simulation_admission') is True and gate.get('full_material_gate')=='PASSED'
    event('FLOWER_MATERIAL_GATE',passed=gate_pass)
    if gate_pass:
        flower_job('deliver','delivery-job-01')
        if flower_job('physics','physics-job-01'):
            if flower_job('render','render-job-01'):
                flower_job('display','paper-display-job-01')
    run([PYTHON,'-B',str(HERE/'noodles_dense_v5_launch_admission.py'),
         '--input',str(HERE/'noodles-dense-v5'),'--job',str(HERE/'noodles-dense-v5-admission-job-01'),
         '--expected-native-sha256',hashlib.sha256(NATIVE.read_bytes()).hexdigest()], 'noodles-exact-admission')
    event('SEQUENCE_FINISHED',optimization_of_heart_or_noodles_granted=False)
    stream.close()


if __name__=='__main__':main()
