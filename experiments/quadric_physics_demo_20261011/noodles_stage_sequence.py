"""Serial input preprocessing after the current Heart production point exits."""
import json
import os
from pathlib import Path
import subprocess
import time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
PYTHON='/home/lin/anaconda3/bin/python3'


def main():
    queue=HERE/'noodles-dense-v5-preprocessing-sequence-01'
    queue.mkdir(exist_ok=False)
    events=(queue/'events.jsonl').open('x')
    def event(status,**data):
        row=dict(status=status,**data)
        events.write(json.dumps(row)+'\n');events.flush();os.fsync(events.fileno())
        print(json.dumps(row),flush=True)
    current=ROOT/'experiments/faster_heart_contact_20261011/rounded-clearance-v3-faster-production-v3/point-job-00'
    event('WAIT_HEART_PRODUCTION_POINT',job=str(current))
    deadline=time.monotonic()+3660
    while not (current/'complete.json').exists():
        if time.monotonic()>deadline:raise RuntimeError('EXISTING_POINT_NOT_EXITED_NO_CONCURRENT_JOB_LAUNCHED')
        time.sleep(2)
    event('HEART_POINT_EXITED',receipt=json.loads((current/'result.json').read_text()))
    previous=None
    for index in range(11):
        job=HERE/('noodles-dense-v5-domain-stage-'+str(index).zfill(2))
        command=[PYTHON,'-B',str(HERE/'noodles_domain_stages.py'),'--input',str(HERE/'noodles-dense-v5'),'--job',str(job)]
        if previous:command+=['--previous',str(previous)]
        event('STAGE_STARTED',index=index,job=str(job))
        with (queue/('stage-'+str(index)+'.log')).open('xb') as stream:
            result=subprocess.run(command,stdout=stream,stderr=subprocess.STDOUT)
        if result.returncode:
            event('STAGE_FAILED',index=index,returncode=result.returncode);return
        previous=job/'output'
        stage=json.loads((previous/'stage.json').read_text())
        event('STAGE_COMPLETED',index=index,excluded_parts=len(stage['excluded_parts']),complete=stage['preprocessing_complete'])
        if stage['preprocessing_complete']:
            event('PREPROCESSING_COMPLETE_GEOMETRY_ADMISSION_PENDING',output=str(previous));break
    events.close()


if __name__=='__main__':main()
