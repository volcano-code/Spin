"""Serial registered production points after current Noodles preprocessing."""
import json
import os
from pathlib import Path
import subprocess
import time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
PYTHON='/home/lin/anaconda3/bin/python3'


def main():
    queue=HERE/'rounded-clearance-v3-production-sequence-01'
    queue.mkdir(exist_ok=False)
    events=(queue/'events.jsonl').open('x')
    def event(status,**data):
        row=dict(status=status,**data)
        events.write(json.dumps(row)+'\n');events.flush();os.fsync(events.fileno())
        print(json.dumps(row),flush=True)
    production=HERE/'rounded-clearance-v3-faster-production-v3'
    base=json.loads((production/'point-job-00/output/point.json').read_text())
    if not base['derivative_defined'] or not base['material_facts']['valid'] or max(base['moment_errors'])>1e-6:
        event('BASE_POINT_BLOCKED');return
    stages=ROOT/'experiments/quadric_physics_demo_20261011/noodles-dense-v5-preprocessing-sequence-01/events.jsonl'
    event('WAIT_NOODLES_INPUT_PREPROCESSING')
    deadline=time.monotonic()+12*3600
    while True:
        lines=stages.read_text().splitlines()
        try:last=json.loads(lines[-1]) if lines else {}
        except json.JSONDecodeError:last={}
        if last.get('status') in ('STAGE_FAILED','PREPROCESSING_COMPLETE_GEOMETRY_ADMISSION_PENDING'):break
        if time.monotonic()>deadline:raise RuntimeError('PREPROCESSING_NOT_TERMINAL_NO_CONCURRENT_POINT_LAUNCHED')
        time.sleep(3)
    event('INPUT_PREPROCESSING_TERMINAL',record=last)
    for index in range(1,81):
        job=production/('point-job-'+str(index).zfill(2))
        command=[PYTHON,'-B',str(HERE/'rounded_faster_production_point.py'),'launch',
                 '--registration',str(production/'registration.json'),'--index',str(index),'--job',str(job)]
        event('POINT_STARTED',index=index,job=str(job))
        with (queue/('point-'+str(index)+'.log')).open('xb') as stream:
            result=subprocess.run(command,stdout=stream,stderr=subprocess.STDOUT)
        if result.returncode:
            event('POINT_FAILED',index=index,returncode=result.returncode);return
        record=json.loads((job/'output/point.json').read_text())
        local=record['derivative_defined'] and record['material_facts']['valid'] and max(record['moment_errors'])<=1e-6
        event('POINT_COMPLETED',index=index,local_mapping_checks=local)
        if not local:
            event('LOCAL_MAPPING_BLOCKER',index=index);return
    command=[PYTHON,'-B',str(HERE/'aggregate_rounded_production.py'),'--directory',str(production),'--output',str(production/'aggregate-01')]
    with (queue/'aggregate.log').open('xb') as stream:
        result=subprocess.run(command,stdout=stream,stderr=subprocess.STDOUT)
    event('AGGREGATE_EXITED',returncode=result.returncode)
    if result.returncode==0:
        report=json.loads((production/'aggregate-01/report.json').read_text())
        event('DERIVATIVE_BATCH_FINISHED',passed=report['derivative_admission'],optimizer_invocations=0)
    events.close()


if __name__=='__main__':main()
