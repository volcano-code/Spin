"""Wait for the active core Heart sequence, then construct the Noodles domain."""
import argparse
import json
import os
from pathlib import Path
import sys
import time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def main(output=None,previous=None):
    sys.path.insert(0,str(ROOT/'python'))
    from make_it_stand.verification.source_job import run_source_job
    output=Path(output) if output is not None else HERE/'noodles-dense-v5-core-preprocessing-v1'
    output.mkdir(exist_ok=False)
    heart=ROOT/'experiments/faster_heart_contact_20261011/rounded-clearance-v3-core-faster-v1'
    with (output/'events.jsonl').open('x') as stream:
        def event(status,**values):
            row=dict(status=status,**values);stream.write(json.dumps(row)+'\n');stream.flush();os.fsync(stream.fileno())
            print(json.dumps(row),flush=True)
        event('WAIT_HEART_CORE_SEQUENCE')
        deadline=time.monotonic()+7300
        while True:
            rows=[json.loads(row) for row in (heart/'events.jsonl').read_text().splitlines()]
            completed=[row for row in rows if row['status']=='CORE_SOLVER_JOB_TERMINAL']
            if any(row['seed']==1 or row['returncode']!=0 or row['status_code']!='PROCESS_EXITED' for row in completed):break
            if time.monotonic()>deadline:raise RuntimeError('HEART_CORE_SEQUENCE_NOT_TERMINAL')
            time.sleep(2)
        case=HERE/'noodles-dense-v5'
        manifest=json.loads((case/'manifest.json').read_text())
        native=ROOT/'build/heart-fixed-facet-20261011-v1/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
        files={'launcher':Path(__file__),'worker':HERE/'noodles_core_preprocess.py','manifest':case/'manifest.json','native':native}
        files['canonicalizer']=HERE/'core_boolean_mesh.py'
        files.update({name:case/row['model_obj'] for name,row in manifest['parts'].items()})
        extra=[]
        if previous is not None:
            previous=Path(previous).resolve();files['core_predecessor']=previous
            extra=['--previous-core',str(previous)]
        event('CORE_DOMAIN_STARTED')
        result=run_source_job(output/'domain-job-01',
            ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
             '/home/lin/anaconda3/bin/python3','-B',str(HERE/'noodles_core_preprocess.py'),
             '--case',str(case),'--output','{job}/output',*extra],files,
             environment={'PYTHONPATH':str(ROOT/'python')+':'+str(native.parent)})
        event('CORE_DOMAIN_TERMINAL',status_code=result.status,returncode=result.returncode)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--directory',type=Path)
    parser.add_argument('--previous-core',type=Path);args=parser.parse_args()
    main(args.directory,args.previous_core)
