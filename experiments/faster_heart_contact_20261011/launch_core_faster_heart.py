"""Serial seeded Faster core runs, with the existing bounded job supervisor."""
import argparse
import json
import os
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def main(output=None):
    sys.path.insert(0,str(ROOT/'python'))
    from make_it_stand.verification.source_job import run_source_job
    output=Path(output) if output is not None else HERE/'rounded-clearance-v3-core-faster-v1'
    output.mkdir(exist_ok=False)
    request=HERE/'rounded-clearance-v3-faster-frozen-v2/request.json'
    native=ROOT/'build/heart-fixed-facet-20261011-v1/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
    files={'launcher':Path(__file__),'worker':HERE/'core_faster_heart.py',
        'journal':HERE/'rounded_solver_journal.py','request':request,'native':native,
        'spin':ROOT/'python/make_it_stand/spin.py'}
    packet=json.loads(request.read_text())
    files.update({'input-'+str(index):Path(path) for index,path in enumerate(packet['input_sources'])})
    with (output/'events.jsonl').open('x') as stream:
        def event(status,**values):
            row=dict(status=status,**values);stream.write(json.dumps(row)+'\n');stream.flush();os.fsync(stream.fileno())
            print(json.dumps(row),flush=True)
        for seed in (0,1):
            event('CORE_SOLVER_STARTED',seed=seed)
            result=run_source_job(output/f'seed-{seed}-job-01',
                ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
                 '/home/lin/anaconda3/bin/python3','-B',str(HERE/'core_faster_heart.py'),
                 '--request',str(request),'--seed',str(seed),'--output','{job}/output'],files,
                environment={'PYTHONPATH':str(ROOT/'python')+':'+str(native.parent)})
            event('CORE_SOLVER_JOB_TERMINAL',seed=seed,status_code=result.status,returncode=result.returncode)
            if result.status!='PROCESS_EXITED' or result.returncode!=0:break


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--directory',type=Path)
    main(parser.parse_args().directory)
