import json
import os
from pathlib import Path
import sys
import argparse
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
native=HERE/'boundary-proof-reuse-build-job-01/output/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
bridge=HERE/'cherchi-inprocess-build-job-01/output/_cherchi_core_boolean.cpython-313-x86_64-linux-gnu.so'
p=argparse.ArgumentParser();p.add_argument('--mode',choices=('yoyo','top'),default='yoyo');p.add_argument('--directory',required=True)
args=p.parse_args()
worker=HERE/('run_small_top_inprocess_solver.py' if args.mode=='top' else 'run_small_inprocess_solver.py')
output=HERE/args.directory;output.mkdir(exist_ok=False)
with (output/'events.jsonl').open('x') as stream:
    for seed in (0,1):
        stream.write(json.dumps(dict(status='STARTED',seed=seed))+'\n');stream.flush();os.fsync(stream.fileno())
        result=run_source_job(output/f'seed-{seed}-job-01',
            ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6','/home/lin/anaconda3/bin/python3','-B',
             str(worker),'--output','{job}/output','--seed',str(seed)],
            dict(worker=worker,launcher=Path(__file__),native=native,bridge=bridge,adapter=HERE/'cherchi_core_map.py',
                configuration=HERE/'same-input-job-03/output/frozen-inputs.json',
                common_worker=HERE/'run_small_inprocess_solver.py',
                differential=ROOT/'python/make_it_stand/verification/cherchi_motion_differential.py',
                spin=ROOT/'python/make_it_stand/spin.py'),
            environment={'PYTHONPATH':str(ROOT/'python')+':'+str(native.parent)+':'+str(bridge.parent)})
        stream.write(json.dumps(dict(status=result.status,returncode=result.returncode,seed=seed))+'\n');stream.flush();os.fsync(stream.fileno())
        if result.status!='PROCESS_EXITED' or result.returncode:break
