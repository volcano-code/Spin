from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
native=HERE/'projective-interval-build-job-02/output/_make_it_stand.cpython-313-x86_64-linux-gnu.so';saved=HERE/'noodles-interval-source-derivative-job-01/output'
files=dict(launcher=Path(__file__),worker=HERE/'probe_interval_algebraic_source_derivative.py',native=native,weights=HERE/'probe_inprocess_derivatives.py',serializer=HERE/'run_small_inprocess_solver.py',config=HERE/'noodles-bound-seed0-production-probe-job-01/output/run/mapping-config.json',manifest=saved/'manifest.json')
for i in range(81):files[f'point{i}']=saved/f'point-{i:03d}.json'
result=run_source_job(HERE/'noodles-interval-algebraic-source-derivative-job-01',
 ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6','/home/lin/anaconda3/bin/python3','-B',str(files['worker']),'{job}/output'],files,
 environment={'PYTHONPATH':str(ROOT/'python')+':'+str(native.parent)+':'+str(HERE)})
print(result.status,result.returncode,flush=True)
