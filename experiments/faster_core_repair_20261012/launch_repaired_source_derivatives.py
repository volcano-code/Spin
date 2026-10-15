from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
worker=HERE/'probe_repaired_heart_source_derivatives.py'
native=HERE/'projective-belt-build-job-04/output/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
request=ROOT/'experiments/faster_heart_contact_20261011/rounded-clearance-v3-faster-frozen-v2/request.json'
result=run_source_job(HERE/'repaired-heart-source-derivative-job-01',
    ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6','/home/lin/anaconda3/bin/python3',
     '-B',str(worker),'--output','{job}/output'],
    dict(launcher=Path(__file__),worker=worker,native=native,request=request,
         states=HERE/'frozen-heart-states.json',weights=HERE/'probe_inprocess_derivatives.py',
         serializer=HERE/'run_small_inprocess_solver.py'),
    environment={'PYTHONPATH':str(ROOT/'python')+':'+str(native.parent)})
print(result.status,result.returncode,flush=True)
