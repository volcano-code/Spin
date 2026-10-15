from pathlib import Path
import json
import sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
native=HERE/'boundary-proof-reuse-build-job-01/output/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
worker=HERE/'replay_saved_source_profiles.py'
request=ROOT/'experiments/faster_heart_contact_20261011/rounded-clearance-v3-faster-frozen-v2/request.json'
files=dict(worker=worker,launcher=Path(__file__),native=native,request=request,
           states=HERE/'frozen-heart-states.json',serialization=HERE/'run_small_inprocess_solver.py')
files.update({f'input_{i}':Path(path) for i,path in enumerate(json.loads(request.read_text())['input_sources'])})
result=run_source_job(HERE/'saved-source-profile-comparison-job-01',
    ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6','/home/lin/anaconda3/bin/python3','-B',str(worker),'--output','{job}/output'],files,
    environment={'PYTHONPATH':str(ROOT/'python')+':'+str(native.parent)})
print(result.status,result.returncode,flush=True)
