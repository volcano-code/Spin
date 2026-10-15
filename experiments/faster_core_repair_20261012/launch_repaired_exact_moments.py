"""New full preregistered point set after the demonstrated exact-policy fix."""
from pathlib import Path
import json
import sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
heart=ROOT/'experiments/faster_heart_contact_20261011'
request=heart/'rounded-clearance-v3-faster-frozen-v2/request.json'
native=HERE/'exact-core-material-build-job-01/output/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
worker=HERE/'probe_repaired_heart_moments.py'
files=dict(launcher=Path(__file__),worker=worker,native=native,request=request,
    states=HERE/'frozen-heart-states.json',integrator=HERE/'independent_tetra_moments.py',
    weights=HERE/'probe_inprocess_derivatives.py',serializer=HERE/'run_small_inprocess_solver.py',
    settings=heart/'core_faster_heart.py')
files.update({f'input_{i}':Path(path) for i,path in enumerate(json.loads(request.read_text())['input_sources'])})
result=run_source_job(HERE/'repaired-heart-exact-production-moments-job-01',
    ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6','/home/lin/anaconda3/bin/python3',
     '-B',str(worker),'--output','{job}/output','--state','heart-first-small-face-rejection',
     '--material-policy','exact'],files,
    environment={'PYTHONPATH':str(ROOT/'python')+':'+str(native.parent)})
print(result.status,result.returncode,flush=True)
