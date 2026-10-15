from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
native=HERE/'boundary-proof-reuse-build-job-01/output/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
bridge=HERE/'cherchi-inprocess-build-job-01/output/_cherchi_core_boolean.cpython-313-x86_64-linux-gnu.so'
worker=HERE/'probe_inprocess_derivatives.py'
run=HERE/'small-top-inprocess-solver-sequence-01/seed-1-job-01/output'
result=run_source_job(HERE/'small-top-derivative-probe-job-02',
    ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6','/home/lin/anaconda3/bin/python3','-B',str(worker),
     '--run',str(run),'--output','{job}/output'],
    dict(worker=worker,launcher=Path(__file__),native=native,bridge=bridge,adapter=HERE/'cherchi_core_map.py',
        independent_integrator=HERE/'independent_tetra_moments.py',
        fixture=HERE/'same-input-job-03/output/frozen-inputs.json',
        configuration=run/'configuration.json',result=run/'numerical-result.json'),
    environment={'PYTHONPATH':str(ROOT/'python')+':'+str(native.parent)+':'+str(bridge.parent)})
print(result.status,result.returncode,flush=True)
