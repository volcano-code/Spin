from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
native=HERE/'newton-capture-build-job-01/output/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
bridge=HERE/'cherchi-inprocess-build-job-01/output/_cherchi_core_boolean.cpython-313-x86_64-linux-gnu.so'
worker=HERE/'audit_captured_newton.py'
baseline=HERE/'small-top-inprocess-solver-sequence-01/seed-1-job-01/output'
saved=HERE/'native-newton-matrix-audit-job-02/output'
result=run_source_job(HERE/'native-newton-matrix-recheck-job-01',
    ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6','/home/lin/anaconda3/bin/python3','-B',str(worker),'--output','{job}/output','--saved',str(saved)],
    dict(worker=worker,launcher=Path(__file__),native=native,bridge=bridge,
         adapter=HERE/'cherchi_core_map.py',symbolic=HERE/'independent_eq14.py',
         serialization=HERE/'run_small_inprocess_solver.py',config=baseline/'configuration.json',events=baseline/'events.jsonl',
         saved_events=saved/'events.jsonl',saved_report=saved/'report.json'),
    environment={'PYTHONPATH':str(ROOT/'python')+':'+str(native.parent)+':'+str(bridge.parent)})
print(result.status,result.returncode,flush=True)
