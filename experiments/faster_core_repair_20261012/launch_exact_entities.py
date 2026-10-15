from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
native=HERE/'boundary-proof-reuse-build-job-01/output/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
worker=HERE/'compare_exact_entities.py';inputs=HERE/'same-input-job-03/output'
files=dict(worker=worker,launcher=Path(__file__),native=native,
    reader=ROOT/'experiments/spin_it_reproduction/strict/run_faster_reference_boolean.py')
for state in ('ellipsoid','one_sheet','two_sheet'):
    for name in ('design.obj','source.obj','cherchi-0.obj','cherchi-0.obj.exact.tsv'):
        files[state+'/'+name]=inputs/state/name
    files[state+'/evaluation']=inputs/(state+'-cgal-evaluation.json')
result=run_source_job(HERE/'exact-entities-job-01',
    ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
     '/home/lin/anaconda3/bin/python3','-B',str(worker),'--inputs',str(inputs),'--output','{job}/output'],files,
    environment={'PYTHONPATH':str(ROOT/'python')+':'+str(native.parent)})
print(result.status,result.returncode,flush=True)
