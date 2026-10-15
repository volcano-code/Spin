from pathlib import Path
import json
import sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
demo=ROOT/'experiments/quadric_physics_demo_20261011';case=demo/'noodles-dense-v5'
manifest=case/'manifest.json';data=json.loads(manifest.read_text())
worker=HERE/'diagnose_noodles_rational71_72.py'
native=HERE/'exact-core-material-build-job-01/output/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
files=dict(launcher=Path(__file__),worker=worker,native=native,manifest=manifest,
    prefix=demo/'noodles-dense-v5-core-preprocessing-v1/domain-job-01/output/domain-070.json',
    witness=HERE/'noodles-rational-authority-diagnosis-job-01/output/vertex-link-witness.json',
    helper=HERE/'noodles_vertex_link_witness.py')
for index in (70,71):
    name=data['construction']['fixed_parts'][index];files[name]=case/data['parts'][name]['model_obj']
result=run_source_job(HERE/'noodles-rational71-72-diagnosis-job-01',
    ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
     '/home/lin/anaconda3/bin/python3','-B',str(worker),'{job}/output'],files,
    environment={'PYTHONPATH':str(ROOT/'python')+':'+str(native.parent)})
print(result.status,result.returncode,flush=True)
