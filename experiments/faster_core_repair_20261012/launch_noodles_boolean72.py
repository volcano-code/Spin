from pathlib import Path
import json
import sys
HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
heart=HERE/'repaired-heart-exact-production-moments-job-01'
if not (heart/'complete.json').exists():raise RuntimeError('SERIAL_HEART_JOB_NOT_TERMINAL')
demo=ROOT/'experiments/quadric_physics_demo_20261011'
manifest=demo/'noodles-dense-v5/manifest.json'
data=json.loads(manifest.read_text());name=data['construction']['fixed_parts'][71]
part=demo/'noodles-dense-v5'/data['parts'][name]['model_obj']
native=HERE/'exact-core-material-build-job-01/output/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
worker=HERE/'capture_noodles_boolean72.py'
result=run_source_job(HERE/'noodles-original-boolean72-capture-job-01',
    ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
     '/home/lin/anaconda3/bin/python3','-B',str(worker),'{job}/output'],
    dict(launcher=Path(__file__),worker=worker,native=native,manifest=manifest,part=part,
        predecessor=demo/'noodles-dense-v5-core-preprocessing-v1/domain-job-01/output/domain-071.json',
        witness=HERE/'noodles-rational-authority-diagnosis-job-01/output/vertex-link-witness.json',
        link_helper=HERE/'noodles_vertex_link_witness.py'),
    environment={'PYTHONPATH':str(ROOT/'python')+':'+str(native.parent)})
print(result.status,result.returncode,flush=True)
