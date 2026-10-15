from pathlib import Path
import json
import sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
case=ROOT/'experiments/quadric_physics_demo_20261011/noodles-dense-v5'
manifest=case/'manifest.json';data=json.loads(manifest.read_text())
if not (HERE/'noodles-rational71-72-diagnosis-job-01/complete.json').exists():
    raise ValueError('SERIAL_DIAGNOSIS_NOT_TERMINAL')
worker=HERE/'continue_noodles_repaired_rational_suffix.py'
native=HERE/'exact-core-material-build-job-01/output/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
files=dict(launcher=Path(__file__),worker=worker,native=native,manifest=manifest,
    start=HERE/'noodles-rational71-72-diagnosis-job-01/output/original-072-rational.json')
files.update({name:case/data['parts'][name]['model_obj'] for name in data['construction']['fixed_parts'][72:]})
result=run_source_job(HERE/'noodles-repaired-rational-suffix-job-01',
    ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
     '/home/lin/anaconda3/bin/python3','-B',str(worker),'{job}/output'],files,
    environment={'PYTHONPATH':str(ROOT/'python')+':'+str(native.parent)})
print(result.status,result.returncode,flush=True)
