from pathlib import Path
import json
import sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
case=ROOT/'experiments/quadric_physics_demo_20261011/noodles-dense-v5'
manifest=case/'manifest.json';data=json.loads(manifest.read_text())
if not (HERE/'noodles-repaired-rational-audit-job-01/complete.json').exists():
    raise ValueError('SERIAL_PREVIOUS_AUDIT_NOT_TERMINAL')
worker=HERE/'construct_noodles_original_rational.py'
native=HERE/'exact-core-material-build-job-01/output/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
files=dict(launcher=Path(__file__),worker=worker,native=native,manifest=manifest,
    encoder=ROOT/'experiments/quadric_physics_demo_20261011/core_boolean_mesh.py')
files.update({name:case/data['parts'][name]['model_obj']
    for name in ['inner-ellipsoid-candidate','outer-ellipsoid']+data['construction']['fixed_parts']})
raise RuntimeError('DO_NOT_RESTART_FULL_SEQUENCE_BEFORE_CONSTRUCTION_PERFORMANCE_REPAIR')
result=run_source_job(HERE/'noodles-all-original-rational-job-03',
    ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
     '/home/lin/anaconda3/bin/python3','-B',str(worker),'{job}/output'],files,
    environment={'PYTHONPATH':str(ROOT/'python')+':'+str(native.parent)})
print(result.status,result.returncode,flush=True)
