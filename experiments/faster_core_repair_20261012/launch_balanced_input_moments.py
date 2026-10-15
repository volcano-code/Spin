from pathlib import Path
import json
import sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
DEMO=ROOT/'experiments/quadric_physics_demo_20261011'
HEART=ROOT/'experiments/faster_heart_contact_20261011'
case=DEMO/'noodles-dense-v5';manifest=case/'manifest.json';data=json.loads(manifest.read_text())
files=dict(launcher=Path(__file__),worker=HERE/'inspect_balanced_input_moments.py',
    capacity=HERE/'saved_input_capacity.py',integration=HEART/'rounded_input.py',
    physical=HEART/'contact_matrix.py',
    common=DEMO/'core_boolean_mesh.py',connectivity=DEMO/'noodles_admission.py',manifest=manifest,
    executable=HERE/'exact-input-union-build-job-01/output/exact-input-boolean',
    construction_report=HERE/'noodles-balanced-full-input-job-01/output/report.json')
for name in ['inner-ellipsoid-candidate','outer-ellipsoid']:
    files[name]=case/data['parts'][name]['model_obj']
for name in ['design','outer','fixed']:
    files[name]=HERE/f'noodles-balanced-full-input-job-01/output/{name}-rational.json'
    files[name+'-audit']=HERE/f'noodles-balanced-{name}-audit-job-01/output/report.json'
    files[name+'-receipt']=HERE/f'noodles-balanced-{name}-audit-job-01/output/receipt.json'
result=run_source_job(HERE/'noodles-balanced-input-moments-job-01',
    ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
     '/home/lin/anaconda3/bin/python3','-B',str(files['worker']),'{job}/output'],files)
print(result.status,result.returncode,flush=True)
