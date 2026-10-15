from pathlib import Path
import json
import sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
case=ROOT/'experiments/quadric_physics_demo_20261011/noodles-dense-v5'
manifest=case/'manifest.json';data=json.loads(manifest.read_text())
files=dict(launcher=Path(__file__),worker=HERE/'pilot_balanced_noodles_input.py',
    balanced=HERE/'balanced_exact_input.py',
    representation=case.parent/'core_boolean_mesh.py',manifest=manifest,
    executable=HERE/'exact-input-union-build-job-01/output/exact-input-boolean',
    build_receipt=HERE/'exact-input-union-build-job-01/output/receipt.json')
for name in ['inner-ellipsoid-candidate','outer-ellipsoid',*data['construction']['fixed_parts'][:8]]:
    files[name]=case/data['parts'][name]['model_obj']
result=run_source_job(HERE/'noodles-balanced-eight-part-pilot-job-01',
    ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
     '/home/lin/anaconda3/bin/python3','-B',str(files['worker']),'{job}/output'],files)
print(result.status,result.returncode,flush=True)
