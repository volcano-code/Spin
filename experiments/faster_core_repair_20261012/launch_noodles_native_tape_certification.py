from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
native=HERE/'cgal-proof-capture-build-job-01/output/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
files=dict(launcher=Path(__file__),worker=HERE/'certify_noodles_native_tapes.py',native=native,
    spin=ROOT/'python/make_it_stand/spin.py',recipes=HERE/'birth_face_recipes.py',chain=HERE/'exact_planar_patches.py',serializer=HERE/'run_small_inprocess_solver.py',
    request=HERE/'noodles-balanced-request-freeze-job-01/output/request.json',
    manifest=HERE/'noodles-algebraic-production-moments-job-01/output/manifest.json',
    numerical_report=HERE/'noodles-algebraic-production-moments-job-01/output/report.json',
    triage=HERE/'noodles-saved-lift-triage-job-01/output/report.json',
    base_recipes=HERE/'noodles-offline-proof-pair-job-01/output/point-000-recipes.json',
    pair_report=HERE/'noodles-offline-proof-pair-job-01/output/report.json')
for i in range(81):
    files[f'source{i}']=HERE/f'noodles-interval-algebraic-source-derivative-job-01/output/point-{i:03d}.json'
    files[f'material{i}']=HERE/f'noodles-algebraic-production-moments-job-01/output/point-{i:03d}-material.json'
for i in (0,1):files[f'capture{i}']=HERE/f'noodles-same-call-proof-capture-job-01/output/point-{i:03d}.json'
result=run_source_job(HERE/'noodles-native-tape-certification-job-01',
    ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6','/home/lin/anaconda3/bin/python3','-B',str(files['worker']),'{job}/output'],files,
    environment={'PYTHONPATH':str(ROOT/'python')+':'+str(native.parent)+':'+str(HERE)})
print(result.status,result.returncode,flush=True)
