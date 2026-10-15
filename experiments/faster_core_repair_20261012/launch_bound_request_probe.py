from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
heart=ROOT/'experiments/faster_heart_contact_20261011'
native=HERE/'exact-core-material-build-job-01/output/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
files=dict(launcher=Path(__file__),worker=HERE/'freeze_balanced_noodles_request.py',
    normalizer=HERE/'bound_faster_input.py',integration=heart/'rounded_input.py',
    native=native,spin=ROOT/'python/make_it_stand/spin.py',
    manifest=HERE/'noodles-balanced-shared-export-job-02/output/manifest.json',
    original_manifest=ROOT/'experiments/quadric_physics_demo_20261011/noodles-dense-v5/manifest.json')
for name in ('outer','design','fixed'):
    files[name]=HERE/f'noodles-balanced-shared-export-job-02/output/{name}-si-candidate.obj'
environment={'PYTHONPATH':str(ROOT/'python')+':'+str(native.parent)+':'+str(HERE)}
freeze=HERE/'noodles-balanced-request-freeze-job-01'
result=run_source_job(freeze,
    ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
     '/home/lin/anaconda3/bin/python3','-B',str(files['worker']),'{job}/output'],files,environment=environment)
if result.returncode!=0:raise RuntimeError('REQUEST_FREEZE_FAILED_NO_PROBE')
files.update(worker=HERE/'probe_bound_noodles_seed0.py',
    request=freeze/'output/request.json',journal=heart/'rounded_solver_journal.py',
    source_motion=ROOT/'src/spin/QuadricSourceMotion.cpp')
result=run_source_job(HERE/'noodles-bound-seed0-production-probe-job-01',
    ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
     '/home/lin/anaconda3/bin/python3','-B',str(files['worker']),str(files['request']),'{job}/output'],
    files,environment=environment)
print(result.status,result.returncode,flush=True)
