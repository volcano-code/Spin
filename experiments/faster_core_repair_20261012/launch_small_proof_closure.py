from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
native=HERE/'cgal-proof-capture-build-job-01/output/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
files=dict(launcher=Path(__file__),worker=HERE/'run_small_proof_closure.py',native=native,
    spin=ROOT/'python/make_it_stand/spin.py',journal=HERE.parent/'faster_heart_contact_20261011/rounded_solver_journal.py',
    independent=HERE/'independent_tetra_moments.py',
    configuration=HERE/'small-top-inprocess-solver-sequence-01/seed-1-job-01/output/configuration.json',
    fixed_fixture=HERE/'same-input-job-03/output/frozen-inputs.json')
for name in ('spin_problem','spin_bundle','spin_common','rational_text'):
    files[name]=ROOT/f'python/make_it_stand/verification/{name}.py'
result=run_source_job(HERE/'small-cgal-proof-closure-job-01',
    ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6','/home/lin/anaconda3/bin/python3','-B',str(files['worker']),'{job}/output'],files,
    environment={'PYTHONPATH':str(ROOT/'python')+':'+str(native.parent)+':'+str(HERE)})
print(result.status,result.returncode,flush=True)
