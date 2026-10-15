from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job

native=HERE/'rational-boundary-binding-build-job-01/output/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
files=dict(launcher=Path(__file__),worker=HERE/'compare_rational_trial_snapshots.py',native=native,
    spin=ROOT/'python/make_it_stand/spin.py')
for name in ('spin_problem','faster_material','spin_common','mass_moments','rational_text'):
    files[name]=ROOT/f'python/make_it_stand/verification/{name}.py'
for state in ('ellipsoid','one_sheet','two_sheet'):
    files[state]=HERE/'projective-formal-pipeline-job-01/output/fixed'/(state+'-material.json')
run=HERE/'core-benchmark-faces-12288-seed-0-job-01/output/run'
for path in run.glob('*.json'):files['run_'+path.name]=path
result=run_source_job(HERE/'rational-trial-snapshots-job-01',
    ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
     '/home/lin/anaconda3/bin/python3','-B',str(files['worker']),'{job}/output'],files,
    environment={'PYTHONPATH':str(ROOT/'python')+':'+str(native.parent)})
print(result.status,result.returncode,flush=True)
