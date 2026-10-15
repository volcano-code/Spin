from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job

if not (HERE/'core-benchmark-faces-12288-seed-0-checkpoint-boundary-gmp-job-01/result.json').exists():
    raise ValueError('PRECEDING_BENCHMARK_NOT_TERMINAL')
native=HERE/'rational-boundary-binding-build-job-01/output/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
generation=sys.argv[1] if len(sys.argv)>1 else '01'
if generation not in ('01','02'):raise ValueError('UNREGISTERED_PROBE_GENERATION')
files=dict(launcher=Path(__file__),worker=HERE/'probe_physical_boundary_costs.py',native=native,
    spin=ROOT/'python/make_it_stand/spin.py',rational=ROOT/'python/make_it_stand/verification/rational_text.py',
    parser=ROOT/'python/make_it_stand/verification/cherchi_inprocess.py')
run=HERE/'core-benchmark-faces-12288-seed-0-job-01/output/run'
for path in run.glob('*.json'):files['saved_'+path.name]=path
result=run_source_job(HERE/('physical-boundary-costs-job-'+generation),
    ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
     '/home/lin/anaconda3/bin/python3','-B',str(files['worker']),'{job}/output'],files,
    environment={'PYTHONPATH':str(ROOT/'python')+':'+str(native.parent)})
print(result.status,result.returncode,flush=True)
