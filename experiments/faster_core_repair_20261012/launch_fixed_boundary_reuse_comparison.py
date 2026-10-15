from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job

native=HERE/'formal-native-core-build-job-06/output/build/python/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
run=HERE/'layered-formal-benchmark-faces-12288-seed-0-job-01/output/run'
files=dict(launcher=Path(__file__),worker=HERE/'compare_fixed_boundary_reuse.py',native=native,
    receipt=HERE/'formal-native-core-build-job-06/output/receipt.json',
    input=HERE/'core-benchmark-inputs-job-01/output/faces-12288.json',config=run/'configuration.json',
    first=run/'evaluation-000001.json',last=run/'evaluation-000073.json',
    spin=ROOT/'python/make_it_stand/spin.py')
result=run_source_job(HERE/'fixed-boundary-reuse-comparison-job-01',
    ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
     '/home/lin/anaconda3/bin/python3','-B',str(files['worker']),'{job}/output'],files,
    environment={'PYTHONPATH':str(ROOT/'python')+':'+str(native.parent)})
print(result.status,result.returncode,flush=True)
