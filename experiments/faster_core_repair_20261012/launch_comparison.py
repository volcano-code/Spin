"""Run one bounded, serial same-input diagnostic batch."""
from pathlib import Path
import argparse
import sys
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job

p=argparse.ArgumentParser()
p.add_argument('--native',type=Path,default=ROOT/'build/heart-fixed-facet-20261011-v1/_make_it_stand.cpython-313-x86_64-linux-gnu.so')
p.add_argument('--job',default='same-input-job-04')
args=p.parse_args()
native=args.native.resolve(strict=True)
probe=ROOT/'build/cherchi-reference-singlethread-20260917/cherchi_probe'
worker=HERE/'compare_core_backends.py'
result=run_source_job(HERE/args.job,
    ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
     '/home/lin/anaconda3/bin/python3','-B',str(worker),'--output','{job}/output','--probe',str(probe)],
    dict(worker=worker,timing=HERE/'phase_timing.py',launcher=Path(__file__),native=native,probe=probe,
         spin=ROOT/'python/make_it_stand/spin.py',
         provenance=ROOT/'experiments/spin_it_reproduction/strict/run_faster_reference_boolean.py',
         differential=ROOT/'python/make_it_stand/verification/cherchi_motion_differential.py'),
    environment={'PYTHONPATH':str(ROOT/'python')+':'+str(native.parent)})
print(result.status,result.returncode,flush=True)
