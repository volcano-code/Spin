from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
module=HERE/'core-failure-tape-build-job-01/output/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
result=run_source_job(HERE/'core-registered-failure-capture-job-01',
 ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
  '/home/lin/anaconda3/bin/python3','-B',str(HERE/'capture_registered_failures.py'),'{job}/output'],
 dict(launcher=Path(__file__),worker=HERE/'capture_registered_failures.py',native=module,
      registration=HERE/'noodles-first-core-failures-registration-01/registration.json',
      config=HERE/'noodles-registered-core-seed-0-job-01/output/run/mapping-config.json',
      spin=ROOT/'python/make_it_stand/spin.py'),
 environment={'PYTHONPATH':str(ROOT/'python')+':'+str(module.parent)})
print(result.status,result.returncode,flush=True)
