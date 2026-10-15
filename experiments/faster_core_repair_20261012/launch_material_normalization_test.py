from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
build=sys.argv[1];job=sys.argv[2]
module=HERE/build/'output/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
result=run_source_job(HERE/job,
 ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
  '/home/lin/anaconda3/bin/python3','-B',str(HERE/'test_registered_material_normalization.py')],
 dict(launcher=Path(__file__),test=HERE/'test_registered_material_normalization.py',native=module,
      config=HERE/'noodles-registered-core-seed-0-job-01/output/run/mapping-config.json',
      registration=HERE/'noodles-first-core-failures-registration-01/registration.json',
      spin=ROOT/'python/make_it_stand/spin.py'),
 environment={'PYTHONPATH':str(ROOT/'python')+':'+str(module.parent)})
print(result.status,result.returncode,flush=True)
