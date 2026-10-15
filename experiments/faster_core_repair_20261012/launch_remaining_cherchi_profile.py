from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
native=HERE/'formal-native-core-build-job-06/output/build/python/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
bridge=HERE/'formal-cherchi-build-job-02/output/build/python/_cherchi_core_boolean.cpython-313-x86_64-linux-gnu.so'
arithmetic=HERE/'formal-rational-arithmetic-build-job-02/output/build/python/_cherchi_rational_moments.cpython-313-x86_64-linux-gnu.so'
run=HERE/'layered-formal-benchmark-faces-12288-seed-0-job-04/output/run'
files=dict(launcher=Path(__file__),worker=HERE/'profile_remaining_cherchi_phases.py',native=native,bridge=bridge,
    arithmetic=arithmetic,config=run/'configuration.json',material=run/'evaluation-000073.json',
    inputs=HERE/'core-benchmark-inputs-job-01/output/faces-12288.json',spin=ROOT/'python/make_it_stand/spin.py')
for path in (ROOT/'python/make_it_stand/verification').glob('*.py'):files['module_'+path.name]=path
result=run_source_job(HERE/'remaining-cherchi-phases-profile-job-01',
    ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6','/home/lin/anaconda3/bin/python3',
     '-B',str(files['worker']),'{job}/output'],files,
    environment={'PYTHONPATH':':'.join(map(str,[ROOT/'python',native.parent,bridge.parent,arithmetic.parent]))})
print(result.status,result.returncode,flush=True)
