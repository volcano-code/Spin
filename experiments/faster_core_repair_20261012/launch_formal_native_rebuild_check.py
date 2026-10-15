from pathlib import Path
import json
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job

build=HERE/'formal-native-core-build-job-04'
if json.loads((build/'result.json').read_text())['status']!='PROCESS_EXITED':
    raise ValueError('FORMAL_NATIVE_BUILD_NOT_SUCCESSFUL_TERMINAL')
native=build/'output/build/python/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
bridge=HERE/'formal-cherchi-build-job-02/output/build/python/_cherchi_core_boolean.cpython-313-x86_64-linux-gnu.so'
arithmetic=HERE/'formal-rational-arithmetic-build-job-01/output/build/python/_cherchi_rational_moments.cpython-313-x86_64-linux-gnu.so'
files=dict(launcher=Path(__file__),worker=HERE/'verify_formal_native_rebuild.py',native=native,bridge=bridge,
    arithmetic=arithmetic,receipt=build/'output/receipt.json',spin=ROOT/'python/make_it_stand/spin.py',
    registered=HERE/'projective-formal-pipeline-job-01/output/registered-inputs.json')
for state in ('ellipsoid','one_sheet','two_sheet'):
    files[state]=HERE/'projective-formal-pipeline-job-01/output/fixed'/(state+'-material.json')
for name in ('cherchi_inprocess','cherchi_motion_differential','spin_problem','mass_moments','rational_text'):
    files[name]=ROOT/f'python/make_it_stand/verification/{name}.py'
result=run_source_job(HERE/'formal-native-rebuild-check-job-01',
    ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
     '/home/lin/anaconda3/bin/python3','-B',str(files['worker']),'{job}/output'],files,
    environment={'PYTHONPATH':':'.join(map(str,[ROOT/'python',native.parent,bridge.parent,arithmetic.parent]))})
print(result.status,result.returncode,flush=True)
