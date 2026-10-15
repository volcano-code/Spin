"""Reuse registered samples only after exact source replay on the formal build."""
from pathlib import Path
import json
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job

if not json.loads((HERE/'formal-native-rebuild-check-job-01/output/report.json').read_text())['passed']:
    raise ValueError('FORMAL_FIXED_MAP_CHECK_NOT_PASSED')
build=HERE/'formal-native-core-build-job-04'
native=build/'output/build/python/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
arithmetic=HERE/'formal-rational-arithmetic-build-job-01/output/build/python/_cherchi_rational_moments.cpython-313-x86_64-linux-gnu.so'
files=dict(launcher=Path(__file__),worker=HERE/'verify_native_saved_derivatives.py',native=native,
    native_receipt=build/'output/receipt.json',arithmetic=arithmetic,
    registered=HERE/'projective-formal-pipeline-job-01/output/registered-inputs.json',
    spin=ROOT/'python/make_it_stand/spin.py')
for name in ('cherchi_inprocess','cherchi_motion_differential','spin_problem','mass_moments','rational_text'):
    files[name]=ROOT/f'python/make_it_stand/verification/{name}.py'
for state in ('ellipsoid','one_sheet','two_sheet'):
    files['certificate_'+state]=HERE/'projective-fixed-branch-job-01/output'/(state+'.json')
    folder=HERE/'projective-formal-pipeline-job-01/output/derivatives'/state
    files['manifest_'+state]=folder/'manifest.json'
    for path in folder.glob('point-*-material.json'):
        files[state+'_'+path.name]=path
result=run_source_job(HERE/'formal-native-saved-derivatives-job-01',
    ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
     '/home/lin/anaconda3/bin/python3','-B',str(files['worker']),'{job}/output'],files,
    environment={'PYTHONPATH':':'.join(map(str,[ROOT/'python',native.parent,arithmetic.parent]))})
print(result.status,result.returncode,flush=True)
