from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / 'python'))
from make_it_stand.verification.source_job import run_source_job

native = HERE / 'single-transport-build-job-01/output/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
engine = HERE / 'formal-rational-arithmetic-build-job-01/output/build/python/_cherchi_rational_moments.cpython-313-x86_64-linux-gnu.so'
files = dict(launcher=Path(__file__), worker=HERE / 'verify_native_saved_derivatives.py', native=native, engine=engine,
    registered=HERE / 'projective-formal-pipeline-job-01/output/registered-inputs.json',
    spin=ROOT / 'python/make_it_stand/spin.py', formal=ROOT / 'python/make_it_stand/verification/cherchi_inprocess.py',
    motion=ROOT / 'python/make_it_stand/verification/cherchi_motion_differential.py')
for state in ('ellipsoid', 'one_sheet', 'two_sheet'):
    files[state + '_certificate'] = HERE / 'projective-fixed-branch-job-01/output' / (state + '.json')
    for path in (HERE / 'projective-formal-pipeline-job-01/output/derivatives' / state).glob('*.json'):
        files[state + '_' + path.name] = path
result = run_source_job(HERE / 'native-saved-derivatives-job-01',
    ['/usr/bin/env', 'LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
     '/home/lin/anaconda3/bin/python3', '-B', str(files['worker']), '{job}/output'], files,
    environment={'PYTHONPATH': str(ROOT / 'python') + ':' + str(native.parent) + ':' + str(engine.parent)})
print(result.status, result.returncode, flush=True)
