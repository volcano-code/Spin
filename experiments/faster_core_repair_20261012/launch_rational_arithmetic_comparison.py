from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / 'python'))
from make_it_stand.verification.source_job import run_source_job

generation = sys.argv[1] if len(sys.argv) > 1 else '01'
if generation not in ('01', '02'):
    raise ValueError('UNREGISTERED_ARITHMETIC_GENERATION')
native = HERE / 'single-transport-build-job-01/output/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
engine = HERE / ('rational-arithmetic-build-job-' + generation) / 'output/_cherchi_rational_moments.cpython-313-x86_64-linux-gnu.so'
files = dict(launcher=Path(__file__), worker=HERE / 'compare_rational_arithmetic.py', native=native, engine=engine,
    config=HERE / 'projective-formal-pipeline-job-01/output/registered-inputs.json',
    physical=HERE / 'core-benchmark-inputs-job-01/output/faces-12288.json',
    spin=ROOT / 'python/make_it_stand/spin.py', formal=ROOT / 'python/make_it_stand/verification/cherchi_inprocess.py',
    motion=ROOT / 'python/make_it_stand/verification/cherchi_motion_differential.py')
for state in ('ellipsoid', 'one_sheet', 'two_sheet'):
    files[state] = HERE / 'projective-formal-pipeline-job-01/output/fixed' / (state + '-material.json')
run = HERE / 'core-benchmark-faces-12288-seed-0-job-01/output/run'
for path in run.glob('*.json'):
    files['run_' + path.name] = path
result = run_source_job(HERE / ('rational-arithmetic-comparison-job-' + generation),
    ['/usr/bin/env', 'LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
     '/home/lin/anaconda3/bin/python3', '-B', str(files['worker']), '{job}/output'], files,
    environment={'PYTHONPATH': str(ROOT / 'python') + ':' + str(native.parent) + ':' + str(engine.parent)})
print(result.status, result.returncode, flush=True)
