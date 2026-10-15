from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / 'python'))
from make_it_stand.verification.source_job import run_source_job

if not (HERE / 'core-benchmark-faces-12288-seed-0-job-01/result.json').exists():
    raise ValueError('PRECEDING_BENCHMARK_NOT_TERMINAL')
generation = sys.argv[1] if len(sys.argv) > 1 else '01'
if generation not in ('01', '02'):
    raise ValueError('UNREGISTERED_ARITHMETIC_GENERATION')
files = dict(launcher=Path(__file__), worker=HERE / 'build_rational_arithmetic.py',
             source=ROOT / 'bindings/cherchi/RationalMoments.cpp',
             tests=ROOT / 'tests/python/test_cherchi_rational_native.py')
result = run_source_job(HERE / ('rational-arithmetic-build-job-' + generation),
    ['/home/lin/anaconda3/bin/python3', '-B', str(files['worker']), '{job}/output'], files)
print(result.status, result.returncode, flush=True)
