from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / 'python'))
from make_it_stand.verification.source_job import run_source_job

files = dict(launcher=Path(__file__), worker=HERE / 'build_formal_rational_arithmetic.py',
    source=ROOT / 'bindings/cherchi/RationalMoments.cpp',
    cmake=ROOT / 'bindings/cherchi/CMakeLists.txt',
    boolean_source=ROOT / 'bindings/cherchi/InProcessBoolean.cpp')
generation = sys.argv[1] if len(sys.argv)>1 else '01'
if generation not in ('01','02'): raise ValueError('UNREGISTERED_FORMAL_ARITHMETIC_BUILD')
result = run_source_job(HERE / ('formal-rational-arithmetic-build-job-'+generation),
    ['/home/lin/anaconda3/bin/python3', '-B', str(files['worker']), '{job}/output'], files)
print(result.status, result.returncode, flush=True)
