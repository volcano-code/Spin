from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / 'python'))
from make_it_stand.verification.source_job import run_source_job

files = dict(launcher=Path(__file__), worker=HERE / 'build_validate_prepared_arithmetic.py',
    builder=HERE / 'build_formal_rational_arithmetic.py',
    source=ROOT / 'bindings/cherchi/RationalMoments.cpp',
    cmake=ROOT / 'bindings/cherchi/CMakeLists.txt',
    boolean_source=ROOT / 'bindings/cherchi/InProcessBoolean.cpp',
    tests=HERE / 'test_prepared_rational_boundary.py',
    legacy_tests=ROOT / 'tests/python/test_cherchi_rational_native.py')
for path in (ROOT / 'python/make_it_stand').rglob('*.py'):
    files[str(path)] = path
result = run_source_job(HERE / 'prepared-arithmetic-build-validation-job-01',
    ['/home/lin/anaconda3/bin/python3', '-B', str(files['worker']), '{job}/output'], files)
print(result.status, result.returncode, flush=True)
