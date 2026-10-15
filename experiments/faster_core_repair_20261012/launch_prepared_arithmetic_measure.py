from pathlib import Path
import json
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / 'python'))
from make_it_stand.verification.source_job import run_source_job

build = HERE / 'prepared-arithmetic-build-validation-job-01'
resource = json.loads((build / 'result.json').read_text())
report = json.loads((build / 'output/report.json').read_text())
if resource['status'] != 'PROCESS_EXITED' or resource['returncode'] != 0 or not report['passed']:
    raise ValueError('CANDIDATE_ARITHMETIC_TESTS_NOT_PASSED')
candidate = build / 'output/candidate/build/python/_cherchi_rational_moments.cpython-313-x86_64-linux-gnu.so'
baseline = HERE / 'formal-rational-arithmetic-build-job-02/output/build/python/_cherchi_rational_moments.cpython-313-x86_64-linux-gnu.so'
files = dict(launcher=Path(__file__), worker=HERE / 'measure_prepared_arithmetic.py',
             candidate=candidate, baseline=baseline, build_result=build / 'result.json',
             validation=build / 'output/report.json', receipt=build / 'output/candidate/receipt.json')
for label in ('small', 'faces-12288'):
    files[label] = HERE / f'layered-formal-benchmark-{label}-seed-0-job-06/output/run/evaluation-000001.json'
result = run_source_job(HERE / 'prepared-arithmetic-comparison-job-01',
    ['/usr/bin/env', 'LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
     '/home/lin/anaconda3/bin/python3', '-B', str(files['worker']), '{job}/output'],
    files, environment={'PYTHONPATH': str(candidate.parent)})
print(result.status, result.returncode, flush=True)
