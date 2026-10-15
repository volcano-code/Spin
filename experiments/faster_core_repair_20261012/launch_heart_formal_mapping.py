"""Independent formal Heart preparation only; never starts a case solve."""
from pathlib import Path
import json
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / 'python'))
from make_it_stand.verification.source_job import run_source_job
from inherited_job_identity import inherited_files

base = HERE / 'heart-cherchi-coordinate-map-job-01'
result = json.loads((base / 'result.json').read_text())
if result['status'] != 'PROCESS_EXITED' or result['returncode'] != 0:
    raise ValueError('COORDINATE_MAP_NOT_COMPLETED')
config = json.loads((base / 'config.json').read_text())
files = inherited_files(config['identity'])
files.update(launcher=Path(__file__), worker=HERE / 'prepare_heart_formal_mapping.py',
    identity_helper=HERE / 'inherited_job_identity.py',
    request_comparison=HERE / 'noodles_formal_case.py',
    coordinate_registration=base / 'output/registration.json',
    coordinate_report=base / 'output/report.json',
    clearance=HERE / 'heart-derived-clearance-transfer-01.json',
    clearance_worker=HERE / 'verify_heart_clearance_transfer.py',
    clearance_inequality=HERE / 'clearance_bound_transfer.py')
for state in ('heart-first-cap-failure', 'heart-first-small-face-rejection'):
    files[state] = base / f'output/{state}-normalized.json'
result = run_source_job(HERE / 'heart-formal-mapping-prepare-job-01',
    ['/usr/bin/env', 'LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
     '/home/lin/anaconda3/bin/python3', '-B', str(files['worker']), '{job}/output'],
    files, environment={'PYTHONPATH': config['environment']['PYTHONPATH']})
print(result.status, result.returncode, flush=True)
