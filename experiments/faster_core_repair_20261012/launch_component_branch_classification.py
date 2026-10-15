from pathlib import Path
import json
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / 'python'))
from make_it_stand.verification.source_job import run_source_job

files = dict(launcher=Path(__file__), worker=HERE / 'classify_component_branch_changes.py',
             registration=HERE / 'noodles-first-core-failures-registration-01/registration.json',
             prior=HERE / 'component-selector-branch-job-01/output/source_side.json',
             replay=HERE / 'cherchi_core_map.py', planar=HERE / 'exact_planar_patches.py',
             formal=ROOT / 'python/make_it_stand/verification/cherchi_inprocess.py',
             motion=ROOT / 'python/make_it_stand/verification/cherchi_motion_differential.py')
registration = json.loads(files['registration'].read_text())
files['config'] = Path(registration['mapping_config'])
for path in (HERE / 'component-selector-derivatives-job-01/output/source_side').glob('*.json'):
    files['saved_' + path.name] = path
result = run_source_job(HERE / 'component-branch-classification-job-01',
    ['/home/lin/anaconda3/bin/python3', '-B', str(files['worker']), '{job}/output'], files,
    environment={'PYTHONPATH': str(ROOT / 'python')})
print(result.status, result.returncode, flush=True)
