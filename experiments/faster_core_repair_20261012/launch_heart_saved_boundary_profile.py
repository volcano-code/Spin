from pathlib import Path
import json
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT/'python'))
from make_it_stand.verification.serial_source_job import run_serial_source_job


def main():
    generation = sys.argv[1] if len(sys.argv) > 1 else '01'
    if generation not in ('01', '02', '03', '04'):
        raise ValueError('UNREGISTERED_HEART_BOUNDARY_PROFILE')
    build = HERE/('formal-native-core-build-job-12' if generation == '04' else 'formal-native-core-build-job-10')
    result = json.loads((build/'result.json').read_text())
    if result['status'] != 'PROCESS_EXITED' or result['returncode'] != 0:
        raise ValueError('NATIVE_BUILD_NOT_COMPLETED')
    receipt = json.loads((build/'output/receipt.json').read_text())
    run = HERE/'heart-guarded-formal-seed-0-job-01/output/run'
    frozen = HERE.parent/'faster_heart_contact_20261011/rounded-clearance-v3-faster-frozen-v2/request.json'
    files = dict(launcher=Path(__file__), worker=HERE/'profile_heart_saved_boundary.py',
        native=Path(receipt['module']), receipt=build/'output/receipt.json',
        config=run/'configuration.json', events=run/'events.jsonl',
        material=run/('evaluation-000002.json' if generation == '04' else 'evaluation-000045.json'),
        frozen=frozen, spin_problem=ROOT/'python/make_it_stand/verification/spin_problem.py',
        rational_text=ROOT/'python/make_it_stand/verification/rational_text.py')
    for name in ('faster_material', 'mass_moments', 'spin_common'):
        files[name] = ROOT/('python/make_it_stand/verification/'+name+'.py')
    if generation == '04':
        for path in (ROOT/'python/make_it_stand/verification').glob('*.py'):
            files['module_'+path.name] = path
    for i, path in enumerate(json.loads(frozen.read_text())['input_sources']):
        files['frozen_input_'+str(i)] = Path(path)
    result = run_serial_source_job(ROOT, HERE/('heart-stalled-boundary-profile-job-'+generation),
        ['/usr/bin/env', 'LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
         '/home/lin/anaconda3/bin/python3', '-B', str(files['worker']), '{job}/output',
         'classify' if generation == '04' else 'compare' if generation == '03' else 'phases'], files,
        environment={'PYTHONPATH': str(ROOT/'python')+':'+str(files['native'].parent)})
    print(result.status, result.returncode, flush=True)
    if result.status != 'PROCESS_EXITED' or result.returncode != 0:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
