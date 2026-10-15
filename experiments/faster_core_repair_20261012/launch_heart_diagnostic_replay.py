from pathlib import Path
import json
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT/'python'))
from make_it_stand.verification.serial_source_job import run_serial_source_job


def main():
    run = HERE/'heart-guarded-formal-seed-0-job-01/output/run'
    native = HERE/'formal-native-core-build-job-10/output/build/python/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
    bridge = HERE/'formal-cherchi-build-job-02/output/build/python/_cherchi_core_boolean.cpython-313-x86_64-linux-gnu.so'
    arithmetic = HERE/'formal-rational-arithmetic-build-job-02/output/build/python/_cherchi_rational_moments.cpython-313-x86_64-linux-gnu.so'
    frozen = HERE.parent/'faster_heart_contact_20261011/rounded-clearance-v3-faster-frozen-v2/request.json'
    diagnostic = HERE/'heart-stalled-prefix-independent-algebra-02.json'
    files = dict(launcher=Path(__file__), worker=HERE/'replay_heart_diagnostic_states.py',
        native=native, bridge=bridge, arithmetic=arithmetic, frozen=frozen, diagnostic=diagnostic,
        configuration=run/'configuration.json', spin=ROOT/'python/make_it_stand/spin.py',
        native_receipt=HERE/'formal-native-core-build-job-10/output/receipt.json')
    for i, path in enumerate(json.loads(frozen.read_text())['input_sources']):
        files['input_'+str(i)] = Path(path)
    for row in json.loads(diagnostic.read_text())['diagnostic_states']:
        files['state_'+str(row['evaluation'])] = run/row['checkpoint']['evaluation_file']
    for path in (ROOT/'python/make_it_stand/verification').glob('*.py'):
        files['module_'+path.name] = path
    result = run_serial_source_job(ROOT, HERE/'heart-four-state-production-replay-job-01',
        ['/usr/bin/env', 'LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
         '/home/lin/anaconda3/bin/python3', '-B', str(files['worker']), '{job}/output'], files,
        environment={'PYTHONPATH': ':'.join(map(str, [ROOT/'python', native.parent, bridge.parent, arithmetic.parent]))})
    print(result.status, result.returncode, flush=True)
    if result.status != 'PROCESS_EXITED' or result.returncode != 0:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
