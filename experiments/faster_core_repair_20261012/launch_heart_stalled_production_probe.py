"""Serial resource job for the registered stalled-state production probe."""
from pathlib import Path
import json
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT/'python'))
from make_it_stand.verification.serial_source_job import run_serial_source_job


def main():
    build = HERE/'formal-native-core-build-job-12/output/build'
    native = build/'python/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
    bridge = HERE/'formal-cherchi-build-job-02/output/build/python/_cherchi_core_boolean.cpython-313-x86_64-linux-gnu.so'
    arithmetic = HERE/'formal-rational-arithmetic-build-job-02/output/build/python/_cherchi_rational_moments.cpython-313-x86_64-linux-gnu.so'
    runtime = build/'compiler-runtime/libstdc++.so.6'
    run = HERE/'heart-guarded-formal-seed-0-job-01/output/run'
    frozen = HERE.parent/'faster_heart_contact_20261011/rounded-clearance-v3-faster-frozen-v2/request.json'
    files = dict(launcher=Path(__file__), worker=HERE/'probe_heart_stalled_production.py',
        native=native, runtime=runtime, native_receipt=build.parent/'receipt.json', bridge=bridge,
        arithmetic=arithmetic, configuration=run/'configuration.json', events=run/'events.jsonl',
        original_material=run/'evaluation-000045.json', frozen=frozen,
        independent_integral=HERE/'independent_tetra_moments.py',
        registered_weights=HERE/'probe_formal_cherchi_derivatives.py', spin=ROOT/'python/make_it_stand/spin.py')
    for i, path in enumerate(json.loads(frozen.read_text())['input_sources']):
        files['input_'+str(i)] = Path(path)
    for path in (ROOT/'python/make_it_stand/verification').glob('*.py'):
        files['module_'+path.name] = path
    result = run_serial_source_job(ROOT, HERE/'heart-stalled-production-81-job-01',
        ['/usr/bin/env', 'LD_PRELOAD='+str(runtime), '/home/lin/anaconda3/bin/python3', '-B',
         str(files['worker']), '{job}/output'], files,
        environment={'PYTHONPATH': ':'.join(map(str, [ROOT/'python', native.parent,
            bridge.parent, arithmetic.parent, HERE]))})
    print(result.status, result.returncode, flush=True)
    if result.status != 'PROCESS_EXITED' or result.returncode != 0:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
