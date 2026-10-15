"""Serial actual-map timing, not a resumed optimization."""
from pathlib import Path
import json
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT/'python'))
from make_it_stand.verification.serial_source_job import run_serial_source_job


def main():
    generation = sys.argv[1] if len(sys.argv) > 1 else '01'
    if generation not in ('01','02','03'):
        raise ValueError('UNREGISTERED_NOODLES_MAP_PROFILE')
    build = HERE/'formal-native-core-build-job-12/output/build'
    native = build/'python/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
    bridge = HERE/'formal-cherchi-build-job-02/output/build/python/_cherchi_core_boolean.cpython-313-x86_64-linux-gnu.so'
    arithmetic_generation = '03' if generation == '03' else '02'
    arithmetic = HERE/f'formal-rational-arithmetic-build-job-{arithmetic_generation}/output/build/python/_cherchi_rational_moments.cpython-313-x86_64-linux-gnu.so'
    runtime = build/'compiler-runtime/libstdc++.so.6'
    run = HERE/'noodles-formal-cherchi-seed-0-job-01/output/run'
    frozen = HERE/'noodles-balanced-request-freeze-job-01/output/request.json'
    files = dict(launcher=Path(__file__), worker=HERE/'profile_noodles_saved_evaluation.py',
        native=native, runtime=runtime, native_receipt=build.parent/'receipt.json',
        bridge=bridge, arithmetic=arithmetic, config=run/'configuration.json', events=run/'events.jsonl',
        frozen=frozen, independent_integral=HERE/'independent_tetra_moments.py',
        spin=ROOT/'python/make_it_stand/spin.py')
    for path in run.glob('evaluation-*.json'):
        files['checkpoint_'+path.name] = path
    for index, name in enumerate(json.loads(frozen.read_text())['input_sources']):
        files['input_'+str(index)] = Path(name)
    for path in (ROOT/'python/make_it_stand/verification').glob('*.py'):
        files['module_'+path.name] = path
    result = run_serial_source_job(ROOT, HERE/('noodles-saved-evaluation-profile-job-'+generation),
        ['/usr/bin/env', 'LD_PRELOAD='+str(runtime), '/home/lin/anaconda3/bin/python3', '-B',
         str(files['worker']), '{job}/output'], files,
        environment={'PYTHONPATH': ':'.join(map(str, [ROOT/'python', native.parent,
            bridge.parent, arithmetic.parent, HERE]))})
    print(result.status, result.returncode, flush=True)
    if result.status != 'PROCESS_EXITED' or result.returncode != 0:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
