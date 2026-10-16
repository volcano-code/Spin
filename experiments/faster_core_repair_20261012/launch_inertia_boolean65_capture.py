"""Resource-limited, serial admission for the actual seed replay capture."""
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT/'python'))
from make_it_stand.verification.serial_source_job import run_serial_source_job


def main():
    build = HERE/'formal-native-core-build-job-13/output/build'
    bridge = HERE/'formal-cherchi-build-job-02/output/build/python'
    arithmetic = HERE/'formal-rational-arithmetic-build-job-03/output/build/python'
    source = HERE/'spin2014-shared-small-top-inertia-only-seed-1-run-01'
    files = dict(worker=HERE/'replay_inertia_boolean65.py', launcher=Path(__file__),
        capture=HERE/'boolean_operand_capture.py', benchmark=HERE/'run_core_benchmark.py',
        native=build/'python/_make_it_stand.cpython-313-x86_64-linux-gnu.so',
        runtime=build/'compiler-runtime/libstdc++.so.6', receipt=build.parent/'receipt.json',
        spin=ROOT/'python/make_it_stand/spin.py')
    for directory in (ROOT/'python/make_it_stand/verification', HERE/'core-benchmark-inputs-job-01/output',
                      source/'output/run', bridge, arithmetic):
        for p in directory.iterdir():
            if p.is_file(): files[str(p.relative_to(ROOT))] = p
    historical = HERE/'historical-inertia-python-01/python'
    files['historical_spin'] = historical/'make_it_stand/spin.py'
    import hashlib, json
    expected = json.loads((source/'config.json').read_text())['identity']
    for key, path in [('spin', files['historical_spin']),
                      ('python/make_it_stand/verification/cherchi_inprocess.py', historical/'make_it_stand/verification/cherchi_inprocess.py'),
                      ('worker', historical.parent/'run_core_benchmark.py')]:
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected[key]['sha256']:
            raise ValueError('HISTORICAL_SOURCE_MISMATCH:'+key)
        files['historical_'+key] = path
    result = run_serial_source_job(ROOT, HERE/'inertia-boolean65-capture-job-07',
        ['/usr/bin/env', 'LD_PRELOAD='+str(files['runtime']), '/home/lin/anaconda3/bin/python', '-B',
         str(files['worker']), '{job}/output', str(source)], files,
        environment={'PYTHONPATH': ':'.join(map(str, [historical, build/'python', bridge, arithmetic, HERE]))})
    print(result.status, result.returncode, flush=True)
    if result.status != 'PROCESS_EXITED' or result.returncode != 0:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
