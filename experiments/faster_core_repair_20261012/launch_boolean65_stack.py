"""Serial 20-second debugger sample on the verified literal operand packet."""
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT/'python'))
from make_it_stand.verification.serial_source_job import run_serial_source_job


def main():
    capture = HERE/'inertia-boolean65-capture-job-07/output/capture'
    bridge = HERE/'formal-cherchi-build-job-02/output/build/python'
    runtime = HERE/'formal-native-core-build-job-13/output/build/compiler-runtime/libstdc++.so.6'
    files = dict(worker=HERE/'replay_saved_boolean.py', debugger=HERE/'sample_boolean_stack.gdb',
        launcher=Path(__file__), packet=capture/'boolean-evaluation65.json',
        prefix=capture/'prefix-comparison.json', runtime=runtime,
        bridge=bridge/'_cherchi_core_boolean.cpython-313-x86_64-linux-gnu.so')
    import json
    if not json.loads(files['prefix'].read_text())['passed']:
        raise ValueError('CAPTURE_PREFIX_NOT_VERIFIED')
    result = run_serial_source_job(ROOT, HERE/'inertia-boolean65-stack-job-04',
        ['/usr/bin/env', 'LD_PRELOAD='+str(runtime), '/usr/bin/gdb', '-q', '-n', '-batch',
         '-x', str(files['debugger']), '--args', '/home/lin/anaconda3/bin/python', '-B',
         str(files['worker']), str(files['packet']), '{job}/output/boolean-diagnostic-output.json'], files,
        environment={'PYTHONPATH':str(bridge)})
    print(result.status, result.returncode, flush=True)
    if result.status != 'PROCESS_EXITED' or result.returncode != 0:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
