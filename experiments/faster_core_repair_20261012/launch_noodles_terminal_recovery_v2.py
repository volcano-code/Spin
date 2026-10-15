"""Recover the saved noodles terminal with the explicit classified/owned path."""
from pathlib import Path
import json
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT/'python'))
from make_it_stand.verification.serial_source_job import run_serial_source_job


def main():
    source = HERE/'noodles-formal-cherchi-seed-0-job-01'
    native = HERE/'formal-native-core-build-job-12/output/build/python/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
    bridge = HERE/'formal-cherchi-build-job-02/output/build/python/_cherchi_core_boolean.cpython-313-x86_64-linux-gnu.so'
    arithmetic = HERE/'formal-rational-arithmetic-build-job-02/output/build/python/_cherchi_rational_moments.cpython-313-x86_64-linux-gnu.so'
    frozen = HERE/'noodles-balanced-request-freeze-job-01/output/request.json'
    files = dict(launcher=Path(__file__), worker=HERE/'run_noodles_terminal_recovery.py',
        native=native, native_receipt=HERE/'formal-native-core-build-job-12/output/receipt.json',
        bridge=bridge, arithmetic=arithmetic, frozen_request=frozen,
        source_job_configuration=source/'config.json', source_kernel=source/'kernel.json',
        source_registration=source/'output/registration.json',
        source_configuration=source/'output/run/configuration.json', source_events=source/'output/run/events.jsonl',
        preparation=HERE/'noodles_formal_case.py', symbolic=HERE/'independent_eq14.py',
        spin=ROOT/'python/make_it_stand/spin.py')
    for path in (source/'output/run').glob('evaluation-*.json'):
        files['accepted_prefix_'+path.name] = path
    for index, path in enumerate(json.loads(frozen.read_text())['input_sources']):
        files['input_'+str(index)] = Path(path)
    for path in (ROOT/'python/make_it_stand/verification').glob('*.py'):
        files['module_'+path.name] = path
    result = run_serial_source_job(ROOT, HERE/'noodles-terminal-recovery-job-02',
        ['/usr/bin/env', 'LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
         '/home/lin/anaconda3/bin/python3', '-B', str(files['worker']), '{job}/output', str(source),
         '--classified-boundary', '--reuse-input-proof', '--owned-terminal'], files,
        environment={'PYTHONPATH':':'.join(map(str, [ROOT/'python', native.parent, bridge.parent,
                                                     arithmetic.parent, HERE]))})
    print(result.status, result.returncode, flush=True)
    if result.status != 'PROCESS_EXITED' or result.returncode != 0:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
