from pathlib import Path
import json
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT/'python'))
from make_it_stand.verification.serial_source_job import run_serial_source_job


def main():
    native = HERE/'formal-native-core-build-job-10/output/build/python/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
    original = HERE/'noodles-formal-cherchi-seed-0-job-01'
    saved = HERE/'noodles-terminal-recovery-job-01/output'
    frozen = HERE/'noodles-balanced-request-freeze-job-01/output/request.json'
    files = dict(launcher=Path(__file__), worker=HERE/'compare_noodles_owned_terminal.py',
        native=native, native_receipt=HERE/'formal-native-core-build-job-10/output/receipt.json',
        frozen=frozen, spin=ROOT/'python/make_it_stand/spin.py')
    for path in (ROOT/'python/make_it_stand/verification').glob('*.py'):
        files['module_'+path.name] = path
    for i, name in enumerate(json.loads(frozen.read_text())['input_sources']):
        files['input_'+str(i)] = Path(name)
    for path in (original/'output/run').glob('*.json*'):
        files['original_'+path.name] = path
    files['original_registration'] = original/'output/registration.json'
    for path in saved.glob('*.json'):
        files['saved_'+path.name] = path
    result = run_serial_source_job(ROOT, HERE/'noodles-owned-terminal-comparison-job-01',
        ['/usr/bin/env', 'LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
         '/home/lin/anaconda3/bin/python3', '-B', str(files['worker']), '{job}/output'], files,
        environment={'PYTHONPATH': str(ROOT/'python')+':'+str(native.parent)})
    print(result.status, result.returncode, flush=True)
    if result.status != 'PROCESS_EXITED' or result.returncode != 0:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
