"""One serial projection audit after the actual input controller terminates."""
import argparse
import json
from pathlib import Path
import sys
import time

from launch_rounded_admission import ROOT, process_start
from rounded_input import publish
from make_it_stand.verification.source_job import run_source_job


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--wait-pid', type=int, required=True)
    args = parser.parse_args()
    case = args.input.resolve()
    start = process_start(args.wait_pid)
    if start is None:
        raise ValueError('EXPECTED_LIVE_INPUT_CONTROLLER')
    publish(case/'projection-schedule.json', dict(schema='serial-rounded-projection-schedule-v1',
        wait_pid=args.wait_pid, process_start=start, automatic_retries=0, optimization_invocations=0))
    print(json.dumps(dict(status='WAITING_FOR_ACTUAL_INPUT_CONTROLLER', pid=args.wait_pid)), flush=True)
    while process_start(args.wait_pid) == start:
        time.sleep(1)
    receipt = case/'admission-controller-terminal.json'
    if not receipt.exists() or json.loads(receipt.read_text())['returncode'] != 0:
        publish(case/'projection-controller-terminal.json', dict(status='NOT_RUN_INPUT_JOB_INCOMPLETE'))
        return 1
    source = case/'input-admission-v1/output'
    if not (source/'admission.json').exists():
        publish(case/'projection-controller-terminal.json', dict(status='NOT_RUN_INPUT_RECORD_MISSING'))
        return 1
    native = ROOT/'build/spin2014-objectives-quadric-20261007/python/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
    worker = Path(__file__).with_name('inspect_rounded_projection.py')
    files = dict(controller=Path(__file__), worker=worker,
        geometry_module=Path(__file__).with_name('rounded_input.py'),
        conversion_module=Path(__file__).with_name('admit_rounded_input.py'),
        admission=source/'admission.json', native=native,
        input_receipt=case/'input-admission-v1/result.json')
    for name in ('outer', 'design', 'fixed'):
        files[name] = source/(name+'.rational.json')
    result = run_source_job(case/'projection-audit-v1', ['/usr/bin/env',
        'LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
        '/home/lin/anaconda3/bin/python3', '-B', str(worker),
        '--admission', str(source), '--output', '{job}/output'], files,
        environment={'PYTHONPATH':str(ROOT/'python')+':'+str(native.parent)})
    publish(case/'projection-controller-terminal.json', dict(status=result.status,
        returncode=result.returncode, scope='PROJECTION_AUDIT_NOT_OPTIMIZATION_OR_SIMULATION',
        job=str(result.directory), unit=result.unit))
    return 0 if result.returncode == 0 else 1


if __name__ == '__main__':
    raise SystemExit(main())
