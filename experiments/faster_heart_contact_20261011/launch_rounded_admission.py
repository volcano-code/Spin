"""One bounded input admission after the existing queue/aggregator have exited.

This launches geometry admission only. It never launches optimization or grants
a production derivative verdict. Process start times prevent PID reuse waits.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
from rounded_input import publish


def process_start(pid):
    try:
        return Path(f'/proc/{pid}/stat').read_text().rsplit(')', 1)[1].split()[19]
    except FileNotFoundError:
        return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--wait-pid', type=int, action='append', default=[])
    args = parser.parse_args()
    case = args.input.resolve()
    job = case/'input-admission-v1'
    if job.exists():
        raise FileExistsError(job)
    native = ROOT/'build/spin2014-objectives-quadric-20261007/python/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
    waiting = {pid: process_start(pid) for pid in args.wait_pid}
    waiting = {pid: start for pid, start in waiting.items() if start is not None}
    module = Path(__file__).with_name('rounded_input.py')
    worker = Path(__file__).with_name('admit_rounded_input.py')
    files = dict(controller=Path(__file__), geometry_module=module, worker=worker,
        spec=case/'spec.json', preparation=case/'preparation.json', native=native,
        legacy_heart_generator=ROOT/'experiments/reconstructed_heart_top/heart_top.py')
    for name in ('heart-body.obj', 'inner-source.obj', 'pin.obj'):
        files[name] = case/name
    frozen = {name: hashlib.sha256(Path(path).read_bytes()).hexdigest() for name, path in files.items()}
    publish(case/'admission-schedule.json', dict(schema='rounded-heart-serial-admission-schedule-v1',
        waiting_processes=waiting, input_hashes=frozen, automatic_retries=0,
        optimization_invocations=0, simulation_invocations=0))
    print(json.dumps(dict(status='WAITING_FOR_EXISTING_QUEUE_AND_AGGREGATOR', pids=waiting)), flush=True)
    while any(process_start(pid) == start for pid, start in waiting.items()):
        time.sleep(1)
    for name, path in files.items():
        if hashlib.sha256(Path(path).read_bytes()).hexdigest() != frozen[name]:
            raise ValueError('SOURCE_CHANGED_WHILE_WAITING: '+name)
    result = run_source_job(job, ['/usr/bin/env',
        'LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
        '/home/lin/anaconda3/bin/python3', '-B', str(worker),
        '--input', str(case), '--output', '{job}/output'], files,
        environment={'PYTHONPATH':str(ROOT/'python')+':'+str(native.parent)})
    publish(case/'admission-controller-terminal.json', dict(status=result.status,
        returncode=result.returncode, job=str(result.directory), unit=result.unit,
        scope='INPUT_ADMISSION_ONLY_NOT_OPTIMIZATION'))
    print(json.dumps(dict(status=result.status, returncode=result.returncode)), flush=True)
    return 0 if result.returncode == 0 else 1


if __name__ == '__main__':
    raise SystemExit(main())
