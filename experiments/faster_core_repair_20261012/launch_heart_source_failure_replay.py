"""Serialize the registered source replay after diagnostic build 14 exits."""
from pathlib import Path
import json
import hashlib
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT/'python'))
from make_it_stand.verification.serial_source_job import run_serial_source_job, reconcile_terminated_source_job


def main():
    generation = sys.argv[1] if len(sys.argv)>1 else '14'
    if generation not in ('14','15','16'):
        raise ValueError('UNREGISTERED_SOURCE_REPLAY_BUILD')
    build_job = HERE/('formal-native-core-build-job-'+generation)
    while not (build_job/'result.json').exists():
        if (build_job/'output/receipt.json').exists():
            observation = reconcile_terminated_source_job(ROOT, build_job)
            receipt = json.loads((build_job/'output/receipt.json').read_text())
            module = Path(receipt['module'])
            if observation['service'].get('ExecMainStatus') != '0' or hashlib.sha256(module.read_bytes()).hexdigest() != receipt['module_sha256']:
                raise RuntimeError('DIAGNOSTIC_BUILD_RECEIPT_INVALID')
            for path, sha in receipt['compiled_sources'].items():
                if hashlib.sha256(Path(path).read_bytes()).hexdigest() != sha:
                    raise RuntimeError('DIAGNOSTIC_BUILD_SOURCE_CHANGED')
            break
        time.sleep(10)
    if (build_job/'result.json').exists():
        result = json.loads((build_job/'result.json').read_text())
        if result['status'] != 'PROCESS_EXITED' or result['returncode'] != 0:
            raise RuntimeError('DIAGNOSTIC_BUILD_FAILED_REPLAY_NOT_STARTED')
    build = build_job/'output/build'
    runtime = build/'compiler-runtime/libstdc++.so.6'
    run = HERE/'heart-certified-formal-seed-1-job-01/output/run'
    files = dict(launcher=Path(__file__), worker=HERE/'replay_heart_source_failure.py',
        registration=HERE/'heart-seed1-source-failure-registration-01.json',
        native=build/'python/_make_it_stand.cpython-313-x86_64-linux-gnu.so',
        runtime=runtime, receipt=build_job/'output/receipt.json',
        configuration=run/'configuration.json', events=run/'events.jsonl')
    while True:
        try:
            result = run_serial_source_job(ROOT, HERE/('heart-source-failure-replay-job-01' if generation=='14'
                                                      else 'heart-source-failure-replay-job-02' if generation=='15'
                                                      else 'heart-source-failure-replay-job-03'),
                ['/usr/bin/env', 'LD_PRELOAD='+str(runtime), '/home/lin/anaconda3/bin/python',
                 '-B', str(files['worker'])], files,
                environment={'PYTHONPATH': str(build/'python')})
            break
        except RuntimeError as error:
            if not str(error).startswith(('PROJECT_JOB_BUSY', 'PROJECT_WORKLOAD_ACTIVE:')):
                raise
            time.sleep(10)
    print(result.status, result.returncode, flush=True)


if __name__ == '__main__':
    main()
