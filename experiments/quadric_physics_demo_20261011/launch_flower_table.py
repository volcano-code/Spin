"""Insert one bounded table export between Heart points; resume its coordinator."""
import json
import os
from pathlib import Path
import signal
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job


def main():
    bundle = HERE/'paper-flower-method-v2/table-scene-v1'
    current = ROOT/'experiments/faster_heart_contact_20261011/rounded-clearance-v3-faster-production-v3/point-job-03'
    coordinator = 1766941
    def checked_signal(value):
        path = Path(f'/proc/{coordinator}/cmdline')
        if path.exists() and b'rounded_production_sequence.py' in path.read_bytes():
            os.kill(coordinator, value)
    try:
        deadline = time.monotonic()+3660
        print('WAIT_CURRENT_HEART_POINT', str(current), flush=True)
        while not (current/'complete.json').exists():
            if time.monotonic() > deadline:
                raise RuntimeError('HEART_POINT_DID_NOT_EXIT')
            time.sleep(2)
        print('HEART_POINT_EXITED_START_TABLE_JOB', flush=True)
        files = {p.name: p for p in bundle.iterdir() if p.is_file()}
        files.update({str(p.relative_to(bundle)): p for p in (bundle/'parts').iterdir() if p.is_file()})
        files['launcher'] = Path(__file__)
        result = run_source_job(bundle/'physics-job-02',
            ['/home/lin/anaconda3/bin/python3', '-B', str(bundle/'run_table_scene.py'),
             '--bundle', str(bundle), '--output', '{job}/output'], files)
        print(json.dumps(dict(status=result.status, returncode=result.returncode)), flush=True)
        if result.status != 'PROCESS_EXITED' or result.returncode != 0:
            raise RuntimeError('TABLE_JOB_FAILED')
    finally:
        checked_signal(signal.SIGCONT)
        print('HEART_COORDINATOR_RESUMED', flush=True)


if __name__ == '__main__':
    main()
