"""Queue one fixed Boolean diagnosis after the existing Heart launcher exits."""
import json
import os
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def main():
    sys.path.insert(0, str(ROOT/'python'))
    from make_it_stand.verification.source_job import run_source_job
    output = HERE/'noodles-dense-v5-core-failed-boolean-capture-v1'
    output.mkdir(exist_ok=False)
    # Identity, not merely PID: never wait on an unrelated reused process.
    process = Path('/proc/872658')
    initial = (process/'stat').read_text().split(') ', 1)[1].split()[19] if process.exists() else None
    with (output/'events.jsonl').open('x') as stream:
        def event(status, **values):
            stream.write(json.dumps(dict(status=status, **values))+'\n')
            stream.flush(); os.fsync(stream.fileno())
        event('WAIT_EXISTING_HEART_CORE_V3')
        started = time.monotonic()
        while process.exists() and initial is not None:
            try:
                if (process/'stat').read_text().split(') ', 1)[1].split()[19] != initial:
                    break
            except FileNotFoundError:
                break
            if time.monotonic()-started > 4000:
                event('QUEUE_WAIT_EXPIRED_NO_JOB_STARTED'); return
            time.sleep(5)
        case = HERE/'noodles-dense-v5'
        manifest = json.loads((case/'manifest.json').read_text())
        part = case/manifest['parts'][manifest['construction']['fixed_parts'][72]]['model_obj']
        predecessor = HERE/'noodles-dense-v5-core-preprocessing-v1/domain-job-01/output/domain-072.json'
        native = ROOT/'build/heart-fixed-facet-20261011-v1/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
        worker = HERE/'capture_noodles_failed_boolean.py'
        event('FIXED_BOOLEAN_CAPTURE_STARTED')
        result = run_source_job(output/'capture-job-01',
            ['/usr/bin/env', 'LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
             '/home/lin/anaconda3/bin/python3', '-B', str(worker), '--case', str(case),
             '--predecessor', str(predecessor), '--output', '{job}/output'],
            dict(launcher=Path(__file__), worker=worker, manifest=case/'manifest.json',
                 predecessor=predecessor, fixed_part=part, native=native),
            environment={'PYTHONPATH':str(ROOT/'python')+':'+str(native.parent)})
        event('FIXED_BOOLEAN_CAPTURE_TERMINAL', status_code=result.status, returncode=result.returncode)


if __name__ == '__main__':
    main()
