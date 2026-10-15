"""Insert one bounded native diagnosis between frozen Heart point jobs."""
from contextlib import contextmanager
import json
import os
from pathlib import Path
import signal
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def process_state(pid):
    for row in Path(f'/proc/{pid}/status').read_text().splitlines():
        if row.startswith('State:'):
            return row.split()[1]
    raise RuntimeError('MISSING_PROCESS_STATE')


def process_identity(pid, marker):
    proc = Path(f'/proc/{pid}')
    if marker.encode() not in (proc/'cmdline').read_bytes().split(b'\0'):
        # Python script paths are absolute; accept exact basename, not a substring.
        names = [Path(x.decode()).name for x in (proc/'cmdline').read_bytes().split(b'\0') if x]
        if marker not in names:
            raise RuntimeError('DISPATCHER_IDENTITY_CHANGED')
    return (proc/'stat').read_text().rsplit(')', 1)[1].split()[19]


@contextmanager
def paused_dispatcher(pid, marker):
    identity = process_identity(pid, marker)
    if process_state(pid) == 'T':
        raise RuntimeError('DISPATCHER_ALREADY_PAUSED_BY_ANOTHER_OWNER')
    os.kill(pid, signal.SIGSTOP)
    try:
        deadline = time.monotonic()+2
        while process_state(pid) != 'T':
            if time.monotonic() >= deadline:
                raise RuntimeError('DISPATCHER_STOP_NOT_CONFIRMED')
            time.sleep(.01)
        yield
    finally:
        if Path(f'/proc/{pid}').exists() and process_identity(pid, marker) == identity:
            os.kill(pid, signal.SIGCONT)


def main():
    sys.path.insert(0, str(ROOT/'python'))
    from make_it_stand.verification.source_job import run_source_job
    output = HERE/'three-dimensional-faster-performance-v1'
    output.mkdir(exist_ok=False)
    with (output/'events.jsonl').open('x') as stream:
        def event(status, **values):
            row = dict(status=status, **values)
            stream.write(json.dumps(row)+'\n'); stream.flush(); os.fsync(stream.fileno())
            print(json.dumps(row), flush=True)
        coordinator = 1766941
        with paused_dispatcher(coordinator, 'rounded_production_sequence.py'):
            queue = ROOT/'experiments/faster_heart_contact_20261011/rounded-clearance-v3-production-sequence-01/events.jsonl'
            rows = [json.loads(row) for row in queue.read_text().splitlines()]
            active = Path(next(row['job'] for row in reversed(rows) if row['status']=='POINT_STARTED'))
            event('WAIT_ACTIVE_HEART_POINT_NO_WORKER_INTERRUPTION', job=str(active), coordinator=coordinator)
            deadline = time.monotonic()+3660
            while not (active/'complete.json').exists():
                if time.monotonic() > deadline:
                    raise RuntimeError('HEART_POINT_NOT_TERMINAL')
                time.sleep(2)
            receipt = json.loads((active/'result.json').read_text())
            event('HEART_POINT_TERMINAL', result=receipt['status'], returncode=receipt['returncode'])
            native = ROOT/'build/spin2014-objectives-quadric-20261007/python/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
            boundary = HERE/'noodles-dense-v5-domain-stage-04/output/design.rational.json'
            profiler = HERE/'native_geometry_profile.py'
            worker = HERE/'profile_noodles_final_region.py'
            event('NATIVE_DIAGNOSTIC_STARTED', budget_seconds=3600,
                  geometry_admission=False, optimization_admission=False)
            result = run_source_job(output/'noodles-native-job-01',
                ['/usr/bin/env', 'LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
                 '/home/lin/anaconda3/bin/python3', '-B', str(profiler),
                 '--output', '{job}/profile', '--interval', '60', '--maximum-samples', '24', '--',
                 '/home/lin/anaconda3/bin/python3', '-B', str(worker),
                 '--source', str(boundary), '--output', '{job}/material-diagnostic'],
                {'launcher':Path(__file__), 'profiler':profiler, 'worker':worker, 'boundary':boundary,
                 'native':native, 'debugger':Path('/usr/bin/gdb'),
                 'snapshot':ROOT/'python/make_it_stand/verification/spin_problem.py'},
                environment={'PYTHONPATH':str(ROOT/'python')+':'+str(native.parent)})
            event('NATIVE_DIAGNOSTIC_TERMINAL', result=result.status, returncode=result.returncode)
        event('HEART_DISPATCHER_RESUMED', coordinator=coordinator)


if __name__ == '__main__':
    main()
