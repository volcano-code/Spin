"""Follow upstream terminal outcomes, then request authoritative serial admission."""
import importlib.util
import json
from pathlib import Path
import time

HERE = Path(__file__).resolve().parent


def wait_record(job):
    path = HERE/job/'result.json'
    while not path.is_file():
        time.sleep(10)
    return json.loads(path.read_text())


def main():
    capture = wait_record('inertia-boolean65-capture-job-01')
    if capture['status'] == 'PROCESS_EXITED' and capture['returncode'] == 0:
        wait_record('inertia-boolean65-stack-job-01')
    else:
        print('UPSTREAM_CAPTURE_FAILED_STACK_NOT_SCHEDULED', flush=True)
    spec = importlib.util.spec_from_file_location('heart_direction_launcher', HERE/'launch_heart_actual_direction.py')
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    while True:
        try:
            module.main()
            print('HEART_DIRECTION_JOB_RETURNED', flush=True)
            return
        except Exception as error:
            if any(code in str(error) for code in ('PROJECT_JOB_BUSY', 'PROJECT_WORKLOAD_ACTIVE')):
                time.sleep(10)
            else:
                raise


if __name__ == '__main__':
    main()
