"""One-time migration inventory and explicit stopped-job reconciliation."""
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT/'python'))
from make_it_stand.verification.serial_source_job import (
    register_legacy_source_jobs, reconcile_terminated_source_job)


def main():
    directories = []
    for path in sorted(HERE.glob('*/config.json')):
        config = json.loads(path.read_text())
        if config.get('schema') == 'source-job-v1':
            directories.append(path.parent)
    register_legacy_source_jobs(ROOT, directories)
    recovered = []
    for directory in directories:
        if not (directory/'result.json').exists():
            record = reconcile_terminated_source_job(ROOT, directory)
            recovered.append(dict(directory=str(directory), unit=record['unit'],
                status=record['status'], scientific_validation=record['scientific_validation']))
    print(json.dumps(dict(registered=len(directories), termination_observations=recovered)), flush=True)


if __name__ == '__main__':
    main()
