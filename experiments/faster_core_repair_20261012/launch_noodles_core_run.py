"""One full registered seed per resource-bounded job; never changes seed/budget."""
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / 'python'))
from make_it_stand.verification.source_job import run_source_job
from noodles_core_run import PATHS, NATIVE, load_evidence, check_evidence

seed = int(sys.argv[1])
evidence = load_evidence(HERE)
check_evidence(evidence, seed)
files = {key: HERE / path for key, path in PATHS.items()}
files.update(launcher=Path(__file__), worker=HERE/'noodles_core_run.py', native=NATIVE,
             spin=ROOT/'python/make_it_stand/spin.py',
             journal=HERE.parent/'faster_heart_contact_20261011/rounded_solver_journal.py',
             build_receipt=NATIVE.parent/'receipt.json')
files.update({f'input_{i}': Path(path) for i, path in enumerate(evidence['frozen']['input_sources'])})
result = run_source_job(HERE/f'noodles-registered-core-seed-{seed}-job-01',
    ['/usr/bin/env', 'LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
     '/home/lin/anaconda3/bin/python3', '-B', str(files['worker']), '{job}/output', str(seed)], files,
    environment={'PYTHONPATH': str(ROOT/'python')+':'+str(NATIVE.parent)+':'+str(HERE)})
print(result.status, result.returncode, flush=True)
