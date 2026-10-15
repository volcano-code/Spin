from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / 'python'))
from make_it_stand.verification.source_job import run_source_job

files = dict(launcher=Path(__file__), worker=HERE / 'prepare_core_benchmark_inputs.py',
    refinement=HERE / 'same_entity_refinement.py', planar=HERE / 'exact_planar_patches.py',
    integrator=HERE / 'independent_tetra_moments.py',
    config=HERE / 'small-top-inprocess-solver-sequence-01/seed-1-job-01/output/configuration.json',
    fixture=HERE / 'same-input-job-03/output/frozen-inputs.json')
result = run_source_job(HERE / 'core-benchmark-inputs-job-01',
    ['/home/lin/anaconda3/bin/python3', '-B', str(files['worker']), '{job}/output'], files)
print(result.status, result.returncode, flush=True)
