from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / 'python'))
from make_it_stand.verification.source_job import run_source_job
files = dict(launcher=Path(__file__), worker=HERE/'classify_noodles_proof_pair.py',
    recipes=HERE/'birth_face_recipes.py', chain=HERE/'exact_planar_patches.py',
    request=HERE/'noodles-balanced-request-freeze-job-01/output/request.json',
    capture_receipt=HERE/'noodles-same-call-proof-capture-job-01/output/report.json')
for i in (0, 1):
    files[f'capture{i}'] = HERE/f'noodles-same-call-proof-capture-job-01/output/point-{i:03d}.json'
    files[f'source{i}'] = HERE/f'noodles-interval-algebraic-source-derivative-job-01/output/point-{i:03d}.json'
    files[f'material{i}'] = HERE/f'noodles-algebraic-production-moments-job-01/output/point-{i:03d}-material.json'
result = run_source_job(HERE/'noodles-offline-proof-pair-job-01',
    ['/home/lin/anaconda3/bin/python3', '-B', str(files['worker']), '{job}/output'], files,
    environment={'PYTHONPATH':str(ROOT/'python')+':'+str(HERE)})
print(result.status, result.returncode, flush=True)
