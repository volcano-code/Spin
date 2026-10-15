from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
state=sys.argv[1]
generation=sys.argv[2] if len(sys.argv)>2 else '01'
if generation not in ('01','02'):raise ValueError('UNREGISTERED_GENERATION')
if state not in ('heart-first-cap-failure','heart-first-small-face-rejection'):raise ValueError('UNREGISTERED_STATE')
inputs=HERE/f'heart-cherchi-{state}-derivative-job-01/output'
files=dict(launcher=Path(__file__),worker=HERE/'certify_heart_cherchi_samples.py',
    classifier=HERE/'saved_connection_certificate.py',replay=HERE/'cherchi_core_map.py',
    planar=HERE/'exact_planar_patches.py',registration=HERE/'heart-cherchi-coordinate-map-job-01/output/registration.json')
for path in inputs.glob('*.json'):files[str(path)]=path
for path in (ROOT/'python/make_it_stand').rglob('*.py'):files[str(path)]=path
result=run_source_job(HERE/f'heart-cherchi-{state}-connection-job-{generation}',
    ['/home/lin/anaconda3/bin/python3','-B',str(files['worker']),'{job}/output',state,
     'classify' if generation=='02' else 'strict'],files,
    environment={'PYTHONPATH':str(ROOT/'python')+':'+str(HERE)})
print(result.status,result.returncode,flush=True)
