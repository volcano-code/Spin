from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
files=dict(launcher=Path(__file__),worker=HERE/'classify_heart_failed_fd_samples.py',
    classifier=HERE/'saved_connection_certificate.py',replay=HERE/'cherchi_core_map.py',
    planar=HERE/'exact_planar_patches.py',registration=HERE/'heart-cherchi-coordinate-map-job-01/output/registration.json')
for directory in ('heart-cherchi-heart-first-small-face-rejection-derivative-job-01/output',
                  'heart-projective-independent-reference-job-01/output'):
    for path in (HERE/directory).glob('*.json'):files[str(path)]=path
for path in (ROOT/'python/make_it_stand').rglob('*.py'):files[str(path)]=path
result=run_source_job(HERE/'heart-small-face-finite-classification-job-01',
    ['/home/lin/anaconda3/bin/python3','-B',str(files['worker']),'{job}/output'],files,
    environment={'PYTHONPATH':str(ROOT/'python')+':'+str(HERE)})
print(result.status,result.returncode,flush=True)
