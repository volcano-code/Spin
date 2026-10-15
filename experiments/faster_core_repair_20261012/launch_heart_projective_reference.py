from pathlib import Path
import sys
import json
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
capture=HERE/'heart-projective-reference-source-job-01/output'
if json.loads((capture/'report.json').read_text()).get('passed') is not True:
    raise ValueError('SOURCE_REPLAY_NOT_VERIFIED')
files=dict(launcher=Path(__file__),worker=HERE/'verify_heart_projective_reference.py',
    projective=HERE/'independent_projective_reference.py',reference=HERE/'independent_directional_reference.py',
    source=capture/'source.json',replay=capture/'report.json',
    material=HERE/'heart-cherchi-heart-first-small-face-rejection-derivative-job-01/output/point-000-material.json',
    registration=HERE/'heart-cherchi-coordinate-map-job-01/output/registration.json',
    normals=HERE/'core-registered-failure-capture-job-01/output/source_side.json')
result=run_source_job(HERE/'heart-projective-independent-reference-job-01',
    ['/home/lin/anaconda3/bin/python3','-B',str(files['worker']),'{job}/output'],files,
    environment={'PYTHONPATH':str(HERE)})
print(result.status,result.returncode,flush=True)
