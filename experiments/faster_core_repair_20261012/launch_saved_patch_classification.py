from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
inputs=HERE/'small-top-derivative-probe-job-02/output'
worker=HERE/'classify_saved_planar_patches.py'
bound=dict(worker=worker,checker=HERE/'exact_planar_patches.py',manifest=inputs/'manifest.json')
bound.update({f'point_{i:03d}':inputs/f'point-{i:03d}-material.json' for i in range(81)})
result=run_source_job(HERE/'small-top-planar-patch-classification-job-01',
    ['/home/lin/anaconda3/bin/python3','-B',str(worker),'--inputs',str(inputs),'--output','{job}/output'],bound,
    environment={'PYTHONPATH':str(ROOT/'python')+':'+str(HERE)})
print(result.status,result.returncode,flush=True)
