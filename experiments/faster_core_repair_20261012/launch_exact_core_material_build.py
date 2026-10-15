from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
worker=HERE/'build_exact_core_material_candidate.py'
result=run_source_job(HERE/'exact-core-material-build-job-01',
    ['/home/lin/anaconda3/bin/python3','-B',str(worker),'{job}/output'],
    dict(worker=worker,launcher=Path(__file__),source=ROOT/'src/spin/QuadricMomentEvaluator.cpp',
         template=ROOT/'build/heart-fixed-facet-20261011-v1/commands.json',
         receipt=HERE/'projective-belt-build-job-04/output/receipt.json'))
print(result.status,result.returncode,flush=True)
