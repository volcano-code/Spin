from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
worker=HERE/'build_projective_belt_candidate.py'
result=run_source_job(HERE/'projective-belt-build-job-04',
    ['/home/lin/anaconda3/bin/python3','-B',str(worker),'{job}/output'],
    dict(launcher=Path(__file__),worker=worker,source=ROOT/'src/spin/QuadricTriangulator.cpp',
         motion=ROOT/'src/spin/QuadricSourceMotion.cpp',header=ROOT/'include/spin/QuadricTriangulator.h',
         template=ROOT/'build/heart-fixed-facet-20261011-v1/commands.json',
         prior_receipt=HERE/'newton-capture-build-job-01/output/receipt.json'))
print(result.status,result.returncode,flush=True)
