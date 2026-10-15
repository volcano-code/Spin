from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
worker=HERE/'build_projective_belt_candidate.py'
prior=HERE/'exact-core-material-build-job-01/output/receipt.json'
suffix=sys.argv[1] if len(sys.argv)==2 else '01'
if suffix not in ('01','02'):raise ValueError('UNREGISTERED_INTERVAL_BUILD_JOB')
result=run_source_job(HERE/('projective-interval-build-job-'+suffix),
    ['/home/lin/anaconda3/bin/python3','-B',str(worker),'{job}/output',str(prior)],
    dict(launcher=Path(__file__),worker=worker,prior_receipt=prior,
         source=ROOT/'src/spin/QuadricTriangulator.cpp',motion=ROOT/'src/spin/QuadricSourceMotion.cpp',
         header=ROOT/'include/spin/QuadricTriangulator.h',
         template=ROOT/'build/heart-fixed-facet-20261011-v1/commands.json',
         regression=HERE/'test_noodles_seed_source.py'))
print(result.status,result.returncode,flush=True)
