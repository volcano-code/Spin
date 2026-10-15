from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
worker=HERE/'build_timing_candidate.py'
result=run_source_job(HERE/'boundary-proof-reuse-build-job-01',
    ['/home/lin/anaconda3/bin/python3','-B',str(worker),'{job}/output',str(HERE/'fixed-input-cache-build-job-01/output')],
    dict(worker=worker,launcher=Path(__file__),commands=ROOT/'build/heart-fixed-facet-20261011-v1/commands.json',
         evaluator=ROOT/'src/spin/QuadricMomentEvaluator.cpp',
         differential=ROOT/'src/spin/PolyhedralMomentDifferential.cpp',
         fullrank=ROOT/'src/spin/FullRankMomentEvaluator.cpp',
         canonical_commands=ROOT/'build/heart-exact-input-20261011-v1/commands.json',
         request_source=ROOT/'src/spin/SpinDesignRequest.cpp',
         reuse_receipt=HERE/'fixed-input-cache-build-job-01/output/receipt.json',
         header=ROOT/'include/spin/PolyhedralMomentDifferential.h'))
print(result.status,result.returncode,flush=True)
