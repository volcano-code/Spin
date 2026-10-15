from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
result=run_source_job(HERE/'cgal-proof-capture-build-job-01',
 ['/home/lin/anaconda3/bin/python3','-B',str(HERE/'build_cgal_proof_capture.py'),'{job}/output'],
 dict(launcher=Path(__file__),worker=HERE/'build_cgal_proof_capture.py',binding=ROOT/'bindings/pybind_module.cpp',prior=HERE/'projective-interval-build-job-02/output/receipt.json',template=ROOT/'build/heart-fixed-facet-20261011-v1/commands.json',regression=HERE/'test_cgal_proof_capture.py'))
print(result.status,result.returncode,flush=True)
