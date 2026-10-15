from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
worker=HERE/'diagnose_noodles_rational_authority.py'
packet=ROOT/'experiments/quadric_physics_demo_20261011/noodles-dense-v5-core-rational-suffix-v1/domain-job-01/output/domain-081-rational.json'
result=run_source_job(HERE/'noodles-rational-authority-diagnosis-job-01',
    ['/home/lin/anaconda3/bin/python3','-B',str(worker),'{job}/output'],
    dict(worker=worker,launcher=Path(__file__),source=HERE/'inspect_rational_boundary.cpp',input=packet,
        template=ROOT/'build/heart-fixed-facet-20261011-v1/commands.json',
        library=ROOT/'build/spin2014-objectives-quadric-20261007/libmake_it_stand_core.a'))
print(result.status,result.returncode,flush=True)
