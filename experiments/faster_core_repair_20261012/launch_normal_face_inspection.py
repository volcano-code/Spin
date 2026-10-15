from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
worker=HERE/'build_and_inspect_normal_faces.py'
result=run_source_job(HERE/'normal-face-inspection-job-01',
    ['/home/lin/anaconda3/bin/python3','-B',str(worker),'{job}/output'],
    dict(launcher=Path(__file__),worker=worker,source=HERE/'inspect_projective_normal_faces.cpp',
         sphere=ROOT/'src/spin/detail/FixedNormalSphere.h',states=HERE/'frozen-heart-states.json',
         library=ROOT/'build/spin2014-objectives-quadric-20261007/libmake_it_stand_core.a',
         template=ROOT/'build/heart-fixed-facet-20261011-v1/commands.json'))
print(result.status,result.returncode,flush=True)
