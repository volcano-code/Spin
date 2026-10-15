from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job

base=HERE/'formal-native-core-build-job-04/output/build'
files=dict(launcher=Path(__file__),worker=HERE/'build_ranked_boundary_candidate.py',
    receipt=base.parent/'receipt.json',commands=base/'compile_commands.json',
    link=base/'CMakeFiles/_make_it_stand.dir/link.txt',source=ROOT/'src/verification/ExactBoundaryAudit.cpp',
    header=ROOT/'include/verification/ExactCoordinateRankBox.h',
    boundary_header=ROOT/'include/verification/ExactBoundaryAudit.h',
    edge=ROOT/'include/verification/ExactSharedEdgeIntersection.h',sha=ROOT/'include/verification/Sha256.h')
for i,path in enumerate(base.rglob('*')):
    if path.is_file() and path.suffix in ('.o','.a'):files['retained_'+str(i)]=path
result=run_source_job(HERE/'ranked-boundary-build-job-01',
    ['/home/lin/anaconda3/bin/python3','-B',str(files['worker']),'{job}/output'],files)
print(result.status,result.returncode,flush=True)
