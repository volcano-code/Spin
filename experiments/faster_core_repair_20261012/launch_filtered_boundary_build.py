from pathlib import Path
import json,sys

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job

prior=HERE/'rational-boundary-binding-build-job-01/output/receipt.json'
files=dict(launcher=Path(__file__),worker=HERE/'build_filtered_boundary.py',prior=prior,
    source=ROOT/'src/verification/ExactBoundaryAudit.cpp',
    audit_header=ROOT/'include/verification/ExactBoundaryAudit.h',
    edge_header=ROOT/'include/verification/ExactSharedEdgeIntersection.h',
    template=ROOT/'build/heart-fixed-facet-20261011-v1/commands.json',
    test=ROOT/'tests/python/test_exact_filtered_boundary.py')
for i,path in enumerate(json.loads(prior.read_text())['commands'][-1]):
    if path.endswith(('.o','.a')) and Path(path).is_file():files['object_'+str(i)]=Path(path)
result=run_source_job(HERE/'filtered-boundary-build-job-01',
    ['/home/lin/anaconda3/bin/python3','-B',str(files['worker']),'{job}/output'],files)
print(result.status,result.returncode,flush=True)
