from pathlib import Path
import json
import sys

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job

prior=HERE/'single-transport-build-job-01/output/receipt.json'
files=dict(launcher=Path(__file__),worker=HERE/'build_rational_boundary_binding.py',prior=prior,
    source=ROOT/'bindings/pybind_module.cpp',audit_header=ROOT/'include/verification/ExactBoundaryAudit.h',
    template=ROOT/'build/heart-fixed-facet-20261011-v1/commands.json',
    test=ROOT/'tests/python/test_rational_boundary_audit.py')
for i,path in enumerate(json.loads(prior.read_text())['commands'][-1]):
    if path.endswith(('.o','.a')) and Path(path).is_file():files['object_'+str(i)]=Path(path)
result=run_source_job(HERE/'rational-boundary-binding-build-job-01',
    ['/home/lin/anaconda3/bin/python3','-B',str(files['worker']),'{job}/output'],files)
print(result.status,result.returncode,flush=True)
