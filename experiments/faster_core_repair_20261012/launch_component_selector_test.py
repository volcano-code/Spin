from pathlib import Path
import sys,json
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
preceding=HERE/'component-selector-derivatives-job-01/result.json'
if not preceding.exists():raise ValueError('PRECEDING_HEAVY_JOB_NOT_TERMINAL')
receipt=HERE/'component-selector-build-job-01/output/receipt.json'
files=dict(launcher=Path(__file__),worker=HERE/'build_component_selector_test.py',receipt=receipt,
    test=ROOT/'tests/quadric_component_selector_regression.cpp',
    template=ROOT/'build/heart-fixed-facet-20261011-v1/commands.json')
for i,path in enumerate(json.loads(receipt.read_text())['commands'][-1]['argv']):
    if path.endswith(('.o','.a')) and Path(path).is_file():files['link_'+str(i)]=Path(path)
result=run_source_job(HERE/'component-selector-contract-job-01',
    ['/home/lin/anaconda3/bin/python3','-B',str(files['worker']),'{job}/output'],files)
print(result.status,result.returncode,flush=True)
