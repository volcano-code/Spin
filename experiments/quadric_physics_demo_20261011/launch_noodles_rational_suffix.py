import json
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
case=HERE/'noodles-dense-v5';manifest=json.loads((case/'manifest.json').read_text())
native=ROOT/'build/heart-fixed-facet-20261011-v1/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
worker=HERE/'noodles_rational_core_suffix.py'
capture=HERE/'noodles-dense-v5-core-failed-boolean-capture-v1/capture-job-02/output'
previous=HERE/'noodles-dense-v5-core-preprocessing-v1/domain-job-01/output/domain-072.json'
files=dict(launcher=Path(__file__),worker=worker,manifest=case/'manifest.json',native=native,
    canonicalizer=HERE/'core_boolean_mesh.py',predecessor=previous,
    capture=capture/'original-rational-output.json',diagnosis=capture/'diagnosis.json')
files.update({name:case/entry['model_obj'] for name,entry in manifest['parts'].items()})
sequence=HERE/'noodles-dense-v5-core-rational-suffix-v1';sequence.mkdir(exist_ok=False)
result=run_source_job(sequence/'domain-job-01',
    ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
     '/home/lin/anaconda3/bin/python3','-B',str(worker),'--case',str(case),
     '--capture',str(capture),'--predecessor',str(previous),'--output','{job}/output'],
    files,environment={'PYTHONPATH':str(ROOT/'python')+':'+str(native.parent)})
print(result.status,result.returncode,flush=True)
