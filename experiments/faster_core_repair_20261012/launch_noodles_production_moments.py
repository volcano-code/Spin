from pathlib import Path
import json,sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
native=HERE/'projective-interval-build-job-02/output/_make_it_stand.cpython-313-x86_64-linux-gnu.so';request=HERE/'noodles-balanced-request-freeze-job-01/output/request.json';source=HERE/'noodles-interval-algebraic-source-derivative-job-01/output'
files=dict(launcher=Path(__file__),worker=HERE/'probe_noodles_production_moments.py',native=native,spin=ROOT/'python/make_it_stand/spin.py',integrator=HERE/'independent_tetra_moments.py',weights=HERE/'probe_inprocess_derivatives.py',serializer=HERE/'run_small_inprocess_solver.py',request=request,fixed=request.parent/'fixed.binary64.json',source_manifest=source/'manifest.json',source_gate=source/'report.json',build_receipt=native.parent/'receipt.json')
for i in range(81):files[f'point{i}']=source/f'point-{i:03d}.json'
for i,path in enumerate(json.loads(request.read_text())['input_sources']):files[f'input{i}']=Path(path)
result=run_source_job(HERE/'noodles-algebraic-production-moments-job-01',
 ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6','/home/lin/anaconda3/bin/python3','-B',str(files['worker']),'{job}/output'],files,
 environment={'PYTHONPATH':str(ROOT/'python')+':'+str(native.parent)+':'+str(HERE)})
print(result.status,result.returncode,flush=True)
