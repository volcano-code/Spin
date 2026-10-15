from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
native=HERE/'projective-interval-build-job-02/output/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
files=dict(launcher=Path(__file__),worker=HERE/'probe_noodles_algebraic_material.py',native=native,spin=ROOT/'python/make_it_stand/spin.py',journal=HERE.parent/'faster_heart_contact_20261011/rounded_solver_journal.py',request=HERE/'noodles-balanced-request-freeze-job-01/output/request.json',source_gate=HERE/'noodles-interval-algebraic-source-derivative-job-01/output/report.json',source_manifest=HERE/'noodles-interval-algebraic-source-derivative-job-01/output/manifest.json',build_receipt=native.parent/'receipt.json')
result=run_source_job(HERE/'noodles-algebraic-first-material-job-01',
 ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6','/home/lin/anaconda3/bin/python3','-B',str(files['worker']),'{job}/output'],files,
 environment={'PYTHONPATH':str(ROOT/'python')+':'+str(native.parent)+':'+str(HERE)})
print(result.status,result.returncode,flush=True)
