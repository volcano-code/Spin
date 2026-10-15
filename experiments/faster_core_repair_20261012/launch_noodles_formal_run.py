"""No fallback or answer initialization; one frozen seed per resource job."""
from pathlib import Path
import hashlib
import json
import sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
from inherited_job_identity import inherited_files
seed=int(sys.argv[1]);repeat=len(sys.argv)==3 and sys.argv[2]=='repeat'
if seed not in (0,1) or len(sys.argv) not in (2,3) or (len(sys.argv)==3 and not repeat):
    raise ValueError('UNREGISTERED_RUN_OR_REPEAT')
job=HERE/'noodles-formal-prepare-job-01'
result=json.loads((job/'result.json').read_text())
if result['status']!='PROCESS_EXITED' or result['returncode']!=0:
    raise ValueError('PREPARATION_RESOURCE_JOB_NOT_SUCCESSFUL')
if not json.loads((job/'output/report.json').read_text())['passed']:
    raise ValueError('PREPARATION_NOT_PASSED')
config=json.loads((job/'config.json').read_text())
files=inherited_files(config['identity'])
files.update(launcher=Path(__file__),worker=HERE/'run_noodles_formal_case.py',
    identity_helper=HERE/'inherited_job_identity.py',
    preparation_helper=HERE/'noodles_formal_case.py',preparation_config=job/'config.json',
    preparation_result=job/'result.json',preparation_report=job/'output/report.json',
    prepared_request=job/'output/prepared-request.json')
suffix='repeat-01' if repeat else 'job-01'
result=run_source_job(HERE/f'noodles-formal-cherchi-seed-{seed}-{suffix}',
    ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
     '/home/lin/anaconda3/bin/python3','-B',str(files['worker']),'{job}/output',str(seed)],files,
    environment={'PYTHONPATH':config['environment']['PYTHONPATH']})
print(result.status,result.returncode,flush=True)
