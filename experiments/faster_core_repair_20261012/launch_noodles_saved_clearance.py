from pathlib import Path
import sys
import json
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
native=HERE/'formal-native-core-build-job-08/output/build/python/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
request=HERE/'noodles-balanced-request-freeze-job-01/output/request.json'
files=dict(launcher=Path(__file__),worker=HERE/'probe_noodles_saved_clearance.py',native=native,request=request)
for name in json.loads(request.read_text())['input_sources']:files[name]=Path(name)
for path in (ROOT/'python/make_it_stand').rglob('*.py'):files[str(path)]=path
result=run_source_job(HERE/'noodles-saved-input-clearance-job-01',
    ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
     '/home/lin/anaconda3/bin/python3','-B',str(files['worker']),'{job}/output'],files,
    environment={'PYTHONPATH':str(ROOT/'python')+':'+str(native.parent)})
print(result.status,result.returncode,flush=True)
