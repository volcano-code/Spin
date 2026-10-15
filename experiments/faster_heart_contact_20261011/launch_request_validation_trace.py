"""One bounded native exception trace of the registered failed base point."""
import json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
HERE=Path(__file__).resolve().parent
registration=HERE/'rounded-clearance-v3-faster-production-v2/registration.json'
data=json.loads(registration.read_text())
native=ROOT/'build/spin2014-objectives-quadric-20261007/python/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
files={str(i):Path(p) for i,p in enumerate(data['sources'])}
files.update(registration=registration,controller=Path(__file__),debugger=Path('/usr/bin/gdb'))
result=run_source_job(HERE/'rounded-clearance-v3-validation-trace-01',
    ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6','/usr/bin/gdb','-q','-batch','-ex','catch throw','-ex','run','-ex','bt 16','-ex','quit',
     '--args','/home/lin/anaconda3/bin/python3','-B',str(HERE/'rounded_faster_production_point.py'),
     'point','--registration',str(registration),'--index','0','--output','{job}/output'],files,
    environment={'PYTHONPATH':str(ROOT/'python')+':'+str(native.parent)})
print(result.status,flush=True)
