from pathlib import Path
import json
import sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
case=ROOT/'experiments/quadric_physics_demo_20261011/noodles-dense-v5'
manifest=case/'manifest.json';data=json.loads(manifest.read_text())
native=HERE/'exact-core-material-build-job-01/output/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
worker=HERE/'diagnose_rational_boolean51_phases.py';name=data['construction']['fixed_parts'][50]
files=dict(launcher=Path(__file__),worker=worker,native=native,manifest=manifest,
    part=case/data['parts'][name]['model_obj'],probe=HERE/'profile_rational_boolean.cpp',
    prefix=HERE/'noodles-all-original-rational-job-02/output/domain-050-rational.json',
    events=HERE/'noodles-all-original-rational-job-02/output/events.jsonl',
    template=ROOT/'build/heart-fixed-facet-20261011-v1/commands.json',
    library=ROOT/'build/spin2014-objectives-quadric-20261007/libmake_it_stand_core.a',
    kernel=ROOT/'src/verification/ExactLibiglBoolean.cpp')
obj=HERE/'noodles-rational51-phase-diagnosis-job-01/output/profile.o'
files.update(precompiled_object=obj,origin=HERE/'noodles-rational51-phase-diagnosis-job-01/config.json',
    link_receipt=HERE/'exact-core-material-build-job-01/output/receipt.json')
for i,path in enumerate(dict.fromkeys(x for x in json.loads(files['link_receipt'].read_text())['commands'][-1]['argv'] if 'libboost_' in x and x.endswith('.a'))):
    files[f'boost_{i}']=Path(path)
result=run_source_job(HERE/'noodles-rational51-phase-diagnosis-job-02',
    ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
     '/home/lin/anaconda3/bin/python3','-B',str(worker),'{job}/output',str(obj)],files,
    environment={'PYTHONPATH':str(ROOT/'python')+':'+str(native.parent)})
print(result.status,result.returncode,flush=True)
