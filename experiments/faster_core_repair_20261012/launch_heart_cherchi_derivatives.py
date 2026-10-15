from pathlib import Path
import json
import sys
HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
state=sys.argv[1]
if state not in ('heart-first-cap-failure','heart-first-small-face-rejection'):raise ValueError('UNREGISTERED_STATE')
native=HERE/'formal-native-core-build-job-08/output/build/python/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
bridge=HERE/'formal-cherchi-build-job-02/output/build/python/_cherchi_core_boolean.cpython-313-x86_64-linux-gnu.so'
arithmetic=HERE/'formal-rational-arithmetic-build-job-02/output/build/python/_cherchi_rational_moments.cpython-313-x86_64-linux-gnu.so'
prior=HERE/'heart-cherchi-coordinate-map-job-01/output'
files=dict(launcher=Path(__file__),worker=HERE/'probe_heart_cherchi_derivatives.py',native=native,bridge=bridge,
    arithmetic=arithmetic,integration=HERE/'independent_tetra_moments.py',weights=HERE/'probe_formal_cherchi_derivatives.py',
    registration=prior/'registration.json',admission=prior/'report.json',base=prior/f'{state}-normalized.json',
    frozen=ROOT/'experiments/faster_heart_contact_20261011/rounded-clearance-v3-faster-frozen-v2/request.json')
for name in json.loads(files['frozen'].read_text())['input_sources']:files[name]=Path(name)
for p in (ROOT/'python/make_it_stand').rglob('*.py'):files[str(p)]=p
result=run_source_job(HERE/f'heart-cherchi-{state}-derivative-job-01',
    ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6','/home/lin/anaconda3/bin/python3',
     '-B',str(files['worker']),'{job}/output',state],files,
    environment={'PYTHONPATH':':'.join(map(str,[ROOT/'python',native.parent,bridge.parent,arithmetic.parent,HERE]))})
print(result.status,result.returncode,flush=True)
