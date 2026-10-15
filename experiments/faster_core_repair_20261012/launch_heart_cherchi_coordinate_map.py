from pathlib import Path
import json
import sys
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
native=HERE/'formal-native-core-build-job-08/output/build/python/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
bridge=HERE/'formal-cherchi-build-job-02/output/build/python/_cherchi_core_boolean.cpython-313-x86_64-linux-gnu.so'
arithmetic=HERE/'formal-rational-arithmetic-build-job-02/output/build/python/_cherchi_rational_moments.cpython-313-x86_64-linux-gnu.so'
heart=ROOT/'experiments/faster_heart_contact_20261011'
files=dict(launcher=Path(__file__),worker=HERE/'probe_heart_cherchi_coordinate_map.py',
    contract=HERE/'case_coordinate_contract.py',states=HERE/'frozen-heart-states.json',
    settings=heart/'core_faster_heart.py',request=heart/'rounded-clearance-v3-faster-frozen-v2/request.json',
    native=native,bridge=bridge,arithmetic=arithmetic)
for filename in json.loads(files['request'].read_text())['input_sources']:
    files[filename]=Path(filename)
for path in (ROOT/'python/make_it_stand').rglob('*.py'):files[str(path)]=path
result=run_source_job(HERE/'heart-cherchi-coordinate-map-job-01',
    ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
     '/home/lin/anaconda3/bin/python3','-B',str(files['worker']),'{job}/output'],files,
    environment={'PYTHONPATH':':'.join(map(str,[ROOT/'python',native.parent,bridge.parent,arithmetic.parent,HERE]))})
print(result.status,result.returncode,flush=True)
