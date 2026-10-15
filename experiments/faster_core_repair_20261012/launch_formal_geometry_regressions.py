from pathlib import Path
import json
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job

if not json.loads((HERE/'formal-native-saved-derivatives-job-01/output/report.json').read_text())['passed']:
    raise ValueError('FORMAL_SAVED_MAP_GATE_NOT_PASSED')
build=HERE/'formal-native-core-build-job-04'
native=build/'output/build/python/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
files=dict(launcher=Path(__file__),worker=HERE/'run_formal_geometry_regressions.py',native=native,
    receipt=build/'output/receipt.json',spin=ROOT/'python/make_it_stand/spin.py')
for name in ('test_saved_heart_projective_source','test_registered_source_side','check_repaired_heart_states',
    'verify_small_face_candidate','test_registered_material_normalization','phase_timing','run_small_inprocess_solver'):
    files[name]=HERE/(name+'.py')
files['heart_states']=HERE/'frozen-heart-states.json'
heart=ROOT/'experiments/faster_heart_contact_20261011'
files['heart_request']=heart/'rounded-clearance-v3-faster-frozen-v2/request.json'
files['heart_settings']=heart/'core_faster_heart.py'
packet=json.loads(files['heart_request'].read_text())
for i,path in enumerate(packet['input_sources']):files['heart_input_'+str(i)]=Path(path)
files['noodles_registration']=HERE/'noodles-first-core-failures-registration-01/registration.json'
registration=json.loads(files['noodles_registration'].read_text())
files['noodles_config']=Path(registration['mapping_config'])
files['source_capture']=HERE/'core-registered-failure-capture-job-01/output/source_side.json'
result=run_source_job(HERE/'formal-registered-geometry-check-job-01',
    ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
     '/home/lin/anaconda3/bin/python3','-B',str(files['worker']),'{job}/output'],files,
    environment={'PYTHONPATH':':'.join(map(str,[ROOT/'python',native.parent,HERE]))})
print(result.status,result.returncode,flush=True)
