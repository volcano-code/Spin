from pathlib import Path
import sys,json
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
native=HERE/'component-selector-build-job-01/output/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
bridge=HERE/'formal-cherchi-build-job-02/output/build/python/_cherchi_core_boolean.cpython-313-x86_64-linux-gnu.so'
registration=HERE/'noodles-first-core-failures-registration-01/registration.json'
files=dict(launcher=Path(__file__),worker=HERE/'probe_component_selector_derivatives.py',native=native,bridge=bridge,
    registration=registration,config=Path(json.loads(registration.read_text())['mapping_config']),
    spin=ROOT/'python/make_it_stand/spin.py',integrator=HERE/'independent_tetra_moments.py',
    probe=HERE/'probe_formal_fixed_derivatives.py',stencil=HERE/'probe_formal_cherchi_derivatives.py')
for name in ('cherchi_inprocess','cherchi_motion_differential','faster_material','spin_problem','mass_moments','rational_text'):
    files[name]=ROOT/f'python/make_it_stand/verification/{name}.py'
result=run_source_job(HERE/'component-selector-derivatives-job-01',
    ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6','/home/lin/anaconda3/bin/python3','-B',str(files['worker']),'{job}/output'],files,
    environment={'PYTHONPATH':str(ROOT/'python')+':'+str(native.parent)+':'+str(bridge.parent)})
print(result.status,result.returncode,flush=True)
