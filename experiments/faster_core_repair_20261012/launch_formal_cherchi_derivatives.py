from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
native=HERE/'core-small-face-build-job-01/output/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
bridge=HERE/'formal-cherchi-build-job-02/output/build/python/_cherchi_core_boolean.cpython-313-x86_64-linux-gnu.so'
files=dict(launcher=Path(__file__),worker=HERE/'probe_formal_cherchi_derivatives.py',native=native,bridge=bridge,
    configuration=HERE/'small-top-inprocess-solver-sequence-01/seed-1-job-01/output/configuration.json',
    fixture=HERE/'same-input-job-03/output/frozen-inputs.json',spin=ROOT/'python/make_it_stand/spin.py')
for name in ('cherchi_inprocess','cherchi_motion_differential','faster_material','spin_problem',
             'spin_bundle','spin_artifact','spin_common','mass_moments','rational_text'):
    files[name]=ROOT/f'python/make_it_stand/verification/{name}.py'
files['independent_integrator']=HERE/'independent_tetra_moments.py'
files['base_configuration']=HERE/'formal-cherchi-small-closure-job-03/output/run/configuration.json'
files['base_result']=HERE/'formal-cherchi-small-closure-job-03/output/run/result.json'
result=run_source_job(HERE/'formal-cherchi-derivative-job-01',
    ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6','/home/lin/anaconda3/bin/python3','-B',str(files['worker']),'--run',str(files['base_result'].parent),'--output','{job}/output'],files,
    environment={'PYTHONPATH':str(ROOT/'python')+':'+str(native.parent)+':'+str(bridge.parent)})
print(result.status,result.returncode,flush=True)

