"""One export after a proven witness defect repair; no new input construction."""
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
demo=ROOT/'experiments/quadric_physics_demo_20261011'
heart=ROOT/'experiments/faster_heart_contact_20261011'
prior=HERE/'noodles-balanced-input-moments-job-01/output'
files=dict(launcher=Path(__file__),worker=HERE/'export_balanced_shared_input.py',
    assembly=HERE/'shared_input_boundary.py',representation=demo/'core_boolean_mesh.py',
    witness_test=demo/'test_edge_erasure_witness.py',diagnosis=HERE/'diagnose_shared_erased_patch.py',
    integration=heart/'rounded_input.py',report=prior/'report.json',capacity=prior/'capacity.json',
    manifest=demo/'noodles-dense-v5/manifest.json',
    native=HERE/'exact-core-material-build-job-01/output/_make_it_stand.cpython-313-x86_64-linux-gnu.so')
for name in ('design','outer','fixed'):
    files[name]=HERE/f'noodles-balanced-full-input-job-01/output/{name}-rational.json'
    files[name+'-moments']=prior/(name+'-moments-exact.json')
result=run_source_job(HERE/'noodles-balanced-shared-export-job-02',
    ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
     '/home/lin/anaconda3/bin/python3','-B',str(files['worker']),'{job}/output','SURVIVING_EDGE_V2'],files,
    environment={'PYTHONPATH':str(ROOT/'python')+':'+str(files['native'].parent)})
print(result.status,result.returncode,flush=True)
