"""Resource-bounded Cherchi FD check; invoke only in a free serial slot."""
from pathlib import Path
import json
import sys
HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
native=HERE/'formal-native-core-build-job-08/output/build/python/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
bridge=HERE/'formal-cherchi-build-job-02/output/build/python/_cherchi_core_boolean.cpython-313-x86_64-linux-gnu.so'
arithmetic=HERE/'formal-rational-arithmetic-build-job-02/output/build/python/_cherchi_rational_moments.cpython-313-x86_64-linux-gnu.so'
fixed=HERE/'noodles-formal-cherchi-fixed-job-01/output/report.json'
record=json.loads(fixed.read_text())
if not record['independent_passed'] or not record['numerical_backend_comparison_passed']:
    raise ValueError('FIXED_STATE_PRECHECK_FAILED')
files=dict(launcher=Path(__file__),worker=HERE/'probe_noodles_formal_cherchi_derivatives.py',
    native=native,bridge=bridge,arithmetic=arithmetic,fixed_state=fixed,
    integration=HERE/'independent_tetra_moments.py',weights=HERE/'probe_formal_cherchi_derivatives.py',
    request=HERE/'noodles-balanced-request-freeze-job-01/output/request.json',
    samples=HERE/'noodles-algebraic-production-moments-job-01/output/manifest.json')
for name in json.loads(files['request'].read_text())['input_sources']:files[name]=Path(name)
for path in (ROOT/'python/make_it_stand').rglob('*.py'): files[str(path)]=path
result=run_source_job(HERE/'noodles-formal-cherchi-derivatives-job-01',
    ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
     '/home/lin/anaconda3/bin/python3','-B',str(files['worker']),'{job}/output'],files,
    environment={'PYTHONPATH':':'.join(map(str,[ROOT/'python',native.parent,bridge.parent,arithmetic.parent,HERE]))})
print(result.status,result.returncode,flush=True)
