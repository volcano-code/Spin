from pathlib import Path
import sys
import json
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
native=HERE/'formal-native-core-build-job-08/output/build/python/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
bridge=HERE/'formal-cherchi-build-job-02/output/build/python/_cherchi_core_boolean.cpython-313-x86_64-linux-gnu.so'
arithmetic=HERE/'formal-rational-arithmetic-build-job-02/output/build/python/_cherchi_rational_moments.cpython-313-x86_64-linux-gnu.so'
files=dict(launcher=Path(__file__),worker=HERE/'noodles_formal_case.py',native=native,bridge=bridge,arithmetic=arithmetic,
    aggregator=HERE/'finite_sample_summary.py',request=HERE/'noodles-balanced-request-freeze-job-01/output/request.json',
    numeric=HERE/'noodles-formal-cherchi-derivatives-job-01/output/report.json',
    connection=HERE/'noodles-formal-cherchi-connection-job-01/output/report.json',
    clearance=HERE/'noodles-saved-input-clearance-job-01/output/report.json',
    base=HERE/'noodles-formal-cherchi-derivatives-job-01/output/point-000-material.json',
    validation=HERE/'layered-formal-validation-job-06/output/report.json')
if not json.loads(files['validation'].read_text())['passed']:raise ValueError('FORMAL_CORE_NOT_VALIDATED')
for name in json.loads(files['request'].read_text())['input_sources']:files[name]=Path(name)
for path in (ROOT/'python/make_it_stand').rglob('*.py'):files[str(path)]=path
result=run_source_job(HERE/'noodles-formal-prepare-job-01',
    ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
     '/home/lin/anaconda3/bin/python3','-B',str(files['worker']),'{job}/output','0'],files,
    environment={'PYTHONPATH':':'.join(map(str,[ROOT/'python',native.parent,bridge.parent,arithmetic.parent,HERE]))})
print(result.status,result.returncode,flush=True)
