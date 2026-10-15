from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
inputs=HERE/'noodles-formal-cherchi-derivatives-job-01/output'
files=dict(launcher=Path(__file__),worker=HERE/'certify_noodles_cherchi_samples.py',
    classifier=HERE/'saved_connection_certificate.py',replay=HERE/'cherchi_core_map.py',
    planar=HERE/'exact_planar_patches.py',
    request=HERE/'noodles-balanced-request-freeze-job-01/output/request.json')
import json
report=json.loads((inputs/'report.json').read_text())
if report.get('numerical_passed') is not True or report.get('point_count')!=81:
    raise ValueError('CHERCHI_DERIVATIVE_GATE_NOT_PASSED')
for name in json.loads(files['request'].read_text())['input_sources']:
    files[name]=Path(name)
for path in inputs.glob('*.json'):files[str(path)]=path
for path in (ROOT/'python/make_it_stand').rglob('*.py'):files[str(path)]=path
result=run_source_job(HERE/'noodles-formal-cherchi-connection-job-01',
    ['/home/lin/anaconda3/bin/python3','-B',str(files['worker']),'{job}/output'],files,
    environment={'PYTHONPATH':str(ROOT/'python')+':'+str(HERE)})
print(result.status,result.returncode,flush=True)
