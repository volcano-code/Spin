from pathlib import Path
import json
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job

files=dict(launcher=Path(__file__),worker=HERE/'compare_filtered_boundary.py',
    before=HERE/'formal-native-core-build-job-04/output/build/python/_make_it_stand.cpython-313-x86_64-linux-gnu.so',
    before_receipt=HERE/'formal-native-core-build-job-04/output/receipt.json',
    after=HERE/'ranked-boundary-build-job-01/output/_make_it_stand.cpython-313-x86_64-linux-gnu.so',
    after_receipt=HERE/'ranked-boundary-build-job-01/output/receipt.json')
for state in ('ellipsoid','one_sheet','two_sheet'):
    files[state]=HERE/'projective-formal-pipeline-job-01/output/fixed'/(state+'-material.json')
run=HERE/'core-benchmark-faces-12288-seed-0-job-01/output/run'
terminal=json.loads((run/'terminal-material-reference.json').read_text())['evaluation_file']
files['initial']=run/'evaluation-000001.json';files['terminal']=run/terminal
files['terminal_reference']=run/'terminal-material-reference.json'
result=run_source_job(HERE/'ranked-boundary-comparison-job-01',
    ['/home/lin/anaconda3/bin/python3','-B',str(files['worker']),'{job}/output','ranked'],files)
print(result.status,result.returncode,flush=True)
