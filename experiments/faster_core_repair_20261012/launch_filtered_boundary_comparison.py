from pathlib import Path
import json
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job

files=dict(launcher=Path(__file__),worker=HERE/'compare_filtered_boundary.py')
for role,job in (('before','rational-boundary-binding-build-job-01'),('after','filtered-boundary-build-job-01')):
    files[role]=HERE/job/'output/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
    files[role+'_receipt']=HERE/job/'output/receipt.json'
for state in ('ellipsoid','one_sheet','two_sheet'):
    files[state]=HERE/'projective-formal-pipeline-job-01/output/fixed'/(state+'-material.json')
run=HERE/'core-benchmark-faces-12288-seed-0-job-01/output/run'
terminal=json.loads((run/'terminal-material-reference.json').read_text())['evaluation_file']
files['initial']=run/'evaluation-000001.json'
files['terminal']=run/terminal
files['terminal_reference']=run/'terminal-material-reference.json'
result=run_source_job(HERE/'filtered-boundary-comparison-job-01',
    ['/home/lin/anaconda3/bin/python3','-B',str(files['worker']),'{job}/output'],files)
print(result.status,result.returncode,flush=True)
