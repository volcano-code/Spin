from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
files=dict(launcher=Path(__file__),worker=HERE/'compare_single_transport.py',
    fixture=HERE/'same-input-job-03/output/frozen-inputs.json',spin=ROOT/'python/make_it_stand/spin.py')
for role,job in (('before','component-selector-build-job-01'),('after','single-transport-build-job-01')):
    files[role]=HERE/job/'output/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
result=run_source_job(HERE/'single-transport-comparison-job-01',
    ['/home/lin/anaconda3/bin/python3','-B',str(files['worker']),'{job}/output'],files)
print(result.status,result.returncode,flush=True)
