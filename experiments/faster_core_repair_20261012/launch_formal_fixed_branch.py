from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
files=dict(launcher=Path(__file__),worker=HERE/'certify_formal_fixed_samples.py',
    fixture=HERE/'same-input-job-03/output/frozen-inputs.json',spin=ROOT/'python/make_it_stand/spin.py')
for name in ('cherchi_inprocess','cherchi_motion_differential'):
    files[name]=ROOT/f'python/make_it_stand/verification/{name}.py'
for name in ('cherchi_core_map.py','exact_planar_patches.py'):
    files[name]=HERE/name
for state in ('ellipsoid','one_sheet','two_sheet'):
    folder=HERE/'formal-fixed-derivatives-job-01/output'/state
    if not (folder/'report.json').exists():raise ValueError('DERIVATIVE_JOB_NOT_FINISHED')
    for path in folder.glob('*.json'):files[state+'_'+path.name]=path
result=run_source_job(HERE/'formal-fixed-branch-job-01',
    ['/home/lin/anaconda3/bin/python3','-B',str(files['worker']),'{job}/output'],files,
    environment={'PYTHONPATH':str(ROOT/'python')})
print(result.status,result.returncode,flush=True)
