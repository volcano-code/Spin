from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
files=dict(launcher=Path(__file__),worker=HERE/'verify_saved_heart_coordinates.py',
    integration=HERE/'independent_tetra_moments.py',contract=HERE/'case_coordinate_contract.py')
for job in ('repaired-heart-production-moments-job-01','repaired-heart-exact-production-moments-job-01'):
    files[job]=HERE/job/'output/point-000-material.json'
result=run_source_job(HERE/'saved-heart-coordinate-equivalence-job-01',
    ['/home/lin/anaconda3/bin/python3','-B',str(files['worker']),'{job}/output'],files)
print(result.status,result.returncode,flush=True)
