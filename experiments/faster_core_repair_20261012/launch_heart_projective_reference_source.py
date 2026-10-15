from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
native=HERE/'formal-native-core-build-job-08/output/build/python/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
files=dict(launcher=Path(__file__),worker=HERE/'capture_heart_projective_reference_source.py',native=native,
    registration=HERE/'heart-cherchi-coordinate-map-job-01/output/registration.json',
    material=HERE/'heart-cherchi-heart-first-small-face-rejection-derivative-job-01/output/point-000-material.json')
for p in (ROOT/'python/make_it_stand').rglob('*.py'):files[str(p)]=p
result=run_source_job(HERE/'heart-projective-reference-source-job-01',
    ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
     '/home/lin/anaconda3/bin/python3','-B',str(files['worker']),'{job}/output'],files,
    environment={'PYTHONPATH':':'.join(map(str,[ROOT/'python',native.parent,HERE]))})
print(result.status,result.returncode,flush=True)
