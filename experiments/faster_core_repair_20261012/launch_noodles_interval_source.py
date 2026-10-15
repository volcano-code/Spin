from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
native=HERE/'projective-interval-build-job-02/output/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
files=dict(launcher=Path(__file__),worker=HERE/'probe_noodles_interval_source.py',native=native,
    weights=HERE/'probe_inprocess_derivatives.py',serializer=HERE/'run_small_inprocess_solver.py',
    map_config=HERE/'noodles-bound-seed0-production-probe-job-01/output/run/mapping-config.json',
    actual_events=HERE/'noodles-bound-seed0-production-probe-job-01/output/run/mapping-events.jsonl',
    triangulator=ROOT/'src/spin/QuadricTriangulator.cpp',motion=ROOT/'src/spin/QuadricSourceMotion.cpp')
result=run_source_job(HERE/'noodles-interval-source-derivative-job-01',
    ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
     '/home/lin/anaconda3/bin/python3','-B',str(files['worker']),'{job}/output'],files,
    environment={'PYTHONPATH':str(ROOT/'python')+':'+str(native.parent)+':'+str(HERE)})
print(result.status,result.returncode,flush=True)
