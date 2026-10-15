"""Serial offline evidence check following the one active native performance probe."""
from pathlib import Path
import json
import os
import sys
import time
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
prior=HERE/'noodles-rational51-phase-diagnosis-job-02'
sequence=HERE/'saved-small-recipe-certification-sequence-01';sequence.mkdir(exist_ok=False)
with (sequence/'events.jsonl').open('x') as stream:
    def event(status,**values):
        stream.write(json.dumps(dict(status=status,**values))+'\n');stream.flush();os.fsync(stream.fileno())
    event('WAIT_ACTIVE_FIXED_BOOLEAN_PERFORMANCE_DIAGNOSIS');deadline=time.monotonic()+3700
    while not (prior/'complete.json').exists():
        if time.monotonic()>deadline:event('WAIT_LIMIT');sys.exit(1)
        time.sleep(5)
    inputs=HERE/'small-top-derivative-probe-job-02/output'
    worker=HERE/'certify_saved_small_recipe_batch.py'
    native=HERE/'boundary-proof-reuse-build-job-01/output/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
    files=dict(launcher=Path(__file__),worker=worker,native=native,
        manifest=inputs/'manifest.json',adapter=HERE/'cherchi_core_map.py',
        chain_checker=HERE/'exact_planar_patches.py',
        base_capture=HERE/'saved-small-provenance-capture-01/point-000.json',
        capture_receipt=HERE/'saved-small-provenance-capture-01/report.json',
        config=HERE/'small-top-inprocess-solver-sequence-01/seed-1-job-01/output/configuration.json')
    files.update({f'point_{i:03d}':inputs/f'point-{i:03d}-material.json' for i in range(81)})
    event('OFFLINE_SAVED_RECIPE_CHAIN_CHECK_STARTED_ZERO_MATERIAL_BOOLEAN')
    result=run_source_job(HERE/'saved-small-recipe-chain-certification-job-01',
        ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
         '/home/lin/anaconda3/bin/python3','-B',str(worker),'{job}/output'],files,
        environment={'PYTHONPATH':str(ROOT/'python')+':'+str(native.parent)})
    event('OFFLINE_RECIPE_CHAIN_CHECK_TERMINAL',status_code=result.status,returncode=result.returncode)
