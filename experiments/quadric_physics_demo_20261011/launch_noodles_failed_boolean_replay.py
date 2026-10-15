"""Replay only the frozen 73rd difference after the serialization repair."""
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / 'python'))
from make_it_stand.verification.source_job import run_source_job

case = HERE / 'noodles-dense-v5'
manifest = json.loads((case / 'manifest.json').read_text())
part = case / manifest['parts'][manifest['construction']['fixed_parts'][72]]['model_obj']
predecessor = HERE / 'noodles-dense-v5-core-preprocessing-v1/domain-job-01/output/domain-072.json'
native = ROOT / 'build/heart-fixed-facet-20261011-v1/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
worker = HERE / 'capture_noodles_failed_boolean.py'
result = run_source_job(
    HERE / 'noodles-dense-v5-core-failed-boolean-capture-v1/capture-job-02',
    ['/usr/bin/env', 'LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
     '/home/lin/anaconda3/bin/python3', '-B', str(worker), '--case', str(case),
     '--predecessor', str(predecessor), '--output', '{job}/output'],
    dict(launcher=Path(__file__), worker=worker, manifest=case / 'manifest.json',
         predecessor=predecessor, fixed_part=part, native=native),
    environment={'PYTHONPATH': str(ROOT / 'python') + ':' + str(native.parent)})
print(result.status, result.returncode, flush=True)
