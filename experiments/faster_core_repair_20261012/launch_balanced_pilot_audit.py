"""Audit the exact eight-part pilot; never starts a complete model solve."""
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
files=dict(launcher=Path(__file__),worker=HERE/'audit_repaired_noodles_rational.py',
    input=HERE/'noodles-balanced-eight-part-pilot-job-01/output/pilot-design.json',
    pilot_report=HERE/'noodles-balanced-eight-part-pilot-job-01/output/report.json',
    pilot_receipt=HERE/'noodles-balanced-eight-part-pilot-job-01/result.json',
    executable=HERE/'noodles-rational-authority-diagnosis-job-01/output/inspect-rational-boundary')
result=run_source_job(HERE/'noodles-balanced-eight-part-audit-job-01',
    ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
     '/home/lin/anaconda3/bin/python3','-B',str(files['worker']),'{job}/output',str(files['input'])],files)
print(result.status,result.returncode,flush=True)
