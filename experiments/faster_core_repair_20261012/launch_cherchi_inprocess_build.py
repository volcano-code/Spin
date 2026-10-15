from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
worker=HERE/'build_cherchi_inprocess.py'
result=run_source_job(HERE/'cherchi-inprocess-build-job-01',
    ['/home/lin/anaconda3/bin/python3','-B',str(worker),'{job}/output'],
    dict(worker=worker,launcher=Path(__file__),source=HERE/'cherchi_inprocess.cpp',
         probe_source=ROOT/'experiments/cherchi_reference_probe/probe.cpp',
         flags=ROOT/'build/cherchi-reference-singlethread-20260917/CMakeFiles/cherchi_probe.dir/flags.make',
         link=ROOT/'build/cherchi-reference-singlethread-20260917/CMakeFiles/cherchi_probe.dir/link.txt'))
print(result.status,result.returncode,flush=True)
