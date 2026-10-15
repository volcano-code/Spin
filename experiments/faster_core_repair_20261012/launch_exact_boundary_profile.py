from pathlib import Path
import json
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job

previous=HERE/'formal-core-benchmark-faces-12288-seed-0-job-01'
if json.loads((previous/'result.json').read_text())['status']!='PROCESS_EXITED':
    raise ValueError('PRECEDING_FORMAL_BENCHMARK_NOT_SUCCESSFUL_TERMINAL')
run=previous/'output/run'
terminal=json.loads((run/'terminal-material-reference.json').read_text())['evaluation_file']
native=HERE/'formal-native-core-build-job-04/output/build/python/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
files=dict(launcher=Path(__file__),worker=HERE/'profile_exact_boundary.py',driver=HERE/'profile_exact_boundary.cpp',
    native=native,commands=HERE/'formal-native-core-build-job-04/output/build/compile_commands.json',
    material=run/terminal,terminal=run/'terminal-material-reference.json',
    source=ROOT/'src/verification/ExactBoundaryAudit.cpp',header=ROOT/'include/verification/ExactBoundaryAudit.h',
    edge=ROOT/'include/verification/ExactSharedEdgeIntersection.h',sha=ROOT/'include/verification/Sha256.h',
    spin=ROOT/'python/make_it_stand/spin.py')
generation=sys.argv[1] if len(sys.argv)>1 else '01'
if generation not in ('01','02','03'):raise ValueError('UNREGISTERED_NATIVE_PROFILE_GENERATION')
result=run_source_job(HERE/('exact-boundary-profile-job-'+generation),
    ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
     '/home/lin/anaconda3/bin/python3','-B',str(files['worker']),'{job}/output',
     'callgraph' if generation=='01' else 'phases'],files,
    environment={'PYTHONPATH':str(ROOT/'python')+':'+str(native.parent)})
print(result.status,result.returncode,flush=True)
