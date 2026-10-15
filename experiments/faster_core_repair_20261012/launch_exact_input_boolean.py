from pathlib import Path
import json
import sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
files=dict(launcher=Path(__file__),worker=HERE/'build_exact_input_boolean.py',
    source=HERE/'exact_input_boolean.cpp',test=HERE/'test_exact_union_recipe.py',
    template=ROOT/'build/heart-fixed-facet-20261011-v1/commands.json',
    receipt=HERE/'exact-core-material-build-job-01/output/receipt.json',
    library=ROOT/'build/spin2014-objectives-quadric-20261007/libmake_it_stand_core.a')
for i,path in enumerate(dict.fromkeys(x for x in json.loads(files['receipt'].read_text())['commands'][-1]['argv'] if 'libboost_' in x and x.endswith('.a'))):
    files[f'boost_{i}']=Path(path)
result=run_source_job(HERE/'exact-input-union-build-job-01',
    ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
     '/home/lin/anaconda3/bin/python3','-B',str(files['worker']),'{job}/output'],files)
print(result.status,result.returncode,flush=True)
