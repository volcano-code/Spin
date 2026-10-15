"""Isolated compile/link/load test of a fixed-input construction candidate."""
import hashlib
import json
from pathlib import Path
import os
import subprocess
import sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
output=Path(sys.argv[1]);output.mkdir(parents=True,exist_ok=False)
source=HERE/'exact_input_boolean.cpp';obj=output/'input.o';exe=output/'exact-input-boolean'
template=ROOT/'build/heart-fixed-facet-20261011-v1/commands.json'
command=list(json.loads(template.read_text())[0])
for flag,value in (('-c',str(source)),('-o',str(obj)),('-MF',str(obj)+'.d'),('-MT',str(obj))):
    command[command.index(flag)+1]=value
library=ROOT/'build/spin2014-objectives-quadric-20261007/libmake_it_stand_core.a'
receipt=HERE/'exact-core-material-build-job-01/output/receipt.json'
boost=list(dict.fromkeys(x for x in json.loads(receipt.read_text())['commands'][-1]['argv'] if 'libboost_' in x and x.endswith('.a')))
if not boost:raise ValueError('PINNED_BOOST_LINK_INPUTS_MISSING')
link=[command[0],str(obj),str(library),'-Wl,--start-group',*boost,'-Wl,--end-group',
    '-L/usr/lib/x86_64-linux-gnu','-lgmp','-lmpfr','-lpthread','-o',str(exe)]
with (output/'build.log').open('x') as stream:
    subprocess.run(command,cwd=ROOT,stdout=stream,stderr=stream,check=True,timeout=1800)
    subprocess.run(link,cwd=ROOT,stdout=stream,stderr=stream,check=True,timeout=300)
test=HERE/'test_exact_union_recipe.py'
with (output/'tests.log').open('x') as stream:
    subprocess.run([sys.executable,'-B',str(test)],cwd=ROOT,
        env=dict(os.environ,EXACT_BOOLEAN_PROBE=str(exe)),stdout=stream,stderr=stream,check=True,timeout=300)
sha=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
with (output/'receipt.json').open('x') as stream:
    json.dump(dict(commands=[command,link],input_sha256={str(p):sha(p) for p in
        [source,template,receipt,library,test,*map(Path,boost)]},executable_sha256=sha(exe),
        role='SMALL_EXACT_SET_SEMANTICS_ONLY',geometry_admission=False,optimization_admission=False),stream,indent=2)
print('Small exact set semantics passed; no model admission',flush=True)
