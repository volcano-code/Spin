"""Retain source/motion repairs, overlay connected exact material proof only."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
output=Path(sys.argv[1]);output.mkdir(parents=True,exist_ok=False)
template=json.loads((ROOT/'build/heart-fixed-facet-20261011-v1/commands.json').read_text())[0]
prior=json.loads((HERE/'projective-belt-build-job-04/output/receipt.json').read_text())
link=list(prior['commands'][-1]['argv'])
source=ROOT/'src/spin/QuadricMomentEvaluator.cpp';obj=output/'QuadricMomentEvaluator.cpp.o'
compile=list(template)
for flag,value in (('-c',str(source)),('-o',str(obj)),('-MF',str(obj)+'.d'),('-MT',str(obj))):
    compile[compile.index(flag)+1]=value
start=time.monotonic()
with (output/'build.log').open('x') as log:
    subprocess.run(compile,cwd=ROOT,stdout=log,stderr=log,check=True,timeout=1800)
    matches=[i for i,x in enumerate(link) if x.endswith('/QuadricMomentEvaluator.cpp.o')]
    if len(matches)!=1:raise ValueError('EXPECTED_ONE_VALUE_EVALUATOR_OBJECT')
    link[matches[0]]=str(obj)
    target=output/'_make_it_stand.cpython-313-x86_64-linux-gnu.so';link[link.index('-o')+1]=str(target)
    subprocess.run(link,cwd=ROOT,stdout=log,stderr=log,check=True,timeout=300)
with (output/'receipt.json').open('x') as stream:
    json.dump(dict(schema='native-overlay-build-v1',production_admission=False,
        source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        native_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),seconds=time.monotonic()-start,
        commands=[dict(argv=compile,object_sha256=hashlib.sha256(obj.read_bytes()).hexdigest()),
                  dict(argv=link,retained_input_hashes={x:hashlib.sha256(Path(x).read_bytes()).hexdigest()
                      for x in link if x.endswith(('.o','.a')) and Path(x).is_file()})]),stream)
