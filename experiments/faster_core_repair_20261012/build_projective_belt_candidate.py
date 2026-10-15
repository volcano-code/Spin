"""Isolated triangulator overlay for the connected finite-normal belt."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
output=Path(sys.argv[1]);output.mkdir(parents=True,exist_ok=False)
template=json.loads((ROOT/'build/heart-fixed-facet-20261011-v1/commands.json').read_text())[0]
prior_path=Path(sys.argv[2]) if len(sys.argv)==3 else HERE/'newton-capture-build-job-01/output/receipt.json'
prior=json.loads(prior_path.read_text())
link=list(prior['commands'][-1]['argv'])
start=time.monotonic()
receipts=[]
with (output/'build.log').open('x') as log:
    for stem in ('QuadricTriangulator','QuadricSourceMotion'):
        source=ROOT/f'src/spin/{stem}.cpp';obj=output/f'{stem}.cpp.o'
        compile=list(template)
        for flag,value in (('-c',str(source)),('-o',str(obj)),('-MF',str(obj)+'.d'),('-MT',str(obj))):
            compile[compile.index(flag)+1]=value
        subprocess.run(compile,cwd=ROOT,stdout=log,stderr=log,check=True,timeout=1800)
        old=[i for i,x in enumerate(link) if x.endswith('/'+stem+'.cpp.o')]
        if len(old)>1:raise ValueError('AMBIGUOUS_SOURCE_OBJECT:'+stem)
        if old:link[old[0]]=str(obj)
        else:link.insert(next(i for i,x in enumerate(link) if x.endswith('.a')),str(obj))
        receipts.append(dict(argv=compile,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                             object_sha256=hashlib.sha256(obj.read_bytes()).hexdigest()))
    target=output/'_make_it_stand.cpython-313-x86_64-linux-gnu.so'
    link[link.index('-o')+1]=str(target)
    subprocess.run(link,cwd=ROOT,stdout=log,stderr=log,check=True,timeout=300)
with (output/'receipt.json').open('x') as stream:
    json.dump(dict(schema='native-overlay-build-v1',production_admission=False,
        native_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),seconds=time.monotonic()-start,
        commands=[*receipts,
                  dict(argv=link,retained_input_hashes={x:hashlib.sha256(Path(x).read_bytes()).hexdigest()
                     for x in link if x.endswith(('.o','.a')) and Path(x).is_file()})]),stream)
