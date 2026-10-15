"""Isolated single-object overlay; do not replace any frozen extension."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time


def main(output):
    here=Path(__file__).resolve().parent;root=here.parents[1]
    output.mkdir(parents=True,exist_ok=False)
    template=json.loads((root/'build/heart-fixed-facet-20261011-v1/commands.json').read_text())[0]
    prior=json.loads((here/'boundary-proof-reuse-build-job-01/output/receipt.json').read_text())
    link=list(prior['commands'][-1]['argv'])
    source=root/'src/spin/QuadricHollowingSolver.cpp';obj=output/'QuadricHollowingSolver.cpp.o'
    command=list(template)
    for flag,value in (('-c',str(source)),('-o',str(obj)),('-MF',str(obj)+'.d'),('-MT',str(obj))):
        command[command.index(flag)+1]=value
    start=time.monotonic()
    with (output/'build.log').open('x') as log:
        subprocess.run(command,cwd=root,stdout=log,stderr=log,check=True,timeout=1800)
        old=[i for i,x in enumerate(link) if x.endswith('/QuadricHollowingSolver.cpp.o')]
        if len(old)>1:raise ValueError('AMBIGUOUS_SOLVER_OBJECT')
        if old:link[old[0]]=str(obj)
        else:link.insert(next(i for i,x in enumerate(link) if x.endswith('.a')),str(obj))
        target=output/'_make_it_stand.cpython-313-x86_64-linux-gnu.so'
        link[link.index('-o')+1]=str(target)
        subprocess.run(link,cwd=root,stdout=log,stderr=log,check=True,timeout=300)
    receipt=dict(schema='native-overlay-build-v1',production_admission=False,
        source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        commands=[dict(argv=command,object_sha256=hashlib.sha256(obj.read_bytes()).hexdigest()),
                  dict(argv=link,retained_input_hashes={x:hashlib.sha256(Path(x).read_bytes()).hexdigest()
                       for x in link if x.endswith(('.o','.a')) and Path(x).is_file()})],
        native_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),seconds=time.monotonic()-start)
    with (output/'receipt.json').open('x') as stream:json.dump(receipt,stream)


if __name__=='__main__':main(Path(sys.argv[1]))
