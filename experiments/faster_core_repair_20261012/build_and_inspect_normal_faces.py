"""Compile a source-only probe, inspect fixed failed states without Boolean."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
output=Path(sys.argv[1]);output.mkdir(parents=True,exist_ok=False)
template=json.loads((ROOT/'build/heart-fixed-facet-20261011-v1/commands.json').read_text())[0]
source=HERE/'inspect_projective_normal_faces.cpp'
obj=output/'inspect.o';exe=output/'inspect-normal-faces'
compile=list(template)
for flag,value in (('-c',str(source)),('-o',str(obj)),('-MF',str(obj)+'.d'),('-MT',str(obj))):
    compile[compile.index(flag)+1]=value
compile.append('-I'+str(ROOT/'src'))
library=ROOT/'build/spin2014-objectives-quadric-20261007/libmake_it_stand_core.a'
link=[compile[0],str(obj),str(library),'-L/usr/lib/x86_64-linux-gnu','-lgmp','-lmpfr','-lpthread','-o',str(exe)]
with (output/'build.log').open('x') as log:
    subprocess.run(compile,cwd=ROOT,stdout=log,stderr=log,check=True)
    subprocess.run(link,cwd=ROOT,stdout=log,stderr=log,check=True)
states=json.loads((HERE/'frozen-heart-states.json').read_text())['records']
report={}
for name in ('heart-first-cap-failure','heart-first-small-face-rejection'):
    result=subprocess.run([str(exe),*[repr(x) for x in states[name]['coefficients']]],
                          capture_output=True,text=True,check=True)
    report[name]=json.loads(result.stdout)
with (output/'report.json').open('x') as stream:json.dump(report,stream,indent=2)
with (output/'receipt.json').open('x') as stream:
    json.dump(dict(source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        library_sha256=hashlib.sha256(library.read_bytes()).hexdigest(),
        binary_sha256=hashlib.sha256(exe.read_bytes()).hexdigest(),commands=[compile,link]),stream)
print({name:len(value['multiple_faces']) for name,value in report.items()},flush=True)
