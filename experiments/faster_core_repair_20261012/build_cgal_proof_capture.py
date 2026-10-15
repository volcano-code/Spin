"""Isolated binding-only build; all actual core objects retained by hash."""
import hashlib,json,subprocess,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=False)
prior=HERE/'projective-interval-build-job-02/output/receipt.json';receipt=json.loads(prior.read_text());link=list(receipt['commands'][-1]['argv'])
compile=list(json.loads((ROOT/'build/heart-fixed-facet-20261011-v1/commands.json').read_text())[0])
source=ROOT/'bindings/pybind_module.cpp';obj=out/'pybind_module.cpp.o';target=out/'_make_it_stand.cpython-313-x86_64-linux-gnu.so'
for flag,value in (('-c',str(source)),('-o',str(obj)),('-MF',str(obj)+'.d'),('-MT',str(obj))):compile[compile.index(flag)+1]=value
compile+=['-D_make_it_stand_EXPORTS','-isystem','/home/lin/anaconda3/include/python3.13','-isystem','/home/lin/anaconda3/lib/python3.13/site-packages/pybind11/include','-fvisibility=hidden']
old=[i for i,x in enumerate(link) if x.endswith('/bindings/pybind_module.cpp.o')]
if len(old)!=1:raise ValueError('AMBIGUOUS_BINDING_OBJECT')
link[old[0]]=str(obj);link[link.index('-o')+1]=str(target);start=time.monotonic()
with (out/'build.log').open('x') as log:
    subprocess.run(compile,cwd=ROOT,stdout=log,stderr=log,check=True,timeout=1800)
    subprocess.run(link,cwd=ROOT,stdout=log,stderr=log,check=True,timeout=300)
with (out/'receipt.json').open('x') as f:json.dump(dict(schema='native-binding-only-proof-overlay-v1',production_admission=False,prior_sha256=hashlib.sha256(prior.read_bytes()).hexdigest(),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),native_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),seconds=time.monotonic()-start,commands=[dict(argv=compile,object_sha256=hashlib.sha256(obj.read_bytes()).hexdigest()),dict(argv=link,retained_input_hashes={x:hashlib.sha256(Path(x).read_bytes()).hexdigest() for x in link if x.endswith(('.o','.a')) and Path(x).is_file()})]),f)
