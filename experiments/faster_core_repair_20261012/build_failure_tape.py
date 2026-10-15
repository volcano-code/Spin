"""Diagnostic-only independent overlay; no existing module is overwritten."""
import hashlib, json, subprocess, sys, time
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=False)
prior=HERE/'cgal-proof-capture-build-job-01/output/receipt.json'
link=list(json.loads(prior.read_text())['commands'][-1]['argv'])
template=json.loads((ROOT/'build/heart-fixed-facet-20261011-v1/commands.json').read_text())[0]
commands=[];start=time.monotonic()
with (out/'build.log').open('x') as log:
    for relative in ('src/spin/QuadricMomentEvaluator.cpp','bindings/pybind_module.cpp'):
        source=ROOT/relative;obj=out/(source.name+'.o');command=list(template)
        for flag,value in (('-c',str(source)),('-o',str(obj)),('-MF',str(obj)+'.d'),('-MT',str(obj))):
            command[command.index(flag)+1]=value
        if relative.startswith('bindings/'):
            command+=['-D_make_it_stand_EXPORTS','-isystem','/home/lin/anaconda3/include/python3.13',
                      '-isystem','/home/lin/anaconda3/lib/python3.13/site-packages/pybind11/include','-fvisibility=hidden']
        subprocess.run(command,cwd=ROOT,stdout=log,stderr=log,check=True,timeout=1800)
        old=[i for i,x in enumerate(link) if x.endswith('/'+source.name+'.o')]
        if len(old)>1:raise ValueError('AMBIGUOUS_OBJECT:'+source.name)
        if old:link[old[0]]=str(obj)
        else:link.insert(next(i for i,x in enumerate(link) if x.endswith('.a')),str(obj))
        commands.append(dict(argv=command,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                             object_sha256=hashlib.sha256(obj.read_bytes()).hexdigest()))
    target=out/'_make_it_stand.cpython-313-x86_64-linux-gnu.so'
    link[link.index('-o')+1]=str(target)
    subprocess.run(link,cwd=ROOT,stdout=log,stderr=log,check=True,timeout=300)
with (out/'receipt.json').open('x') as stream:
    json.dump(dict(schema='native-failure-tape-overlay-v1',production_admission=False,
                   prior_sha256=hashlib.sha256(prior.read_bytes()).hexdigest(),
                   native_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),seconds=time.monotonic()-start,
                   commands=commands+[dict(argv=link,retained_input_hashes={x:hashlib.sha256(Path(x).read_bytes()).hexdigest()
                      for x in link if x.endswith(('.o','.a')) and Path(x).is_file()})]),stream)
