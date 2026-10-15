"""Link a focused selector contract test against the exact candidate objects."""
import hashlib,json,subprocess,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=False)
receipt=HERE/'component-selector-build-job-01/output/receipt.json'
prior=json.loads(receipt.read_text())
source=ROOT/'tests/quadric_component_selector_regression.cpp'
obj=out/'selector_test.o';binary=out/'selector_test'
command=list(json.loads((ROOT/'build/heart-fixed-facet-20261011-v1/commands.json').read_text())[0])
for flag,value in (('-c',str(source)),('-o',str(obj)),('-MF',str(obj)+'.d'),('-MT',str(obj))):
    command[command.index(flag)+1]=value
link=[x for x in prior['commands'][-1]['argv'] if x!='-shared' and not (x.endswith('.o') and ('/bindings/' in x or x.endswith('/pybind_module.cpp.o')))]
link[link.index('-o')+1]=str(binary)
link.insert(next(i for i,x in enumerate(link) if x.endswith('.o')),str(obj))
with (out/'build.log').open('x') as log:
    subprocess.run(command,cwd=ROOT,stdout=log,stderr=log,check=True)
    subprocess.run(link,cwd=ROOT,stdout=log,stderr=log,check=True)
run=subprocess.run([str(binary)],capture_output=True,text=True)
(out/'report.json').write_text(json.dumps(dict(returncode=run.returncode,stdout=run.stdout,stderr=run.stderr,
    source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
    prior_sha256=hashlib.sha256(receipt.read_bytes()).hexdigest(),commands=[command,link],core_admission=False)))
if run.returncode:raise SystemExit(run.returncode)
