"""Independent exact-audit predicate variant; all other objects are retained."""
import hashlib,json,subprocess,sys
from pathlib import Path

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
output=Path(sys.argv[1]);output.mkdir(parents=True,exist_ok=False)
prior=HERE/'rational-boundary-binding-build-job-01/output/receipt.json'
link=list(json.loads(prior.read_text())['commands'][-1])
template=json.loads((ROOT/'build/heart-fixed-facet-20261011-v1/commands.json').read_text())[0]
source=ROOT/'src/verification/ExactBoundaryAudit.cpp';obj=output/'ExactBoundaryAudit.cpp.o'
command=list(template)
for flag,value in (('-c',str(source)),('-o',str(obj)),('-MF',str(obj)+'.d'),('-MT',str(obj))):
    command[command.index(flag)+1]=value
command.append('-DMAKE_IT_STAND_FILTERED_BOUNDARY_KERNEL=1')
if any(x.endswith('/ExactBoundaryAudit.cpp.o') for x in link):raise ValueError('UNEXPECTED_EXPLICIT_AUDIT_OBJECT')
link.insert(next(i for i,x in enumerate(link) if x.endswith('.a')),str(obj))
target=output/'_make_it_stand.cpython-313-x86_64-linux-gnu.so';link[link.index('-o')+1]=str(target)
with (output/'build.log').open('x') as log:
    subprocess.run(command,cwd=ROOT,stdout=log,stderr=log,check=True)
    subprocess.run(link,cwd=ROOT,stdout=log,stderr=log,check=True)
(output/'receipt.json').write_text(json.dumps(dict(commands=[command,link],
    source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
    prior_sha256=hashlib.sha256(prior.read_bytes()).hexdigest(),
    native_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),
    retained_object_sha256={x:hashlib.sha256(Path(x).read_bytes()).hexdigest()
        for x in link if x.endswith(('.o','.a')) and Path(x).is_file()},core_admission=False)))
