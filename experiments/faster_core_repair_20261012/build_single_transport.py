"""Independent overlay, enabled only for proven identical transport objects."""
import hashlib,json,subprocess,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=False)
transport_path=HERE/'fixed-design-cache-build-job-01/output/receipt.json'
transport=json.loads(transport_path.read_text())
rows=[r for r in transport['commands'] if '-c' in r['argv']]
def semantic(argv):
    result=[];i=0
    while i<len(argv):
        if argv[i] in ('-o','-MF','-MT'):i+=2;continue
        if not argv[i].startswith('-Ddifferentiate'):result.append(argv[i])
        i+=1
    return result
if len(rows)!=2 or rows[0]['source_sha256']!=rows[1]['source_sha256'] or semantic(rows[0]['argv'])!=semantic(rows[1]['argv']):
    raise ValueError('TRANSPORT_IMPLEMENTATIONS_NOT_IDENTICAL')
for row in rows:
    obj=Path(row['argv'][row['argv'].index('-o')+1])
    if hashlib.sha256(obj.read_bytes()).hexdigest()!=row['object_sha256']:raise ValueError('TRANSPORT_OBJECT_CHANGED')
prior=HERE/'component-selector-build-job-01/output/receipt.json'
link=list(json.loads(prior.read_text())['commands'][-1]['argv'])
templates=json.loads((ROOT/'build/heart-exact-input-20261011-v1/commands.json').read_text())
command=list(next(c for c in templates if '-c' in c and c[c.index('-c')+1].endswith('/FullRankMomentEvaluator.cpp')))
source=ROOT/'src/spin/FullRankMomentEvaluator.cpp';obj=out/'FullRankMomentEvaluator.cpp.o'
for flag,value in (('-c',str(source)),('-o',str(obj)),('-MF',str(obj)+'.d'),('-MT',str(obj))):command[command.index(flag)+1]=value
command.append('-DMAKE_IT_STAND_IDENTICAL_TRANSPORT_OVERLAY=1')
for row in rows:
    path=row['argv'][row['argv'].index('-o')+1]
    if path not in link:raise ValueError('VERIFIED_TRANSPORT_NOT_LINKED')
indices=[i for i,x in enumerate(link) if x.endswith('/FullRankMomentEvaluator.cpp.o')]
if len(indices)!=1:raise ValueError('AMBIGUOUS_FULLRANK_OBJECT')
link[indices[0]]=str(obj);target=out/'_make_it_stand.cpython-313-x86_64-linux-gnu.so';link[link.index('-o')+1]=str(target)
with (out/'build.log').open('x') as log:
    subprocess.run(command,cwd=ROOT,stdout=log,stderr=log,check=True)
    subprocess.run(link,cwd=ROOT,stdout=log,stderr=log,check=True)
(out/'receipt.json').write_text(json.dumps(dict(commands=[command,link],
    source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
    transport_receipt_sha256=hashlib.sha256(transport_path.read_bytes()).hexdigest(),
    native_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),core_admission=False)))
