"""Rebuild only the audited bbox unit, retaining verified ordinary source objects."""
import hashlib
import json
from pathlib import Path
import shlex
import subprocess
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
output=Path(sys.argv[1]);output.mkdir(parents=True,exist_ok=False)
variant=sys.argv[2] if len(sys.argv)>2 else 'ranked'
if variant not in ('ranked','sharedvertex'):raise ValueError('UNREGISTERED_BOUNDARY_CANDIDATE')
header_name='ExactCoordinateRankBox.h' if variant=='ranked' else 'ExactSharedVertexIntersection.h'
definition='MAKE_IT_STAND_RANKED_BOUNDARY_BOXES' if variant=='ranked' else 'MAKE_IT_STAND_SHARED_VERTEX_CONE_PREDICATE'
base=HERE/'formal-native-core-build-job-04/output/build'
receipt=HERE/'formal-native-core-build-job-04/output/receipt.json'
prior=json.loads(receipt.read_text())
rows=json.loads((base/'compile_commands.json').read_text())
command=shlex.split(next(r['command'] for r in rows if Path(r['file']).name=='ExactBoundaryAudit.cpp'))
obj=output/'ExactBoundaryAudit.cpp.o'
command[command.index('-o')+1]=str(obj)
command.append('-D'+definition+'=1')
link=shlex.split((base/'CMakeFiles/_make_it_stand.dir/link.txt').read_text())
target=output/'_make_it_stand.cpython-313-x86_64-linux-gnu.so'
link[link.index('-o')+1]=str(target)
link.insert(next(i for i,x in enumerate(link) if x.endswith('.a')),str(obj))
retained={}
for token in link:
    path=Path(token)
    if not path.is_absolute():path=base/path
    if path.is_file() and path.suffix in ('.o','.a') and path!=obj:
        retained[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
with (output/'build.log').open('x') as log:
    subprocess.run(command,cwd=base,stdout=log,stderr=log,check=True)
    subprocess.run(link,cwd=base,stdout=log,stderr=log,check=True)
(output/'receipt.json').write_text(json.dumps(dict(commands=[command,link],link_directory=str(base),
    prior_receipt_sha256=hashlib.sha256(receipt.read_bytes()).hexdigest(),
    retained_object_sha256=retained,native_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),
    source_sha256=hashlib.sha256((ROOT/'src/verification/ExactBoundaryAudit.cpp').read_bytes()).hexdigest(),
    header_sha256=hashlib.sha256((ROOT/'include/verification'/header_name).read_bytes()).hexdigest(),variant=variant,
    core_admission=False)))
