import json
from pathlib import Path
import shlex
import subprocess
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
output=Path(sys.argv[1]);output.mkdir(parents=True,exist_ok=False)
rows=json.loads((HERE/'formal-native-core-build-job-04/output/build/compile_commands.json').read_text())
command=shlex.split(next(r['command'] for r in rows if Path(r['file']).name=='ExactBoundaryAudit.cpp'))
command[command.index('-c')+1]=str(ROOT/'tests/exact_shared_vertex_predicate_regression.cpp')
command[command.index('-o')+1]=str(output/'test.o')
filtered=len(sys.argv)>2 and sys.argv[2]=='filtered'
if filtered:command.append('-DTEST_FILTERED_SHARED_VERTEX_KERNEL=1')
commands=[command,['/usr/bin/c++','-o',str(output/'vertex-test'),str(output/'test.o'),'-lgmp','-lmpfr']]
with (output/'build.log').open('x') as log:
    for item in commands:subprocess.run(item,cwd=ROOT,stdout=log,stderr=log,check=True)
with (output/'test.log').open('x') as log:
    subprocess.run([str(output/'vertex-test')],cwd=ROOT,stdout=log,stderr=log,check=True)
(output/'report.json').write_text(json.dumps(dict(passed=True,commands=commands,filtered=filtered,core_admission=False,
    scope='Shared-vertex exact predicate versus independent construction; not solver admission')))
