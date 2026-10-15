"""Compile an isolated mesh-Boolean utility from the existing pinned CGAL/libigl."""
from pathlib import Path
import shlex,subprocess,json,hashlib
root=Path(__file__).resolve().parents[2]
out=root/'build/heart-stem-contact-tools'; out.mkdir(exist_ok=True)
if (out/'mesh-boolean').exists(): raise FileExistsError('do not overwrite an experiment dependency')
flags=(root/'build/task-7v1-viewer/CMakeFiles/make_it_stand_viewer.dir/flags.make').read_text()
includes=shlex.split(next(x.split(' = ',1)[1] for x in flags.splitlines() if x.startswith('CXX_INCLUDES')))
source=Path(__file__).with_name('boolean_cli.cpp')
command=['/usr/bin/c++','-O1','-fmax-errors=2','-DNDEBUG','-DCGAL_HEADER_ONLY=1','-DCGAL_USE_CORE=1','-std=c++17','-frounding-math',*includes,
         str(source),'/home/lin/anaconda3/lib/libboost_thread.a','-pthread','-lgmp','-lmpfr','-o',str(out/'mesh-boolean')]
result=subprocess.run(command,capture_output=True,text=True)
(out/'compile.log').write_text(result.stdout+result.stderr)
if result.returncode:
    print('\n'.join(x for x in result.stderr.splitlines() if 'error:' in x))
    raise SystemExit(result.returncode)
(out/'build.json').write_text(json.dumps(dict(command=command,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
    binary_sha256=hashlib.sha256((out/'mesh-boolean').read_bytes()).hexdigest()),indent=2)+'\n')
