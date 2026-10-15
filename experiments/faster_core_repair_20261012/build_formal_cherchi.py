import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[2]
out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=False)
commands=[['/usr/bin/cmake','-S',str(ROOT/'bindings/cherchi'),'-B',str(out/'build'),
           '-DCMAKE_BUILD_TYPE=Release','-DPython_EXECUTABLE='+sys.executable],
          ['/usr/bin/cmake','--build',str(out/'build'),'--parallel','4']]
with (out/'build.log').open('x') as stream:
    for command in commands:
        subprocess.run(command,stdout=stream,stderr=stream,check=True,timeout=1800)
modules=list((out/'build/python').glob('_cherchi_core_boolean*.so'))
if len(modules)!=1:raise ValueError('EXPECTED_ONE_NATIVE_CHERCHI_MODULE')
observed={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
          for p in (ROOT/'bindings/cherchi').glob('*') if p.is_file()}
with (out/'receipt.json').open('x') as stream:
    json.dump(dict(commands=commands,module=str(modules[0]),
        sha256=hashlib.sha256(modules[0].read_bytes()).hexdigest(),sources=observed,
        production_admission=False),stream,sort_keys=True)
