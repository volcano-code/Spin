import hashlib
import json
from pathlib import Path
import shlex
import subprocess
import sys
import sysconfig
import pybind11

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent


def main(output):
    output=Path(output);output.mkdir(parents=True,exist_ok=False)
    reference=ROOT/'build/cherchi-reference-singlethread-20260917'
    flags={}
    for row in (reference/'CMakeFiles/cherchi_probe.dir/flags.make').read_text().splitlines():
        if row.startswith('CXX_') and ' = ' in row:
            key,value=row.split(' = ',1);flags[key]=shlex.split(value)
    suffix=sysconfig.get_config_var('EXT_SUFFIX')
    target=output/('_cherchi_core_boolean'+suffix)
    libraries=[(reference/token).resolve() for token in shlex.split((reference/'CMakeFiles/cherchi_probe.dir/link.txt').read_text())
        if token.endswith('.a') or '.so' in token and not token.startswith('-')]
    command=['/usr/bin/c++',*flags['CXX_DEFINES'],*flags['CXX_INCLUDES'],*flags['CXX_FLAGS'],
        '-shared','-fPIC','-I'+sysconfig.get_path('include'),'-I'+pybind11.get_include(),
        str(HERE/'cherchi_inprocess.cpp'),'-o',str(target),*[str(p) for p in libraries],
        '-Wl,-rpath,'+str(libraries[0].parent)]
    with (output/'build.log').open('x') as log:
        subprocess.run(command,cwd=ROOT,stdout=log,stderr=log,check=True,timeout=1800)
    with (output/'receipt.json').open('x') as stream:
        json.dump(dict(command=command,native_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),
            source_sha256=hashlib.sha256((HERE/'cherchi_inprocess.cpp').read_bytes()).hexdigest(),
            libraries={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in libraries},production_admission=False),stream)


if __name__=='__main__':main(sys.argv[1])
