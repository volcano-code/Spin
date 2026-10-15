"""Independent GMP arithmetic module; no replacement of live libraries."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import sysconfig
import pybind11

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
output = Path(sys.argv[1])
output.mkdir(parents=True, exist_ok=False)
source = ROOT / 'bindings/cherchi/RationalMoments.cpp'
target = output / ('_cherchi_rational_moments' + sysconfig.get_config_var('EXT_SUFFIX'))
command = ['/usr/bin/c++', '-O3', '-std=c++17', '-shared', '-fPIC',
    '-I' + sysconfig.get_path('include'), '-I' + pybind11.get_include(),
    str(source), '-o', str(target), '-lgmpxx', '-lgmp']
with (output / 'build.log').open('x') as stream:
    subprocess.run(command, cwd=ROOT, stdout=stream, stderr=stream, check=True)
(output / 'receipt.json').write_text(json.dumps(dict(command=command,
    source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
    module_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),
    core_admission=False)))
