"""Build the production CMake arithmetic target in an independent directory."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
output = Path(sys.argv[1])
output.mkdir(parents=True, exist_ok=False)
commands = [
    ['/usr/bin/cmake', '-S', str(ROOT / 'bindings/cherchi'), '-B', str(output / 'build'),
     '-DCMAKE_BUILD_TYPE=Release', '-DPython_EXECUTABLE=' + sys.executable,
     '-DMAKE_IT_STAND_CHERCHI_RATIONAL_ARITHMETIC=ON'],
    ['/usr/bin/cmake', '--build', str(output / 'build'), '--target', '_cherchi_rational_moments', '--parallel', '4']]
with (output / 'build.log').open('x') as stream:
    for command in commands:
        subprocess.run(command, cwd=ROOT, stdout=stream, stderr=stream, check=True)
modules = list((output / 'build/python').glob('_cherchi_rational_moments*.so'))
if len(modules) != 1:
    raise ValueError('EXPECTED_ONE_ARITHMETIC_MODULE')
(output / 'receipt.json').write_text(json.dumps(dict(commands=commands, module=str(modules[0]),
    module_sha256=hashlib.sha256(modules[0].read_bytes()).hexdigest(),
    source_sha256=hashlib.sha256((ROOT / 'bindings/cherchi/RationalMoments.cpp').read_bytes()).hexdigest(),
    cmake_sha256=hashlib.sha256((ROOT / 'bindings/cherchi/CMakeLists.txt').read_bytes()).hexdigest(),
    core_admission=False)))
