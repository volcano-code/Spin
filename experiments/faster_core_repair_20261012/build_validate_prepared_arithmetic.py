"""Independent candidate build and affected arithmetic tests, never installation."""
from pathlib import Path
import json
import os
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
output = Path(sys.argv[1])
output.mkdir(parents=True, exist_ok=False)
subprocess.run([sys.executable, '-B', str(HERE / 'build_formal_rational_arithmetic.py'),
                str(output / 'candidate')], cwd=ROOT, check=True)
environment = dict(os.environ)
environment['LD_PRELOAD'] = '/usr/lib/x86_64-linux-gnu/libstdc++.so.6'
environment['PYTHONPATH'] = ':'.join(map(str, [output / 'candidate/build/python', ROOT / 'python', HERE]))
command = [sys.executable, '-B', '-m', 'pytest', '-q',
           str(HERE / 'test_prepared_rational_boundary.py'),
           str(ROOT / 'tests/python/test_cherchi_rational_native.py')]
with (output / 'tests.log').open('x') as stream:
    result = subprocess.run(command, cwd=ROOT, env=environment,
                            stdout=stream, stderr=subprocess.STDOUT)
(output / 'report.json').write_text(json.dumps({
    'passed': result.returncode == 0, 'test_command': command,
    'test_returncode': result.returncode, 'production_integrated': False,
    'performance_passed': False, 'full_core_admission': False}, indent=2))
raise SystemExit(result.returncode)
