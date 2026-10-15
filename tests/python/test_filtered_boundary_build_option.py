"""The explicit CMake choice must reach the actual audit compile command."""
import json
from pathlib import Path
import subprocess
import sys
import pytest

ROOT = Path(__file__).resolve().parents[2]


@pytest.mark.parametrize('filtered', [False, True])
def test_cmake_filtered_choice_reaches_only_boundary_translation_unit(tmp_path, filtered):
    command = ['/usr/bin/cmake', '-S', str(ROOT), '-B', str(tmp_path/'build'),
               '-DMAKE_IT_STAND_BUILD_VIEWER=OFF', '-DMAKE_IT_STAND_BUILD_REACHABILITY=OFF',
               '-DMAKE_IT_STAND_BUILD_EXACT_GEOMETRY_CHECKS=ON',
               '-DMAKE_IT_STAND_BUILD_SPIN_IT_FASTER=OFF', '-DMAKE_IT_STAND_BUILD_PYTHON=OFF',
               '-DCMAKE_EXPORT_COMPILE_COMMANDS=ON',
               '-DBoost_INCLUDE_DIR='+str(Path(sys.base_prefix)/'include'),
               '-DMAKE_IT_STAND_FILTERED_BOUNDARY_KERNEL='+('ON' if filtered else 'OFF')]
    configured = subprocess.run(command, cwd=ROOT, capture_output=True, text=True, timeout=120)
    assert configured.returncode == 0, configured.stdout+'\n'+configured.stderr
    rows = json.loads((tmp_path/'build/compile_commands.json').read_text())
    audit = [r for r in rows if Path(r['file']).name == 'ExactBoundaryAudit.cpp']
    assert len(audit) == 1
    flag = '-DMAKE_IT_STAND_FILTERED_BOUNDARY_KERNEL=1'
    assert (flag in audit[0]['command']) == filtered, 'requested exact filter not bound to audit object'
    assert all(flag not in r['command'] for r in rows if Path(r['file']).name != 'ExactBoundaryAudit.cpp')


def test_filtered_choice_without_exact_backend_is_rejected(tmp_path):
    configured = subprocess.run(['/usr/bin/cmake','-S',str(ROOT),'-B',str(tmp_path/'build'),
        '-DMAKE_IT_STAND_BUILD_VIEWER=OFF','-DMAKE_IT_STAND_BUILD_REACHABILITY=OFF',
        '-DMAKE_IT_STAND_BUILD_EXACT_GEOMETRY_CHECKS=OFF','-DMAKE_IT_STAND_BUILD_SPIN_IT_FASTER=OFF',
        '-DMAKE_IT_STAND_FILTERED_BOUNDARY_KERNEL=ON'],
        cwd=ROOT,capture_output=True,text=True,timeout=120)
    assert configured.returncode != 0, 'filtered selection silently built without its exact backend'
    assert 'requires MAKE_IT_STAND_BUILD_EXACT_GEOMETRY_CHECKS' in ' '.join(
        (configured.stdout+configured.stderr).split())
