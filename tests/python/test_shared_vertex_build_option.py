"""Formal explicit selection of the exact shared-vertex cone predicate."""
import json
from pathlib import Path
import subprocess
import sys
import pytest

ROOT=Path(__file__).resolve().parents[2]


@pytest.mark.parametrize('enabled',[False,True])
def test_shared_vertex_selection_reaches_only_audit_unit(tmp_path,enabled):
    command=['/usr/bin/cmake','-S',str(ROOT),'-B',str(tmp_path/'build'),
        '-DMAKE_IT_STAND_BUILD_VIEWER=OFF','-DMAKE_IT_STAND_BUILD_REACHABILITY=OFF',
        '-DMAKE_IT_STAND_BUILD_SPIN_IT_FASTER=OFF','-DMAKE_IT_STAND_BUILD_PYTHON=OFF',
        '-DMAKE_IT_STAND_BUILD_EXACT_GEOMETRY_CHECKS=ON','-DCMAKE_EXPORT_COMPILE_COMMANDS=ON',
        '-DBoost_INCLUDE_DIR='+str(Path(sys.base_prefix)/'include'),
        '-DMAKE_IT_STAND_SHARED_VERTEX_CONE_PREDICATE='+('ON' if enabled else 'OFF')]
    result=subprocess.run(command,cwd=ROOT,capture_output=True,text=True,timeout=120)
    assert result.returncode==0,result.stdout+'\n'+result.stderr
    rows=json.loads((tmp_path/'build/compile_commands.json').read_text())
    audit=[r for r in rows if Path(r['file']).name=='ExactBoundaryAudit.cpp']
    assert len(audit)==1
    flag='-DMAKE_IT_STAND_SHARED_VERTEX_CONE_PREDICATE=1'
    assert (flag in audit[0]['command'])==enabled,'exact shared-vertex selection not bound to audit object'
    assert all(flag not in r['command'] for r in rows if Path(r['file']).name!='ExactBoundaryAudit.cpp')


def test_shared_vertex_selection_without_exact_checks_is_rejected(tmp_path):
    result=subprocess.run(['/usr/bin/cmake','-S',str(ROOT),'-B',str(tmp_path/'build'),
        '-DMAKE_IT_STAND_BUILD_VIEWER=OFF','-DMAKE_IT_STAND_BUILD_REACHABILITY=OFF',
        '-DMAKE_IT_STAND_BUILD_SPIN_IT_FASTER=OFF','-DMAKE_IT_STAND_BUILD_EXACT_GEOMETRY_CHECKS=OFF',
        '-DMAKE_IT_STAND_SHARED_VERTEX_CONE_PREDICATE=ON'],cwd=ROOT,capture_output=True,text=True,timeout=120)
    assert result.returncode!=0,'shared-vertex selection silently ignored without exact backend'
    assert 'requires MAKE_IT_STAND_BUILD_EXACT_GEOMETRY_CHECKS' in ' '.join((result.stdout+result.stderr).split())
