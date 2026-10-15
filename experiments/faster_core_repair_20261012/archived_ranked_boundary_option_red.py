"""Archived RED option test: rank candidate was not promoted after no speed benefit."""
import json
from pathlib import Path
import subprocess
import sys
import pytest

ROOT=Path(__file__).resolve().parents[2]


@pytest.mark.parametrize('ranked,filtered',[(False,False),(True,False),(True,True)])
def test_ranked_bbox_choice_is_bound_to_exact_audit_only(tmp_path,ranked,filtered):
    command=['/usr/bin/cmake','-S',str(ROOT),'-B',str(tmp_path/'build'),
        '-DMAKE_IT_STAND_BUILD_VIEWER=OFF','-DMAKE_IT_STAND_BUILD_REACHABILITY=OFF',
        '-DMAKE_IT_STAND_BUILD_SPIN_IT_FASTER=OFF','-DMAKE_IT_STAND_BUILD_PYTHON=OFF',
        '-DMAKE_IT_STAND_BUILD_EXACT_GEOMETRY_CHECKS=ON','-DCMAKE_EXPORT_COMPILE_COMMANDS=ON',
        '-DBoost_INCLUDE_DIR='+str(Path(sys.base_prefix)/'include'),
        '-DMAKE_IT_STAND_RANKED_BOUNDARY_BOXES='+('ON' if ranked else 'OFF'),
        '-DMAKE_IT_STAND_FILTERED_BOUNDARY_KERNEL='+('ON' if filtered else 'OFF')]
    result=subprocess.run(command,cwd=ROOT,capture_output=True,text=True,timeout=120)
    assert result.returncode==0,result.stdout+'\n'+result.stderr
    rows=json.loads((tmp_path/'build/compile_commands.json').read_text())
    audit=[r for r in rows if Path(r['file']).name=='ExactBoundaryAudit.cpp']
    assert len(audit)==1
    flag='-DMAKE_IT_STAND_RANKED_BOUNDARY_BOXES=1'
    assert (flag in audit[0]['command'])==ranked, 'requested exact rank bbox choice is not in compiled unit'
    assert all(flag not in r['command'] for r in rows if Path(r['file']).name!='ExactBoundaryAudit.cpp')


def test_ranked_bbox_without_exact_backend_is_rejected(tmp_path):
    result=subprocess.run(['/usr/bin/cmake','-S',str(ROOT),'-B',str(tmp_path/'build'),
        '-DMAKE_IT_STAND_BUILD_VIEWER=OFF','-DMAKE_IT_STAND_BUILD_REACHABILITY=OFF',
        '-DMAKE_IT_STAND_BUILD_SPIN_IT_FASTER=OFF','-DMAKE_IT_STAND_BUILD_EXACT_GEOMETRY_CHECKS=OFF',
        '-DMAKE_IT_STAND_RANKED_BOUNDARY_BOXES=ON'],cwd=ROOT,capture_output=True,text=True,timeout=120)
    assert result.returncode!=0,'exact rank bbox choice silently ignored without its backend'
    assert 'requires MAKE_IT_STAND_BUILD_EXACT_GEOMETRY_CHECKS' in ' '.join((result.stdout+result.stderr).split())
