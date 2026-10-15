"""A linked extension must resolve all declared canonical transport overloads."""
import os
import subprocess
import sys


def test_native_extension_loads_and_exposes_core_evaluator():
    result=subprocess.run([sys.executable,'-c',
        'import _make_it_stand as n; '
        'assert callable(n.evaluate_spin_it_faster_fullrank_quadric); '
        'assert callable(n.evaluate_spin_it_faster_source_motion)'],
        env=dict(os.environ),capture_output=True,text=True,timeout=30)
    assert result.returncode==0, result.stdout+'\n'+result.stderr
