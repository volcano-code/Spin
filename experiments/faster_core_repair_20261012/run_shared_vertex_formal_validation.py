"""Ordinary-source predicate variant: fixed maps and relevant regressions."""
import json
from pathlib import Path
import subprocess
import sys
from verify_formal_native_rebuild import main as verify_maps

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def main(output):
    output.mkdir(parents=True,exist_ok=False)
    verify_maps(output/'fixed')
    names=('test_formal_native_module_load','test_canonical_transport_observer','test_filtered_boundary_build_option',
        'test_shared_vertex_build_option','test_exact_filtered_boundary','test_exact_boundary_audit',
        'test_rational_boundary_audit','test_rational_solver_snapshot','test_cherchi_rational_native',
        'test_cherchi_native_arithmetic_adapter','test_cherchi_inprocess_core','test_cherchi_rational_boundary_profile',
        'test_cherchi_run_journal','test_cherchi_formal_run_evidence','test_faster_material_bridge')
    command=[sys.executable,'-m','pytest',*[str(ROOT/'tests/python'/(name+'.py')) for name in names],'-q']
    with (output/'tests.log').open('x') as log:
        subprocess.run(command,cwd=ROOT,stdout=log,stderr=log,check=True)
    fixed=json.loads((output/'fixed/report.json').read_text())
    (output/'report.json').write_text(json.dumps(dict(passed=fixed['passed'],command=command,
        core_admission=False,scope='Ordinary-source shared-vertex variant fixed maps and relevant core regressions; full solve pending')))


if __name__=='__main__':main(Path(sys.argv[1]))
