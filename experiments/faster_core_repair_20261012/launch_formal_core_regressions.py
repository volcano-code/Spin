"""Relevant formal core regressions under the same serialized resource limits."""
from pathlib import Path
import json
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job

if not json.loads((HERE/'formal-core-benchmark-small-seed-0-job-01/output/report.json').read_text())['material_closed_loop_passed']:
    raise ValueError('FORMAL_SMALL_CLOSURE_NOT_PASSED')
build=HERE/'formal-native-core-build-job-04'
native=build/'output/build/python/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
bridge=HERE/'formal-cherchi-build-job-02/output/build/python/_cherchi_core_boolean.cpython-313-x86_64-linux-gnu.so'
arithmetic=HERE/'formal-rational-arithmetic-build-job-01/output/build/python/_cherchi_rational_moments.cpython-313-x86_64-linux-gnu.so'
names=('test_formal_native_module_load','test_canonical_transport_observer','test_filtered_boundary_build_option',
       'test_exact_filtered_boundary','test_exact_boundary_audit','test_rational_boundary_audit',
       'test_rational_solver_snapshot','test_cherchi_rational_native','test_cherchi_native_arithmetic_adapter',
       'test_cherchi_inprocess_core','test_cherchi_rational_boundary_profile','test_cherchi_run_journal',
       'test_cherchi_formal_run_evidence','test_faster_material_bridge')
tests=[ROOT/'tests/python'/(name+'.py') for name in names]
files=dict(launcher=Path(__file__),native=native,bridge=bridge,arithmetic=arithmetic,
    receipt=build/'output/receipt.json',spin=ROOT/'python/make_it_stand/spin.py',
    verification=HERE/'formal-core-benchmark-small-seed-0-job-01/output/run/record-verification.json')
for path in tests:files[path.stem]=path
for path in (ROOT/'python/make_it_stand/verification').glob('*.py'):
    files['verification_'+path.name]=path
generation=sys.argv[1] if len(sys.argv)>1 else '01'
if generation not in ('01','02','03'):
    raise ValueError('UNREGISTERED_FORMAL_REGRESSION_GENERATION')
result=run_source_job(HERE/('formal-core-regressions-job-'+generation),
    ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
     'MAKE_IT_STAND_FORMAL_EVIDENCE_RUN='+str(HERE/'formal-core-benchmark-small-seed-0-job-01/output/run'),
     '/home/lin/anaconda3/bin/python3','-m','pytest',*[str(p) for p in tests],'-q'],files,
    environment={'PYTHONPATH':':'.join(map(str,[ROOT/'python',native.parent,bridge.parent,arithmetic.parent]))})
print(result.status,result.returncode,flush=True)
