"""Benchmark only the successful ordinary source build, explicitly configured."""
from pathlib import Path
import json
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job

scale,seed=sys.argv[1],int(sys.argv[2])
repeat=len(sys.argv)>3 and sys.argv[3]=='repeat'
if len(sys.argv)>4 or (len(sys.argv)>3 and not repeat) or scale not in ('small','faces-12288','faces-98304') or seed not in range(10):
    raise ValueError('UNREGISTERED_FORMAL_BENCHMARK')
if not json.loads((HERE/'formal-native-rebuild-check-job-01/output/report.json').read_text())['passed']:
    raise ValueError('FORMAL_FIXED_MAP_CHECK_NOT_PASSED')
for gate in ('formal-native-saved-derivatives-job-01','formal-registered-geometry-check-job-01'):
    if not json.loads((HERE/gate/'output/report.json').read_text())['passed']:
        raise ValueError('FORMAL_CORE_GATE_NOT_PASSED:'+gate)
if scale!='small' or seed!=0 or repeat:
    regression=json.loads((HERE/'formal-core-regressions-job-03/result.json').read_text())
    if regression['status']!='PROCESS_EXITED' or regression['returncode']!=0:
        raise ValueError('FORMAL_CORE_REGRESSIONS_NOT_PASSED')
build=HERE/'formal-native-core-build-job-04'
native=build/'output/build/python/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
bridge=HERE/'formal-cherchi-build-job-02/output/build/python/_cherchi_core_boolean.cpython-313-x86_64-linux-gnu.so'
arithmetic=HERE/'formal-rational-arithmetic-build-job-01/output/build/python/_cherchi_rational_moments.cpython-313-x86_64-linux-gnu.so'
files=dict(launcher=Path(__file__),worker=HERE/'run_core_benchmark.py',native=native,bridge=bridge,
    arithmetic_module=arithmetic,native_receipt=build/'output/receipt.json',
    fixed_map_check=HERE/'formal-native-rebuild-check-job-01/output/report.json',
    saved_map_check=HERE/'formal-native-saved-derivatives-job-01/output/report.json',
    geometry_check=HERE/'formal-registered-geometry-check-job-01/output/report.json',
    arithmetic_receipt=HERE/'formal-rational-arithmetic-build-job-01/output/receipt.json',
    config=HERE/'core-benchmark-inputs-job-01/output'/(scale+'.json'),
    manifest=HERE/'core-benchmark-inputs-job-01/output/manifest.json',spin=ROOT/'python/make_it_stand/spin.py')
if scale!='small' or seed!=0 or repeat:
    files['regression_result']=HERE/'formal-core-regressions-job-03/result.json'
    files['regression_log']=HERE/'formal-core-regressions-job-03/stdout.log'
for name in ('cherchi_inprocess','cherchi_run_evidence','cherchi_run_journal','cherchi_motion_differential',
    'faster_material','spin_problem','spin_bundle','spin_artifact','spin_common','mass_moments','rational_text'):
    files[name]=ROOT/f'python/make_it_stand/verification/{name}.py'
suffix='repeat-01' if repeat else 'job-01'
result=run_source_job(HERE/f'formal-core-benchmark-{scale}-seed-{seed}-{suffix}',
    ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
     '/home/lin/anaconda3/bin/python3','-B',str(files['worker']),'{job}/output',scale,str(seed),
     'NATIVE_GMP_SURFACE_V1','EXACT_RATIONAL_BOUNDARY_V1','EVALUATION_CHECKPOINT_V1'],files,
    environment={'PYTHONPATH':':'.join(map(str,[ROOT/'python',native.parent,bridge.parent,arithmetic.parent]))})
print(result.status,result.returncode,flush=True)
