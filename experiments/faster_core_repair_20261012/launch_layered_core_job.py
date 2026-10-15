"""Serial resource-supervised layered validation and unchanged seeded benchmarks."""
from pathlib import Path
import json
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job

mode = sys.argv[1]
build = HERE/'formal-native-core-build-job-05'
native = build/'output/build/python/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
bridge = HERE/'formal-cherchi-build-job-02/output/build/python/_cherchi_core_boolean.cpython-313-x86_64-linux-gnu.so'
arithmetic = HERE/'formal-rational-arithmetic-build-job-01/output/build/python/_cherchi_rational_moments.cpython-313-x86_64-linux-gnu.so'
files = dict(launcher=Path(__file__), native=native, bridge=bridge, arithmetic=arithmetic,
    receipt=build/'output/receipt.json', spin=ROOT/'python/make_it_stand/spin.py')
for path in (ROOT/'python/make_it_stand/verification').glob('*.py'):
    files['module_'+path.name] = path
if mode == 'validate':
    generation = sys.argv[2] if len(sys.argv) == 3 else '01'
    if len(sys.argv) not in (2,3) or generation not in ('01','02','03','04','05','06'):
        raise ValueError('INVALID_VALIDATION_ARGUMENTS')
    files.update(worker=HERE/'run_layered_validation.py', maps=HERE/'verify_formal_native_rebuild.py',
        registered=HERE/'projective-formal-pipeline-job-01/output/registered-inputs.json')
    for state in ('ellipsoid','one_sheet','two_sheet'):
        files[state]=HERE/'projective-formal-pipeline-job-01/output/fixed'/(state+'-material.json')
    for path in (ROOT/'tests/python').glob('test_*.py'): files['test_'+path.name]=path
    job=HERE/f'layered-formal-validation-job-{generation}'
    args=[]
elif mode in ('benchmark','audit'):
    scale, seed = sys.argv[2], int(sys.argv[3])
    generation = sys.argv[4] if len(sys.argv)>4 else '01'
    repeat = len(sys.argv)==6 and sys.argv[5]=='repeat'
    if (len(sys.argv) not in (4,5,6) or generation not in ('01','02','03','04','05','06') or
        (len(sys.argv)==6 and not repeat) or
        scale not in ('small','faces-12288','faces-98304') or seed not in range(10)):
        raise ValueError('UNREGISTERED_LAYERED_BENCHMARK')
    validation=HERE/f'layered-formal-validation-job-{generation}/output/report.json'
    if not json.loads(validation.read_text())['passed']: raise ValueError('LAYERED_VALIDATION_NOT_PASSED')
    files.update(worker=HERE/'run_core_benchmark.py', validation=validation,
        config=HERE/'core-benchmark-inputs-job-01/output'/(scale+'.json'),
        manifest=HERE/'core-benchmark-inputs-job-01/output/manifest.json')
    suffix = f'repeat-{generation}' if repeat else f'job-{generation}'
    profile = 'LAYERED_SAME_CALL_V1' if mode=='benchmark' else 'INDEPENDENT_EACH_EVALUATION_V1'
    prefix = 'layered-formal-benchmark' if mode=='benchmark' else 'layered-independent-audit'
    job=HERE/f'{prefix}-{scale}-seed-{seed}-{suffix}'
    boundary='EXACT_RATIONAL_FIXED_SUBSET_V1' if generation in ('04','05','06') else 'EXACT_RATIONAL_BOUNDARY_V1'
    args=[scale, str(seed), 'NATIVE_GMP_SURFACE_V1', boundary,
          'EVALUATION_CHECKPOINT_V1', profile]
else: raise ValueError('UNKNOWN_LAYERED_JOB')
if generation in ('03','04','05','06'):
    arithmetic=HERE/'formal-rational-arithmetic-build-job-02/output/build/python/_cherchi_rational_moments.cpython-313-x86_64-linux-gnu.so'
    files['arithmetic']=arithmetic
    files['arithmetic_receipt']=HERE/'formal-rational-arithmetic-build-job-02/output/receipt.json'
    files['arithmetic_source']=ROOT/'bindings/cherchi/RationalMoments.cpp'
if generation in ('04','05'):
    native=HERE/'formal-native-core-build-job-06/output/build/python/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
    files['native']=native
    files['receipt']=HERE/'formal-native-core-build-job-06/output/receipt.json'
    files['boundary_cpp']=ROOT/'src/verification/ExactBoundaryAudit.cpp'
    files['boundary_header']=ROOT/'include/verification/ExactBoundaryAudit.h'
    files['bindings_cpp']=ROOT/'bindings/pybind_module.cpp'
if generation=='06':
    build=HERE/'formal-native-core-build-job-08'
    native=build/'output/build/python/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
    files.update(native=native,receipt=build/'output/receipt.json',
        candidate_comparison=HERE/'boundary-phase-measurement-job-02/output/report.json')
    if not json.loads(files['candidate_comparison'].read_text())['passed']:
        raise ValueError('DIAGONAL_BOUNDARY_COMPARISON_NOT_PASSED')
result=run_source_job(job,
    ['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
     'MAKE_IT_STAND_FORMAL_EVIDENCE_RUN='+str(HERE/'shared-vertex-formal-benchmark-small-seed-0-job-01/output/run'),
     '/home/lin/anaconda3/bin/python3','-B',str(files['worker']),'{job}/output',*args], files,
    environment={'PYTHONPATH':':'.join(map(str,[ROOT/'python',native.parent,bridge.parent,arithmetic.parent,HERE]))})
print(result.status,result.returncode,flush=True)
