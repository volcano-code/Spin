from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / 'python'))
from make_it_stand.verification.source_job import run_source_job

scale, seed = sys.argv[1], int(sys.argv[2])
flags=sys.argv[3:]
if len(set(flags))!=len(flags) or any(flag not in ('repeat','gmp','boundary','checkpoint') for flag in flags):
    raise ValueError('UNREGISTERED_BENCHMARK_FLAGS')
suffix = 'repeat-01' if 'repeat' in flags else 'job-01'
arithmetic_backend='NATIVE_GMP_SURFACE_V1' if 'gmp' in flags else 'PYTHON_FRACTION_V1'
if 'gmp' in flags:suffix='gmp-'+suffix
boundary_profile='EXACT_RATIONAL_BOUNDARY_V1' if 'boundary' in flags else 'EXACT_NEF_REGION_V1'
if 'boundary' in flags:suffix='boundary-'+suffix
durability_profile='EVALUATION_CHECKPOINT_V1' if 'checkpoint' in flags else 'PER_EVENT_FSYNC_V1'
if 'checkpoint' in flags:suffix='checkpoint-'+suffix
if scale not in ('small', 'faces-12288', 'faces-98304') or seed not in range(10):
    raise ValueError('UNREGISTERED_SCALE_OR_SEED')
native = HERE / 'single-transport-build-job-01/output/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
native_receipt=HERE/'single-transport-build-job-01/output/receipt.json'
if 'boundary' in flags:
    native=HERE/'rational-boundary-binding-build-job-01/output/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
    native_receipt=HERE/'rational-boundary-binding-build-job-01/output/receipt.json'
bridge = HERE / 'formal-cherchi-build-job-02/output/build/python/_cherchi_core_boolean.cpython-313-x86_64-linux-gnu.so'
files = dict(launcher=Path(__file__), worker=HERE / 'run_core_benchmark.py', native=native, bridge=bridge,
    native_receipt=native_receipt,
    config=HERE / 'core-benchmark-inputs-job-01/output' / (scale + '.json'),
    manifest=HERE / 'core-benchmark-inputs-job-01/output/manifest.json',
    spin=ROOT / 'python/make_it_stand/spin.py')
for name in ('cherchi_inprocess', 'cherchi_run_evidence', 'cherchi_run_journal', 'cherchi_motion_differential', 'faster_material',
             'spin_problem', 'spin_bundle', 'spin_artifact', 'spin_common', 'mass_moments', 'rational_text'):
    files[name] = ROOT / f'python/make_it_stand/verification/{name}.py'
engine=None
if 'gmp' in flags:
    engine=HERE/'formal-rational-arithmetic-build-job-01/output/build/python/_cherchi_rational_moments.cpython-313-x86_64-linux-gnu.so'
    files['arithmetic_module']=engine
    files['arithmetic_receipt']=HERE/'formal-rational-arithmetic-build-job-01/output/receipt.json'
result = run_source_job(HERE / f'core-benchmark-{scale}-seed-{seed}-{suffix}',
    ['/usr/bin/env', 'LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
     '/home/lin/anaconda3/bin/python3', '-B', str(files['worker']), '{job}/output', scale, str(seed),arithmetic_backend,boundary_profile,durability_profile], files,
    environment={'PYTHONPATH': str(ROOT / 'python') + ':' + str(native.parent) + ':' + str(bridge.parent)
        +(':'+str(engine.parent) if engine is not None else '')})
print(result.status, result.returncode, flush=True)
