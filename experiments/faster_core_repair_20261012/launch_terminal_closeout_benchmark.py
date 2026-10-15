"""Serial seeded benchmark for durable terminals and one-time preparation."""
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT/'python'))
from make_it_stand.verification.serial_source_job import run_serial_source_job


def main():
    args = list(sys.argv[1:])
    prepared_arithmetic = '--prepared-arithmetic' in args
    if prepared_arithmetic:
        args.remove('--prepared-arithmetic')
    separated = '--separated-roots' in args
    if separated:
        args.remove('--separated-roots')
    classified = '--classified-roots' in args
    if classified:
        args.remove('--classified-roots')
    if separated and classified:
        raise ValueError('SELECT_ONE_BOUNDARY_CLASSIFICATION_PROFILE')
    reuse_input = '--reuse-input-proof' in args
    if reuse_input:
        args.remove('--reuse-input-proof')
    owned_terminal = '--owned-terminal' in args
    if owned_terminal:
        args.remove('--owned-terminal')
    if len(args) not in (3, 4):
        raise ValueError('UNREGISTERED_CLOSEOUT_BENCHMARK')
    scale, seed, generation = args[0], int(args[1]), args[2]
    repeat = len(args) == 4 and args[3] == 'repeat'
    if ((len(args) == 4 and not repeat) or
            scale not in ('small', 'faces-12288', 'faces-98304') or seed not in range(10) or
            not generation.isdigit() or len(generation) != 2):
        raise ValueError('UNREGISTERED_CLOSEOUT_BENCHMARK')
    native_generation = '12' if owned_terminal or classified else '10' if separated else '08'
    native = HERE/f'formal-native-core-build-job-{native_generation}/output/build/python/_make_it_stand.cpython-313-x86_64-linux-gnu.so'
    bridge = HERE/'formal-cherchi-build-job-02/output/build/python/_cherchi_core_boolean.cpython-313-x86_64-linux-gnu.so'
    arithmetic_generation = '03' if prepared_arithmetic else '02'
    arithmetic = HERE/f'formal-rational-arithmetic-build-job-{arithmetic_generation}/output/build/python/_cherchi_rational_moments.cpython-313-x86_64-linux-gnu.so'
    runtime = native.parent.parent/'compiler-runtime/libstdc++.so.6'
    if not runtime.is_file():
        raise ValueError('MATCHING_NATIVE_COMPILER_RUNTIME_MISSING')
    files = dict(launcher=Path(__file__), worker=HERE/'run_core_benchmark.py',
        native=native, runtime=runtime, bridge=bridge, arithmetic=arithmetic, spin=ROOT/'python/make_it_stand/spin.py',
        configuration=HERE/'core-benchmark-inputs-job-01/output'/(scale+'.json'),
        input_manifest=HERE/'core-benchmark-inputs-job-01/output/manifest.json')
    for path in (ROOT/'python/make_it_stand/verification').glob('*.py'):
        files['module_'+path.name] = path
    suffix = ('repeat-' if repeat else 'job-')+generation
    prefix = 'classified-closeout' if classified else 'separated-closeout' if separated else 'terminal-closeout'
    if reuse_input:
        prefix = 'prepared-'+prefix
    if owned_terminal:
        prefix = 'owned-'+prefix
    if prepared_arithmetic:
        prefix = 'arithmetic-prepared-'+prefix
    job = HERE/f'{prefix}-{scale}-seed-{seed}-{suffix}'
    if separated or owned_terminal:
        files['native_receipt'] = HERE/f'formal-native-core-build-job-{native_generation}/output/receipt.json'
    result = run_serial_source_job(ROOT, job,
        ['/usr/bin/env', 'LD_PRELOAD='+str(runtime),
         '/home/lin/anaconda3/bin/python3', '-B', str(files['worker']), '{job}/output', scale, str(seed),
         'NATIVE_GMP_SURFACE_V1',
         'EXACT_RATIONAL_FIXED_CLASSIFIED_V1' if classified else
         'EXACT_RATIONAL_FIXED_SEPARATED_ROOTS_V1' if separated else 'EXACT_RATIONAL_FIXED_SUBSET_V1',
         'EVALUATION_CHECKPOINT_V1', 'LAYERED_SAME_CALL_V1',
         'CONSTRUCTOR_PROOF_REUSE_V1' if reuse_input else 'RECHECK_INPUT_V1',
         'OWNED_INTERIOR_DIFFERENCE_V1' if owned_terminal else 'INDEXED_BOUNDARY_REBUILD_V1'], files,
        environment={'PYTHONPATH':':'.join(map(str, [ROOT/'python', native.parent, bridge.parent, arithmetic.parent]))})
    print(result.status, result.returncode, flush=True)
    if result.status != 'PROCESS_EXITED' or result.returncode != 0:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
