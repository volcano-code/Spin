"""One genuine original-objective seeded Cherchi run and independent seal."""
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT/'python'))
from make_it_stand.verification.serial_source_job import run_serial_source_job


def main():
    selector = sys.argv[1] if len(sys.argv)>1 else 'TOP_LOWERING_ONLY'
    label = sys.argv[2] if len(sys.argv)>2 else 'repeat-01'
    if selector not in ('TOP_INERTIA_ONLY','TOP_LOWERING_ONLY','TOP_MASS_AND_LOWERING'):
        raise ValueError('unsupported spinning-top smoke objective')
    if label not in ('run-01','repeat-01'): raise ValueError('unregistered smoke run')
    build = HERE/'formal-native-core-build-job-13/output/build'
    bridge = HERE/'formal-cherchi-build-job-02/output/build/python'
    arithmetic = HERE/'formal-rational-arithmetic-build-job-03/output/build/python'
    runtime = build/'compiler-runtime/libstdc++.so.6'
    files = dict(worker=HERE/'run_core_benchmark.py', launcher=Path(__file__),
        native=build/'python/_make_it_stand.cpython-313-x86_64-linux-gnu.so',
        receipt=build.parent/'receipt.json', runtime=runtime, spin=ROOT/'python/make_it_stand/spin.py')
    for directory in (ROOT/'python/make_it_stand/verification', HERE/'core-benchmark-inputs-job-01/output'):
        for path in directory.glob('*'):
            if path.is_file(): files[str(path.relative_to(ROOT))] = path
    for directory in (bridge, arithmetic):
        for path in directory.glob('*.so'): files[path.name] = path
    name = 'spin2014-shared-small-'+selector.lower().replace('_','-')+'-seed-1-'+label
    result = run_serial_source_job(ROOT, HERE/name,
        ['/usr/bin/env', 'LD_PRELOAD='+str(runtime), '/home/lin/anaconda3/bin/python', '-B',
         str(files['worker']), '{job}/output', 'small', '1', 'NATIVE_GMP_SURFACE_V1',
         'EXACT_RATIONAL_FIXED_CLASSIFIED_V1', 'EVALUATION_CHECKPOINT_V1', 'LAYERED_SAME_CALL_V1',
         'CONSTRUCTOR_PROOF_REUSE_V1', 'OWNED_INTERIOR_DIFFERENCE_V1', selector], files,
        environment={'PYTHONPATH': ':'.join(map(str, [ROOT/'python', build/'python', bridge, arithmetic, HERE]))})
    print(result.status, result.returncode, flush=True)
    if result.status != 'PROCESS_EXITED' or result.returncode != 0: raise SystemExit(1)


if __name__ == '__main__': main()
