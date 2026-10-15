"""Registered Heart seeds and repeats; project serial admission is mandatory."""
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT/'python'))
from make_it_stand.verification.serial_source_job import run_serial_source_job


def main():
    seed = int(sys.argv[1])
    repeat = len(sys.argv) == 3 and sys.argv[2] == 'repeat'
    if seed not in (0, 1) or len(sys.argv) not in (2, 3) or (len(sys.argv) == 3 and not repeat):
        raise ValueError('UNREGISTERED_HEART_RUN')
    admission_path = HERE/'heart-guarded-solve-admission-01/report.json'
    admission = json.loads(admission_path.read_text())
    if not admission['guarded_solve_admission']:
        raise ValueError('HEART_NOT_ADMITTED_FOR_GUARDED_ATTEMPT')
    files = dict(launcher=Path(__file__), worker=HERE/'run_heart_guarded_formal_case.py',
        admission=admission_path, runtime_result=admission_path.parent/'runtime-test-result.json',
        coordinate_contract=HERE/'case_coordinate_contract.py', solver_protocol=HERE/'run_noodles_formal_case.py',
        coordinate_registration=HERE/'heart-cherchi-coordinate-map-job-01/output/registration.json')
    for key, record in admission['bindings'].items():
        path = Path(record['path'])
        if hashlib.sha256(path.read_bytes()).hexdigest() != record['sha256']:
            raise ValueError('HEART_ADMISSION_DEPENDENCY_CHANGED:'+key)
        files['admission_'+key] = path
    native = Path(admission['bindings']['native']['path'])
    bridge = HERE/'formal-cherchi-build-job-02/output/build/python/_cherchi_core_boolean.cpython-313-x86_64-linux-gnu.so'
    arithmetic = HERE/'formal-rational-arithmetic-build-job-02/output/build/python/_cherchi_rational_moments.cpython-313-x86_64-linux-gnu.so'
    prepared_config = json.loads(Path(admission['bindings']['preparation_resource_config']['path']).read_text())
    for path in (native, bridge, arithmetic):
        references = [record for record in prepared_config['identity'].values()
                      if Path(record['path']).resolve() == path.resolve()]
        if not references or any(hashlib.sha256(path.read_bytes()).hexdigest() != record['sha256']
                                 for record in references):
            raise ValueError('HEART_PREPARED_NATIVE_BACKEND_CHANGED:'+str(path))
    files.update(bridge=bridge, arithmetic=arithmetic)
    for path in (ROOT/'python/make_it_stand/verification').glob('*.py'):
        files['module_'+path.name] = path
    suffix = 'repeat-01' if repeat else 'job-01'
    result = run_serial_source_job(ROOT, HERE/f'heart-guarded-formal-seed-{seed}-{suffix}',
        ['/usr/bin/env', 'LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
         '/home/lin/anaconda3/bin/python3', '-B', str(files['worker']), '{job}/output', str(seed)], files,
        environment={'PYTHONPATH':':'.join(map(str, [ROOT/'python', native.parent, bridge.parent, arithmetic.parent, HERE]))})
    print(result.status, result.returncode, flush=True)
    if result.status != 'PROCESS_EXITED' or result.returncode != 0:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
