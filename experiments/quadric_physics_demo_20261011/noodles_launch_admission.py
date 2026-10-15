"""Launch exactly one resource-controlled geometry-only Noodles admission."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
from noodles_admission import publish,validate_manifest


def main(*,input_version='v1'):
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--job', type=Path, required=True)
    parser.add_argument('--native', type=Path, default=ROOT/'build/spin2014-objectives-quadric-20261007/python/_make_it_stand.cpython-313-x86_64-linux-gnu.so')
    parser.add_argument('--expected-native-sha256', required=True)
    args = parser.parse_args()
    case, native = args.input.resolve(), args.native.resolve()
    if hashlib.sha256(native.read_bytes()).hexdigest() != args.expected_native_sha256:
        raise ValueError('PINNED_NATIVE_IDENTITY_MISMATCH')
    manifest = json.loads((case/'manifest.json').read_text())
    validate_manifest(manifest,input_version)
    workers={'v1':'noodles_admission.py','dense-v2':'noodles_dense_admission.py','dense-v3':'noodles_dense_v3_admission.py','dense-v4':'noodles_dense_v4_admission.py','dense-v5':'noodles_dense_v5_admission.py'}
    generators={'v1':'noodles_input.py','dense-v2':'noodles_dense_input.py','dense-v3':'noodles_dense_v3_input.py','dense-v4':'noodles_dense_v4_input.py','dense-v5':'noodles_dense_v5_input.py'}
    worker=HERE/workers[input_version]
    generator=HERE/generators[input_version]
    heart = ROOT/'experiments/faster_heart_contact_20261011'
    files = dict(controller=Path(__file__), worker=worker,
        generator=generator, manifest=case/'manifest.json', native=native,
        exact_construction_helper=heart/'admit_rounded_input.py', independent_moments_helper=heart/'rounded_input.py',
        independent_heart_moments=ROOT/'experiments/reconstructed_heart_top/heart_top.py',
        problem=ROOT/'python/make_it_stand/verification/spin_problem.py',
        common=ROOT/'python/make_it_stand/verification/spin_common.py',
        libstdcpp=Path('/usr/lib/x86_64-linux-gnu/libstdc++.so.6').resolve(),
        python=Path('/home/lin/anaconda3/bin/python3').resolve())
    if input_version=='dense-v2':
        files.update(dense_launch_adapter=HERE/'noodles_dense_launch_admission.py',
            shared_admission=HERE/'noodles_admission.py',tube_mesh_helper=HERE/'noodles_input.py')
    elif input_version=='dense-v3':
        files.update(dense_launch_adapter=HERE/'noodles_dense_v3_launch_admission.py',
            shared_admission=HERE/'noodles_admission.py',tube_mesh_helper=HERE/'noodles_input.py')
    elif input_version in ('dense-v4','dense-v5'):
        files.update(dense_launch_adapter=HERE/('noodles_dense_v4_launch_admission.py' if input_version=='dense-v4' else 'noodles_dense_v5_launch_admission.py'),
            shared_admission=HERE/'noodles_admission.py',tube_mesh_helper=HERE/'noodles_input.py',
            inherited_geometry_generator=HERE/'noodles_dense_v3_input.py')
    for name, entry in manifest['parts'].items():
        for suffix in ('model', 'si'): files[name+'-'+suffix] = case/entry[suffix+'_obj']
    result = run_source_job(args.job.resolve(), ['/usr/bin/env',
        'LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
        '/home/lin/anaconda3/bin/python3', '-B', str(worker),
        '--input', str(case), '--output', '{job}/output'], files,
        maximum_wall_ns=3600_000_000_000, maximum_memory_bytes=24*1024**3,
        environment={'PYTHONPATH':str(ROOT/'python')+':'+str(native.parent)})
    publish(result.directory/'noodles-controller-terminal.json', dict(status=result.status,
        returncode=result.returncode, unit=result.unit, job=str(result.directory),
        scope='SOURCE_GEOMETRY_ADMISSION_ONLY', optimization_invocations=0, simulation_invocations=0))
    print(json.dumps(dict(status=result.status, returncode=result.returncode)), flush=True)
    return 0 if result.status=='PROCESS_EXITED' and result.returncode==0 else 1


if __name__ == '__main__': raise SystemExit(main())
