"""Serially scheduled native geometry gate; no optimization or simulation."""
import argparse
import json
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
from rational_request_view import sha


def main():
    parser=argparse.ArgumentParser()
    for name in ('admission','candidate','source-binding','pin','job','native'):
        parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--expected-native-sha256',required=True)
    parser.add_argument('--pin-policy',choices=('ORIGINAL_PIN_STRICT_V1','BOUNDED_PIN_INTERSECTION_V1'),default='ORIGINAL_PIN_STRICT_V1')
    parser.add_argument('--pin-companion',type=Path)
    args=parser.parse_args()
    paths={name:getattr(args,name).resolve() for name in ('admission','candidate','source_binding','pin','job','native')}
    if sha(paths['native'])!=args.expected_native_sha256:raise ValueError('NATIVE_IDENTITY_MISMATCH')
    binding=json.loads(paths['source_binding'].read_text())
    if (binding.get('source_identity_verified')is not True or
        binding.get('admission_sha256')!=sha(paths['admission']/'admission.json') or
        binding['frozen_source']['pin.obj']!=sha(paths['pin'])):
        raise ValueError('SOURCE_BINDING_MISMATCH')
    files=dict(controller=Path(__file__),worker=HERE/'admit_bounded_request.py',
        request_view=HERE/'rational_request_view.py',
        independent_moments=HERE/'rounded_input.py',
        independent_heart_dependency=ROOT/'experiments/reconstructed_heart_top/heart_top.py',
        problem=ROOT/'python/make_it_stand/verification/spin_problem.py',
        common=ROOT/'python/make_it_stand/verification/spin_common.py',
        rational_text=ROOT/'python/make_it_stand/verification/rational_text.py',
        source_binding_checker=HERE/'clearance_source_binding.py',
        source_binding=paths['source_binding'],pin=paths['pin'],native=paths['native'],
        original_admission=paths['admission']/'admission.json',candidate=paths['candidate']/'candidate.json',
        runtime=Path(sys.executable).resolve(),libstdcpp=Path('/usr/lib/x86_64-linux-gnu/libstdc++.so.6').resolve())
    candidate=json.loads((paths['candidate']/'candidate.json').read_text())
    for name,want in candidate['files'].items():
        path=paths['candidate']/name
        if sha(path)!=want:raise ValueError('STALE_CANDIDATE: '+name)
        files['candidate-'+name]=path
    record=json.loads((paths['admission']/'admission.json').read_text())
    for name in ('outer','fixed','design'):
        path=paths['admission']/(name+'.rational.json')
        if sha(path)!=record['files'][path.name]:raise ValueError('STALE_AUTHORITY: '+name)
        files['authority-'+name]=path
    if args.pin_policy=='ORIGINAL_PIN_STRICT_V1' and args.pin_companion is not None:
        raise ValueError('EXPLICIT_PIN_POLICY_REQUIRED')
    if args.pin_policy=='BOUNDED_PIN_INTERSECTION_V1':
        if args.pin_companion is None:raise ValueError('PIN_COMPANION_REQUIRED')
        companion=args.pin_companion.resolve()
        evidence=json.loads(companion.read_text())
        if evidence.get('schema')!='bounded-pin-companion-admission-v1' or \
                evidence.get('status')!='PIN_CORRESPONDENCE_CERTIFIED':
            raise ValueError('PIN_COMPANION_NOT_CERTIFIED')
        files['bounded_pin_companion']=companion
        for name,want in evidence.get('files',{}).items():
            path=companion.parent/name
            if path.parent.resolve()!=companion.parent.resolve() or path.is_symlink() or sha(path)!=want:
                raise ValueError('PIN_COMPANION_FILE_IDENTITY_MISMATCH')
            files['bounded_pin_attachment-'+name]=path
        for name in ('bounded_pin_attachment.py','bounded_pin_correspondence.py','pin_role_witness.py'):
            files['bounded-pin-checker-'+name]=HERE/name
    import fractions,gmpy2
    files['fractions_runtime']=Path(fractions.__file__)
    files['gmpy2_package']=Path(gmpy2.__file__)
    for index,path in enumerate(sorted(Path(gmpy2.__file__).parent.glob('*.so'))):files[f'gmpy2_binary_{index}']=path
    libraries={line.split()[-1] for line in Path('/proc/self/maps').read_text().splitlines()
               if '/' in line and any(name in line for name in ('libgmp','libmpfr','libmpc'))}
    for index,path in enumerate(sorted(libraries)):files[f'exact_arithmetic_library_{index}']=Path(path)
    argv=['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6',
        'OMP_NUM_THREADS=1','OPENBLAS_NUM_THREADS=1','MKL_NUM_THREADS=1',str(Path(sys.executable).resolve()),
        '-B',str(HERE/'admit_bounded_request.py')]
    for name in ('admission','candidate','source_binding','pin'):
        argv+=['--'+name.replace('_','-'),str(paths[name])]
    argv+=['--output','{job}/output','--expected-native-sha256',args.expected_native_sha256]
    argv+=['--pin-policy',args.pin_policy]
    if args.pin_companion is not None:argv+=['--pin-companion',str(args.pin_companion.resolve())]
    result=run_source_job(paths['job'],argv,files,maximum_wall_ns=3600_000_000_000,
        maximum_memory_bytes=24*1024**3,
        environment={'PYTHONPATH':str(ROOT/'python')+':'+str(paths['native'].parent)})
    print(json.dumps(dict(status=result.status,returncode=result.returncode,
        job=str(result.directory),unit=result.unit,scope='DERIVED_REQUEST_GEOMETRY_ONLY',
        optimization_invocations=0,simulation_invocations=0)))
    return 0 if result.status=='PROCESS_EXITED' and result.returncode==0 else 1


if __name__=='__main__':raise SystemExit(main())
