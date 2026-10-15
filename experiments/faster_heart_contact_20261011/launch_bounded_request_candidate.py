"""One bounded Python-only candidate job; explicitly invoked, no native call."""
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
    parser.add_argument('--admission',type=Path,required=True)
    parser.add_argument('--source-binding',type=Path,required=True)
    parser.add_argument('--job',type=Path,required=True)
    args=parser.parse_args()
    admission=args.admission.resolve(); binding=args.source_binding.resolve()
    b=json.loads(binding.read_text())
    if b.get('source_identity_verified') is not True or b.get('admission_sha256')!=sha(admission/'admission.json'):
        raise ValueError('SOURCE_BINDING_NOT_FOR_ACTUAL_ADMISSION')
    if b.get('status')!='INPUT_ADMITTED_FEASIBILITY_UNKNOWN':
        raise ValueError('INPUT_NOT_ADMITTED_UNKNOWN_FEASIBILITY')
    import fractions,json as json_module,gmpy2
    files=dict(controller=Path(__file__),worker=HERE/'bounded_request_candidate.py',
        request_view_checker=HERE/'rational_request_view.py',
        independent_moments=HERE/'rounded_input.py',
        heart_generator_dependency=ROOT/'experiments/reconstructed_heart_top/heart_top.py',
        admission=admission/'admission.json',source_binding=binding,
        source_binding_checker=HERE/'clearance_source_binding.py',
        runtime=Path(sys.executable).resolve(),fractions=Path(fractions.__file__),
        json_runtime=Path(json_module.__file__),gmpy2_package=Path(gmpy2.__file__))
    for index,path in enumerate(sorted(Path(gmpy2.__file__).parent.glob('*.so'))):
        files[f'gmpy2_binary_{index}']=path
    loaded_libraries={line.split()[-1] for line in Path('/proc/self/maps').read_text().splitlines()
                      if '/' in line and any(name in line for name in ('libgmp','libmpfr','libmpc'))}
    for index,path in enumerate(sorted(loaded_libraries)):
        files[f'exact_arithmetic_library_{index}']=Path(path)
    for name in ('outer','fixed','design'):
        files[name]=admission/(name+'.rational.json')
    result=run_source_job(args.job.resolve(),['/usr/bin/env','OMP_NUM_THREADS=1',
        'OPENBLAS_NUM_THREADS=1','MKL_NUM_THREADS=1',str(Path(sys.executable).resolve()),
        '-B',str(HERE/'bounded_request_candidate.py'),'--admission',str(admission),
        '--output','{job}/output'],files,maximum_wall_ns=3600_000_000_000,
        maximum_memory_bytes=24*1024**3,environment={'PYTHONPATH':str(ROOT/'python')})
    print(json.dumps(dict(status=result.status,returncode=result.returncode,
        scope='PYTHON_ONLY_DERIVED_CANDIDATE_NOT_NATIVE_CERTIFICATION_OR_OPTIMIZATION',
        job=str(result.directory),unit=result.unit)))
    return 0 if result.status=='PROCESS_EXITED' and result.returncode==0 else 1


if __name__=='__main__':raise SystemExit(main())
