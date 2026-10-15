"""Frozen explicit companion certificate job; no solver or simulation."""
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
    parser.add_argument('--failed-job',type=Path,required=True);parser.add_argument('--job',type=Path,required=True)
    parser.add_argument('--expected-native-sha256',required=True)
    args=parser.parse_args();config_path=args.failed_job.resolve()/'config.json'
    config=json.loads(config_path.read_text());identities=config['identity']
    native=Path(identities['native']['path'])
    if sha(native)!=args.expected_native_sha256:raise ValueError('NATIVE_IDENTITY_MISMATCH')
    files=dict(controller=Path(__file__),worker=HERE/'build_bounded_pin.py',
        algorithm=HERE/'bounded_pin_correspondence.py',pin_loader=HERE/'pin_role_witness.py',
        helpers=HERE/'admit_bounded_request.py',failed_config=config_path,
        request_view=HERE/'rational_request_view.py',moments=HERE/'rounded_input.py',
        heart_dependency=ROOT/'experiments/reconstructed_heart_top/heart_top.py',
        problem=ROOT/'python/make_it_stand/verification/spin_problem.py',
        common=ROOT/'python/make_it_stand/verification/spin_common.py',
        rational_text=ROOT/'python/make_it_stand/verification/rational_text.py',
        native=native,runtime=Path(sys.executable).resolve(),
        libstdcpp=Path('/usr/lib/x86_64-linux-gnu/libstdc++.so.6').resolve())
    for key in ('authority-outer','authority-fixed','candidate-fixed.rational.json','candidate-design.rational.json',
                'pin','source_binding','original_admission','candidate'):
        info=identities[key];path=Path(info['path'])
        if sha(path)!=info['sha256']:raise ValueError('STALE_PIN_COMPANION_SOURCE: '+key)
        files['input-'+key]=path
    import fractions,gmpy2,numpy
    for name,module in [('fractions',fractions),('gmpy2',gmpy2),('numpy',numpy)]:files['runtime-'+name]=Path(module.__file__)
    for index,path in enumerate(sorted(Path(gmpy2.__file__).parent.glob('*.so'))):files[f'gmpy2_binary_{index}']=path
    for index,path in enumerate(sorted({line.split()[-1] for line in Path('/proc/self/maps').read_text().splitlines()
        if '/' in line and any(name in line for name in ('libgmp','libmpfr','libmpc'))})):
        files[f'exact_arithmetic_library_{index}']=Path(path)
    argv=['/usr/bin/env','LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6','OMP_NUM_THREADS=1',
        'OPENBLAS_NUM_THREADS=1','MKL_NUM_THREADS=1',str(Path(sys.executable).resolve()),'-B',str(HERE/'build_bounded_pin.py'),
        '--failed-job',str(args.failed_job.resolve()),'--output','{job}/output',
        '--expected-native-sha256',args.expected_native_sha256]
    result=run_source_job(args.job.resolve(),argv,files,maximum_wall_ns=3600_000_000_000,
        maximum_memory_bytes=24*1024**3,
        environment={'PYTHONPATH':str(ROOT/'python')+':'+str(native.parent)})
    print(json.dumps(dict(status=result.status,returncode=result.returncode,job=str(result.directory),
        optimization_invocations=0,simulation_invocations=0)))
    return 0 if result.status=='PROCESS_EXITED' and result.returncode==0 else 1


if __name__=='__main__':raise SystemExit(main())
