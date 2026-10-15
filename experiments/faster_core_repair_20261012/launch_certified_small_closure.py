"""Resource admission for the certified integral formal closed loop."""
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.serial_source_job import run_serial_source_job

if __name__=='__main__':
    build=HERE/'formal-native-core-build-job-12/output/build'
    bridge=HERE/'formal-cherchi-build-job-02/output/build/python'
    arithmetic=HERE/'formal-rational-arithmetic-build-job-05/output/build/python'
    files=dict(launcher=Path(__file__),worker=HERE/'run_certified_small_closure.py',
        benchmark=HERE/'run_core_benchmark.py',native=build/'python/_make_it_stand.cpython-313-x86_64-linux-gnu.so',
        runtime=build/'compiler-runtime/libstdc++.so.6',spin=ROOT/'python/make_it_stand/spin.py')
    for directory in (ROOT/'python/make_it_stand/verification',bridge,arithmetic,HERE/'core-benchmark-inputs-job-01/output'):
        for path in directory.iterdir():
            if path.is_file():files[str(path.relative_to(ROOT))]=path
    result=run_serial_source_job(ROOT,HERE/('certified-integral-small-closure-repeat-01' if '--repeat' in sys.argv else
                                           'certified-integral-small-closure-job-02'),
        ['/usr/bin/env','LD_PRELOAD='+str(files['runtime']),'/home/lin/anaconda3/bin/python','-B',
         str(files['worker']),'{job}/output'],files,
        environment={'PYTHONPATH':':'.join(map(str,[ROOT/'python',build/'python',bridge,arithmetic,HERE]))})
    print(result.status,result.returncode,flush=True)
    if result.status!='PROCESS_EXITED' or result.returncode!=0:raise SystemExit(1)
