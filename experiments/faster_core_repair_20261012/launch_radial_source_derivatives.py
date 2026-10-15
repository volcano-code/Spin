"""Run actual radial source validation after independent build17."""
from pathlib import Path
import json
import sys
import time
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.serial_source_job import run_serial_source_job


def main():
    job=HERE/'formal-native-core-build-job-17'
    while not (job/'result.json').exists():time.sleep(10)
    result=json.loads((job/'result.json').read_text())
    if result['status']!='PROCESS_EXITED' or result['returncode']!=0:
        raise RuntimeError('RADIAL_SOURCE_BUILD_FAILED')
    build=job/'output/build';run=HERE/'heart-certified-formal-seed-1-job-01/output/run'
    files=dict(launcher=Path(__file__),worker=HERE/'probe_radial_source_derivatives.py',
        native=build/'python/_make_it_stand.cpython-313-x86_64-linux-gnu.so',
        runtime=build/'compiler-runtime/libstdc++.so.6',receipt=job/'output/receipt.json',
        configuration=run/'configuration.json',
        registration=HERE/'heart-seed1-source-failure-registration-01.json')
    while True:
        try:
            result=run_serial_source_job(ROOT,HERE/'heart-radial-source-derivatives-job-01',
                ['/usr/bin/env','LD_PRELOAD='+str(files['runtime']),'/home/lin/anaconda3/bin/python',
                 '-B',str(files['worker']),'{job}/output'],files,
                environment={'PYTHONPATH':str(build/'python')})
            print(result.status,result.returncode,flush=True);return
        except RuntimeError as error:
            if not str(error).startswith(('PROJECT_JOB_BUSY','PROJECT_WORKLOAD_ACTIVE:')):raise
            time.sleep(10)


if __name__=='__main__':main()
