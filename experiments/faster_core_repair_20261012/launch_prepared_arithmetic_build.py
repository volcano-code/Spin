"""Independent serial build; never replaces the frozen arithmetic02 module."""
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.serial_source_job import run_serial_source_job


def main():
    files = dict(launcher=Path(__file__),worker=HERE/'build_formal_rational_arithmetic.py',
        source=ROOT/'bindings/cherchi/RationalMoments.cpp',cmake=ROOT/'bindings/cherchi/CMakeLists.txt',
        boolean_source=ROOT/'bindings/cherchi/InProcessBoolean.cpp')
    result = run_serial_source_job(ROOT,HERE/'formal-rational-arithmetic-build-job-03',
        ['/home/lin/anaconda3/bin/python3','-B',str(files['worker']),'{job}/output'],files)
    print(result.status,result.returncode,flush=True)
    if result.status != 'PROCESS_EXITED' or result.returncode != 0:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
