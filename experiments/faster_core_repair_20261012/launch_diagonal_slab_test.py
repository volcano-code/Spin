from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
generation=sys.argv[1]
if generation not in ('01','02'):raise ValueError('UNREGISTERED_TEST_GENERATION')
files=dict(launcher=Path(__file__),worker=HERE/'build_diagonal_slab_test.py',
    test=ROOT/'tests/conservative_diagonal_slabs_regression.cpp',
    commands=HERE/'formal-native-core-build-job-07/output/build/compile_commands.json')
header=ROOT/'include/verification/ConservativeDiagonalSlabs.h'
if header.exists():files['header']=header
result=run_source_job(HERE/f'diagonal-slab-test-job-{generation}',
    ['/home/lin/anaconda3/bin/python3','-B',str(files['worker']),'{job}/output'],files)
print(result.status,result.returncode,flush=True)
