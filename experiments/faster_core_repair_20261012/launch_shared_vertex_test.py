from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job

generation=sys.argv[1] if len(sys.argv)>1 else '01'
if generation not in ('01','02','03'):raise ValueError('UNREGISTERED_VERTEX_UNIT_GENERATION')
files=dict(launcher=Path(__file__),worker=HERE/'build_shared_vertex_test.py',
    test=ROOT/'tests/exact_shared_vertex_predicate_regression.cpp',
    commands=HERE/'formal-native-core-build-job-04/output/build/compile_commands.json')
header=ROOT/'include/verification/ExactSharedVertexIntersection.h'
if header.exists():files['header']=header
result=run_source_job(HERE/('shared-vertex-predicate-test-job-'+generation),
    ['/home/lin/anaconda3/bin/python3','-B',str(files['worker']),'{job}/output',
     'filtered' if generation=='03' else 'unfiltered'],files)
print(result.status,result.returncode,flush=True)
