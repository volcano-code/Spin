from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job
source=ROOT/'external/reference-sources-20260917/sources/InteractiveAndRobustMeshBooleans-7bd6c2697695f555f9a15d7e154e084d841316f0'
files=dict(launcher=Path(__file__),worker=HERE/'build_formal_cherchi.py',
           cmake=ROOT/'bindings/cherchi/CMakeLists.txt',bridge=ROOT/'bindings/cherchi/InProcessBoolean.cpp')
for index,path in enumerate(sorted(source.rglob('*'))):
    if path.is_file() and (path.suffix in ('.h','.hpp','.cpp','.c','.cmake') or path.name=='CMakeLists.txt'):
        files[f'dependency_{index}']=path
result=run_source_job(HERE/'formal-cherchi-build-job-02',
    ['/home/lin/anaconda3/bin/python3','-B',str(files['worker']),'{job}/output'],files)
print(result.status,result.returncode,flush=True)
