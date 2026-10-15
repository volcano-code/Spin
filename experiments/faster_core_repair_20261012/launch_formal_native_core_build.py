from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job

files=dict(launcher=Path(__file__),worker=HERE/'build_formal_native_core.py',cmake=ROOT/'CMakeLists.txt')
files['boost_version']=Path(sys.base_prefix)/'include/boost/version.hpp'
for folder in ('src','include','bindings','cmake'):
    for path in (ROOT/folder).rglob('*'):
        if path.is_file() and path.suffix in ('.cpp','.h','.hpp','.cmake','.txt'):
            files['source_'+str(path.relative_to(ROOT))]=path
for path in ('libigl/CMakeLists.txt','libigl/cmake/LibiglDownloadExternal.cmake'):
    files[path]=ROOT/path
generation=sys.argv[1] if len(sys.argv)>1 else '01'
if generation not in ('01','02','03','04','05','06','07','08'):
    raise ValueError('UNREGISTERED_FORMAL_BUILD_GENERATION')
result=run_source_job(HERE/('formal-native-core-build-job-'+generation),
    ['/home/lin/anaconda3/bin/python3','-B',str(files['worker']),'{job}/output',
     'sharedvertex-slabs' if generation=='08' else
     ('sharedvertex' if generation in ('05','06','07') else 'base')],files)
print(result.status,result.returncode,flush=True)
