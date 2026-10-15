"""Whole-job launch of the existing reconstructed teapot protocol.

Does not change science settings, admit an original paper asset, resume an
arbitrary native checkpoint, or equate process exit with a successful seal.
"""
import argparse
import hashlib
import importlib.machinery
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'python'))
from make_it_stand.verification.source_job import run_source_job


def verify_native_module(expected, loaded):
    if Path(expected).resolve(strict=True) != Path(loaded).resolve(strict=True):
        raise ValueError('LOADED_NATIVE_MISMATCH')


def prepare(input_dir,strategy,native_module,ld_preload=None):
    if strategy not in {'inertia','lowering','mass_lowering'}:raise ValueError('UNSUPPORTED_STRATEGY')
    inputs=Path(input_dir).resolve(strict=True);binary=Path(native_module).resolve(strict=True)
    if not binary.is_file() or binary.name not in {
        '_make_it_stand'+suffix for suffix in importlib.machinery.EXTENSION_SUFFIXES
    }:
        raise ValueError('EXPLICIT_NATIVE_MODULE_REQUIRED')
    manifest=inputs/'manifest.json';data=json.loads(manifest.read_text())
    geometry=data['geometry_hashes']
    if not {str(inputs/(name+'.obj')) for name in ('outer','design','fixed')}<=set(geometry):
        raise ValueError('INCOMPLETE_GEOMETRY_MANIFEST')
    declared={**data['source_files'],**geometry}
    for name,expected in declared.items():
        with Path(name).open('rb') as stream:actual=hashlib.file_digest(stream,'sha256').hexdigest()
        if actual!=expected:raise ValueError('STALE_TEAPOT_INPUT: '+name)
    runner=Path(__file__).resolve().with_name('run.py')
    identity={'manifest':manifest,'native':binary,'adapter':Path(__file__).resolve(),
              **{name:Path(name) for name in declared}}
    for folder in ['include','src','bindings','python/make_it_stand','experiments/teapot_spin2014','experiments/heart_stem_contact']:
        for path in (ROOT/folder).rglob('*'):
            if path.is_file() and path.suffix in {'.py','.h','.hpp','.cpp'}:identity[str(path)]=path
    identity['cmake']=ROOT/'CMakeLists.txt'
    argv=[sys.executable,str(runner),'--input',str(inputs),'--output','{job}/payload','--strategy',strategy,
          '--expected-native-module',str(binary)]
    if ld_preload is not None:
        library=Path(ld_preload).resolve(strict=True)
        if not library.is_file():raise ValueError('PRELOAD_FILE_REQUIRED')
        identity['explicit_preload']=library
        argv=['/usr/bin/env','LD_PRELOAD='+str(library),*argv]
    return dict(argv=argv,identity_files=identity,
        environment={'PYTHONPATH':str(ROOT/'python')+':'+str(binary.parent)},
        scope='RECONSTRUCTED_EXISTING_TEAPOT_PROTOCOL_NOT_PAPER_ORIGINAL')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--strategy',choices=['inertia','lowering','mass_lowering'],required=True)
    p.add_argument('--native-module',type=Path,required=True)
    p.add_argument('--ld-preload',type=Path)
    p.add_argument('--wall-seconds',type=int,default=3600)
    a=p.parse_args()
    if not 1<=a.wall_seconds<=3600:p.error('wall-seconds must be in [1,3600]')
    request=prepare(a.input,a.strategy,a.native_module,a.ld_preload)
    result=run_source_job(a.output,request['argv'],request['identity_files'],
        environment=request['environment'],maximum_wall_ns=a.wall_seconds*1_000_000_000)
    print(json.dumps(dict(status=result.status,returncode=result.returncode,scope=request['scope'],
                         scientific_validation='READ_EXISTING_RUNNER_SEAL_AND_INDEPENDENT_EVIDENCE')))
    raise SystemExit(0 if result.status=='PROCESS_EXITED' else 2)
