"""Fresh ordinary CMake core build, not an overlay of historical objects."""
import hashlib
import json
from pathlib import Path
import shlex
import subprocess
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(output,profile='base'):
    output.mkdir(parents=True,exist_ok=False)
    build=output/'build'
    configure=['/usr/bin/cmake','-S',str(ROOT),'-B',str(build),
        '-DCMAKE_BUILD_TYPE=Release','-DCMAKE_CXX_FLAGS_RELEASE=-O2 -DNDEBUG',
        '-DCMAKE_INTERPROCEDURAL_OPTIMIZATION=OFF','-DCMAKE_EXPORT_COMPILE_COMMANDS=ON',
        '-DPython_EXECUTABLE='+sys.executable,
        '-DPython3_EXECUTABLE='+sys.executable,
        '-DBoost_INCLUDE_DIR='+str(Path(sys.base_prefix)/'include'),
        '-DMAKE_IT_STAND_BUILD_VIEWER=OFF','-DMAKE_IT_STAND_BUILD_REACHABILITY=OFF',
        '-DMAKE_IT_STAND_BUILD_SPIN_IT_FASTER=ON','-DMAKE_IT_STAND_BUILD_SPIN_IT_2014=OFF',
        '-DMAKE_IT_STAND_BUILD_EXACT_GEOMETRY_CHECKS=ON','-DMAKE_IT_STAND_BUILD_PYTHON=ON',
        '-DMAKE_IT_STAND_FILTERED_BOUNDARY_KERNEL=ON','-DBUILD_TESTING=OFF']
    compile_command=['/usr/bin/cmake','--build',str(build),'--target','_make_it_stand','--parallel','2']
    if profile == 'sharedvertex-slabs-spin2014':
        configure.remove('-DMAKE_IT_STAND_BUILD_SPIN_IT_2014=OFF')
        configure.append('-DMAKE_IT_STAND_BUILD_SPIN_IT_2014=ON')
        profile = 'sharedvertex-slabs'
    if profile in ('sharedvertex','sharedvertex-slabs'):configure.append('-DMAKE_IT_STAND_SHARED_VERTEX_CONE_PREDICATE=ON')
    elif profile!='base':raise ValueError('UNREGISTERED_FORMAL_NATIVE_PROFILE')
    if profile=='sharedvertex-slabs':configure.append('-DMAKE_IT_STAND_DIAGONAL_BOUNDARY_SLABS=ON')
    with (output/'build.log').open('x') as log:
        subprocess.run(configure,cwd=ROOT,stdout=log,stderr=log,check=True)
        rows=json.loads((build/'compile_commands.json').read_text())
        sources={r['file']:sha(Path(r['file'])) for r in rows}
        audit=[r for r in rows if Path(r['file']).name=='ExactBoundaryAudit.cpp']
        if len(audit)!=1 or '-DMAKE_IT_STAND_FILTERED_BOUNDARY_KERNEL=1' not in audit[0]['command']:
            raise ValueError('FORMAL_AUDIT_FILTER_CONFIGURATION_MISSING')
        if profile in ('sharedvertex','sharedvertex-slabs') and '-DMAKE_IT_STAND_SHARED_VERTEX_CONE_PREDICATE=1' not in audit[0]['command']:
            raise ValueError('FORMAL_SHARED_VERTEX_CONFIGURATION_MISSING')
        if profile=='sharedvertex-slabs' and '-DMAKE_IT_STAND_DIAGONAL_BOUNDARY_SLABS=1' not in audit[0]['command']:
            raise ValueError('FORMAL_DIAGONAL_SLAB_CONFIGURATION_MISSING')
        # The standard canonical implementation is deliberately preserved.
        if any('MAKE_IT_STAND_IDENTICAL_TRANSPORT_OVERLAY' in r['command'] for r in rows):
            raise ValueError('UNPROVEN_OVERLAY_MACRO_IN_FORMAL_BUILD')
        subprocess.run(compile_command,cwd=ROOT,stdout=log,stderr=log,check=True)
    for path,digest in sources.items():
        if sha(Path(path))!=digest:
            raise ValueError('COMPILED_SOURCE_CHANGED_DURING_BUILD')
    modules=list((build/'python').glob('_make_it_stand*.so'))
    if len(modules)!=1:
        raise ValueError('ONE_FORMAL_NATIVE_MODULE_REQUIRED')
    dependencies={}
    objects={}
    for path in build.rglob('*.o.d'):
        text=path.read_text().replace('\\\n',' ')
        for token in shlex.split(text.split(':',1)[1]):
            dep=Path(token)
            if dep.is_file():dependencies[str(dep.resolve())]=sha(dep)
        obj=Path(str(path)[:-2])
        if obj.is_file():objects[str(obj)]=sha(obj)
    link_files={str(p):p.read_text() for p in build.rglob('link.txt')}
    loaded_script=('import _make_it_stand as n,json; '
        'print(json.dumps(dict(module=n.__file__, libraries=sorted(set('
        'line.split()[-1] for line in open("/proc/self/maps") '
        'if len(line.split())>=6 and line.split()[-1].startswith("/") and ".so" in line.split()[-1])))))')
    import os
    env=dict(os.environ,LD_PRELOAD='/usr/lib/x86_64-linux-gnu/libstdc++.so.6',PYTHONPATH=str(build/'python'))
    loaded=json.loads(subprocess.check_output([sys.executable,'-c',loaded_script],env=env,text=True))
    if Path(loaded['module']).resolve()!=modules[0].resolve():
        raise ValueError('FORMAL_NATIVE_LOAD_IDENTITY_MISMATCH')
    loaded['library_sha256']={p:sha(Path(p)) for p in loaded['libraries'] if Path(p).is_file()}
    (output/'receipt.json').write_text(json.dumps(dict(commands=[configure,compile_command],
        module=str(modules[0]),module_sha256=sha(modules[0]),compiled_sources=sources,
        actual_dependency_sha256=dependencies,object_sha256=objects,actual_links=link_files,
        actual_loaded=loaded,cmake_sha256=sha(ROOT/'CMakeLists.txt'),
        canonical_transport_profile='STANDARD_CANONICAL_IMPLEMENTATION_NO_OVERLAY',
        boundary_predicate_profile=profile,
        core_admission=False),sort_keys=True))


if __name__=='__main__':main(Path(sys.argv[1]),sys.argv[2] if len(sys.argv)>2 else 'base')
