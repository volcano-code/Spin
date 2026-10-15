"""Native-only termination-policy comparison on a frozen reconstructed input.

Copies every numeric setting from an existing attempt. The ONLY treatment is
termination_policy. Does not load an old result as a witness or grant a seal.
Run under prlimit/timeout; never reuse an output directory.
"""
import argparse
import json
from pathlib import Path
import resource
import subprocess
import time
import traceback
import numpy as np
import trimesh
import _make_it_stand as native
from teapot import ROOT, check_files, file_hashes, write, write_mesh


def execute(reference, output, policy):
    output.mkdir(exist_ok=False,parents=True)
    start=time.monotonic()
    cfg=json.loads(reference.read_text())
    inputs=Path(cfg['input'])
    manifest=json.loads((inputs/'manifest.json').read_text())
    solver={**cfg['solver'],'termination_policy':policy}
    # Record effective code, not reports (which would create circular hashes).
    paths=[Path(native.__file__),Path(__file__),Path(__file__).with_name('teapot.py'),
           reference,inputs/'manifest.json']
    for directory in ('include','src','bindings','python/make_it_stand'):
        paths.extend(p for p in (ROOT/directory).rglob('*') if p.is_file() and
                     p.suffix in ('.h','.hpp','.cpp','.py'))
    paths.extend([ROOT/'CMakeLists.txt'])
    build=Path(native.__file__).resolve().parent.parent
    for name in ('CMakeCache.txt','compile_commands.json','build.ninja','libmake_it_stand_core.a'):
        if (build/name).is_file():paths.append(build/name)
    files=file_hashes(paths);files.update(manifest['geometry_hashes'])
    write(output/'source.json',dict(files=files,
        commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        source_diff=subprocess.check_output(['git','diff','--','src/spin2014','include/spin2014',
            'bindings/SourceSpinMaterialBindings.cpp'],cwd=ROOT,text=True),
        scope='SOURCE_AND_LOADED_BINARY_OBSERVATION_NOT_FROZEN_RELEASE_ATTESTATION'))
    write(output/'config.json',dict(reference_config=str(reference.resolve()),
        reference_config_hash=file_hashes([reference]),source_settings=cfg['geometry'],
        solver_settings=solver,boundary_settings=cfg['boundary'],density=1.,
        minimum_outer_wall=cfg['minimum_wall'],asset_evidence='RECONSTRUCTED',
        original_paper_case=False,only_treatment='termination_policy',
        claim_scope='NATIVE_OPTIMIZATION_AND_BINARY_CHECK_NOT_SEALED_OR_PHYSICALLY_VALIDATED'))
    try:
        check_files(files)
        meshes=[trimesh.load(inputs/(name+'.obj'),force='mesh',process=False)
                for name in ('outer','design')]
        args=[a for m in meshes for a in (np.asarray(m.vertices),np.asarray(m.faces))]
        print('START',cfg['strategy'],policy,flush=True)
        result=native.solve_source_spin2014_material(*args,source_settings=cfg['geometry'],
            solver_settings=solver,boundary_settings=cfg['boundary'],density=1.,
            minimum_outer_wall=cfg['minimum_wall'])
        facts=result.diagnostics
        write(output/'native-run.json',facts)
        if result.has_geometry:
            for name in ('outer','material','cavity'):
                write_mesh(output/(name+'.obj'),trimesh.Trimesh(*getattr(result,name),process=False))
        check_files(files)
        s=facts['solve']
        write(output/'summary.json',dict(status=s['status'],passes=len(s['passes']),
            binary_check=s['binary_check'],terminal_relaxation_kkt_applies=s['terminal_relaxation_kkt_applies'],
            trace_hash=s['trace_hash'],output_hash=s['output_hash'],
            native_geometry_available=result.has_geometry,independent_python_seal=False,
            physical_validation='NOT_RUN',original_paper_case=False,source_unchanged=True,
            elapsed_seconds=time.monotonic()-start,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
        write(output/'outputs.json',file_hashes(p for p in output.iterdir() if p.is_file()))
        print('RESULT',s['status'],'passes',len(s['passes']),'binary',s['binary_check'],flush=True)
        return 0 if s['status']=='BINARY_MOMENT_VALIDATED' else 2
    except Exception as exc:
        write(output/'failure.json',dict(error=str(exc),traceback=traceback.format_exc(),
            elapsed_seconds=time.monotonic()-start))
        raise


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--reference-config',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--policy',choices=['PAPER_POST_TRANSITION_V1','UNCHANGED_PARTITION_KKT_V1'],
                        default='PAPER_POST_TRANSITION_V1')
    args=parser.parse_args()
    raise SystemExit(execute(args.reference_config,args.output,args.policy))
