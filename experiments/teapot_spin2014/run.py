"""Three separately executed native Spin2014 objectives on one locked teapot.

    Use a new output directory for every attempt. Not level-9/paper-original.
"""
import argparse
from dataclasses import asdict
from pathlib import Path
import json
import resource
import subprocess
import time
import traceback
import numpy as np
import trimesh
from teapot import ROOT, STRATEGIES, strategy_settings, write_mesh, write, file_hashes, check_files


def execute(input_dir, output_dir, strategy, expected_native_module=None):
    from make_it_stand.verification import source_spin_material as b
    from make_it_stand.verification.spin_common import native
    from make_it_stand.verification.spin_problem import SpinArtifactProblem
    from make_it_stand.verification.spin_bundle import save_spin_artifact
    if expected_native_module is not None:
        from supervised import verify_native_module
        verify_native_module(expected_native_module, native().__file__)
    input_dir=input_dir.resolve(); out=output_dir.resolve(); out.mkdir(exist_ok=False)
    start=time.monotonic()
    manifest=json.loads((input_dir/'manifest.json').read_text())
    sources=file_hashes([Path(__file__),Path(__file__).with_name('teapot.py'),Path(native().__file__),Path(b.__file__),input_dir/'manifest.json'])
    sources.update(manifest['geometry_hashes']); sources.update(manifest['source_files'])
    write(out/'source.json',dict(files=sources,head=subprocess.check_output(['git','rev-parse','HEAD'],text=True,cwd=ROOT).strip(),
        scope='LOADED_BINARY_AND_DECLARED_FILES_NOT_FULL_BUILD_ATTESTATION'))
    try:
        check_files(sources)
        mesh=lambda name: trimesh.load(input_dir/(name+'.obj'),force='mesh',process=False)
        arrays=lambda m:(np.asarray(m.vertices),np.asarray(m.faces))
        outer,design,fixed=map(mesh,['outer','design','fixed'])
        problem=SpinArtifactProblem(original=arrays(outer),fixed_regions=[('fixed-shell-and-shaft',arrays(fixed),1.)],
            design_regions=[('body-interior',arrays(design),1.)],length_unit='25mm',mass_unit='18.75g',
            target={'mode':'SPINNING_TOP','axis':[0,0,1],'origin':[0,0,0],'support_points':[]},minimum_outer_wall=.015)
        if problem.input_digest!=manifest['problem_digest']: raise ValueError('INPUT_IDENTITY_MISMATCH')
        geometry=b.SourceGeometrySettings(root_origin=(-1.5,-1.5,-.05),root_side=3.,maximum_depth=6)
        solver=strategy_settings(strategy); boundary=b.SourceBoundarySettings(maximum_input_cuboids=12000)
        write(out/'config.json',dict(strategy=strategy,geometry=asdict(geometry),solver=asdict(solver),boundary=asdict(boundary),
            input=str(input_dir),asset_evidence='RECONSTRUCTED',original_paper_case=False,
            resolution_scope='SOURCE_DEPTH_6_OPTIMIZATION_MAX_5_NOT_PAPER_LEVEL_9',
            meters_per_model_unit=.025,density_kg_m3=1200.,minimum_wall=.015,
            physical_protocol=dict(duration=8.,tilt_degrees=2.,azimuths=8,
                hz=[10.,20.,35.,45.],dt=1/2000,refined_dt=1/4000,friction=.8)))
        print('START',strategy,solver.objective,flush=True)
        result=b.solve_and_seal_source_spin2014(problem,source_settings=geometry,solver_settings=solver,boundary_settings=boundary)
        write(out/'native-run.json',result.identity);write(out/'evidence.json',result.evidence)
        for name,arrays_ in [('outer',result.outer),('material',result.material),('cavity',result.cavity)]:
            if arrays_ is not None: write_mesh(out/(name+'.obj'),trimesh.Trimesh(*arrays_,process=False))
        if result.artifact is not None: save_spin_artifact(result.artifact,out/'sealed')
        check_files(sources)
        write(out/'summary.json',dict(status=result.status,sealed=result.artifact is not None,seal_error=result.seal_error,
            strategy=strategy,problem_digest=problem.input_digest,scientific_digest=result.scientific_digest,
            source_unchanged=True,original_paper_case=False,elapsed_seconds=time.monotonic()-start,
            max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
        print('RESULT',strategy,result.status,'SEALED',result.artifact is not None,result.seal_error,flush=True)
        return 0 if result.artifact is not None else 2
    except Exception as e:
        write(out/'failure.json',dict(error=str(e),traceback=traceback.format_exc(),elapsed_seconds=time.monotonic()-start))
        raise


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--strategy',choices=STRATEGIES,required=True)
    p.add_argument('--expected-native-module',type=Path)
    a=p.parse_args();raise SystemExit(execute(a.input,a.output,a.strategy,a.expected_native_module))
