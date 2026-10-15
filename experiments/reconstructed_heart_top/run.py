"""Run a fixed reconstructed heart through the existing source2014 solver.
Use a new output directory; never overwrite an attempt. No paper-case credit.
"""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import subprocess
import time
import traceback
import numpy as np
import trimesh
from heart_top import heart_mesh, integrate


def write(path,data):
    path.write_text(json.dumps(data,sort_keys=True,indent=2,allow_nan=False)+'\n')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args(); out=args.output; out.mkdir(parents=True,exist_ok=False)
    start=time.monotonic()
    try:
        from make_it_stand.verification import source_spin_material as b
        from make_it_stand.verification.spin_problem import SpinArtifactProblem
        from make_it_stand.verification.spin_bundle import save_spin_artifact
        from make_it_stand.verification.spin_common import native
        v,f=heart_mesh(); center=np.array([0.,0.,.95])
        dv=center+.72*(v-center)
        shell=(np.vstack((v,dv)),np.vstack((f,f[:,::-1]+len(v))))
        for name,mesh in [('outer',(v,f)),('design',(dv,f)),('fixed-shell',shell)]:
            trimesh.Trimesh(*mesh,process=False).export(out/(name+'.obj'))
        geometry=b.SourceGeometrySettings(root_origin=(-1.5,-1.5,-.25),root_side=3.,maximum_depth=5)
        solver=b.SourceSolverSettings(objective='TOP_MASS_AND_LOWERING',initial_level=3,
            maximum_optimization_level=4,maximum_adaptive_passes=10,maximum_leaf_count=4000,
            initial_fill=.5,initial_phi=.2,laplacian_weight=.001,
            slqp=b.AuditedSlqpSettings(maximum_iterations=160))
        boundary=b.SourceBoundarySettings(maximum_input_cuboids=4000)
        config=dict(asset='RECONSTRUCTED_SELF_BUILT_HEART_V1',author_original=False,
            profile_samples=48,profile_layers=12,design_scale=.72,design_center=center.tolist(),
            source=asdict(geometry),solver=asdict(solver),boundary=asdict(boundary),
            minimum_outer_wall=.015,solver_density=1.,meters_per_model_unit=.025,
            physical_density_kg_per_m3=1200.,physical_protocol='IDEAL_FIXED_PIVOT_GRAVITY',
            rotation_hz=[2.,6.,15.],tilt_degrees=2.,azimuth_count=8,duration_seconds=8.,
            max_step_seconds=1/240,refined_max_step_seconds=1/480)
        write(out/'config.json',config)
        root=Path(__file__).resolve().parents[2]
        paths=[Path(__file__),Path(__file__).with_name('heart_top.py'),Path(native().__file__),
               root/'python/make_it_stand/verification/source_spin_material.py']
        identities={str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
        write(out/'source.json',dict(head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip(),
            files=identities,scope='LOADED_BINARY_AND_DECLARED_FILES_NOT_FULL_BUILD_ATTESTATION'))
        print('Building locked heart input',flush=True)
        problem=SpinArtifactProblem(original=(v,f),fixed_regions=[('shell',shell,1.)],
            design_regions=[('design',(dv,f),1.)],length_unit='25mm',mass_unit='18.75g',
            target={'mode':'SPINNING_TOP','axis':[0,0,1],'origin':[0,0,0],'support_points':[]},
            minimum_outer_wall=.015)
        write(out/'input.json',problem.identity)
        print('Running actual nonzero-objective source2014 optimization',flush=True)
        r=b.solve_and_seal_source_spin2014(problem,source_settings=geometry,solver_settings=solver,boundary_settings=boundary)
        write(out/'native-run.json',r.identity); write(out/'evidence.json',r.evidence)
        for name,mesh in [('material',r.material),('cavity',r.cavity)]:
            if mesh is not None and len(mesh[0]): trimesh.Trimesh(*mesh,process=False).export(out/(name+'.obj'))
        if r.artifact is not None: save_spin_artifact(r.artifact,out/'sealed')
        unchanged=all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in identities.items())
        write(out/'summary.json',dict(status=r.status,sealed=r.artifact is not None,seal_error=r.seal_error,
            elapsed_seconds=time.monotonic()-start,source_unchanged=unchanged,
            problem_digest=problem.input_digest,outer_moments=integrate(v,f).tolist(),
            scientific_digest=r.scientific_digest,original_paper_case=False))
        print(json.dumps(dict(status=r.status,sealed=r.artifact is not None,seal_error=r.seal_error)),flush=True)
        return 0 if r.artifact is not None and unchanged else 2
    except Exception as e:
        write(out/'failure.json',dict(error=str(e),traceback=traceback.format_exc(),elapsed_seconds=time.monotonic()-start))
        raise


if __name__=='__main__': raise SystemExit(main())
