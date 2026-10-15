"""Fresh optimization of the complete shaft-equipped heart. Never reuse old carving."""
import argparse,hashlib,json,time,traceback
from dataclasses import asdict
from pathlib import Path
import numpy as np
import trimesh
from stem_geometry import make_input,write_mesh,ROOT
from run import write


def main():
    p=argparse.ArgumentParser(); p.add_argument('--output',type=Path,required=True)
    p.add_argument('--geometry-depth',type=int,default=5)
    p.add_argument('--optimization-level',type=int,default=4)
    a=p.parse_args()
    out=a.output; out.mkdir(exist_ok=False); start=time.monotonic()
    try:
        from make_it_stand.verification import source_spin_material as b
        from make_it_stand.verification.spin_problem import SpinArtifactProblem
        from make_it_stand.verification.spin_bundle import save_spin_artifact
        from make_it_stand.verification.spin_common import native
        sources=[Path(__file__),Path(__file__).with_name('stem_geometry.py'),
            ROOT/'experiments/reconstructed_heart_top/heart_top.py',Path(native().__file__),
            ROOT/'build/heart-stem-contact-tools/mesh-boolean',Path(b.__file__)]
        hashes={str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
        write(out/'source.json',dict(files=hashes,scope='DECLARED_FILES_AND_LOADED_BINARY_NOT_FULL_BUILD_ATTESTATION'))
        print('Constructing actual CGAL heart/shaft union and fixed-material partition',flush=True)
        shapes=make_input(out)
        # Coordinate change declared BEFORE solving: bottom tip is the origin.
        for shape in shapes.values(): shape.vertices=np.asarray(shape.vertices)+[0,0,.16]
        for name,shape in shapes.items(): write_mesh(out/(name+'.obj'),shape)
        geom=b.SourceGeometrySettings(root_origin=(-1.5,-1.5,-.1),root_side=3.,maximum_depth=a.geometry_depth)
        solver=b.SourceSolverSettings(objective='TOP_MASS_AND_LOWERING',initial_level=3,
            maximum_optimization_level=a.optimization_level,maximum_adaptive_passes=10,maximum_leaf_count=4000,
            initial_fill=.5,initial_phi=.2,laplacian_weight=.001,
            slqp=b.AuditedSlqpSettings(maximum_iterations=160))
        boundary=b.SourceBoundarySettings(maximum_input_cuboids=4000)
        write(out/'config.json',dict(asset='RECONSTRUCTED_HEART_WITH_FIXED_SHAFT_V1',
            meters_per_model_unit=.025,density_kg_m3=1200.,translation_model=[0,0,.16],
            shaft_radius_model=.065,tip_radius_model=.025,shaft_top_before_translation=2.25,
            geometry=asdict(geom),solver=asdict(solver),boundary=asdict(boundary),minimum_wall=.015,
            contact_protocol=dict(hz=[6,35,45],azimuths=8,tilt_degrees=2.,duration_seconds=8.,
                nominal_time_step=1/2000,refined_time_step=1/4000,lateral_friction=.8,
                rolling_friction=0.,spinning_friction=0.,linear_damping=0.,angular_damping=0.,
                disclosure='Zero phenomenological rolling/spin friction and damping; uncalibrated dry sliding contact, no maintained torque')))
        mesh=lambda name:(np.asarray(shapes[name].vertices),np.asarray(shapes[name].faces))
        print('Validating complete input partition',flush=True)
        problem=SpinArtifactProblem(original=mesh('outer'),fixed_regions=[('shaft-and-shell',mesh('fixed'),1.)],
            design_regions=[('design-excluding-shaft',mesh('design'),1.)],length_unit='25mm',mass_unit='18.75g',
            target={'mode':'SPINNING_TOP','axis':[0,0,1],'origin':[0,0,0],'support_points':[]},minimum_outer_wall=.015)
        write(out/'input.json',problem.identity)
        print('Running fresh source2014 optimization including shaft mass',flush=True)
        r=b.solve_and_seal_source_spin2014(problem,source_settings=geom,solver_settings=solver,boundary_settings=boundary)
        write(out/'native-run.json',r.identity); write(out/'evidence.json',r.evidence)
        for name,mesh in [('material',r.material),('cavity',r.cavity)]:
            if mesh is not None: write_mesh(out/(name+'.obj'),trimesh.Trimesh(*mesh,process=False))
        if r.artifact is not None:
            save_spin_artifact(r.artifact,out/'sealed')
            ov,of=r.artifact.outer.arrays(); mv,mf=r.artifact.final_material.arrays(); cv,cf=r.artifact.newly_carved.arrays()
            np.savez(out/'sealed-meshes.npz',outer_v=ov,outer_f=of,material_v=mv,material_f=mf,cavity_v=cv,cavity_f=cf)
        unchanged=all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in hashes.items())
        write(out/'summary.json',dict(status=r.status,sealed=r.artifact is not None,seal_error=r.seal_error,
            problem_digest=problem.input_digest,elapsed_seconds=time.monotonic()-start,source_unchanged=unchanged,
            scientific_digest=r.scientific_digest,author_original=False))
        print('RESULT',r.status,'sealed=',r.artifact is not None,'error=',r.seal_error,flush=True)
        return 0 if r.artifact is not None and unchanged else 2
    except Exception as e:
        write(out/'failure.json',dict(error=str(e),traceback=traceback.format_exc(),elapsed_seconds=time.monotonic()-start))
        raise


if __name__=='__main__': raise SystemExit(main())
