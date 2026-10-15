"""Independently reintegrate sealed material, then run the fixed physical matrix."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import trimesh
from heart_top import integrate, properties, simulate
from run import write


def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('run',type=Path)
    p.add_argument('--output-name',default='validation')
    a=p.parse_args(); root=a.run; dest=root/a.output_name; dest.mkdir(exist_ok=False)
    summary=json.loads((root/'summary.json').read_text()); config=json.loads((root/'config.json').read_text())
    if not summary['sealed'] or not summary['source_unchanged']: raise ValueError('not a current sealed result')
    from make_it_stand.verification.spin_bundle import load_spin_artifact
    print('Independently reloading exact sealed heart',flush=True)
    artifact=load_spin_artifact(root/'sealed',expected_problem_digest=summary['problem_digest'])
    original=artifact.problem.original.arrays(); outer=artifact.outer.arrays()
    if not all(np.array_equal(a,b) for a,b in zip(original,outer)): raise ValueError('external shape changed')
    material=artifact.final_material.arrays(); cavity=artifact.newly_carved.arrays()
    arrays=dict(outer_v=original[0],outer_f=original[1],material_v=material[0],material_f=material[1],
                cavity_v=cavity[0],cavity_f=cavity[1])
    np.savez(dest/'actual_meshes.npz',**arrays)
    s0=integrate(*original); s1=integrate(*material); removed=integrate(*cavity)
    error=float(np.max(np.abs(s0-s1-removed))/max(np.abs(s0).max(),1e-30))
    if error>1e-10: raise ValueError('independent partition moments inconsistent')
    native_run=json.loads((root/'native-run.json').read_text())
    native_s=np.array(native_run['native_run']['solve']['binary_moments'],float)
    native_error=float(np.max(np.abs(s1-native_s))/max(np.abs(s1).max(),1e-30))
    if native_error>1e-6: raise ValueError('actual final mesh differs from solver moments')
    model_records=[]; runs=[]
    powers=np.array([3,4,4,4,5,5,5,5,5,5])
    for label,s,mesh in [('original',s0,original),('carved',s1,material)]:
        mass,c,I=properties(s*config['physical_density_kg_per_m3']*config['meters_per_model_unit']**powers)
        tm=trimesh.Trimesh(*mesh,process=False)
        if not tm.is_watertight or not tm.is_winding_consistent: raise ValueError('invalid material mesh')
        model_records.append(dict(model=label,mass_kg=mass,center_m=c.tolist(),pivot_inertia_kg_m2=I.tolist(),
            boundary_components=len(tm.split(only_watertight=False)),geometry_moments=s.tolist()))
        for hz in config['rotation_hz']:
            for direction in range(config['azimuth_count']):
                print(f'{label}: {hz:g}Hz direction{direction}',flush=True)
                q=simulate(mass,c,I,[0,0,2*np.pi*hz],duration=config['duration_seconds'],
                    tilt_degrees=config['tilt_degrees'],azimuth=direction*2*np.pi/config['azimuth_count'],
                    max_step=config['max_step_seconds'])
                record=dict(model=label,hz=hz,direction=direction,**q)
                if q['relative_energy_drift']>1e-6 or q['relative_vertical_momentum_drift']>1e-6:
                    raise ValueError('physical invariant check failed')
                # All directions recorded, not only the most favorable. This is
                # a finite-time observation, NOT proof of indefinite stability.
                record['within_15deg_for_8s']=q['max_tilt_degrees']<=15
                write(dest/f'{label}-{hz:g}hz-dir{direction}.json',record)
                if direction==0:
                    fine=simulate(mass,c,I,[0,0,2*np.pi*hz],duration=config['duration_seconds'],
                        tilt_degrees=config['tilt_degrees'],azimuth=0,max_step=config['refined_max_step_seconds'])
                    angle_error=float(np.max(np.abs(np.array(q['tilt_degrees'])-fine['tilt_degrees'])))
                    record['step_refinement_max_angle_difference_degrees']=angle_error
                    write(dest/f'{label}-{hz:g}hz-refined.json',fine)
                    record['step_refinement_passed']=angle_error<=1e-3
                    if angle_error>1e-3:
                        print('NUMERICAL_REFINEMENT_FAILED: preserved; continuing remaining scenarios',flush=True)
                runs.append({k:v for k,v in record.items() if k not in ('time','quaternions','tilt_degrees')})
    report=dict(scope='RECONSTRUCTED_HEART_IDEAL_FIXED_PIVOT_NOT_CONTACT_NOT_PAPER_ORIGINAL',
        mass_partition_error=error,native_moment_relative_error=native_error,outer_unchanged=True,
        meshes_sha256=hashlib.sha256((dest/'actual_meshes.npz').read_bytes()).hexdigest(),
        artifact_scientific_digest=artifact.scientific_digest,models=model_records,scenarios=runs,
        tested_duration_seconds=8,gravity_m_s2=9.81,contact_simulation=False,physical_print_validation=False,
        moment_order='V,x,y,z,xx,yy,zz,yz,xz,xy',
        evaluator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        dynamics_sha256=hashlib.sha256(Path(__file__).with_name('heart_top.py').read_bytes()).hexdigest(),
        all_step_refinements_passed=all(r.get('step_refinement_passed',True) for r in runs))
    write(dest/'report.json',report)
    print(f"Recorded {len(runs)} scenarios; all step refinements passed: {report['all_step_refinements_passed']}",flush=True)


if __name__=='__main__': main()
