"""Actual admitted flower dynamics, with the whole preregistered scene matrix."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np


def read_json(path):
    return json.loads(Path(path).read_text())


def scene_record(result, *, scene_id, model, variant, phase_index, dt_seconds,
                 trajectory_sha256, actual_steps):
    """Merge worker diagnostics without allowing them to replace scene identity.

    Flexible workers already report dt_seconds. Verify matching duplicates,
    retain every non-identity diagnostic, and publish each identity field once.
    """
    identity=dict(id=scene_id,model=model,variant=variant,phase_index=phase_index,
                  dt_seconds=dt_seconds,trajectory_sha256=trajectory_sha256,actual_steps=actual_steps)
    for key,value in identity.items():
        if key in result and result[key]!=value:
            raise ValueError('SCENE_IDENTITY_MISMATCH: '+key)
    if 'dt' in result and result['dt']!=dt_seconds:
        raise ValueError('SCENE_IDENTITY_MISMATCH: dt')
    return {**result,**identity}


def validate_clamp_mass_agreement(record, clamps):
    actual = sum(x['mass_kg'] for x in clamps.values())
    expected = record['clamp']['mass_kg']
    # Only account for binary64 summation/integration roundoff. This is not
    # the geometric 1e-6 gate and never changes the support's authoritative mass.
    bound = 64*np.finfo(float).eps*max(abs(actual), abs(expected))
    if not math.isfinite(actual+expected) or min(actual,expected)<=0 or abs(actual-expected)>bound:
        raise ValueError('CLAMP_MASS_IDENTITY_MISMATCH')
    return dict(admitted_mass_kg=expected, independently_integrated_mass_kg=actual,
                absolute_difference_kg=abs(actual-expected), arithmetic_roundoff_bound_kg=bound)


def compare_time_step_trajectories(coarse, fine, model):
    """Compare actual samples on the common grid without interpolating motion."""
    t0,t1=np.asarray(coarse['t']),np.asarray(fine['t'])
    if t0.ndim!=1 or t1.ndim!=1 or len(t0)<2 or len(t1)<len(t0):
        raise ValueError('INVALID_TIMESTEP_TRAJECTORY')
    if not np.isfinite(t0).all() or not np.isfinite(t1).all() or np.any(np.diff(t0)<=0) or np.any(np.diff(t1)<=0):
        raise ValueError('INVALID_TIMESTEP_TRAJECTORY')
    indices=np.searchsorted(t1,t0)
    indices=np.clip(indices,0,len(t1)-1)
    preceding=np.maximum(indices-1,0)
    indices=np.where(abs(t1[preceding]-t0)<abs(t1[indices]-t0),preceding,indices)
    time_error=float(np.max(abs(t1[indices]-t0)))
    if time_error>16*np.finfo(float).eps*max(1.,abs(float(t0[-1]))):
        raise ValueError('TIMESTEP_COMMON_GRID_MISMATCH')
    energy_key='energy_J' if model=='flexible' else 'energy'
    keys=['support_xy_m','angle','omega',energy_key,'reaction']
    if model=='flexible': keys+=['relative_energy_balance_error']
    for key in keys:
        if len(coarse[key])!=len(t0) or len(fine[key])!=len(t1) or not np.isfinite(coarse[key]).all() or not np.isfinite(fine[key]).all():
            raise ValueError('INVALID_TIMESTEP_TRAJECTORY:'+key)
    position_difference=np.linalg.norm(coarse['support_xy_m']-fine['support_xy_m'][indices],axis=1)
    result=dict(scope='ALIGNED_ACTUAL_SAMPLES_NO_AUTOMATIC_CONVERGENCE_LABEL',common_samples=len(t0),
        maximum_sample_time_difference_seconds=time_error,
        maximum_support_position_difference_m=float(position_difference.max()),
        maximum_angle_difference_rad=float(np.max(abs(coarse['angle']-fine['angle'][indices]))),
        maximum_angular_speed_difference_rad_s=float(np.max(abs(coarse['omega']-fine['omega'][indices]))),
        maximum_energy_difference_J=float(np.max(abs(coarse[energy_key]-fine[energy_key][indices]))),
        maximum_reaction_force_difference_N=float(np.max(np.linalg.norm(coarse['reaction'][:,:3]-fine['reaction'][indices,:3],axis=1))))
    if model=='flexible':
        errors=[np.asarray(x['relative_energy_balance_error']) for x in (coarse,fine)]
        result.update(maximum_energy_balance_error_coarse=float(errors[0].max()),
            maximum_energy_balance_error_fine=float(errors[1].max()),
            maximum_energy_balance_error_difference=float(np.max(abs(errors[0]-errors[1][indices]))))
        for key,field in [('support_velocity_m_s','maximum_support_velocity_difference_m_s'),
                          ('damping_loss_J','maximum_damping_loss_difference_J')]:
            if key in coarse and key in fine:
                difference=coarse[key]-fine[key][indices]
                result[field]=float(np.max(np.linalg.norm(difference,axis=1) if difference.ndim==2 else abs(difference)))
    else:
        for label,data in [('coarse',coarse),('fine',fine)]:
            e0=float(data['energy'][0])
            if e0<=0:raise ValueError('NONPOSITIVE_INITIAL_ENERGY')
            result['maximum_relative_energy_drift_'+label]=float(np.max(abs(data['energy']/e0-1)))
    return result


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _checked_properties(value):
    mass=float(value['mass_kg']);center=np.asarray(value['center_m'],float)
    inertia=np.asarray(value['inertia_kg_m2'],float)
    if mass<=0 or not math.isfinite(mass) or center.shape!=(3,) or inertia.shape!=(3,3) or \
            not np.isfinite(center).all() or not np.isfinite(inertia).all() or \
            not np.allclose(inertia,inertia.T,rtol=1e-12,atol=0.) or np.linalg.eigvalsh(inertia).min()<=0:
        raise ValueError('INVALID_SAVED_PHYSICS_PROPERTIES')
    return dict(mass_kg=mass,center_m=center.tolist(),inertia_kg_m2=inertia.tolist())


def _saved_balance(folder, case, optimized):
    """A separate application gate using saved independent M and first moments.

    The planar dogleg problem targets the origin; no alternate target is inferred
    from the retained object. Normalization uses frozen input geometry only.
    """
    from paper_flower_fixture import exact_extrusion_moments
    frozen=read_json(folder/'input.json')
    manifest=read_json(folder/'flower-cases.json')['cases'][0]
    problem=manifest['planar_problem']
    if case.get('id')!=manifest['id'] or case.get('planar_problem')!=problem or \
            case.get('solver_family')!='PLANAR_DOGLEG_V1' or case.get('mode')!='FIDGET_SPINNER':
        raise ValueError('SAVED_OPTIMIZATION_APPLICATION_MISMATCH')
    if any(name in problem for name in ('target_origin','target_center','target_axes')):
        raise ValueError('UNSUPPORTED_PLANAR_TARGET_CONTRACT')
    length=float(frozen['meters_per_model_unit']);density=float(frozen['pla_density_kg_m3'])
    points=np.asarray([p for loop in problem['design_region_loops'] for p in loop],float)
    span=float(np.linalg.norm(np.ptp(points,axis=0)))*length
    entities=[dict(loops=problem['design_region_loops'],z_interval=problem['extrusion_interval'],
                   volume_density=problem['extrusion_volume_density']),*problem['fixed_extrusion_layers']]
    initial_mass=float(sum(exact_extrusion_moments(e)[0] for e in entities))*density*length**3
    if not math.isfinite(span+initial_mass) or min(span,initial_mass,length,density)<=0:
        raise ValueError('INVALID_FROZEN_BALANCE_NORMALIZATION')
    first=np.asarray(optimized['center_m'][:2])*optimized['mass_kg']
    residual=first/(initial_mass*span)
    success=next(a for a in case['attempts'] if a['seed']==0 and a['status']=='CONVERGED')
    circle=success['solver_result'];radius=float(circle['radius'])
    if not math.isfinite(radius) or radius<=0 or not math.isclose(radius,math.hypot(*circle['center']),rel_tol=1e-14,abs_tol=0.):
        raise ValueError('SAVED_CIRCLE_TARGET_AXIS_MISMATCH')
    if not np.isfinite(residual).all() or float(np.linalg.norm(residual))>1e-6:
        raise ValueError('SAVED_MATERIAL_BALANCE_CONSTRAINT_FAILED')
    return dict(scope='SAVED_OBJ_FIRST_MOMENT_BALANCE_NOT_DOGLEG_OPTIMALITY',
                target_axis_xy_model=[0.,0.],initial_input_mass_kg=initial_mass,
                input_design_diagonal_m=span,independent_saved_first_moments_kg_m=first.tolist(),
                normalized_first_moment_residual=residual.tolist(),tolerance=1e-6,passed=True)


def admit(output, audit_directory, folder, *, saved_material_check=None, assembly_reference=None):
    output,audit_directory,folder=map(Path,(output,audit_directory,folder))
    record = read_json(output/'output-preparation.json')
    shared=record.get('derived_geometry_policy')=='BOUNDED_SHARED_RATIONAL_ARRANGEMENT_V1'
    if shared and saved_material_check is None:
        raise ValueError('SAVED_MATERIAL_CHECK_REQUIRED_FOR_SHARED_POLICY')
    if saved_material_check is not None:
        saved_material_check=Path(saved_material_check)
        audit=read_json(saved_material_check)
        if audit.get('schema')!='paper-flower-saved-independent-check-v1':
            raise ValueError('SAVED_MATERIAL_CHECK_SCHEMA_MISMATCH')
        required=[folder/'input.json',folder/'flower-cases.json',folder/'input-sha256.json',
                  folder/'solve-job-01/output/application.json',output/'output-preparation.json']
        if assembly_reference is not None:required.append(Path(assembly_reference))
        sources={Path(path).resolve():sha for path,sha in audit.get('sources',{}).items()}
        if any(path.resolve() not in sources or sources[path.resolve()]!=_sha(path) for path in required) or \
                any(not path.is_file() or _sha(path)!=sha for path,sha in sources.items()):
            raise ValueError('SAVED_CHECK_SOURCE_IDENTITY_MISMATCH')
        if not record.get('files'):
            raise ValueError('SAVED_CHECK_ACTUAL_OBJ_REQUIRED')
        partitions=audit.get('partitions',{})
        subsets=audit.get('bearing_subsets',{})
        errors=audit.get('moment_errors',{})
        if set(partitions)!=set(('editable','middle','original')) or \
                any(any(row.get(key) is not True for key in ('completed','complete_partition','no_positive_volume_overlap')) for row in partitions.values()) or \
                set(subsets)!=set(('outer-race','bearing-seat')) or \
                any(row.get('completed') is not True or row.get('a_subset_b') is not True for row in subsets.values()) or \
                not isinstance(audit.get('material_pairs'),list) or \
                any(row['relation'].get('completed') is not True or row['relation'].get('positive_volume_overlap') is not False for row in audit['material_pairs']):
            raise ValueError('SAVED_CHECK_GEOMETRY_PROOF_NOT_PASSED')
        if set(errors)!=set(('original','optimized','authority')) or \
                any(np.asarray(row).shape!=(10,) or not np.isfinite(row).all() or \
                    np.min(row)<0 or np.max(row)>1e-6 for row in errors.values()):
            raise ValueError('SAVED_CHECK_INDEPENDENT_MOMENT_PROOF_NOT_PASSED')
        physical=['retained',*record['rotor_layers'],*record['clamp_layers']]
        expected_pairs={tuple(sorted((a,b))) for i,a in enumerate(physical) for b in physical[i+1:]}
        actual_pairs=[tuple(sorted((row['a'],row['b']))) for row in audit['material_pairs']]
        if len(physical)!=len(set(physical)) or any(name not in record['files'] for name in physical) or \
                len(actual_pairs)!=len(set(actual_pairs)) or set(actual_pairs)!=expected_pairs:
            raise ValueError('SAVED_CHECK_OVERLAP_PROOF_INCOMPLETE')
    else:
        audit = read_json(audit_directory/'exact-audit.json')
        if audit.get('preparation_sha256') != _sha(output/'output-preparation.json'):
            raise ValueError('MATERIAL_AUDIT_IDENTITY_MISMATCH')
    if audit.get('full_material_gate') != 'PASSED' or audit.get('simulation_admission') is not True:
        raise ValueError('MATERIAL_NOT_ADMITTED_FOR_PHYSICS')
    for info in record['files'].values():
        path=output/info['file']
        if path.is_symlink() or path.parent.resolve()!=output.resolve() or _sha(path) != info['sha256']:
            raise ValueError('STALE_PHYSICS_MESH')
    if record.get('input_sha256')!=hashlib.sha256((folder/'input.json').read_bytes()).hexdigest():
        raise ValueError('MATERIAL_INPUT_IDENTITY_MISMATCH')
    application_paths=[folder/f'solve-job-0{i}/output/application.json' for i in (1,2)]
    if record.get('application_sha256')!=hashlib.sha256(application_paths[0].read_bytes()).hexdigest():
        raise ValueError('MATERIAL_APPLICATION_IDENTITY_MISMATCH')
    applications = [read_json(path) for path in application_paths]
    cases = [a['cases'][0] for a in applications]
    if cases[0]['scientific_hash'] != cases[1]['scientific_hash']:
        raise ValueError('INDEPENDENT_SOLVER_REPEAT_MISMATCH')
    if saved_material_check is not None:
        selected=next(a for a in cases[0]['attempts'] if a['seed']==0 and a['status']=='CONVERGED')['solver_result']
        if audit.get('selected_seed')!=0 or record.get('selected_seed')!=0 or \
                audit.get('circle_center')!=record.get('circle_center') or \
                audit.get('circle_radius')!=record.get('circle_radius') or \
                audit.get('circle_center')!=selected['center'] or audit.get('circle_radius')!=selected['radius']:
            raise ValueError('SAVED_CHECK_CIRCLE_IDENTITY_MISMATCH')
    for i in (1, 2):
        check_path=folder/f'check-job-0{i}/output/scientific.json'
        checked_file=read_json(check_path)
        if saved_material_check is not None and (checked_file.get('schema')!='saved-planar-material-moment-audit-v1' or \
                checked_file.get('input',{}).get('sha256')!=_sha(application_paths[i-1]) or \
                Path(checked_file.get('input',{}).get('path','')).resolve()!=application_paths[i-1].resolve()):
            raise ValueError('INDEPENDENT_APPLICATION_CHECK_SOURCE_MISMATCH')
        checked = checked_file['cases'][0]['attempts']
        case = cases[i-1]
        success = next(a for a in case['attempts'] if a['seed'] == 0 and a['status'] == 'CONVERGED')
        proof = next(a for a in checked if a['run'] == success['run'] and a['restart'] == success['restart'])
        if not proof['moments_match'] or not proof['balance_constraints_pass']:
            raise ValueError('SAVED_MATERIAL_INDEPENDENT_CHECK_FAILED')
    if saved_material_check is not None:
        record={**record, 'original':_checked_properties(audit['original']),
                'optimized':_checked_properties(audit['derived_optimized']),
                'clamp':_checked_properties(audit['clamp'])}
        audit={**audit, 'independent_saved_balance':_saved_balance(folder,cases[0],record['optimized']),
               'consumed_saved_check':dict(path=str(saved_material_check.resolve()),sha256=_sha(saved_material_check),
                    material_admission='SAVED_INDEPENDENT_OBJ_GATE',optimization_admission='TWO_SAVED_APPLICATION_CHECKS_AND_INDEPENDENT_OBJ_BALANCE',
                    application_sha256=[_sha(path) for path in application_paths],
                    application_check_sha256=[_sha(folder/f'check-job-0{i}/output/scientific.json') for i in (1,2)])}
    return record, audit


def run(folder, material, material_audit, grip, out, *, saved_material_check=None):
    import trimesh
    from bearing_physics import simulate as ideal
    from flexible_bearing import support_protocol, simulate as flexible
    record, audit = admit(material, material_audit, folder,saved_material_check=saved_material_check,
                          assembly_reference=grip/'assembly.json')
    # Measure the frozen clamp's geometry. These parts translate with the
    # support, not with the spinning flower. Rolling elements are not modeled.
    assembly = read_json(grip/'assembly.json')
    length = assembly['meters_per_model_unit']
    densities = dict(upper_cap=1200., lower_cap=1200., shaft=7800., inner_race=7800.)
    clamps = {}
    for name, density in densities.items():
        mesh = trimesh.load(grip/(name+'.obj'), force='mesh', process=False)
        if not mesh.is_watertight or not mesh.is_winding_consistent or mesh.volume <= 0:
            raise ValueError('INVALID_CLAMP_GEOMETRY')
        clamps[name] = dict(density_kg_m3=density, mass_kg=float(mesh.volume*length**3*density),
            source_sha256=hashlib.sha256((grip/(name+'.obj')).read_bytes()).hexdigest())
    legacy_clamp_check=validate_clamp_mass_agreement(record,clamps)
    admitted_clamps={}
    for name,density in densities.items():
        info=record['files'][name]
        mesh=trimesh.load(material/info['file'],force='mesh',process=False)
        admitted_clamps[name]=dict(density_kg_m3=density,mass_kg=float(mesh.volume*density),source_sha256=info['sha256'])
    admitted_clamp_check=validate_clamp_mass_agreement(record,admitted_clamps)
    protocol = support_protocol(record['original']['mass_kg'], record['clamp']['mass_kg'])
    out.mkdir(parents=True, exist_ok=False)
    rows = []
    with (out/'events.jsonl').open('x') as events:
        for model in ('ideal', 'flexible'):
            for variant in ('original', 'optimized'):
                props = record[variant]
                for phase_index in range(8):
                    for dt in (.0005, .00025):
                        name = f'{model}-{variant}-phase{phase_index}-dt{dt}'
                        kwargs = dict(dt=dt, duration=8., hz=35., phase=phase_index*math.pi/4)
                        if model == 'flexible':
                            kwargs['protocol'] = protocol
                        result = (ideal if model == 'ideal' else flexible)(props['mass_kg'],
                            props['center_m'], np.array(props['inertia_kg_m2']), **kwargs)
                        trajectory = result.pop('records')
                        if model == 'flexible':
                            fields = ('t', 'angle', 'omega', 'energy_J', 'damping_loss_J', 'relative_energy_balance_error')
                            arrays = {field: np.array([r[field] for r in trajectory]) for field in fields}
                            arrays['support_xy_m'] = np.array([r['support_xy_m'] for r in trajectory])
                            arrays['support_velocity_m_s'] = np.array([r['support_velocity_m_s'] for r in trajectory])
                        else:
                            arrays = {field: np.array([r[field] for r in trajectory]) for field in ('t', 'angle', 'omega', 'energy')}
                            arrays['support_xy_m'] = np.zeros((len(trajectory), 2))
                        arrays['reaction'] = np.array([r['reaction'] for r in trajectory])
                        arrays['motor_torque'] = np.array([r['motor_torque'] for r in trajectory])
                        np.savez_compressed(out/(name+'.npz'), **arrays)
                        row = scene_record(result,scene_id=name, model=model, variant=variant, phase_index=phase_index,
                            dt_seconds=dt, trajectory_sha256=hashlib.sha256((out/(name+'.npz')).read_bytes()).hexdigest(),
                            actual_steps=len(trajectory)-1)
                        if model == 'ideal':
                            actual_force = float(np.linalg.norm(arrays['reaction'][1:, :2], axis=1).mean())
                            analytic_force = props['mass_kg']*np.linalg.norm(props['center_m'][:2])*(2*math.pi*35)**2
                            row.update(mean_horizontal_force_N=actual_force, analytic_horizontal_force_N=float(analytic_force),
                                analytic_force_absolute_error_N=abs(actual_force-analytic_force))
                        rows.append(row)
                        events.write(json.dumps({k:row[k] for k in ('id', 'actual_steps', 'trajectory_sha256')})+'\n')
                        events.flush()
                        print(name, flush=True)
    comparisons = []
    for variant in ('original', 'optimized'):
        for phase_index in range(8):
            for model in ('ideal','flexible'):
                pair = [r for r in rows if r['model']==model and r['variant']==variant and r['phase_index']==phase_index]
                pair.sort(key=lambda row:row.get('dt_seconds',row.get('dt')),reverse=True)
                with np.load(out/(pair[0]['id']+'.npz')) as coarse, np.load(out/(pair[1]['id']+'.npz')) as fine:
                    comparison=compare_time_step_trajectories(coarse,fine,model)
                if model=='flexible':
                    peaks = [r['maximum_support_displacement_m'] for r in pair]
                    comparison['peak_displacement_relative_difference']=abs(peaks[0]-peaks[1])/max(max(peaks),1e-12)
                comparisons.append(dict(variant=variant,phase_index=phase_index,model=model,**comparison))
    summary = dict(schema='admitted-flower-actual-physics-matrix-v1', classification='RECONSTRUCTED_METHOD_VALIDATION',
        original=record['original'], optimized=record['optimized'], support_protocol=protocol,
        fixed_clamp_geometry=admitted_clamps, independent_legacy_clamp_geometry=clamps,
        clamp_mass_checks=dict(legacy=legacy_clamp_check,admitted_si=admitted_clamp_check),
        models=['ideal', 'flexible'], initial_hz=35., duration_seconds=8.,
        phase_count=8, time_steps_seconds=[.0005, .00025], expected_scenes=64, recorded_scenes=len(rows),
        sources=dict(material_preparation_sha256=hashlib.sha256((material/'output-preparation.json').read_bytes()).hexdigest(),
                     **({'saved_material_check_sha256':audit['consumed_saved_check']['sha256'],
                         'saved_material_check_path':audit['consumed_saved_check']['path']}
                        if saved_material_check is not None else
                        {'exact_material_audit_sha256':_sha(material_audit/'exact-audit.json')})),
        material_admission=audit.get('consumed_saved_check',dict(material_admission='LEGACY_EXACT_AUDIT')),
        independent_saved_balance=audit.get('independent_saved_balance'),
        dt_comparisons=comparisons, scenes=rows, motor_drive=False,
        video_selection='phase0; fine timestep; no ranking or best-scene selection',
        real_hand_calibration=False, stopping_time_prediction=False,
        scientific_conclusions='See complete timestep, energy and vibration records; no automatic stability label')
    (out/'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    return summary


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    for name in ('folder', 'material', 'material-audit', 'grip', 'output'):
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--saved-material-check',type=Path,
                        help='Actual saved-independent-check.json; required for shared derived policy')
    args=parser.parse_args()
    run(args.folder, args.material, args.material_audit, args.grip, args.output,
        saved_material_check=args.saved_material_check)
