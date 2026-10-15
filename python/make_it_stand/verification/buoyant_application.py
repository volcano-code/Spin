"""Prescribed-waterline method dispatch using measured geometry and saved A→B.

This bounded entry supports a sealed outer displacement hull and one explicit
initial plane per stage. It is not the author's random Table 1 protocol, does
not infer flooding semantics, and cannot grant manufacturing or paper approval.
"""
from dataclasses import replace
from pathlib import Path

import numpy as np

from ..spin import SpinItFasterRequest, compute_hydrostatic_geometry, solve_spin_it_faster_buoyant_plane
from .faster_material import (seal_buoyant_plane_result_as_artifact,
    solve_floating_island_b_buoyant_from_sealed_a, recompute_sealed_buoyant_stability_at_pose,
    _raw_moments, _json_value)
from .spin_problem import MeshSnapshot, SpinArtifactProblem
from .spin_bundle import save_spin_artifact
from .spin_common import canonical, parse, digest, require, scalar

FAMILY = 'BUOYANT_PLANE_GEOMETRY_V1'
PROTOCOL = 'LOCKED_EXPLICIT_BUOYANT_PLANES_V1'


def _mesh(value):
    return MeshSnapshot((value['vertices'], value['faces']))


def locked_problem(case):
    require(isinstance(case.get('buoyant_problem'), dict),
            'buoyant_problem with actual displacement and fixed geometry is required')
    p = parse(canonical(case['buoyant_problem']))
    require(p['displacement']['semantics'] == 'SEALED_OUTER_HULL',
            'unsupported displacement semantics; flooding is not inferred')
    require(_mesh(p['outer']).identity == _mesh(p['displacement']).identity,
            'sealed displacement hull differs from actual outer boundary')
    _mesh(p['design'])
    scalar(p['design']['density'], 'design density', positive=True)
    scalar(p['liquid_density'], 'liquid density', positive=True)
    require(p['fixed_regions'], 'actual nonzero fixed material is required')
    ids = [r['id'] for r in p['fixed_regions']]
    require(len(ids) == len(set(ids)), 'duplicate fixed material ID')
    for r in p['fixed_regions']:
        _mesh(r)
        scalar(r['density'], 'fixed density', positive=True)
    axes = np.asarray(p['target_axes'], dtype=float)
    origin = np.asarray(p['target_origin'], dtype=float)
    require(axes.shape == (3,3) and origin.shape == (3,) and
            np.isfinite(axes).all() and np.isfinite(origin).all() and
            np.allclose(axes.T@axes, np.eye(3), rtol=0., atol=1e-10) and
            abs(np.linalg.det(axes)-1.) <= 1e-10, 'invalid target frame')
    controls = p['controls']
    for key in ('maximum_iterations', 'maximum_line_search_steps'):
        require(type(controls[key]) is int and controls[key] > 0, 'invalid '+key)
    require(0 < scalar(controls['constraint_tolerance'], 'constraint tolerance') <= 1e-10,
            'prescribed-waterline tolerance relaxed')
    for stage in [p] + ([p['stage_b']] if 'stage_b' in p else []):
        for key in ('origin', 'normal'):
            value = np.asarray(stage['waterline'][key], dtype=float)
            require(value.shape == (3,) and np.isfinite(value).all(), 'invalid waterline '+key)
        normal = np.asarray(stage['waterline']['normal'], dtype=float)
        require(np.linalg.norm(normal) > 0 and np.allclose(
            axes.T@(normal/np.linalg.norm(normal)), (0,0,1), rtol=0., atol=1e-10),
            'waterline normal must align with target vertical')
        initial = np.asarray(stage['initial_parameters'], dtype=float)
        require(initial.shape == (3,) and np.isfinite(initial).all(), 'invalid initial plane')
    if 'stage_b' in p:
        require(set(p['stage_b']) == {'design_from', 'waterline', 'initial_parameters'} and
                p['stage_b']['design_from'] == 'A_CAVITY',
                'B must consume A cavity with the same fluid, frame and design density')
    return p


def buoyant_case_plan(case, starts, max_restarts, seeds):
    require(type(starts) is int and starts > 0 and len(seeds) == starts,
            'invalid buoyant run count')
    require(max_restarts == 1, 'explicit buoyant plane protocol has one attempt; random restarts not implemented')
    require(len(set(seeds)) == starts, 'duplicate run labels')
    p = locked_problem(case)
    return {'id': case['id'], 'mode': case['mode'], 'solver_family': FAMILY,
            'initialization_protocol': PROTOCOL, 'buoyant_problem': p,
            'buoyant_problem_sha256': digest(p),
            'runs': [{'run': i, 'seed': seed, 'seed_role': 'REPEAT_LABEL_NOT_RANDOM_INITIALIZATION'}
                     for i, seed in enumerate(seeds)]}


def _base_request(p, design):
    v, f = design.arrays()
    return SpinItFasterRequest(p['outer']['vertices'], p['outer']['faces'], v.tolist(), f.tolist(),
        material_density=p['design']['density'], target_origin=p['target_origin'],
        target_axes=p['target_axes'], minimum_wall_thickness=p['minimum_outer_wall'])


def _a_request(p):
    request = _base_request(p, _mesh(p['design']))
    axes, origin = np.asarray(p['target_axes']), np.asarray(p['target_origin'])
    fixed = np.zeros(10)
    volume = 0.
    for r in p['fixed_regions']:
        mesh = _mesh(r)
        v, f = mesh.arrays()
        fixed += _raw_moments(mesh, r['density'])
        volume += float(_raw_moments(mesh, 1.)[0])
    d = p['displacement']
    hydro = compute_hydrostatic_geometry((np.asarray(d['vertices'])-origin)@axes, d['faces'],
        waterline_origin=(np.asarray(p['waterline']['origin'])-origin)@axes,
        waterline_normal=(0,0,1))
    return replace(request, mode='BUOYANT_STABILITY', fixed_shell_moments=tuple(fixed),
        fixed_shell_volume=volume, fixed_shell_density=float(fixed[0]/volume),
        liquid_density=p['liquid_density'], displaced_liquid_volume=hydro.displaced_volume,
        center_of_buoyancy_target=hydro.center_of_buoyancy, waterline_area=hydro.waterline_area,
        waterline_central_second_moment=hydro.waterline_central_second_moment)


def _a_problem(p):
    return SpinArtifactProblem(original=_mesh(p['outer']),
        fixed_regions=[(r['id'], _mesh(r), r['density']) for r in p['fixed_regions']],
        design_regions=[('design:A', _mesh(p['design']), p['design']['density'])],
        length_unit=p['length_unit'], mass_unit=p['mass_unit'],
        target={'mode': 'BUOYANT_STABILITY', 'axis': p['waterline']['normal'],
                'origin': p['waterline']['origin'], 'support_points': []},
        minimum_outer_wall=p['minimum_outer_wall'])


def _stage_record(name, solved, p, waterline, initial):
    result = _json_value(solved)
    return {'stage': name, 'status': result['status'], 'solver_result': result,
            'trace': result['trace'], 'function_evaluations': len(result['trace']),
            'accepted_steps': result['accepted_steps'], 'total_raw_moments': result['total_raw_moments'],
            'liquid_density': p['liquid_density'], 'waterline': waterline,
            'initial_parameters': initial, 'artifact_scientific_digest': None,
            'relative_bundle': None, 'seal_error': None}


def _save_stage(stage, sealed, root, relative):
    require(sealed.artifact is not None, 'no material produced: '+str(sealed.identity))
    stage['attempted_bundle'] = relative
    artifact = save_spin_artifact(sealed.artifact, root/relative)
    stage.update(relative_bundle=relative, artifact_scientific_digest=artifact.scientific_digest,
                 problem_digest=artifact.problem.input_digest, seal_identity=sealed.identity)
    return artifact


def run_buoyant_case(case, starts, max_restarts, bundle_root, seeds):
    plan = buoyant_case_plan(case, starts, max_restarts, seeds)
    require(bundle_root is not None, 'buoyant dispatch requires a bundle_root for reloadable material')
    root = Path(bundle_root).resolve()
    root.mkdir(parents=True, exist_ok=True)
    p = plan['buoyant_problem']
    problem = _a_problem(p)
    request = _a_request(p)
    runs = []
    for label in plan['runs']:
        # Explicit fixed initial planes: repetitions are never represented as random statistics.
        solved = solve_spin_it_faster_buoyant_plane(request,
            initial_parameters=p['initial_parameters'], **p['controls'])
        a = _stage_record('A', solved, p, p['waterline'], p['initial_parameters'])
        row = {**label, 'stages': [a], 'status': solved['status']}
        runs.append(row)
        if solved['status'] != 'CONVERGED':
            continue
        try:
            sealed = seal_buoyant_plane_result_as_artifact(problem, request, solved,
                displacement_vertices=p['displacement']['vertices'], displacement_faces=p['displacement']['faces'],
                waterline_origin=p['waterline']['origin'], waterline_normal=p['waterline']['normal'])
            prefix = digest({'id': case['id'], 'problem': p})[:16]+'-run-'+str(label['run'])
            artifact_a = _save_stage(a, sealed, root, prefix+'-A')
        except (ValueError, RuntimeError, OSError) as error:
            a['seal_error'] = {'type': type(error).__name__, 'message': str(error)}
            row['status'] = 'MATERIAL_VALIDATION_FAILED'
            continue
        if 'stage_b' in p:
            b_input = p['stage_b']
            design_b = artifact_a.design_results[0][2]
            try:
                result_b = solve_floating_island_b_buoyant_from_sealed_a(
                    _base_request(p, design_b), artifact_a, region_prefix='A',
                    source_to_request_rotation=np.eye(3), source_to_request_translation=(0,0,0),
                    source_length_to_request_unit=1., reference_density=1.,
                    source_length_unit=p['length_unit'], target_length_unit=p['length_unit'],
                    source_mass_unit=p['mass_unit'], target_mass_unit=p['mass_unit'],
                    displacement_vertices=p['displacement']['vertices'], displacement_faces=p['displacement']['faces'],
                    waterline_origin=b_input['waterline']['origin'], waterline_normal=b_input['waterline']['normal'],
                    liquid_density=p['liquid_density'], initial_parameters=b_input['initial_parameters'],
                    capture_seal_failure=True, **p['controls'])
                b = _stage_record('B', result_b.identity['solver_result'], p,
                                  b_input['waterline'], b_input['initial_parameters'])
                b['source_a_scientific_digest'] = artifact_a.scientific_digest
                b['composition_identity'] = result_b.identity
                b['seal_error'] = result_b.identity.get('seal_failure')
                row['stages'].append(b)
                row['status'] = result_b.status
                if result_b.status != 'CONVERGED':
                    continue
                artifact_b = _save_stage(b, result_b, root, prefix+'-B')
                axes = np.asarray(p['target_axes'])
                reverse = axes @ np.diag([1.,-1.,-1.]) @ axes.T
                b['pose_checks'] = [recompute_sealed_buoyant_stability_at_pose(artifact_b,
                    displacement_vertices=p['displacement']['vertices'], displacement_faces=p['displacement']['faces'],
                    rotation=rotation, translation=(0,0,0), reference_origin=b_input['waterline']['origin'],
                    waterline_normal=b_input['waterline']['normal'], liquid_density=p['liquid_density'])
                    for rotation in (np.eye(3), reverse)]
            except (ValueError, RuntimeError, OSError) as error:
                row['status'] = 'B_CHAIN_FAILED'
                row['failure'] = {'type': type(error).__name__, 'message': str(error)}
                if len(row['stages']) == 2:
                    row['stages'][1]['seal_error'] = row['failure']
    payload = {**{k:v for k,v in plan.items() if k != 'runs'}, 'schema': 'buoyant-table-method-v1',
        'asset_class': case.get('asset_class'), 'runs': runs,
        'passed': all(r['status'] == 'CONVERGED' and
                      all(s['artifact_scientific_digest'] for s in r['stages']) for r in runs),
        'scope': 'PRESCRIBED_WATERLINE_METHOD_NOT_MANUFACTURING_OR_ORIGINAL_PAPER_APPROVAL',
        'formal_paper_admission': 'NOT_GRANTED_BY_METHOD_RUNNER'}
    payload['scientific_hash'] = digest(payload)
    return payload
