"""Reload buoyant material, reintegrate labelled regions, and recheck equilibrium.

The mass checker is independent signed tetrahedron integration. Hydrostatic
geometry uses the existing clipping kernel, not an independent fluid simulator.
No optimizer is called. Passing is scoped method evidence, not manufacturing,
random Table 1 statistics, author-asset admission or physical validation.
"""
from pathlib import Path

import numpy as np

from ..spin import compute_hydrostatic_geometry
from .mass_moments import TriangleMesh, compute_mass_properties
from .spin_bundle import load_spin_artifact
from .spin_problem import MeshSnapshot
from .spin_common import digest, require
from .buoyant_application import FAMILY, PROTOCOL, locked_problem, _a_problem
from .faster_material import recompute_sealed_buoyant_stability_at_pose


def _moments(mesh, density, origin, axes):
    if mesh.explicit_empty:
        return np.zeros(10)
    v, f = mesh.arrays()
    p = compute_mass_properties(TriangleMesh(((v-origin)@axes).tolist(), f.tolist()), density)
    q = p.second_moment
    return np.array([p.mass, *p.first_moment, q[0][0], q[1][1], q[2][2],
                     q[1][2], q[0][2], q[0][1]])


def _same_region(actual, mesh, density, label):
    require(actual.mesh.identity == mesh.identity and actual.density == density,
            label+' geometry or density differs from input')


def _check_trace(s, initial, controls):
    result = s['solver_result']
    trace = result['trace']
    require(s['trace'] == trace and s['function_evaluations'] == len(trace),
            'evaluation ledger/count mismatch')
    require(s['status'] == result['status'] and s['accepted_steps'] == result['accepted_steps'],
            'solver terminal status/steps mismatch')
    require(isinstance(trace, list), 'invalid evaluation trace')
    if not trace:
        require(result['status'] == 'INVALID_REQUEST' and bool(result['diagnostic']) and
                result['accepted_steps'] == 0 and s['artifact_scientific_digest'] is None,
                'missing initial evaluation')
        return
    accepted = 0
    for i, event in enumerate(trace):
        require(event['evaluation'] == i, 'missing or reordered evaluation')
        require(type(event['accepted']) is bool and type(event['evaluation_succeeded']) is bool,
                'invalid evaluation state')
        require(type(event['iteration']) is int and 0 <= event['iteration'] <= controls['maximum_iterations'],
                'iteration outside registered budget')
        if event['accepted']:
            require(event['evaluation_succeeded'], 'failed evaluation accepted')
            accepted += 1
        if not event['evaluation_succeeded']:
            require(bool(event['diagnostic']), 'failure diagnostic omitted')
    require(accepted == result['accepted_steps'], 'accepted step ledger differs')
    require(len(trace) <= 1+controls['maximum_iterations']*controls['maximum_line_search_steps'],
            'evaluation ledger exceeds registered budget')
    require(trace[0]['iteration'] == 0 and trace[0]['accepted'] is False and
            np.allclose(trace[0]['plane_coefficients'], [*initial, -1.], rtol=1e-14, atol=0.),
            'initial evaluation differs from frozen initial plane')
    if result['status'] == 'CONVERGED':
        terminal = next((event for event in reversed(trace) if event['accepted']), trace[0])
        require(terminal['evaluation_succeeded'] and
                terminal['plane_coefficients'] == result['plane_coefficients'] and
                terminal['residual_norm'] <= controls['constraint_tolerance'],
                'converged terminal plane/residual differs from accepted trajectory')
        require(all(result.get(key) is True for key in ('equilibrium_constraints_verified',
                    'eq7_numerical_conditions_verified', 'stable_at_prescribed_waterline')),
                'converged solver lacks equilibrium/optimality/stability verification')


def audit_saved_buoyant_case(case, *, bundle_root):
    require(case.get('schema') == 'buoyant-table-method-v1' and case.get('solver_family') == FAMILY and
            case.get('initialization_protocol') == PROTOCOL, 'unsupported buoyant method record')
    p = locked_problem(case)
    require(digest(p) == case.get('buoyant_problem_sha256'), 'buoyant input identity mismatch')
    require(digest({k:v for k,v in case.items() if k != 'scientific_hash'}) == case.get('scientific_hash'),
            'buoyant run scientific identity mismatch')
    root = Path(bundle_root).resolve()
    axes, origin = np.asarray(p['target_axes']), np.asarray(p['target_origin'])
    d = p['displacement']
    outer = MeshSnapshot((p['outer']['vertices'], p['outer']['faces']))
    length = float(np.max(np.ptp((np.asarray(p['outer']['vertices'])-origin)@axes, axis=0)))
    scales = p['design']['density']*length**np.array([3,4,4,4,5,5,5,5,5,5])
    require(np.isfinite(scales).all() and (scales > 0).all(), 'invalid input normalization')
    runs = case.get('runs')
    require(isinstance(runs, list) and bool(runs), 'empty buoyant run collection')
    require([r['run'] for r in runs] == list(range(len(runs))) and
            len({r['seed'] for r in runs}) == len(runs), 'duplicate/missing run identity')
    qualified = 0
    maximum_error = 0.
    checked = []
    for run in runs:
        stages = run['stages']
        expected = ['A','B'] if 'stage_b' in p else ['A']
        require([s['stage'] for s in stages] == expected[:len(stages)] and bool(stages),
                'missing, duplicate or reordered stage')
        a = None
        stage_checks = []
        for index, s in enumerate(stages):
            inp = p if index == 0 else p['stage_b']
            require(s['waterline'] == inp['waterline'] and s['initial_parameters'] == inp['initial_parameters'] and
                    s['liquid_density'] == p['liquid_density'], 'stage differs from locked fluid/waterline/initial')
            _check_trace(s, inp['initial_parameters'], p['controls'])
            if s['artifact_scientific_digest'] is None:
                require(s['relative_bundle'] is None and
                        (s['status'] != 'CONVERGED' or s.get('seal_error')), 'missing converged material')
                stage_checks.append({'stage': s['stage'], 'qualified': False})
                continue
            path = root/s['relative_bundle']
            require(not Path(s['relative_bundle']).is_absolute() and root in path.resolve().parents,
                    'material bundle path escapes evidence root')
            expected_problem = _a_problem(p).input_digest if index == 0 else s['problem_digest']
            require(s['problem_digest'] == expected_problem, 'material problem binding mismatch')
            artifact = load_spin_artifact(path, expected_problem_digest=expected_problem)
            require(artifact.scientific_digest == s['artifact_scientific_digest'], 'sealed material digest mismatch')
            seal = s['seal_identity'] if index == 0 else s['seal_identity']['b_seal']
            require(seal['solver_result_digest'] == digest(s['solver_result']) and
                    seal['plane_coefficients_target'] == s['solver_result']['plane_coefficients'] and
                    seal['artifact_scientific_digest'] == artifact.scientific_digest,
                    'solver ledger/terminal plane differs from sealed material binding')
            problem = artifact.problem
            require(artifact.outer.identity == outer.identity and problem.length_unit == p['length_unit'] and
                    problem.mass_unit == p['mass_unit'] and problem.minimum_outer_wall == p['minimum_outer_wall'],
                    'saved material outer/units/wall differ from frozen input')
            if index == 1:
                require(a is not None and s['source_a_scientific_digest'] == a.scientific_digest,
                        'B source A binding mismatch')
                require(len(problem.design_regions) == 1, 'unexpected B design regions')
                _same_region(problem.design_regions[0], a.design_results[0][2], p['design']['density'], 'B design')
                expected_fixed = { 'A:'+r.region_id:(r.mesh,r.density) for r in a.problem.fixed_regions
                                   if not r.mesh.explicit_empty }
                adensity = {r.region_id:r.density for r in a.problem.design_regions}
                expected_fixed.update({'A:'+i:(m,adensity[i]) for i,m,_ in a.design_results if not m.explicit_empty})
                require({r.region_id for r in problem.fixed_regions} == set(expected_fixed), 'B lost A material region')
                for r in problem.fixed_regions:
                    _same_region(r, *expected_fixed[r.region_id], 'B fixed A')
            total = np.zeros(10)
            for r in problem.fixed_regions:
                total += _moments(r.mesh, r.density, origin, axes)
            densities = {r.region_id:r.density for r in problem.design_regions}
            for i, retained, _ in artifact.design_results:
                total += _moments(retained, densities[i], origin, axes)
            for claim in (s['total_raw_moments'], s['solver_result']['total_raw_moments']):
                claim = np.asarray(claim, dtype=float)
                require(claim.shape == (10,) and np.isfinite(claim).all(), 'invalid saved mass moments')
                error = float(np.max(np.abs(total-claim)/scales))
                maximum_error = max(maximum_error, error)
                require(error <= 1e-6, 'saved material does not reproduce solver moments')
            hydro = compute_hydrostatic_geometry((np.asarray(d['vertices'])-origin)@axes, d['faces'],
                waterline_origin=(np.asarray(inp['waterline']['origin'])-origin)@axes, waterline_normal=(0,0,1))
            residual = np.array([total[0]-p['liquid_density']*hydro.displaced_volume,
                total[1]-hydro.center_of_buoyancy[0]*total[0], total[2]-hydro.center_of_buoyancy[1]*total[0]])
            normalized = residual/scales[:3]
            margin = float(np.min(np.linalg.eigvalsh(hydro.waterline_central_second_moment))/hydro.displaced_volume -
                           (total[3]/total[0]-hydro.center_of_buoyancy[2]))
            require(np.max(np.abs(normalized)) <= 1e-6 and margin > 0,
                    'independent prescribed-waterline equilibrium or stability failed')
            if index == 1:
                require(len(s.get('pose_checks', [])) == 2, 'missing nominal/reverse free-waterline checks')
                reverse = axes @ np.diag([1.,-1.,-1.]) @ axes.T
                for saved, rotation in zip(s['pose_checks'], (np.eye(3), reverse)):
                    actual = recompute_sealed_buoyant_stability_at_pose(artifact,
                        displacement_vertices=d['vertices'], displacement_faces=d['faces'], rotation=rotation,
                        translation=(0,0,0), reference_origin=inp['waterline']['origin'],
                        waterline_normal=inp['waterline']['normal'], liquid_density=p['liquid_density'])
                    require(digest(actual) == digest(saved), 'pose check does not reproduce actual geometry')
            ok = s['status'] == 'CONVERGED'
            stage_checks.append({'stage': s['stage'], 'qualified': ok, 'residual': normalized.tolist(),
                                 'minimum_metacentric_margin': margin})
            if index == 0:
                a = artifact
        ok = len(stages) == len(expected) and all(s['qualified'] for s in stage_checks)
        require((run['status'] == 'CONVERGED') == ok, 'run terminal state conflicts with independent material checks')
        qualified += int(ok)
        checked.append({'run':run['run'], 'qualified':ok, 'stages':stage_checks})
    require(type(case.get('passed')) is bool and case['passed'] == (qualified == len(runs)),
            'reported pass conflicts with independent method checks')
    return {'schema':'saved-buoyant-method-audit-v1', 'qualified_runs':qualified, 'total_runs':len(runs),
            'maximum_normalized_moment_error':maximum_error, 'runs':checked,
            'scope':'INDEPENDENT_MASS_MOMENTS_AND_RECOMPUTED_HYDROSTATICS_NOT_MANUFACTURING_OR_PAPER_APPROVAL'}
