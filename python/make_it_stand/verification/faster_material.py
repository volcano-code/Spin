"""Locked hollow-only Faster candidates/runs and independently sealed material.

Native coordinates are input-normalized model coordinates; quadric coefficients
are native p (matrix-entry convention) in the derived input-normalized target
frame. This module never accepts caller-supplied moments or repaired geometry.
The immutable result binds the run to its material seal. Saving just its artifact
with spin_bundle persists only the material seal, not this run binding.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import math

import numpy as np

from ..spin import (SpinEvaluatorConfig, SpinItFasterRequest, SpinSolverConfig,
                    _finite_tuple, _native_module, evaluate_spin_it_faster_quadric,
                    integrate_raw_mass_moments, solve_spin_it_faster)
from .spin_artifact import SealedSpinArtifact, _recheck_problem
from .spin_common import FRAME, ORDER, SCOPE, canonical, digest, parse, require
from .spin_problem import MeshSnapshot, SpinArtifactProblem


@dataclass(frozen=True, init=False)
class FasterMaterialResult:
    """Run binding kept separate from the independently reloadable material seal."""
    status: str
    artifact: SealedSpinArtifact | None
    seal_error: str | None
    scientific_digest: str
    _identity: str

    def __init__(self):
        raise TypeError('use evaluate_and_seal_faster_candidate or solve_and_seal_faster')

    @property
    def identity(self):
        return parse(self._identity)


def _result(status, identity, artifact=None, seal_error=None):
    identity = {**identity, 'status': status, 'seal_error': seal_error,
                'artifact_scientific_digest': None if artifact is None else artifact.scientific_digest}
    result = object.__new__(FasterMaterialResult)
    for name, value in dict(status=status, artifact=artifact, seal_error=seal_error,
                            scientific_digest=digest(identity), _identity=canonical(identity)).items():
        object.__setattr__(result, name, value)
    return result


def _transform_moments(moments, axes, translation):
    """Push x to A*x+t; mass unchanged, full symmetric second moments retained."""
    m, lx, ly, lz, xx, yy, zz, yz, xz, xy = moments
    a, t = np.asarray(axes), np.asarray(translation)
    first = a @ np.array([lx, ly, lz])
    second = a @ np.array([[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]]) @ a.T
    second += np.outer(first, t) + np.outer(t, first) + m*np.outer(t, t)
    first += m*t
    return [float(m), *first.tolist(), *np.diag(second).tolist(),
            float(second[1,2]), float(second[0,2]), float(second[0,1])]


def _target_axes(axis):
    # Scaling first avoids overflow for finite very large direction components.
    z = np.asarray(axis, dtype=float)
    z = z / np.max(np.abs(z))
    z = z / np.linalg.norm(z)
    seed = np.eye(3)[int(np.argmin(np.abs(z)))]
    x = seed - z*np.dot(seed, z)
    x /= np.linalg.norm(x)
    return np.column_stack((x, np.cross(z, x), z))


def _prepare(problem, evaluator_settings, solver_settings):
    n = _native_module()
    if not all(hasattr(n, name) for name in ('SpinEvaluatorSettings', 'evaluate_spin_it_faster_quadric',
                                            'solve_spin_it_faster', 'integrate_raw_mass_moments')):
        raise RuntimeError('strict Faster material bridge requires native Faster evaluator support')
    require(isinstance(problem, SpinArtifactProblem), 'expected locked SpinArtifactProblem')
    problem = _recheck_problem(problem)
    require(problem.validation_scope == SCOPE, 'NOT_SUPPORTED: material scope')
    require(len(problem.design_regions) == 1, 'NOT_SUPPORTED: exactly one design region required')
    require(isinstance(evaluator_settings, SpinEvaluatorConfig) and
            evaluator_settings.derivative_method == 'DISCRETE_ELLIPSOID',
            'strict bridge requires explicit DISCRETE_ELLIPSOID settings')
    require(isinstance(solver_settings, SpinSolverConfig) and
            solver_settings.equation_profile == 'DUAL_MONOMIAL_RAW_V1',
            'strict bridge requires explicit DUAL_MONOMIAL_RAW_V1 settings')
    target = problem.target
    require(target['mode'] != 'BUOYANT_STABILITY', 'NOT_SUPPORTED: buoyancy requires actual hydrostatics')
    if target['mode'] == 'STATIC_BALANCE':
        require(target['support_points'] == [target['origin']],
                'NOT_SUPPORTED: static target needs exactly one support point equal to target origin')
    else:
        require(not target['support_points'], 'NOT_SUPPORTED: rotational target with support points')
    norm = problem.normalization
    origin, length, density = np.array(norm['origin']), norm['L'], norm['R']
    def normalized_mesh(mesh):
        v, f = mesh.arrays()
        return ((v-origin)/length).tolist(), f.tolist()
    fixed = []
    for region in problem.fixed_regions:
        v, f = normalized_mesh(region.mesh)
        fixed.append(integrate_raw_mass_moments(v, f, region.density/density))
    moments = tuple(math.fsum(row[i] for row in fixed) for i in range(10))
    outer_v, outer_f = normalized_mesh(problem.original)
    design_v, design_f = normalized_mesh(problem.design_regions[0].mesh)
    axes = _target_axes(target['axis'])
    target_origin = (np.array(target['origin']) - origin)/length
    request = SpinItFasterRequest(outer_v, outer_f, design_v, design_f,
        fixed_shell_moments=moments, target_origin=target_origin.tolist(), target_axes=axes.tolist(),
        mode=target['mode'], material_density=problem.design_regions[0].density/density,
        original_units_per_normalized_unit=length, minimum_wall_thickness=problem.minimum_outer_wall/length,
        settings=solver_settings)
    identity = {'schema': 'faster-locked-material-run-v1', 'problem_digest': problem.input_digest,
                'scope': SCOPE, 'project_status': 'PARTIAL_BLOCKED',
                'solver_target_feasibility': 'UNKNOWN', 'physical_validation': 'NOT_IN_SCOPE',
                'evaluator_settings': asdict(evaluator_settings), 'solver_settings': asdict(solver_settings),
                'native_request_digest': digest(asdict(request)),
                'coefficient_frame': 'INPUT_NORMALIZED_TARGET_V1',
                'coefficient_convention': 'NATIVE_MATRIX_ENTRY_P',
                'normalization': norm, 'target_axes': axes.tolist(), 'target_origin': target_origin.tolist(),
                'fixed_moments_input_normalized': list(moments)}
    return problem, request, identity


def _seal(problem, request, retained, carved, outer, total_moments):
    # Compare actual native outer output against the exact request arrays. The
    # separately sealed raw outer is the locked original, never a roundtrip mesh.
    for actual, expected, dtype in zip(outer, (request.outer_vertices, request.outer_faces), ('<f8', '<i4')):
        require(np.asarray(actual, dtype=dtype).shape == np.asarray(expected, dtype=dtype).shape and
                np.asarray(actual, dtype=dtype).tobytes() == np.asarray(expected, dtype=dtype).tobytes(),
                'native output outer bytes differ from normalized locked request')
    norm = problem.normalization
    def model_mesh(mesh):
        v, f = mesh
        if not len(v) and not len(f):
            return MeshSnapshot(None)
        return MeshSnapshot((np.asarray(v)*norm['L'] + np.array(norm['origin']), np.asarray(f, dtype=np.int32)))
    kept, removed = model_mesh(retained), model_mesh(carved)
    if removed.explicit_empty:
        material = problem.original
    else:
        ov, of = problem.original.arrays()
        cv, cf = removed.arrays()
        material = MeshSnapshot((np.vstack((ov, cv)), np.vstack((of, cf[:, ::-1]+len(ov)))))
    claim = _transform_moments(total_moments, request.target_axes, request.target_origin)
    return SealedSpinArtifact.create(problem, outer=problem.original,
        design_results=[(problem.design_regions[0].region_id, kept, removed)],
        final_material=material, newly_carved=removed,
        solver_claim={'frame': FRAME, 'order': list(ORDER), 'moments': claim})


def evaluate_and_seal_faster_candidate(problem, coefficients, *, evaluator_settings, solver_settings):
    """Evaluate p and -p exactly once each; seal material without optimizing."""
    problem, request, identity = _prepare(problem, evaluator_settings, solver_settings)
    coefficients = _finite_tuple(coefficients, 10, 'candidate coefficients')
    retained = evaluate_spin_it_faster_quadric(request, coefficients, evaluator_settings=evaluator_settings)
    require(retained.derivative_defined, 'unsupported discrete derivative: '+retained.derivative_diagnostic)
    carved = evaluate_spin_it_faster_quadric(request, tuple(-x for x in coefficients),
                                            evaluator_settings=evaluator_settings)
    require(carved.derivative_defined, 'unsupported discrete complement derivative: '+carved.derivative_diagnostic)
    require(retained.request_hash == carved.request_hash, 'native evaluation request identity mismatch')
    artifact = _seal(problem, request, (retained.material_vertices, retained.material_faces),
                     (carved.material_vertices, carved.material_faces),
                     (retained.outer_vertices, retained.outer_faces), retained.total_raw_moments)
    identity.update(operation='CANDIDATE', coefficients=list(coefficients), native_evaluation_count=2,
                    native_evaluations=[asdict(retained), asdict(carved)])
    return _result('NOT_OPTIMIZED_CANDIDATE', identity, artifact)


def solve_and_seal_faster(problem, *, initial_quadrics, evaluator_settings, solver_settings):
    """Preserve an actual native run; only converged exported geometry may seal."""
    problem, request, identity = _prepare(problem, evaluator_settings, solver_settings)
    initial = tuple(_finite_tuple(row, 10, 'initial quadric') for row in initial_quadrics)
    require(bool(initial), 'strict bridge requires an explicit nonempty initial sequence')
    run = solve_spin_it_faster(request, initial_quadrics=initial, evaluator_settings=evaluator_settings)
    identity.update(operation='SOLVE', initial_quadrics=initial, native_run=run.diagnostics,
                    derivative_metadata=run.derivative_metadata)
    if run.status != 'CONVERGED':
        return _result(run.status, identity)
    try:
        artifact = _seal(problem, request, (run.material_vertices, run.material_faces),
            (run.cavity_vertices, run.cavity_faces), (run.outer_vertices, run.outer_faces),
            run.final_raw_moments)
    except (ValueError, RuntimeError) as exc:
        return _result(run.status, identity, seal_error=str(exc))
    return _result(run.status, identity, artifact)
