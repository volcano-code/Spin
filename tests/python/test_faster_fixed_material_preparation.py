"""Fixed-material moment integration must reuse authenticated input geometry."""

import numpy as np

from make_it_stand import spin
from make_it_stand.verification import faster_material
from make_it_stand.verification.spin_problem import SpinArtifactProblem
from test_exact_material_region import box, combine


def _problem():
    outer = box((-2.0,) * 3, (2.0,) * 3)
    design = box((-1.0,) * 3, (1.0,) * 3)
    fixed = combine(outer, (design[0], design[1][:, ::-1]))
    return SpinArtifactProblem(
        original=outer,
        fixed_regions=[('shell', fixed, 1.3)],
        design_regions=[('design', design, 2.0)],
        length_unit='m', mass_unit='kg', minimum_outer_wall=0.25,
        target={'mode': 'SPINNING_TOP', 'axis': [1.0, 2.0, 3.0],
                'origin': [0.25, -0.5, 0.3], 'support_points': []})


def test_fixed_moment_integration_does_not_reconstruct_an_exact_region(monkeypatch):
    problem = _problem()
    evaluator = spin.SpinEvaluatorConfig(
        derivative_method='DISCRETE_FULL_RANK_QUADRIC', icosphere_subdivisions=0)
    solver = spin.SpinSolverConfig(equation_profile='DUAL_MONOMIAL_RAW_V1')
    _, request, _ = faster_material._prepare(problem, evaluator, solver)

    class NoReconstruction:
        def __init__(self, *args, **kwargs):
            raise AssertionError('fixed moment integration rebuilt an exact material region')

    monkeypatch.setattr(faster_material, 'MeshSnapshot', NoReconstruction)
    prepared = faster_material._prepare_external_fixed_material(problem, request)
    axes = np.asarray(request.target_axes)
    origin = np.asarray(request.target_origin)
    expected = faster_material._transform_moments(
        request.fixed_shell_moments, axes.T, -axes.T @ origin)
    np.testing.assert_allclose(prepared.moments, expected, rtol=0.0, atol=1e-10)
