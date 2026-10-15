"""The explicit terminal strategy keeps full material auditing and roundtrip."""
import pytest

from make_it_stand import spin
from make_it_stand.verification import faster_material as fm
from make_it_stand.verification.spin_bundle import save_spin_artifact
from test_faster_material_bridge import problem
from test_exact_material_region import box, combine


def fixture():
    p = problem(density=1., fixed_density=1.)
    evaluator = spin.SpinEvaluatorConfig(derivative_method='DISCRETE_FULL_RANK_QUADRIC')
    solver = spin.SpinSolverConfig(equation_profile='DUAL_MONOMIAL_RAW_V1')
    p, request, _ = fm._prepare(p, evaluator, solver)
    design = box((-.25,)*3, (.25,)*3)
    carved = box((-.125,)*3, (.125,)*3)
    retained = combine(design, (carved[0], carved[1][:, ::-1]))
    total = [63/64, 0, 0, 0, 341/4096, 341/4096, 341/4096, 0, 0, 0]
    return p, request, retained, carved, total


def test_owned_terminal_seal_preserves_physics_and_reloads_actual_operation(tmp_path):
    p, request, retained, carved, total = fixture()
    args = (p, request, retained, carved, (request.outer_vertices, request.outer_faces), total)
    baseline = fm._seal(*args)
    owned = fm._seal(*args, terminal_region_profile='OWNED_INTERIOR_DIFFERENCE_V1', reuse_problem_validation=True)
    assert owned.final_material.exact_region().facts['algorithm'] == 'strict-interior-nef-difference-shared-boundary-v1'
    for key in ('physical', 'normalized', 'clearance'):
        assert baseline.facts[key] == owned.facts[key]
    loaded = save_spin_artifact(owned, tmp_path/'bundle')
    assert loaded.facts == owned.facts


def test_terminal_region_profile_cannot_be_an_unregistered_strategy():
    p, request, retained, carved, total = fixture()
    with pytest.raises(ValueError, match='terminal region profile'):
        fm._seal(p, request, retained, carved, (request.outer_vertices, request.outer_faces), total,
                 terminal_region_profile='TRUST_BOUNDARY')
