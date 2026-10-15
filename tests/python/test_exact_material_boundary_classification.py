"""Exact shell nesting without constructing a Nef for every Newton trial."""
from fractions import Fraction

import pytest

from make_it_stand.verification.spin_common import native
import make_it_stand.verification.spin_problem as spin_problem
from test_exact_material_region import box, combine, l_prism


def classify(mesh, proof=None):
    operation = getattr(native(), 'audit_exact_rational_material_boundary', None)
    assert callable(operation), 'exact point-location material boundary classification missing'
    vertices, faces = mesh
    rows = [[f'{q.numerator}/{q.denominator}' for q in (Fraction(float(x)) for x in row)] for row in vertices]
    return operation(rows, faces, proof)


@pytest.mark.parametrize('nested', [False, True])
def test_valid_classification_matches_independent_nef_hierarchy(nested):
    outer = box((0,)*3, (8,)*3)
    if nested:
        cavity = box((1,)*3, (7,)*3)
        island = box((2,)*3, (6,)*3)
        mesh = combine(outer, (cavity[0], cavity[1][:, ::-1]), island)
    else:
        mesh = combine(l_prism(), box((1.5, 1.5, .2), (2.5, 2.5, .8)))
    result = classify(mesh)
    reference = native().make_exact_material_region(*mesh).facts
    assert result['completed'] and result['classification_valid'] and reference['valid']
    for key in ('parents', 'depths', 'original_signs', 'material_components', 'bounded_void_components'):
        assert result[key] == reference[key]
    assert result['scope'] == 'exact-embedded-boundary-material-classification-only'


def test_wrong_nesting_sign_is_rejected():
    mesh = combine(box((0,)*3, (4,)*3), box((1,)*3, (2,)*3))
    result = classify(mesh)
    assert result['completed'] and not result['classification_valid']
    assert 'orientation' in result['error']


def test_touching_or_crossing_shells_are_rejected_before_nesting_queries():
    for offset in (.5, 1.):
        result = classify(combine(box(), box((offset, 0, 0), (offset+1, 1, 1))))
        assert not result['classification_valid']
        assert not result['boundary']['boundary_valid']


def test_fixed_proof_is_owned_and_cannot_be_a_report_dictionary():
    operation = getattr(native(), 'audit_exact_rational_material_boundary', None)
    assert callable(operation), 'exact material boundary classifier missing'
    vertices, faces = box()
    rows = [[f'{int(x)}/1' for x in row] for row in vertices]
    proof = native().verify_fixed_rational_boundary(rows, faces)
    result = operation(rows, faces, proof)
    assert result['classification_valid']
    assert result['fixed_proof_identity'] == proof.reference_identity
    with pytest.raises((TypeError, RuntimeError)):
        operation(rows, faces, dict(boundary_valid=True))


def test_negative_root_and_extremely_small_positive_solid_are_distinguished():
    vertices, faces = box((0,)*3, (1e-30,)*3)
    result = classify((vertices, faces))
    assert result['classification_valid'] and result['material_components'] == 1
    reversed_result = classify((vertices, faces[:, ::-1]))
    assert not reversed_result['classification_valid']


def test_formal_classified_snapshot_uses_one_complete_boundary_audit(monkeypatch):
    vertices, faces = box()
    rows = [[f'{int(x)}/1' for x in row] for row in vertices]
    real = native()
    calls = []
    class ClassifierOnly:
        def audit_exact_rational_material_boundary(self, *args):
            calls.append(args)
            return real.audit_exact_rational_material_boundary(*args)
        def audit_exact_rational_triangle_boundary(self, *args):
            raise AssertionError('formal classified path repeated its boundary audit')
    monkeypatch.setattr(spin_problem, 'native', lambda: ClassifierOnly())
    snapshot, mode = spin_problem.rational_solver_candidate_snapshot(
        rows, faces, allow_signed_classification=True)
    assert mode == 'EXACT_RATIONAL_SIGNED_BOUNDARY_CLASSIFICATION'
    assert len(calls) == 1
    assert snapshot.classification_proof['parents'] == [-1]
