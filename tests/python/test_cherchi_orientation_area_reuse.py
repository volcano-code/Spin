"""Reuse exact original face area without changing emitted/coherent volume signs."""
from fractions import Fraction
import numpy as np
import pytest
from make_it_stand.verification import cherchi_motion_differential as motion


@pytest.mark.parametrize('offset,width', [(Fraction(0),Fraction(1)), (Fraction(2**50),Fraction(1,2**60))])
def test_exact_orientation_reuses_each_area_and_preserves_local_flip_audit(monkeypatch, offset, width):
    points = [(offset,offset,offset), (offset+width,offset,offset),
              (offset,offset+width,offset), (offset,offset,offset+width)]
    vertices = np.asarray([[float(x) for x in row] for row in points])
    expected = np.asarray([[1,2,3],[0,3,2],[0,1,3],[0,2,1]],dtype=np.int64)
    emitted = expected.copy(); emitted[2] = emitted[2][[0,2,1]]
    calls = []
    original = motion._rational_area_vector
    def observe(*args):
        calls.append(True)
        return original(*args)
    monkeypatch.setattr(motion,'_rational_area_vector',observe)
    faces, audit = motion.normalize_closed_component_orientation(vertices, emitted, rational_vertices=points)
    np.testing.assert_array_equal(faces, expected)
    volume = width**3/6
    component = audit['components'][0]
    assert component['exact_signed_volume'] == f'{volume.numerator}/{volume.denominator}'
    assert component['exact_emitted_signed_volume'] == f'{volume.numerator}/{volume.denominator}'
    assert component['relative_flip_count'] == 1
    assert len(calls) == 4, 'exact triangle areas are recomputed for emitted/coherent volumes'
