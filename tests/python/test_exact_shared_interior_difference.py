"""Reuse owned exact operands without rebuilding their boundary components."""
import numpy as np
import pytest

from make_it_stand import spin
from test_spin_it_faster import cube


def factory():
    native = spin._native_module()
    operation = getattr(native, 'make_exact_interior_difference_region', None)
    assert callable(operation), 'missing exact shared-operand interior difference construction'
    return native, operation


@pytest.mark.parametrize('two_cavities', [False, True])
def test_shared_difference_equals_independently_reconstructed_boundary(two_cavities):
    native, operation = factory()
    outer = native.make_exact_material_region(*cube(-4., 4.))
    vertices, faces = map(np.asarray, cube(-2., -1.))
    if two_cavities:
        extra_v, extra_f = map(np.asarray, cube(1., 2.))
        faces = np.vstack((faces, extra_f+len(vertices)))
        vertices = np.vstack((vertices, extra_v))
    carved = native.make_exact_material_region(vertices, faces)
    result = operation(outer, carved)
    assert result.facts['completed'] and result.facts['valid']
    assert result.facts['material_components'] == 1
    assert result.facts['bounded_void_components'] == (2 if two_cavities else 1)
    boundary = result.rational_boundary
    assert result.matches_rational_boundary(boundary['vertices'], boundary['faces'])
    independent = native.make_exact_material_region_from_rationals(boundary['vertices'], boundary['faces'])
    comparison = native.relate_exact_material_regions(result, independent)
    assert comparison['completed'] and comparison['equal']
    assert comparison['positive_volume_overlap'] is True
    assert comparison['symmetric_difference_volume_exact'] in ('0', '0/1')


def test_shared_difference_rejects_noninterior_or_unowned_operands():
    native, operation = factory()
    outer = native.make_exact_material_region(*cube(-1., 1.))
    outside = native.make_exact_material_region(*cube(2., 3.))
    touching = native.make_exact_material_region(*cube(0., 1.))
    with pytest.raises(ValueError):
        operation(outer, outside)
    with pytest.raises(ValueError):
        operation(outer, touching)
    with pytest.raises(TypeError):
        operation(dict(valid=True), outside)


def test_exported_boundary_is_a_copy_and_cannot_change_its_exact_proof():
    native, operation = factory()
    result = operation(native.make_exact_material_region(*cube(-2., 2.)),
                       native.make_exact_material_region(*cube(-1., 1.)))
    boundary = result.rational_boundary
    boundary['vertices'][0][0] = '999/1'
    assert not result.matches_rational_boundary(boundary['vertices'], boundary['faces'])
    fresh = result.rational_boundary
    assert fresh['vertices'][0][0] != '999/1'
    assert result.matches_rational_boundary(fresh['vertices'], fresh['faces'])
