"""Final snapshots may export an owned exact region without rebuilding it."""
from unittest.mock import patch

import numpy as np
import pytest

from make_it_stand.verification.spin_problem import MeshSnapshot, SpinArtifactProblem
from make_it_stand.verification.spin_artifact import SealedSpinArtifact
from make_it_stand.verification.spin_common import native, FRAME, ORDER
from test_exact_material_region import box, combine


def exported(region):
    factory = getattr(MeshSnapshot, 'from_exact_region', None)
    assert callable(factory), 'owned exact snapshot export is missing'
    return factory(region)


def hollow_region():
    n = native()
    return n.make_exact_interior_difference_region(n.make_exact_material_region(*box()),
        n.make_exact_material_region(*box((.25,)*3, (.75,)*3)))


def test_export_retains_actual_owned_region_without_new_nef_construction():
    region = hollow_region()
    with patch.object(native(), 'make_exact_material_region_from_rationals',
                      side_effect=AssertionError('unexpected second exact construction')):
        result = exported(region)
    assert result.exact_region() is region
    assert region.matches_rational_boundary(result.rational_vertices, result.arrays()[1])
    derived = result.arrays()[0]
    derived[0, 0] = 999
    assert result.arrays()[0][0, 0] != 999


def test_exact_export_preserves_nonbinary_rational_coordinates():
    vertices, faces = box()
    rows = [['1/3' if x == 0 and axis == 0 else f'{int(x)}/1'
             for axis, x in enumerate(row)] for row in vertices]
    region = native().make_exact_material_region_from_rationals(rows, faces)
    result = exported(region)
    assert result.rational_vertices == tuple(map(tuple, rows))
    assert result.arrays()[0][0, 0] == 1/3


def test_self_reported_or_invalid_regions_are_rejected():
    with pytest.raises((TypeError, ValueError)):
        exported(dict(completed=True, valid=True))
    invalid = native().make_exact_material_region(*combine(box(), box((.2,)*3, (.8,)*3)))
    assert not invalid.facts['valid']
    with pytest.raises(ValueError):
        exported(invalid)


def test_owned_export_still_passes_independent_complete_material_seal():
    outer, cavity = box(), box((.25,)*3, (.75,)*3)
    retained = combine(outer, (cavity[0], cavity[1][:, ::-1]))
    problem = SpinArtifactProblem(original=outer, fixed_regions=[], design_regions=[('d', outer, 1.)],
        length_unit='m', mass_unit='kg', minimum_outer_wall=.25,
        target=dict(mode='SPINNING_TOP', axis=[0., 0., 1.], origin=[0., 0., 0.], support_points=[]))
    claim = dict(frame=FRAME, order=list(ORDER), moments=[7/8, 0, 0, 0, 31/384, 31/384, 31/384, 0, 0, 0])
    def seal(material):
        return SealedSpinArtifact.create(problem, outer=outer, design_results=[('d', retained, cavity)],
            final_material=material, newly_carved=cavity, solver_claim=claim)
    baseline = seal(retained)
    candidate = seal(exported(hollow_region()))
    for key in ('physical', 'normalized', 'clearance', 'material_components', 'bounded_void_components'):
        assert baseline.facts[key] == candidate.facts[key]
    relation = native().relate_exact_material_regions(baseline.final_material.exact_region(),
                                                     candidate.final_material.exact_region())
    assert relation['completed'] and relation['equal']
