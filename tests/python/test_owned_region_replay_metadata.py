"""Replay the actual owned construction, preserving operand coordinate policy."""
import gc

import numpy as np

from make_it_stand.verification.spin_common import native
from test_exact_material_region import box


def test_owned_difference_exports_real_immutable_operand_bindings():
    n = native()
    outer = n.make_exact_material_region(*box())
    carved = n.make_exact_material_region(*box((.25,)*3, (.75,)*3))
    region = n.make_exact_interior_difference_region(outer, carved)
    assert hasattr(region, 'construction_operands'), 'owned replay metadata is missing'
    operands = region.construction_operands
    assert [r.facts['region_digest'] for r in operands] == [outer.facts['region_digest'], carved.facts['region_digest']]
    operands.clear()
    del outer, carved
    gc.collect()
    assert len(region.construction_operands) == 2
    assert region.boundary_coordinate_policy == 'RATIONAL_INDEXED'


def test_binary_operand_export_preserves_original_binary64_bits():
    vertices, faces = box()
    vertices[0, 0] = -0.
    region = native().make_exact_material_region(vertices, faces)
    assert hasattr(region, 'boundary_coordinate_policy'), 'exact coordinate policy is missing'
    assert region.boundary_coordinate_policy == 'BINARY64_INDEXED'
    actual = region.binary64_boundary_vertices
    assert actual.tobytes() == vertices.tobytes()
    assert np.signbit(actual[0, 0])
    actual[0, 0] = 999
    assert np.signbit(region.binary64_boundary_vertices[0, 0])
    assert region.construction_operands == []
