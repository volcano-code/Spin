"""Rational and binary64 entity distinctions must survive input inspection."""
import unittest
import numpy as np
import trimesh
import _make_it_stand as native
from admit_rounded_input import exact_rows
from inspect_rounded_projection import inspect_region


class ProjectionTests(unittest.TestCase):
    def test_exactly_representable_box_matches(self):
        mesh = trimesh.creation.box(extents=[2, 2, 2])
        result = inspect_region(native, exact_rows(mesh.vertices), mesh.faces.tolist())
        self.assertTrue(result['binary64_region_valid'])
        self.assertTrue(result['represented_sets_equal'])
        self.assertEqual(result['maximum_coordinate_rounding_exact'], '0')

    def test_valid_float_projection_is_not_exact_entity_equality(self):
        mesh = trimesh.creation.box(extents=[2, 2, 2])
        rows = exact_rows(mesh.vertices)
        rows = [['-1/1' if row[0]=='-1/1' else '1/3', *row[1:]] for row in rows]
        result = inspect_region(native, rows, mesh.faces.tolist())
        self.assertTrue(result['binary64_region_valid'])
        self.assertFalse(result['represented_sets_equal'])
        self.assertGreater(result['moment_errors'][0], 0)
        self.assertNotEqual(result['maximum_coordinate_rounding_exact'], '0')


if __name__ == '__main__':
    unittest.main()
