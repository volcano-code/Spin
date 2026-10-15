"""Exact multi-cavity separation, not a complete wall-thickness certificate.

Analytic expectations catch missed material-island boundaries, strict instead
of inclusive thresholds, component truncation and accidental seal admission.
"""
import importlib
from fractions import Fraction
import unittest

import numpy as np

from test_exact_material_region import box, combine


class InterCavityClearanceTests(unittest.TestCase):
    def setUp(self):
        name = 'make_it_stand.verification.inter_cavity_clearance'
        self.assertIsNotNone(importlib.util.find_spec(name),
                             'integrated inter-cavity checker is missing')
        self.audit = importlib.import_module(name).audit_inter_cavity_clearance

    def two_voids(self):
        return combine(box((-2, -2, -2), (4, 3, 3)),
                       box((0, 0, 0), (1, 1, 1), True),
                       box((1.25, 0, 0), (2.25, 1, 1), True))

    def check(self, mesh, threshold=.25, maximum_pairs=32):
        result = self.audit(mesh, minimum_clearance=threshold,
                            maximum_pairs=maximum_pairs)
        self.assertEqual(result['complete_wall_thickness'], 'UNKNOWN')
        self.assertEqual(result['manufacturing_admission'], 'NOT_GRANTED')
        self.assertEqual(result['physical_admission'], 'NOT_GRANTED')
        return result

    def test_two_voids_in_one_material_threshold_equality(self):
        mesh = self.two_voids()
        r = self.check(mesh)
        self.assertTrue(r['completed'])
        self.assertTrue(r['at_least_minimum'])
        self.assertEqual(r['squared_distance_exact'], '1/16')
        self.assertEqual((r['cavity_count'], r['pairs_required'],
                          r['pairs_checked']), (2, 1, 1))
        self.assertEqual(r['witness_cavities'], [1, 2])
        self.assertTrue(12 <= r['witness_faces'][0] < 24)
        self.assertTrue(24 <= r['witness_faces'][1] < 36)
        self.assertFalse(self.check(mesh, np.nextafter(.25, np.inf))['at_least_minimum'])
        self.assertEqual(r, self.check(mesh))

    def test_nested_cavity_in_material_island_includes_positive_shell(self):
        mesh = combine(box((0, 0, 0), (8, 8, 8)),
                       box((1, 1, 1), (7, 7, 7), True),
                       box((2, 2, 2), (6, 6, 6)),
                       box((2.25, 2.25, 2.25), (5.75, 5.75, 5.75), True))
        r = self.check(mesh, .5)
        self.assertEqual(r['cavity_count'], 2)
        self.assertEqual(r['cavities'][0]['boundary_shells'], [1, 2])
        self.assertEqual(r['cavities'][1]['boundary_shells'], [3])
        self.assertEqual(r['squared_distance_exact'], '1/16')
        self.assertFalse(r['at_least_minimum'])
        # Closest boundary of the large void is the material island, NOT its
        # negative outer shell. Testing negative shells alone gives 25/16.
        self.assertTrue(24 <= r['witness_faces'][0] < 36)
        self.assertTrue(36 <= r['witness_faces'][1] < 48)

    def test_non_axis_aligned_exact_metric(self):
        v, f = self.two_voids()
        v[:, 0] += .5 * v[:, 1]
        r = self.check((v, f), .22)
        self.assertEqual(Fraction(r['squared_distance_exact']), Fraction(1, 20))
        self.assertTrue(r['at_least_minimum'])
        self.assertFalse(self.check((v, f), .23)['at_least_minimum'])

    def test_exterior_air_gap_does_not_become_wall_failure(self):
        mesh = combine(box((0, 0, 0), (2, 2, 2)),
                       box((.125, .125, .125), (1.875, 1.875, 1.875), True),
                       box((2.125, 0, 0), (4.125, 2, 2)),
                       box((2.25, .125, .125), (4, 1.875, 1.875), True))
        r = self.check(mesh, .5)
        self.assertTrue(r['completed'])
        self.assertEqual(r['squared_distance_exact'], '9/64')
        self.assertFalse(r['at_least_minimum'])
        self.assertEqual(r['complete_wall_thickness'], 'UNKNOWN')

    def test_reordered_faces_and_vertices_preserve_distance_and_witness(self):
        from test_exact_partition_clearance import triangle_distance
        v, f = self.two_voids()
        # Independent Fraction oracle checks the returned original-index faces.
        rng = np.random.default_rng(20260921)
        order = rng.permutation(len(v))
        inverse = np.argsort(order)
        v, f = v[order], inverse[f[rng.permutation(len(f))]]
        r = self.check((v, f))
        self.assertTrue(r['completed'])
        self.assertEqual(r['squared_distance_exact'], '1/16')
        triangles = [tuple(tuple(Fraction(float(x)) for x in p) for p in v[f[i]])
                     for i in r['witness_faces']]
        self.assertEqual(triangle_distance(*triangles), Fraction(1, 16))

    def test_zero_one_and_island_in_one_void_are_not_applicable(self):
        meshes = [None, box(),
                  combine(box((0, 0, 0), (8, 8, 8)),
                          box((1, 1, 1), (7, 7, 7), True)),
                  combine(box((0, 0, 0), (8, 8, 8)),
                          box((1, 1, 1), (7, 7, 7), True),
                          box((2, 2, 2), (6, 6, 6)))]
        for mesh, count in zip(meshes, [0, 0, 1, 1]):
            with self.subTest(count=count, empty=mesh is None):
                r = self.check(mesh)
                self.assertTrue(r['completed'])
                self.assertEqual(r['state'], 'NOT_APPLICABLE')
                self.assertEqual(r['cavity_count'], count)
                for key in ('at_least_minimum', 'squared_distance_exact',
                            'witness_faces', 'witness_cavities'):
                    self.assertIsNone(r[key])

    def test_budget_exhaustion_cannot_report_partial_minimum(self):
        mesh = combine(box((-1, -1, -1), (7, 2, 2)),
                       box((0, 0, 0), (1, 1, 1), True),
                       box((2, 0, 0), (3, 1, 1), True),
                       box((3.25, 0, 0), (4.25, 1, 1), True))
        r = self.check(mesh, .5, maximum_pairs=2)
        self.assertFalse(r['completed'])
        self.assertEqual(r['state'], 'RESOURCE_LIMIT')
        self.assertEqual((r['pairs_required'], r['pairs_checked']), (3, 0))
        self.assertIsNone(r['at_least_minimum'])
        self.assertIsNone(r['squared_distance_exact'])
        self.assertIsNone(r['witness_faces'])
        full = self.check(mesh, .5, maximum_pairs=3)
        self.assertTrue(full['completed'])
        self.assertEqual(full['pairs_checked'], 3)
        self.assertEqual(full['squared_distance_exact'], '1/16')
        self.assertEqual(full['witness_cavities'], [2, 3])

    def test_input_and_threshold_change_identity_and_no_mutation(self):
        v, f = self.two_voids()
        before_v, before_f = v.copy(), f.copy()
        original = self.check((v, f))
        changed_threshold = self.check((v, f), .5)
        self.assertNotEqual(original['audit_digest'], changed_threshold['audit_digest'])
        np.testing.assert_array_equal(v, before_v)
        np.testing.assert_array_equal(f, before_f)
        v[16:, 0] += .25
        changed = self.check((v, f))
        self.assertEqual(changed['squared_distance_exact'], '1/4')
        self.assertNotEqual(original['audit_digest'], changed['audit_digest'])

    def test_invalid_input_threshold_and_budget_rejected(self):
        for threshold in (-1, float('nan'), float('inf'), True, '0.1'):
            with self.subTest(threshold=threshold), self.assertRaises((ValueError, TypeError)):
                self.check(self.two_voids(), threshold)
        for cap in (-1, 0, True, 1.5):
            with self.subTest(cap=cap), self.assertRaises((ValueError, TypeError)):
                self.check(self.two_voids(), maximum_pairs=cap)
        v, f = box()
        for mesh in ((v, f[:, ::-1]), (v, f[:-1])):
            with self.assertRaises(ValueError):
                self.check(mesh)


if __name__ == '__main__':
    unittest.main()
