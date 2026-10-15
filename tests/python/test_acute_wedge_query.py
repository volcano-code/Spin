"""Exact moving-center counterexamples; never a general wall-thickness pass."""
from fractions import Fraction as Q
import importlib
import unittest
import numpy as np

from experiments.teapot_spin2014.segment_material_query import ClosedBoundaryQueries
from experiments.teapot_spin2014.boundary_cap_query import BoundaryCapQuery
from experiments.teapot_spin2014.near_feature_witness import add, mul
from test_boundary_cap_query import box


def prism():
    vertices = np.array([[0, 0, 0], [1, 0, 0], [1, .5, 0],
                         [0, 0, 1], [1, 0, 1], [1, .5, 1]], dtype=float)
    faces = np.array([[0, 2, 1], [3, 4, 5], [0, 1, 4], [0, 4, 3],
                      [1, 2, 5], [1, 5, 4], [2, 0, 3], [2, 3, 5]])
    return vertices, faces


class AcuteWedgeTests(unittest.TestCase):
    def setUp(self):
        name = 'experiments.teapot_spin2014.acute_wedge_query'
        self.assertIsNotNone(importlib.util.find_spec(name), 'moving-center acute wedge checker missing')
        self.Query = importlib.import_module(name).AcuteWedgeQuery

    def query(self, mesh):
        return self.Query(BoundaryCapQuery(ClosedBoundaryQueries(*mesh)))

    def test_material_ball_cut_has_continuum_witness_not_only_vertex_sample(self):
        query = self.query(prism())
        self.assertIn([0, 3], query.convex_acute_edges())
        row = query.certificate([0, 3], Q(1, 10000), maximum_halvings=100)
        self.assertEqual(row['state'], 'CONDITIONAL_MATERIAL_BALL_CUT_INFIMUM_ZERO')
        self.assertEqual(row['operator'], 'MATERIAL_MINUS_OPEN_BALL')
        self.assertEqual(row['material_ball_cut_radius_infimum_exact'], '0')
        self.assertEqual(row['complete_wall_thickness'], 'UNKNOWN')
        self.assertFalse(row['wall_metric_adopted'])
        self.assertEqual(row['wall_thickness_admission'], 'NOT_GRANTED')
        self.assertEqual(row['local_cut_euler_change'], -1)
        maximum = Q(row['maximum_scale_exact'])
        midpoint = tuple(map(Q, row['edge_midpoint_exact']))
        direction = tuple(map(Q, row['source_direction_exact']))
        target_direction = tuple(map(Q, row['target_direction_exact']))
        coefficient = Q(row['radius_squared_coefficient_exact'])
        # Independent full-boundary filtration/segment checks at several scales.
        # The certificate's inequalities, not this finite list, cover ALL scales.
        for factor in (1, 2, 7, 100):
            t = maximum / factor
            center = add(midpoint, mul(direction, t))
            target = add(midpoint, mul(target_direction, t))
            radius2 = t*t*coefficient
            cap = query.caps.audit(center, radius2)
            self.assertFalse(cap['threshold_cap_single_disk'])
            self.assertEqual([c['euler'] for c in cap['critical_intervals'][-1]['components']], [1, 1])
            self.assertEqual(query.caps.boundary.audit_segment(center, target)['state'],
                             'OPEN_SEGMENT_PARITY_INSIDE')
            self.assertLess(radius2, Q(1, 10000))

    def test_right_angle_and_coplanar_edges_are_not_acute(self):
        query = self.Query(box())
        self.assertEqual(query.convex_acute_edges(), [])
        r = query.certificate([0, 1], Q(1, 100), maximum_halvings=100)
        self.assertEqual(r['state'], 'NOT_CONVEX_ACUTE')
        self.assertEqual(r['complete_wall_thickness'], 'UNKNOWN')

    def test_inward_cavity_corner_is_not_positive_material_wedge(self):
        v, f = prism()
        query = self.query((v, f[:, ::-1]))
        self.assertEqual(query.convex_acute_edges(), [])
        self.assertEqual(query.certificate([0, 3], Q(1, 100), maximum_halvings=100)['state'],
                         'NOT_CONVEX_ACUTE')

    def test_halving_limit_cannot_publish_a_continuum_certificate(self):
        row = self.query(prism()).certificate([0, 3], Q(1, 10**40), maximum_halvings=1)
        self.assertEqual(row['state'], 'RESOURCE_LIMIT')
        self.assertNotIn('maximum_scale_exact', row)
        self.assertEqual(row['wall_thickness_admission'], 'NOT_GRANTED')

    def test_source_triangle_subdivision_preserves_acute_edge(self):
        v, f = prism()
        # Dyadic interior splits preserve the exact original triangle plane.
        vertices = v.tolist()
        faces = []
        for a, b, c in f:
            p = len(vertices)
            vertices.append((v[a]/2 + v[b]/4 + v[c]/4).tolist())
            faces.extend([[int(a), int(b), p], [int(b), int(c), p], [int(c), int(a), p]])
        query = self.query((np.array(vertices), np.array(faces)))
        self.assertIn([0, 3], query.convex_acute_edges())
        self.assertEqual(query.certificate([0, 3], Q(1, 10000), maximum_halvings=100)['state'],
                         'CONDITIONAL_MATERIAL_BALL_CUT_INFIMUM_ZERO')

    def test_invalid_edge_budget_and_radius_rejected(self):
        query = self.query(prism())
        for edge in ([0, 0], [0, 5], [True, 3], [0, 3, 4]):
            with self.assertRaises(ValueError):
                query.certificate(edge, Q(1), maximum_halvings=100)
        for radius in (0, -1, True, float('inf'), float('nan')):
            with self.assertRaises((ValueError, OverflowError)):
                query.certificate([0, 3], radius, maximum_halvings=100)
        for cap in (0, -1, True, 1.5):
            with self.assertRaises(ValueError):
                query.certificate([0, 3], Q(1), maximum_halvings=cap)


if __name__ == '__main__':
    unittest.main()
