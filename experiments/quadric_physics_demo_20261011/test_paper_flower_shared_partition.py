import unittest
from fractions import Fraction as F


class SharedArrangementTests(unittest.TestCase):
    def test_full_middle_source_uses_same_noded_boundary_as_design(self):
        import paper_flower_shared_partition as shared
        from shapely.geometry import Polygon
        source=Polygon([[0,0],[1,3],[0,3]])
        cutter=Polygon([[-1,1],[2,1],[2,2],[-1,2]])
        result=shared.shared_partition(source,source,cutter,.02)
        self.assertTrue(result['rim'].is_empty)
        self.assertEqual(result['original'].symmetric_difference(result['design']).area,0)
        self.assertTrue(result['certificate']['source_and_design_share_noded_boundary'])

    def test_rounded_intersection_is_shared_by_both_sides_not_original_edge(self):
        # A rational crossing generally cannot remain exactly on a stored
        # binary64 edge. Both sides must adopt the same derived boundary.
        import paper_flower_shared_partition as shared
        crossing=shared.segment_crossing((F(0),F(0)),(F(1),F(3)),
                                         (F(0),F(1)),(F(1),F(1)))
        self.assertEqual(crossing,(F(1,3),F(1)))
        stored=tuple(F(float(v)) for v in crossing)
        self.assertNotEqual(stored[1]-3*stored[0],0)
        from shapely.geometry import Polygon
        original=Polygon([[-1,-1],[2,-1],[2,4],[-1,4]])
        design=Polygon([[0,0],[1,3],[0,3]])
        cutter=Polygon([[-1,1],[2,1],[2,2],[-1,2]])
        result=shared.shared_partition(original,design,cutter,1.)
        self.assertEqual(result['design'].symmetric_difference(result['retained'].union(result['cavity'])).area,0)
        self.assertEqual(result['original'].symmetric_difference(result['design'].union(result['rim'])).area,0)
        self.assertEqual(result['retained'].intersection(result['rim']).area,0)

    def test_curve_bound_covers_entire_chord_not_samples(self):
        import paper_flower_shared_partition as shared
        bound=shared.circle_chord_bound([(F(1),F(0)),(F(0),F(1)),
                                        (F(-1),F(0)),(F(0),F(-1))],(F(0),F(0)),F(1))
        self.assertEqual(bound,F(1,2))

    def test_saved_bound_checker_rejects_changed_authority_circle(self):
        import paper_flower_saved_check as checker
        with self.assertRaisesRegex(ValueError,'CERTIFICATE_CIRCLE_CHANGED'):
            checker.verify_shared_certificate(dict(schema='shared-rational-arrangement-v2',
                circle_center_model=[0,0],circle_radius_model=1),{},{},dict(center=[0,0],radius=2))

    def test_source_edge_indices_use_authenticated_ring_order_not_manifest_start(self):
        import paper_flower_saved_check as checker
        frozen=[[[0.,0.],[1.,0.],[1.,1.],[0.,1.]]]
        # Same source boundary, cyclic start moved by GEOS; edge0 now means
        # manifest edge1. These are not interchangeable local indices.
        ordered=[[["1","0"],["1","1"],["0","1"],["0","0"]]]
        actual=checker.authenticate_source_loop_order(ordered,frozen)
        self.assertEqual(actual[0][0],(F(1),F(0)))
        damaged=[[list(p) for p in ordered[0]]];damaged[0][0]=["2","0"]
        with self.assertRaisesRegex(ValueError,'SOURCE_ORDER_BOUNDARY_CHANGED'):
            checker.authenticate_source_loop_order(damaged,frozen)


if __name__=='__main__':unittest.main()
