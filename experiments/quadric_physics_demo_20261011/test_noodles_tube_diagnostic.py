import unittest
from fractions import Fraction as Q
from diagnose_noodles_tube import segment_triangle


class ExactTubeWitnessTests(unittest.TestCase):
    def test_closed_v4_tube003_still_has_exact_self_intersection(self):
        from pathlib import Path
        from diagnose_noodles_tube import diagnose
        source=Path(__file__).parent/'noodles-dense-v4/tube-003-fixed-model.obj'
        record=diagnose(source)
        self.assertEqual(record['status'],'SELF_INTERSECTION_WITNESSED')
        pairs={tuple(w['faces']) for w in record['strict_exact_edge_triangle_witnesses']}
        self.assertIn((456,504),pairs)
        self.assertFalse(record['complete_admission'])

    def test_strict_crossing_has_exact_barycentric_witness(self):
        a,b,c=(Q(0),Q(0),Q(0)),(Q(1),Q(0),Q(0)),(Q(0),Q(1),Q(0))
        p,q=(Q(1,4),Q(1,4),Q(-1)),(Q(1,4),Q(1,4),Q(1))
        hit=segment_triangle(p,q,a,b,c)
        self.assertEqual(hit['point'],['1/4','1/4','0'])
        self.assertEqual(hit['segment_parameter'],'1/2')

    def test_touch_or_outside_is_not_positive_crossing(self):
        a,b,c=(Q(0),Q(0),Q(0)),(Q(1),Q(0),Q(0)),(Q(0),Q(1),Q(0))
        self.assertIsNone(segment_triangle((Q(0),Q(0),Q(-1)),(Q(0),Q(0),Q(1)),a,b,c))
        self.assertIsNone(segment_triangle((Q(2),Q(2),Q(-1)),(Q(2),Q(2),Q(1)),a,b,c))


if __name__=='__main__':unittest.main()
