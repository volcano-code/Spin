import unittest
from noodles_vertex_link_witness import link_components


class VertexLinkTests(unittest.TestCase):
    def test_closed_single_star(self):
        result = link_components([[0, 1, 2], [0, 2, 3], [0, 3, 1]], 0)
        self.assertTrue(result['is_single_cycle'])

    def test_two_touching_closed_stars_are_not_manifold(self):
        result = link_components([[0, 1, 2], [0, 2, 3], [0, 3, 1],
                                  [0, 4, 5], [0, 5, 6], [0, 6, 4]], 0)
        self.assertFalse(result['is_single_cycle'])
        self.assertEqual(len(result['components']), 2)

    def test_open_star_is_rejected(self):
        self.assertFalse(link_components([[0, 1, 2]], 0)['is_single_cycle'])
