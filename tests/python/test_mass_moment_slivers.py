"""Nonzero source triangles must not disappear during checker normalization."""
import unittest
import numpy as np
from make_it_stand.verification.mass_moments import (
    TriangleMesh, GeometryValidationError, compute_mass_properties,
)


def tetra(split=False, epsilon=1e-20):
    vertices = [(0., 0., 0.), (1., 0., 0.), (0., 2., 0.), (0., 0., 3.)]
    faces = [(0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3)]
    if split:
        vertices.append((0., epsilon, 0.))
        faces = [(0, 4, 1), (4, 2, 1), (0, 1, 3), (0, 3, 4), (4, 3, 2), (1, 2, 3)]
    return TriangleMesh(tuple(vertices), tuple(faces))


class SliverTests(unittest.TestCase):
    def test_nonzero_edge_subdivision_preserves_all_moments(self):
        expected = compute_mass_properties(tetra())
        for epsilon in [1e-18, 1e-20, 1e-100]:
            actual = compute_mass_properties(tetra(True, epsilon))
            self.assertAlmostEqual(actual.mass, expected.mass, places=14)
            for name in ['first_moment', 'second_moment', 'center_of_mass', 'central_inertia']:
                np.testing.assert_allclose(getattr(actual, name), getattr(expected, name), atol=1e-14, rtol=1e-14)

    def test_truly_degenerate_subdivision_still_rejected(self):
        with self.assertRaisesRegex(GeometryValidationError, 'degenerate'):
            compute_mass_properties(tetra(True, 0.))

    def test_artifact_affine_frame_does_not_round_source_vertices_first(self):
        from make_it_stand.verification.spin_artifact import _integrate
        class Mesh:
            def __init__(self, value): self.value = value
            def arrays(self): return np.array(self.value.vertices), np.array(self.value.faces)
        frame = dict(origin=[.5, 1., 1.5], L=3., R=2.)
        geo, physical = _integrate(Mesh(tetra(True)), 4., frame)
        expected_geo, expected = _integrate(Mesh(tetra()), 4., frame)
        self.assertAlmostEqual(geo['volume'], expected_geo['volume'], places=14)
        for key in ['moments', 'center_of_mass', 'central_inertia']:
            np.testing.assert_allclose(physical[key], expected[key], atol=1e-14, rtol=1e-14)


if __name__ == '__main__': unittest.main()
