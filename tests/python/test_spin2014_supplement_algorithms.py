"""Tests for an independent transcription of Spin-It 2014 Supplement Algorithms 1--2.

These tests deliberately do not call the production C++ kernel while evaluating
the reference.  The comparison to the extension is only the final test oracle.
"""

from __future__ import annotations

import importlib
import unittest

import numpy as np

from test_independent_mass_moments import TETRAHEDRON, box_mesh


class Spin2014SupplementAlgorithmTests(unittest.TestCase):
    def setUp(self) -> None:
        self.reference = importlib.import_module(
            "make_it_stand.verification.spin2014_supplement_algorithms"
        )

    def test_algorithm_one_hand_derived_asymmetric_tetrahedron(self) -> None:
        values = self.reference.integrate_surface_moments(
            TETRAHEDRON.vertices, TETRAHEDRON.faces, density=1.0
        )
        np.testing.assert_allclose(
            values,
            [1.0, 0.25, 0.5, 0.75, 0.1, 0.4, 0.9, 0.3, 0.15, 0.1],
            rtol=0.0,
            atol=1e-14,
        )
        self.assertEqual(tuple(self.reference.ORDER),
                         ("M", "LX", "LY", "LZ", "QXX", "QYY", "QZZ", "QYZ", "QXZ", "QXY"))

    def test_algorithm_two_matches_centered_finite_differences(self) -> None:
        vertices = np.asarray(TETRAHEDRON.vertices, dtype=float)
        analytic = self.reference.differentiate_surface_moments(
            vertices, TETRAHEDRON.faces, density=2.7
        )
        self.assertEqual(analytic.shape, (10, 12))
        h = 1e-6
        finite = np.zeros_like(analytic)
        for vertex in range(len(vertices)):
            for axis in range(3):
                plus, minus = vertices.copy(), vertices.copy()
                plus[vertex, axis] += h
                minus[vertex, axis] -= h
                finite[:, 3 * vertex + axis] = (
                    self.reference.integrate_surface_moments(plus, TETRAHEDRON.faces, 2.7)
                    - self.reference.integrate_surface_moments(minus, TETRAHEDRON.faces, 2.7)
                ) / (2 * h)
        np.testing.assert_allclose(analytic, finite, rtol=2e-8, atol=2e-9)

    def test_reference_and_native_agree_without_shared_calls(self) -> None:
        native_module = importlib.import_module("_make_it_stand")
        mesh = box_mesh((-2.0, -1.0, -3.0), (1.5, 2.0, 0.75))
        reference_moments = self.reference.integrate_surface_moments(
            mesh.vertices, mesh.faces, density=1.7
        )
        reference_jacobian = self.reference.differentiate_surface_moments(
            mesh.vertices, mesh.faces, density=1.7
        )
        native = native_module.differentiate_surface_raw_moments(
            mesh.vertices, mesh.faces, 1.7
        )
        np.testing.assert_allclose(reference_moments, native["moments"], rtol=0.0, atol=2e-13)
        np.testing.assert_allclose(reference_jacobian, native["jacobian"], rtol=0.0, atol=4e-13)


if __name__ == "__main__":
    unittest.main()
