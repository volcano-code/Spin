"""Quasi-static fixed-support table/gravity tilt oracle tests."""

from __future__ import annotations

import unittest

from make_it_stand.stability_oracle import (
    evaluate_quasi_static_tilt,
    evaluate_tilt_azimuths,
)


class StabilityOracleTests(unittest.TestCase):
    def setUp(self) -> None:
        self.support = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 0.0, 1.0],
            [0.0, 0.0, 1.0],
        ]
        self.center = [0.5, 1.0, 0.5]
        self.gravity = [0.0, -1.0, 0.0]

    def test_square_critical_tilt_matches_analytic_angle(self) -> None:
        nominal = evaluate_quasi_static_tilt(
            self.support,
            self.center,
            self.gravity,
            tilt_degrees=0.0,
            azimuth_degrees=0.0,
        )
        self.assertTrue(nominal.inside_support)
        self.assertAlmostEqual(nominal.signed_margin, 0.5, places=12)
        self.assertAlmostEqual(nominal.critical_tilt_degrees, 26.565051177, places=8)

        self.assertTrue(evaluate_quasi_static_tilt(
            self.support,
            self.center,
            self.gravity,
            tilt_degrees=26.5,
            azimuth_degrees=0.0,
        ).inside_support)
        self.assertFalse(evaluate_quasi_static_tilt(
            self.support,
            self.center,
            self.gravity,
            tilt_degrees=26.6,
            azimuth_degrees=0.0,
        ).inside_support)

    def test_eight_two_degree_directions_preserve_fixed_support(self) -> None:
        results = evaluate_tilt_azimuths(
            self.support,
            self.center,
            self.gravity,
            tilt_degrees=2.0,
            azimuths_degrees=range(0, 360, 45),
        )
        self.assertEqual(len(results), 8)
        self.assertTrue(all(item.inside_support for item in results))
        self.assertTrue(all(item.protocol == "quasi_static_table_tilt" for item in results))

    def test_nonplanar_support_heights_are_projected_not_reseated(self) -> None:
        support = [point[:] for point in self.support]
        support[1][1] = 0.02
        result = evaluate_quasi_static_tilt(
            support,
            self.center,
            self.gravity,
            tilt_degrees=2.0,
            azimuth_degrees=90.0,
        )
        self.assertTrue(result.inside_support)
        self.assertEqual(result.contact_policy, "fixed_projected_support_polygon")


if __name__ == "__main__":
    unittest.main()
