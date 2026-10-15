"""Independent checks for Spin-It 2014 Eq. (1)'s strict stability threshold."""

from __future__ import annotations

import importlib
import unittest

from make_it_stand.verification.spin2014_givens_reference import stability_threshold_squared


class Spin2014StabilityThresholdTests(unittest.TestCase):
    def setUp(self) -> None:
        self.native = importlib.import_module("_make_it_stand")

    def test_native_matches_published_scalar_formula(self) -> None:
        expected = stability_threshold_squared(4.0, 9.81, 0.5, 10.0, 7.0)
        actual = self.native.spin_it_2014_stability_threshold_squared(4.0, 9.81, 0.5, 10.0, 7.0)
        self.assertAlmostEqual(actual, expected, places=14)

    def test_nonphysical_or_degenerate_inputs_are_rejected(self) -> None:
        for values in ((0.0, 9.81, 0.5, 10.0, 7.0),
                       (4.0, -9.81, 0.5, 10.0, 7.0),
                       (4.0, 9.81, -0.5, 10.0, 7.0),
                       (4.0, 9.81, 0.5, 7.0, 7.0),
                       (4.0, 9.81, 0.5, 6.0, 7.0)):
            with self.subTest(values=values), self.assertRaises(ValueError):
                self.native.spin_it_2014_stability_threshold_squared(*values)


if __name__ == "__main__":
    unittest.main()
