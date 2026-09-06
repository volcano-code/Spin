"""Hand-derived represented-LP checks, independent of a native solver."""
from dataclasses import FrozenInstanceError, replace
from fractions import Fraction as Q
import importlib
import importlib.util
import math
import unittest


class ExactFeasibilityTests(unittest.TestCase):
    def setUp(self):
        name = "make_it_stand.verification.feasibility"
        self.assertIsNotNone(importlib.util.find_spec(name), "exact checker is missing")
        self.f = importlib.import_module(name)

    def snapshot(self, x=3.0, *, cells=True, fixed=2.0, polygon=None):
        return self.f.MomentSnapshot(
            1.0, (fixed, 0.0),
            (self.f.MomentCell(7, 1.0, (x, 0.0)),) if cells else (),
            tuple(polygon or ((-1., -1.), (1., -1.), (1., 1.), (-1., 1.))),
            1e-9,
        )

    def test_right_edge_impossible_and_scaled_certificate(self):
        for weights, margin, column_sum in [
            ((0, 1, 0, 0), Q(2), Q(4)),
            ((0, 3, 0, 0), Q(6), Q(12)),
            ((0, .1, 0, 0), Q(3602879701896397, 18014398509481984),
             Q(3602879701896397, 9007199254740992)),
        ]:
            report = self.f.verify_snapshot(self.snapshot(), weights=weights)
            self.assertEqual(report.status, "INFEASIBLE_CERTIFIED")
            self.assertEqual(report.scope, "REPRESENTED_DISCRETE_MOMENT_LP")
            self.assertEqual(report.certificate.margin, margin)
            self.assertEqual(report.certificate.column_sums, (column_sum,))
            self.assertEqual(report.certificate.lower_bound, 0)
            self.assertEqual(report.certificate.weighted_rhs, -margin)

    def test_boundary_witness_never_promotes_geometry(self):
        report = self.f.verify_snapshot(self.snapshot(0.), weights=(0, 1, 0, 0), retained_fractions=(1.,))
        self.assertEqual(report.status, "UNKNOWN")
        self.assertFalse(report.certificate.valid)
        self.assertEqual(report.certificate.margin, 0)
        self.assertTrue(report.witness.valid)
        self.assertEqual(report.witness.mass, 2)
        self.assertEqual(report.witness.center, (Q(1), Q(0)))
        self.assertEqual(report.geometry_validation, "NOT_VERIFIED")
        self.assertEqual(report.continuous_geometry, "UNKNOWN")
        self.assertEqual(report.physical_validation, "NOT_IN_SCOPE")
        self.assertIn("FEASIBLE_WITNESS", self.f.STATUSES)
        with self.assertRaises(FrozenInstanceError):
            report.status = "FEASIBLE_WITNESS"

    def test_no_cells_outside_and_inside(self):
        outside = self.f.verify_snapshot(self.snapshot(cells=False), weights=(0, 1, 0, 0), retained_fractions=())
        inside = self.f.verify_snapshot(self.snapshot(cells=False, fixed=0.), weights=(0, 1, 0, 0), retained_fractions=())
        self.assertEqual(outside.status, "INFEASIBLE_CERTIFIED")
        self.assertFalse(outside.witness.valid)
        self.assertEqual(inside.status, "UNKNOWN")
        self.assertTrue(inside.witness.valid)

    def test_oblique_rows_include_xy_moments(self):
        snapshot = replace(self.snapshot(polygon=((0., 0.), (2., 0.), (0., 2.))),
                           fixed_first_moment=(1., 2.),
                           cells=(self.f.MomentCell(7, 2., (3., 5.)),))
        system = self.f.build_system(snapshot)
        self.assertEqual(system.C, ((Q(-10),), (Q(8),), (Q(-6),)))
        self.assertEqual(system.d, (Q(4), Q(-2), Q(2)))

    def test_binary_decimal_and_subtolerance_margin(self):
        snapshot = replace(self.snapshot(cells=False, fixed=math.nextafter(.1, math.inf)),
                           polygon=((-1., -1.), (.1, -1.), (.1, 1.), (-1., 1.)))
        report = self.f.verify_snapshot(snapshot, weights=(0, 1., 0, 0))
        self.assertEqual(report.certificate.margin, Q(1, 36028797018963968))
        self.assertEqual(report.status, "INFEASIBLE_CERTIFIED")
        self.assertNotEqual(Q(.1), Q(1, 10))

    def test_candidate_rejection_and_invalid_witness(self):
        for weights in [(), (0, 1), (0, -1, 0, 0), (0, -5e-324, 0, 0), (0, 0, 0, 0),
                        (0, float('nan'), 0, 0), (0, float('inf'), 0, 0),
                        (0, True, 0, 0), 'bad', None]:
            with self.subTest(weights=weights):
                report = self.f.verify_snapshot(self.snapshot(), weights=weights)
                self.assertEqual(report.status, "UNKNOWN")
                self.assertFalse(report.certificate.valid)
        for retained in [(), (-.1,), (1.1,), (True,), (float('nan'),), ('1',)]:
            self.assertFalse(self.f.verify_snapshot(self.snapshot(0.), retained_fractions=retained).witness.valid)

    def test_exact_validation_rejects_malformed_input(self):
        snapshot = self.snapshot()
        bad = [replace(snapshot, fixed_mass=x) for x in (True, float('nan'), float('inf'), 0, -1, '1')]
        bad += [replace(snapshot, tolerance=x) for x in (True, 0, float('nan'))]
        bad += [replace(snapshot, polygon=p) for p in (
            (), ((0, 0), (1, 0)), ((0, 0), (1, 0), (2, 0)),
            ((0, 0), (1, 1), (0, 1), (1, 0)),
            ((0, 0), (3, 3), (0, 2), (2, 0)),
            ((0, 0), (1, 0), (1, 1), (0, 0)),
            ((0, 0), (0, 1), (1, 1), (1, 0)),
            ((0, 0), (1, 0), (float('inf'), 1)),
        )]
        bad += [replace(snapshot, cells=(c,)) for c in (
            self.f.MomentCell(-1, 1., (0., 0.)), self.f.MomentCell(True, 1., (0., 0.)),
            self.f.MomentCell(0, 0., (1., 0.)), self.f.MomentCell(0, -1., (0., 0.)),
            self.f.MomentCell(0, 1., (float('nan'), 0.)),
        )]
        bad.append(replace(snapshot, cells=snapshot.cells * 2))
        bad += [replace(snapshot, fixed_first_moment=()), {}, None, True]
        for item in bad:
            with self.subTest(item=item):
                report = self.f.verify_snapshot(item, weights=(0, 1, 0, 0))
                self.assertEqual(report.status, "UNKNOWN")
                self.assertIn("INVALID_INPUT", report.reason)

    def test_collinear_redundant_vertices_are_validated(self):
        snapshot = self.snapshot(polygon=((-1, -1), (0, -1), (1, -1), (1, 1), (-1, 1)))
        self.assertEqual(self.f.verify_snapshot(snapshot, weights=(0, 0, 1, 0, 0)).status, "INFEASIBLE_CERTIFIED")
        zero_cell = replace(self.snapshot(), cells=(self.f.MomentCell(7, 0., (0., 0.)),))
        self.assertEqual(self.f.verify_snapshot(zero_cell, weights=(0, 1, 0, 0)).status, "INFEASIBLE_CERTIFIED")

    def test_hash_binds_cell_order_ids_target_and_not_caller_hash(self):
        base = replace(self.snapshot(), cells=(self.f.MomentCell(7, 1., (3., 0.)), self.f.MomentCell(8, 2., (4., 0.))))
        first = self.f.verify_snapshot(base, weights=(0, 1, 0, 0))
        self.assertEqual(first.to_dict(), self.f.verify_snapshot(base, weights=(0, 1, 0, 0)).to_dict())
        for mutated in (replace(base, cells=base.cells[::-1]),
                        replace(base, cells=(replace(base.cells[0], cell_id=9), base.cells[1])),
                        replace(base, polygon=((-1, -1), (2, -1), (2, 1), (-1, 1)))):
            report = self.f.verify_snapshot(mutated, weights=(0, 1, 0, 0), expected_input_hash=first.input_hash)
            self.assertEqual(report.status, "UNKNOWN")
            self.assertIn("INPUT_IDENTITY_MISMATCH", report.reason)
            self.assertNotEqual(first.input_hash, report.input_hash)

    def test_primal_and_certificate_cannot_both_validate(self):
        for x in (0., math.nextafter(0., 1.), 1., 3.):
            for fixed in (1., 2., math.nextafter(2., 3.)):
                for retained in (0., .5, 1.):
                    report = self.f.verify_snapshot(self.snapshot(x, fixed=fixed), weights=(0, 1, 0, 0), retained_fractions=(retained,))
                    self.assertFalse(report.witness.valid and report.certificate.valid)


if __name__ == '__main__':
    unittest.main()
