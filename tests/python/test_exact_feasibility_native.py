"""Real backend integration. A missing backend is an error, never a skip."""
import importlib
import math
import unittest
from unittest.mock import patch
from types import SimpleNamespace

import _make_it_stand as native
import numpy as np
from test_reachability_bindings import cube


class ExactNativeFeasibilityTests(unittest.TestCase):
    def setUp(self):
        self.assertTrue(hasattr(native, 'CarvingReachabilityInput'), 'native typed snapshot missing')
        self.f = importlib.import_module('make_it_stand.verification.feasibility')

    def typed(self, x=3., fixed=2.):
        return native.CarvingReachabilityInput(1., [fixed, 0.],
            [native.ReachabilityCell(7, 1., [x, 0.])],
            [[-1., -1.], [1., -1.], [1., 1.], [-1., 1.]], 1e-9)

    def model(self, degrees=5.):
        config = native.BalanceModelConfig()
        config.voxel_resolution = 8
        model = native.BalanceModel(*cube(), config)
        model.add_standing_constraint_polygon([[0, 0, 0], [1, 0, 0], [0, 0, 1]], [0, -1, 0], math.radians(degrees))
        return model

    def test_real_typed_native_candidate_and_exact_certificate(self):
        data = self.typed()
        result = native.analyze_carving_reachability_input(data)
        self.assertEqual(result.status, 'PROVEN_UNREACHABLE')
        report = self.f.check_native_feasibility(data)
        self.assertEqual(report.status, 'INFEASIBLE_CERTIFIED')
        self.assertGreater(report.certificate.margin, 0)
        self.assertEqual(report.native_input_hash, result.input_hash)
        self.assertEqual(report.native_input_hash, data.input_hash)
        self.assertRegex(report.native_identity.sha256, r'^[a-f0-9]{64}$')
        self.assertRegex(report.checker_sha256, r'^[a-f0-9]{64}$')
        self.assertEqual(report.build_source_attestation, 'NOT_VERIFIED')
        self.assertEqual(report.to_dict(), self.f.check_native_feasibility(data).to_dict())

    def test_real_typed_boundary_native_result_cannot_grant_geometry(self):
        report = self.f.check_native_feasibility(self.typed(0.))
        self.assertEqual(report.status, 'UNKNOWN')
        self.assertTrue(report.witness.valid)
        self.assertEqual(report.witness.center, (1, 0))

    def test_model_snapshot_same_as_analyzed_and_is_read_only(self):
        model = self.model()
        snapshot = model.carving_reachability_snapshot(0)
        self.assertEqual(snapshot.preparation_status, 'READY')
        self.assertEqual(snapshot.input_hash, model.analyze_carving_reachability(0).input_hash)
        self.assertGreater(len(snapshot.input.variable_cells), 0)
        with self.assertRaises(AttributeError):
            snapshot.input.fixed_mass = 0.
        polygon = snapshot.input.safe_polygon_2d
        polygon[0, 0] += 5
        self.assertNotEqual(polygon[0, 0], snapshot.input.safe_polygon_2d[0, 0])
        report = self.f.check_native_feasibility(model)
        self.assertEqual(report.native_input_hash, snapshot.input_hash)
        self.assertNotEqual(report.status, 'FEASIBLE_WITNESS')
        self.assertEqual(report.native_status, 'RELAXATION_REACHABLE')
        self.assertEqual(report.native_pre_input_hash, report.native_post_input_hash)

    def test_degenerate_model_is_unknown_with_actual_reason(self):
        report = self.f.check_native_feasibility(self.model(60.))
        self.assertEqual(report.status, 'UNKNOWN')
        self.assertEqual(report.native_status, 'INDETERMINATE')
        self.assertIn('DEGENERATE_OR_UNATTAINABLE_SAFE_REGION', report.reason)

    def test_stale_native_result_and_mutation_fail_closed(self):
        model = self.model()
        actual = native.BalanceModel.analyze_carving_reachability
        def mutate(value, index=0):
            result = actual(value, index)
            value.clear_constraints()
            value.add_standing_constraint_polygon([[0, 0, 0], [1, 0, 0], [1, 0, 1], [0, 0, 1]], [0, -1, 0], 0.)
            return result
        with patch.object(native.BalanceModel, 'analyze_carving_reachability', mutate):
            report = self.f.check_native_feasibility(model)
        self.assertEqual(report.status, 'UNKNOWN')
        self.assertIn('INPUT_IDENTITY_MISMATCH', report.reason)
        self.assertNotEqual(report.native_pre_input_hash, report.native_post_input_hash)
        with patch.object(native, 'analyze_carving_reachability_input', lambda data: native_result):
            native_result = actual(self.model(), 0)
            report = self.f.check_native_feasibility(self.typed())
        self.assertEqual(report.status, 'UNKNOWN')
        self.assertIn('INPUT_IDENTITY_MISMATCH', report.reason)

    def test_no_old_json_or_flags_can_grant_native_proof(self):
        for old in ({'status': 'PROVEN_UNREACHABLE', 'certificate': {'valid': True}}, None, True):
            report = self.f.check_native_feasibility(old)
            self.assertEqual(report.status, 'UNKNOWN')
            self.assertIn('INVALID_NATIVE_SUBJECT', report.reason)

    def test_native_tolerance_claim_requires_strict_margin(self):
        # HiGHS accepts this subtolerance outside point as already safe.
        data = self.typed(3., fixed=math.nextafter(1., 2.))
        report = self.f.check_native_feasibility(data)
        self.assertNotEqual(report.status, 'FEASIBLE_WITNESS')
        self.assertFalse(report.witness.valid)

    def test_empty_certificate_solver_statistics_cannot_certify(self):
        # Corrupt only the observed status on a real freshly solved native
        # result. The candidate itself is not a mock proof generator.
        with patch.object(native.CarvingReachabilityResult, 'certificate_solver_status', ''):
            report = self.f.check_native_feasibility(self.typed())
        self.assertEqual(report.status, 'UNKNOWN')
        self.assertIn('empty native certificate solver status', report.reason)

    def test_legacy_true_flag_and_numeric_failure_do_not_grant_certificate(self):
        empty = SimpleNamespace(weights=np.array([]), valid=True,
                                minimum_weight=0., normalization_residual=0.,
                                stationarity_infinity_norm=0., contradiction=-1.)
        for status, solver in [('PROVEN_UNREACHABLE', 'Infeasible'), ('INDETERMINATE', 'API_ERROR')]:
            with patch.object(native.CarvingReachabilityResult, 'certificate', empty), \
                 patch.object(native.CarvingReachabilityResult, 'status', status), \
                 patch.object(native.CarvingReachabilityResult, 'primal_solver_status', solver):
                report = self.f.check_native_feasibility(self.typed())
            self.assertEqual(report.status, 'UNKNOWN')
            self.assertEqual(report.native_status, status)
            self.assertEqual(report.native_primal_status, solver)
            self.assertEqual(report.certificate.reason, 'MISSING_CANDIDATE')

    def test_malformed_native_statistics_fail_closed(self):
        for attribute, value in [('fixed_mass', float('nan')), ('final_mass', float('inf')),
                                 ('input_cell_count', True), ('primal_solver_status', ''),
                                 ('certificate_system_nonzero_count', -1)]:
            with self.subTest(attribute=attribute), patch.object(native.CarvingReachabilityResult, attribute, value):
                report = self.f.check_native_feasibility(self.typed())
                self.assertEqual(report.status, 'UNKNOWN')
                self.assertIn('NATIVE_CHECK_FAILURE', report.reason)

    def test_missing_backend_is_explicit_failure(self):
        with patch.dict('sys.modules', {'_make_it_stand': None}):
            report = self.f.check_native_feasibility(self.typed())
        self.assertEqual(report.status, 'UNKNOWN')
        self.assertIn('MISSING_NATIVE_BACKEND', report.reason)

    def test_real_native_no_variable_and_oblique_inputs(self):
        square = [[-1., -1.], [1., -1.], [1., 1.], [-1., 1.]]
        for x in (0., 2.):
            data = native.CarvingReachabilityInput(1., [x, 0.], [], square)
            report = self.f.check_native_feasibility(data)
            self.assertEqual(report.witness.valid, x == 0.)
            self.assertNotEqual(report.status, 'FEASIBLE_WITNESS')
            self.assertEqual(report.native_input_hash, data.input_hash)
        data = native.CarvingReachabilityInput(1., [1., 2.],
            [native.ReachabilityCell(7, 2., [3., 5.])], [[0., 0.], [2., 0.], [0., 2.]])
        report = self.f.check_native_feasibility(data)
        self.assertEqual(report.status, 'INFEASIBLE_CERTIFIED')
        self.assertGreater(report.certificate.margin, 0)


if __name__ == '__main__':
    unittest.main()
