"""Strict full-domain KKT witness; no changed weights or success-selected starts."""
import copy
from pathlib import Path
import tempfile
import unittest

import numpy as np
from test_faster_static_table_dispatch import method_case
from test_spin_it_faster_table1_runner import load_runner
from make_it_stand.spin import solve_spin_it_faster_static_plane_family


class StaticFullDomainTests(unittest.TestCase):
    def test_full_state_cannot_use_small_interior_residual_as_optimality(self):
        runner = load_runner()
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); case=method_case(root,"LOWER_CENTER_OF_MASS")
            case["material_density"]=1e12
            result=solve_spin_it_faster_static_plane_family(runner.build_request(case,asset_root=root),
                initial_coefficients=(2.,.1,0.,-.5))
            self.assertLess(result["normalized_stationarity_residual"],1e-10)
            self.assertNotEqual(result["status"],"CONVERGED")

    def test_failed_full_state_without_witness_remains_auditable_not_qualified(self):
        from make_it_stand.verification.static_material import audit_static_case_material
        runner=load_runner()
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); case=method_case(root,"LOWER_CENTER_OF_MASS")
            actual=runner.run_case(case,1,1,asset_root=root)
            result=solve_spin_it_faster_static_plane_family(runner.build_request(case,asset_root=root),
                initial_coefficients=(2.,.1,0.,-.5))
            self.assertNotEqual(result["status"],"CONVERGED")
            actual["attempts"][0].update(status=result["status"],solver_result=runner._planar_json_value(result),
                                         function_evaluations=len(result["trace"]))
            checked=audit_static_case_material(actual)["attempts"][0]
            self.assertTrue(checked["moments_match"])
            self.assertFalse(checked["terminal_numerics_qualified"])

    def test_typed_full_record_cannot_strip_witness_and_use_legacy_interior_gate(self):
        from make_it_stand.verification.static_material import audit_static_case_material
        runner = load_runner()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            forged = runner.run_case(method_case(root,"PARETO_TRADEOFF",.5),1,2,asset_root=root)
            row = forged["attempts"][0]; result = row["solver_result"]
            plane = result["boundary_optimality"]["plane_coefficients"]
            result["plane_coefficients"] = plane
            result["trace"][-1]["plane_coefficients"] = plane
            result["boundary_optimality"] = None
            with self.subTest(checker="terminal"):
                self.assertFalse(runner.static_plane_attempt_acceptance(row)["passed"])
            with self.subTest(checker="saved_geometry"):
                with self.assertRaisesRegex(ValueError,"witness"):
                    audit_static_case_material(forged)

    def test_frozen_pareto_full_material_has_checked_boundary_optimality(self):
        runner = load_runner()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            case = method_case(root, "PARETO_TRADEOFF", .5)
            actual = runner.run_case(case, 2, 8, asset_root=root)
            self.assertTrue(actual["passed"], str(actual["summary"]))
            self.assertEqual(len(actual["attempts"]), 2)
            for row in actual["attempts"]:
                result = row["solver_result"]
                self.assertGreater(result["accepted_steps"], 0)
                witness = result.get("boundary_optimality")
                self.assertIsInstance(witness, dict)
                self.assertEqual(witness["kind"], "STRICT_FULL_DOMAIN_ZERO_EQUALITY_MULTIPLIERS_V1")
                self.assertGreater(witness["minimum_gradient_polynomial"], 0.)
                # Actual total mass, no arbitrary regularization or truncation.
                self.assertAlmostEqual(result["total_raw_moments"][0], 108.8, places=10)
                self.assertEqual(result["material_region"], "STRICT_FULL")
                self.assertFalse(np.allclose(result["plane_coefficients"], witness["plane_coefficients"]))

    def test_saved_witness_reintegrates_and_rejects_forged_or_missing_evidence(self):
        from make_it_stand.verification.static_material import audit_static_case_material
        runner = load_runner()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            actual = runner.run_case(method_case(root,"PARETO_TRADEOFF",.5),1,2,asset_root=root)
            checked = audit_static_case_material(actual)["attempts"][0]
            self.assertTrue(checked["terminal_numerics_qualified"])
            self.assertTrue(checked["full_domain_gradient_sign_verified"])
            self.assertGreater(checked["stationarity_residual"],1e-10)
            for field,value in (("minimum_gradient_polynomial", -1.),
                                ("plane_coefficients", [1.,0.,0.,0.]),
                                ("objective_measure", "GEOMETRIC_VOLUME")):
                forged = copy.deepcopy(actual)
                row = forged["attempts"][0]
                row["solver_result"]["boundary_optimality"][field] = value
                self.assertFalse(runner.static_plane_attempt_acceptance(row)["passed"])
                with self.assertRaises(ValueError): audit_static_case_material(forged)
            for field,value in (("material_region", "UNCLASSIFIED"), ("boundary_optimality",None)):
                forged = copy.deepcopy(actual)
                forged["attempts"][0]["solver_result"][field]=value
                self.assertFalse(runner.static_plane_attempt_acceptance(forged["attempts"][0])["passed"])
                with self.assertRaises(ValueError): audit_static_case_material(forged)
            forged=copy.deepcopy(actual)
            forged["attempts"][0]["design_target_vertices"] = np.asarray(
                forged["attempts"][0]["design_target_vertices"]).tolist()
            point=list(forged["attempts"][0]["design_target_vertices"][0]); point[2] += .125
            forged["attempts"][0]["design_target_vertices"][0]=point
            with self.assertRaisesRegex(ValueError,"design"):
                audit_static_case_material(forged)

    def test_feasible_full_state_without_optimality_is_not_promoted(self):
        runner = load_runner()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for objective, weight in (("LOWER_CENTER_OF_MASS",0.), ("PARETO_TRADEOFF",.001)):
                request = runner.build_request(method_case(root,objective,weight), asset_root=root)
                result = solve_spin_it_faster_static_plane_family(request,
                    initial_coefficients=(2.,.1,0.,-.5))
                self.assertNotEqual(result["status"], "CONVERGED")
                self.assertIsNone(result.get("boundary_optimality"))
            case = method_case(root,"PARETO_TRADEOFF",.5)
            case["target_origin"] = [.25,0.,0.]
            request = runner.build_request(case,asset_root=root)
            result = solve_spin_it_faster_static_plane_family(request,
                initial_coefficients=(2.,.1,0.,-.5))
            self.assertNotEqual(result["status"], "CONVERGED")
            self.assertIsNone(result.get("boundary_optimality"))

    def test_initial_full_state_witness_is_not_an_accepted_optimizer_update(self):
        runner = load_runner()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            request = runner.build_request(method_case(root,"PARETO_TRADEOFF",.5),asset_root=root)
            result = solve_spin_it_faster_static_plane_family(request,
                initial_coefficients=(2.,.1,0.,-.5))
            self.assertEqual(result["status"], "CONVERGED")
            self.assertEqual(result["accepted_steps"], 0)
            self.assertEqual(len(result["trace"]), 1)
            self.assertFalse(result["trace"][0]["accepted"])


if __name__ == "__main__": unittest.main()
