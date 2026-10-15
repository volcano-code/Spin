from __future__ import annotations

import copy
import unittest

from make_it_stand.verification import faster_theorem_exact


class FasterExactTheoremInstanceTests(unittest.TestCase):
    def test_exact_reference_instances_check_both_theorem_two_branches(self):
        certificate = faster_theorem_exact.reference_certificate()
        result = faster_theorem_exact.verify_certificate(certificate)

        self.assertEqual(result["status"], "EXACT_INSTANCE_CONDITIONS_VERIFIED")
        self.assertEqual(result["general_theorem_proof"], "NOT_CLAIMED")
        cases = {case["case_id"]: case for case in result["cases"]}
        self.assertEqual(set(cases), {
            "affine_cut_with_equality_constraint",
            "flat_objective_p_zero_branch",
        })
        cut = cases["affine_cut_with_equality_constraint"]
        self.assertEqual(cut["raw_moments"], [
            "3/4", "3/8", "3/8", "9/32", "1/4", "1/4", "9/64", "9/64", "9/64", "3/16",
        ])
        self.assertEqual(cut["theorem_1"], "EXACT_HYPOTHESES_AND_STRICT_SIGN_CONCLUSION")
        self.assertEqual(cut["theorem_2_branch"], "P_NONZERO_STRICT_LOCAL_OPTIMUM")
        self.assertEqual(cut["zero_set_measure"], "ZERO")
        self.assertEqual(cases["flat_objective_p_zero_branch"]["theorem_2_branch"],
                         "P_ZERO_UNBOUNDED_FIRST_ORDER_STATIONARITY")

    def test_exact_checker_rejects_a_wrong_material_cut(self):
        certificate = copy.deepcopy(faster_theorem_exact.reference_certificate())
        certificate["instances"][0]["material_box"]["maximum"][2] = "1/2"
        with self.assertRaisesRegex(ValueError, "material box does not equal the exact beta-positive region"):
            faster_theorem_exact.verify_certificate(certificate)

    def test_exact_checker_rejects_a_false_eq7_witness(self):
        certificate = copy.deepcopy(faster_theorem_exact.reference_certificate())
        certificate["instances"][0]["objective_gradient"][0] = "0"
        with self.assertRaisesRegex(ValueError, "Eq. 7 stationarity is not exact"):
            faster_theorem_exact.verify_certificate(certificate)

    def test_exact_checker_rejects_a_moment_vector_not_integrated_from_geometry(self):
        certificate = copy.deepcopy(faster_theorem_exact.reference_certificate())
        certificate["instances"][0]["raw_moments"][0] = "1/2"
        with self.assertRaisesRegex(ValueError, "raw moments do not match exact box integration"):
            faster_theorem_exact.verify_certificate(certificate)

    def test_exact_checker_rejects_a_non_unit_quadric_coefficient_vector(self):
        certificate = copy.deepcopy(faster_theorem_exact.reference_certificate())
        certificate["instances"][0]["beta_coefficients"][0] = "1"
        with self.assertRaisesRegex(ValueError, "not exactly unit norm"):
            faster_theorem_exact.verify_certificate(certificate)

    def test_exact_checker_rejects_wrong_p_zero_moments(self):
        certificate = copy.deepcopy(faster_theorem_exact.reference_certificate())
        certificate["instances"][1]["raw_moments"][0] = "1/3"
        with self.assertRaisesRegex(ValueError, "moments do not match exact constant characteristic"):
            faster_theorem_exact.verify_certificate(certificate)

    def test_exact_checker_rejects_missing_or_duplicate_theorem_branches(self):
        certificate = copy.deepcopy(faster_theorem_exact.reference_certificate())
        certificate["instances"].pop()
        with self.assertRaisesRegex(ValueError, "must contain both declared theorem instances"):
            faster_theorem_exact.verify_certificate(certificate)

        certificate = copy.deepcopy(faster_theorem_exact.reference_certificate())
        certificate["instances"][1]["case_id"] = certificate["instances"][0]["case_id"]
        with self.assertRaisesRegex(ValueError, "IDs are incomplete or duplicated"):
            faster_theorem_exact.verify_certificate(certificate)


if __name__ == "__main__":
    unittest.main()
