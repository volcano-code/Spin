"""Exact, pure-Python checks for the two published Faster Section 6.6 cases."""

from __future__ import annotations

import copy
import importlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
SPEC = ROOT / "experiments/spin_it_reproduction/sources/faster_section_6_6_spec.json"
PDF = ROOT / "build/strict-paper-sources/sif-final.pdf"
RUNNER = ROOT / "experiments/spin_it_reproduction/strict/run_faster_section_6_6.py"


class FasterSection66ExactTests(unittest.TestCase):
    def setUp(self):
        self.assertIsNotNone(
            importlib.util.find_spec("make_it_stand.verification.faster_section_6_6"),
            "exact Faster Section 6.6 checker is missing",
        )
        self.api = importlib.import_module("make_it_stand.verification.faster_section_6_6")
        self.spec = self.api.load_spec(SPEC)
        self.witnesses = self.api.published_witness_definitions()

    def verify(self, spec=None, witnesses=None):
        return self.api.verify_published_section_6_6(
            self.spec if spec is None else spec,
            self.witnesses if witnesses is None else witnesses,
        )

    def test_exact_computations_prove_both_published_outcomes(self):
        result = self.verify()
        self.assertEqual(result["project_status"], "PARTIAL_BLOCKED")
        self.assertEqual(result["author_bitwise_exact"], "NOT_VERIFIED")
        self.assertEqual(result["physical_validation"], "NOT_IN_SCOPE")
        self.assertEqual(result["case_count"], 2)
        by_id = {case["source_case_id"]: case for case in result["cases"]}

        nonunique = by_id["faster.section6_6.nonunique"]
        self.assertEqual(nonunique["catalogue_id"], "faster.section6_6.example.1")
        self.assertEqual(nonunique["scientific_outcome"], "VERIFIED_ANALYTIC_PROPERTIES")
        self.assertEqual(nonunique["volume_polynomial_in_tau"], ["1/2", "0"])
        self.assertEqual(nonunique["rectangle_corner_values"], ["1", "0", "0", "1"])
        self.assertEqual(nonunique["parameter_derivative_l2_squared"], "1/3")
        self.assertEqual(nonunique["binary_half_interval_integrals"], ["1/2", "1/2"])
        self.assertEqual(nonunique["constant_nonbinary_integral"], "1/2")
        self.assertEqual(nonunique["shifted_binary_overlap_polynomial"], ["1/2", "-1"])
        self.assertEqual(nonunique["shifted_binary_symmetric_difference_polynomial"], ["0", "2"])
        self.assertTrue(nonunique["uncountable_affine_family_modulo_ae"])
        self.assertTrue(nonunique["uncountable_binary_family_modulo_ae"])

        full = by_id["faster.section6_6.full_boundary"]
        self.assertEqual(full["catalogue_id"], "faster.section6_6.example.2")
        self.assertEqual(full["scientific_outcome"], "VERIFIED_ANALYTIC_PROPERTIES")
        self.assertEqual(full["domain_measure"], "1/2")
        self.assertEqual(full["witness_integral"], "1/2")
        self.assertEqual(full["objective"], "0")
        self.assertEqual(full["objective_derivative"], "0")
        self.assertEqual(full["p"], "0")
        self.assertEqual(full["uniqueness"], "UNIQUE_ALMOST_EVERYWHERE")
        self.assertTrue(full["nonnegative_zero_integral_lemma_applies"])
        self.assertEqual(full["deficit_range_from_box"], ["0", "1"])
        self.assertFalse(result["claims_all_p_zero_problems_classified"])
        self.assertEqual(result["native_extensions_required"], [])

    def test_corrupted_witness_coefficients_and_parameter_ranges_fail(self):
        mutations = []
        changed = copy.deepcopy(self.witnesses)
        changed["affine_family"]["tau_coefficient"] = ["-1", "3"]
        mutations.append(changed)
        changed = copy.deepcopy(self.witnesses)
        changed["affine_family"]["tau_upper"] = "0"
        mutations.append(changed)
        changed = copy.deepcopy(self.witnesses)
        changed["binary_interval_family"]["parameter_upper"] = "0"
        mutations.append(changed)
        changed = copy.deepcopy(self.witnesses)
        changed["full_domain_witness"]["constant"] = "0"
        mutations.append(changed)
        for mutation in mutations:
            with self.subTest(mutation=mutation), self.assertRaises(ValueError):
                self.verify(witnesses=mutation)

    def test_affine_interval_and_zero_integral_conditions_are_derived(self):
        interval = self.api.derive_shifted_binary_interval_identity("2/5")
        self.assertEqual(interval["interval_measure_polynomial"], ["2/5", "0"])
        self.assertEqual(interval["overlap_polynomial"], ["2/5", "-1"])
        self.assertEqual(interval["symmetric_difference_polynomial"], ["0", "2"])

        invalid_box = self.api.derive_full_domain_uniqueness_conditions(
            "0", "2", "1/2", "1/2"
        )
        self.assertEqual(invalid_box["deficit_range_from_box"], ["-1", "1"])
        self.assertFalse(invalid_box["box_makes_deficit_nonnegative"])
        self.assertTrue(invalid_box["deficit_integral_at_any_optimum_is_zero"])
        self.assertFalse(invalid_box["nonnegative_zero_integral_lemma_applies"])

    def test_mutated_published_definitions_missing_and_duplicate_cases_fail(self):
        mutations = []
        changed = copy.deepcopy(self.spec)
        changed["objective"] = "(v(chi) - 1/3)^2"
        mutations.append(changed)
        changed = copy.deepcopy(self.spec)
        changed["constraints"] = ["chi is binary"]
        mutations.append(changed)
        changed = copy.deepcopy(self.spec)
        changed["cases"][0]["domain"]["upper"] = 2
        mutations.append(changed)
        changed = copy.deepcopy(self.spec)
        changed["cases"][0]["domain"]["upper"] = 1.0
        mutations.append(changed)
        changed = copy.deepcopy(self.spec)
        changed["cases"][0]["domain"]["lower"] = 0.0
        mutations.append(changed)
        changed = copy.deepcopy(self.spec)
        changed["cases"][0]["domain"]["dimension"] = True
        mutations.append(changed)
        changed = copy.deepcopy(self.spec)
        changed["cases"][0]["case_id"] = "renamed"
        mutations.append(changed)
        changed = copy.deepcopy(self.spec)
        changed["cases"] = changed["cases"][:1]
        mutations.append(changed)
        changed = copy.deepcopy(self.spec)
        changed["cases"] = [changed["cases"][0], copy.deepcopy(changed["cases"][0])]
        mutations.append(changed)
        changed = copy.deepcopy(self.spec)
        changed["passed"] = True
        changed["objective"] = "(v(chi) - 1/3)^2"
        mutations.append(changed)
        for mutation in mutations:
            with self.subTest(mutation=mutation), self.assertRaises(ValueError):
                self.verify(spec=mutation)

    def test_nonfinite_and_duplicate_key_json_fail_closed(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "bad.json"
            path.write_text('{"schema":"x","value":NaN}', encoding="utf-8")
            with self.assertRaises(ValueError):
                self.api.load_spec(path)
            path.write_text('{"schema":"x","schema":"y"}', encoding="utf-8")
            with self.assertRaises(ValueError):
                self.api.load_spec(path)


class FasterSection66RunnerTests(unittest.TestCase):
    def test_consistently_rehashed_wrong_science_is_not_current(self):
        api = importlib.import_module("make_it_stand.verification.faster_section_6_6")
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "run"
            result = self.run_cli(SPEC, PDF, output)
            self.assertEqual(result.returncode, 0, result.stderr)
            result_path = output / "result.json"
            evidence_path = output / "evidence.json"
            scientific = json.loads(result_path.read_text())
            evidence = json.loads(evidence_path.read_text())
            scientific["scientific_data"]["cases"][0]["objective"] = "1"
            digest = api.scientific_sha256(scientific["scientific_data"])
            scientific["scientific_data_sha256"] = digest
            result_path.write_text(api.canonical_json(scientific))
            evidence["outputs"]["result"] = api.file_record(result_path)
            evidence["scientific_data_sha256"] = digest
            evidence_path.write_text(api.canonical_json(evidence))
            assessment = api.check_evidence_freshness(evidence_path)
            self.assertEqual(assessment["status"], "STALE_RESULT")
            self.assertTrue(any("recomputed" in x for x in assessment["mismatches"]))

    def run_cli(self, spec, paper, output):
        environment = dict(os.environ)
        environment["PYTHONPATH"] = str(ROOT / "python")
        return subprocess.run(
            [sys.executable, str(RUNNER), "--spec", str(spec), "--paper", str(paper),
             "--output-directory", str(output)],
            cwd=ROOT,
            env=environment,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

    def test_cli_twice_is_deterministic_exclusive_and_pure_python(self):
        with tempfile.TemporaryDirectory() as temporary:
            first = Path(temporary) / "first"
            second = Path(temporary) / "second"
            run1 = self.run_cli(SPEC, PDF, first)
            run2 = self.run_cli(SPEC, PDF, second)
            self.assertEqual(run1.returncode, 0, run1.stderr)
            self.assertEqual(run2.returncode, 0, run2.stderr)
            first_result = json.loads((first / "result.json").read_text(encoding="utf-8"))
            second_result = json.loads((second / "result.json").read_text(encoding="utf-8"))
            self.assertEqual(first_result["scientific_data_sha256"], second_result["scientific_data_sha256"])
            self.assertEqual(first_result["scientific_data"], second_result["scientific_data"])
            evidence = json.loads((first / "evidence.json").read_text(encoding="utf-8"))
            self.assertEqual(evidence["runtime"]["native_extensions_required"], [])
            self.assertEqual(evidence["build_source_attestation"], "NOT_APPLICABLE_PURE_PYTHON_FORMULA_PROOF")
            self.assertEqual(evidence["author_bitwise_exact"], "NOT_VERIFIED")
            self.assertEqual(evidence["physical_validation"], "NOT_IN_SCOPE")
            overwrite = self.run_cli(SPEC, PDF, first)
            self.assertNotEqual(overwrite.returncode, 0)
            self.assertIn("already exists", overwrite.stderr)

    def test_cli_rejects_bad_source_hash_and_freshness_rechecks_files(self):
        self.assertIsNotNone(
            importlib.util.find_spec("make_it_stand.verification.faster_section_6_6"),
            "exact Faster Section 6.6 checker is missing",
        )
        api = importlib.import_module("make_it_stand.verification.faster_section_6_6")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            bad_pdf = root / "sif-final.pdf"
            shutil.copyfile(PDF, bad_pdf)
            with bad_pdf.open("ab") as stream:
                stream.write(b"corruption")
            rejected = self.run_cli(SPEC, bad_pdf, root / "bad-output")
            self.assertNotEqual(rejected.returncode, 0)
            self.assertIn("PDF SHA256 mismatch", rejected.stderr)

            output = root / "valid-output"
            accepted = self.run_cli(SPEC, PDF, output)
            self.assertEqual(accepted.returncode, 0, accepted.stderr)
            evidence_path = output / "evidence.json"
            self.assertEqual(api.check_evidence_freshness(evidence_path)["status"], "CURRENT")
            original = json.loads(evidence_path.read_text(encoding="utf-8"))
            for label in ("checker", "spec", "paper_pdf"):
                changed = copy.deepcopy(original)
                changed["inputs"][label]["sha256"] = "0" * 64
                evidence_path.write_text(json.dumps(changed), encoding="utf-8")
                with self.subTest(label=label):
                    assessment = api.check_evidence_freshness(evidence_path)
                    self.assertEqual(assessment["status"], "STALE_RESULT")
                    self.assertTrue(any(label in item for item in assessment["mismatches"]))
            evidence_path.write_text(json.dumps(original), encoding="utf-8")

            for label, unrelated in (("checker", SPEC), ("runner", PDF)):
                changed = copy.deepcopy(original)
                changed["inputs"][label] = api.file_record(unrelated)
                changed["source_content_sha256"] = api.scientific_sha256(changed["inputs"])
                evidence_path.write_text(json.dumps(changed), encoding="utf-8")
                with self.subTest(substituted=label):
                    assessment = api.check_evidence_freshness(evidence_path)
                    self.assertEqual(assessment["status"], "STALE_RESULT")
                    self.assertTrue(any(label in item for item in assessment["mismatches"]))

            changed = copy.deepcopy(original)
            changed["source_content_sha256"] = "0" * 64
            evidence_path.write_text(json.dumps(changed), encoding="utf-8")
            assessment = api.check_evidence_freshness(evidence_path)
            self.assertEqual(assessment["status"], "STALE_RESULT")
            self.assertTrue(any("source_content_sha256" in item for item in assessment["mismatches"]))


if __name__ == "__main__":
    unittest.main()
