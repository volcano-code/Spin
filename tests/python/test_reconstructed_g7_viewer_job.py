import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiments/teapot_spin2014"))

from reconstructed_g7_viewer_admission import validate_viewer_admission_result  # noqa: E402


class ReconstructedG7ViewerAdmissionResultTests(unittest.TestCase):
    def _summary(self):
        return {
            "artifact_digest": "a" * 64,
            "problem_digest": "p" * 64,
            "geometry_status": "MESH_MATERIAL_VALIDATED",
            "support_display": "TARGET_POINTS_ONLY_NO_CERTIFIED_HULL",
            "disclosures": {
                "project_status": "PARTIAL_BLOCKED",
                "physical_validation": "NOT_IN_SCOPE",
                "solver_target_feasibility": "UNKNOWN",
            },
            "loaded_paths": {
                "native": "/build/recorded/_make_it_stand.cpython-313-x86_64-linux-gnu.so",
                "package": "/project/python/make_it_stand/__init__.py",
                "adapter": "/project/python/make_it_stand/verification/spin_view.py",
            },
            "layers": {"original": {"faces": 100}, "outer": {"faces": 90},
                       "material": {"faces": 80}, "cavity": {"faces": 10}},
        }

    def test_accepts_only_expected_bundle_summary_and_exact_loaded_extension(self):
        summary = self._summary()
        text = json.dumps(summary)
        accepted = validate_viewer_admission_result(
            0, text, "", expected_artifact="a" * 64, expected_problem="p" * 64,
            expected_extension="/build/recorded/_make_it_stand.cpython-313-x86_64-linux-gnu.so")
        self.assertEqual(accepted["artifact_digest"], "a" * 64)
        for returncode, mutated, expected_error in (
            (1, summary, "NATIVE_VIEWER_CHECK_FAILED"),
            (0, {**summary, "artifact_digest": "b" * 64}, "NATIVE_VIEWER_ARTIFACT_DIGEST_MISMATCH"),
            (0, {**summary, "problem_digest": "q" * 64}, "NATIVE_VIEWER_PROBLEM_DIGEST_MISMATCH"),
            (0, {**summary, "loaded_paths": {**summary["loaded_paths"], "native": "/other/module.so"}},
             "NATIVE_VIEWER_EXTENSION_PATH_MISMATCH"),
        ):
            with self.subTest(expected_error=expected_error):
                with self.assertRaisesRegex(ValueError, expected_error):
                    validate_viewer_admission_result(returncode, json.dumps(mutated), "",
                        expected_artifact="a" * 64, expected_problem="p" * 64,
                        expected_extension="/build/recorded/_make_it_stand.cpython-313-x86_64-linux-gnu.so")


if __name__ == "__main__":
    unittest.main()
