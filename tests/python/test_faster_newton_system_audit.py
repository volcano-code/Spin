import unittest
import json
from pathlib import Path
import tempfile

from make_it_stand.verification.newton_system_audit import (
    audit_newton_event_journal,
    audit_newton_system_payload,
)


class NewtonSystemAuditTests(unittest.TestCase):
    def test_exactly_full_rank_badly_scaled_system_is_solved_and_checked_in_original_units(self):
        payload = {
            "schema": "spin-newton-system-payload-v2",
            "matrix": [[1.0e12, 0.0], [0.0, 1.0e-6]],
            "right_hand_side": [2.0e12, -3.0e-6],
            "multiplier_mu": 1.0,
        }

        report = audit_newton_system_payload(payload)

        self.assertEqual(report["exact_binary64_rank"], 2)
        self.assertEqual(report["scaled_numerical_rank"], 2)
        self.assertEqual(report["interpretation"], "SCALE_SENSITIVE_FULL_RANK")
        self.assertEqual(report["solution"], [2.0, -3.0])
        self.assertLessEqual(report["original_system_backward_error"], 1e-15)

    def test_exactly_singular_system_is_not_repaired_by_scaling(self):
        payload = {
            "schema": "spin-newton-system-payload-v2",
            "matrix": [[1.0, 1.0], [2.0, 2.0]],
            "right_hand_side": [1.0, 2.0],
            "multiplier_mu": 0.5,
        }

        report = audit_newton_system_payload(payload)

        self.assertEqual(report["exact_binary64_rank"], 1)
        self.assertEqual(report["interpretation"], "EXACTLY_RANK_DEFICIENT")
        self.assertIsNone(report["solution"])

    def test_audit_rejects_nonfinite_or_non_square_newton_payloads(self):
        with self.assertRaisesRegex(ValueError, "square"):
            audit_newton_system_payload({
                "schema": "spin-newton-system-payload-v2",
                "matrix": [[1.0, 2.0]],
                "right_hand_side": [1.0],
            })
        with self.assertRaisesRegex(ValueError, "finite"):
            audit_newton_system_payload({
                "schema": "spin-newton-system-payload-v2",
                "matrix": [[float("nan")]],
                "right_hand_side": [1.0],
            })

    def test_event_journal_replay_keeps_event_identity_and_hash(self):
        payload = {
            "schema": "spin-newton-system-payload-v2",
            "matrix": [[1.0e12, 0.0], [0.0, 1.0e-6]],
            "right_hand_side": [2.0e12, -3.0e-6],
            "multiplier_mu": 1.0,
        }
        with tempfile.TemporaryDirectory() as temporary:
            journal = Path(temporary) / "events.jsonl"
            journal.write_text("\n".join((
                json.dumps({"schema": "solver-event-journal-v1"}),
                json.dumps({"record_type": "solver_event", "record": {
                    "restart": 1, "iteration": 4, "evaluation_index": 12,
                    "reason": "bordered_system_rank_deficient",
                    "newton_system_payload": payload,
                }}),
                json.dumps({"record_type": "geometry_evaluation_finished", "record": {}}),
            )) + "\n", encoding="utf-8")

            replay = audit_newton_event_journal(journal)

        self.assertEqual(replay["payload_count"], 1)
        self.assertEqual(len(replay["events"]), 1)
        self.assertEqual(replay["events"][0]["evaluation_index"], 12)
        self.assertEqual(replay["events"][0]["restart"], 1)
        self.assertEqual(replay["events"][0]["audit"]["exact_binary64_rank"], 2)
        self.assertRegex(replay["journal_sha256"], r"^[0-9a-f]{64}$")
        self.assertRegex(replay["auditor"]["source_sha256"], r"^[0-9a-f]{64}$")
        self.assertTrue(replay["auditor"]["python_version"])
        self.assertTrue(replay["auditor"]["numpy_version"])


if __name__ == "__main__":
    unittest.main()
