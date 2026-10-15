import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
STRICT_RUNNERS = ROOT / "experiments" / "spin_it_reproduction" / "strict"
sys.path.insert(0, str(STRICT_RUNNERS))

from run_faster_fullrank_end_to_end import (  # noqa: E402
    parse_unknown_seed_diagnostics,
    parse_unknown_seed_records,
    require_unknown_seed_replay,
)


def record(seed: int, status: str, solver_diagnostic: str = "not-converged",
           restarts: int = 1, evaluations: int = 60,
           restart_evaluations: str = "60") -> str:
    digest = f"{seed:064x}"
    return (
        f"seed={seed} generated_start={digest} status={status} "
        f"evaluations={evaluations} trace=3 restarts={restarts} "
        f"restart_evaluations={restart_evaluations} residual=0.25 "
        f"final_verified=0 input_hash={digest} event_hash={digest} "
        f"trace_hash={digest} output_hash={digest} "
        f"initial_subtype=ELLIPSOID initial_source_motion_defined=1 "
        f"initial_source_diagnostic_hex={b'regular'.hex()} "
        f"first_derivative_failure_evaluation=4 first_derivative_failure_iteration=2 "
        f"first_derivative_failure_subtype=ONE_SHEET_HYPERBOLOID "
        f"first_derivative_failure_step=0.125 "
        f"first_derivative_failure_hex={b'unmatched root'.hex()} "
        f"first_derivative_failure_q=0,1,2,3,4,5,6,7,8,9 "
        f"source_witness_vertex=7 source_witness_nearest_normal=3 "
        f"source_witness_nearest_distance=0.25 source_witness_radius_delta=0.5 "
        f"source_witness_q_residual=0.75 source_witness_point=1,2,3 "
        f"solver_diagnostic_hex={solver_diagnostic.encode().hex()} "
        f"verification_diagnostic_hex={b'not-audited'.hex()}"
    )


def transcript(statuses=None, *, changed_event_seed=None, seeds=range(5)):
    statuses = statuses or {seed: "RESOURCE_LIMIT_REACHED" for seed in range(5)}
    lines = []
    for seed in seeds:
        digest = f"{seed:064x}"
        lines.extend([
            f"initial-start seed={seed} restart=0 hash={digest} coefficients=1,0,0,0,0,0,0,0,0,0",
            f"initial-start seed={seed} restart=1 hash={seed + 10:064x} coefficients=2,0,0,0,0,0,0,0,0,0",
            f"solver-event seed={seed} restart=0 iteration=0 evaluation=0 type=attempt_start outcome=started reason=none detail_hex= coefficients=1,0,0,0,0,0,0,0,0,0",
        ])
        if changed_event_seed == seed:
            lines.append(
                f"solver-event seed={seed} restart=0 iteration=1 evaluation=1 type=evaluation outcome=failed reason=changed detail_hex=626164 coefficients=1,0,0,0,0,0,0,0,0,0"
            )
        lines.append(
            f"solver-event seed={seed} restart=0 iteration=3 evaluation=60 type=attempt_terminal outcome=terminated reason=budget detail_hex= coefficients=none"
        )
        lines.append(record(seed, statuses[seed]))
    if tuple(seeds) == tuple(range(5)):
        lines.append("[PASS] unknown-seed-sweep")
    return "\n".join(lines)


class FasterUnknownSeedEvidenceTests(unittest.TestCase):
    def test_preserves_resource_limited_and_geometry_failed_seeds(self):
        output = "\n".join([
            record(0, "RESOURCE_LIMIT_REACHED"),
            record(1, "GEOMETRY_EVALUATION_FAILED", "Boolean audit failed"),
            record(2, "LINE_SEARCH_EXHAUSTED"),
            record(3, "OPTIMIZATION_DID_NOT_CONVERGE"),
            record(4, "RESOURCE_LIMIT_REACHED"),
            "[PASS] unknown-seed-sweep",
        ])
        rows = parse_unknown_seed_records(output)
        self.assertEqual([0, 1, 2, 3, 4], [row["seed"] for row in rows])
        self.assertEqual("RESOURCE_LIMIT_REACHED", rows[0]["status"])
        self.assertEqual("GEOMETRY_EVALUATION_FAILED", rows[1]["status"])
        self.assertEqual("Boolean audit failed", rows[1]["solver_diagnostic"])
        self.assertEqual("not-audited", rows[1]["verification_diagnostic"])
        self.assertEqual("ELLIPSOID", rows[0]["initial_subtype"])
        self.assertTrue(rows[0]["initial_source_motion_defined"])
        self.assertEqual("regular", rows[0]["initial_source_diagnostic"])
        self.assertEqual(4, rows[0]["first_derivative_failure_evaluation"])
        self.assertEqual("ONE_SHEET_HYPERBOLOID",
                         rows[0]["first_derivative_failure_subtype"])
        self.assertEqual("unmatched root", rows[0]["first_derivative_failure"])
        self.assertEqual("0,1,2,3,4,5,6,7,8,9",
                         rows[0]["first_derivative_failure_q"])
        self.assertEqual(7, rows[0]["source_witness_vertex"])
        self.assertEqual(3, rows[0]["source_witness_nearest_normal"])
        self.assertEqual("1,2,3", rows[0]["source_witness_point"])

    def test_rejects_missing_or_duplicate_seed_records(self):
        incomplete = "\n".join(record(seed, "RESOURCE_LIMIT_REACHED") for seed in range(4))
        duplicated = "\n".join(
            [record(seed, "RESOURCE_LIMIT_REACHED") for seed in range(4)]
            + [record(3, "RESOURCE_LIMIT_REACHED")]
        )
        with self.assertRaises(ValueError):
            parse_unknown_seed_records(incomplete)
        with self.assertRaises(ValueError):
            parse_unknown_seed_records(duplicated)

    def test_rejects_restart_count_above_preregistered_per_seed_limit(self):
        output = "\n".join([
            record(0, "RESOURCE_LIMIT_REACHED", restarts=3),
            *[record(seed, "RESOURCE_LIMIT_REACHED") for seed in range(1, 5)],
            "[PASS] unknown-seed-sweep",
        ])
        with self.assertRaisesRegex(ValueError, "restart count exceeds frozen limit"):
            parse_unknown_seed_records(output)

    def test_rejects_per_restart_evaluation_budget_overrun(self):
        output = "\n".join([
            record(0, "RESOURCE_LIMIT_REACHED", restarts=2,
                   evaluations=121, restart_evaluations="60,61"),
            *[record(seed, "RESOURCE_LIMIT_REACHED") for seed in range(1, 5)],
            "[PASS] unknown-seed-sweep",
        ])
        with self.assertRaisesRegex(ValueError, "per-restart budget"):
            parse_unknown_seed_records(output)

    def test_requires_independent_process_scientific_replay_to_match(self):
        first = transcript()
        require_unknown_seed_replay(first, first)
        changed = transcript(changed_event_seed=3)
        with self.assertRaises(ValueError):
            require_unknown_seed_replay(first, changed)

    def test_rejects_missing_or_duplicate_restart_start_evidence(self):
        first = transcript()
        missing = first.replace(
            f"initial-start seed=2 restart=1 hash={12:064x} coefficients=2,0,0,0,0,0,0,0,0,0\n", ""
        )
        with self.assertRaisesRegex(ValueError, "exactly two starts"):
            require_unknown_seed_replay(first, missing)

    def test_rejects_duplicate_attempt_start_event(self):
        first = transcript()
        duplicate = first.replace(
            f"solver-event seed=1 restart=0 iteration=0 evaluation=0 type=attempt_start outcome=started reason=none detail_hex= coefficients=1,0,0,0,0,0,0,0,0,0\n",
            f"solver-event seed=1 restart=0 iteration=0 evaluation=0 type=attempt_start outcome=started reason=none detail_hex= coefficients=1,0,0,0,0,0,0,0,0,0\n"
            f"solver-event seed=1 restart=0 iteration=0 evaluation=0 type=attempt_start outcome=started reason=none detail_hex= coefficients=1,0,0,0,0,0,0,0,0,0\n",
        )
        with self.assertRaisesRegex(ValueError, "exactly one start"):
            parse_unknown_seed_diagnostics(duplicate)

    def test_requires_one_terminal_event_for_each_complete_attempt(self):
        first = transcript()
        missing = first.replace(
            "solver-event seed=1 restart=0 iteration=3 evaluation=60 type=attempt_terminal outcome=terminated reason=budget detail_hex= coefficients=none\n",
            "",
        )
        with self.assertRaisesRegex(ValueError, "exactly one terminal"):
            parse_unknown_seed_diagnostics(missing)

    def test_rejects_terminal_evaluation_count_inconsistent_with_budget_summary(self):
        first = transcript()
        inconsistent = first.replace(
            "solver-event seed=1 restart=0 iteration=3 evaluation=60 type=attempt_terminal outcome=terminated reason=budget detail_hex= coefficients=none",
            "solver-event seed=1 restart=0 iteration=3 evaluation=59 type=attempt_terminal outcome=terminated reason=budget detail_hex= coefficients=none",
        )
        with self.assertRaisesRegex(ValueError, "terminal evaluation index"):
            parse_unknown_seed_diagnostics(inconsistent)

    def test_attempt_start_coefficients_match_frozen_restart(self):
        first = transcript()
        inconsistent = first.replace(
            "solver-event seed=1 restart=0 iteration=0 evaluation=0 type=attempt_start outcome=started reason=none detail_hex= coefficients=1,0,0,0,0,0,0,0,0,0",
            "solver-event seed=1 restart=0 iteration=0 evaluation=0 type=attempt_start outcome=started reason=none detail_hex= coefficients=2,0,0,0,0,0,0,0,0,0",
        )
        with self.assertRaisesRegex(ValueError, "does not match frozen start"):
            parse_unknown_seed_diagnostics(inconsistent)

    def test_accepts_one_seed_diagnostic_without_marking_a_full_sweep(self):
        from run_faster_fullrank_end_to_end import parse_unknown_seed_diagnostics

        parsed = parse_unknown_seed_diagnostics(
            transcript(seeds=(2,)), expected_seeds=(2,), require_completion=False
        )
        self.assertEqual(2, len(parsed["initial_starts"]))
        self.assertEqual(1, sum(
            event["type"] == "attempt_start" for event in parsed["solver_events"]
        ))


if __name__ == "__main__":
    unittest.main()
