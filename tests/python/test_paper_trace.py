"""Unified paper/released-beta iteration trace contracts."""

from __future__ import annotations

import unittest

from make_it_stand.paper_trace import TraceProfile, build_iteration_trace_v2


def _fixture(profile: TraceProfile) -> dict:
    beta = profile is TraceProfile.RELEASED_BETA_PROFILE
    return {
        "configuration": {
            "freeze_targets_after_initial_refresh": beta,
            "released_beta_outer_loop": beta,
            "line_search": {
                "required_relative_decrease": 0.03 if beta else 0.01,
                "released_beta_semantics": beta,
            },
            "initial_mu": 0.75,
            "shape_scale_lambda": 20.0,
        },
        "problem": {
            "gravity": [0.0, -1.0, 0.0],
            "standing_constraint": {
                "gamma_radians": 0.08726646259971647,
                "support_polygon": [
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [1.0, 0.0, 1.0],
                    [0.0, 0.0, 1.0],
                ],
            },
        },
        "r10_events": [
            {
                "ordinal": 0,
                "kind": "initial",
                "block": -1,
                "iteration": -1,
                "all_standing_feasible": False,
                "all_targets_reached": False,
                "target_generation": 1,
                "target_snapshot_hash": "target-1",
                "alpha_hash": "alpha-0",
                "geometry_hash": "geometry-0",
                "projected_target_energy": 0.4,
            },
            {
                "ordinal": 1,
                "kind": "carving_iteration",
                "block": 0,
                "iteration": 0,
                "all_standing_feasible": True,
                "all_targets_reached": False,
                "target_generation": 1,
                "target_snapshot_hash": "target-1",
                "alpha_hash": "alpha-1",
                "geometry_hash": "geometry-1",
                "projected_target_energy": 0.1,
            },
            {
                "ordinal": 2,
                "kind": "accepted_deformation",
                "block": 0,
                "iteration": 0,
                "all_standing_feasible": True,
                "all_targets_reached": True,
                "target_generation": 1 if beta else 2,
                "target_snapshot_hash": "target-1" if beta else "target-2",
                "alpha_hash": "alpha-1",
                "geometry_hash": "geometry-2",
                "projected_target_energy": 0.0,
            },
        ],
        "carving_iterations": [{
            "block": 0,
            "iteration": 0,
            "frozen_energy_before": 0.4,
            "frozen_energy_after": 0.1,
            "refreshed_energy_after": 0.1,
            "target_generation_before": 1,
            "target_generation_after": 1,
            "alpha_hash_before": "alpha-0",
            "alpha_hash_after": "alpha-1",
            "stop_reason": "none",
        }],
        "deformation_iterations": [{
            "block": 0,
            "iteration": 0,
            "decision": "accepted",
            "accepted": True,
            "rolled_back": False,
            "energy_before": 0.1,
            "trial_energy": 0.05,
            "energy_after": 0.05,
            "refreshed_energy_after": 0.0,
            "attempted_step": 1.0,
            "step_after": 1.0,
            "mu_before": 0.75,
            "mu_after": 0.75,
            "target_generation_before": 1,
            "target_generation_after": 1 if beta else 2,
            "alpha_hash_before": "alpha-1",
            "alpha_hash_after": "alpha-1",
            "geometry_hash_before": "geometry-1",
            "geometry_hash_after": "geometry-2",
            "trial_balance_calls": 1 if beta else 0,
            "rollback_balance_calls": 0,
            "trial_carving_iterations": 1 if beta else 0,
            "rollback_carving_iterations": 0,
        }],
    }


class PaperTraceTests(unittest.TestCase):
    def test_trace_separates_event_and_algorithm_counters(self) -> None:
        trace = build_iteration_trace_v2(
            _fixture(TraceProfile.RELEASED_BETA_PROFILE),
            TraceProfile.RELEASED_BETA_PROFILE,
        )
        self.assertEqual(trace["schema_version"], 2)
        self.assertEqual(trace["profile"], "RELEASED_BETA_PROFILE")
        self.assertEqual(trace["events"][0]["standing_state"], "OUTSIDE")
        self.assertEqual(trace["events"][1]["standing_state"], "FEASIBLE")
        self.assertEqual(trace["events"][2]["standing_state"], "SAFE")
        self.assertEqual(trace["events"][2]["beta_optimize_call"], 0)
        self.assertEqual(trace["carving_trials"][0]["carving_call"], 0)
        self.assertEqual(trace["deformation_trials"][0]["deformation_trial"], 0)
        self.assertEqual(trace["deformation_trials"][0]["carving_calls"], 1)
        self.assertRegex(trace["support_polygon_sha256"], r"^[0-9a-f]{64}$")
        self.assertRegex(trace["safe_region_sha256"], r"^[0-9a-f]{64}$")

    def test_profile_contracts_reject_mixed_semantics(self) -> None:
        mixed = _fixture(TraceProfile.PAPER_PROFILE)
        mixed["configuration"]["line_search"]["required_relative_decrease"] = 0.03
        with self.assertRaisesRegex(ValueError, "PAPER_PROFILE"):
            build_iteration_trace_v2(mixed, TraceProfile.PAPER_PROFILE)

        invalid_beta = _fixture(TraceProfile.RELEASED_BETA_PROFILE)
        invalid_beta["deformation_iterations"][0]["trial_carving_iterations"] = 2
        with self.assertRaisesRegex(ValueError, "at most one carving"):
            build_iteration_trace_v2(
                invalid_beta,
                TraceProfile.RELEASED_BETA_PROFILE,
            )

    def test_paper_dynamic_target_generation_is_preserved(self) -> None:
        trace = build_iteration_trace_v2(
            _fixture(TraceProfile.PAPER_PROFILE),
            TraceProfile.PAPER_PROFILE,
        )
        self.assertEqual(trace["events"][-1]["target_generation"], 2)
        self.assertEqual(trace["deformation_trials"][0]["target_generation_after"], 2)
        self.assertEqual(trace["deformation_trials"][0]["carving_calls"], 0)


if __name__ == "__main__":
    unittest.main()
