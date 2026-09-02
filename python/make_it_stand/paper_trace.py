"""Unified, profile-checked trace projection for paper simulation evidence."""

from __future__ import annotations

from enum import Enum
import hashlib
import json
import math
from typing import Any, Mapping


class TraceProfile(str, Enum):
    PAPER_PROFILE = "PAPER_PROFILE"
    RELEASED_BETA_PROFILE = "RELEASED_BETA_PROFILE"


def _hash(value: Any) -> str:
    payload = json.dumps(
        value,
        allow_nan=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _require_mapping(value: Any, name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be an object")
    return value


def _require_list(value: Any, name: str) -> list[Any]:
    if not isinstance(value, list):
        raise ValueError(f"{name} must be an array")
    return value


def _standing_state(event: Mapping[str, Any]) -> str:
    if bool(event.get("all_targets_reached")):
        return "SAFE"
    if bool(event.get("all_standing_feasible")):
        return "FEASIBLE"
    return "OUTSIDE"


def _validate_profile(
    diagnostics: Mapping[str, Any],
    profile: TraceProfile,
) -> Mapping[str, Any]:
    config = _require_mapping(diagnostics.get("configuration"), "configuration")
    line_search = _require_mapping(config.get("line_search"), "line_search")
    beta = profile is TraceProfile.RELEASED_BETA_PROFILE
    required = 0.03 if beta else 0.01
    actual = float(line_search.get("required_relative_decrease", math.nan))
    checks = {
        "required_relative_decrease": math.isclose(
            actual,
            required,
            rel_tol=0.0,
            abs_tol=1.0e-15,
        ),
        "freeze_targets_after_initial_refresh": (
            config.get("freeze_targets_after_initial_refresh") is beta
        ),
        "released_beta_outer_loop": config.get("released_beta_outer_loop") is beta,
        "released_beta_semantics": (
            line_search.get("released_beta_semantics") is beta
        ),
    }
    failed = [name for name, accepted in checks.items() if not accepted]
    if failed:
        raise ValueError(
            f"{profile.value} configuration mismatch: {', '.join(failed)}"
        )
    return config


def build_iteration_trace_v2(
    diagnostics: Mapping[str, Any],
    profile: TraceProfile | str,
) -> dict[str, Any]:
    """Project native diagnostics without conflating its three counters."""

    selected = TraceProfile(profile)
    config = _validate_profile(diagnostics, selected)
    problem = _require_mapping(diagnostics.get("problem"), "problem")
    standing = _require_mapping(
        problem.get("standing_constraint"),
        "problem.standing_constraint",
    )
    support_polygon = _require_list(
        standing.get("support_polygon"),
        "problem.standing_constraint.support_polygon",
    )
    safe_region_contract = {
        "support_polygon": support_polygon,
        "gravity": problem.get("gravity"),
        "gamma_radians": standing.get("gamma_radians"),
    }

    events = []
    event_lookup: dict[tuple[str, int, int], Mapping[str, Any]] = {}
    for raw in _require_list(diagnostics.get("r10_events"), "r10_events"):
        event = _require_mapping(raw, "r10 event")
        kind = str(event.get("kind"))
        carving_event = kind in {"carving", "carving_iteration"}
        deformation_event = kind in {"deformation", "accepted_deformation"}
        block = int(event.get("block", -1))
        iteration = int(event.get("iteration", -1))
        event_lookup[(kind, block, iteration)] = event
        if carving_event:
            event_lookup[("carving", block, iteration)] = event
        if deformation_event:
            event_lookup[("deformation", block, iteration)] = event
        events.append({
            "committed_event": int(event.get("ordinal", len(events))),
            "kind": kind,
            "paper_outer_iteration": block if block >= 0 else None,
            "beta_optimize_call": (
                iteration
                if selected is TraceProfile.RELEASED_BETA_PROFILE
                and deformation_event
                and iteration >= 0
                else None
            ),
            "carving_call": (
                iteration if carving_event and iteration >= 0 else None
            ),
            "deformation_trial": (
                iteration if deformation_event and iteration >= 0 else None
            ),
            "standing_state": _standing_state(event),
            "target_generation": event.get("target_generation"),
            "target_snapshot_hash": event.get("target_snapshot_hash"),
            "projected_target_energy": event.get("projected_target_energy"),
            "occupancy_hash": event.get("alpha_hash"),
            "geometry_hash": event.get("geometry_hash"),
        })

    carving_trials = []
    for index, raw in enumerate(
        _require_list(diagnostics.get("carving_iterations"), "carving_iterations")
    ):
        trial = _require_mapping(raw, "carving iteration")
        block = int(trial.get("block", -1))
        iteration = int(trial.get("iteration", index))
        carving_trials.append({
            "paper_outer_iteration": block if block >= 0 else None,
            "carving_call": iteration,
            "energy_before_same_context": trial.get("frozen_energy_before"),
            "energy_after_same_context": trial.get("frozen_energy_after"),
            "energy_after_refreshed_context": trial.get("refreshed_energy_after"),
            "target_generation_before": trial.get("target_generation_before"),
            "target_generation_after": trial.get("target_generation_after"),
            "occupancy_hash_before": trial.get("alpha_hash_before"),
            "occupancy_hash_after": trial.get("alpha_hash_after"),
            "stop_reason": trial.get("stop_reason"),
        })

    deformation_trials = []
    for index, raw in enumerate(
        _require_list(
            diagnostics.get("deformation_iterations"),
            "deformation_iterations",
        )
    ):
        trial = _require_mapping(raw, "deformation iteration")
        block = int(trial.get("block", -1))
        iteration = int(trial.get("iteration", index))
        trial_carving = int(trial.get("trial_carving_iterations", 0))
        rollback_carving = int(trial.get("rollback_carving_iterations", 0))
        if selected is TraceProfile.RELEASED_BETA_PROFILE:
            if trial_carving > 1 or rollback_carving > 1:
                raise ValueError(
                    "RELEASED_BETA_PROFILE balance call may perform at most one carving pass"
                )
            if trial.get("target_generation_before") != trial.get(
                "target_generation_after"
            ):
                raise ValueError("RELEASED_BETA_PROFILE target must remain frozen")
        elif trial_carving or rollback_carving:
            raise ValueError("PAPER_PROFILE fixed-alpha deformation cannot recarve")
        matched_event = event_lookup.get(("deformation", block, iteration), {})
        deformation_trials.append({
            "paper_outer_iteration": block if block >= 0 else None,
            "beta_optimize_call": (
                iteration
                if selected is TraceProfile.RELEASED_BETA_PROFILE
                else None
            ),
            "deformation_trial": iteration,
            "decision": trial.get("decision"),
            "accepted": trial.get("accepted"),
            "rolled_back": trial.get("rolled_back"),
            "standing_state_after": _standing_state(matched_event),
            "energy_before_same_context": trial.get("energy_before"),
            "trial_energy_same_context": trial.get("trial_energy"),
            "energy_after_same_context": trial.get("energy_after"),
            "energy_after_refreshed_context": trial.get("refreshed_energy_after"),
            "step_attempted": trial.get("attempted_step"),
            "step_after": trial.get("step_after"),
            "mu_before": trial.get("mu_before"),
            "mu_after": trial.get("mu_after"),
            "lambda": config.get("shape_scale_lambda"),
            "target_generation_before": trial.get("target_generation_before"),
            "target_generation_after": trial.get("target_generation_after"),
            "occupancy_hash_before": trial.get("alpha_hash_before"),
            "occupancy_hash_after": trial.get("alpha_hash_after"),
            "geometry_hash_before": trial.get("geometry_hash_before"),
            "geometry_hash_after": trial.get("geometry_hash_after"),
            "balance_calls": trial.get("trial_balance_calls", 0),
            "carving_calls": trial_carving,
            "rollback_balance_calls": trial.get("rollback_balance_calls", 0),
            "rollback_carving_calls": rollback_carving,
        })

    return {
        "schema": "make-it-stand-iteration-trace-v2",
        "schema_version": 2,
        "profile": selected.value,
        "profile_contract_sha256": _hash(config),
        "support_polygon_sha256": _hash(support_polygon),
        "safe_region_sha256": _hash(safe_region_contract),
        "events": events,
        "carving_trials": carving_trials,
        "deformation_trials": deformation_trials,
    }


__all__ = ["TraceProfile", "build_iteration_trace_v2"]
