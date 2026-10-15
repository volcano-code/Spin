"""Post-hoc comparison of independent Faster reconstructed-method runs.

This audit is intentionally separate from immutable run evidence. It compares
the scientific core while permitting only a small, named set of runtime fields
and files to differ. It never rewrites either run's evidence.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import re
from typing import Any


_VOLATILE_EVENT_KEYS = frozenset({
    "evaluator_wall_seconds",
    "residual_evaluation_seconds",
    "derivative_evaluation_seconds",
})
_RUNTIME_METADATA_FILES = frozenset({
    "begin.json",
    "result.json",
    "evaluation-events.jsonl",
    "solver-events.jsonl",
})


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _without_volatile_event_fields(value: Any) -> Any:
    if isinstance(value, dict):
        return {
            key: _without_volatile_event_fields(item)
            for key, item in value.items()
            if key not in _VOLATILE_EVENT_KEYS
        }
    if isinstance(value, list):
        return [_without_volatile_event_fields(item) for item in value]
    return value


def _output_inventory(evidence: dict, root: Path) -> tuple[dict[str, dict], list[str]]:
    outputs = evidence.get("scientific", {}).get("outputs")
    if not isinstance(outputs, list):
        return {}, ["scientific.outputs is not a list"]
    inventory: dict[str, dict] = {}
    errors: list[str] = []
    resolved_root = root.resolve()
    for entry in outputs:
        if not isinstance(entry, dict) or not isinstance(entry.get("path"), str):
            errors.append("malformed output entry")
            continue
        path = Path(entry["path"]).resolve()
        try:
            relative = path.relative_to(resolved_root).as_posix()
        except ValueError:
            errors.append(f"output escapes declared root: {entry['path']}")
            continue
        if relative in inventory:
            errors.append(f"duplicate output path: {relative}")
            continue
        if not path.is_file():
            errors.append(f"missing output: {relative}")
            continue
        actual_size = path.stat().st_size
        actual_hash = sha256_file(path)
        if entry.get("size") != actual_size:
            errors.append(f"output size mismatch: {relative}")
        if entry.get("sha256") != actual_hash:
            errors.append(f"output SHA-256 mismatch: {relative}")
        inventory[relative] = {
            "sha256": actual_hash,
            "size": actual_size,
        }
    return inventory, errors


def _external_audit_ok(result: dict) -> bool:
    audit = result.get("external_terminal_geometry_audit")
    return bool(
        isinstance(audit, dict)
        and audit.get("status") == "PASSED"
        and audit.get("geometry_verified") is True
        and audit.get("material_domain_status") == "VERIFIED_SUBSET"
        and isinstance(audit.get("terminal_material_sha256"), str)
    )


def _sealed_identity(result: dict) -> tuple:
    sealed = result.get("sealed_artifact")
    if not isinstance(sealed, dict):
        return (None, None)
    return sealed.get("scientific_digest"), sealed.get("problem_digest")


def _build_identity(evidence: dict) -> tuple | None:
    observation = evidence.get("begin", {}).get("observation")
    if not isinstance(observation, dict):
        return None
    source = observation.get("source")
    runtime = observation.get("runtime")
    dependencies = observation.get("dependencies")
    if not isinstance(source, dict) or not isinstance(runtime, dict) or not isinstance(dependencies, dict):
        return None
    source_hash = source.get("content_sha256")
    modules = runtime.get("native_modules")
    roots = dependencies.get("roots")
    if not isinstance(source_hash, str) or re.fullmatch(r"[0-9a-f]{64}", source_hash) is None:
        return None
    if not isinstance(modules, dict) or not modules or not isinstance(roots, dict) or not roots:
        return None
    for module in modules.values():
        if not isinstance(module, dict) or not isinstance(module.get("path"), str) or \
                not isinstance(module.get("sha256"), str) or \
                re.fullmatch(r"[0-9a-f]{64}", module["sha256"]) is None or \
                not isinstance(module.get("size"), int):
            return None
    return source_hash, modules, roots


def _inventory_digest(inventory: dict[str, dict]) -> str:
    rows = [
        [name, entry["sha256"], entry["size"]]
        for name, entry in sorted(inventory.items())
    ]
    return hashlib.sha256(
        json.dumps(rows, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    ).hexdigest()


def _result_matches_file(result: dict, inventory: dict[str, dict], root: Path) -> bool:
    if "result.json" not in inventory or not isinstance(result, dict):
        return False
    try:
        return json.loads((root / "result.json").read_text(encoding="utf-8")) == result
    except (OSError, UnicodeError, json.JSONDecodeError):
        return False


def _normalize_run_paths(value: Any, root: Path) -> Any:
    if isinstance(value, dict):
        return {key: _normalize_run_paths(item, root) for key, item in value.items()}
    if isinstance(value, list):
        return [_normalize_run_paths(item, root) for item in value]
    if isinstance(value, str):
        prefix = str(root.resolve())
        if value == prefix:
            return "<RUN>"
        if value.startswith(prefix + "/"):
            return "<RUN>/" + value[len(prefix) + 1:]
    return value


def _checkpoint_scientific_state(root: Path, result: dict, evidence: dict) -> dict:
    """Validate the mixed checkpoint and bind its scientific fields to the run.

    Only wall-time values, declared run-root paths, contextual Git observations,
    and hashes of separately validated runtime journals are excluded. Unknown
    payload fields remain in the signature, so they cannot silently disappear.
    This is a comparison check, not a geometry or execution attestation.
    """
    def pairs(items):
        out = {}
        for key, value in items:
            if key in out:
                raise ValueError("duplicate checkpoint field")
            out[key] = value
        return out
    checkpoint = json.loads((root / "terminal-audit-checkpoint.json").read_text(),
                            object_pairs_hook=pairs)
    if not isinstance(checkpoint, dict):
        raise ValueError("checkpoint must be an object")
    claimed = checkpoint.pop("checkpoint_digest", None)
    encoded = json.dumps(checkpoint, sort_keys=True, separators=(",", ":"),
                         ensure_ascii=True, allow_nan=False).encode("ascii")
    if claimed != hashlib.sha256(encoded).hexdigest():
        raise ValueError("checkpoint digest mismatch")
    if (checkpoint.get("schema") != "faster-terminal-audit-checkpoint-v1" or
            checkpoint.get("status") != "PENDING_TERMINAL_AUDIT" or
            checkpoint.get("complete_run_evidence") is not False):
        raise ValueError("unsupported checkpoint state")
    payload = checkpoint["payload"]
    if not isinstance(payload, dict):
        raise ValueError("checkpoint payload must be an object")
    context = payload["run_finalization"]
    if not isinstance(context, dict) or not isinstance(context.get("native_result"), dict):
        raise ValueError("checkpoint finalization/native context must be objects")
    if (context.get("schema") != "faster-run-finalization-context-v1" or
            Path(context["output_directory"]).resolve() != root.resolve()):
        raise ValueError("checkpoint run directory mismatch")
    expected_status = result.get("solver_status", result["status"])
    if payload["solver_status"] != expected_status or context["native_result"]["status"] != expected_status:
        raise ValueError("checkpoint solver status differs from completed result")
    for key in ("trace", "evaluations", "evaluation_map_binding_audit"):
        if payload[key] != result[key]:
            raise ValueError("checkpoint scientific field differs from completed result: " + key)
    for key, value in context["native_result"].items():
        if key != "status" and key in result and value != result[key]:
            raise ValueError("checkpoint native record differs from completed result: " + key)
    if (context["locked_problem_digest"] != result["locked_problem_digest"] or
            context["classification"] != result["classification"] or
            context["expected_artifact_digest"] != _sealed_identity(result)[0] or
            context["begin"] != evidence["begin"]):
        raise ValueError("checkpoint problem/artifact/source binding mismatch")
    begin_integrity = context["begin"].get("integrity_sha256")
    if begin_integrity is not None:
        begin_bytes = json.dumps({k: v for k, v in context["begin"].items()
                                  if k != "integrity_sha256"}, sort_keys=True,
                                 separators=(",", ":"), ensure_ascii=True,
                                 allow_nan=False).encode("ascii")
        if begin_integrity != hashlib.sha256(begin_bytes).hexdigest():
            raise ValueError("checkpoint begin integrity mismatch")
    bridge_digest = context.get("bridge_scientific_digest")
    if bridge_digest is not None and (
            not isinstance(bridge_digest, str) or
            re.fullmatch(r"[0-9a-f]{64}", bridge_digest) is None or
            bridge_digest != result.get("bridge_scientific_digest")):
        raise ValueError("checkpoint bridge digest differs from completed result")
    if ("bridge_identity" in context) != ("bridge_identity_attachment" in result):
        raise ValueError("checkpoint/result bridge identity presence mismatch")
    if "bridge_identity" in context:
        bridge = context["bridge_identity"]
        required = {"operation", "native_run", "status", "seal_error",
                    "problem_digest", "artifact_scientific_digest"}
        if not isinstance(bridge, dict) or not required <= bridge.keys():
            raise ValueError("checkpoint bridge identity incomplete")
        encoded_bridge = json.dumps(bridge, sort_keys=True, separators=(",", ":"),
                                    ensure_ascii=True, allow_nan=False).encode("ascii")
        expected_attachment = {"schema": "cherchi-bridge-identity-attachment-v1",
            "status": "BOUND_IDENTITY", "identity": bridge, "digest": bridge_digest}
        if (hashlib.sha256(encoded_bridge).hexdigest() != bridge_digest or
                bridge["operation"] != "SOLVE_CHERCHI_EXTERNAL_GEOMETRY" or
                bridge["native_run"] != context["native_result"] or
                bridge["status"] != context.get("bridge_status") or
                bridge["seal_error"] != context.get("bridge_seal_error") or
                bridge["problem_digest"] != context["locked_problem_digest"] or
                bridge["artifact_scientific_digest"] != context["expected_artifact_digest"] or
                result.get("bridge_identity_attachment") != expected_attachment):
            raise ValueError("checkpoint bridge identity is not bound to native/problem/material/result")
    dependencies = checkpoint["dependencies"]
    identities = checkpoint["checker_identity"]
    if (not isinstance(identities, list) or not identities or
            any(not isinstance(path, str) or not path for path in identities) or
            len(set(identities)) != len(identities) or
            not isinstance(dependencies, list) or not dependencies or
            any(not isinstance(row, dict) or not isinstance(row.get("path"), str) or
                not isinstance(row.get("sha256"), str) for row in dependencies)):
        raise ValueError("malformed checkpoint checker dependencies")
    dependency_paths = [row["path"] for row in dependencies]
    if (len(set(dependency_paths)) != len(dependency_paths) or
            not set(identities).issubset(dependency_paths)):
        raise ValueError("checkpoint checker dependencies incomplete")
    comparable_dependencies = []
    for row in dependencies:
        path = Path(row["path"])
        if not path.is_file() or sha256_file(path) != row["sha256"]:
            raise ValueError("checkpoint dependency changed: " + str(path))
        try:
            relative = path.resolve().relative_to(root.resolve()).as_posix()
        except ValueError:
            relative = None
        if relative not in _RUNTIME_METADATA_FILES:
            comparable_dependencies.append(row)
    comparable = {**checkpoint, "dependencies": comparable_dependencies}
    # Make a by-value copy: validation above must remain against original bytes.
    comparable = json.loads(json.dumps(comparable, allow_nan=False))
    begin = comparable["payload"]["run_finalization"]["begin"]
    # Validate original seals before excluding their derived hashes. This begin
    # seal includes the diagnostic git observation removed immediately below.
    begin.pop("integrity_sha256", None)
    begin["observation"]["source"].pop("git", None)
    for dependency in begin["observation"]["dependencies"]["roots"].values():
        dependency.pop("git_head", None)
    comparable["payload"]["run_finalization"].pop("solver_wall_seconds")
    # This native digest includes duration-bearing event payloads. It was
    # bound to result.json above; the full events remain in the signature and
    # are compared after removing only the named volatile duration fields.
    comparable["payload"]["run_finalization"]["native_result"].pop("event_hash", None)
    # Native timing rows remain bound to the saved result above. Ignore only
    # elapsed seconds, preserving evaluation identities, statuses and counts.
    native = comparable["payload"]["run_finalization"]["native_result"]
    timings = native.get("evaluation_timings")
    if timings is not None:
        if not isinstance(timings, list):
            raise ValueError("native evaluation timings must be a list")
        for row in timings:
            if not isinstance(row, dict):
                raise ValueError("native evaluation timing must be an object")
            seconds = row.get("wall_seconds")
            if type(seconds) not in (int, float) or not math.isfinite(seconds) or seconds < 0:
                raise ValueError("invalid native evaluation duration")
            row.pop("wall_seconds")
    # This duplicate native copy was bound before normalization. Preserve all
    # other bridge fields, including unknown future scientific fields.
    bridge = comparable["payload"]["run_finalization"].get("bridge_identity")
    if bridge is not None:
        bridge["native_run"] = json.loads(json.dumps(native, allow_nan=False))
    # The bridge aggregate seals native_run including the duration rows. Bind
    # its original value above; compare actual native content and material
    # independently, not this duration-bearing aggregate across processes.
    comparable["payload"]["run_finalization"].pop("bridge_scientific_digest", None)
    comparable = _normalize_run_paths(_without_volatile_event_fields(comparable), root)
    return comparable


def compare_repeat_records(
    first_result,
    repeat_result,
    first_evidence,
    repeat_evidence,
    *,
    first_output_root,
    repeat_output_root,
):
    """Compare two run records and validate their evidence-listed outputs."""
    first_outputs, first_errors = _output_inventory(first_evidence, Path(first_output_root))
    repeat_outputs, repeat_errors = _output_inventory(repeat_evidence, Path(repeat_output_root))

    first_events = first_result.get("events")
    repeat_events = repeat_result.get("events")
    same_inventory = set(first_outputs) == set(repeat_outputs)
    shared_names = sorted(set(first_outputs) & set(repeat_outputs))
    different_names = [
        name for name in shared_names
        if first_outputs[name]["sha256"] != repeat_outputs[name]["sha256"]
    ]
    different_metadata = [name for name in different_names if Path(name).name in _RUNTIME_METADATA_FILES]
    different_scientific = [name for name in different_names if Path(name).name not in _RUNTIME_METADATA_FILES]
    checkpoint_name = "terminal-audit-checkpoint.json"
    checkpoint_required = checkpoint_name in first_outputs or checkpoint_name in repeat_outputs
    checkpoint_equal = not checkpoint_required
    checkpoint_errors = []
    if checkpoint_required:
        try:
            first_checkpoint = _checkpoint_scientific_state(Path(first_output_root), first_result, first_evidence)
            repeat_checkpoint = _checkpoint_scientific_state(Path(repeat_output_root), repeat_result, repeat_evidence)
            checkpoint_equal = first_checkpoint == repeat_checkpoint
        except (ValueError, TypeError, KeyError, OSError) as error:
            checkpoint_errors.append(str(error))
        # The mixed-role checkpoint has its own scientific-state gate. It is
        # not added to the general runtime filename allowlist.
        different_scientific = [name for name in different_scientific if name != checkpoint_name]
    first_binding = first_result.get("evaluation_map_binding_audit")
    repeat_binding = repeat_result.get("evaluation_map_binding_audit")
    first_audit = first_result.get("external_terminal_geometry_audit") or {}
    repeat_audit = repeat_result.get("external_terminal_geometry_audit") or {}
    first_sealed = _sealed_identity(first_result)
    repeat_sealed = _sealed_identity(repeat_result)
    first_build = _build_identity(first_evidence)
    repeat_build = _build_identity(repeat_evidence)

    checks = {
        "result_records_match_verified_files": (
            _result_matches_file(first_result, first_outputs, Path(first_output_root))
            and _result_matches_file(repeat_result, repeat_outputs, Path(repeat_output_root))
        ),
        "source_binary_and_dependencies_match": (
            first_build is not None and first_build == repeat_build
        ),
        "checkpoint_scientific_states_equal_and_bound": checkpoint_equal,
        "evidence_outputs_current_and_confined": not first_errors and not repeat_errors,
        "output_inventory_matches": same_inventory,
        "scientific_output_bytes_equal": same_inventory and not different_scientific,
        "solver_events_equal_without_timing": (
            isinstance(first_events, list)
            and isinstance(repeat_events, list)
            and _without_volatile_event_fields(first_events)
            == _without_volatile_event_fields(repeat_events)
        ),
        "evidence_scientific_events_equal": (
            first_evidence.get("scientific", {}).get("events")
            == repeat_evidence.get("scientific", {}).get("events")
        ),
        "trace_and_trace_hash_equal": (
            first_result.get("trace") == repeat_result.get("trace")
            and first_result.get("trace_hash") == repeat_result.get("trace_hash")
        ),
        "locked_problem_and_solver_counts_equal": all(
            first_result.get(key) == repeat_result.get(key)
            for key in (
                "locked_problem_digest",
                "function_evaluations",
                "restart_count",
                "restart_function_evaluations",
            )
        ),
        "run_terminal_states_equal_and_converged": (
            first_result.get("status") == repeat_result.get("status") == "CONVERGED"
            and first_result.get("state") == repeat_result.get("state")
        ),
        "evaluation_map_bindings_equal_and_consistent": (
            first_binding == repeat_binding
            and isinstance(first_binding, dict)
            and first_binding.get("status") == "CONSISTENT"
        ),
        "external_geometry_audits_pass_and_match": (
            _external_audit_ok(first_result)
            and _external_audit_ok(repeat_result)
            and first_audit.get("terminal_material_sha256")
            == repeat_audit.get("terminal_material_sha256")
        ),
        "sealed_geometry_and_moments_identical": (
            first_sealed == repeat_sealed
            and first_sealed[0] is not None
            and first_sealed[1] is not None
        ),
        "classification_and_initialization_match": (
            first_result.get("classification") == repeat_result.get("classification")
            and first_result.get("initialization_mode") == repeat_result.get("initialization_mode")
            and first_result.get("deterministic_seed") == repeat_result.get("deterministic_seed")
        ),
    }
    report = {
        "schema": "faster-independent-repeat-comparison-v2",
        "status": "SCIENTIFIC_CORE_IDENTICAL" if all(checks.values()) else "MISMATCH",
        "checks": checks,
        "errors": {"first": first_errors, "repeat": repeat_errors},
        "checkpoint_comparison_errors": checkpoint_errors,
        "volatile_event_fields_ignored": sorted(_VOLATILE_EVENT_KEYS),
        "checkpoint_runtime_normalization": {
            "derived_bound_digests": ["begin.integrity_sha256", "native_result.event_hash",
                                      "bridge_scientific_digest"],
            "durations": ["solver_wall_seconds", "native_result.evaluation_timings[].wall_seconds"],
            "scope": "validated_original_bindings_before_comparing_nonruntime_content",
        },
        "outputs": {
            "first_count": len(first_outputs),
            "repeat_count": len(repeat_outputs),
            "first_verified_inventory_sha256": _inventory_digest(first_outputs),
            "repeat_verified_inventory_sha256": _inventory_digest(repeat_outputs),
            "identical_scientific_file_count": len(shared_names) - len(different_names),
            "different_runtime_metadata_files": [Path(name).name for name in different_metadata],
            "different_runtime_metadata_paths": different_metadata,
            "different_scientific_files": different_scientific,
            "different_mixed_checkpoint_files": [name for name in different_names if name == checkpoint_name],
            "missing_from_first": sorted(set(repeat_outputs) - set(first_outputs)),
            "missing_from_repeat": sorted(set(first_outputs) - set(repeat_outputs)),
            "runtime_metadata_allowlist": sorted(_RUNTIME_METADATA_FILES),
        },
        "comparison_scope": "posthoc_scientific_repeat_audit_not_paper_completion_or_native_determinism_digest",
    }
    canonical = json.dumps(report, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    report["comparison_digest"] = hashlib.sha256(canonical).hexdigest()
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--first-result", type=Path, required=True)
    parser.add_argument("--repeat-result", type=Path, required=True)
    parser.add_argument("--first-evidence", type=Path, required=True)
    parser.add_argument("--repeat-evidence", type=Path, required=True)
    parser.add_argument("--first-output-root", type=Path, required=True)
    parser.add_argument("--repeat-output-root", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = compare_repeat_records(
        json.loads(args.first_result.read_text()),
        json.loads(args.repeat_result.read_text()),
        json.loads(args.first_evidence.read_text()),
        json.loads(args.repeat_evidence.read_text()),
        first_output_root=args.first_output_root,
        repeat_output_root=args.repeat_output_root,
    )
    encoded = json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.output:
        args.output.write_text(encoded)
    else:
        print(encoded, end="")
    return 0 if report["status"] == "SCIENTIFIC_CORE_IDENTICAL" else 1


if __name__ == "__main__":
    raise SystemExit(main())
