#!/usr/bin/env python3
"""Run the license-safe synthetic R10 ablation and robustness matrix.

The driver creates each immutable manifest before loading the model, then
emits one canonical record per repeat and a canonical aggregate.  The only
geometry distributed with the protocol is an original CC0 fixture.  The
Make It Stand beta archive is neither read nor copied by this program.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Any, Callable, Mapping, Sequence

from make_it_stand.model import BalanceModel, ModelConfig, RunConfig
from make_it_stand.r10 import (
    CANONICALIZATION,
    MANIFEST_SCHEMA,
    RunRecorder,
    aggregate_run_records,
    bind_artifact,
    validate_manifest,
    write_canonical_json,
)


@dataclass(frozen=True)
class ExperimentCase:
    experiment_id: str
    variant: str
    resolution: int
    handle_layout: str
    enable_carving: bool
    enable_deformation: bool
    enable_laplacian: bool
    constraint: str = "standing"


@dataclass(frozen=True)
class RepositoryState:
    commit_sha: str
    dirty_worktree: bool


METRIC_DEFINITIONS: tuple[dict[str, str], ...] = (
    {
        "name": "projected_target_energy",
        "unit": "normalized_objective",
        "direction": "lower_is_better",
        "definition": "Final projected balance-target energy reported by the R10 diagnostics.",
        "source": "local_measurement",
    },
    {
        "name": "removed_mass_fraction",
        "unit": "fraction",
        "direction": "higher_is_better",
        "definition": "Initial mass minus final mass, divided by initial mass.",
        "source": "local_measurement",
    },
    {
        "name": "deformation_energy",
        "unit": "normalized_objective",
        "direction": "lower_is_better",
        "definition": "Final unscaled Laplacian shape energy reported by the R10 diagnostics.",
        "source": "local_measurement",
    },
    {
        "name": "normalized_rms_displacement",
        "unit": "normalized_length",
        "direction": "lower_is_better",
        "definition": "RMS source-vertex displacement divided by source bounding-box diagonal.",
        "source": "local_measurement",
    },
    {
        "name": "normalized_maximum_displacement",
        "unit": "normalized_length",
        "direction": "lower_is_better",
        "definition": "Maximum source-vertex displacement divided by source bounding-box diagonal.",
        "source": "local_measurement",
    },
    {
        "name": "committed_event_count",
        "unit": "events",
        "direction": "descriptive",
        "definition": "Count of committed R10 geometry states, including the initial state.",
        "source": "local_measurement",
    },
    {
        "name": "carving_iteration_count",
        "unit": "iterations",
        "direction": "descriptive",
        "definition": "Count of recorded carving iterations under the modern controller semantics.",
        "source": "local_measurement",
    },
    {
        "name": "deformation_iteration_count",
        "unit": "iterations",
        "direction": "descriptive",
        "definition": "Count of recorded deformation line-search iterations.",
        "source": "local_measurement",
    },
    {
        "name": "algorithm_success",
        "unit": "indicator",
        "direction": "higher_is_better",
        "definition": "One when the native result classifies the optimization outcome as successful.",
        "source": "local_measurement",
    },
    {
        "name": "targets_reached",
        "unit": "indicator",
        "direction": "higher_is_better",
        "definition": "One when all declared final balance targets are reached.",
        "source": "local_measurement",
    },
    {
        "name": "manufacturing_ready",
        "unit": "indicator",
        "direction": "higher_is_better",
        "definition": "One when the sealed manufacturing audit reports ready.",
        "source": "local_measurement",
    },
    {
        "name": "measured_minimum_wall_clearance",
        "unit": "model_length",
        "direction": "higher_is_better",
        "definition": "Minimum audited outer-to-inner clearance when a cavity exists and the audit completes.",
        "source": "local_measurement",
    },
)

_METRIC_SEMANTICS = {
    item["name"]: item["definition"] for item in METRIC_DEFINITIONS
}


def default_cases() -> tuple[ExperimentCase, ...]:
    """Return the fixed R10 matrix; no case depends on a third-party asset."""

    return (
        ExperimentCase("synthetic-full-r16-h0", "full", 16, "low", True, True, True),
        ExperimentCase(
            "synthetic-carving-only-r16-h0",
            "carving-only",
            16,
            "low",
            True,
            False,
            True,
        ),
        ExperimentCase(
            "synthetic-deformation-only-r16-h0",
            "deformation-only",
            16,
            "low",
            False,
            True,
            True,
        ),
        ExperimentCase(
            "synthetic-no-laplacian-r16-h0",
            "no-laplacian",
            16,
            "low",
            True,
            True,
            False,
        ),
        ExperimentCase("synthetic-full-r12-h0", "resolution-r12", 12, "low", True, True, True),
        ExperimentCase("synthetic-full-r20-h0", "resolution-r20", 20, "low", True, True, True),
        ExperimentCase("synthetic-full-r16-h1", "handle-high", 16, "high", True, True, True),
        ExperimentCase(
            "synthetic-suspension-r16-h0",
            "suspension",
            16,
            "low",
            True,
            True,
            True,
            constraint="suspension",
        ),
    )


def repository_state(repository: Path) -> RepositoryState:
    """Read Git provenance once, before the evidence directory is created."""

    commit = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=repository,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    status = subprocess.run(
        ["git", "status", "--porcelain", "--untracked-files=no"],
        cwd=repository,
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    return RepositoryState(commit, bool(status.strip()))


def _copy_fixture(repository: Path, experiment_root: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    asset_directory = repository / "experiments" / "r10" / "assets"
    input_directory = experiment_root / "inputs"
    input_directory.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(asset_directory / "synthetic_overhang.obj", input_directory / "synthetic_overhang.obj")
    shutil.copyfile(asset_directory / "LICENSE.txt", input_directory / "LICENSE.txt")
    model = bind_artifact(
        experiment_root,
        "inputs/synthetic_overhang.obj",
        artifact_id="synthetic-overhang",
        role="model_input",
        media_type="model/obj",
        description="Original watertight CC0 overhang fixture.",
    )
    license_binding = bind_artifact(
        experiment_root,
        "inputs/LICENSE.txt",
        artifact_id="synthetic-overhang-license",
        role="asset_license",
        media_type="text/plain",
    )
    return model, license_binding


def build_manifest(
    case: ExperimentCase,
    *,
    experiment_root: Path,
    repository: Path,
    provenance: RepositoryState,
    repeats: int,
) -> dict[str, Any]:
    """Materialize and validate a manifest before native execution."""

    model_binding, license_binding = _copy_fixture(repository, experiment_root)
    manifest: dict[str, Any] = {
        "schema": MANIFEST_SCHEMA,
        "schema_version": 1,
        "canonicalization": CANONICALIZATION,
        "experiment_id": case.experiment_id,
        "title": f"License-safe synthetic R10 study: {case.variant}",
        "claim": {
            "level": "local_robustness",
            "statement": "This run measures the modern implementation on one original synthetic fixture.",
            "exact_reproduction_claim": False,
            "limitations": [
                "A single synthetic fixture cannot establish object-class generality.",
                "Controller events are not the beta-2013 iteration unit.",
                "The result is not an exact reproduction of a published artifact.",
                "This computational matrix does not constitute physical print validation.",
            ],
        },
        "paper_reference": {
            "citation_key": "prevost2013makeitstand",
            "scope": "Algorithmic ablation context only; no paper value is pooled with local measurements.",
            "reported_metrics": [],
        },
        "protocol": {
            "variant": case.variant,
            "repeats": repeats,
            "seed_policy": {
                "kind": "fixed_sequence",
                "first_seed": 41000,
                "use": "audit identifiers and timing repeats; the optimizer is deterministic",
            },
            "preprocessing": {
                "fixture": "synthetic-overhang-v1",
                "normalization": "native facade normalization",
                "voxel_resolution": case.resolution,
                "minimum_wall_thickness_model_units": 0.06,
            },
            "algorithm": {
                "carving": case.enable_carving,
                "deformation": case.enable_deformation,
                "laplacian": case.enable_laplacian,
                "constraint": case.constraint,
                "handle_layout": case.handle_layout,
                "shape_scale_lambda": 1.0,
                "maximum_carving_blocks": 8,
                "maximum_deformation_iterations_per_block": 200,
            },
            "metric_definitions": list(METRIC_DEFINITIONS),
        },
        "provenance": {
            "repository_commit_sha": provenance.commit_sha,
            "dirty_worktree": provenance.dirty_worktree,
            "implementation_profile": "release-headless-r10-synthetic-v1",
            "resources": [license_binding],
            "deviations": [
                "The fixture is repository-authored rather than a model from the paper.",
                "The adaptive local voxel grid differs from the publication fixed cubic grid.",
                "The modern guarded carving and manufacturing audits postdate the beta implementation.",
            ],
        },
        "inputs": [model_binding],
    }
    validate_manifest(manifest)
    write_canonical_json(experiment_root / "manifest.json", manifest)
    return manifest


def _configure_model(model: BalanceModel, case: ExperimentCase) -> None:
    model.set_gravity((0.0, -1.0, 0.0))
    if case.constraint == "standing":
        model.set_support_polygon(
            (
                (0.0, 0.0, 0.0),
                (0.39, 0.0, 0.0),
                (0.39, 0.0, 0.35),
                (0.0, 0.0, 0.35),
            ),
            0.0,
        )
    elif case.constraint == "suspension":
        model.add_suspension_constraint((0.12, 0.94, 0.175), math.radians(5.0))
    else:
        raise ValueError(f"unsupported constraint profile: {case.constraint}")

    pivots = {
        "low": (0.18, 0.35, 0.175),
        "high": (0.18, 0.82, 0.175),
    }
    model.add_handle(pivots[case.handle_layout])


def _number(value: Any, field: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"diagnostics field {field!r} is not numeric")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"diagnostics field {field!r} is not finite")
    return result


def extract_metrics(result: Any, diagnostics: Mapping[str, Any]) -> list[tuple[str, float]]:
    """Extract only measured values; absent wall clearance stays absent."""

    r10 = diagnostics["r10_metrics"]
    targets = diagnostics["targets"]
    manufacturing = diagnostics["manufacturing"]
    metrics = [
        ("projected_target_energy", _number(r10["final_projected_target_energy"], "final_projected_target_energy")),
        ("removed_mass_fraction", _number(r10["removed_mass_fraction"], "removed_mass_fraction")),
        ("deformation_energy", _number(r10["final_unscaled_shape_energy"], "final_unscaled_shape_energy")),
        ("normalized_rms_displacement", _number(r10["normalized_rms_vertex_displacement"], "normalized_rms_vertex_displacement")),
        ("normalized_maximum_displacement", _number(r10["normalized_maximum_vertex_displacement"], "normalized_maximum_vertex_displacement")),
        ("committed_event_count", _number(r10["committed_event_count"], "committed_event_count")),
        ("carving_iteration_count", float(len(diagnostics["carving_iterations"]))),
        ("deformation_iteration_count", float(len(diagnostics["deformation_iterations"]))),
        ("algorithm_success", 1.0 if result.succeeded else 0.0),
        ("targets_reached", 1.0 if targets["all_reached"] else 0.0),
        ("manufacturing_ready", 1.0 if result.manufacturing_ready else 0.0),
    ]
    wall = manufacturing.get("wall_clearance", {})
    clearance = wall.get("measured_minimum_clearance")
    if clearance is not None:
        metrics.append(("measured_minimum_wall_clearance", _number(clearance, "measured_minimum_clearance")))
    return metrics


def run_case(
    case: ExperimentCase,
    *,
    output_root: Path,
    repository: Path,
    provenance: RepositoryState,
    repeats: int,
    build: Mapping[str, Any],
    model_factory: Callable[..., Any] = BalanceModel,
) -> dict[str, Any]:
    """Run one case and return its canonical aggregate document."""

    experiment_root = output_root / case.experiment_id
    manifest = build_manifest(
        case,
        experiment_root=experiment_root,
        repository=repository,
        provenance=provenance,
        repeats=repeats,
    )
    records: list[dict[str, Any]] = []
    for repeat_index in range(repeats):
        run_id = f"{case.experiment_id}-r{repeat_index:02d}"
        recorder = RunRecorder(
            manifest,
            run_id=run_id,
            repeat_index=repeat_index,
            seed=41000 + repeat_index,
            build=build,
            timing_scope="model-load-through-artifact-export",
            warmup_runs=0,
        )
        try:
            model_path = experiment_root / "inputs" / "synthetic_overhang.obj"
            with recorder.phase("model_setup"):
                model = model_factory(
                    os.fspath(model_path),
                    config=ModelConfig(
                        path=os.fspath(model_path),
                        voxel_resolution=case.resolution,
                        minimum_wall_thickness=0.06,
                    ),
                )
                _configure_model(model, case)
            run_config = RunConfig(
                maximum_carving_blocks=8,
                maximum_deformation_iterations_per_block=200,
                shape_scale_lambda=1.0,
                enable_carving=case.enable_carving,
                enable_deformation=case.enable_deformation,
                enable_laplacian=case.enable_laplacian,
            )
            with recorder.phase("optimization"):
                result = model.optimize(run_config)
            diagnostics = result.diagnostics
            for name, value in extract_metrics(result, diagnostics):
                recorder.add_local_metric(name, value, semantics=_METRIC_SEMANTICS[name])
            recorder.add_note(f"native optimizer status: {result.status}")
            run_root = experiment_root / "runs" / f"r{repeat_index:02d}"
            run_root.mkdir(parents=True, exist_ok=True)
            with recorder.phase("artifact_export"):
                result.save_diagnostics(run_root / "diagnostics.json")
                result.save_outer_mesh(run_root / "outer.obj")
                if result.manufacturing_ready:
                    result.save_printable_mesh(run_root / "printable.stl")
            for relative, artifact_id, role, media_type in (
                (f"runs/r{repeat_index:02d}/diagnostics.json", f"diagnostics-r{repeat_index:02d}", "diagnostics", "application/json"),
                (f"runs/r{repeat_index:02d}/outer.obj", f"outer-r{repeat_index:02d}", "result_mesh", "model/obj"),
            ):
                recorder.add_artifact(
                    bind_artifact(
                        experiment_root,
                        relative,
                        artifact_id=artifact_id,
                        role=role,
                        media_type=media_type,
                    )
                )
            if result.manufacturing_ready:
                recorder.add_artifact(
                    bind_artifact(
                        experiment_root,
                        f"runs/r{repeat_index:02d}/printable.stl",
                        artifact_id=f"printable-r{repeat_index:02d}",
                        role="manufacturing_mesh",
                        media_type="model/stl",
                    )
                )
            record = recorder.finish(status="succeeded")
        except Exception as error:  # keep the batch auditable without leaking paths
            recorder.add_note(f"execution failed with {type(error).__name__}")
            record = recorder.finish(status="failed")
            print(f"{case.experiment_id} repeat {repeat_index} failed: {error}", file=sys.stderr)
        write_canonical_json(
            experiment_root / "records" / f"r{repeat_index:02d}.json",
            record,
        )
        records.append(record)
    aggregate = aggregate_run_records(manifest, records)
    write_canonical_json(experiment_root / "aggregate.json", aggregate)
    return aggregate


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--case", action="append", default=[], help="experiment_id to run; repeatable")
    parser.add_argument("--repeats", type=int, default=1)
    parser.add_argument("--build-type", required=True)
    parser.add_argument("--compiler", required=True)
    parser.add_argument("--compiler-version", required=True)
    parser.add_argument("--build-flag", action="append", default=[])
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    arguments = _parser().parse_args(argv)
    if arguments.repeats < 1:
        raise SystemExit("--repeats must be at least one")
    repository = Path(__file__).resolve().parents[1]
    # A build-tree extension normally lives beside this script's package root
    # (for example build/python/python).  Make explicit PYTHONPATH entries work
    # even when the source package itself was imported first.
    import make_it_stand

    for entry in sys.path:
        candidate = Path(entry) / "make_it_stand"
        if candidate.is_dir():
            text = os.fspath(candidate)
            if text not in make_it_stand.__path__:
                make_it_stand.__path__.append(text)
    provenance = repository_state(repository)
    selected = set(arguments.case)
    cases = tuple(case for case in default_cases() if not selected or case.experiment_id in selected)
    unknown = selected - {case.experiment_id for case in cases}
    if unknown:
        raise SystemExit(f"unknown --case value(s): {', '.join(sorted(unknown))}")
    build = {
        "build_type": arguments.build_type,
        "compiler": arguments.compiler,
        "compiler_version": arguments.compiler_version,
        "commit_sha": provenance.commit_sha,
        "dirty_worktree": provenance.dirty_worktree,
        "flags": list(arguments.build_flag),
    }
    summaries = []
    for case in cases:
        aggregate = run_case(
            case,
            output_root=arguments.output_root,
            repository=repository,
            provenance=provenance,
            repeats=arguments.repeats,
            build=build,
        )
        summaries.append(
            {
                "experiment_id": case.experiment_id,
                "successful_repeat_count": aggregate["successful_repeat_count"],
                "repeat_count": aggregate["repeat_count"],
            }
        )
    print(json.dumps(summaries, sort_keys=True))
    return 0 if all(item["successful_repeat_count"] == item["repeat_count"] for item in summaries) else 2


if __name__ == "__main__":
    raise SystemExit(main())
