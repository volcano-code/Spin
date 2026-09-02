#!/usr/bin/env python3
"""Run the three auditable alternative-pose reachability classifications."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
from typing import Any, Mapping

from make_it_stand.io import read_json
from make_it_stand.model import ModelConfig
from make_it_stand.official_beta import load_released_beta_archive
from make_it_stand.paper_compat import PaperCompatibilityModel
from make_it_stand.reachability import ReachabilityResult, run_manifest


DEFAULT_BETA_ROOT = Path(
    "/home/lin/xvxinyao/reference-only/make-it-stand/"
    "beta-extracted/make-it-stand_src_beta"
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _configure_alternative_constraint(native: Any, case: Mapping[str, Any]) -> None:
    native.clear_constraints()
    native.add_standing_constraint_from_mesh(
        case["gravity"],
        float(case["contact_threshold"]),
        float(case["gamma_radians"]),
        float(case["minimum_stability_inset"]),
        bool(case["support_plane_at_band_top"]),
    )


def main() -> int:
    script_directory = Path(__file__).resolve().parent
    repository = script_directory.parents[2]
    parser = argparse.ArgumentParser()
    parser.add_argument("--repository", type=Path, default=repository)
    parser.add_argument("--beta-root", type=Path, default=DEFAULT_BETA_ROOT)
    parser.add_argument("--manifest", type=Path, default=script_directory / "cases.json")
    parser.add_argument("--output", type=Path, default=script_directory / "results")
    arguments = parser.parse_args()

    repository = arguments.repository.expanduser().resolve()
    beta_root = arguments.beta_root.expanduser().resolve()
    manifest = read_json(arguments.manifest)

    def evaluate(case: Mapping[str, Any]) -> ReachabilityResult:
        kind = case["source_kind"]
        source = Path(str(case["source"]))
        if kind == "RELEASED_BETA_ARCHIVE":
            archive = load_released_beta_archive(
                beta_root / source,
                beta_root=beta_root,
            )
            model = archive.build_model()
            _configure_alternative_constraint(model.native, case)
            return ReachabilityResult.from_native(
                model.native.analyze_carving_reachability(0)
            )
        if kind == "RECONSTRUCTED_OBJ":
            mesh = repository / source
            expected_hash = str(case["source_sha256"])
            observed_hash = _sha256(mesh)
            if observed_hash != expected_hash:
                raise RuntimeError(
                    f"{case['id']}: reconstructed asset hash mismatch: "
                    f"expected {expected_hash}, observed {observed_hash}"
                )
            model = PaperCompatibilityModel(
                mesh,
                config=ModelConfig(
                    path=str(mesh),
                    voxel_resolution=int(case["voxel_resolution"]),
                ),
            )
            _configure_alternative_constraint(model.native, case)
            return ReachabilityResult.from_native(
                model.native.analyze_carving_reachability(0)
            )
        raise ValueError(f"{case['id']}: unsupported source_kind {kind!r}")

    summary = run_manifest(manifest, arguments.output, evaluate)
    for case in summary["cases"]:
        print(
            f"{case['id']}: {case['result']['status']} -> "
            f"{case['classification']}"
        )
    return 0 if summary["classification_complete"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
