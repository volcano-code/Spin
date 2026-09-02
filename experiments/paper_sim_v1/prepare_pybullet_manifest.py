#!/usr/bin/env python3
"""Create a PyBullet manifest only from current, post-flatten sealed artifacts."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

from make_it_stand.paper_evidence import (
    audit_manifest_freshness,
    capture_execution_fingerprint,
)
from make_it_stand.sealed_artifact import load_sealed_pybullet_entry


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case-run", action="append", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    repository_root = Path(__file__).resolve().parents[2]
    cases = []
    for case_root in arguments.case_run:
        case_root = case_root.resolve()
        manifest_path = case_root / "manifest.json"
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        audited = audit_manifest_freshness(manifest, repository_root)
        if audited["freshness"]["status"] != "CURRENT":
            raise ValueError(
                f"{manifest.get('case', case_root.name)}: optimizer result is "
                f"STALE_RESULT ({', '.join(audited['freshness']['mismatches'])})"
            )
        run = manifest["runs"][0]
        if run["completion"] != "stable" or not run["succeeded"]:
            raise ValueError(f"{manifest['case']}: optimizer result is not stable")
        sealed_path = case_root / "repeat-1" / "sealed" / "sealed-artifact.json"
        expected_sealed_hash = run["output_sha256"].get("sealed_artifact")
        if not expected_sealed_hash or _sha256(sealed_path) != expected_sealed_hash:
            raise ValueError(f"{manifest['case']}: sealed artifact hash changed")
        sealed_payload = json.loads(sealed_path.read_text(encoding="utf-8"))
        entry = load_sealed_pybullet_entry(sealed_path)
        entry["name"] = manifest["case"]
        entry["source_profile"] = manifest["profile"]
        entry["evidence_class"] = manifest["classification"]
        entry["manufacturing_ready"] = bool(
            sealed_payload["manufacturing_ready"]
        )
        entry["sealing_failures"] = list(sealed_payload["failures"])
        cases.append(entry)
    payload = {
        "schema": "make-it-stand-paper-sim-pybullet-input-v2",
        "schema_version": 2,
        "input_contract": "current-post-flatten-sealed-artifact-only",
        "cases": cases,
    }
    payload["freshness"] = {"status": "CURRENT", "mismatches": []}
    payload["execution_fingerprint"] = capture_execution_fingerprint(
        repository_root,
        source_paths=[
            repository_root / "python" / "make_it_stand" / "sealed_artifact.py",
            Path(__file__).resolve(),
            *(root.resolve() / "manifest.json" for root in arguments.case_run),
            *(
                root.resolve() / "repeat-1" / "sealed" / "sealed-artifact.json"
                for root in arguments.case_run
            ),
        ],
        executable_paths=[Path(sys.executable).resolve()],
        configuration={"case_names": [item["name"] for item in cases]},
    )
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = arguments.output.with_suffix(arguments.output.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(arguments.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
