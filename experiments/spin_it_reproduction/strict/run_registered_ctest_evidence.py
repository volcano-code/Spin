#!/usr/bin/env python3
"""Capture complete CTest evidence; inventory-only output never grants verification."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import time

from make_it_stand.verification.run_evidence import begin_run, check_freshness, finish_run
from make_it_stand.verification.spin_common import canonical, parse


ROOT = Path(__file__).resolve().parents[3]
SCOPE = "REGISTERED_CTEST_COMPLETE_EXECUTION_GATE_ONLY_NOT_PAPER_CASE_EVIDENCE"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write(path: Path, value: object) -> None:
    with path.open("x", encoding="ascii") as stream:
        stream.write(canonical(value) + "\n")


def _inventory_ids(inventory_text: str) -> list[int]:
    entries = [int(item) for item in re.findall(r"^\s*Test\s+#(\d+): .+$", inventory_text, flags=re.MULTILINE)]
    total_match = re.search(r"^Total Tests: (\d+)\s*$", inventory_text, flags=re.MULTILINE)
    if total_match is None:
        raise ValueError("CTest inventory lacks Total Tests")
    total = int(total_match.group(1))
    if entries != list(range(1, total + 1)):
        raise ValueError("CTest inventory does not enumerate every registered test")
    return entries


def parse_complete_ctest_output(inventory_text: str, execution_text: str) -> dict[str, object]:
    """Require every numbered inventory entry to have one passing CTest receipt."""
    inventory_ids = _inventory_ids(inventory_text)
    total = len(inventory_ids)
    if "Not Run" in execution_text:
        raise ValueError("CTest execution contains Not Run")
    receipts = re.findall(
        r"^\s*(\d+)/(\d+)\s+Test\s+#(\d+):.*\s+(Passed|Failed)\b.*$",
        execution_text,
        flags=re.MULTILINE,
    )
    if len(receipts) != total:
        raise ValueError("CTest execution receipt count differs from inventory")
    executed_ids = [int(test_id) for _, _, test_id, _ in receipts]
    if executed_ids != inventory_ids or any(int(denominator) != total for _, denominator, _, _ in receipts):
        raise ValueError("CTest execution does not cover the full inventory")
    if any(status != "Passed" for _, _, _, status in receipts):
        raise ValueError("CTest execution has a non-passing receipt")
    summary = re.search(r"^100% tests passed, 0 tests failed out of (\d+)\s*$", execution_text,
                        flags=re.MULTILINE)
    if summary is None or int(summary.group(1)) != total:
        raise ValueError("CTest execution summary disagrees with inventory")
    return {"inventory_total": total, "executed_test_ids": executed_ids}


def run(*, build_directory: Path, output_directory: Path, manifest_only: bool) -> dict[str, object]:
    """Record the CTest inventory and, unless requested otherwise, one full execution."""
    build = build_directory.resolve(strict=True)
    output = output_directory.absolute()
    output.mkdir(parents=False, exist_ok=False)
    ctest_raw = shutil.which("ctest")
    if ctest_raw is None:
        raise RuntimeError("CTest executable is unavailable")
    ctest = Path(ctest_raw).resolve(strict=True)
    cache = (build / "CMakeCache.txt").resolve(strict=True)
    testfiles = sorted(path.resolve(strict=True) for path in build.rglob("CTestTestfile.cmake"))
    if not testfiles:
        raise RuntimeError("CTest test-file inventory is unavailable")
    started, tick, begin, events = datetime.now(timezone.utc).isoformat(), time.perf_counter(), None, []
    try:
        inputs = parse(canonical({
            "schema": "registered-ctest-inputs-v2", "scope": SCOPE,
            "build_directory": str(build), "ctest": {"path": str(ctest), "sha256": _sha256(ctest)},
            "manifest_command": [str(ctest), "--test-dir", str(build), "-N"],
            # The isolated HiGHS dependency contract is a deliberately cold
            # build and can exceed its per-test property on this host. Keep
            # the registered run sequential, but give CTest the explicit
            # evidence-run budget needed for a complete receipt.
            "execution_command": [str(ctest), "--test-dir", str(build), "--timeout", "420", "--output-on-failure", "-j1"],
            "mode": "MANIFEST_ONLY" if manifest_only else "COMPLETE_EXECUTION",
        }))
        input_path = output / "inputs.json"
        _write(input_path, inputs)
        begin = begin_run(
            repository_root=ROOT,
            runner_files=[Path(__file__)],
            checker_files=[
                ROOT / "python/make_it_stand/verification/run_evidence.py",
                ROOT / "CMakeLists.txt",
            ],
            asset_files=[input_path, ctest],
            build_configuration_files=[cache, *testfiles],
            parameters=inputs,
            seed_protocol={"seed": 0, "protocol": "deterministic sequential registered CTest execution"},
            dependency_roots={"cmake_modules": ROOT / "cmake"}, native_modules=[], pure_python=True,
        )
        _write(output / "begin.json", begin.as_dict())
        manifest = subprocess.run(inputs["manifest_command"], capture_output=True, text=True, check=False)
        inventory = output / "ctest-inventory.txt"
        inventory.write_text(manifest.stdout + manifest.stderr, encoding="utf-8")
        if manifest.returncode:
            raise RuntimeError("CTest inventory returned nonzero")
        inventory_ids = _inventory_ids(inventory.read_text(encoding="utf-8"))
        scientific = {
            "schema": "registered-ctest-scientific-v2", "scope": SCOPE,
            "mode": inputs["mode"], "inventory": {
                "path": str(inventory), "sha256": _sha256(inventory), "total": len(inventory_ids),
            },
        }
        if manifest_only:
            events.append({"iteration": 0, "status": "ok", "data": {
                "kind": "REGISTERED_CTEST_INVENTORY", "mode": "MANIFEST_ONLY",
                "inventory_total": len(inventory_ids), "verification_effect": "NONE_REQUIRES_COMPLETE_EXECUTION",
            }})
        else:
            completed = subprocess.run(inputs["execution_command"], capture_output=True, text=True, check=False)
            transcript = output / "ctest-execution.txt"
            transcript.write_text(completed.stdout + completed.stderr, encoding="utf-8")
            try:
                if completed.returncode:
                    raise RuntimeError("CTest execution returned nonzero")
                execution_summary = parse_complete_ctest_output(
                    inventory.read_text(encoding="utf-8"), transcript.read_text(encoding="utf-8"),
                )
            except (RuntimeError, ValueError) as exc:
                events.append({"iteration": 0, "status": "failed", "data": {
                    "kind": "REGISTERED_CTEST_EXECUTION", "mode": "COMPLETE_EXECUTION", "reason": str(exc),
                }})
                scientific["execution"] = {"path": str(transcript), "sha256": _sha256(transcript),
                                           "returncode": completed.returncode, "accepted": False, "reason": str(exc)}
            else:
                events.append({"iteration": 0, "status": "ok", "data": {
                    "kind": "REGISTERED_CTEST_EXECUTION", "mode": "COMPLETE_EXECUTION",
                    "inventory_total": execution_summary["inventory_total"],
                    "executed_test_ids_sha256": hashlib.sha256(
                        canonical(execution_summary["executed_test_ids"]).encode("ascii")
                    ).hexdigest(),
                }})
                scientific["execution"] = {"path": str(transcript), "sha256": _sha256(transcript),
                                           "returncode": completed.returncode, "accepted": True,
                                           "summary": execution_summary}
        _write(output / "scientific.json", scientific)
        evidence = finish_run(
            begin, events=events,
            output_files=sorted(path for path in output.rglob("*")
                                if path.is_file() and path.name not in {"begin.json", "evidence.json"}),
            measurements={"started_at": started, "elapsed_seconds": time.perf_counter() - tick},
        )
        _write(output / "evidence.json", evidence.as_dict())
        if check_freshness(evidence).status != "CURRENT":
            raise RuntimeError("new CTest evidence is unexpectedly stale")
        if not manifest_only and not scientific["execution"]["accepted"]:
            raise RuntimeError("CTest execution is incomplete, failed, or contains Not Run")
        return evidence.as_dict()
    except Exception:
        if begin is not None and not (output / "evidence.json").exists():
            _write(output / "failure.json", {"schema": "registered-ctest-failure-v1", "scope": SCOPE})
        raise


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-directory", type=Path, required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    parser.add_argument("--dry-run-manifest-only", action="store_true")
    args = parser.parse_args()
    print(json.dumps(run(build_directory=args.build_directory, output_directory=args.output_directory,
                         manifest_only=args.dry_run_manifest_only), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
