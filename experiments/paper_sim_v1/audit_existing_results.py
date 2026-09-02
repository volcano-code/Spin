#!/usr/bin/env python3
"""Re-evaluate experiment manifests against the current source/executables."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from make_it_stand.io import write_json_atomic
from make_it_stand.paper_evidence import audit_manifest_freshness


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", nargs="+", type=Path)
    parser.add_argument("--write", action="store_true")
    arguments = parser.parse_args()
    repository_root = Path(__file__).resolve().parents[2]
    stale = 0
    summaries = []
    for raw_path in arguments.manifest:
        path = raw_path.expanduser().resolve()
        payload = json.loads(path.read_text(encoding="utf-8"))
        audited = audit_manifest_freshness(payload, repository_root)
        if audited["freshness"]["status"] == "STALE_RESULT":
            stale += 1
        if arguments.write:
            write_json_atomic(path, audited)
        summaries.append({
            "manifest": str(path),
            "freshness": audited["freshness"],
        })
    print(json.dumps(summaries, indent=2, sort_keys=True))
    return 2 if stale else 0


if __name__ == "__main__":
    raise SystemExit(main())
