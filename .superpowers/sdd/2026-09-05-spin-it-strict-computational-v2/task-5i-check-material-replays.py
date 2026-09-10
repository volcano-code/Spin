"""Reload two actual synthetic runner outputs; never optimizes or grants release."""
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

from make_it_stand.verification.run_evidence import PaperRunEvidenceV2, check_freshness
from make_it_stand.verification.spin_bundle import load_spin_artifact
from make_it_stand.verification.spin_common import native


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


native()
roots = [Path(value) for value in sys.argv[1:]]
require(len(roots) == 2 and roots[0].resolve() != roots[1].resolve(), "two distinct actual runs required")
artifacts = []
for root in roots:
    summary = json.loads((root / "summary.json").read_text())
    require(summary["classification"] == "SYNTHETIC_REGRESSION", "not a paper-case checker")
    require(summary["project_status"] == "PARTIAL_BLOCKED", "inappropriate project claim")
    scientific = json.loads((root / "scientific.json").read_text())
    identity = scientific["identity"]
    freshness = check_freshness(PaperRunEvidenceV2.from_json((root / "evidence.json").read_text()))
    require(freshness.status == "CURRENT", str(freshness))
    artifact = load_spin_artifact(root / "sealed", expected_problem_digest=identity["problem_digest"])
    require(artifact.scientific_digest == identity["artifact_scientific_digest"], "seal identity mismatch")
    facts = artifact.facts
    # Literal independent box-minus-72-cuboid integrals from original fixture.
    world = [503/64, 0, 0, 2003/256, 32687/12288, 32687/12288, 128363/12288, 0, 0, 0]
    np.testing.assert_allclose(facts["physical"]["moments"], world, rtol=0, atol=1e-12)
    np.testing.assert_allclose(facts["physical"]["central_inertia"],
                              np.diag([16432021/3090432]*2+[32687/6144]), rtol=0, atol=1e-12)
    require(facts["bounded_void_components"] == facts["newly_carved_components"] == 1, "actual cavity count")
    require(facts["clearance"]["squared_distance_exact"] == "25/64", "actual wall clearance")
    for kind, array in zip(("vertices", "faces"), artifact.outer.arrays()):
        original = np.load(root / f"outer_{kind}.npy", allow_pickle=False)
        require(original.dtype == array.dtype and original.shape == array.shape and
                original.tobytes() == array.tobytes(), "outer bytes changed")
    run = identity["native_run"]
    require(run["geometry"]["union"]["input_cuboids"] == 72, "native carved cuboid count")
    require(sum(len(p["trace"]) for p in run["solve"]["passes"]) == 39, "complete ordered trace")
    artifacts.append(artifact)
    print(json.dumps({"directory": str(root), "freshness": freshness.status,
                      "artifact_digest": artifact.scientific_digest,
                      "scope": "SYNTHETIC_REGRESSION_GEOMETRY_ONLY"}, sort_keys=True))
for name in ("inputs.json", "scientific.json", "native-run.json"):
    first, second = [(root / name).read_bytes() for root in roots]
    require(first == second, f"scientific replay mismatch: {name}")
    print(name, hashlib.sha256(first).hexdigest(), "BYTE_IDENTICAL")
require(artifacts[0].scientific_digest == artifacts[1].scientific_digest, "different sealed materials")
print("PASS: fresh repeated native runs, independent reload/integration, literal moments/inertia, outer bytes and actual cavity")
