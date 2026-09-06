from __future__ import annotations

import copy
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
STRICT = ROOT / "experiments/spin_it_reproduction/strict"
AUDIT = STRICT / "audit.py"


def load_audit():
    if not AUDIT.is_file():
        raise AssertionError("strict audit CLI has not been implemented")
    spec = importlib.util.spec_from_file_location("spin_it_strict_audit", AUDIT)
    if spec is None or spec.loader is None:
        raise AssertionError("cannot load strict audit CLI")
    module = importlib.util.module_from_spec(spec)
    sys.path.insert(0, str(STRICT))
    try:
        spec.loader.exec_module(module)
    finally:
        sys.path.pop(0)
    return module


def exact_asset(root: Path) -> dict[str, object]:
    asset = root / "author-original.obj"
    asset.write_bytes(b"strict-fixture\n")
    return {
        "asset_id": "fixture-author-object",
        "asset_class": "AUTHOR_ORIGINAL",
        "path": asset.name,
        "sha256": hashlib.sha256(asset.read_bytes()).hexdigest(),
        "source_reference": "author delivery record fixture-001",
        "required_critical_parameters": ["density"],
        "parameter_evidence": [
            {
                "parameter": "density",
                "status": "AUTHOR_PROVIDED",
                "evidence_description": "Author supplied density with fixture-001.",
                "value": 1.25,
            }
        ],
    }


def blocked_payload(audit, root: Path) -> dict[str, object]:
    catalogue = audit.load_catalogue()
    requirements = []
    for row in catalogue["requirements"]:
        status = (
            "NOT_IN_SCOPE"
            if row["requirement_kind"] == "PHYSICAL_MEASUREMENT"
            else "NOT_IMPLEMENTED"
        )
        requirements.append(
            {
                "requirement_id": row["requirement_id"],
                "status": status,
                "evidence_description": "No independent checker result exists in Task 1.",
                "details": {},
            }
        )
    return {
        "schema": "spin-it-strict-evidence-v1",
        "schema_version": 1,
        "requirements": requirements,
        "assets": [exact_asset(root)],
    }


class SpinItStrictCatalogueTests(unittest.TestCase):
    def test_fixed_catalogue_has_individual_required_claims(self) -> None:
        audit = load_audit()
        catalogue = audit.load_catalogue()
        ids = {row["requirement_id"] for row in catalogue["requirements"]}
        self.assertTrue(
            {
                "spin2014.eq.1",
                "spin2014.eq.8",
                "spin2014.case.teapot.inertia_only",
                "spin2014.case.teapot.lowering_only",
                "spin2014.case.teapot.mass_and_lowering",
                "faster.eq14.newton",
                "faster.table1.floating_island_cat",
                "faster.figure.13",
                "shared.real_desktop_visualization",
                "physical.2014.measurements",
                "physical.faster.measurements",
            }.issubset(ids)
        )
        self.assertEqual(len(ids), len(catalogue["requirements"]))
        for row in catalogue["requirements"]:
            self.assertEqual(
                set(row),
                {
                    "requirement_id",
                    "paper_locator",
                    "requirement_kind",
                    "expected_outcome",
                    "dependencies",
                    "required_evidence",
                },
            )
            self.assertTrue(row["required_evidence"].strip())

    def test_catalogue_rejects_missing_duplicate_and_unknown_ids(self) -> None:
        audit = load_audit()
        rows = audit.load_catalogue()["requirements"]
        mutations = {
            "missing": rows[:-1],
            "duplicate": rows + [copy.deepcopy(rows[0])],
            "unknown": rows[:-1]
            + [{**copy.deepcopy(rows[-1]), "requirement_id": "caller.invented"}],
        }
        for label, mutated in mutations.items():
            with self.subTest(label=label):
                with self.assertRaisesRegex(ValueError, label):
                    audit.validate_catalogue(mutated)

    def test_physical_measurements_cannot_be_dropped(self) -> None:
        audit = load_audit()
        rows = [
            row
            for row in audit.load_catalogue()["requirements"]
            if row["requirement_id"] != "physical.2014.measurements"
        ]
        with self.assertRaisesRegex(ValueError, "physical.2014.measurements"):
            audit.validate_catalogue(rows)

    def test_2014_catalogue_uses_paper_sections_and_named_table1_models(self) -> None:
        audit = load_audit()
        rows = {
            row["requirement_id"]: row
            for row in audit.load_catalogue()["requirements"]
        }
        expected_locators = {
            "spin2014.fill.continuous": "Spin-It (2014), Sections 5.2-5.3",
            "spin2014.fill.binary": "Spin-It (2014), Section 5.3",
            "spin2014.octree.boundary_subtrees": "Spin-It (2014), Section 5.4",
            "spin2014.octree.split_merge": "Spin-It (2014), Section 5.3",
            "spin2014.regularizer.laplacian": "Spin-It (2014), Section 5.3",
            "spin2014.solver.slqp": "Spin-It (2014), Section 5.4",
            "spin2014.solver.active_set": "Spin-It (2014), Section 5.4",
            "spin2014.solver.lbfgs": "Spin-It (2014), Section 5.4",
            "spin2014.cage.complete_octree": "Spin-It (2014), Section 6 and Figure 4",
            "spin2014.cage.bbw_interpolation": "Spin-It (2014), Section 6 and Figure 4",
            "spin2014.deformation.joint": "Spin-It (2014), Section 6.1",
            "spin2014.material.dual_density": "Spin-It (2014), Section 7",
            "spin2014.case.teapot.inertia_only": "Spin-It (2014), Section 8 rotational-stability comparison: inertia-only strategy",
            "spin2014.case.teapot.lowering_only": "Spin-It (2014), Section 8 rotational-stability comparison: lowering-only strategy",
            "spin2014.case.teapot.mass_and_lowering": "Spin-It (2014), Section 8 rotational-stability comparison: mass-and-lowering strategy",
            "spin2014.table1.family.c1_cut1": "Spin-It (2014), Table 1: Cuboid, Ellipsoid, Heart, and T-Rex (1 void, 1 cut)",
            "spin2014.table1.family.c1_cut3": "Spin-It (2014), Table 1: Armadillo (hand) (1 void, 3 cuts)",
            "spin2014.table1.family.c1_cut4": "Spin-It (2014), Table 1: Woven Ring (1 void, 4 cuts)",
            "spin2014.table1.family.c1_cut5": "Spin-It (2014), Table 1: Dancing Couple (1 void, 5 cuts)",
            "spin2014.table1.family.c2_cut1": "Spin-It (2014), Table 1: Teapot (2 voids, 1 cut)",
            "spin2014.table1.family.c3_cut2": "Spin-It (2014), Table 1: Elephant (3 voids, 2 cuts)",
            "spin2014.table1.family.c5_cut2": "Spin-It (2014), Table 1: Armadillo (shell) (5 voids, 2 cuts)",
        }
        self.assertEqual(
            {requirement_id: rows[requirement_id]["paper_locator"] for requirement_id in expected_locators},
            expected_locators,
        )

    def test_degenerate_and_hollowing_cases_do_not_require_generic_success(self) -> None:
        audit = load_audit()
        rows = {
            row["requirement_id"]: row
            for row in audit.load_catalogue()["requirements"]
        }
        self.assertEqual(
            rows["spin2014.case.dancing_couple.hollowing"]["expected_outcome"],
            "MATCH_HOLLOWING_INSUFFICIENCY",
        )
        self.assertEqual(
            rows["spin2014.case.dancing_couple.joint_deformation"]["expected_outcome"],
            "MATCH_JOINT_DEFORMATION_IMPROVEMENT",
        )
        self.assertEqual(
            rows["faster.section6_6.example.1"]["expected_outcome"],
            "NONUNIQUE_ZERO_SENSITIVITY",
        )
        self.assertEqual(
            rows["faster.section6_6.example.2"]["expected_outcome"],
            "FULL_DOMAIN_ZERO_SENSITIVITY",
        )


class SpinItStrictAssetPreflightTests(unittest.TestCase):
    def test_exact_specification_asset_can_be_ready(self) -> None:
        audit = load_audit()
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            asset = exact_asset(root)
            asset["asset_class"] = "PUBLISHED_EXACT_SPECIFICATION"
            asset["parameter_evidence"][0]["status"] = "EXACTLY_DERIVED"
            result = audit.preflight_assets([asset], root)
        self.assertEqual(result["status"], "READY")
        self.assertTrue(result["ready"])
        self.assertIs(result.get("source_authenticity_verified"), False)
        self.assertIs(result.get("strict_requirement_verification"), False)

    def test_asset_preflight_rejects_untrusted_or_malformed_inputs(self) -> None:
        audit = load_audit()
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            baseline = exact_asset(root)
            mutations = {}

            bad_digest = copy.deepcopy(baseline)
            bad_digest["sha256"] = "0" * 64
            mutations["digest"] = bad_digest

            source_only = copy.deepcopy(baseline)
            source_only["asset_class"] = "SOURCE_ONLY"
            mutations["asset class"] = source_only

            reconstructed = copy.deepcopy(baseline)
            reconstructed["asset_class"] = "RECONSTRUCTED_ASSET"
            mutations["reconstructed"] = reconstructed

            unknown_parameter = copy.deepcopy(baseline)
            unknown_parameter["parameter_evidence"][0]["status"] = "UNKNOWN_CRITICAL"
            mutations["UNKNOWN_CRITICAL"] = unknown_parameter

            undeclared_parameter = copy.deepcopy(baseline)
            undeclared_parameter["parameter_evidence"].append(
                {
                    "parameter": "invented",
                    "status": "PUBLISHED",
                    "evidence_description": "Invented extra parameter.",
                    "value": 2.0,
                }
            )
            mutations["unknown parameter"] = undeclared_parameter

            nan_parameter = copy.deepcopy(baseline)
            nan_parameter["parameter_evidence"][0]["value"] = math.nan
            mutations["nonfinite"] = nan_parameter

            for message, asset in mutations.items():
                with self.subTest(message=message):
                    with self.assertRaisesRegex(ValueError, message):
                        audit.preflight_assets([asset], root)


class SpinItStrictAuditTests(unittest.TestCase):
    def test_ready_assets_and_invented_passes_still_cannot_publish(self) -> None:
        audit = load_audit()
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            payload = blocked_payload(audit, root)
            result = audit.audit(payload, asset_root=root)
            self.assertEqual(result["asset_preflight"]["status"], "READY")
            self.assertEqual(result["status"], "PARTIAL_BLOCKED")
            self.assertFalse(result["public_complete_reproduction"])
            self.assertFalse(result["formula_coverage"])
            self.assertFalse(result["experiment_coverage"])
            self.assertIn("independent requirement checkers", " ".join(result["blockers"]))

            invented = copy.deepcopy(payload)
            invented["requirements"][0]["passed"] = True
            with self.assertRaisesRegex(ValueError, "passed"):
                audit.audit(invented, asset_root=root)

            invented = copy.deepcopy(payload)
            invented["requirements"][0]["status"] = "VERIFIED"
            with self.assertRaisesRegex(ValueError, "VERIFIED"):
                audit.audit(invented, asset_root=root)

    def test_requirement_evidence_rejects_bad_sets_empty_text_and_nan(self) -> None:
        audit = load_audit()
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            baseline = blocked_payload(audit, root)
            mutations = {}

            missing = copy.deepcopy(baseline)
            missing["requirements"].pop()
            mutations["missing"] = missing

            duplicate = copy.deepcopy(baseline)
            duplicate["requirements"].append(copy.deepcopy(duplicate["requirements"][0]))
            mutations["duplicate"] = duplicate

            unknown = copy.deepcopy(baseline)
            unknown["requirements"][-1]["requirement_id"] = "caller.invented"
            mutations["unknown"] = unknown

            empty = copy.deepcopy(baseline)
            empty["requirements"][0]["evidence_description"] = "   "
            mutations["empty evidence"] = empty

            nonfinite = copy.deepcopy(baseline)
            nonfinite["requirements"][0]["details"] = {"residual": math.inf}
            mutations["nonfinite"] = nonfinite

            for message, payload in mutations.items():
                with self.subTest(message=message):
                    with self.assertRaisesRegex(ValueError, message):
                        audit.audit(payload, asset_root=root)

    def test_cli_exits_nonzero_and_writes_detailed_blockers(self) -> None:
        audit = load_audit()
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            evidence = root / "evidence.json"
            output = root / "audit.json"
            evidence.write_text(
                json.dumps(blocked_payload(audit, root), allow_nan=False),
                encoding="utf-8",
            )
            completed = subprocess.run(
                [
                    sys.executable,
                    str(AUDIT),
                    "--evidence", str(evidence),
                    "--asset-root", str(root),
                    "--output", str(output),
                ],
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(completed.returncode, 2)
            result = json.loads(output.read_text(encoding="utf-8"))
            self.assertEqual(result["status"], "PARTIAL_BLOCKED")
            self.assertGreater(len(result["blockers"]), 1)
            self.assertIn("PARTIAL_BLOCKED", completed.stdout)
            self.assertIn("independent requirement checkers", completed.stdout)


if __name__ == "__main__":
    unittest.main()
