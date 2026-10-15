from __future__ import annotations

import copy
import hashlib
import importlib.util
import json
import math
import inspect
import os
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


def case_asset_readiness_payload(audit, *, exact_ids=()) -> dict[str, object]:
    exact_ids = set(exact_ids)
    rows = []
    for row in audit.load_catalogue()["requirements"]:
        if row["requirement_kind"] != "ORIGINAL_EXPERIMENT":
            continue
        exact = row["requirement_id"] in exact_ids
        rows.append({
            "requirement_id": row["requirement_id"],
            "asset_status": "PUBLISHED_EXACT_SPECIFICATION" if exact else "MISSING",
            "missing_artifacts": [] if exact else ["authenticated original geometry"],
            "missing_critical_parameters": [] if exact else ["target frame"],
            "source_summary": "Synthetic strict-audit readiness fixture.",
        })
    source = ROOT / "experiments/spin_it_reproduction/sources/strict_asset_audit.md"
    return {
        "schema": "spin-it-case-asset-readiness-v1",
        "schema_version": 1,
        "source_audit_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "cases": rows,
    }


class CaseReadinessSnapshotTests(unittest.TestCase):
    def fixture(self, root):
        audit = load_audit()
        payload = case_asset_readiness_payload(audit)
        source = b"Frozen acquisition evidence; no authenticated original input.\n"
        digest = hashlib.sha256(source).hexdigest()
        relative = f"snapshots/{digest}.md"
        snapshot = root/relative
        snapshot.parent.mkdir()
        snapshot.write_bytes(source)
        payload.update(schema="spin-it-case-asset-readiness-v2", schema_version=2,
                       source_audit_path=relative, source_audit_sha256=digest)
        path = root/"readiness.json"
        path.write_text(json.dumps(payload))
        return audit, path, snapshot, payload

    def test_snapshot_binding_survives_unrelated_report_append_without_case_promotion(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            audit, path, snapshot, payload = self.fixture(root)
            try:
                rows = audit.load_case_asset_readiness(path)
            except ValueError as error:
                self.fail(f"explicit frozen source binding is unsupported: {error}")
            (root/"strict_asset_audit.md").write_text("New report referring to inventory hash.\n")
            record = audit.load_case_asset_readiness_record(path)
            self.assertEqual(record["cases"], rows)
            self.assertTrue(all(row["asset_status"] == "MISSING" for row in rows.values()))
            self.assertEqual(record["source_audit"]["path"], str(snapshot.resolve()))
            self.assertEqual(record["source_audit"]["sha256"], payload["source_audit_sha256"])
            self.assertEqual(record["source_audit"]["scope"], "FROZEN_SOURCE_SNAPSHOT_NOT_ASSET_AUTHENTICITY")
            result = audit.audit(blocked_payload(audit, root), asset_root=root,
                                 case_asset_readiness=path, paper="faster")
            self.assertEqual(result["case_asset_readiness"]["source_audit"], record["source_audit"])
            self.assertEqual(result["status"], "PARTIAL_BLOCKED")

    def test_snapshot_tampering_is_not_accepted(self):
        with tempfile.TemporaryDirectory() as directory:
            audit, path, snapshot, _ = self.fixture(Path(directory))
            snapshot.write_text("Changed acquisition evidence\n")
            with self.assertRaisesRegex(ValueError, "SHA-256 mismatch"):
                audit.load_case_asset_readiness(path)

    def test_snapshot_reference_cannot_escape_or_impersonate_another_version(self):
        with tempfile.TemporaryDirectory() as directory:
            audit, path, snapshot, payload = self.fixture(Path(directory))
            for reference in (str(snapshot.resolve()), "../strict_asset_audit.md", "strict_asset_audit.md"):
                payload["source_audit_path"] = reference
                path.write_text(json.dumps(payload))
                with self.subTest(reference=reference), self.assertRaisesRegex(ValueError, "snapshot path"):
                    audit.load_case_asset_readiness(path)
            payload["source_audit_path"] = f"snapshots/{payload['source_audit_sha256']}.md"
            payload["schema_version"] = 1
            path.write_text(json.dumps(payload))
            with self.assertRaisesRegex(ValueError, "version"):
                audit.load_case_asset_readiness(path)

    def test_snapshot_symlink_escape_and_missing_file_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory, tempfile.TemporaryDirectory() as outside:
            audit, path, snapshot, _ = self.fixture(Path(directory))
            external = Path(outside)/"source.md"
            external.write_bytes(snapshot.read_bytes())
            snapshot.unlink()  # Only this test's private temporary file.
            with self.assertRaises(FileNotFoundError):
                audit.load_case_asset_readiness(path)
            snapshot.symlink_to(external)
            with self.assertRaisesRegex(ValueError, "snapshot path"):
                audit.load_case_asset_readiness(path)


class SpinItStrictCatalogueTests(unittest.TestCase):
    def test_paper_scope_is_fixed_and_includes_shared_requirements(self) -> None:
        audit = load_audit()
        faster = set(audit.paper_scope_requirement_ids("faster"))
        spin = set(audit.paper_scope_requirement_ids("spin2014"))
        all_ids = set(audit.paper_scope_requirement_ids("all"))
        shared = {item for item in all_ids if item.startswith("shared.")}
        self.assertTrue(shared)
        self.assertTrue({item for item in faster if item.startswith("faster.")})
        self.assertTrue({item for item in spin if item.startswith("spin2014.")})
        self.assertTrue(shared <= faster)
        self.assertTrue(shared <= spin)
        self.assertIn("physical.faster.measurements", faster)
        self.assertNotIn("physical.2014.measurements", faster)
        self.assertIn("physical.2014.measurements", spin)
        self.assertEqual(faster | spin, all_ids)
        with self.assertRaisesRegex(ValueError, "paper scope"):
            audit.paper_scope_requirement_ids("invented")

    def test_paper_completion_is_derived_from_audited_requirement_rows(self) -> None:
        audit = load_audit()
        rows = [
            {"requirement_id": requirement_id, "status": "VERIFIED"}
            for requirement_id in audit.paper_scope_requirement_ids("faster")
        ]
        kinds = {row["requirement_id"]: row["requirement_kind"]
                 for row in audit.load_catalogue()["requirements"]}
        readiness = {
            requirement_id: {"asset_status": "PUBLISHED_EXACT_SPECIFICATION"}
            for requirement_id in audit.paper_scope_requirement_ids("faster")
            if kinds[requirement_id] == "ORIGINAL_EXPERIMENT"
        }
        complete = audit.derive_paper_status("faster", rows, readiness)
        self.assertEqual(complete["status"], "VERIFIED_COMPUTATIONAL_REPRODUCTION")
        self.assertTrue(complete["public_complete_reproduction"])

        blocked = audit.derive_paper_status("faster", rows[:-1], readiness)
        self.assertEqual(blocked["status"], "PARTIAL_BLOCKED")
        self.assertFalse(blocked["public_complete_reproduction"])
        self.assertTrue(blocked["blockers"])

    def test_author_original_readiness_can_satisfy_the_input_gate_only_after_checks(self) -> None:
        audit = load_audit()
        rows = [
            {"requirement_id": requirement_id, "status": "VERIFIED"}
            for requirement_id in audit.paper_scope_requirement_ids("faster")
        ]
        kinds = {row["requirement_id"]: row["requirement_kind"]
                 for row in audit.load_catalogue()["requirements"]}
        readiness = {
            requirement_id: {"asset_status": "AUTHOR_ORIGINAL"}
            for requirement_id in audit.paper_scope_requirement_ids("faster")
            if kinds[requirement_id] == "ORIGINAL_EXPERIMENT"
        }
        result = audit.derive_paper_status("faster", rows, readiness)
        self.assertEqual(result["status"], "VERIFIED_COMPUTATIONAL_REPRODUCTION")

        rows[0]["status"] = "BLOCKED"
        rows[0]["evidence_description"] = "Independent check absent."
        blocked = audit.derive_paper_status("faster", rows, readiness)
        self.assertEqual(blocked["status"], "PARTIAL_BLOCKED")

    def test_audit_can_report_faster_without_hiding_spin_requirements(self) -> None:
        audit = load_audit()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            payload = blocked_payload(audit, root)
            faster = audit.audit(payload, asset_root=root, paper="faster")
            all_papers = audit.audit(payload, asset_root=root)
        self.assertEqual(faster["paper_scope"], "faster")
        self.assertEqual(faster["status"], "PARTIAL_BLOCKED")
        self.assertFalse(faster["public_complete_reproduction"])
        self.assertTrue(all(row["requirement_id"] in {
            item["requirement_id"] for item in faster["requirements"]
        } for row in all_papers["requirements"]))
        self.assertFalse(any(item.startswith("spin2014.") for item in faster["blockers"]))
        self.assertEqual(all_papers["paper_scope"], "all")
        self.assertEqual(all_papers["status"], "PARTIAL_BLOCKED")

    def test_cli_accepts_fixed_faster_scope_and_keeps_blocked_exit_nonzero(self) -> None:
        audit = load_audit()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            evidence_path = root / "evidence.json"
            evidence_path.write_text(json.dumps(blocked_payload(audit, root)), encoding="utf-8")
            output_path = root / "audit.json"
            completed = subprocess.run([
                sys.executable, str(AUDIT), "--evidence", str(evidence_path),
                "--asset-root", str(root), "--output", str(output_path),
                "--paper", "faster",
            ], capture_output=True, text=True, check=False)
            self.assertEqual(completed.returncode, 2)
            self.assertTrue(output_path.is_file(), completed.stderr + completed.stdout)
            result = json.loads(output_path.read_text(encoding="utf-8"))
        self.assertEqual(result["paper_scope"], "faster")
        self.assertEqual(result["status"], "PARTIAL_BLOCKED")
        self.assertFalse(any(item.startswith("spin2014.") for item in result["blockers"]))

    def test_missing_fixed_shell_requirement_cannot_be_silently_omitted(self) -> None:
        audit = load_audit()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            payload = blocked_payload(audit, root)
            payload['requirements'] = [row for row in payload['requirements']
                if row['requirement_id'] != 'spin2014.hollowing.fixed_shell']
            with self.assertRaisesRegex(ValueError, 'spin2014.hollowing.fixed_shell'):
                audit.audit(payload, asset_root=root)

    def test_missing_minimum_wall_requirement_cannot_be_silently_omitted(self) -> None:
        """Section 5.1's fabrication thickness rule is not fixed-shell clearance."""
        audit = load_audit()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            payload = blocked_payload(audit, root)
            payload['requirements'] = [
                row for row in payload['requirements']
                if row['requirement_id'] != 'spin2014.hollowing.minimum_wall'
            ]
            with self.assertRaisesRegex(ValueError, 'spin2014.hollowing.minimum_wall'):
                audit.audit(payload, asset_root=root)

    def test_missing_deformation_regularizers_cannot_be_silently_omitted(self) -> None:
        """Section 6.1's two deformation regularizers are not jointness alone."""
        audit = load_audit()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            payload = blocked_payload(audit, root)
            payload['requirements'] = [
                row for row in payload['requirements']
                if row['requirement_id'] != 'spin2014.deformation.regularizers'
            ]
            with self.assertRaisesRegex(
                ValueError, 'spin2014.deformation.regularizers'
            ):
                audit.audit(payload, asset_root=root)

    def test_missing_faster_application_contracts_cannot_be_silently_omitted(self) -> None:
        """Table 1 cases do not replace Section 4's application contracts."""
        audit = load_audit()
        requirement_ids = (
            'faster.application.static_stability',
            'faster.application.buoyant_stability',
            'faster.application.fidget_spinner',
            'faster.application.yoyo',
            'faster.application.spinning_top',
        )
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for requirement_id in requirement_ids:
                with self.subTest(requirement_id=requirement_id):
                    payload = blocked_payload(audit, root)
                    payload['requirements'] = [
                        row for row in payload['requirements']
                        if row['requirement_id'] != requirement_id
                    ]
                    with self.assertRaisesRegex(ValueError, requirement_id):
                        audit.audit(payload, asset_root=root)

    def test_fixed_shell_handwritten_pass_does_not_replace_independent_check(self) -> None:
        audit = load_audit()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            payload = blocked_payload(audit, root)
            rows = [r for r in payload['requirements']
                    if r['requirement_id'] == 'spin2014.hollowing.fixed_shell']
            self.assertEqual(len(rows), 1, 'fixed shell must be an explicit requirement')
            rows[0]['status'] = 'VERIFIED'
            with self.assertRaisesRegex(ValueError, 'caller status cannot replace'):
                audit.audit(payload, asset_root=root)

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

    def test_caller_mutation_cannot_change_canonical_catalogue(self) -> None:
        audit = load_audit()
        exposed = getattr(audit, "REQUIREMENTS", ())
        if exposed:
            first = exposed[0]
            dependent = next(row for row in exposed if row["dependencies"])
            original_outcome = first["expected_outcome"]
            original_evidence = first["required_evidence"]
            try:
                first["expected_outcome"] = "CALLER_ALTERED"
                first["required_evidence"] = "Caller replacement."
                dependent["dependencies"].append("caller.invented")
                fresh = audit.load_catalogue()["requirements"]
                self.assertEqual(fresh[0]["expected_outcome"], original_outcome)
                self.assertEqual(fresh[0]["required_evidence"], original_evidence)
                self.assertNotIn(
                    "caller.invented",
                    next(
                        row["dependencies"]
                        for row in fresh
                        if row["requirement_id"] == dependent["requirement_id"]
                    ),
                )
            finally:
                first["expected_outcome"] = original_outcome
                first["required_evidence"] = original_evidence
                dependent["dependencies"].remove("caller.invented")

        first_copy = audit.load_catalogue()
        first_copy["requirements"][0]["required_evidence"] = "Caller replacement."
        dependent_copy = next(
            row for row in first_copy["requirements"] if row["dependencies"]
        )
        dependent_copy["dependencies"].append("caller.invented")
        fresh_copy = audit.load_catalogue()["requirements"]
        self.assertNotEqual(fresh_copy[0]["required_evidence"], "Caller replacement.")
        self.assertNotIn(
            "caller.invented",
            next(
                row["dependencies"]
                for row in fresh_copy
                if row["requirement_id"] == dependent_copy["requirement_id"]
            ),
        )

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

    def test_faster_theorems_use_authenticated_numbered_locators(self) -> None:
        audit = load_audit()
        rows = {
            row["requirement_id"]: row
            for row in audit.load_catalogue()["requirements"]
        }
        self.assertEqual(
            rows["faster.theorem.quadric_reduction"]["paper_locator"],
            "Spin-It Faster, Theorem 1, Section 3.3",
        )
        self.assertEqual(
            rows["faster.theorem.optimality"]["paper_locator"],
            "Spin-It Faster, Theorem 2, Section 3.4",
        )
        self.assertEqual(
            rows["faster.theorem.degenerate_conditions"]["paper_locator"],
            "Spin-It Faster, Theorem 2, Section 3.4 (p = 0 branch)",
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

    def test_asset_preflight_rejects_empty_or_boolean_parameter_values(self) -> None:
        audit = load_audit()
        invalid_values = {
            "null": None,
            "boolean true": True,
            "boolean false": False,
            "empty text": "",
            "whitespace text": "   ",
            "empty array": [],
            "empty object": {},
            "nested null": {"published": None},
            "nested boolean": ["published", False],
            "nested empty text": {"published": "  "},
            "nested empty array": {"published": []},
            "nested empty object": ["published", {}],
        }
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            baseline = exact_asset(root)
            for label, invalid_value in invalid_values.items():
                with self.subTest(label=label):
                    asset = copy.deepcopy(baseline)
                    asset["parameter_evidence"][0]["value"] = invalid_value
                    with self.assertRaisesRegex(ValueError, "parameter value"):
                        audit.preflight_assets([asset], root)

            zero = copy.deepcopy(baseline)
            zero["parameter_evidence"][0]["value"] = 0
            self.assertTrue(audit.preflight_assets([zero], root)["ready"])

            structured = copy.deepcopy(baseline)
            structured["parameter_evidence"][0]["value"] = {
                "components": [0, "published zero"]
            }
            self.assertTrue(audit.preflight_assets([structured], root)["ready"])

            missing = copy.deepcopy(baseline)
            del missing["parameter_evidence"][0]["value"]
            with self.assertRaisesRegex(ValueError, "missing fields.*value"):
                audit.preflight_assets([missing], root)


class SpinItStrictAuditTests(unittest.TestCase):
    def test_author_original_case_readiness_is_eligible_but_cannot_verify_itself(self) -> None:
        audit = load_audit()
        requirement_id = "faster.table1.heart_spinning_top"
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            readiness = case_asset_readiness_payload(audit, exact_ids={requirement_id})
            case = next(row for row in readiness["cases"]
                        if row["requirement_id"] == requirement_id)
            case["asset_status"] = "AUTHOR_ORIGINAL"
            readiness_path = root / "case-readiness.json"
            readiness_path.write_text(json.dumps(readiness), encoding="utf-8")
            loaded = audit.load_case_asset_readiness(readiness_path)
            result = audit.audit(blocked_payload(audit, root), asset_root=root,
                                 case_asset_readiness=readiness_path, paper="faster")
        self.assertEqual(loaded[requirement_id]["asset_status"], "AUTHOR_ORIGINAL")
        self.assertEqual(result["requirement_status_index"][requirement_id]["original_experiment"],
                         "ELIGIBLE_INPUT_AWAITING_VERIFICATION")
        self.assertEqual(result["status"], "PARTIAL_BLOCKED")

    def test_author_original_case_cannot_claim_readiness_with_missing_parameters(self) -> None:
        audit = load_audit()
        requirement_id = "faster.table1.heart_spinning_top"
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            readiness = case_asset_readiness_payload(audit, exact_ids={requirement_id})
            case = next(row for row in readiness["cases"]
                        if row["requirement_id"] == requirement_id)
            case["asset_status"] = "AUTHOR_ORIGINAL"
            case["missing_critical_parameters"] = ["target axis"]
            readiness_path = root / "case-readiness.json"
            readiness_path.write_text(json.dumps(readiness), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "cannot list missing evidence"):
                audit.load_case_asset_readiness(readiness_path)

    def test_case_asset_readiness_is_complete_and_controls_case_eligibility(self) -> None:
        audit = load_audit()
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            readiness_path = root / "case-readiness.json"
            readiness_path.write_text(json.dumps(case_asset_readiness_payload(
                audit, exact_ids={"faster.section6_6.example.1", "faster.section6_6.example.2"}
            )), encoding="utf-8")

            readiness = audit.load_case_asset_readiness(readiness_path)
            readiness_digest = hashlib.sha256(readiness_path.read_bytes()).hexdigest()
            self.assertEqual(
                readiness["faster.section6_6.example.1"]["asset_status"],
                "PUBLISHED_EXACT_SPECIFICATION",
            )
            result = audit.audit(blocked_payload(audit, root), asset_root=root,
                                 case_asset_readiness=readiness_path)

        index = result["requirement_status_index"]
        self.assertEqual(
            result["case_asset_readiness"]["source"]["sha256"],
            readiness_digest,
        )
        self.assertEqual(
            index["faster.section6_6.example.1"]["original_experiment"],
            "ELIGIBLE_INPUT_AWAITING_VERIFICATION",
        )
        self.assertEqual(
            index["spin2014.case.teapot.inertia_only"]["original_experiment"],
            "NOT_ELIGIBLE_MISSING_ORIGINAL_EVIDENCE",
        )

    def test_case_asset_readiness_rejects_missing_original_case_row(self) -> None:
        audit = load_audit()
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            readiness = case_asset_readiness_payload(audit)
            readiness["cases"] = [row for row in readiness["cases"]
                                  if row["requirement_id"] != "spin2014.case.ellipsoid"]
            readiness_path = root / "case-readiness.json"
            readiness_path.write_text(json.dumps(readiness), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "spin2014.case.ellipsoid"):
                audit.load_case_asset_readiness(readiness_path)

    def test_audit_emits_four_column_requirement_status_index(self) -> None:
        audit = load_audit()
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            result = audit.audit(blocked_payload(audit, root), asset_root=root)

        index = result["requirement_status_index"]
        self.assertEqual(set(index), {
            row["requirement_id"] for row in audit.load_catalogue()["requirements"]
        })
        row = index["spin2014.case.teapot.inertia_only"]
        self.assertEqual(row["implementation"], "UNASSESSED")
        self.assertEqual(row["independent_verification"], "NOT_VERIFIED")
        self.assertEqual(row["original_experiment"], "NOT_ELIGIBLE_MISSING_ORIGINAL_EVIDENCE")
        self.assertIn("asset_evidence", row["blocker_dependencies"])

        physical = index["physical.2014.measurements"]
        self.assertEqual(physical["implementation"], "NOT_APPLICABLE")
        self.assertEqual(physical["independent_verification"], "NOT_IN_SCOPE")
        self.assertEqual(physical["original_experiment"], "NOT_IN_SCOPE")

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


class IndependentCheckerIntegrationTests(unittest.TestCase):
    def test_faster_extrusion_checker_audits_the_reported_arc_mesh_convergence(self):
        audit = load_audit()
        extrusion = {
            "schema": "faster-planar-component-result-v1",
            "classification": "SYNTHETIC_REGRESSION",
            "component": "extrusion",
            "planar_boundary_loop_count": 2,
            "surface_shell_count": 2,
            "extrusion_height": 1.5,
            "density": 2.4,
            "planar_mass": 73.8,
            "extruded_mass": 73.80000000000001,
            "maximum_scaled_absolute_moment_error": 5.736544440489255e-16,
            "exact_arc_extrusion_mesh_convergence_error": 3.9218280490659897e-7,
            "maximum_scaled_objective_constraint_error": 1.7351470958732605e-16,
            "planar_objective": -556.125,
            "spatial_objective": -556.125,
            "planar_constraints": [40.949999999999996, 112.725],
            "spatial_constraints": [40.95, 112.72499999999998],
            "scaled_absolute_error_tolerance": 2e-12,
        }
        self.assertEqual(
            audit._validate_faster_planar_component_result("extrusion", extrusion),
            extrusion,
        )
        negative_error = dict(extrusion, exact_arc_extrusion_mesh_convergence_error=-1e-12)
        with self.assertRaisesRegex(ValueError, "arc-extrusion mesh convergence"):
            audit._validate_faster_planar_component_result("extrusion", negative_error)
        excessive_error = dict(extrusion, exact_arc_extrusion_mesh_convergence_error=1.1e-6)
        with self.assertRaisesRegex(ValueError, "arc-extrusion mesh convergence"):
            audit._validate_faster_planar_component_result("extrusion", excessive_error)

    def test_faster_planar_component_contract_rejects_out_of_contract_numbers(self):
        audit = load_audit()
        dogleg = {
            "schema": "faster-planar-component-result-v1",
            "classification": "SYNTHETIC_REGRESSION",
            "component": "dogleg",
            "status": "CONVERGED",
            "solved_center": [1.0000000001, 0.0],
            "known_root_center": [1.0, 0.0],
            "center_error": 1e-10,
            "normalized_residual": 8e-12,
            "residual_tolerance": 1e-10,
            "iterations": 5,
            "function_evaluations": 26,
            "function_evaluation_budget": 60,
            "exact_root_residual": 4e-18,
            "zero_radius_initial_converged": True,
            "budget_exhaustion_distinguished": True,
            "topology_transition_detected": True,
        }
        self.assertEqual(
            audit._validate_faster_planar_component_result("dogleg", dogleg),
            dogleg,
        )
        bad_residual = dict(dogleg, normalized_residual=2e-10)
        with self.assertRaisesRegex(ValueError, "dogleg residual"):
            audit._validate_faster_planar_component_result("dogleg", bad_residual)
        bad_evaluations = dict(dogleg, function_evaluations=61)
        with self.assertRaisesRegex(ValueError, "function-evaluation budget"):
            audit._validate_faster_planar_component_result("dogleg", bad_evaluations)

    def test_source_synthetic_repeats_verify_only_scoped_fill_slqp_and_fixed_shell_methods(self):
        """A real native synthetic trace may establish method behavior only.

        It must not grant an original model, provenance, or adaptive
        split/merge requirement merely because the small fixture seals.  The
        fixed-shell *method* is different: two independently reloadable seals
        whose exact partitions and outer-wall clearance are recomputed may
        establish that narrowly scoped synthetic capability.
        """

        audit = load_audit()
        import pybind11
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            payload = blocked_payload(audit, root)
            paths = []
            for number in range(2):
                output = root / str(number)
                completed = subprocess.run([
                    sys.executable, str(STRICT / "run_source_spin2014_material.py"),
                    "--output-directory", str(output),
                    "--build-directory", str(ROOT / "build/strict-algorithms-exact-current-20260922"),
                    "--dependency-root", "pybind11=" + pybind11.get_include(),
                ], capture_output=True, text=True,
                    env={**os.environ,
                         "PYTHONPATH": os.pathsep.join((str(ROOT / "python"),
                                                        str(ROOT / "build/strict-algorithms-exact-current-20260922/python"))),
                         "LD_PRELOAD": str(ROOT / "build/strict-algorithms-exact-current-20260922/compiler-runtime/libstdc++.so.6"),
                         "PYTHONDONTWRITEBYTECODE": "1"})
                self.assertEqual(completed.returncode, 0, completed.stderr)
                paths.append(str(output / "evidence.json"))
            transition_runner = STRICT / "run_source_spin2014_adaptive_transitions.py"
            transition_paths = []
            for number in range(2):
                output = root / ("transition-" + str(number))
                completed = subprocess.run([
                    sys.executable, str(transition_runner),
                    "--output-directory", str(output),
                    "--build-directory", str(ROOT / "build/strict-algorithms-exact-current-20260922"),
                    "--dependency-root", "pybind11=" + pybind11.get_include(),
                ], capture_output=True, text=True,
                    env={**os.environ,
                         "PYTHONPATH": os.pathsep.join((str(ROOT / "python"),
                                                        str(ROOT / "build/strict-algorithms-exact-current-20260922/python"))),
                         "LD_PRELOAD": str(ROOT / "build/strict-algorithms-exact-current-20260922/compiler-runtime/libstdc++.so.6"),
                         "PYTHONDONTWRITEBYTECODE": "1"})
                self.assertEqual(completed.returncode, 0, completed.stderr)
                transition_paths.append(str(output / "evidence.json"))
            cage_paths = []
            cage_runner = STRICT / "run_source_spin2014_cage.py"
            for number in range(2):
                output = root / ("cage-" + str(number))
                completed = subprocess.run([
                    sys.executable, str(cage_runner),
                    "--output-directory", str(output),
                    "--build-directory", str(ROOT / "build/strict-algorithms-exact-current-20260922"),
                    "--dependency-root", "pybind11=" + pybind11.get_include(),
                ], capture_output=True, text=True,
                    env={**os.environ,
                         "PYTHONPATH": os.pathsep.join((str(ROOT / "python"),
                                                        str(ROOT / "build/strict-algorithms-exact-current-20260922/python"))),
                         "LD_PRELOAD": str(ROOT / "build/strict-algorithms-exact-current-20260922/compiler-runtime/libstdc++.so.6"),
                         "PYTHONDONTWRITEBYTECODE": "1"})
                self.assertEqual(completed.returncode, 0, completed.stderr)
                cage_paths.append(str(output / "evidence.json"))
            deformation_paths = []
            for number in range(2):
                output = root / ("deformation-" + str(number))
                completed = subprocess.run([sys.executable, str(STRICT / "run_source_spin2014_deformation.py"),
                    "--output-directory", str(output), "--build-directory", str(ROOT / "build/strict-algorithms-exact-current-20260922"),
                    "--dependency-root", "pybind11=" + pybind11.get_include()], capture_output=True, text=True,
                    env={**os.environ, "PYTHONPATH": os.pathsep.join((str(ROOT / "python"),str(ROOT / "build/strict-algorithms-exact-current-20260922/python"))),"LD_PRELOAD":str(ROOT / "build/strict-algorithms-exact-current-20260922/compiler-runtime/libstdc++.so.6"),"PYTHONDONTWRITEBYTECODE":"1"})
                self.assertEqual(completed.returncode, 0, completed.stderr)
                deformation_paths.append(str(output / "evidence.json"))
            dual_material_paths = []
            dual_runner = STRICT / "run_spin2014_dual_material.py"
            for number in range(2):
                output = root / ("dual-material-" + str(number))
                completed = subprocess.run([
                    sys.executable, str(dual_runner),
                    "--output-directory", str(output),
                    "--build-directory", str(ROOT / "build/strict-algorithms-exact-current-20260922"),
                    "--dependency-root", "pybind11=" + pybind11.get_include(),
                ], capture_output=True, text=True,
                    env={**os.environ,
                         "PYTHONPATH": os.pathsep.join((str(ROOT / "python"),
                                                        str(ROOT / "build/strict-algorithms-exact-current-20260922/python"))),
                         "LD_PRELOAD": str(ROOT / "build/strict-algorithms-exact-current-20260922/compiler-runtime/libstdc++.so.6"),
                         "PYTHONDONTWRITEBYTECODE": "1"})
                self.assertEqual(completed.returncode, 0, completed.stderr)
                self.assertIn(
                    "independent_label_moments=validated",
                    (output / "dual-material.stdout").read_text(encoding="utf-8"),
                    "dual-material evidence must include the independent per-label cuboid integral",
                )
                dual_material_paths.append(str(output / "evidence.json"))
            for evidence_path in (
                paths + transition_paths + cage_paths + deformation_paths
                + dual_material_paths
            ):
                with Path(evidence_path).open("r", encoding="utf-8") as stream:
                    request = json.load(stream)["begin"]["request"]
                self.assertEqual(
                    request.get("source_policy"),
                    "make-it-stand-explicit-effective-source-v2",
                    "method evidence must bind only its declared effective inputs",
                )
            result = audit.audit(payload, asset_root=root,
                                 independent_evidence={
                                     "spin2014.source_method": paths,
                                     "spin2014.source_adaptive_transitions": transition_paths,
                                     "spin2014.source_cage": cage_paths,
                                     "spin2014.source_deformation": deformation_paths,
                                     "spin2014.dual_material": dual_material_paths,
                                 })
            verified = {row["requirement_id"] for row in result["requirements"]
                        if row["status"] == "VERIFIED"}
            self.assertEqual(verified, {
                "spin2014.fill.continuous", "spin2014.fill.binary",
                "spin2014.solver.slqp", "spin2014.solver.active_set",
                "spin2014.solver.lbfgs",
                "spin2014.hollowing.fixed_shell",
                "spin2014.octree.boundary_subtrees", "spin2014.octree.split_merge",
                "spin2014.regularizer.laplacian",
                "spin2014.cage.complete_octree", "spin2014.cage.bbw_interpolation",
                "spin2014.deformation.regularizers",
                "spin2014.deformation.joint",
                "spin2014.material.dual_density",
            })
            self.assertEqual(result["requirement_summary"]["verified"], 14)
            self.assertNotIn("shared.sealed_geometry", verified)
            self.assertNotIn("shared.provenance", verified)
            with (root / "0/native-run.json").open("a", encoding="utf-8") as stream:
                stream.write(" ")
            with self.assertRaisesRegex(ValueError, "STALE_RESULT"):
                audit.audit(payload, asset_root=root,
                            independent_evidence={"spin2014.source_method": paths})

    def test_published_spin2014_givens_repeats_verify_only_scoped_equations(self):
        audit = load_audit()
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            payload = blocked_payload(audit, root)
            paths = []
            for number in range(2):
                output = root / str(number)
                completed = subprocess.run([
                    sys.executable, str(STRICT / "run_spin2014_givens_equations.py"),
                    "--paper-text", str(ROOT / "build/strict-paper-sources/spin2014.txt"),
                    "--native-module-directory", str(ROOT / "build/strict-current-tree-20260921/python"),
                    "--output-directory", str(output)], capture_output=True, text=True,
                    env={**os.environ, "PYTHONPATH": str(ROOT / "python"),
                         "LD_PRELOAD": "/lib/x86_64-linux-gnu/libstdc++.so.6",
                         "PYTHONDONTWRITEBYTECODE": "1"})
                self.assertEqual(completed.returncode, 0, completed.stderr)
                paths.append(str(output / "evidence.json"))
            result = audit.audit(payload, asset_root=root,
                                 independent_evidence={"spin2014.givens_equations": paths})
            scientific = json.loads(
                (root / "0" / "result.json").read_text(encoding="utf-8")
            )["scientific_data"]
            objective_checks = {
                check["objective"] for check in scientific["checks"]
            }
            self.assertEqual(objective_checks, {
                "YOYO", "TOP_INERTIA_ONLY", "TOP_LOWERING_ONLY",
                "TOP_MASS_AND_LOWERING", "STATIC_LOWERING",
            })
            verified = {row["requirement_id"] for row in result["requirements"]
                        if row["status"] == "VERIFIED"}
            self.assertEqual(verified, {
                "spin2014.eq.1", "spin2014.eq.2", "spin2014.eq.3", "spin2014.eq.5",
                "spin2014.eq.6", "spin2014.eq.7", "spin2014.eq.8",
                "spin2014.objective.static", "spin2014.objective.variants",
            })
            self.assertEqual(result["requirement_summary"]["verified"], 9)
            self.assertNotIn("spin2014.eq.4", verified)
            with (root / "0/result.json").open("a", encoding="utf-8") as stream:
                stream.write(" ")
            with self.assertRaisesRegex(ValueError, "STALE_RESULT"):
                audit.audit(payload, asset_root=root,
                            independent_evidence={"spin2014.givens_equations": paths})

    def test_spin2014_givens_evidence_ignores_unrelated_document_observation(self):
        """A documentation-only workspace change is not a scientific input."""
        runner_path = STRICT / "run_spin2014_givens_equations.py"
        spec = importlib.util.spec_from_file_location("spin2014_givens_runner", runner_path)
        self.assertIsNotNone(spec)
        self.assertIsNotNone(spec.loader)
        runner = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(runner)
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / "givens"
            completed = subprocess.run([
                sys.executable, str(runner_path),
                "--paper-text", str(ROOT / "build/strict-paper-sources/spin2014.txt"),
                "--native-module-directory", str(ROOT / "build/strict-current-tree-20260921/python"),
                "--output-directory", str(output)], capture_output=True, text=True,
                env={**os.environ, "PYTHONPATH": str(ROOT / "python"),
                     "LD_PRELOAD": "/lib/x86_64-linux-gnu/libstdc++.so.6",
                     "PYTHONDONTWRITEBYTECODE": "1"})
            self.assertEqual(completed.returncode, 0, completed.stderr)
            unrelated = ROOT / "docs" / ".strict-evidence-unrelated-observation.md"
            self.assertFalse(unrelated.exists())
            unrelated.write_text("temporary non-scientific observation\n", encoding="utf-8")
            try:
                assessment = runner.check_evidence_freshness(output / "evidence.json")
            finally:
                unrelated.unlink(missing_ok=True)
            self.assertEqual(assessment["status"], "CURRENT", assessment["mismatches"])

    def test_published_spin2014_supplement_repeats_verify_eq4_and_two_registered_algorithms(self):
        audit = load_audit()
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            payload = blocked_payload(audit, root)
            paths = []
            for number in range(2):
                output = root / str(number)
                completed = subprocess.run([
                    sys.executable, str(STRICT / "run_spin2014_supplement_algorithms.py"),
                    "--paper", str(ROOT / "build/strict-paper-sources/SpinItSupMat.pdf"),
                    "--paper-text", str(ROOT / "build/strict-paper-sources/spin2014.txt"),
                    "--native-module-directory", str(ROOT / "build/strict-current-tree-20260921/python"),
                    "--output-directory", str(output)], capture_output=True, text=True,
                    env={**os.environ, "PYTHONPATH": str(ROOT / "python"),
                         "LD_PRELOAD": "/lib/x86_64-linux-gnu/libstdc++.so.6",
                         "PYTHONDONTWRITEBYTECODE": "1"})
                self.assertEqual(completed.returncode, 0, completed.stderr)
                paths.append(str(output / "evidence.json"))
            result = audit.audit(payload, asset_root=root,
                                 independent_evidence={"spin2014.supplement_algorithms": paths})
            verified = {row["requirement_id"] for row in result["requirements"]
                        if row["status"] == "VERIFIED"}
            self.assertEqual(verified, {
                "spin2014.eq.4", "spin2014.supplement.algorithm.1",
                "spin2014.supplement.algorithm.2"
            })
            self.assertEqual(result["requirement_summary"]["verified"], 3)
            self.assertIn(
                "Only registered independent checkers can verify their explicitly scoped requirements.",
                result["limitations"],
            )
            with (root / "0/result.json").open("a", encoding="utf-8") as stream:
                stream.write(" ")
            with self.assertRaisesRegex(ValueError, "STALE_RESULT"):
                audit.audit(payload, asset_root=root,
                            independent_evidence={"spin2014.supplement_algorithms": paths})

    def test_real_exact_repeats_only_verify_two_registered_requirements(self):
        audit = load_audit()
        self.assertIn("independent_evidence", inspect.signature(audit.audit).parameters)
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            payload = blocked_payload(audit, root)
            paths = []
            for number in range(2):
                output = root / str(number)
                completed = subprocess.run([
                    sys.executable, str(STRICT / "run_faster_section_6_6.py"),
                    "--spec", str(ROOT / "experiments/spin_it_reproduction/sources/faster_section_6_6_spec.json"),
                    "--paper", str(ROOT / "build/strict-paper-sources/sif-final.pdf"),
                    "--output-directory", str(output)], capture_output=True, text=True,
                    env={**os.environ, "PYTHONPATH": str(ROOT / "python"), "PYTHONDONTWRITEBYTECODE": "1"})
                self.assertEqual(completed.returncode, 0, completed.stderr)
                paths.append(str(output / "evidence.json"))
            result = audit.audit(payload, asset_root=root,
                                 independent_evidence={"faster.section6_6": paths})
            verified = {r["requirement_id"] for r in result["requirements"] if r["status"] == "VERIFIED"}
            self.assertEqual(verified, {
                "faster.section6_6.example.1", "faster.section6_6.example.2",
                "faster.theorem.degenerate_conditions",
            })
            self.assertEqual(result["requirement_summary"]["verified"], 3)
            self.assertEqual(result["status"], "PARTIAL_BLOCKED")
            self.assertFalse(result["public_complete_reproduction"])
            (root / "payload.json").write_text(json.dumps(payload))
            (root / "checks.json").write_text(json.dumps({"faster.section6_6": paths}))
            cli = subprocess.run([sys.executable, str(AUDIT),
                "--evidence", str(root / "payload.json"), "--asset-root", str(root),
                "--independent-evidence", str(root / "checks.json"),
                "--output", str(root / "audit.json")], capture_output=True, text=True)
            self.assertEqual(cli.returncode, 2, cli.stderr)
            self.assertEqual(json.loads((root / "audit.json").read_text())["requirement_summary"]["verified"], 3)
            for invalid in ({"unknown": paths}, {"faster.section6_6": paths[:1]},
                            {"faster.section6_6": True},
                            {"faster.section6_6": [paths[0], paths[0]]}):
                with self.subTest(invalid=invalid), self.assertRaises(ValueError):
                    audit.audit(payload, asset_root=root, independent_evidence=invalid)
            with (root / "0/result.json").open("a") as stream:
                stream.write(" ")
            with self.assertRaisesRegex(ValueError, "STALE_RESULT"):
                audit.audit(payload, asset_root=root,
                            independent_evidence={"faster.section6_6": paths})

    def test_faster_discrete_moment_repeats_verify_only_scoped_methods(self):
        """The native fixed-topology derivative gate must be evidence-bound.

        This intentionally exercises the future strict runner rather than
        allowing a passing CTest result to be promoted by a caller boolean.
        It may grant only the derivative and inverse-Gauss construction
        methods, never an author case or an independent Boolean backend.
        """
        audit = load_audit()
        import pybind11
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            payload = blocked_payload(audit, root)
            paths = []
            for number in range(2):
                output = root / str(number)
                completed = subprocess.run([
                    sys.executable, str(STRICT / "run_faster_discrete_moments.py"),
                    "--build-directory", str(ROOT / "build/strict-algorithms-exact-current-20260922"),
                    "--output-directory", str(output),
                    "--dependency-root", "pybind11=" + pybind11.get_include(),
                ], capture_output=True, text=True,
                    env={**os.environ, "PYTHONPATH": os.pathsep.join((str(ROOT / "python"),
                                                                         str(ROOT / "build/strict-algorithms-exact-current-20260922/python"))),
                         "LD_PRELOAD": "/lib/x86_64-linux-gnu/libstdc++.so.6",
                         "PYTHONDONTWRITEBYTECODE": "1"})
                self.assertEqual(completed.returncode, 0, completed.stderr)
                paths.append(str(output / "evidence.json"))
            discrete_evidence = json.loads(Path(paths[0]).read_text(encoding="utf-8"))
            request = discrete_evidence["begin"]["request"]
            self.assertEqual(
                request.get("source_policy"),
                "make-it-stand-explicit-effective-source-v2",
            )
            self.assertNotIn(
                "tests/python/test_spin_it_sources.py",
                discrete_evidence["begin"]["observation"]["source"]["selected_paths"],
            )
            result = audit.audit(payload, asset_root=root,
                                 independent_evidence={"faster.discrete_moments": paths})
            verified = {row["requirement_id"] for row in result["requirements"]
                        if row["status"] == "VERIFIED"}
            self.assertEqual(verified, {
                "faster.discrete_moment_map.differentiable",
                "faster.geometry.inverse_gauss_map",
            })
            self.assertEqual(result["requirement_summary"]["verified"], 2)
            self.assertNotIn("faster.geometry.boolean_independent", verified)
            with (root / "0" / "discrete-moments.stdout").open("a", encoding="utf-8") as stream:
                stream.write(" ")
            with self.assertRaisesRegex(ValueError, "STALE_RESULT"):
                audit.audit(payload, asset_root=root,
                            independent_evidence={"faster.discrete_moments": paths})

    def test_reference_boolean_paraboloid_operation_follows_positive_region(self):
        """A frozen positive-region paraboloid must not be audited as subtraction."""
        audit = load_audit()
        self.assertEqual(
            audit._reference_boolean_expected_operation("paraboloid"),
            "intersection",
        )
        self.assertEqual(
            audit._reference_boolean_expected_operation("two_sheet_hyperboloid"),
            "subtraction",
        )
        with self.assertRaisesRegex(ValueError, "unknown reference Boolean fixture"):
            audit._reference_boolean_expected_operation("unknown")

    def test_reference_boolean_admission_checks_input_operands_not_just_backend_agreement(self):
        """Two agreeing outputs must not pass if both violate the exact input Boolean."""
        audit = load_audit()
        native_dir = ROOT / "build/strict-algorithms-exact-current-20260922/python"
        self.assertTrue(native_dir.is_dir())
        sys.path.insert(0, str(native_dir))
        try:
            import _make_it_stand as native
        finally:
            sys.path.pop(0)

        def box(lo, hi):
            x0, y0, z0 = lo
            x1, y1, z1 = hi
            vertices = ((x0, y0, z0), (x1, y0, z0), (x1, y1, z0), (x0, y1, z0),
                        (x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1))
            faces = ((0, 2, 1), (0, 3, 2), (4, 5, 6), (4, 6, 7),
                     (0, 1, 5), (0, 5, 4), (1, 2, 6), (1, 6, 5),
                     (2, 3, 7), (2, 7, 6), (3, 0, 4), (3, 4, 7))
            return native.make_exact_material_region(vertices, faces)

        design = box((0, 0, 0), (2, 2, 2))
        cutter = box((1, 1, 1), (3, 3, 3))
        valid_intersection = box((1, 1, 1), (2, 2, 2))
        invalid_output = design
        intersection = audit._reference_boolean_operand_conditions(
            native, "intersection", design, cutter, valid_intersection)
        self.assertEqual(intersection["status"], "PASS")
        for operation in ("intersection", "subtraction"):
            invalid = audit._reference_boolean_operand_conditions(
                native, operation, design, cutter, invalid_output)
            self.assertEqual(invalid["status"], "FAIL")
            self.assertTrue(invalid["subset_of_design"])
            self.assertFalse(invalid["quadric_operation_condition"])

    def test_faster_reference_boolean_repeats_record_exact_solid_mismatch(self):
        """Matching moments are insufficient when exact material solids differ.

        The fixture is synthetic and cannot validate an author model or a
        complete Faster application.
        """
        audit = load_audit()
        import pybind11
        probe = ROOT / "build/cherchi-reference-singlethread-20260917/cherchi_probe"
        self.assertTrue(probe.is_file(), "pinned Cherchi probe must be built for this integration test")
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            payload = blocked_payload(audit, root)
            paths = []
            for number in range(2):
                output = root / str(number)
                completed = subprocess.run([
                    sys.executable, str(STRICT / "run_faster_reference_boolean.py"),
                    "--build-directory", str(ROOT / "build/strict-algorithms-exact-current-20260922"),
                    "--cherchi-probe", str(probe),
                    "--output-directory", str(output),
                    "--dependency-root", "pybind11=" + pybind11.get_include(),
                ], capture_output=True, text=True,
                    env={**os.environ, "PYTHONPATH": os.pathsep.join((str(ROOT / "python"),
                                                                         str(ROOT / "build/strict-algorithms-exact-current-20260922/python"))),
                         "LD_PRELOAD": str(ROOT / "build/strict-algorithms-exact-current-20260922/compiler-runtime/libstdc++.so.6"),
                         "PYTHONDONTWRITEBYTECODE": "1"})
                self.assertEqual(completed.returncode, 0, completed.stderr)
                paths.append(str(output / "evidence.json"))
            scientific = json.loads((root / "0" / "scientific.json").read_text(encoding="utf-8"))
            fixtures = {row["name"]: row for row in scientific["fixtures"]}
            self.assertEqual(set(fixtures), {
                "ellipsoid", "one_sheet_hyperboloid", "two_sheet_hyperboloid", "plane", "paraboloid",
            })
            self.assertEqual(
                {name: row["operation"] for name, row in fixtures.items()},
                {"ellipsoid": "subtraction", "one_sheet_hyperboloid": "subtraction",
                 # For these frozen coefficients the two bounded sheets are
                 # the q<0 lobes; the paraboloid is q=x^2+y^2-z, whose capped
                 # solid is the q>0 side. These signs determine minus/intersect.
                 "two_sheet_hyperboloid": "subtraction", "plane": "intersection",
                 "paraboloid": "intersection"},
            )
            # These are necessary conditions on the exact represented output
            # meshes, not a sufficient proof of the requested Boolean solid.
            expected_operand_checks = {
                "ellipsoid": {"cherchi": (True, True), "cgal": (True, True)},
                "one_sheet_hyperboloid": {"cherchi": (False, False), "cgal": (True, False)},
                "two_sheet_hyperboloid": {"cherchi": (False, False), "cgal": (True, False)},
                "plane": {"cherchi": (False, True), "cgal": (True, True)},
                "paraboloid": {"cherchi": (True, False), "cgal": (True, False)},
            }
            expected_export_outside_box = {
                "ellipsoid": 0,
                "one_sheet_hyperboloid": 11,
                "two_sheet_hyperboloid": 7,
                "plane": 1,
                "paraboloid": 0,
            }
            expected_exact_construction_equal = {
                "ellipsoid": True,
                "one_sheet_hyperboloid": False,
                "two_sheet_hyperboloid": False,
                "plane": True,
                "paraboloid": False,
            }
            expected_cgal_matches_exact_input = {
                "ellipsoid": True, "one_sheet_hyperboloid": False,
                "two_sheet_hyperboloid": False, "plane": True,
                "paraboloid": False,
            }
            native_dir = ROOT / "build/strict-algorithms-exact-current-20260922/python"
            sys.path.insert(0, str(native_dir))
            try:
                import _make_it_stand as native
            finally:
                sys.path.pop(0)
            for row in fixtures.values():
                checks = row["geometry_comparison"]["operand_necessary_conditions"]
                box_audit = row["face_source_provenance"]["fixed_box_vertex_audit"]
                self.assertEqual(box_audit["scope"], "AXIS_ALIGNED_FIXED_DESIGN_BOX_VERTICES_ONLY")
                self.assertEqual(box_audit["exact_rational_outside_count"], 0)
                self.assertEqual(box_audit["exported_binary64_outside_count"],
                                 expected_export_outside_box[row["name"]])
                rational_relation = row["geometry_comparison"]["exact_construction_relation"]
                self.assertEqual(rational_relation["scope"],
                                 "cherchi-exact-construction-vs-cgal-exported-polyhedron-only")
                self.assertTrue(rational_relation["construction_region_valid"])
                self.assertTrue(rational_relation["completed"])
                self.assertEqual(rational_relation["equal"],
                                 expected_exact_construction_equal[row["name"]])
                oracle = row["geometry_comparison"]["exact_input_boolean_oracle"]
                self.assertEqual(oracle["scope"],
                                 "exact-Boolean-of-serialized-input-polyhedra-only")
                self.assertTrue(oracle["cherchi_construction"]["completed"])
                self.assertTrue(oracle["cherchi_construction"]["equal"])
                self.assertTrue(oracle["cgal_export"]["completed"])
                self.assertEqual(oracle["cgal_export"]["equal"],
                                 expected_cgal_matches_exact_input[row["name"]])
                self.assertIn("cgal_construction", oracle)
                self.assertTrue(oracle["cgal_construction"]["completed"], oracle)
                internal = oracle["cgal_cherchi_exact_construction_relation"]
                self.assertTrue(internal["completed"], internal)
                self.assertTrue(internal["equal"], internal)
                self.assertEqual(internal["symmetric_difference_volume_exact"], "0/1")
                rounding = oracle["cgal_output_rounding_moments"]
                self.assertEqual(rounding["scope"],
                    "ten-geometric-moments-only;not-exact-solid-equivalence")
                self.assertTrue(rounding["within_1e_6_numeric_gate"], rounding)
                self.assertGreaterEqual(rounding["maximum_absolute_normalized_difference"], 0)
                if row["name"] in ("one_sheet_hyperboloid", "two_sheet_hyperboloid", "paraboloid"):
                    self.assertGreater(rounding["maximum_absolute_normalized_difference"], 0)
                exact_output_path = root / "0" / f"cgal-{row['name']}-exact-construction.json"
                self.assertTrue(exact_output_path.is_file())
                exact_output = json.loads(exact_output_path.read_text(encoding="ascii"))
                self.assertEqual(exact_output["scope"],
                                 "serialized-binary64-input-polyhedra-before-output-rounding")
                self.assertEqual(len(exact_output["birth_faces"]), len(exact_output["faces"]))
                self.assertEqual(oracle, audit._reference_boolean_input_oracle(
                    native, root / "0", row["name"], row["operation"]))
                cgal_vertices, cgal_faces = audit._obj_mesh_arrays(
                    root / "0" / f"cgal-material-{row['name']}.obj")
                cherchi_vertices, cherchi_faces = audit._obj_mesh_arrays(
                    root / "0" / f"cherchi-{row['name']}-run-1.obj")
                independent_rational_relation = audit._reference_boolean_exact_construction_relation(
                    native, root / "0" / f"cherchi-{row['name']}-provenance.obj.exact.tsv",
                    cherchi_vertices, cherchi_faces,
                    native.make_exact_material_region(cgal_vertices, cgal_faces))
                self.assertEqual(independent_rational_relation, rational_relation)
                quadric_vertices, quadric_faces = audit._obj_mesh_arrays(
                    root / "0" / f"quadric-{row['name']}.obj")
                actual_input_digest = native.make_exact_material_region(
                    quadric_vertices, quadric_faces).facts["region_digest"]
                for backend, expected in expected_operand_checks[row["name"]].items():
                    self.assertEqual(checks[backend]["quadric_digest"], actual_input_digest,
                                     "the check must bind the serialized operand read by Cherchi")
                    self.assertTrue(checks[backend]["completed"])
                    self.assertEqual(
                        (checks[backend]["subset_of_design"],
                         checks[backend]["quadric_operation_condition"]),
                        expected,
                    )
                    self.assertEqual(checks[backend]["status"],
                                     "PASS" if all(expected) else "FAIL")
                self.assertEqual(len(row["cherchi_runs"]), 2)
                self.assertLessEqual(row["cherchi_runs"][0]["maximum_absolute_error"], 1e-12)
                provenance = row["face_source_provenance"]
                construction = row["vertex_construction_provenance"]
                self.assertEqual(provenance["geometry_sha256"], row["cherchi_runs"][0]["sha256"])
                if row["name"] == "ellipsoid":
                    self.assertEqual(provenance["status"], "DIAGNOSTIC_ONLY")
                    self.assertEqual(sum(provenance["status_counts"].values()),
                                     row["cherchi_runs"][0]["face_count"])
                else:
                    self.assertEqual(provenance["status"], "APPROXIMATE_OUTPUT_COORDINATES")
                    self.assertGreater(provenance["approximation_disagreement_count"], 0)
                    self.assertEqual(provenance["source_records"], 0)
                self.assertEqual(construction["status"], "DIAGNOSTIC_ONLY")
                self.assertEqual(construction["source_identity_status"], "UNIQUE")
                self.assertEqual(construction["vertex_count"],
                                 row["cherchi_runs"][0]["vertex_count"])
                self.assertEqual(sum(construction["point_type_counts"].values()),
                                 construction["vertex_count"])
                self.assertTrue(construction["all_operands_are_exact_input_vertices"])
                self.assertEqual(construction["surface_plane_source_counts"].get("AMBIGUOUS", 0), 0)
                self.assertEqual(construction["surface_plane_source_counts"].get("UNMATCHED", 0), 0)
            self.assertEqual(fixtures["ellipsoid"]["geometry_comparison"]["status"],
                             "SAME_EXACT_POLYHEDRAL_SOLID")
            self.assertEqual(fixtures["one_sheet_hyperboloid"]["geometry_comparison"]["status"],
                             "DIFFERENT_EXACT_POLYHEDRAL_SOLIDS")
            plane_volume = fixtures["plane"]["geometry_comparison"]["symmetric_difference_volume_comparison"]
            self.assertEqual(plane_volume["status"], "WITHIN_INPUT_SCALE_VOLUME_THRESHOLD_NOT_EXACT_MATCH")
            self.assertFalse(plane_volume["exact_set_equal"])
            curved_volume = fixtures["one_sheet_hyperboloid"]["geometry_comparison"]["symmetric_difference_volume_comparison"]
            self.assertEqual(curved_volume["status"], "WITHIN_INPUT_SCALE_VOLUME_THRESHOLD_NOT_EXACT_MATCH")
            self.assertFalse(curved_volume["exact_set_equal"])
            result = audit.audit(payload, asset_root=root,
                                 independent_evidence={"faster.reference_boolean": paths})
            verified = {row["requirement_id"] for row in result["requirements"]
                        if row["status"] == "VERIFIED"}
            self.assertIn("faster.geometry.boolean_independent", verified)
            self.assertNotIn("faster.geometry.exported_material_boolean", verified)
            self.assertIn("faster.geometry.exported_material_boolean",
                          {row["requirement_id"] for row in result["requirements"]})
            from make_it_stand.verification.run_evidence import PaperRunEvidenceV2, check_freshness
            with (root / "0" / "cherchi-ellipsoid-run-1.obj").open("a", encoding="utf-8") as stream:
                stream.write(" ")
            evidence = PaperRunEvidenceV2.from_json(Path(paths[0]).read_text(encoding="utf-8"))
            self.assertEqual(check_freshness(evidence).status, "STALE_RESULT")

    def test_vertex_construction_provenance_does_not_hide_duplicate_source_ids(self):
        """Exact coordinates cannot make duplicated source identities unique."""
        runner_path = STRICT / "run_faster_reference_boolean.py"
        specification = importlib.util.spec_from_file_location(
            "reference_boolean_runner_duplicate_source", runner_path)
        self.assertIsNotNone(specification)
        self.assertIsNotNone(specification.loader)
        runner = importlib.util.module_from_spec(specification)
        specification.loader.exec_module(runner)
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            source = root / "source.obj"
            source.write_text(
                "v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 3\n", encoding="ascii")
            header = ["vertex_index", "point_type", "operand_count", "parameter_t"] + [
                f"operand{index}_{axis}" for index in range(9) for axis in "xyz"
            ]
            rows = []
            for index, point in enumerate((
                ("0/1", "0/1", "0/1"),
                ("0b1/0b1", "0/1", "0/1"),
                ("0/1", "0b1/0b1", "0/1"),
            )):
                row = [str(index), "EXPLICIT3D", "1", "0", *point]
                row.extend(["-"] * 24)
                rows.append("\t".join(row))
            sidecar = root / "construction.tsv"
            sidecar.write_text("\t".join(header) + "\n" + "\n".join(rows) + "\n",
                               encoding="ascii")

            summary = runner._summarize_vertex_construction_provenance(
                (source, source), source, sidecar)

            self.assertEqual(summary["source_identity_status"], "AMBIGUOUS")
            self.assertTrue(summary["all_operands_are_exact_input_vertices"])
            self.assertEqual(summary["operand_source_counts"],
                             {"UNIQUE": 0, "AMBIGUOUS": 3, "UNMATCHED": 0})

    def test_lnc_provenance_retains_parameter_and_binds_topology_identity(self):
        runner_path = STRICT / "run_faster_reference_boolean.py"
        specification = importlib.util.spec_from_file_location(
            "reference_boolean_runner_lnc_parameter", runner_path)
        self.assertIsNotNone(specification)
        self.assertIsNotNone(specification.loader)
        runner = importlib.util.module_from_spec(specification)
        specification.loader.exec_module(runner)
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            source, output, sidecar = root / "source.obj", root / "output.obj", root / "construction.tsv"
            source.write_text(
                "v 0 0 0\nv 1 0 0\nv 0 1 0\nv 0 0 1\n"
                "f 1 3 2\nf 1 2 4\nf 2 3 4\nf 3 1 4\n", encoding="ascii")
            output.write_text(
                "v 0 0 0\nv 0.25 0 0\nv 0 1 0\nv 0 0 1\n"
                "f 1 3 2\nf 1 2 4\nf 2 3 4\nf 3 1 4\n", encoding="ascii")
            header = ["vertex_index", "point_type", "operand_count", "parameter_t"] + [
                f"operand{index}_{axis}" for index in range(9) for axis in "xyz"
            ]

            def row(index: int, kind: str, count: int, parameter: str, operands: list[tuple[str, str, str]]) -> str:
                fields = [str(index), kind, str(count), parameter]
                flattened = [coordinate for point in operands for coordinate in point]
                fields.extend(flattened + ["-"] * (27 - len(flattened)))
                return "\t".join(fields)

            rows = [
                row(0, "EXPLICIT3D", 1, "0", [("0/1", "0/1", "0/1")]),
                row(1, "LNC", 2, "0.25", [("0/1", "0/1", "0/1"), ("1/1", "0/1", "0/1")]),
                row(2, "EXPLICIT3D", 1, "0", [("0/1", "1/1", "0/1")]),
                row(3, "EXPLICIT3D", 1, "0", [("0/1", "0/1", "1/1")]),
            ]
            sidecar.write_text("\t".join(header) + "\n" + "\n".join(rows) + "\n", encoding="ascii")
            summary = runner._summarize_vertex_construction_provenance((source,), output, sidecar)
            self.assertEqual(summary["records"][1]["parameter_t"], 0.25)
            faces = ((0, 2, 1), (0, 1, 3), (1, 2, 3), (2, 0, 3))
            first = runner._construction_topology_signature(summary, faces)
            changed = dict(summary)
            changed["records"] = [dict(record) for record in summary["records"]]
            changed["records"][1]["parameter_t"] = 0.5
            second = runner._construction_topology_signature(changed, faces)
            self.assertNotEqual(first, second)

    def test_reference_boolean_runner_can_compare_actual_oriented_triangle_boundaries(self):
        runner_path = STRICT / "run_faster_reference_boolean.py"
        specification = importlib.util.spec_from_file_location("reference_boolean_runner_geometry", runner_path)
        self.assertIsNotNone(specification)
        self.assertIsNotNone(specification.loader)
        runner = importlib.util.module_from_spec(specification)
        specification.loader.exec_module(runner)
        self.assertTrue(hasattr(runner, "_triangle_boundary_signature"),
                        "moment equality must be accompanied by an actual mesh-boundary comparison")
        vertices = [(0., 0., 0.), (1., 0., 0.), (0., 1., 0.), (0., 0., 1.)]
        faces = [(0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3)]
        reference = runner._triangle_boundary_signature(vertices, faces)
        reordered = runner._triangle_boundary_signature(vertices, [faces[2], faces[0], faces[3], faces[1]])
        self.assertEqual(reference, reordered, "face ordering must not affect geometric boundary identity")
        changed = list(vertices)
        changed[3] = (0., 0., 1.01)
        self.assertNotEqual(reference, runner._triangle_boundary_signature(changed, faces))
        flipped = list(faces)
        flipped[0] = tuple(reversed(flipped[0]))
        self.assertNotEqual(reference, runner._triangle_boundary_signature(vertices, flipped),
                            "orientation changes must not be hidden by geometry matching")

    def test_reference_boolean_runner_distinguishes_exact_solid_mismatch_from_unknown(self):
        runner_path = STRICT / "run_faster_reference_boolean.py"
        specification = importlib.util.spec_from_file_location("reference_boolean_runner_relation", runner_path)
        self.assertIsNotNone(specification)
        self.assertIsNotNone(specification.loader)
        runner = importlib.util.module_from_spec(specification)
        specification.loader.exec_module(runner)
        self.assertTrue(hasattr(runner, "_exact_solid_relation_status"),
                        "triangle signature mismatch must be classified by an exact solid relation")
        self.assertEqual(runner._exact_solid_relation_status({
            "completed": True, "equal": True, "a_subset_b": True, "b_subset_a": True,
        }), "SAME_EXACT_POLYHEDRAL_SOLID")
        self.assertEqual(runner._exact_solid_relation_status({
            "completed": True, "equal": False, "a_subset_b": False, "b_subset_a": False,
        }), "DIFFERENT_EXACT_POLYHEDRAL_SOLIDS")
        self.assertEqual(runner._exact_solid_relation_status({
            "completed": False, "equal": False, "a_subset_b": False, "b_subset_a": False,
        }), "UNKNOWN")
        self.assertTrue(hasattr(runner, "_symmetric_difference_volume_comparison"))
        within = runner._symmetric_difference_volume_comparison({
            "completed": True, "equal": False,
            "symmetric_difference_volume_exact": f"1/{2**60}",
        }, [(0., 0., 0.), (1., 0., 0.), (0., 1., 0.)],
           [(0., 0., 0.), (1., 0., 0.), (0., 1., 0.)])
        self.assertEqual(within["status"], "WITHIN_INPUT_SCALE_VOLUME_THRESHOLD_NOT_EXACT_MATCH")
        self.assertFalse(within["exact_set_equal"], "volume tolerance must not upgrade exact set inequality")
        outside = runner._symmetric_difference_volume_comparison({
            "completed": True, "equal": False, "symmetric_difference_volume_exact": "1/1",
        }, [(0., 0., 0.), (1., 0., 0.), (0., 1., 0.)],
           [(0., 0., 0.), (1., 0., 0.), (0., 1., 0.)])
        self.assertEqual(outside["status"], "OUTSIDE_INPUT_SCALE_VOLUME_THRESHOLD")
        audit = load_audit()
        independently_checked = audit._exact_symmetric_difference_volume_comparison({
            "completed": True, "equal": False,
            "symmetric_difference_volume_exact": f"1/{2**60}",
        }, [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]],
           [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]])
        self.assertEqual(independently_checked, within)

    def test_faster_stratum_matrix_repeats_verify_only_geometry_scope(self):
        audit = load_audit()
        import pybind11
        executable = ROOT / "build/m1-stratum-matrix-20260922/quadric_stratum_matrix"
        self.assertTrue(executable.is_file())
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            payload = blocked_payload(audit, root)
            paths = []
            for number in range(2):
                output = root / str(number)
                completed = subprocess.run([
                    sys.executable, str(STRICT / "run_faster_stratum_matrix.py"),
                    "--executable", str(executable), "--output-directory", str(output),
                    "--dependency-root", "pybind11=" + pybind11.get_include(),
                ], capture_output=True, text=True,
                    env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"})
                self.assertEqual(completed.returncode, 0, completed.stderr)
                paths.append(str(output / "evidence.json"))
            result = audit.audit(payload, asset_root=root,
                                 independent_evidence={"faster.stratum_matrix": paths})
            verified = {row["requirement_id"] for row in result["requirements"] if row["status"] == "VERIFIED"}
            self.assertEqual(verified, {"faster.geometry.ideal_points", "faster.geometry.caps",
                                        "faster.geometry.degenerate_types"})
            with (root / "0" / "stratum.json").open("a", encoding="utf-8") as stream:
                stream.write(" ")
            with self.assertRaisesRegex(ValueError, "STALE_RESULT"):
                audit.audit(payload, asset_root=root,
                            independent_evidence={"faster.stratum_matrix": paths})

    def test_faster_eq14_newton_repeats_verify_only_solver_scope(self):
        """Equation 14 evidence must execute all protected native modes twice.

        The fixtures are intentionally synthetic.  They establish the strict
        spherical Newton/restart protocol, not success on an author model or
        a global-optimality claim.
        """
        audit = load_audit()
        import pybind11
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            payload = blocked_payload(audit, root)
            paths = []
            for number in range(2):
                output = root / str(number)
                completed = subprocess.run([
                    sys.executable, str(STRICT / "run_faster_eq14_newton.py"),
                    "--build-directory", str(ROOT / "build/strict-current-tree-20260921"),
                    "--output-directory", str(output),
                    "--dependency-root", "pybind11=" + pybind11.get_include(),
                ], capture_output=True, text=True,
                    env={**os.environ, "PYTHONPATH": os.pathsep.join((str(ROOT / "python"),
                                                                         str(ROOT / "build/strict-current-tree-20260921/python"))),
                         "LD_PRELOAD": "/lib/x86_64-linux-gnu/libstdc++.so.6",
                         "PYTHONDONTWRITEBYTECODE": "1"})
                self.assertEqual(completed.returncode, 0, completed.stderr)
                paths.append(str(output / "evidence.json"))
            result = audit.audit(payload, asset_root=root,
                                 independent_evidence={"faster.eq14_newton": paths})
            verified = {row["requirement_id"] for row in result["requirements"]
                        if row["status"] == "VERIFIED"}
            self.assertEqual(verified, {
                "faster.eq14.newton", "faster.eq14.geodesic",
                "faster.eq14.multiplier", "faster.eq14.restart",
            })
            self.assertEqual(result["requirement_summary"]["verified"], 4)
            self.assertNotIn("faster.theorem.optimality", verified)
            with (root / "0" / "newton-default.stdout").open("a", encoding="utf-8") as stream:
                stream.write(" ")
            with self.assertRaisesRegex(ValueError, "STALE_RESULT"):
                audit.audit(payload, asset_root=root,
                            independent_evidence={"faster.eq14_newton": paths})

    def test_faster_moment_foundations_repeats_verify_only_raw_and_duality_scope(self):
        """Raw moments and coefficient pairing do not prove either theorem.

        The runner must execute the native raw-moment and literal cross-term
        regressions twice, bind their transcripts, and grant only the two
        directly checked catalogue requirements.
        """
        audit = load_audit()
        import pybind11
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            payload = blocked_payload(audit, root)
            paths = []
            for number in range(2):
                output = root / str(number)
                completed = subprocess.run([
                    sys.executable, str(STRICT / "run_faster_moment_foundations.py"),
                    "--build-directory", str(ROOT / "build/strict-current-tree-20260921"),
                    "--output-directory", str(output),
                    "--dependency-root", "pybind11=" + pybind11.get_include(),
                ], capture_output=True, text=True,
                    env={**os.environ, "PYTHONPATH": os.pathsep.join((str(ROOT / "python"),
                                                                         str(ROOT / "build/strict-current-tree-20260921/python"))),
                         "LD_PRELOAD": "/lib/x86_64-linux-gnu/libstdc++.so.6",
                         "PYTHONDONTWRITEBYTECODE": "1"})
                self.assertEqual(completed.returncode, 0, completed.stderr)
                paths.append(str(output / "evidence.json"))
            result = audit.audit(payload, asset_root=root,
                                 independent_evidence={"faster.moment_foundations": paths})
            verified = {row["requirement_id"] for row in result["requirements"]
                        if row["status"] == "VERIFIED"}
            self.assertEqual(verified, {"faster.moments.raw", "faster.moments.coefficient_duality"})
            self.assertNotIn("faster.theorem.quadric_reduction", verified)
            self.assertNotIn("faster.theorem.optimality", verified)
            self.assertNotIn("faster.theorem.degenerate_conditions", verified)
            with (root / "0" / "dual-monomial.stdout").open("a", encoding="utf-8") as stream:
                stream.write(" ")
            with self.assertRaisesRegex(ValueError, "STALE_RESULT"):
                audit.audit(payload, asset_root=root,
                            independent_evidence={"faster.moment_foundations": paths})

    def test_faster_theorem_numerical_conditions_are_supporting_not_theorem_verification(self):
        """A finite Eq. (7) witness cannot substitute for either theorem proof.

        The evidence binds both the small condition checker and the canonical
        material replay.  It must remain visible in the audit while leaving
        every theorem catalogue requirement blocked.
        """
        audit = load_audit()
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            payload = blocked_payload(audit, root)
            paths = []
            for number in range(2):
                output = root / str(number)
                completed = subprocess.run([
                    sys.executable, str(STRICT / "run_faster_theorem_numerical_conditions.py"),
                    "--build-directory", str(ROOT / "build/strict-current-tree-20260921"),
                    "--canonical-build-directory", str(ROOT / "build/m1-canonical-theorem-20260922"),
                    "--output-directory", str(output),
                    "--dependency-root", "pybind11=" + __import__("pybind11").get_include(),
                ], capture_output=True, text=True,
                    env={**os.environ, "PYTHONPATH": os.pathsep.join((str(ROOT / "python"),
                                                                         str(ROOT / "build/strict-current-tree-20260921/python"))),
                         "LD_PRELOAD": "/lib/x86_64-linux-gnu/libstdc++.so.6",
                         "PYTHONDONTWRITEBYTECODE": "1"})
                self.assertEqual(completed.returncode, 0, completed.stderr)
                paths.append(str(output / "evidence.json"))
                recorded = json.loads((output / "evidence.json").read_text())
                self.assertEqual(
                    recorded["begin"]["observation"]["source"]["policy"],
                    "make-it-stand-explicit-effective-source-v2",
                )
                self.assertTrue({
                    "experiments/spin_it_reproduction/strict/run_faster_theorem_numerical_conditions.py",
                    "python/make_it_stand/verification/run_evidence.py",
                    "python/make_it_stand/verification/spin_common.py",
                }.issubset(recorded["begin"]["observation"]["source"]["selected_paths"]))
                scientific = json.loads((output / "scientific.json").read_text())
                exact_instances = json.loads((output / "exact-instances.json").read_text())
                self.assertEqual(exact_instances["status"], "EXACT_INSTANCE_CONDITIONS_VERIFIED")
                self.assertEqual(exact_instances["general_theorem_proof"], "NOT_CLAIMED")
                self.assertEqual(scientific["exact_instance_coverage"], "EXACT_INSTANCE_CONDITIONS_VERIFIED")
                theorem_command = next(
                    command for command in scientific["commands"]
                    if command["name"] == "theorem-conditions"
                )
                self.assertIn(
                    "[PASS] finite signed variation certificate proves a strict discrete increase",
                    theorem_command["required_markers"],
                )
            result = audit.audit(
                payload, asset_root=root,
                supporting_evidence={"faster.theorem_numerical_conditions": paths},
            )
            verified = {row["requirement_id"] for row in result["requirements"]
                        if row["status"] == "VERIFIED"}
            self.assertEqual(verified, set())
            self.assertEqual(result["requirement_summary"]["verified"], 0)
            support = result["supporting_evidence"]["faster.theorem_numerical_conditions"]
            self.assertEqual(
                support["scope"],
                "SYNTHETIC_EXACT_RATIONAL_THEOREM_INSTANCES_AND_CANONICAL_EQ7_"
                "MATERIAL_METHODS_ONLY_NOT_GENERAL_THEOREM_PROOF_OR_AUTHOR_CASE",
            )
            self.assertNotIn("faster.theorem.quadric_reduction", verified)
            self.assertNotIn("faster.theorem.optimality", verified)
            self.assertNotIn("faster.theorem.degenerate_conditions", verified)
            with (root / "0" / "exact-instances.json").open("a", encoding="utf-8") as stream:
                stream.write(" ")
            with self.assertRaisesRegex(ValueError, "STALE_RESULT"):
                audit.audit(
                    payload, asset_root=root,
                    supporting_evidence={"faster.theorem_numerical_conditions": paths},
                )
            with (root / "0" / "canonical-replay.stdout").open("a", encoding="utf-8") as stream:
                stream.write(" ")
            with self.assertRaisesRegex(ValueError, "STALE_RESULT"):
                audit.audit(
                    payload, asset_root=root,
                    supporting_evidence={"faster.theorem_numerical_conditions": paths},
                )

    def test_faster_hydrostatic_geometry_repeats_verify_only_method_scope(self):
        audit = load_audit()
        import pybind11
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            payload = blocked_payload(audit, root)
            paths = []
            for number in range(2):
                output = root / str(number)
                completed = subprocess.run([
                    sys.executable, str(STRICT / "run_faster_hydrostatic_geometry.py"),
                    "--build-directory", str(ROOT / "build/strict-current-tree-20260921"),
                    "--output-directory", str(output),
                    "--dependency-root", "pybind11=" + pybind11.get_include(),
                ], capture_output=True, text=True,
                    env={**os.environ, "PYTHONPATH": os.pathsep.join((str(ROOT / "python"),
                                                                         str(ROOT / "build/strict-current-tree-20260921/python"))),
                         "LD_PRELOAD": "/lib/x86_64-linux-gnu/libstdc++.so.6",
                         "PYTHONDONTWRITEBYTECODE": "1"})
                self.assertEqual(completed.returncode, 0, completed.stderr)
                paths.append(str(output / "evidence.json"))
            result = audit.audit(payload, asset_root=root,
                                 independent_evidence={"faster.hydrostatic_geometry": paths})
            verified = {row["requirement_id"] for row in result["requirements"]
                        if row["status"] == "VERIFIED"}
            self.assertEqual(verified, {"faster.hydrostatics.actual"})
            self.assertNotIn("faster.floating_island.chained_assets", verified)
            hydro = next(row for row in result["requirements"]
                         if row["requirement_id"] == "faster.hydrostatics.actual")
            self.assertEqual(
                hydro["details"]["scope"],
                "SYNTHETIC_CLOSED_DISPLACEMENT_BOUNDARY_ACTUAL_WATERLINE_AND_VOLUME_EQUILIBRIUM_METHODS_ONLY_NOT_AUTHOR_FLOATING_CASE",
            )
            with (root / "0" / "hydrostatic.stdout").open("a", encoding="utf-8") as stream:
                stream.write(" ")
            with self.assertRaisesRegex(ValueError, "STALE_RESULT"):
                audit.audit(payload, asset_root=root,
                            independent_evidence={"faster.hydrostatic_geometry": paths})

    def test_faster_planar_methods_repeats_verify_only_method_scope(self):
        """A real five-binary planar suite proves methods, never author cases.

        This catches a broken or omitted protected planar runner, a missing
        polygon/disk or extrusion regression, and any accidental promotion of
        the synthetic fixtures to the three original-experiment rows.
        """
        audit = load_audit()
        import pybind11
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            payload = blocked_payload(audit, root)
            paths = []
            for number in range(2):
                output = root / str(number)
                completed = subprocess.run([
                    sys.executable, str(STRICT / "run_faster_planar_methods.py"),
                    "--build-directory", str(ROOT / "build/strict-algorithms-exact-current-20260922"),
                    "--output-directory", str(output),
                    "--dependency-root", "pybind11=" + pybind11.get_include(),
                ], capture_output=True, text=True,
                    env={**os.environ, "PYTHONPATH": os.pathsep.join((str(ROOT / "python"),
                                                                         str(ROOT / "build/strict-algorithms-exact-current-20260922/python"))),
                         "LD_PRELOAD": "/lib/x86_64-linux-gnu/libstdc++.so.6",
                         "PYTHONDONTWRITEBYTECODE": "1"})
                self.assertEqual(completed.returncode, 0, completed.stderr)
                paths.append(str(output / "evidence.json"))
            result = audit.audit(payload, asset_root=root,
                                 independent_evidence={"faster.planar_methods": paths})
            verified = {row["requirement_id"] for row in result["requirements"]
                        if row["status"] == "VERIFIED"}
            self.assertEqual(verified, {
                "faster.2d.methods",
                "faster.2d.dogleg",
                "faster.2d.polygon_circle",
                "faster.2d.extrusion",
                "faster.application.fidget_spinner",
            })
            self.assertNotIn("faster.table1.flower_fidget_spinner", verified)
            for requirement_id, expected_scope in (
                ("faster.2d.dogleg",
                 "SYNTHETIC_2D_DOGLEG_SOLVER_CONTRACT_ONLY_NOT_AUTHOR_FIDGET_CASE"),
                ("faster.2d.polygon_circle",
                 "SYNTHETIC_POLYGON_DISK_MOMENT_INTEGRATION_CONTRACT_ONLY_NOT_AUTHOR_FIDGET_CASE"),
                ("faster.2d.extrusion",
                 "SYNTHETIC_PLANAR_TO_3D_EXTRUSION_CROSSCHECK_ONLY_NOT_AUTHOR_FIDGET_CASE"),
            ):
                row = next(item for item in result["requirements"]
                           if item["requirement_id"] == requirement_id)
                self.assertEqual(row["details"]["scope"], expected_scope)
                self.assertRegex(row["details"]["scientific_data_sha256"], r"^[0-9a-f]{64}$")
            with (root / "0" / "extrusion.stdout").open("a", encoding="utf-8") as stream:
                stream.write(" ")
            with self.assertRaisesRegex(ValueError, "STALE_RESULT"):
                audit.audit(payload, asset_root=root,
                            independent_evidence={"faster.planar_methods": paths})

    def test_faster_static_application_contract_repeats_verify_only_method_scope(self):
        """Static Section 4.2 evidence verifies its contract, never Table 1."""
        audit = load_audit()
        import pybind11
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            payload = blocked_payload(audit, root)
            paths = []
            for number in range(2):
                output = root / str(number)
                completed = subprocess.run([
                    sys.executable, str(STRICT / "run_faster_static_stability_contract.py"),
                    "--build-directory", str(ROOT / "build/strict-current-tree-20260921"),
                    "--output-directory", str(output),
                    "--dependency-root", "pybind11=" + pybind11.get_include(),
                ], capture_output=True, text=True,
                    env={**os.environ, "PYTHONPATH": os.pathsep.join((str(ROOT / "python"),
                                                                         str(ROOT / "build/strict-current-tree-20260921/python"))),
                         "LD_PRELOAD": "/lib/x86_64-linux-gnu/libstdc++.so.6",
                         "PYTHONDONTWRITEBYTECODE": "1"})
                self.assertEqual(completed.returncode, 0, completed.stderr)
                paths.append(str(output / "evidence.json"))
            result = audit.audit(payload, asset_root=root,
                                 independent_evidence={"faster.static_application_contract": paths})
            verified = {row["requirement_id"] for row in result["requirements"]
                        if row["status"] == "VERIFIED"}
            self.assertEqual(verified, {"faster.application.static_stability"})
            self.assertNotIn("faster.table1.static_tilted_vase", verified)
            self.assertNotIn("faster.application.buoyant_stability", verified)
            with (root / "0" / "static-contract.stdout").open("a", encoding="utf-8") as stream:
                stream.write(" ")
            with self.assertRaisesRegex(ValueError, "STALE_RESULT"):
                audit.audit(payload, asset_root=root,
                            independent_evidence={"faster.static_application_contract": paths})

    def test_faster_buoyant_application_contract_repeats_verify_only_method_scope(self):
        """Section 4.3 validation cannot promote a floating author case."""
        audit = load_audit()
        import pybind11
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            payload = blocked_payload(audit, root)
            buoyant_paths = []
            hydrostatic_paths = []
            for number in range(2):
                hydrostatic_output = root / f"hydrostatic-{number}"
                hydrostatic = subprocess.run([
                    sys.executable, str(STRICT / "run_faster_hydrostatic_geometry.py"),
                    "--build-directory", str(ROOT / "build/strict-current-tree-20260921"),
                    "--output-directory", str(hydrostatic_output),
                    "--dependency-root", "pybind11=" + pybind11.get_include(),
                ], capture_output=True, text=True,
                    env={**os.environ, "PYTHONPATH": os.pathsep.join((str(ROOT / "python"),
                                                                         str(ROOT / "build/strict-current-tree-20260921/python"))),
                         "LD_PRELOAD": "/lib/x86_64-linux-gnu/libstdc++.so.6",
                         "PYTHONDONTWRITEBYTECODE": "1"})
                self.assertEqual(hydrostatic.returncode, 0, hydrostatic.stderr)
                hydrostatic_paths.append(str(hydrostatic_output / "evidence.json"))
                output = root / str(number)
                completed = subprocess.run([
                    sys.executable, str(STRICT / "run_faster_buoyant_stability_contract.py"),
                    "--build-directory", str(ROOT / "build/strict-current-tree-20260921"),
                    "--output-directory", str(output),
                    "--dependency-root", "pybind11=" + pybind11.get_include(),
                ], capture_output=True, text=True,
                    env={**os.environ, "PYTHONPATH": os.pathsep.join((str(ROOT / "python"),
                                                                         str(ROOT / "build/strict-current-tree-20260921/python"))),
                         "LD_PRELOAD": "/lib/x86_64-linux-gnu/libstdc++.so.6",
                         "PYTHONDONTWRITEBYTECODE": "1"})
                self.assertEqual(completed.returncode, 0, completed.stderr)
                scientific = json.loads((output / "scientific.json").read_text(encoding="ascii"))
                required_markers = scientific["command"]["required_markers"]
                self.assertEqual(scientific["command"]["name"], "buoyant-plane-balance-contract")
                self.assertEqual(Path(scientific["command"]["binary"]).name,
                                 "buoyant-plane-balance-regression")
                self.assertIn(
                    "[PASS] prescribed-waterline plane solve and finite-precision Eq. 7 audit passed",
                    required_markers,
                )
                buoyant_paths.append(str(output / "evidence.json"))
            result = audit.audit(payload, asset_root=root,
                                 independent_evidence={
                                     "faster.hydrostatic_geometry": hydrostatic_paths,
                                     "faster.buoyant_application_contract": buoyant_paths,
                                 })
            verified = {row["requirement_id"] for row in result["requirements"]
                        if row["status"] == "VERIFIED"}
            self.assertEqual(verified, {
                "faster.application.buoyant_stability",
                "faster.hydrostatics.actual",
            })
            self.assertNotIn("faster.table1.floating_tea_infuser", verified)
            self.assertNotIn("faster.floating_island.chained_assets", verified)
            with (root / "0" / "buoyant-contract.stdout").open("a", encoding="utf-8") as stream:
                stream.write(" ")
            with self.assertRaisesRegex(ValueError, "STALE_RESULT"):
                audit.audit(payload, asset_root=root,
                            independent_evidence={
                                "faster.hydrostatic_geometry": hydrostatic_paths,
                                "faster.buoyant_application_contract": buoyant_paths,
                            })

    def test_faster_yoyo_application_contract_repeats_verify_only_method_scope(self):
        """Section 4.5 validation cannot promote the Jellyfish Table 1 case."""
        audit = load_audit()
        import pybind11
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            payload = blocked_payload(audit, root)
            paths = []
            for number in range(2):
                output = root / str(number)
                completed = subprocess.run([
                    sys.executable, str(STRICT / "run_faster_yoyo_contract.py"),
                    "--build-directory", str(ROOT / "build/strict-current-tree-20260921"),
                    "--output-directory", str(output),
                    "--dependency-root", "pybind11=" + pybind11.get_include(),
                ], capture_output=True, text=True,
                    env={**os.environ, "PYTHONPATH": os.pathsep.join((str(ROOT / "python"),
                                                                         str(ROOT / "build/strict-current-tree-20260921/python"))),
                         "LD_PRELOAD": "/lib/x86_64-linux-gnu/libstdc++.so.6",
                         "PYTHONDONTWRITEBYTECODE": "1"})
                self.assertEqual(completed.returncode, 0, completed.stderr)
                paths.append(str(output / "evidence.json"))
            result = audit.audit(payload, asset_root=root,
                                 independent_evidence={"faster.yoyo_application_contract": paths})
            verified = {row["requirement_id"] for row in result["requirements"]
                        if row["status"] == "VERIFIED"}
            self.assertEqual(verified, {"faster.application.yoyo"})
            self.assertNotIn("faster.table1.jellyfish_yoyo", verified)
            with (root / "0" / "yoyo-contract.stdout").open("a", encoding="utf-8") as stream:
                stream.write(" ")
            with self.assertRaisesRegex(ValueError, "STALE_RESULT"):
                audit.audit(payload, asset_root=root,
                            independent_evidence={"faster.yoyo_application_contract": paths})

    def test_faster_spinning_top_application_contract_repeats_verify_only_method_scope(self):
        """Section 4.6 validation cannot promote Heart or Noodles Table 1 cases."""
        audit = load_audit()
        import pybind11
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            payload = blocked_payload(audit, root)
            paths = []
            for number in range(2):
                output = root / str(number)
                completed = subprocess.run([
                    sys.executable, str(STRICT / "run_faster_spinning_top_contract.py"),
                    "--build-directory", str(ROOT / "build/strict-current-tree-20260921"),
                    "--output-directory", str(output),
                    "--dependency-root", "pybind11=" + pybind11.get_include(),
                ], capture_output=True, text=True,
                    env={**os.environ, "PYTHONPATH": os.pathsep.join((str(ROOT / "python"),
                                                                         str(ROOT / "build/strict-current-tree-20260921/python"))),
                         "LD_PRELOAD": "/lib/x86_64-linux-gnu/libstdc++.so.6",
                         "PYTHONDONTWRITEBYTECODE": "1"})
                self.assertEqual(completed.returncode, 0, completed.stderr)
                paths.append(str(output / "evidence.json"))
            result = audit.audit(payload, asset_root=root,
                                 independent_evidence={"faster.spinning_top_application_contract": paths})
            verified = {row["requirement_id"] for row in result["requirements"]
                        if row["status"] == "VERIFIED"}
            self.assertEqual(verified, {"faster.application.spinning_top"})
            self.assertNotIn("faster.table1.heart_top", verified)
            self.assertNotIn("faster.table1.noodles_top", verified)
            with (root / "0" / "spinning-top-contract.stdout").open("a", encoding="utf-8") as stream:
                stream.write(" ")
            with self.assertRaisesRegex(ValueError, "STALE_RESULT"):
                audit.audit(payload, asset_root=root,
                            independent_evidence={"faster.spinning_top_application_contract": paths})

    def test_registered_ctest_evidence_requires_two_complete_no_skip_runs(self):
        """A hand-written pass count cannot verify the registered test gate."""
        audit = load_audit()
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            payload = blocked_payload(audit, root)
            paths = []
            for number in range(2):
                output = root / str(number)
                completed = subprocess.run([
                    sys.executable, str(STRICT / "run_registered_ctest_evidence.py"),
                    "--build-directory", str(ROOT / "build/strict-current-tree-20260921"),
                    "--output-directory", str(output), "--dry-run-manifest-only",
                ], capture_output=True, text=True,
                    env={**os.environ, "PYTHONPATH": str(ROOT / "python"), "PYTHONDONTWRITEBYTECODE": "1"})
                self.assertEqual(completed.returncode, 0, completed.stderr)
                paths.append(str(output / "evidence.json"))
            self.assertTrue(all(Path(path).is_file() for path in paths))
            from make_it_stand.verification.run_evidence import PaperRunEvidenceV2
            for path in paths:
                observation = PaperRunEvidenceV2.from_json(Path(path).read_text(encoding="utf-8"))
                self.assertEqual(
                    observation.as_dict()["scientific"]["events"][0]["data"]["mode"],
                    "MANIFEST_ONLY",
                )
            with self.assertRaisesRegex(ValueError, "COMPLETE_EXECUTION"):
                audit.audit(payload, asset_root=root,
                            independent_evidence={"shared.registered_tests": paths})

    def test_registered_ctest_parser_rejects_execution_count_mismatch(self):
        """A pass banner cannot hide an execution shorter than its inventory."""
        runner_path = STRICT / "run_registered_ctest_evidence.py"
        specification = importlib.util.spec_from_file_location("registered_ctest_runner", runner_path)
        self.assertIsNotNone(specification)
        self.assertIsNotNone(specification.loader)
        runner = importlib.util.module_from_spec(specification)
        specification.loader.exec_module(runner)
        inventory = "  Test #1: alpha\n  Test #2: beta\n\nTotal Tests: 2\n"
        truncated_execution = (
            "    Start 1: alpha\n1/1 Test #1: alpha ................ Passed\n"
            "100% tests passed, 0 tests failed out of 1\n"
        )
        with self.assertRaisesRegex(ValueError, "inventory"):
            runner.parse_complete_ctest_output(inventory, truncated_execution)

    def test_pinned_cherchi_probe_serializes_one_sheet_geometry_deterministically(self):
        """Reference-Boolean evidence must not depend on container iteration order."""
        runner_path = STRICT / "run_faster_reference_boolean.py"
        specification = importlib.util.spec_from_file_location("reference_boolean_runner", runner_path)
        self.assertIsNotNone(specification)
        self.assertIsNotNone(specification.loader)
        runner = importlib.util.module_from_spec(specification)
        specification.loader.exec_module(runner)
        probe = ROOT / "build/cherchi-reference-20260917/cherchi_probe"
        self.assertTrue(probe.is_file())
        request, fixtures = runner._fixtures()
        name, _, evaluation = next(row for row in fixtures if row[0] == "one_sheet_hyperboloid")
        operation = "intersection" if evaluation.quadric_solid_is_positive_region else "subtraction"
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            design = root / "design.obj"
            quadric = root / "quadric.obj"
            runner._write_obj(design, request.design_vertices, request.design_faces)
            runner._write_obj(quadric, evaluation.interface_vertices, evaluation.interface_faces)
            outputs = []
            for index in range(2):
                output = root / f"output-{index}.obj"
                completed = subprocess.run(
                    [str(probe), operation, str(design), str(quadric), str(output)],
                    capture_output=True, text=True, check=False,
                )
                self.assertEqual(completed.returncode, 0, completed.stderr)
                outputs.append(output.read_bytes())
            self.assertEqual(outputs[0], outputs[1], name)

    def test_faster_cad_methods_repeats_verify_only_method_scope(self):
        """Pinned OCCT CAD transport proves a method, never the author CAD case."""
        audit = load_audit()
        probe = ROOT / "build/occt-step-roundtrip-20260922a/occt_step_roundtrip"
        occt_source = ROOT / "experiments/occt_step_roundtrip"
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            payload = blocked_payload(audit, root)
            paths = []
            for number in range(2):
                output = root / str(number)
                completed = subprocess.run([
                    sys.executable, str(STRICT / "run_faster_cad_methods.py"),
                    "--probe", str(probe), "--output-directory", str(output),
                    "--dependency-root", "occt-probe-source=" + str(occt_source),
                ], capture_output=True, text=True,
                    env={**os.environ, "PYTHONPATH": os.pathsep.join((str(ROOT / "python"),
                                                                         str(ROOT / "build/strict-current-tree-20260921/python"))),
                         "LD_PRELOAD": "/lib/x86_64-linux-gnu/libstdc++.so.6",
                         "PYTHONDONTWRITEBYTECODE": "1"})
                self.assertEqual(completed.returncode, 0, completed.stderr)
                paths.append(str(output / "evidence.json"))
            result = audit.audit(payload, asset_root=root,
                                 independent_evidence={"faster.cad_methods": paths})
            verified = {row["requirement_id"] for row in result["requirements"]
                        if row["status"] == "VERIFIED"}
            self.assertEqual(verified, {"faster.cad.methods"})
            self.assertNotIn("faster.cad.roundtrip", verified)
            with (root / "0" / "step.step").open("ab") as stream:
                stream.write(b" ")
            with self.assertRaisesRegex(ValueError, "STALE_RESULT"):
                audit.audit(payload, asset_root=root,
                            independent_evidence={"faster.cad_methods": paths})


if __name__ == "__main__":
    unittest.main()
