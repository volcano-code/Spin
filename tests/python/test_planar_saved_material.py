"""Independent saved line/arc extrusion checks, not native moment self-comparison."""
import copy
import importlib
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np


class PlanarSavedMaterialTests(unittest.TestCase):
    def checker(self):
        spec = importlib.util.find_spec("make_it_stand.verification.planar_material")
        self.assertIsNotNone(spec, "saved planar entity independent checker is not implemented")
        return importlib.import_module("make_it_stand.verification.planar_material")

    def test_translated_triangle_prism_density_and_non_diagonal_moments(self):
        checker = self.checker()
        entity = {"geometry_representation": "exact_polygonal_prism", "role": "fixed_material",
                  "loops": [[[3., -2.], [5., -2.], [3., -1.]]],
                  "z_interval": [-.5, 1.5], "volume_density": 3.}
        actual = checker.integrate_planar_entity(entity, origin=(3., -2.), length_scale=5.)
        # Unit-area triangle translated by (3,-2), height 2, density 3.
        np.testing.assert_allclose(actual["raw_mass_moments"],
            [6., 22., -10., 3., 82., 17., 3.5, -5., 11., -37.], rtol=1e-12, atol=1e-12)

    def test_offset_annulus_prism_retains_arc_winding_and_hole(self):
        checker = self.checker()
        def circle(radius, sweep):
            return [{"type": "CIRCULAR_ARC", "center": [1., -2.], "radius": radius,
                     "start_angle": 0., "sweep": sweep,
                     "start": [1.+radius, -2.], "end": [1.+radius, -2.]}]
        entity = {"geometry_representation": "exact_line_arc_profile_swept_between_planes",
                  "role": "editable_design", "boundary_cycles": [circle(2., 2*math.pi),
                     circle(1., -2*math.pi)], "z_interval": [-1., 2.], "volume_density": 2.}
        actual = checker.integrate_planar_entity(entity, origin=(1., -2.), length_scale=4.)
        np.testing.assert_allclose(actual["raw_mass_moments"],
            np.array([18., 18., -36., 9., 40.5, 94.5, 18., -18., 9., -36.])*math.pi,
            rtol=1e-12, atol=1e-11)
        self.assertEqual(actual["quadrature_orders"], [32, 64, 128])
        self.assertFalse(actual["full_topology_verified"])
        altered = copy.deepcopy(entity)
        altered["boundary_cycles"][1][0]["end"][0] += .01
        with self.assertRaisesRegex(ValueError, "endpoint|closure"):
            checker.integrate_planar_entity(altered, origin=(1., -2.), length_scale=4.)

    def test_saved_actual_layered_result_is_reintegrated_and_tampering_rejected(self):
        checker = self.checker()
        from test_spin_it_faster_table1_runner import load_runner, planar_method_case
        runner = load_runner()
        case = runner.run_case(planar_method_case(), starts=1, max_restarts=2)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"application.json"
            path.write_text(json.dumps({"cases": [case]}))
            audit = checker.audit_saved_planar_application(path)
            self.assertEqual(audit["scope"], "SAVED_PLANAR_MOMENTS_AND_BALANCE_NOT_FULL_MATERIAL_OR_PAPER_APPROVAL")
            row = audit["cases"][0]["attempts"][0]
            self.assertTrue(row["moments_match"])
            self.assertTrue(row["balance_constraints_pass"])
            self.assertFalse(audit["manufacturing_approved"])
            self.assertLess(row["maximum_normalized_moment_error"], 1e-10)
            for mutation in ("density", "missing_layer", "claimed_moments", "arc", "total"):
                altered = copy.deepcopy(case)
                result = altered["attempts"][0]["solver_result"]
                if mutation == "density":
                    result["extruded_material_entities"][1]["volume_density"] *= 2
                elif mutation == "missing_layer":
                    result["extruded_material_entities"].pop()
                elif mutation == "claimed_moments":
                    result["extruded_material_entities"][0]["raw_mass_moments"][0] += .01
                elif mutation == "arc":
                    result["extruded_material_entities"][0]["boundary_cycles"][1][0]["radius"] *= .9
                else:
                    result["extruded_complete_material_raw_moments"][0] += .01
                path.write_text(json.dumps({"cases": [altered]}))
                with self.subTest(mutation=mutation), self.assertRaises(ValueError):
                    checker.audit_saved_planar_application(path)

    def test_initial_geometry_sets_normalization_and_failed_run_never_grants_balance(self):
        checker = self.checker()
        from test_spin_it_faster_table1_runner import load_runner, planar_method_case
        case = load_runner().run_case(planar_method_case(maximum_function_evaluations=1),
                                     starts=1, max_restarts=1)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"application.json"
            path.write_text(json.dumps({"cases": [case]}))
            row = checker.audit_saved_planar_application(path)["cases"][0]
            self.assertAlmostEqual(row["normalization"]["length"], math.sqrt(72))
            self.assertAlmostEqual(row["normalization"]["initial_mass"], 39.96)
            self.assertTrue(row["attempts"][0]["moments_match"])
            self.assertFalse(row["attempts"][0]["balance_constraints_pass"])

    def test_fixed_layer_cannot_switch_to_unbound_curve_geometry(self):
        from test_spin_it_faster_table1_runner import load_runner, planar_method_case
        case = load_runner().run_case(planar_method_case(), starts=1, max_restarts=2)
        result = case["attempts"][0]["solver_result"]
        entity = result["extruded_material_entities"][1]
        center = math.pi/4
        rectangle = [[center-.2, -.05], [center+.2, -.05],
                     [center+.2, .05], [center-.2, .05]]
        entity["geometry_representation"] = "exact_line_arc_profile_swept_between_planes"
        entity["boundary_cycles"] = [[{"type": "LINE_SEGMENT", "start": a, "end": b}
            for a, b in zip(rectangle, rectangle[1:]+rectangle[:1])]]
        # Different rectangle: same area and centroid, changed second moments.
        # Leave original input loops in place to expose the representation bypass.
        for raw in (entity["raw_mass_moments"], result["extruded_complete_material_raw_moments"]):
            raw[4] += .02
            raw[5] -= .005
        with self.assertRaisesRegex(ValueError, "fixed material representation"):
            self.checker().audit_planar_case_material(case)

    def test_repeat_identity_binds_saved_material_not_only_its_moments(self):
        from test_spin_it_faster_table1_runner import load_runner, planar_method_case
        case = load_runner().run_case(planar_method_case(), starts=1, max_restarts=2)
        entities = case["attempts"][0]["solver_result"]["extruded_material_entities"]
        encoded = json.dumps(entities, sort_keys=True, separators=(",", ":"),
                             ensure_ascii=True, allow_nan=False).encode("ascii")
        report = self.checker().audit_planar_case_material(case)
        self.assertEqual(report["attempts"][0].get("material_sha256"),
                         hashlib.sha256(encoded).hexdigest())
        altered = copy.deepcopy(case)
        cycle = altered["attempts"][0]["solver_result"]["extruded_material_entities"][0]["boundary_cycles"][0]
        cycle[:] = cycle[1:]+cycle[:1]
        changed = self.checker().audit_planar_case_material(altered)
        np.testing.assert_allclose(report["attempts"][0]["independent_total_raw_mass_moments"],
                                   changed["attempts"][0]["independent_total_raw_mass_moments"], atol=1e-12)
        self.assertNotEqual(report["attempts"][0]["material_sha256"],
                            changed["attempts"][0]["material_sha256"])

    def test_saved_checker_cli_enters_existing_auditor_as_non_promotional_support(self):
        from test_spin_it_faster_table1_runner import load_runner, planar_method_case
        from test_spin_it_strict_audit import load_audit, blocked_payload
        root = Path(__file__).resolve().parents[2]
        cli = root/"experiments/spin_it_reproduction/faster/check_saved_planar_material.py"
        self.assertTrue(cli.is_file(), "actual saved-geometry checker CLI is missing")
        audit = load_audit()
        environment = {**os.environ, "PYTHONPATH": str(root/"python")}
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            evidence_paths = []
            for index in (1, 2):
                case = load_runner().run_case(planar_method_case(), starts=1, max_restarts=2)
                application = directory/f"application-{index}.json"
                application.write_text(json.dumps({"cases": [case]}))
                output = directory/f"check-{index}"
                result = subprocess.run([sys.executable, str(cli), "--application", str(application),
                    "--output", str(output)], env=environment, capture_output=True, text=True, timeout=60)
                self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
                evidence_paths.append(str(output/"evidence.json"))
            baseline = blocked_payload(audit, directory)
            actual = audit.audit(baseline, asset_root=directory,
                supporting_evidence={"faster.planar_saved_material": evidence_paths})
            self.assertEqual(actual["status"], "PARTIAL_BLOCKED")
            self.assertEqual(actual["supporting_evidence"]["faster.planar_saved_material"]["verification_effect"],
                             "SUPPORTING_ONLY_NO_CATALOGUE_REQUIREMENT_VERIFIED")
            with self.assertRaises(ValueError):
                audit._supporting_checks({"faster.planar_saved_material": [evidence_paths[0]]*2})
            changed = json.loads((directory/"application-2.json").read_text())
            changed["cases"][0]["attempts"][0]["solver_result"]["extruded_material_entities"].pop()
            (directory/"application-2.json").write_text(json.dumps(changed))
            with self.assertRaisesRegex(ValueError, "STALE_RESULT"):
                audit._supporting_checks({"faster.planar_saved_material": evidence_paths})


if __name__ == "__main__":
    unittest.main()
