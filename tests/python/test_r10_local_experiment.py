"""Contract tests for the license-safe R10 local experiment driver."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "r10_local_experiment",
    REPOSITORY / "python" / "r10_local_experiment.py",
)
assert SPEC is not None and SPEC.loader is not None
DRIVER = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = DRIVER
SPEC.loader.exec_module(DRIVER)


class _FakeResult:
    succeeded = True
    manufacturing_ready = False
    status = "stable_after_carving"

    def __init__(self) -> None:
        self.diagnostics = {
            "r10_metrics": {
                "final_projected_target_energy": 0.125,
                "final_cavity_volume": 0.0,
                "final_cavity_volume_fraction": 0.0,
                "final_outer_volume": 2.0,
                "material_mass_change_fraction": 0.2,
                "final_unscaled_shape_energy": 0.75,
                "normalized_rms_vertex_displacement": 0.02,
                "normalized_maximum_vertex_displacement": 0.04,
                "maximum_user_handle_scale_deviation": 0.03,
                "committed_event_count": 3,
            },
            "targets": {"all_reached": True},
            "manufacturing": {
                "outer_boundary": {
                    "empty": False,
                    "matrix_shape_valid": True,
                    "finite_vertices": True,
                    "face_indices_valid": True,
                    "nondegenerate_faces": True,
                    "unreferenced_vertices": 0,
                    "edge_manifold": True,
                    "vertex_manifold": True,
                    "closed": True,
                    "consistently_oriented": True,
                    "component_orientations_valid": True,
                    "self_intersection_completed": True,
                    "self_intersection_free": True,
                    "signed_volume": 2.0,
                },
                "inner_boundary": {"empty": True},
                "wall_clearance": {"measured_minimum_clearance": None},
            },
            "carving_iterations": [{}, {}],
            "deformation_iterations": [{}],
        }


class R10LocalExperimentContractTests(unittest.TestCase):
    def test_matrix_covers_required_algorithm_and_robustness_axes(self) -> None:
        cases = DRIVER.default_cases()
        identifiers = {case.experiment_id for case in cases}
        self.assertEqual(len(identifiers), len(cases))
        self.assertTrue(any(case.enable_carving and case.enable_deformation and case.enable_laplacian for case in cases))
        self.assertTrue(any(case.enable_carving and not case.enable_deformation for case in cases))
        self.assertTrue(any(not case.enable_carving and case.enable_deformation for case in cases))
        self.assertTrue(any(not case.enable_laplacian for case in cases))
        self.assertTrue(
            any(
                not case.enable_carving
                and case.enable_deformation
                and not case.enable_laplacian
                for case in cases
            )
        )
        self.assertGreaterEqual(len({case.resolution for case in cases}), 3)
        self.assertGreaterEqual(len({case.handle_layout for case in cases}), 2)
        self.assertTrue(any(case.constraint == "suspension" for case in cases))

    def test_manifests_are_canonical_license_safe_and_path_portable(self) -> None:
        from make_it_stand.r10 import read_canonical_json, validate_manifest

        provenance = DRIVER.RepositoryState("a" * 40, False)
        case = DRIVER.default_cases()[0]
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "evidence" / case.experiment_id
            manifest = DRIVER.build_manifest(
                case,
                experiment_root=root,
                repository=REPOSITORY,
                provenance=provenance,
                repeats=2,
            )
            validate_manifest(manifest)
            self.assertEqual(read_canonical_json(root / "manifest.json"), manifest)
            serialized = (root / "manifest.json").read_text(encoding="utf-8")
            self.assertNotIn(str(REPOSITORY), serialized)
            self.assertNotIn("make-it-stand_src_beta", serialized)
            self.assertFalse(manifest["claim"]["exact_reproduction_claim"])
            self.assertEqual(manifest["paper_reference"]["reported_metrics"], [])
            self.assertTrue((root / "inputs" / "LICENSE.txt").is_file())

    def test_metric_extraction_never_invents_absent_wall_measurement(self) -> None:
        metrics = dict(DRIVER.extract_metrics(_FakeResult(), _FakeResult().diagnostics))
        self.assertEqual(metrics["projected_target_energy"], 0.125)
        self.assertEqual(metrics["final_cavity_volume_fraction"], 0.0)
        self.assertEqual(metrics["material_mass_change_fraction"], 0.2)
        self.assertEqual(metrics["carving_iteration_count"], 2.0)
        self.assertEqual(metrics["deformation_iteration_count"], 1.0)
        self.assertEqual(metrics["algorithm_success"], 1.0)
        self.assertNotIn("measured_minimum_wall_clearance", metrics)

    def test_invalid_or_null_diagnostic_measurements_fail_closed(self) -> None:
        result = _FakeResult()
        result.diagnostics["r10_metrics"]["material_mass_change_fraction"] = None
        with self.assertRaisesRegex(ValueError, "not numeric"):
            DRIVER.extract_metrics(result, result.diagnostics)

    def test_empty_inner_boundary_rejects_nonzero_cavity(self) -> None:
        result = _FakeResult()
        result.diagnostics["r10_metrics"]["final_cavity_volume"] = 0.5
        result.diagnostics["r10_metrics"]["final_cavity_volume_fraction"] = 0.25
        with self.assertRaisesRegex(ValueError, "empty inner boundary"):
            DRIVER.extract_metrics(result, result.diagnostics)

    def test_valid_boundary_rejects_missing_cavity_metrics(self) -> None:
        result = _FakeResult()
        result.diagnostics["r10_metrics"]["final_cavity_volume_fraction"] = None
        with self.assertRaisesRegex(ValueError, "cavity metrics are missing"):
            DRIVER.extract_metrics(result, result.diagnostics)

    def test_nonempty_boundary_requires_strict_audit_and_consistent_volume(self) -> None:
        result = _FakeResult()
        result.diagnostics["r10_metrics"].update(
            final_cavity_volume=0.5,
            final_cavity_volume_fraction=0.25,
            final_outer_volume=2.0,
        )
        result.diagnostics["manufacturing"]["inner_boundary"] = {
            "empty": False,
            "matrix_shape_valid": True,
            "finite_vertices": True,
            "face_indices_valid": True,
            "nondegenerate_faces": True,
            "unreferenced_vertices": 1,
            "edge_manifold": True,
            "vertex_manifold": True,
            "closed": True,
            "consistently_oriented": True,
            "component_orientations_valid": True,
            "self_intersection_completed": True,
            "self_intersection_free": True,
            "signed_volume": -0.5,
        }
        metrics = dict(DRIVER.extract_metrics(result, result.diagnostics))
        self.assertNotIn("final_cavity_volume_fraction", metrics)

        result.diagnostics["manufacturing"]["inner_boundary"]["unreferenced_vertices"] = 0
        result.diagnostics["manufacturing"]["inner_boundary"]["signed_volume"] = -0.4
        with self.assertRaisesRegex(ValueError, "signed volume"):
            DRIVER.extract_metrics(result, result.diagnostics)


if __name__ == "__main__":
    unittest.main()
