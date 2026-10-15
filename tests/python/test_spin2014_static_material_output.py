"""Actual saved-plane material output, without a manufacturing/seal claim."""
from dataclasses import replace
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from make_it_stand.spin import (SpinItFasterRequest, Spin2014QuadricObjectiveConfig,
    solve_spin2014_static_quadric_hollowing)
from test_spin2014_quadric_interface import box


class CanonicalBoundaryAuthentication(unittest.TestCase):
    vertices = (("0/1", "0/1", "0/1"), ("1/3", "0/1", "0/1"),
                ("0/1", "1/1", "0/1"), ("0/1", "0/1", "1/1"))
    faces = ((0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3))

    def test_reindexing_authenticates_actual_region_once(self):
        # Break caught: building the same expensive Nef twice merely to
        # canonicalize indexing. The spy forwards every actual native call.
        from make_it_stand.verification.faster_material import (
            _canonical_exact_boolean_boundary, _raw_moments)
        from make_it_stand.verification.spin_common import native
        backend = native()
        with patch.object(backend, "make_exact_material_region_from_rationals",
                          wraps=backend.make_exact_material_region_from_rationals) as observed:
            result = _canonical_exact_boolean_boundary(self.vertices, self.faces)
        self.assertTrue(result.exact_region().facts["valid"])
        self.assertAlmostEqual(_raw_moments(result, 1.)[0], 1/18, places=15)
        self.assertEqual(result.rational_vertices,
            (("0/1", "0/1", "0/1"), ("0/1", "0/1", "1/1"),
             ("0/1", "1/1", "0/1"), ("1/3", "0/1", "0/1")))
        self.assertEqual(observed.call_count, 1)

    def test_invalid_faces_still_rejected_before_indexing(self):
        from make_it_stand.verification.faster_material import _canonical_exact_boolean_boundary
        cases = [((0, 2, -1),), ((0, 2, 4),), ((0, 2, True),),
                 ((0, 2, 1.0),), ((0, 2),),
                 np.asarray(self.faces, dtype=float), np.asarray(self.faces, dtype=bool),
                 ((0, 0, 1),), self.faces[:-1]]
        for faces in cases:
            with self.subTest(faces=str(faces)):
                with self.assertRaises((ValueError, TypeError)):
                    _canonical_exact_boolean_boundary(self.vertices, faces)


class StaticMaterialOutput(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        root = Path(__file__).resolve().parents[2]
        spec = importlib.util.spec_from_file_location("static_material_frozen_input",
            root / "experiments/spin2014_objectives_quadric_v1/run.py")
        runner = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(runner)
        request, objective, _, _ = runner.load_frozen_static(
            "static-faceted-shell-lowering", "STATIC_LOWERING", 0)
        cls.static_run = solve_spin2014_static_quadric_hollowing(
            request, objective, fixed_material_density=1.)

    def test_same_retained_complement_forms_actual_closed_material(self):
        # Break caught: reevaluating the opposite plane, which creates a
        # different rounded interface; or replacing the original outer mesh.
        from make_it_stand.verification.faster_material import prepare_spin2014_static_material_output
        output = prepare_spin2014_static_material_output(self.static_run)
        report = output.audit
        self.assertTrue(report["design_partition"]["complete_partition"])
        self.assertTrue(report["design_partition"]["no_positive_volume_overlap"])
        self.assertTrue(report["whole_partition"]["complete_partition"])
        self.assertTrue(report["material_partition"]["complete_partition"])
        self.assertTrue(report["independent_moments_verified"])
        self.assertFalse(report["manufacturing_approved"])
        self.assertFalse(report["sealed_artifact"])
        self.assertGreater(report["cavity_geometric_moments"][0], 0.)
        self.assertAlmostEqual(report["total_moments_target"][0], 8.93933099899741, places=10)
        self.assertLess(report["normalized_moment_error"], 1e-6)
        self.assertIsNotNone(output.cavity.rational_vertices)
        vertices, faces = output.outer.arrays()
        np.testing.assert_array_equal(vertices, self.static_run.measured_request.outer_vertices)
        np.testing.assert_array_equal(faces, self.static_run.measured_request.outer_faces)
        self.assertTrue(output.material.exact_region().facts["valid"])

    def test_exact_complement_serialization_is_invariant_to_indexing(self):
        # Break caught: Boolean traversal/vertex numbering leaking into the
        # persisted scientific identity. No coordinates/facets may change.
        from make_it_stand.verification.faster_material import _canonical_exact_boolean_boundary, _raw_moments
        vertices = (("0/1","0/1","0/1"),("1/3","0/1","0/1"),
                    ("0/1","1/1","0/1"),("0/1","0/1","1/1"))
        faces = ((0,2,1),(0,1,3),(0,3,2),(1,2,3))
        one = _canonical_exact_boolean_boundary(vertices,faces)
        # Fixed reindexing 0->2, 1->0, 2->3, 3->1; reorder faces and
        # rotate each face cyclically, never reverse its orientation.
        two = _canonical_exact_boolean_boundary((vertices[1],vertices[3],vertices[0],vertices[2]),
            ((3,1,0),(1,3,2),(0,1,2),(3,0,2)))
        self.assertEqual(one.identity,two.identity)
        self.assertAlmostEqual(_raw_moments(one,1.)[0],1/18,places=15)
        self.assertEqual(one.vertex_count,4)
        self.assertEqual(one.face_count,4)

    def test_actual_distinct_densities_and_empty_cavity(self):
        # Break caught: density multiplied twice, omitted shell, or treating
        # a full retained design as an invalid empty triangulation.
        from make_it_stand.verification.faster_material import prepare_spin2014_static_material_output
        outer, faces = box(0., 1.)
        outer = tuple((x, y, 4*z) for x, y, z in outer)
        design, inner = box(.1, .9)
        request = SpinItFasterRequest(outer, faces, design, inner, mode="STATIC_BALANCE",
            material_density=3., target_origin=(.5, .5, 0.))
        objective = Spin2014QuadricObjectiveConfig(objective="STATIC_LOWERING", inertia_weight=0.)
        run = solve_spin2014_static_quadric_hollowing(request, objective, fixed_material_density=2.)
        output = prepare_spin2014_static_material_output(run)
        self.assertTrue(output.cavity.explicit_empty)
        self.assertEqual(output.material.identity, output.outer.identity)
        self.assertAlmostEqual(output.audit["fixed_moments_model"][0], 6.976, places=12)
        self.assertAlmostEqual(output.audit["retained_moments_model"][0], 1.536, places=12)
        self.assertAlmostEqual(output.audit["total_moments_target"][0], 8.512, places=12)
        self.assertAlmostEqual(output.audit["total_moments_target"][3], 16.256, places=11)

    def test_forged_density_moments_and_missing_retained_material_rejected(self):
        from make_it_stand.verification.faster_material import prepare_spin2014_static_material_output
        with self.assertRaisesRegex(ValueError, "moment|density"):
            prepare_spin2014_static_material_output(replace(self.static_run, fixed_material_density=2.))
        candidate = dict(self.static_run.candidate)
        candidate["total_raw_moments"] = [99.] + list(candidate["total_raw_moments"])[1:]
        with self.assertRaisesRegex(ValueError, "moment"):
            prepare_spin2014_static_material_output(replace(self.static_run, candidate=candidate))
        candidate = dict(self.static_run.candidate)
        candidate["variable_material_vertices"] = []
        candidate["variable_material_faces"] = []
        with self.assertRaises(ValueError):
            prepare_spin2014_static_material_output(replace(self.static_run, candidate=candidate))

    def test_unqualified_candidate_cannot_become_valid_material_output(self):
        from make_it_stand.verification.faster_material import prepare_spin2014_static_material_output
        candidate = dict(self.static_run.candidate)
        candidate["status"] = "LINE_SEARCH_FAILED"
        with self.assertRaisesRegex(ValueError, "terminal|converged"):
            prepare_spin2014_static_material_output(replace(self.static_run, candidate=candidate))
        square = dict(self.static_run.squared_height)
        square["squared_height_conditions_verified"] = False
        with self.assertRaisesRegex(ValueError, "square|squared"):
            prepare_spin2014_static_material_output(replace(self.static_run, squared_height=square))

    def test_saved_geometry_reloads_and_cannot_become_manufacturing_seal(self):
        # Break caught: output remaining memory-only, or reload accepting an
        # altered cavity/false audit instead of the actual source run/material.
        from make_it_stand.verification.faster_material import prepare_spin2014_static_material_output
        from make_it_stand.verification.spin_bundle import (
            save_spin2014_static_material_output, load_spin2014_static_material_output)
        output = prepare_spin2014_static_material_output(self.static_run)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)/"material"
            saved = save_spin2014_static_material_output(output, self.static_run, root)
            self.assertEqual(saved.scientific_digest, output.scientific_digest)
            reloaded = load_spin2014_static_material_output(root,
                expected_source_run_digest=output.identity["source_run_sha256"])
            self.assertEqual(reloaded.cavity.identity, output.cavity.identity)
            self.assertFalse(reloaded.audit["sealed_artifact"])
            with self.assertRaises(FileExistsError):
                save_spin2014_static_material_output(output, self.static_run, root)
            with self.assertRaisesRegex(ValueError, "binding|source"):
                load_spin2014_static_material_output(root, expected_source_run_digest="0"*64)
            manifest = json.loads((root/"manifest.json").read_text())
            manifest["identity"]["audit"]["manufacturing_approved"] = True
            (root/"manifest.json").write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError, "identity|audit"):
                load_spin2014_static_material_output(root,
                    expected_source_run_digest=output.identity["source_run_sha256"])


if __name__ == "__main__":
    unittest.main()
