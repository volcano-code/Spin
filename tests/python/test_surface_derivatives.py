"""Independent numerical contracts for published surface Algorithms 1 and 2."""
import importlib
import json
import unittest
from unittest.mock import patch

import numpy as np

from make_it_stand.verification.mass_moments import DensityRegion, TriangleMesh
from test_independent_mass_moments import TETRAHEDRON, box_mesh, combine_meshes, transformed_mesh


ORDER = ("M", "LX", "LY", "LZ", "QXX", "QYY", "QZZ", "QYZ", "QXZ", "QXY")
ALGORITHM = "spin-it-2014-supplement-algorithms-1-2-surface-polynomial"


class PublishedSurfaceDerivativeTests(unittest.TestCase):
    def setUp(self):
        self.native = importlib.import_module("_make_it_stand")

    def analytic(self, mesh=TETRAHEDRON, density=1.0):
        self.assertTrue(hasattr(self.native, "differentiate_surface_raw_moments"),
                        "the actual published surface differential binding is missing")
        return self.native.differentiate_surface_raw_moments(mesh.vertices, mesh.faces, density)

    def verifier(self):
        name = "make_it_stand.verification.surface_derivatives"
        self.assertIsNotNone(importlib.util.find_spec(name), "independent derivative checker is missing")
        return importlib.import_module(name)

    def test_hand_derived_asymmetric_tetrahedron_and_identity(self):
        out = self.analytic()
        self.assertEqual(tuple(out["order"]), ORDER)
        self.assertEqual(out["algorithm"], ALGORITHM)
        np.testing.assert_allclose(out["moments"], [1, .25, .5, .75, .1, .4, .9, .3, .15, .1], atol=1e-14)
        j = out["jacobian"]
        self.assertEqual(j.shape, (10, 12))
        np.testing.assert_allclose(j[0], [-1, -.5, -1/3, 1, 0, 0, 0, .5, 0, 0, 0, 1/3], atol=1e-14)
        np.testing.assert_allclose(j[9, :3], [0, 0, -1/30], atol=1e-14)
        np.testing.assert_allclose(j[8, :3], [0, -3/40, 0], atol=1e-14)
        np.testing.assert_allclose(j[7, :3], [-.3, 0, 0], atol=1e-14)

    def test_every_entry_three_steps_geometries_and_density_regions(self):
        checker = self.verifier()
        outer = box_mesh((-2, -2, -2), (2, 2, 2))
        cavity = box_mesh((-.8, -.5, -.4), (-.2, .3, .7), inward=True)
        cavity2 = box_mesh((.3, -.7, -.8), (1.1, .2, -.2), inward=True)
        rotated = transformed_mesh(TETRAHEDRON, lambda v: (-v[1]+4, v[0]-3, v[2]+2))
        cases = {
            "tetrahedron": (DensityRegion(TETRAHEDRON, 1),),
            "rotated_translated": (DensityRegion(rotated, 2.7),),
            "cube": (DensityRegion(outer, 1),),
            "eccentric_cavity": (DensityRegion(combine_meshes(outer, cavity), 1.8),),
            "multiple_cavities": (DensityRegion(combine_meshes(outer, cavity, cavity2), .7),),
            "density_regions": (DensityRegion(TETRAHEDRON, 2), DensityRegion(rotated, 7)),
            "thin": (DensityRegion(box_mesh((0, 0, 0), (1, .7, .001)), 3),),
        }
        for name, regions in cases.items():
            with self.subTest(name=name):
                result = checker.verify_surface_derivatives(regions)
                self.assertEqual(result["status"], "consistent", json.dumps({
                    "moment_error": result["moment_max_absolute_error"],
                    "steps": [{k: s[k] for k in ("h", "max_absolute_error", "max_smooth_relative_error", "max_near_zero_absolute_error")} for s in result["steps"]],
                    "trends": [t["counts"] for t in result["trends"]],
                }))
                self.assertEqual(len(result["steps"]), 3)
                columns = 3 * sum(len(r.mesh.vertices) for r in regions)
                for step in result["steps"]:
                    self.assertEqual(np.asarray(step["absolute_errors"]).shape, (10, columns))
                    self.assertTrue(step["all_entries_consistent"])
                self.assertEqual(len(result["trends"]), 2)
                self.assertTrue(result["identity"]["native_module"]["sha256"])
                json.dumps(result, allow_nan=False)
                print("SURFACE_MEASUREMENT", name, json.dumps({
                    "moment_max_absolute_error": result["moment_max_absolute_error"],
                    "steps": [{k: s[k] for k in ("h", "max_absolute_error", "max_smooth_relative_error", "max_near_zero_absolute_error")} for s in result["steps"]],
                    "trends": [t["counts"] for t in result["trends"]],
                }))

    def test_physical_scale_density_translation_covariance(self):
        base = self.analytic()
        degrees = np.array([3, 4, 4, 4, 5, 5, 5, 5, 5, 5])
        for scale, density in ((1e-4, .03), (17, 9), (1e4, 1e-3)):
            mesh = transformed_mesh(TETRAHEDRON, lambda v: tuple(scale*x for x in v))
            out = self.analytic(mesh, density)
            np.testing.assert_allclose(out["moments"] / (density*scale**degrees), base["moments"], rtol=1e-12, atol=1e-13)
            np.testing.assert_allclose(out["jacobian"] / (density*scale**(degrees-1))[:, None], base["jacobian"], rtol=1e-12, atol=1e-13)
        q = transformed_mesh(TETRAHEDRON, lambda v: (v[0]+4, v[1]-3, v[2]+2))
        out = self.analytic(q)
        # Independent binomial translation map for raw moments/Jacobian rows.
        shift = np.array([4., -3., 2.])
        translation = np.eye(10)
        translation[1:4, 0] = shift
        for row, (a, b) in enumerate(((0,0),(1,1),(2,2),(1,2),(0,2),(0,1)), 4):
            translation[row, 0] = shift[a]*shift[b]
            translation[row, 1+b] += shift[a]
            translation[row, 1+a] += shift[b]
        np.testing.assert_allclose(out["moments"], translation @ base["moments"], rtol=1e-12, atol=1e-13)
        np.testing.assert_allclose(out["jacobian"], translation @ base["jacobian"], rtol=1e-12, atol=1e-13)
        self.assertTrue(np.all(np.abs(out["moments"][-3:]) > 0))
        j = out["jacobian"]
        sums = j.reshape(10, 4, 3).sum(axis=1)
        expected = np.zeros((10, 3)); expected[1:4] = np.eye(3)
        for row, (a, b) in enumerate(((0,0),(1,1),(2,2),(1,2),(0,2),(0,1)), 4):
            expected[row, a] += out["moments"][1+b]
            expected[row, b] += out["moments"][1+a]
        np.testing.assert_allclose(sums, expected, atol=2e-13)
        np.testing.assert_allclose(j @ np.asarray(q.vertices).ravel(), degrees*out["moments"], atol=3e-12)

    def test_inward_components_keep_signed_contributions(self):
        inward = TriangleMesh(TETRAHEDRON.vertices, tuple((a,c,b) for a,b,c in TETRAHEDRON.faces))
        positive, negative = self.analytic(), self.analytic(inward)
        np.testing.assert_allclose(negative["moments"], -positive["moments"], atol=1e-14)
        np.testing.assert_allclose(negative["jacobian"], -positive["jacobian"], atol=1e-14)

    def test_invalid_native_geometry_and_density_fail(self):
        self.analytic()
        for v, f, rho in (([], [], 1), (TETRAHEDRON.vertices, TETRAHEDRON.faces[:-1], 1),
                           (np.empty((0, 3)), np.empty((0, 3), dtype=int), 1),
                           ([(0,0,0)]*4, TETRAHEDRON.faces, 1),
                           ([(float("nan"),0,0)]*4, TETRAHEDRON.faces, 1),
                           (TETRAHEDRON.vertices, np.asarray(TETRAHEDRON.faces, float), 1),
                           (TETRAHEDRON.vertices, TETRAHEDRON.faces, float("nan")),
                           (TETRAHEDRON.vertices, ((0, 2, 9), *TETRAHEDRON.faces[1:]), 1),
                           (TETRAHEDRON.vertices, ((0, 1, 2), *TETRAHEDRON.faces[1:]), 1),
                           (TETRAHEDRON.vertices, TETRAHEDRON.faces, 0)):
            with self.subTest(v=v, rho=rho), self.assertRaises((ValueError, TypeError)):
                self.native.differentiate_surface_raw_moments(v, f, rho)

    def test_corrupted_native_data_and_absent_backend_never_consistent(self):
        checker = self.verifier(); self.analytic()
        original = self.native.differentiate_surface_raw_moments
        def corrupt(kind):
            def call(*args):
                out = original(*args)
                if kind == "cross": out["jacobian"][9, 2] += .01
                if kind == "density": out["jacobian"] *= 2
                if kind == "nan": out["jacobian"][0, 0] = float("nan")
                if kind == "empty": out["jacobian"] = np.zeros((10, 0))
                if kind == "order": out["order"] = list(reversed(out["order"]))
                return out
            return call
        for kind in ("cross", "density", "nan", "empty", "order"):
            with self.subTest(kind=kind), patch.object(self.native, "differentiate_surface_raw_moments", corrupt(kind)):
                if kind in ("nan", "empty", "order"):
                    with self.assertRaisesRegex(ValueError, "native"):
                        checker.verify_surface_derivatives((DensityRegion(TETRAHEDRON, 2),))
                else:
                    result = checker.verify_surface_derivatives((DensityRegion(TETRAHEDRON, 2),))
                    self.assertEqual(result["status"], "inconsistent")
        with patch.object(checker, "_native_module", side_effect=RuntimeError("missing backend")):
            with self.assertRaises(RuntimeError):
                checker.verify_surface_derivatives((DensityRegion(TETRAHEDRON, 1),))
        with self.assertRaises(ValueError): checker.verify_surface_derivatives(())

    def test_one_lucky_step_does_not_override_other_steps(self):
        checker = self.verifier(); self.analytic()
        original = checker._independent_vector
        calls = 0
        def corrupt(mesh, density):
            nonlocal calls
            calls += 1
            out = original(mesh, density)
            # Baseline call, then 24 calls per step on the 4-vertex mesh.
            if 26 <= calls <= 49: out[9] += (.01 if calls % 2 == 0 else -.01)
            return out
        with patch.object(checker, "_independent_vector", corrupt):
            result = checker.verify_surface_derivatives((DensityRegion(TETRAHEDRON, 1),))
        self.assertTrue(result["steps"][0]["all_entries_consistent"])
        self.assertFalse(result["steps"][1]["all_entries_consistent"])
        self.assertEqual(result["status"], "inconsistent")

    def test_trend_assessment_rejects_constant_above_floor_errors(self):
        checker = self.verifier()
        classify = checker.assess_error_trend
        self.assertEqual(classify(4e-7, 1e-7, 1e-10, 2e-10), "SECOND_ORDER")
        self.assertEqual(classify(1e-7, 1e-7, 1e-10, 2e-10), "INCONSISTENT")
        self.assertEqual(classify(0, 0, 1e-10, 2e-10), "NOISE_LIMITED")
        self.assertEqual(classify(4e-10, 1e-10, 1e-10, 2e-10), "FLOOR_TRANSITION")
        with self.assertRaises(ValueError): classify(float("nan"), 1, 1, 1)


if __name__ == "__main__": unittest.main()
