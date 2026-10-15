from __future__ import annotations

import json
import hashlib
import math
import os
from pathlib import Path
import re
import subprocess
import tempfile
import unittest

import numpy as np
import make_it_stand
import make_it_stand.spin as spin_api

from make_it_stand.spin import (
    SpinItFasterRequest,
    solve_spin_it_faster_buoyant_plane,
    with_measured_buoyancy_geometry,
)
from make_it_stand.verification.faster_material import (
    quadric_coefficients_to_cad_frame,
    evaluate_and_seal_faster_candidate,
)
from make_it_stand.verification.spin_problem import SpinArtifactProblem


ROOT = Path(__file__).resolve().parents[2]
PROBE = Path(os.environ.get(
    "MAKE_IT_STAND_OCCT_STEP_PROBE",
    ROOT / "build" / "occt-step-roundtrip-20260922a" / "occt_step_roundtrip",
))


class OcctStepRoundtripTests(unittest.TestCase):
    @staticmethod
    def _write_box_obj(path):
        path.write_text(
            "v -2 -2 -2\nv 2 -2 -2\nv 2 2 -2\nv -2 2 -2\n"
            "v -2 -2 2\nv 2 -2 2\nv 2 2 2\nv -2 2 2\n"
            "f 1 3 2\nf 1 4 3\nf 5 6 7\nf 5 7 8\n"
            "f 1 2 6\nf 1 6 5\nf 4 8 7\nf 4 7 3\n"
            "f 1 5 8\nf 1 8 4\nf 2 3 7\nf 2 7 6\n",
            encoding="ascii")

    @staticmethod
    def _write_oriented_box_obj(path, center, axes, half_extents):
        signs = np.asarray(((-1., -1., -1.), (1., -1., -1.),
                            (1., 1., -1.), (-1., 1., -1.),
                            (-1., -1., 1.), (1., -1., 1.),
                            (1., 1., 1.), (-1., 1., 1.)))
        vertices = np.asarray(center) + (axes @ (signs * np.asarray(half_extents)).T).T
        faces = ((1,3,2),(1,4,3),(5,6,7),(5,7,8),(1,2,6),(1,6,5),
                 (4,8,7),(4,7,3),(1,5,8),(1,8,4),(2,3,7),(2,7,6))
        lines = [*(f"v {x:.17g} {y:.17g} {z:.17g}" for x, y, z in vertices),
                 *("f " + " ".join(str(index) for index in face) for face in faces)]
        path.write_text("\n".join(lines) + "\n", encoding="ascii")

    @staticmethod
    def _raw_from_center_covariance(volume, center, covariance):
        second = covariance + np.outer(center, center)
        return np.asarray((volume, *(volume * center),
                          *(volume * np.diag(second)),
                          volume * second[1, 2], volume * second[0, 2],
                          volume * second[0, 1]))

    def test_layered_planar_solver_result_exports_exact_multisolid_step(self):
        """Exact solved arcs and density layers must survive STEP with their moments."""
        export = getattr(make_it_stand, "export_spin_it_faster_planar_material_step", None)
        self.assertTrue(callable(export), "the exact 2D result-to-CAD handoff is missing")
        pi = math.pi
        outer = ((-3,-3),(3,-3),(3,3),(-3,3))
        patch = ((pi/4-0.1,-0.1),(pi/4+0.1,-0.1),
                 (pi/4+0.1,0.1),(pi/4-0.1,0.1))
        result = spin_api.solve_spin_it_faster_planar_fidget(
            design_region_loops=(outer, tuple(reversed(patch))),
            fixed_material_loops=(), fixed_areal_density=None,
            design_areal_density=1.0, initial_center=(0.75,0.2),
            extrusion_interval=(-0.7,0.8), extrusion_volume_density=2.0/3.0,
            fixed_extrusion_layers=(
                {"loops": (patch,), "z_interval": (-0.1,0.0), "volume_density": 500.0},
                {"loops": (patch,), "z_interval": (0.0,0.1), "volume_density": 500.0},
            ),
        )
        self.assertEqual(result["status"], "CONVERGED")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            brep, step = root / "layered.brep", root / "layered.step"
            cad = export(result, occt_executable=PROBE, brep_path=brep, step_path=step)
            self.assertEqual(cad["status"], "PLANAR_LAYERED_STEP_BREP_ROUNDTRIP_VALIDATED")
            self.assertEqual(cad["material_solid_count"], 3)
            self.assertEqual(cad["entity_roles"], ["editable_design", "fixed_material", "fixed_material"])
            np.testing.assert_allclose(cad["entity_densities"], (2.0/3.0,500.0,500.0), rtol=0, atol=0)
            entity_ids = cad["material_entity_ids"]
            self.assertEqual(len(entity_ids), 3)
            self.assertEqual(len(set(entity_ids)), 3)
            self.assertTrue(all(re.fullmatch(r"[0-9a-f]{64}", identity)
                                for identity in entity_ids))
            identities = cad["material_identity_sidecar"]
            self.assertEqual([record["entity_id"] for record in identities], entity_ids)
            self.assertEqual([record["role"] for record in identities], cad["entity_roles"])
            np.testing.assert_array_equal(
                [record["volume_density"] for record in identities], cad["entity_densities"])
            self.assertEqual([record["solid_count"] for record in identities], [1, 1, 1])
            self.assertEqual(cad["step_material_identity_scope"],
                             "geometry-density-bound-ordered-sidecar-v1")
            manifest_path = Path(cad["material_identity_manifest_path"])
            self.assertEqual(manifest_path, step.with_suffix(step.suffix + ".identity.json"))
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            self.assertEqual(manifest["schema"], "PLANAR_STEP_MATERIAL_IDENTITY_V1")
            self.assertEqual(cad["material_identity_manifest_sha256"],
                             hashlib.sha256(manifest_path.read_bytes()).hexdigest())
            self.assertEqual(manifest["serialized_geometry_sha256"],
                             cad["serialized_geometry_sha256"])
            self.assertEqual(manifest["files"]["step"]["sha256"],
                             hashlib.sha256(step.read_bytes()).hexdigest())
            self.assertEqual(manifest["files"]["brep"]["sha256"],
                             hashlib.sha256(brep.read_bytes()).hexdigest())
            self.assertEqual([record["entity_id"] for record in manifest["entities"]],
                             entity_ids)
            expected = np.asarray(result["extruded_complete_material_raw_moments"])
            for stage in ("brep", "step"):
                observed = np.asarray(cad["raw_mass_moments"][stage])
                scale = max(1.0, float(np.max(np.abs(expected))))
                self.assertLessEqual(float(np.max(np.abs(observed-expected)))/scale, 1e-6)
            np.testing.assert_allclose(cad["raw_mass_moments"]["brep"],
                                       cad["raw_mass_moments"]["step"],
                                       rtol=5e-8, atol=1e-8)
            self.assertTrue(brep.is_file() and brep.stat().st_size > 0)
            self.assertTrue(step.is_file() and step.stat().st_size > 0)

            wrong_density = dict(result)
            wrong_entities = [dict(entity) for entity in result["extruded_material_entities"]]
            wrong_entities[1]["volume_density"] = 250.0
            wrong_density["extruded_material_entities"] = tuple(wrong_entities)
            bad_root = root / "wrong-density"
            bad_root.mkdir()
            with self.assertRaisesRegex(RuntimeError, "entity mass moments differ"):
                export(wrong_density, occt_executable=PROBE,
                       brep_path=bad_root / "wrong.brep", step_path=bad_root / "wrong.step")

    def test_rotated_elliptic_cylinder_quadric_clips_exactly_through_step(self):
        """The elliptic-cylinder caps must be outside and preserve clipped moments."""
        self.assertTrue(PROBE.is_file(), "fixed-OCCT roundtrip probe is not built")
        axis = np.asarray((1., 2., 3.)) / math.sqrt(14.)
        first = np.cross(axis, np.asarray((0., 0., 1.)))
        first /= np.linalg.norm(first)
        second = np.cross(axis, first)
        axes = np.column_stack((first, second, axis))
        center_local = np.asarray((0.2, -0.3, 0.1))
        center = 1.7 * center_local
        half_extents = np.asarray((2., 2., 1.5))
        scale = 1.7
        canonical_axis_point = center - axis * np.dot(center, axis)
        radii = scale * np.asarray((0.5, 2. / 3.))
        quadratic = 4. * np.outer(first, first) + 2.25 * np.outer(second, second)
        linear = -2. * quadratic @ center
        coefficients = np.asarray((center @ quadratic @ center - scale**2,
                                   *linear, quadratic[0, 0], quadratic[1, 1],
                                   quadratic[2, 2], quadratic[1, 2],
                                   quadratic[0, 2], quadratic[0, 1]))
        coefficients /= np.linalg.norm(coefficients)

        box_extents = scale * half_extents
        box_volume = 8. * np.prod(box_extents)
        box_covariance = axes @ np.diag(box_extents**2 / 3.) @ axes.T
        box_raw = self._raw_from_center_covariance(box_volume, center, box_covariance)
        cylinder_length = 2. * box_extents[2]
        cylinder_volume = math.pi * np.prod(radii) * cylinder_length
        cylinder_covariance = axes @ np.diag((radii[0]**2/4., radii[1]**2/4.,
                                              cylinder_length**2/12.)) @ axes.T
        cylinder_raw = self._raw_from_center_covariance(
            cylinder_volume, center, cylinder_covariance)
        expected = box_raw - cylinder_raw

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            obj, brep, step = root / "domain.obj", root / "cylinder.brep", root / "cylinder.step"
            self._write_oriented_box_obj(obj, center_local, axes, half_extents)
            completed = subprocess.run(
                [str(PROBE), "--quadric-cylinder", str(obj), str(scale),
                 *(format(x, ".17g") for x in coefficients), str(brep), str(step)],
                check=False, capture_output=True, text=True)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            result = json.loads(completed.stdout)
            self.assertEqual(result["classification"], "SYNTHETIC_ANALYTIC_QUADRIC_CAD_METHOD")
            self.assertEqual(result["status"], "ANALYTIC_ELLIPTIC_CYLINDER_QUADRIC_STEP_BREP_VALIDATED")
            np.testing.assert_allclose(result["cylinder_center"], canonical_axis_point,
                                       rtol=0., atol=2e-10)
            np.testing.assert_allclose(result["cylinder_axis"], axis, rtol=0., atol=2e-10)
            np.testing.assert_allclose(np.sort(result["cylinder_radii"]), np.sort(radii),
                                       rtol=0., atol=2e-10)
            self.assertGreater(result["carrier_height"], cylinder_length)
            for stage in ("brep", "step"):
                np.testing.assert_allclose(result["raw_moments"][stage], expected,
                                           rtol=3e-8, atol=3e-8)

    def test_cylinder_cad_adapter_rejects_paraboloid_axis_linear_term(self):
        self.assertTrue(PROBE.is_file(), "fixed-OCCT roundtrip probe is not built")
        axis = np.asarray((1., 2., 3.)) / math.sqrt(14.)
        first = np.cross(axis, np.asarray((0., 0., 1.)))
        first /= np.linalg.norm(first)
        second = np.cross(axis, first)
        quadratic = 4. * np.outer(first, first) + 2.25 * np.outer(second, second)
        center = np.asarray((0.2, -0.3, 0.1))
        linear = -2. * quadratic @ center + 0.25 * axis
        coefficients = np.asarray((center @ quadratic @ center - 1., *linear,
                                   quadratic[0, 0], quadratic[1, 1], quadratic[2, 2],
                                   quadratic[1, 2], quadratic[0, 2], quadratic[0, 1]))
        coefficients /= np.linalg.norm(coefficients)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            obj, brep, step = root / "domain.obj", root / "invalid.brep", root / "invalid.step"
            self._write_box_obj(obj)
            completed = subprocess.run(
                [str(PROBE), "--quadric-cylinder", str(obj), "1.0",
                 *(format(x, ".17g") for x in coefficients), str(brep), str(step)],
                check=False, capture_output=True, text=True)
            self.assertNotEqual(completed.returncode, 0)
            self.assertIn("linear component along its cylinder axis", completed.stderr)
            self.assertFalse(brep.exists())
            self.assertFalse(step.exists())

    def test_rotated_elliptic_paraboloid_boolean_matches_closed_form_segment_moments(self):
        """A paraboloid carrier must clip to the finite analytic segment in the domain."""
        self.assertTrue(PROBE.is_file(), "fixed-OCCT roundtrip probe is not built")
        axis = np.asarray((1., 2., 3.)) / math.sqrt(14.)
        first = np.cross(axis, np.asarray((0., 0., 1.)))
        first /= np.linalg.norm(first)
        second = np.cross(axis, first)
        axes = np.column_stack((first, second, axis))
        half_extents = np.asarray((2., 2., 1.2))
        lambda_x, lambda_y, axial = 4., 2.25, 1.
        quadratic = lambda_x*np.outer(first, first) + lambda_y*np.outer(second, second)
        linear = -axial * axis
        coefficients = np.asarray((0., *linear, quadratic[0, 0], quadratic[1, 1],
                                   quadratic[2, 2], quadratic[1, 2],
                                   quadratic[0, 2], quadratic[0, 1]))
        coefficients /= np.linalg.norm(coefficients)

        height = 1.2
        scale_x = math.sqrt(axial / lambda_x)
        scale_y = math.sqrt(axial / lambda_y)
        volume = scale_x * scale_y * math.pi * height**2 / 2.
        first_local = np.asarray((0., 0., scale_x*scale_y*math.pi*height**3/3.))
        second_local = np.diag((
            scale_x**3 * scale_y * math.pi * height**3 / 12.,
            scale_x * scale_y**3 * math.pi * height**3 / 12.,
            scale_x * scale_y * math.pi * height**4 / 4.))
        first_world = axes @ first_local
        second_world = axes @ second_local @ axes.T
        cavity_raw = np.asarray((volume, *first_world,
                                 second_world[0, 0], second_world[1, 1], second_world[2, 2],
                                 second_world[1, 2], second_world[0, 2], second_world[0, 1]))
        outer_volume = 8. * np.prod(half_extents)
        outer_covariance = axes @ np.diag(half_extents**2/3.) @ axes.T
        outer_raw = self._raw_from_center_covariance(
            outer_volume, np.zeros(3), outer_covariance)
        expected = outer_raw - cavity_raw

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            obj, brep, step = root / "domain.obj", root / "paraboloid.brep", root / "paraboloid.step"
            self._write_oriented_box_obj(obj, np.zeros(3), axes, half_extents)
            completed = subprocess.run(
                [str(PROBE), "--quadric-paraboloid", str(obj), "1.0",
                 *(format(x, ".17g") for x in coefficients), str(brep), str(step)],
                check=False, capture_output=True, text=True)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            result = json.loads(completed.stdout)
            self.assertEqual(result["classification"], "SYNTHETIC_ANALYTIC_QUADRIC_CAD_METHOD")
            self.assertEqual(result["status"], "ANALYTIC_ELLIPTIC_PARABOLOID_QUADRIC_STEP_BREP_VALIDATED")
            np.testing.assert_allclose(result["paraboloid_axis"], axis, rtol=0., atol=2e-10)
            for stage in ("brep", "step"):
                np.testing.assert_allclose(result["raw_moments"][stage], expected,
                                           rtol=5e-8, atol=5e-8)

    def test_paraboloid_with_positive_axial_linear_term_reverses_carrier_axis(self):
        self.assertTrue(PROBE.is_file(), "fixed-OCCT roundtrip probe is not built")
        axis = np.asarray((1., 2., 3.)) / math.sqrt(14.)
        first = np.cross(axis, np.asarray((0., 0., 1.)))
        first /= np.linalg.norm(first)
        second = np.cross(axis, first)
        axes = np.column_stack((first, second, axis))
        half_extents = np.asarray((2., 2., 1.2))
        lambda_x, lambda_y, height = 4., 2.25, 1.2
        quadratic = lambda_x*np.outer(first, first) + lambda_y*np.outer(second, second)
        coefficients = np.asarray((0., *axis, quadratic[0, 0], quadratic[1, 1],
                                   quadratic[2, 2], quadratic[1, 2],
                                   quadratic[0, 2], quadratic[0, 1]))
        coefficients /= np.linalg.norm(coefficients)
        scale_x = math.sqrt(1. / lambda_x)
        scale_y = math.sqrt(1. / lambda_y)
        volume = scale_x*scale_y*math.pi*height**2/2.
        first_local = np.asarray((0., 0., -scale_x*scale_y*math.pi*height**3/3.))
        second_local = np.diag((scale_x**3*scale_y*math.pi*height**3/12.,
                                scale_x*scale_y**3*math.pi*height**3/12.,
                                scale_x*scale_y*math.pi*height**4/4.))
        cavity_raw = np.asarray((volume, *(axes @ first_local),
                                 *((axes @ second_local @ axes.T)[i, i] for i in range(3)),
                                 (axes @ second_local @ axes.T)[1, 2],
                                 (axes @ second_local @ axes.T)[0, 2],
                                 (axes @ second_local @ axes.T)[0, 1]))
        outer_volume = 8. * np.prod(half_extents)
        outer_covariance = axes @ np.diag(half_extents**2/3.) @ axes.T
        expected = self._raw_from_center_covariance(
            outer_volume, np.zeros(3), outer_covariance) - cavity_raw

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            obj, brep, step = root / "domain.obj", root / "reverse-paraboloid.brep", root / "reverse-paraboloid.step"
            self._write_oriented_box_obj(obj, np.zeros(3), axes, half_extents)
            completed = subprocess.run(
                [str(PROBE), "--quadric-paraboloid", str(obj), "1.0",
                 *(format(x, ".17g") for x in coefficients), str(brep), str(step)],
                check=False, capture_output=True, text=True)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            result = json.loads(completed.stdout)
            np.testing.assert_allclose(result["paraboloid_axis"], -axis,
                                       rtol=0., atol=2e-10)
            for stage in ("brep", "step"):
                np.testing.assert_allclose(result["raw_moments"][stage], expected,
                                           rtol=5e-8, atol=5e-8)

    def test_rotated_hyperbolic_paraboloid_boolean_matches_closed_form_moments(self):
        """A saddle carrier must clip the exact q<=0 region and survive STEP."""
        self.assertTrue(PROBE.is_file(), "fixed-OCCT roundtrip probe is not built")
        axis = np.asarray((1., 2., 3.)) / math.sqrt(14.)
        first = np.cross(axis, np.asarray((0., 0., 1.)))
        first /= np.linalg.norm(first)
        second = np.cross(axis, first)
        axes = np.column_stack((first, second, axis))
        half_xy, half_z = 1.0, 2.0
        quadratic = np.outer(first, first) - np.outer(second, second)
        vertex = axes @ np.asarray((0.23, -0.17, 0.11))
        linear = axis - 2.0 * quadratic @ vertex
        constant = float(vertex @ quadratic @ vertex - axis @ vertex)
        coefficients = np.asarray((constant, *linear, quadratic[0, 0], quadratic[1, 1],
                                   quadratic[2, 2], quadratic[1, 2],
                                   quadratic[0, 2], quadratic[0, 1]))
        coefficients /= np.linalg.norm(coefficients)

        cavity_volume = 4.0 * half_xy**2 * half_z
        cavity_first_local = np.asarray((
            0.0, 0.0, 16.0*half_xy**6/45.0 - 2.0*half_xy**2*half_z**2))
        cavity_second_local = np.diag((
            -16.0*half_xy**6/45.0 + 4.0*half_z*half_xy**4/3.0,
            16.0*half_xy**6/45.0 + 4.0*half_z*half_xy**4/3.0,
            4.0*half_xy**2*half_z**3/3.0))
        cavity_local_first_world = axes @ cavity_first_local
        cavity_first = cavity_volume * vertex + cavity_local_first_world
        cavity_second = (
            axes @ cavity_second_local @ axes.T +
            cavity_volume * np.outer(vertex, vertex) +
            np.outer(vertex, cavity_local_first_world) +
            np.outer(cavity_local_first_world, vertex))
        cavity_raw = np.asarray((cavity_volume, *cavity_first,
                                 cavity_second[0, 0], cavity_second[1, 1],
                                 cavity_second[2, 2], cavity_second[1, 2],
                                 cavity_second[0, 2], cavity_second[0, 1]))
        outer_volume = 8.0 * half_xy**2 * half_z
        outer_covariance = axes @ np.diag((half_xy**2/3.0, half_xy**2/3.0,
                                           half_z**2/3.0)) @ axes.T
        outer_raw = self._raw_from_center_covariance(
            outer_volume, vertex, outer_covariance)
        expected = outer_raw - cavity_raw

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            obj = root / "domain.obj"
            brep = root / "saddle.brep"
            step = root / "saddle.step"
            self._write_oriented_box_obj(
                obj, vertex, axes, (half_xy, half_xy, half_z))
            completed = subprocess.run(
                [str(PROBE), "--quadric-hyperbolic-paraboloid", str(obj), "1.0",
                 *(format(x, ".17g") for x in coefficients), str(brep), str(step)],
                check=False, capture_output=True, text=True)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            result = json.loads(completed.stdout)
            self.assertEqual(result["classification"],
                             "SYNTHETIC_ANALYTIC_QUADRIC_CAD_METHOD")
            self.assertEqual(result["status"],
                             "ANALYTIC_HYPERBOLIC_PARABOLOID_QUADRIC_STEP_BREP_VALIDATED")
            self.assertEqual(result["surface_representation"],
                             "exact_tensor_product_quadratic_bezier")
            np.testing.assert_allclose(result["vertex"], vertex, rtol=0., atol=2e-10)
            np.testing.assert_allclose(result["axis"], axis, rtol=0., atol=2e-10)
            for stage in ("brep", "step"):
                np.testing.assert_allclose(result["raw_moments"][stage], expected,
                                           rtol=8e-8, atol=8e-8)
            repeat = subprocess.run(
                [str(PROBE), "--quadric-hyperbolic-paraboloid", str(obj), "1.0",
                 *(format(x, ".17g") for x in coefficients),
                 str(root / "saddle-repeat.brep"), str(root / "saddle-repeat.step")],
                check=False, capture_output=True, text=True)
            self.assertEqual(repeat.returncode, 0, repeat.stderr)
            self.assertEqual(repeat.stdout, completed.stdout)

            # Reversing only the axial linear coefficient switches q<=0 to
            # the region above the saddle. The canonical frame must reverse
            # its axis rather than reusing the previous carrier orientation.
            negative_linear = -axis - 2.0 * quadratic @ vertex
            negative_constant = float(vertex @ quadratic @ vertex + axis @ vertex)
            negative_coefficients = np.asarray((negative_constant, *negative_linear,
                                                quadratic[0, 0], quadratic[1, 1],
                                                quadratic[2, 2], quadratic[1, 2],
                                                quadratic[0, 2], quadratic[0, 1]))
            negative_coefficients /= np.linalg.norm(negative_coefficients)
            upper_first_local = cavity_first_local.copy()
            upper_first_local[2] *= -1.0
            # In the caller's frame the sign-reversed axial coefficient gives
            # w=u^2-v^2; reflecting the canonical axis changes which physical
            # side is retained but does not exchange the principal directions.
            upper_second_local = cavity_second_local.copy()
            upper_local_first_world = axes @ upper_first_local
            upper_first = cavity_volume * vertex + upper_local_first_world
            upper_second = (
                axes @ upper_second_local @ axes.T +
                cavity_volume * np.outer(vertex, vertex) +
                np.outer(vertex, upper_local_first_world) +
                np.outer(upper_local_first_world, vertex))
            upper_raw = np.asarray((cavity_volume, *upper_first,
                                    upper_second[0, 0], upper_second[1, 1],
                                    upper_second[2, 2], upper_second[1, 2],
                                    upper_second[0, 2], upper_second[0, 1]))
            upper_expected = outer_raw - upper_raw
            upper_brep, upper_step = root / "saddle-upper.brep", root / "saddle-upper.step"
            upper_completed = subprocess.run(
                [str(PROBE), "--quadric-hyperbolic-paraboloid", str(obj), "1.0",
                 *(format(x, ".17g") for x in negative_coefficients),
                 str(upper_brep), str(upper_step)],
                check=False, capture_output=True, text=True)
            self.assertEqual(upper_completed.returncode, 0, upper_completed.stderr)
            upper_result = json.loads(upper_completed.stdout)
            np.testing.assert_allclose(upper_result["axis"], -axis, rtol=0., atol=2e-10)
            for stage in ("brep", "step"):
                np.testing.assert_allclose(upper_result["raw_moments"][stage],
                                           upper_expected, rtol=8e-8, atol=8e-8)

    def test_hyperbolic_paraboloid_cad_rejects_rank_two_cylinder_without_axial_linear_term(self):
        self.assertTrue(PROBE.is_file(), "fixed-OCCT roundtrip probe is not built")
        quadratic = np.diag((1., -1., 0.))
        vertex = np.asarray((0.2, -0.3, 0.4))
        coefficients = np.asarray((vertex @ quadratic @ vertex,
                                   *(-2.0 * quadratic @ vertex),
                                   quadratic[0, 0], quadratic[1, 1], quadratic[2, 2],
                                   0., 0., 0.))
        coefficients /= np.linalg.norm(coefficients)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            obj, brep, step = root / "domain.obj", root / "invalid.brep", root / "invalid.step"
            self._write_box_obj(obj)
            completed = subprocess.run(
                [str(PROBE), "--quadric-hyperbolic-paraboloid", str(obj), "1.0",
                 *(format(x, ".17g") for x in coefficients), str(brep), str(step)],
                check=False, capture_output=True, text=True)
            self.assertNotEqual(completed.returncode, 0)
            self.assertIn("nonzero linear component along its null axis", completed.stderr)
            self.assertFalse(brep.exists())
            self.assertFalse(step.exists())

    def test_hyperbolic_paraboloid_cad_refuses_to_overwrite_existing_evidence(self):
        self.assertTrue(PROBE.is_file(), "fixed-OCCT roundtrip probe is not built")
        coefficients = np.asarray((0., 0., 0., 1., 1., -1., 0., 0., 0., 0.))
        coefficients /= np.linalg.norm(coefficients)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            obj, brep, step = root / "domain.obj", root / "existing.brep", root / "new.step"
            self._write_box_obj(obj)
            brep.write_bytes(b"preserve-existing-cad-evidence")
            completed = subprocess.run(
                [str(PROBE), "--quadric-hyperbolic-paraboloid", str(obj), "1.0",
                 *(format(x, ".17g") for x in coefficients), str(brep), str(step)],
                check=False, capture_output=True, text=True)
            self.assertNotEqual(completed.returncode, 0)
            self.assertIn("refusing to overwrite", completed.stderr)
            self.assertEqual(brep.read_bytes(), b"preserve-existing-cad-evidence")
            self.assertFalse(step.exists())

    def test_rotated_one_sheet_hyperboloid_quadric_clips_to_analytic_domain_moments(self):
        """A hyperboloid carrier must preserve the finite one-sheet sublevel volume."""
        self.assertTrue(PROBE.is_file(), "fixed-OCCT roundtrip probe is not built")
        axis = np.asarray((1., 2., 3.)) / math.sqrt(14.)
        first = np.cross(axis, np.asarray((0., 0., 1.)))
        first /= np.linalg.norm(first)
        second = np.cross(axis, first)
        axes = np.column_stack((first, second, axis))
        radius_x, radius_y, axial_scale, half_height = 0.7, 0.55, 0.8, 0.7
        quadratic = (np.outer(first, first)/radius_x**2 +
                     np.outer(second, second)/radius_y**2 -
                     np.outer(axis, axis)/axial_scale**2)
        coefficients = np.asarray((-1., 0., 0., 0., quadratic[0, 0],
                                   quadratic[1, 1], quadratic[2, 2],
                                   quadratic[1, 2], quadratic[0, 2], quadratic[0, 1]))
        coefficients /= np.linalg.norm(coefficients)

        axial_integral_0 = 2.*half_height + 2.*half_height**3/(3.*axial_scale**2)
        axial_integral_2 = (2.*half_height + 4.*half_height**3/(3.*axial_scale**2) +
                            2.*half_height**5/(5.*axial_scale**4))
        axial_integral_z2 = (2.*half_height**3/3. +
                             2.*half_height**5/(5.*axial_scale**2))
        cavity_volume = math.pi*radius_x*radius_y*axial_integral_0
        cavity_local_second = np.diag((
            math.pi*radius_x**3*radius_y*axial_integral_2/4.,
            math.pi*radius_x*radius_y**3*axial_integral_2/4.,
            math.pi*radius_x*radius_y*axial_integral_z2))
        cavity_second = axes @ cavity_local_second @ axes.T
        cavity_raw = np.asarray((cavity_volume, 0., 0., 0.,
                                 cavity_second[0, 0], cavity_second[1, 1],
                                 cavity_second[2, 2], cavity_second[1, 2],
                                 cavity_second[0, 2], cavity_second[0, 1]))
        half_extents = np.asarray((2., 2., half_height))
        outer_volume = 8.*np.prod(half_extents)
        outer_second = outer_volume * axes @ np.diag(half_extents**2/3.) @ axes.T
        outer_raw = np.asarray((outer_volume, 0., 0., 0., outer_second[0, 0],
                                outer_second[1, 1], outer_second[2, 2],
                                outer_second[1, 2], outer_second[0, 2], outer_second[0, 1]))
        expected = outer_raw - cavity_raw

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            obj, brep, step = root / "domain.obj", root / "hyperboloid.brep", root / "hyperboloid.step"
            self._write_oriented_box_obj(obj, np.zeros(3), axes, half_extents)
            completed = subprocess.run(
                [str(PROBE), "--quadric-hyperboloid-one-sheet", str(obj), "1.0",
                 *(format(x, ".17g") for x in coefficients), str(brep), str(step)],
                check=False, capture_output=True, text=True)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            result = json.loads(completed.stdout)
            self.assertEqual(result["classification"], "SYNTHETIC_ANALYTIC_QUADRIC_CAD_METHOD")
            self.assertEqual(result["status"], "ANALYTIC_ONE_SHEET_HYPERBOLOID_QUADRIC_STEP_BREP_VALIDATED")
            self.assertAlmostEqual(abs(np.dot(result["hyperboloid_axis"], axis)), 1.0, places=10)
            for stage in ("brep", "step"):
                np.testing.assert_allclose(result["raw_moments"][stage], expected,
                                           rtol=5e-8, atol=5e-8)

    def test_one_sheet_hyperboloid_path_rejects_other_full_rank_signatures(self):
        """A superficially similar quadric must not be exported with wrong topology."""
        self.assertTrue(PROBE.is_file(), "fixed-OCCT roundtrip probe is not built")
        coefficients = np.asarray((-1., 0., 0., 0., 1., -1., -1., 0., 0., 0.))
        coefficients /= np.linalg.norm(coefficients)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            obj, brep, step = root / "domain.obj", root / "unsupported.obj.brep", root / "unsupported.obj.step"
            self._write_box_obj(obj)
            completed = subprocess.run(
                [str(PROBE), "--quadric-hyperboloid-one-sheet", str(obj), "1.0",
                 *(format(x, ".17g") for x in coefficients), str(brep), str(step)],
                check=False, capture_output=True, text=True)
            self.assertNotEqual(completed.returncode, 0)
            self.assertIn("one-sheet hyperboloid signature (+,+,-)", completed.stderr)
            self.assertFalse(brep.exists())
            self.assertFalse(step.exists())

    def test_near_cone_positive_level_preserves_analytic_material_moments(self):
        """A strictly positive small neck is not a zero-level cone.

        The old absolute center-level cutoff rejects this normalized, full
        rank input before building geometry. Integration expectations below
        come from disks r(z)^2=level+z^2, not from the CAD implementation.
        """
        level, h = 1e-14, .5
        coefficients = np.asarray((-level, 0., 0., 0., 1., 1., -1., 0., 0., 0.))
        coefficients /= np.linalg.norm(coefficients)
        cavity = np.asarray((2*math.pi*(level*h+h**3/3), 0., 0., 0.,
            math.pi/2*(level**2*h+2*level*h**3/3+h**5/5),
            math.pi/2*(level**2*h+2*level*h**3/3+h**5/5),
            2*math.pi*(level*h**3/3+h**5/5), 0., 0., 0.))
        outer = np.asarray((16.,0.,0.,0.,64./3,64./3,4./3,0.,0.,0.))
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            obj,brep,step=root/'domain.obj',root/'material.brep',root/'material.step'
            self._write_oriented_box_obj(obj,np.zeros(3),np.eye(3),(2.,2.,h))
            completed=subprocess.run([str(PROBE),'--quadric-hyperboloid-one-sheet',
                str(obj),'1.0',*(format(x,'.17g') for x in coefficients),str(brep),str(step)],
                capture_output=True,text=True,timeout=60)
            self.assertEqual(completed.returncode,0,completed.stderr)
            result=json.loads(completed.stdout)
            for stage in ('brep','step'):
                np.testing.assert_allclose(result['raw_moments'][stage],outer-cavity,
                                           rtol=5e-8,atol=5e-8)


    def test_invalid_nearcone_retains_subshape_witness_without_approving_cad(self):
        # Missing structured rejection evidence must fail independently of
        # the still-red positive CAD acceptance test above.
        q = np.asarray((-1e-14,0.,0.,0.,1.,1.,-1.,0.,0.,0.))
        q /= np.linalg.norm(q)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            obj, brep, step = root/'domain.obj', root/'out.brep', root/'out.step'
            self._write_oriented_box_obj(obj,np.zeros(3),np.eye(3),(2.,2.,.5))
            witness = root/'rejected'
            env = dict(os.environ, MAKE_IT_STAND_OCCT_FAILURE_DIR=str(witness))
            command = [str(PROBE),'--quadric-hyperboloid-one-sheet',str(obj),'1.0',
                       *(format(v,'.17g') for v in q),str(brep),str(step)]
            result = subprocess.run(command,env=env,capture_output=True,text=True,timeout=60)
            self.assertNotEqual(result.returncode,0)
            self.assertFalse(brep.exists())
            self.assertFalse(step.exists())
            stages = [json.loads(line[len('OCCT_CARRIER_STAGE '):])
                      for line in result.stderr.splitlines()
                      if line.startswith('OCCT_CARRIER_STAGE ')]
            self.assertEqual([item['stage'] for item in stages],
                             ['meridian','revolution','affine'])
            for item in stages:
                self.assertGreater(item['subshape_count'],0)
                self.assertEqual(sum(item['type_counts']),item['subshape_count'])
            report_path = witness/'topology.json'
            self.assertTrue(report_path.is_file(),result.stderr)
            report = json.loads(report_path.read_text())
            self.assertEqual(report['schema'],'occt-rejected-topology-v1')
            self.assertIs(report['valid'],False)
            self.assertEqual(report['stage'],'one-sheet hyperboloid Boolean result')
            self.assertTrue(report['defects'])
            self.assertGreater(report['subshape_count'],0)
            for defect in report['defects']:
                self.assertGreater(defect['subshape_index'],0)
                self.assertLessEqual(defect['subshape_index'],report['subshape_count'])
                self.assertNotEqual(defect['status_code'],0)
                self.assertTrue(defect['status'])
                self.assertGreaterEqual(defect['context_index'],0)
                self.assertLessEqual(defect['context_index'],report['subshape_count'])
            self.assertGreater((witness/'invalid.brep').stat().st_size,0)
            before = {p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                      for p in witness.iterdir() if p.is_file()}
            repeated = subprocess.run(command,env=env,capture_output=True,text=True,timeout=60)
            self.assertNotEqual(repeated.returncode,0)
            self.assertIn('refusing to overwrite',repeated.stderr)
            self.assertEqual(before,{p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                             for p in witness.iterdir() if p.is_file()})

    def test_typed_invalid_material_retains_entity_topology_witness(self):
        q=np.asarray((-1e-14,0.,0.,0.,1.,1.,-1.,0.,0.,0.))
        q/=np.linalg.norm(q)
        # Identity tokens are not restricted to hexadecimal by this protocol.
        # The rejection JSON must escape them instead of becoming unreadable.
        identity='entity\x01"\\suffix'
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            obj=root/'domain.obj'
            self._write_oriented_box_obj(obj,np.zeros(3),np.eye(3),(2.,2.,.5))
            protocol=root/'input.txt'
            protocol.write_text('FASTER_3D_QUADRIC_MATERIAL_V2\n0\n'+
                ' '.join(format(v,'.17g') for v in q)+'\nONE_SHEET_HYPERBOLOID\n1\n'+
                f'ENTITY design {identity} 1 {json.dumps(str(obj))}\n')
            witness=root/'rejected'
            run=subprocess.run([str(PROBE),'--faster-3d-quadric-material',str(protocol),
                str(root/'out.brep'),str(root/'out.step')],capture_output=True,text=True,timeout=60,
                env=dict(os.environ,MAKE_IT_STAND_OCCT_FAILURE_DIR=str(witness)))
            self.assertNotEqual(run.returncode,0)
            self.assertTrue((witness/'topology.json').is_file(),run.stderr)
            report=json.loads((witness/'topology.json').read_text())
            self.assertEqual(report['stage'],'Faster material entity '+identity)
            self.assertTrue(report['defects'])
            self.assertIs(report['valid'],False)
            self.assertFalse((root/'out.brep').exists())
            self.assertFalse((root/'out.step').exists())

    def test_saved_nearcone_meridian_does_not_silently_drop_curved_edge(self):
        # Saved polynomial, not a new optimized input or an approximate cone.
        q = (-4.8194380240747099e-38,2.0338229586518573e-19,
             -1.3585749154241787e-18,7.0991444041997789e-20,
             -.82631939729131509,.39824375457455563,.39824385193305739,
             4.153467832522168e-20,-2.1685150389145536e-18,1.6741243732157981e-17)
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            obj=root/'domain.obj'
            self._write_oriented_box_obj(obj,np.zeros(3),np.eye(3),(.25,.25,.25))
            run=subprocess.run([str(PROBE),'--quadric-hyperboloid-one-sheet',str(obj),'1',
                *(format(v,'.17g') for v in q),str(root/'out.brep'),str(root/'out.step')],
                env=dict(os.environ,MAKE_IT_STAND_OCCT_FAILURE_DIR=str(root/'witness')),
                capture_output=True,text=True,timeout=60)
            stages=[json.loads(line[len('OCCT_CARRIER_STAGE '):]) for line in run.stderr.splitlines()
                    if line.startswith('OCCT_CARRIER_STAGE ')]
            self.assertTrue(stages,run.stderr)
            meridian=stages[0]
            self.assertEqual(meridian['stage'],'meridian')
            self.assertEqual(meridian['type_counts'][6],4,
                             'axis + two caps + hyperbola must all remain in the wire')
            self.assertIs(meridian['valid'],True)
            # This assertion covers the input face only. Full CAD has its
            # own unwaived acceptance test and may still fail downstream.

    def test_saved_current_terminal_cannot_approve_wrong_cad_material(self):
        # Frozen current terminal's design entity (fixed shell not needed here),
        # not a fitted cone or an easier new problem. Provenance: sprint
        # cherchi-cad/current-terminal-carrier-02/request.json.
        # Analytic carrier accuracy is separate from solver convergence. A
        # conversion rejection is safe, but cannot count as successful CAD.
        q=(-1.1479438025969792e-28,5.703740330169272e-20,
           -1.689149773316858e-20,-3.975533057052151e-20,
           -.826319397291315,.3982437545745472,.3982438519330661,
           1.27852423201439e-21,-5.729176781573e-20,1.0566499509012466e-19)
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);obj=root/'domain.obj'
            self._write_oriented_box_obj(obj,np.zeros(3),np.eye(3),(.25,.25,.25))
            protocol=root/'input.txt';brep=root/'out.brep';step=root/'out.step'
            protocol.write_text('FASTER_3D_QUADRIC_MATERIAL_V2\n0\n'+
                ' '.join(format(v,'.17g') for v in q)+'\nONE_SHEET_HYPERBOLOID\n1\n'+
                f'ENTITY design frozen-current-terminal 1 {json.dumps(str(obj))}\n')
            witness=root/'rejected'
            run=subprocess.run([str(PROBE),'--faster-3d-quadric-material',str(protocol),
                str(brep),str(step)],capture_output=True,text=True,timeout=60,
                env=dict(os.environ,MAKE_IT_STAND_OCCT_FAILURE_DIR=str(witness)))
            if run.returncode==0:
                result=json.loads(run.stdout)
                # This low-level probe constructs the ANALYTIC polynomial,
                # not the solver's discretized boundary. Independent full-q
                # integration: nearcone-diagnostic/independent-quadrature-01/
                # quadrature.json; two tolerances agree within 1e-12 and a
                # sphere control meets 1e-10. The full exporter separately
                # enforces agreement with the actual sealed material.
                self.assertLessEqual(abs(result['entity_mass_moments'][0][0]-
                                         .06640229393839628),1e-6,
                    'valid B-rep topology cannot approve the uncarved full domain')
            else:
                self.assertFalse(brep.exists());self.assertFalse(step.exists())
                self.assertTrue((witness/'topology.json').is_file(),run.stderr)
                report=json.loads((witness/'topology.json').read_text())
                self.assertIs(report['valid'],False)
                self.assertTrue(report['defects'])
                self.assertIn(report['stage'],('one-sheet hyperboloid NURBS conversion',
                    'Faster material entity frozen-current-terminal'))

    def test_rotated_two_sheet_hyperboloid_preserves_both_material_lobes(self):
        """The q>=0 material must be both analytic lobes, not the center band."""
        self.assertTrue(PROBE.is_file(), "fixed-OCCT roundtrip probe is not built")
        axis = np.asarray((1., 2., 3.)) / math.sqrt(14.)
        first = np.cross(axis, np.asarray((0., 0., 1.)))
        first /= np.linalg.norm(first)
        second = np.cross(axis, first)
        axes = np.column_stack((axis, first, second))
        axial_scale, radial_x_scale, radial_y_scale = 0.5, 0.25, 0.2
        quadratic = (np.outer(axis, axis)/axial_scale**2 -
                     np.outer(first, first)/radial_x_scale**2 -
                     np.outer(second, second)/radial_y_scale**2)
        coefficients = np.asarray((-1., 0., 0., 0., quadratic[0, 0],
                                   quadratic[1, 1], quadratic[2, 2],
                                   quadratic[1, 2], quadratic[0, 2], quadratic[0, 1]))
        coefficients /= np.linalg.norm(coefficients)

        normalized_cap = 2.0 / axial_scale
        volume_integral = normalized_cap**3/3. - normalized_cap + 2./3.
        axial_second_integral = (normalized_cap**5/5. - normalized_cap**3/3. + 2./15.)
        radial_second_integral = (normalized_cap**5/5. - 2.*normalized_cap**3/3. +
                                  normalized_cap - 8./15.)
        jacobian = axial_scale*radial_x_scale*radial_y_scale
        material_volume = 2.*math.pi*jacobian*volume_integral
        local_second = np.diag((
            2.*math.pi*jacobian*axial_scale**2*axial_second_integral,
            math.pi/2.*jacobian*radial_x_scale**2*radial_second_integral,
            math.pi/2.*jacobian*radial_y_scale**2*radial_second_integral,
        ))
        world_second = axes @ local_second @ axes.T
        expected = np.asarray((material_volume, 0., 0., 0.,
                               world_second[0, 0], world_second[1, 1],
                               world_second[2, 2], world_second[1, 2],
                               world_second[0, 2], world_second[0, 1]))

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            obj, brep, step = root / "domain.obj", root / "two-sheet.brep", root / "two-sheet.step"
            self._write_oriented_box_obj(obj, np.zeros(3), axes, (2., 2., 2.))
            completed = subprocess.run(
                [str(PROBE), "--quadric-hyperboloid-two-sheet", str(obj), "1.0",
                 *(format(x, ".17g") for x in coefficients), str(brep), str(step)],
                check=False, capture_output=True, text=True)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            result = json.loads(completed.stdout)
            self.assertEqual(result["classification"], "SYNTHETIC_ANALYTIC_QUADRIC_CAD_METHOD")
            self.assertEqual(result["status"], "ANALYTIC_TWO_SHEET_HYPERBOLOID_QUADRIC_STEP_BREP_VALIDATED")
            self.assertEqual(result["material_solid_count"], 2)
            self.assertAlmostEqual(abs(np.dot(result["hyperboloid_axis"], axis)), 1.0, places=10)
            brep_moments = np.asarray(result["raw_moments"]["brep"])
            step_moments = np.asarray(result["raw_moments"]["step"])
            np.testing.assert_allclose(step_moments, brep_moments, rtol=5e-8, atol=1e-8)
            for stage in ("brep", "step"):
                np.testing.assert_allclose(result["raw_moments"][stage], expected,
                                           rtol=1e-6, atol=1e-8)

    def test_two_sheet_hyperboloid_path_rejects_one_sheet_material_signature(self):
        """The two-lobe CAD branch must fail closed on the opposite topology."""
        self.assertTrue(PROBE.is_file(), "fixed-OCCT roundtrip probe is not built")
        # q>=0 is the connected outside of a one-sheet hyperboloid, not two
        # axial lobes. The two-sheet adapter must not relabel it as such.
        coefficients = np.asarray((-1., 0., 0., 0., -1., 1., 1., 0., 0., 0.))
        coefficients /= np.linalg.norm(coefficients)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            obj, brep, step = root / "domain.obj", root / "wrong-signature.brep", root / "wrong-signature.step"
            self._write_box_obj(obj)
            completed = subprocess.run(
                [str(PROBE), "--quadric-hyperboloid-two-sheet", str(obj), "1.0",
                 *(format(x, ".17g") for x in coefficients), str(brep), str(step)],
                check=False, capture_output=True, text=True)
            self.assertNotEqual(completed.returncode, 0)
            self.assertIn("hyperboloid signature", completed.stderr)
            self.assertFalse(brep.exists())
            self.assertFalse(step.exists())

    def test_solver_quadric_sphere_coefficients_drive_exact_occt_boolean_and_step(self) -> None:
        """A sphere-valued solver quadric must remain analytic through STEP."""
        self.assertTrue(PROBE.is_file(), "fixed-OCCT roundtrip probe is not built")
        center = np.asarray((0.25, -0.125, 0.5), dtype=np.float64)
        radius = 0.5
        coefficients = np.asarray((center @ center - radius**2,
                                   -2. * center[0], -2. * center[1], -2. * center[2],
                                   1., 1., 1., 0., 0., 0.), dtype=np.float64)
        coefficients /= np.linalg.norm(coefficients)
        target_axes = np.asarray(((0., -1., 0.), (1., 0., 0.), (0., 0., 1.)))
        target_origin = np.asarray((0.2, -0.1, 0.05))
        coordinate_scale = 2.5
        cad_coefficients = quadric_coefficients_to_cad_frame(
            coefficients, target_axes=target_axes, target_origin=target_origin,
            original_units_per_normalized_unit=coordinate_scale)
        cad_center = coordinate_scale * (target_origin + target_axes @ center)
        cad_radius = coordinate_scale * radius
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            obj, brep, step = root / "domain.obj", root / "sphere-cut.brep", root / "sphere-cut.step"
            self._write_box_obj(obj)
            completed = subprocess.run(
                [str(PROBE), "--quadric-sphere", str(obj), str(coordinate_scale),
                 *(format(x, ".17g") for x in cad_coefficients), str(brep), str(step)],
                check=False, capture_output=True, text=True)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            result = json.loads(completed.stdout)
            self.assertEqual(result["classification"], "SYNTHETIC_ANALYTIC_QUADRIC_CAD_METHOD")
            self.assertEqual(result["status"], "ANALYTIC_SPHERE_QUADRIC_STEP_BREP_VALIDATED")
            self.assertTrue(brep.is_file() and brep.stat().st_size > 0)
            self.assertTrue(step.is_file() and step.stat().st_size > 0)

            sphere_volume = 4. * math.pi * cad_radius**3 / 3.
            side_scale = coordinate_scale
            outer_raw = np.asarray((64. * side_scale**3, 0., 0., 0.,
                                    256./3. * side_scale**5,
                                    256./3. * side_scale**5,
                                    256./3. * side_scale**5, 0., 0., 0.))
            sphere_raw = np.asarray((sphere_volume,
                                     *(sphere_volume * cad_center),
                                     *(sphere_volume * (cad_radius**2 / 5. + cad_center**2)),
                                     sphere_volume * cad_center[1] * cad_center[2],
                                     sphere_volume * cad_center[0] * cad_center[2],
                                     sphere_volume * cad_center[0] * cad_center[1]))
            expected = outer_raw - sphere_raw
            for stage in ("brep", "step"):
                np.testing.assert_allclose(result["raw_moments"][stage], expected,
                                           rtol=2e-9, atol=2e-9)
            np.testing.assert_allclose(result["sphere_center"], cad_center, atol=1e-12, rtol=0.)
            self.assertAlmostEqual(result["sphere_radius"], cad_radius, places=12)

    def test_quadric_cad_rejects_nonspherical_coefficients_instead_of_mesh_approximating(self) -> None:
        self.assertTrue(PROBE.is_file(), "fixed-OCCT roundtrip probe is not built")
        coefficients = np.asarray((-.25, 0., 0., 0., 1., 2., 1., 0., 0., 0.), dtype=float)
        coefficients /= np.linalg.norm(coefficients)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            obj, brep, step = root / "domain.obj", root / "invalid.brep", root / "invalid.step"
            self._write_box_obj(obj)
            completed = subprocess.run(
                [str(PROBE), "--quadric-sphere", str(obj), "1.0",
                 *(format(x, ".17g") for x in coefficients), str(brep), str(step)],
                check=False, capture_output=True, text=True)
            self.assertNotEqual(completed.returncode, 0)
            self.assertIn("is not an isotropic sphere quadric", completed.stderr)
            self.assertFalse(brep.exists())
            self.assertFalse(step.exists())

    def test_solver_ellipsoid_quadric_preserves_cross_terms_through_exact_cad_step(self) -> None:
        self.assertTrue(PROBE.is_file(), "fixed-OCCT roundtrip probe is not built")
        center = np.asarray((0.2, -0.15, 0.1), dtype=np.float64)
        radii = np.asarray((0.6, 0.45, 0.35), dtype=np.float64)
        angle = 0.4
        principal_axes = np.asarray(((math.cos(angle), -math.sin(angle), 0.),
                                     (math.sin(angle), math.cos(angle), 0.),
                                     (0., 0., 1.)))
        quadratic = principal_axes @ np.diag(1. / radii**2) @ principal_axes.T
        linear = -2. * quadratic @ center
        coefficients = np.asarray((center @ quadratic @ center - 1., *linear,
                                   quadratic[0, 0], quadratic[1, 1], quadratic[2, 2],
                                   quadratic[1, 2], quadratic[0, 2], quadratic[0, 1]))
        coefficients /= np.linalg.norm(coefficients)
        target_axes = np.asarray(((0., -1., 0.), (1., 0., 0.), (0., 0., 1.)))
        target_origin = np.asarray((0.2, -0.1, 0.05))
        scale = 2.5
        cad_coefficients = quadric_coefficients_to_cad_frame(
            coefficients, target_axes=target_axes, target_origin=target_origin,
            original_units_per_normalized_unit=scale)
        cad_center = scale * (target_origin + target_axes @ center)
        affine_axes = scale * target_axes @ principal_axes @ np.diag(radii)
        volume = 4. * math.pi * np.prod(radii) * scale**3 / 3.
        central_second = affine_axes @ affine_axes.T / 5.
        expected_sphere_raw = np.asarray((
            volume, *(volume * cad_center),
            *(volume * (np.diag(central_second) + cad_center**2)),
            volume * (central_second[1, 2] + cad_center[1] * cad_center[2]),
            volume * (central_second[0, 2] + cad_center[0] * cad_center[2]),
            volume * (central_second[0, 1] + cad_center[0] * cad_center[1])))
        outer_raw = np.asarray((64. * scale**3, 0., 0., 0.,
                                256./3. * scale**5, 256./3. * scale**5,
                                256./3. * scale**5, 0., 0., 0.))
        expected = outer_raw - expected_sphere_raw

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            obj, brep, step = root / "domain.obj", root / "ellipsoid.brep", root / "ellipsoid.step"
            self._write_box_obj(obj)
            completed = subprocess.run(
                [str(PROBE), "--quadric-ellipsoid", str(obj), str(scale),
                 *(format(x, ".17g") for x in cad_coefficients), str(brep), str(step)],
                check=False, capture_output=True, text=True)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            result = json.loads(completed.stdout)
            self.assertEqual(result["classification"], "SYNTHETIC_ANALYTIC_QUADRIC_CAD_METHOD")
            self.assertEqual(result["status"], "ANALYTIC_ELLIPSOID_QUADRIC_STEP_BREP_VALIDATED")
            np.testing.assert_allclose(result["ellipsoid_center"], cad_center, rtol=0., atol=1e-11)
            self.assertAlmostEqual(result["ellipsoid_volume"], volume, places=10)
            self.assertLessEqual(abs(result["emitted_ellipsoid_volume"] - volume),
                                 volume * 3e-8)
            self.assertTrue(brep.is_file() and step.is_file())
            for stage in ("brep", "step"):
                np.testing.assert_allclose(result["raw_moments"][stage], expected,
                                           rtol=3e-8, atol=3e-8)

    def test_box_minus_sphere_brep_and_step_roundtrip_match_analytic_oracle(self) -> None:
        """Changing a CAD export, import, or moment check must fail this test."""
        self.assertTrue(PROBE.is_file(), "fixed-OCCT roundtrip probe is not built")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            brep, step = root / "material.brep", root / "material.step"
            completed = subprocess.run(
                [str(PROBE), str(brep), str(step)], check=False,
                capture_output=True, text=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            result = json.loads(completed.stdout)
            self.assertEqual(result["classification"], "SYNTHETIC_REGRESSION")
            self.assertEqual(result["status"], "SYNTHETIC_OCCT_STEP_BREP_ROUNDTRIP_VALIDATED")
            self.assertTrue(brep.is_file() and brep.stat().st_size > 0)
            self.assertTrue(step.is_file() and step.stat().st_size > 0)
            expected_volume = 8.0 - math.pi / 6.0
            self.assertAlmostEqual(result["analytic_volume"], expected_volume, places=14)
            for key in ("original", "brep", "step"):
                self.assertAlmostEqual(result["moments"][key]["volume"], expected_volume, places=10)
                self.assertEqual(result["moments"][key]["center"], [1.0, 1.0, 1.0])

    def test_solver_output_mesh_roundtrips_through_pinned_occt_with_all_moments(self) -> None:
        """The CAD bridge must consume a real native solver boundary, not a canned primitive."""
        self.assertTrue(PROBE.is_file(), "fixed-OCCT roundtrip probe is not built")

        def box(lower, upper):
            x0, y0, z0 = lower
            x1, y1, z1 = upper
            return (
                ((x0,y0,z0),(x1,y0,z0),(x1,y1,z0),(x0,y1,z0),
                 (x0,y0,z1),(x1,y0,z1),(x1,y1,z1),(x0,y1,z1)),
                ((0,2,1),(0,3,2),(4,5,6),(4,6,7),(0,1,5),(0,5,4),
                 (3,7,6),(3,6,2),(0,4,7),(0,7,3),(1,2,6),(1,6,5)),
            )

        def combine(*meshes):
            vertices, faces, offset = [], [], 0
            for mesh_vertices, mesh_faces in meshes:
                vertices.extend(mesh_vertices)
                faces.extend(tuple(index + offset for index in face) for face in mesh_faces)
                offset += len(mesh_vertices)
            return tuple(vertices), tuple(faces)

        outer_v, outer_f = box((-3,-3,-2),(4,3,4))
        design_v, design_f = box((-1,-1,-1),(2,1,1))
        fixed = combine(box((-1,1,-1),(2,1.125,1)), box((-1,-1.125,-1),(2,-1,1)))
        displacement = box((-1,-2,-1),(2,2,1))
        request = with_measured_buoyancy_geometry(
            SpinItFasterRequest(
                outer_v, outer_f, design_v, design_f, material_density=1.0,
                original_units_per_normalized_unit=2.5,
            ),
            *displacement, *fixed, waterline_origin=(0,0,0),
            waterline_normal=(0,0,1), liquid_density=0.625,
        )
        solved = solve_spin_it_faster_buoyant_plane(
            request, initial_parameters=(0.1,0.2,-0.1),
            maximum_iterations=40, constraint_tolerance=1e-11,
        )
        self.assertEqual(solved["status"], "CONVERGED")
        vertices = np.asarray(solved["variable_material_vertices"], dtype=np.float64)
        faces = np.asarray(solved["variable_material_faces"], dtype=np.int64)
        self.assertGreater(len(faces), 0)

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            obj, brep, step = root / "solver-material.obj", root / "solver-material.brep", root / "solver-material.step"
            with obj.open("w", encoding="ascii") as stream:
                for vertex in vertices:
                    stream.write("v %.17g %.17g %.17g\n" % tuple(vertex))
                for face in faces:
                    stream.write("f %d %d %d\n" % tuple(face + 1))
            completed = subprocess.run(
                [str(PROBE), "--mesh", str(obj),
                 str(request.original_units_per_normalized_unit), str(brep), str(step)],
                check=False, capture_output=True, text=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            result = json.loads(completed.stdout)
            self.assertEqual(result["status"], "SOLVER_MESH_OCCT_STEP_BREP_ROUNDTRIP_VALIDATED")
            self.assertEqual(result["coordinate_scale"], 2.5)
            self.assertTrue(brep.is_file() and brep.stat().st_size > 0)
            self.assertTrue(step.is_file() and step.stat().st_size > 0)
            scale = request.original_units_per_normalized_unit
            scale_factors = np.asarray(
                (scale**3, scale**4, scale**4, scale**4,
                 scale**5, scale**5, scale**5, scale**5, scale**5, scale**5),
                dtype=np.float64,
            )
            expected = np.asarray(solved["variable_raw_moments"], dtype=np.float64) * scale_factors
            for stage in ("original", "brep", "step"):
                actual = np.asarray(result["raw_moments"][stage], dtype=np.float64)
                np.testing.assert_allclose(actual, expected, rtol=1e-8, atol=1e-9)

    def test_mesh_path_rejects_open_surface_instead_of_fabricating_a_solid(self) -> None:
        self.assertTrue(PROBE.is_file(), "fixed-OCCT roundtrip probe is not built")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            obj, brep, step = root / "open.obj", root / "open.brep", root / "open.step"
            obj.write_text(
                "v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 3\n",
                encoding="ascii",
            )
            completed = subprocess.run(
                [str(PROBE), "--mesh", str(obj), "1.0", str(brep), str(step)],
                check=False, capture_output=True, text=True,
            )
            self.assertNotEqual(completed.returncode, 0)
            self.assertIn("closed shell", completed.stderr)
            self.assertFalse(brep.exists())
            self.assertFalse(step.exists())

    def test_mesh_path_rejects_nonpositive_or_nonfinite_unit_scale(self) -> None:
        self.assertTrue(PROBE.is_file(), "fixed-OCCT roundtrip probe is not built")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            obj = root / "closed.obj"
            obj.write_text(
                "v 0 0 0\nv 1 0 0\nv 1 1 0\nv 0 1 0\n"
                "v 0 0 1\nv 1 0 1\nv 1 1 1\nv 0 1 1\n"
                "f 1 3 2\nf 1 4 3\nf 5 6 7\nf 5 7 8\n"
                "f 1 2 6\nf 1 6 5\nf 4 8 7\nf 4 7 3\n"
                "f 1 5 8\nf 1 8 4\nf 2 3 7\nf 2 7 6\n",
                encoding="ascii",
            )
            for index, scale in enumerate(("0", "-1", "nan", "inf")):
                brep, step = root / f"bad-{index}.brep", root / f"bad-{index}.step"
                completed = subprocess.run(
                    [str(PROBE), "--mesh", str(obj), scale, str(brep), str(step)],
                    check=False, capture_output=True, text=True,
                )
                self.assertNotEqual(completed.returncode, 0)
                self.assertIn("coordinate scale", completed.stderr)
                self.assertFalse(brep.exists())
                self.assertFalse(step.exists())

    def test_sealed_faster_result_exports_locked_3d_material_to_step(self) -> None:
        """The 3D CAD handoff must be bound to a sealed native result."""
        export = getattr(make_it_stand,
                         "export_spin_it_faster_quadric_material_step", None)
        self.assertTrue(callable(export), "sealed 3D Faster-to-OCCT handoff is missing")
        from test_exact_material_region import box
        vertices, faces = box((-2.,) * 3, (2.,) * 3)
        design = box((-1.,) * 3, (1.,) * 3)
        fixed_regions = [
            ("shell_x_minus", box((-2., -2., -2.), (-1., 2., 2.)), 2.0),
            ("shell_x_plus", box((1., -2., -2.), (2., 2., 2.)), 2.0),
            ("shell_y_minus", box((-1., -2., -2.), (1., -1., 2.)), 2.0),
            ("shell_y_plus", box((-1., 1., -2.), (1., 2., 2.)), 2.0),
            ("shell_z_minus", box((-1., -1., -2.), (1., 1., -1.)), 2.0),
            ("shell_z_plus", box((-1., -1., 1.), (1., 1., 2.)), 2.0),
        ]
        problem = SpinArtifactProblem(
            original=(vertices, faces), fixed_regions=fixed_regions,
            design_regions=[("design", design, 2.0)],
            length_unit="m", mass_unit="kg",
            target={"mode": "YO_YO", "axis": [0., 0., 1.],
                    "origin": [0., 0., 0.], "support_points": []},
            minimum_outer_wall=0.1,
        )
        # Cut a small analytic ellipsoid from the actual normalized design.
        radius = 0.01
        center = np.zeros(3)
        coefficients = np.asarray((center @ center - radius**2,
                                   *(-2. * center), 1., 1., 1.,
                                   0., 0., 0.), dtype=np.float64)
        coefficients /= np.linalg.norm(coefficients)
        candidate = evaluate_and_seal_faster_candidate(
            problem, coefficients,
            evaluator_settings=spin_api.SpinEvaluatorConfig(
                derivative_method="DISCRETE_ELLIPSOID", icosphere_subdivisions=3),
            solver_settings=spin_api.SpinSolverConfig(
                equation_profile="DUAL_MONOMIAL_RAW_V1", maximum_iterations=1,
                maximum_restarts=1, maximum_function_evaluations_per_restart=1),
        )
        self.assertEqual(candidate.status, "NOT_OPTIMIZED_CANDIDATE")
        self.assertIsNotNone(candidate.artifact)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            cad = export(candidate, occt_executable=PROBE,
                         brep_path=root / "sealed.brep", step_path=root / "sealed.step")
            self.assertTrue((root / "sealed.brep").is_file())
            self.assertTrue((root / "sealed.step").is_file())
            sidecar = json.loads((root / "sealed.step.identity.json").read_text(
                encoding="ascii"))
            self.assertEqual(sidecar["cad_input_digest"], cad["cad_input_digest"])
            self.assertRegex(sidecar["occt_executable_sha256"], r"^[0-9a-f]{64}$")
            self.assertRegex(sidecar["brep_sha256"], r"^[0-9a-f]{64}$")
            self.assertRegex(sidecar["step_sha256"], r"^[0-9a-f]{64}$")
        self.assertEqual(cad["classification"], "METHOD_VALIDATION_ONLY")
        self.assertEqual(cad["status"], "FASTER_3D_MATERIAL_STEP_BREP_VALIDATED")
        self.assertEqual(cad["occt_version"], "7.9.3")
        self.assertEqual(np.asarray(cad["cad_source_mass_moments"]).shape, (10,))
        self.assertEqual(cad["material_region_ids"], [
            "design", "shell_x_minus", "shell_x_plus", "shell_y_minus",
            "shell_y_plus", "shell_z_minus", "shell_z_plus"])
        self.assertEqual(cad["entity_roles"], ["design", "fixed", "fixed", "fixed",
                                                "fixed", "fixed", "fixed"])
        np.testing.assert_array_equal(cad["entity_densities"], np.ones(7))
        self.assertLessEqual(cad["maximum_normalized_material_moment_error"], 1e-6)
        tampered = object.__new__(type(candidate))
        for name in ("status", "artifact", "seal_error", "scientific_digest"):
            object.__setattr__(tampered, name, getattr(candidate, name))
        altered = candidate.identity
        altered["coefficients"][0] += 1e-3
        object.__setattr__(tampered, "_identity", json.dumps(
            altered, sort_keys=True, separators=(",", ":")))
        with self.assertRaisesRegex(ValueError, "identity digest mismatch"):
            export(tampered, occt_executable=PROBE,
                   brep_path=Path("unused.brep"), step_path=Path("unused.step"))

    def test_reduced_stratum_cylinder_candidate_seals_after_exact_partition(self) -> None:
        """A sealed reduced-stratum cylinder candidate round-trips through CAD."""
        export = getattr(make_it_stand,
                         "export_spin_it_faster_quadric_material_step", None)
        self.assertTrue(callable(export))
        from test_exact_material_region import box
        vertices, faces = box((-2.,) * 3, (2.,) * 3)
        design = box((-1.,) * 3, (1.,) * 3)
        fixed_regions = [
            ("shell_x_minus", box((-2., -2., -2.), (-1., 2., 2.)), 2.0),
            ("shell_x_plus", box((1., -2., -2.), (2., 2., 2.)), 2.0),
            ("shell_y_minus", box((-1., -2., -2.), (1., -1., 2.)), 2.0),
            ("shell_y_plus", box((-1., 1., -2.), (1., 2., 2.)), 2.0),
            ("shell_z_minus", box((-1., -1., -2.), (1., 1., -1.)), 2.0),
            ("shell_z_plus", box((-1., -1., 1.), (1., 1., 2.)), 2.0),
        ]
        problem = SpinArtifactProblem(
            original=(vertices, faces), fixed_regions=fixed_regions,
            design_regions=[("design", design, 2.0)], length_unit="m", mass_unit="kg",
            target={"mode": "YO_YO", "axis": [0., 0., 1.],
                    "origin": [0., 0., 0.], "support_points": []},
            minimum_outer_wall=0.1)
        # q=x^2+2y^2-r^2 is an anisotropic elliptic cylinder, not a full-rank
        # ellipsoid. Sign-invariant triangulation must make p/-p a valid partition.
        radius = 0.005
        coefficients = np.asarray((-radius**2, 0., 0., 0., 1., 2., 0., 0., 0., 0.),
                                  dtype=np.float64)
        coefficients /= np.linalg.norm(coefficients)
        with tempfile.TemporaryDirectory() as directory:
            candidate = evaluate_and_seal_faster_candidate(
                problem, coefficients,
                evaluator_settings=spin_api.SpinEvaluatorConfig(
                    derivative_method="DISCRETE_REDUCED_STRATUM",
                    cylinder_angular_resolution=96,
                    truncation_radius=2.,
                    cap_radius_policy="GROW_TO_CERTIFIED_V1"),
                solver_settings=spin_api.SpinSolverConfig(
                    equation_profile="DUAL_MONOMIAL_RAW_V1", maximum_iterations=1,
                    maximum_restarts=1, maximum_function_evaluations_per_restart=1))
            self.assertEqual(candidate.status, "NOT_OPTIMIZED_CANDIDATE")
            self.assertIsNotNone(candidate.artifact)
            retained_evaluation = candidate.identity["native_evaluations"][0]
            self.assertEqual(retained_evaluation["interface_type"], "ELLIPTIC_CYLINDER")
            self.assertFalse(retained_evaluation["derivative_defined"])
            cad = export(candidate, occt_executable=PROBE,
                         brep_path=Path(directory) / "cylinder.brep",
                         step_path=Path(directory) / "cylinder.step")
            self.assertEqual(cad["quadric_type"], "ELLIPTIC_CYLINDER")
            self.assertLessEqual(cad["maximum_normalized_material_moment_error"], 1e-6)
            self.assertTrue((Path(directory) / "cylinder.brep").is_file())
            self.assertTrue((Path(directory) / "cylinder.step").is_file())

if __name__ == "__main__":
    unittest.main()
