"""Analytic motion and moment transport for retained Cherchi constructions."""

import unittest

from test_rational_cherchi_transport import RationalBoundaryTransportTests

import numpy as np

import _make_it_stand as native
from make_it_stand.verification.cherchi_motion_differential import (
    CherchiDifferentialUnavailable,
    differentiate_cherchi_point,
    differentiate_polyhedral_raw_moments,
    normalize_closed_component_orientation,
    recover_construction_vertex_motion,
)
from make_it_stand.verification.mass_moments import TriangleMesh, integrate_geometric_moments


def line_plane_intersection(points):
    line_a, line_b, plane_a, plane_b, plane_c = points
    direction = line_b - line_a
    normal = np.cross(plane_b - plane_a, plane_c - plane_a)
    parameter = np.dot(normal, plane_a - line_a) / np.dot(normal, direction)
    return line_a + parameter * direction


class CherchiPointMotionTests(unittest.TestCase):
    def test_lnc_fixed_parameter_motion_matches_finite_difference(self):
        points = np.asarray(((0.4, -0.2, 1.1), (-0.7, 0.8, -0.3)))
        velocities = np.asarray((
            ((0.2, -0.1), (0.3, 0.0), (-0.4, 0.2)),
            ((-0.1, 0.5), (0.0, -0.2), (0.6, 0.1)),
        ))
        parameter = 0.375
        point, velocity = differentiate_cherchi_point(
            "LNC", points, velocities, parameter_t=parameter,
        )
        self.assertEqual(velocity.shape, (3, 2))
        np.testing.assert_allclose(point, (1 - parameter) * points[0] + parameter * points[1])
        np.testing.assert_allclose(
            velocity, (1 - parameter) * velocities[0] + parameter * velocities[1],
            rtol=0, atol=1e-15,
        )
        for step in (1e-4, 3e-5, 1e-5):
            plus = (1 - parameter) * (points[0] + step * velocities[0, :, 0]) + \
                parameter * (points[1] + step * velocities[1, :, 0])
            minus = (1 - parameter) * (points[0] - step * velocities[0, :, 0]) + \
                parameter * (points[1] - step * velocities[1, :, 0])
            np.testing.assert_allclose(velocity[:, 0], (plus - minus) / (2 * step), atol=2e-12)

    def test_lnc_motion_rejects_invalid_parameter(self):
        points = np.asarray(((0., 0., 0.), (1., 2., 3.)))
        for parameter in (True, "0.5", np.nan, -0.1, 1.1):
            with self.subTest(parameter=parameter), self.assertRaises(CherchiDifferentialUnavailable):
                differentiate_cherchi_point(
                    "LNC", points, np.zeros((2, 3, 1)), parameter_t=parameter,
                )

    def test_lpi_motion_matches_finite_difference_of_same_five_point_construction(self):
        points = np.asarray((
            (-1.0, 0.2, 1.5), (1.2, -0.1, -0.5),
            (-0.5, -0.8, 0.0), (0.9, -0.6, 0.1), (-0.2, 1.1, -0.2),
        ))
        velocities = np.asarray((
            ((0.1, -0.2), (0.3, 0.0), (-0.1, 0.4)),
            ((-0.2, 0.1), (0.0, 0.3), (0.2, -0.1)),
            ((0.0, 0.1), (-0.1, 0.0), (0.0, 0.2)),
            ((0.2, -0.1), (0.0, 0.2), (0.1, 0.0)),
            ((-0.1, 0.0), (0.2, -0.1), (0.0, 0.1)),
        ))
        expected_point, velocity = differentiate_cherchi_point("LPI", points, velocities)
        self.assertEqual(velocity.shape, (3, 2))
        for step in (1e-4, 3e-5, 1e-5):
            plus = line_plane_intersection(points + step * velocities[:, :, 0])
            minus = line_plane_intersection(points - step * velocities[:, :, 0])
            finite_difference = (plus - minus) / (2.0 * step)
            np.testing.assert_allclose(expected_point, line_plane_intersection(points), atol=1e-14)
            self.assertLess(np.linalg.norm(velocity[:, 0] - finite_difference), 1e-8)

    def test_parallel_lpi_is_rejected_instead_of_returning_nonfinite_velocity(self):
        points = np.asarray(((0., 0., 1.), (1., 0., 1.),
                             (0., 0., 0.), (1., 0., 0.), (0., 1., 0.)))
        velocities = np.zeros((5, 3, 1))
        with self.assertRaises(CherchiDifferentialUnavailable):
            differentiate_cherchi_point("LPI", points, velocities)

    def test_tpi_motion_matches_finite_difference_of_three_moving_planes(self):
        points = np.asarray((
            (1., 0., 0.), (1., 1., 0.), (1., 0., 1.),
            (0., 2., 0.), (1., 2., 0.), (0., 2., 1.),
            (0., 0., 3.), (1., 0., 3.), (0., 1., 3.),
        ))
        rng = np.random.default_rng(17)
        velocities = rng.normal(size=(9, 3, 2)) * 0.1
        point, jacobian = differentiate_cherchi_point("TPI", points, velocities)
        np.testing.assert_allclose(point, (1., 2., 3.), atol=1e-14)
        self.assertEqual(jacobian.shape, (3, 2))
        for step in (1e-4, 3e-5, 1e-5):
            plus = differentiate_cherchi_point(
                "TPI", points + step * velocities[:, :, 0], np.zeros_like(velocities)
            )[0]
            minus = differentiate_cherchi_point(
                "TPI", points - step * velocities[:, :, 0], np.zeros_like(velocities)
            )[0]
            self.assertLess(np.linalg.norm(jacobian[:, 0] - (plus - minus) / (2 * step)), 2e-8)

    def test_construction_recovery_requires_unique_sources_and_matches_output_point(self):
        points = np.asarray((
            (-1.0, 0.2, 1.5), (1.2, -0.1, -0.5),
            (-0.5, -0.8, 0.0), (0.9, -0.6, 0.1), (-0.2, 1.1, -0.2),
        ))
        velocities = np.zeros((5, 3, 2))
        velocities[0, 0, :] = (0.1, -0.2)
        velocities[2, 1, :] = (-0.1, 0.3)
        point, expected = differentiate_cherchi_point("LPI", points, velocities)
        record = {
            "vertex_index": 0, "point_type": "LPI",
            "operands": [{"status": "UNIQUE", "source_vertices": [
                {"input_object": 1, "input_vertex": i}
            ]} for i in range(5)],
            "operand_coordinates": points.tolist(),
        }
        recovered = recover_construction_vertex_motion(
            [record], np.asarray(((0., 0., 0.),)), points, velocities, point[None, :]
        )
        np.testing.assert_allclose(recovered[0], expected, atol=1e-14)
        record["operands"][0]["status"] = "AMBIGUOUS"
        with self.assertRaises(CherchiDifferentialUnavailable):
            recover_construction_vertex_motion(
                [record], np.asarray(((0., 0., 0.),)), points, velocities, point[None, :]
            )

    def test_lnc_recovery_uses_retained_fixed_parameter_and_unique_sources(self):
        design = np.asarray(((0.2, 0.4, -0.1),))
        source = np.asarray(((1.2, -0.3, 0.5),))
        motion = np.asarray((((0.3,), (-0.2,), (0.1,)),))
        parameter = 0.25
        point = (1 - parameter) * design[0] + parameter * source[0]
        record = {
            "vertex_index": 0, "point_type": "LNC", "parameter_t": parameter,
            "operands": [
                {"status": "UNIQUE", "source_vertices": [{"input_object": 0, "input_vertex": 0}]},
                {"status": "UNIQUE", "source_vertices": [{"input_object": 1, "input_vertex": 0}]},
            ],
            "operand_coordinates": [design[0].tolist(), source[0].tolist()],
        }
        recovered = recover_construction_vertex_motion(
            [record], design, source, motion, point[None, :],
        )
        np.testing.assert_allclose(recovered[0, :, 0], parameter * motion[0, :, 0], atol=1e-15)
        record["parameter_t"] = None
        with self.assertRaises(CherchiDifferentialUnavailable):
            recover_construction_vertex_motion([record], design, source, motion, point[None, :])


class PolyhedralMomentMotionTests(unittest.TestCase):
    @staticmethod
    def independent_raw_moments(vertices, faces):
        moments = integrate_geometric_moments(TriangleMesh(
            tuple(tuple(float(value) for value in row) for row in vertices),
            tuple(tuple(int(value) for value in face) for face in faces),
        ))
        q = moments.second_moment
        return np.asarray((
            moments.volume, *moments.first_moment,
            q[0][0], q[1][1], q[2][2], q[1][2], q[0][2], q[0][1],
        ))

    def test_uniform_translation_differential_matches_independent_integrator(self):
        vertices = np.asarray(((0., 0., 0.), (1.2, 0., 0.),
                               (0., 0.8, 0.), (0., 0., 1.1)))
        faces = np.asarray(((0, 2, 1), (0, 1, 3), (1, 2, 3), (2, 0, 3)), dtype=np.int32)
        velocities = np.repeat(np.asarray((((0.3,), (-0.2,), (0.1,)),)), 4, axis=0)
        density = 1.7
        jacobian = differentiate_polyhedral_raw_moments(vertices, faces, velocities, density)
        self.assertEqual(jacobian.shape, (10, 1))
        step = 1e-5
        offset = step * velocities[:, :, 0]
        plus = self.independent_raw_moments(vertices + offset, faces)
        minus = self.independent_raw_moments(vertices - offset, faces)
        finite_difference = density * (plus - minus) / (2.0 * step)
        np.testing.assert_allclose(jacobian[:, 0], finite_difference, rtol=2e-8, atol=2e-9)

    def test_open_boolean_boundary_reports_edge_witness_counts(self):
        tetra_vertices = np.asarray((
            (0., 0., 0.), (1., 0., 0.), (0., 1., 0.), (0., 0., 1.),
        ))
        tetra_faces = np.asarray((
            (0, 2, 1), (0, 1, 3), (1, 2, 3), (2, 0, 3),
        ), dtype=np.int64)
        vertices = np.vstack((tetra_vertices, tetra_vertices + 3.0))
        open_faces = np.vstack((tetra_faces, tetra_faces[:-1] + 4))
        motion = np.zeros((len(vertices), 3, 1))

        with self.assertRaisesRegex(
            CherchiDifferentialUnavailable,
            r"boundary_edges=3; nonmanifold_edges=0; inconsistent_edges=0; sample_edges=",
        ):
            differentiate_polyhedral_raw_moments(vertices, open_faces, motion)

    def test_orientation_normalization_repairs_local_winding_but_preserves_component_signs(self):
        # Outer tetrahedron is outward oriented; nested inner tetrahedron is
        # inward oriented to represent a cavity. Flip a few faces as seen in
        # the pinned Boolean backend's emitted boundary.
        def outward_tetra(points):
            faces = np.asarray(((0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3)), dtype=np.int64)
            center = points.mean(axis=0)
            for index, face in enumerate(faces):
                a, b, c = points[face]
                if np.dot(np.cross(b - a, c - a), (a + b + c) / 3.0 - center) < 0:
                    faces[index, [1, 2]] = faces[index, [2, 1]]
            return faces

        outer = np.asarray(((0., 0., 0.), (2., 0., 0.), (0., 2., 0.), (0., 0., 2.)))
        inner = 0.2 + 0.3 * np.asarray(((0., 0., 0.), (2., 0., 0.), (0., 2., 0.), (0., 0., 2.)))
        vertices = np.vstack((outer, inner))
        outer_faces = outward_tetra(outer)
        inner_faces = outward_tetra(inner)[:, (0, 2, 1)] + 4
        expected_faces = np.vstack((outer_faces, inner_faces))
        emitted_faces = expected_faces.copy()
        # Choose faces whose individual signed-volume contributions do not
        # overwhelm the emitted component-wide orientation evidence.
        emitted_faces[[1, 5], :] = emitted_faces[[1, 5]][:, (0, 2, 1)]

        normalized, audit = normalize_closed_component_orientation(vertices, emitted_faces)
        np.testing.assert_array_equal(normalized, expected_faces)
        self.assertEqual(audit["component_count"], 2)
        self.assertEqual(audit["flipped_face_count"], 2)
        self.assertEqual([c["orientation_sign"] for c in audit["components"]], [1, -1])
        np.testing.assert_allclose(
            self.independent_raw_moments(vertices, normalized),
            self.independent_raw_moments(vertices, expected_faces), rtol=0, atol=1e-13,
        )

    def test_orientation_normalization_rejects_open_and_nonmanifold_boundaries(self):
        tetra = np.asarray(((0., 0., 0.), (1., 0., 0.), (0., 1., 0.), (0., 0., 1.)))
        faces = np.asarray(((0, 2, 1), (0, 1, 3), (1, 2, 3), (2, 0, 3)), dtype=np.int64)
        two_tetra = np.vstack((tetra, tetra + 3.0))
        one_open_component = np.vstack((faces, faces[:-1] + 4))
        with self.assertRaisesRegex(
            CherchiDifferentialUnavailable, r"boundary_edges=3.*nonmanifold_edges=0",
        ):
            normalize_closed_component_orientation(two_tetra, one_open_component)
        with self.assertRaisesRegex(
            CherchiDifferentialUnavailable, r"boundary_edges=0.*nonmanifold_edges=3",
        ):
            normalize_closed_component_orientation(tetra, np.vstack((faces, faces[0])))

    def test_exact_nonzero_sliver_face_is_not_misclassified_as_degenerate(self):
        # This closed unit cube has one top-face triangle with nonzero area
        # 1e-15, below the former scale-relative floating threshold. The mesh
        # remains a valid oriented boundary with unit volume.
        vertices = np.asarray((
            (0., 0., 0.), (1., 0., 0.), (1., 1., 0.), (0., 1., 0.),
            (0., 0., 1.), (1., 0., 1.), (1., 1., 1.), (0., 1., 1.),
            (.5, 1e-15, 1.),
        ))
        faces = np.asarray((
            (0, 2, 1), (0, 3, 2),
            (4, 5, 8), (5, 6, 8), (6, 4, 8), (4, 6, 7),
            (0, 1, 5), (0, 5, 4),
            (1, 2, 6), (1, 6, 5),
            (2, 3, 7), (2, 7, 6),
            (3, 0, 4), (3, 4, 7),
        ), dtype=np.int64)

        oriented, audit = normalize_closed_component_orientation(vertices, faces)
        np.testing.assert_array_equal(oriented, faces)
        self.assertEqual(audit["component_count"], 1)
        self.assertEqual(audit["components"][0]["orientation_sign"], 1)

        motion = np.zeros((len(vertices), 3, 1))
        motion[8, 2, 0] = 1.0
        jacobian = differentiate_polyhedral_raw_moments(vertices, oriented, motion)
        self.assertTrue(np.isfinite(jacobian).all())
        for step in (1e-4, 3e-5, 1e-5):
            plus = self.independent_raw_moments(vertices + step * motion[:, :, 0], oriented)
            minus = self.independent_raw_moments(vertices - step * motion[:, :, 0], oriented)
            finite_difference = (plus - minus) / (2.0 * step)
            np.testing.assert_allclose(jacobian[:, 0], finite_difference, rtol=2e-7, atol=2e-8)

    def test_exactly_collinear_face_remains_rejected(self):
        vertices = np.asarray((
            (0., 0., 0.), (1., 0., 0.), (1., 1., 0.), (0., 1., 0.),
            (0., 0., 1.), (1., 0., 1.), (1., 1., 1.), (0., 1., 1.),
            (.5, 0., 1.),
        ))
        faces = np.asarray((
            (0, 2, 1), (0, 3, 2),
            (4, 5, 8), (5, 6, 8), (6, 4, 8), (4, 6, 7),
            (0, 1, 5), (0, 5, 4),
            (1, 2, 6), (1, 6, 5),
            (2, 3, 7), (2, 7, 6),
            (3, 0, 4), (3, 4, 7),
        ), dtype=np.int64)
        with self.assertRaisesRegex(CherchiDifferentialUnavailable, "zero-area face"):
            normalize_closed_component_orientation(vertices, faces)

    def test_small_closed_component_uses_its_own_volume_scale(self):
        cube_faces = np.asarray((
            (0, 2, 1), (0, 3, 2), (4, 5, 6), (4, 6, 7),
            (0, 1, 5), (0, 5, 4), (1, 2, 6), (1, 6, 5),
            (2, 3, 7), (2, 7, 6), (3, 0, 4), (3, 4, 7),
        ), dtype=np.int64)
        vertices = []
        faces = []
        centers = (-0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75)
        half_side = 0.03
        for x in centers:
            for y in centers:
                for z in centers:
                    if x == y == z == 0.0:
                        continue
                    base = len(vertices)
                    vertices.extend((
                        (x-half_side, y-half_side, z-half_side),
                        (x+half_side, y-half_side, z-half_side),
                        (x+half_side, y+half_side, z-half_side),
                        (x-half_side, y+half_side, z-half_side),
                        (x-half_side, y-half_side, z+half_side),
                        (x+half_side, y-half_side, z+half_side),
                        (x+half_side, y+half_side, z+half_side),
                        (x-half_side, y+half_side, z+half_side),
                    ))
                    faces.extend(cube_faces + base)

        # This is the exact scale of the four-face component rejected in the
        # seed-0 Boolean output: nonzero volume, but below the old whole-mesh
        # triangle-count-dependent volume floor.
        tiny_base = len(vertices)
        vertices.extend(((0., 0., 0.), (1e-3, 0., 0.), (0., 1e-3, 0.), (0., 0., 1e-3)))
        faces.extend(np.asarray(((0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3))) + tiny_base)

        normalized, audit = normalize_closed_component_orientation(
            np.asarray(vertices, dtype=np.float64), np.asarray(faces, dtype=np.int64),
        )
        np.testing.assert_array_equal(normalized, np.asarray(faces, dtype=np.int64))
        self.assertEqual(audit["component_count"], 343)
        tiny = audit["components"][-1]
        self.assertEqual(tiny["face_count"], 4)
        self.assertGreater(tiny["emitted_signed_volume"], 1e-10)
        self.assertEqual(tiny["orientation_sign"], 1)


class NativeQuadricSourceMotionTests(unittest.TestCase):
    def test_paper_projective_config_drives_actual_source_and_provenance(self):
        from make_it_stand import spin
        config=spin.SpinEvaluatorConfig(paper_projective_hyperboloid=True,
            derivative_method="DISCRETE_FULL_RANK_QUADRIC",
            icosphere_subdivisions=2,maximum_icosphere_subdivisions=2,
            truncation_radius=16.,maximum_radius_doublings=8)
        settings=spin._native_evaluator_settings(native,config).triangulation
        q=np.array((-1.,.1,-.12,.08,1.,1.3,-.9,.03,-.04,.02))
        domain=np.array(((-.2,-.2,-.2),(.2,.2,.2)))
        result=native.evaluate_spin_it_faster_source_motion(q,domain,settings)
        self.assertTrue(result['defined'],result['diagnostic'])
        self.assertTrue(result['paper_projective_hyperboloid'])
        rays=np.asarray(result['source_projective_rays'])
        self.assertEqual(rays.shape,(len(result['vertices']),4))
        self.assertGreater(np.count_nonzero(rays[:,0]>=0),0)
        self.assertEqual(len(result['source_normal_vertices']),len(rays))
        self.assertEqual(len(result['source_normal_edges']),len(rays))
        settings.paper_projective_hyperboloid=False
        legacy=native.evaluate_spin_it_faster_source_motion(q,domain,settings)
        self.assertTrue(legacy['defined'],legacy['diagnostic'])
        self.assertFalse(legacy['paper_projective_hyperboloid'])
        self.assertFalse(np.array_equal(result['vertices'],legacy['vertices']) and
                         np.array_equal(result['faces'],legacy['faces']),
                         'map tags changed without changing the actual source geometry')
        self.assertNotEqual(result['triangulation_hash'],legacy['triangulation_hash'])
        self.assertNotEqual(result['scientific_hash'],legacy['scientific_hash'])
        self.assertFalse(spin.SpinEvaluatorConfig().paper_projective_hyperboloid)
        for bad in (1,0,None,'true'):
            with self.assertRaises(ValueError):
                spin.SpinEvaluatorConfig(paper_projective_hyperboloid=bad)

    def test_native_motion_returns_the_source_mesh_and_full_tangent_columns(self):
        coefficients = np.asarray((-1., 0., 0., 0., 1., 1., 1., 0., 0., 0.))
        domain = np.asarray(((-2., -2., -2.), (2., 2., 2.)))
        settings = native.SpinTriangulationSettings()
        settings.icosphere_subdivisions = 1
        result = native.evaluate_spin_it_faster_source_motion(coefficients, domain, settings)
        self.assertTrue(result["defined"], result["diagnostic"])
        vertices = np.asarray(result["vertices"])
        faces = np.asarray(result["faces"])
        velocity = np.asarray(result["vertex_velocity"])
        self.assertEqual(vertices.shape[1], 3)
        self.assertEqual(faces.shape[1], 3)
        self.assertEqual(velocity.shape, (3 * len(vertices), 10))
        self.assertTrue(np.isfinite(velocity).all())
        self.assertEqual(len(result["provenance"]), len(vertices))
        self.assertTrue(result["scientific_hash"])
        self.assertTrue(result["side_selection_defined"], result["side_selection_diagnostic"])
        self.assertFalse(result["interior_is_positive"])
        self.assertEqual(
            result["side_selection_method"],
            "interior-centroid-winding-or-oriented-interface-v1",
        )


if __name__ == "__main__":
    unittest.main()
