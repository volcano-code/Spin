"""Paper-coordinate and actual support-flattening regression tests."""

from __future__ import annotations

import importlib.util
from dataclasses import replace
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from make_it_stand.model import ModelConfig, RunConfig
from make_it_stand.paper_compat import (
    PaperCompatibilityModel,
    PaperNormalization,
    PAPER_RUN_CONFIG,
    RELEASED_BETA_RUN_CONFIG,
    SupportFlattening,
    released_beta_run_config,
    paper_run_config,
)
from make_it_stand.sealed_artifact import (
    SealedBalanceArtifact,
    _repair_structured_cavity,
    _structured_boundary_issues,
    _structured_cavity_boundary,
    _structured_occupancy_hash,
    _triangulate_simple_polygon,
    load_sealed_pybullet_entry,
)


def _native_available() -> bool:
    for name in ("_make_it_stand", "make_it_stand._make_it_stand"):
        try:
            if importlib.util.find_spec(name) is not None:
                return True
        except (ImportError, ModuleNotFoundError, ValueError):
            continue
    return False


def _sloped_box() -> tuple[np.ndarray, np.ndarray]:
    vertices = np.asarray(
        [
            [0.0, 10.0, 0.0],
            [2.0, 10.1, 0.0],
            [2.0, 11.1, 0.0],
            [0.0, 11.0, 0.0],
            [0.0, 10.0, 1.0],
            [2.0, 10.1, 1.0],
            [2.0, 11.1, 1.0],
            [0.0, 11.0, 1.0],
        ],
        dtype=np.float64,
    )
    faces = np.asarray(
        [
            [0, 2, 1], [0, 3, 2],
            [4, 5, 6], [4, 6, 7],
            [0, 1, 5], [0, 5, 4],
            [3, 7, 6], [3, 6, 2],
            [0, 4, 7], [0, 7, 3],
            [1, 2, 6], [1, 6, 5],
        ],
        dtype=np.int32,
    )
    return vertices, faces


class PaperTransformTests(unittest.TestCase):
    def test_cavity_repair_allows_one_neutral_topology_transition(self) -> None:
        dimensions = (2, 2, 2)
        void = {
            (0, 0, 1),
            (0, 1, 0),
            (0, 1, 1),
            (1, 0, 0),
            (1, 0, 1),
            (1, 1, 0),
        }
        repaired, audit = _repair_structured_cavity(dimensions, void)
        bad_edges, bad_vertices = _structured_boundary_issues(
            _structured_cavity_boundary(dimensions, repaired)
        )
        self.assertEqual(bad_edges, [])
        self.assertEqual(bad_vertices, [])
        self.assertEqual(len(audit["refilled_cells_xyz"]), 2)

    def test_planar_ear_clipping_preserves_concave_polygon_boundary(self) -> None:
        vertices = np.asarray(
            [
                [0.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],
                [2.0, 2.0, 0.0],
                [1.0, 0.75, 0.0],
                [0.0, 2.0, 0.0],
            ],
            dtype=np.float64,
        )
        loop = (0, 1, 2, 3, 4)
        triangles = _triangulate_simple_polygon(
            vertices,
            loop,
            (0.0, 0.0, 1.0),
        )

        self.assertEqual(triangles.shape, (3, 3))
        counts: dict[tuple[int, int], int] = {}
        for face in triangles:
            for first, second in (
                (int(face[0]), int(face[1])),
                (int(face[1]), int(face[2])),
                (int(face[2]), int(face[0])),
            ):
                edge = tuple(sorted((first, second)))
                counts[edge] = counts.get(edge, 0) + 1
        boundary = {edge for edge, count in counts.items() if count == 1}
        expected = {
            tuple(sorted((loop[index], loop[(index + 1) % len(loop)])))
            for index in range(len(loop))
        }
        self.assertEqual(boundary, expected)
        self.assertTrue(all(count in (1, 2) for count in counts.values()))

    def test_released_beta_profile_is_explicit_and_modern_defaults_stay_modern(
        self,
    ) -> None:
        modern = RunConfig()
        paper = paper_run_config(maximum_carving_blocks=9)
        beta = released_beta_run_config(maximum_carving_blocks=7)

        self.assertEqual(paper, replace(PAPER_RUN_CONFIG, maximum_carving_blocks=9))
        self.assertEqual(paper.initial_mu, 0.75)
        self.assertEqual(paper.shape_scale_lambda, 20.0)
        self.assertEqual(paper.required_relative_decrease, 0.01)
        self.assertFalse(paper.freeze_targets_after_initial_refresh)
        self.assertFalse(paper.released_beta_semantics)
        self.assertFalse(paper.released_beta_outer_loop)

        self.assertEqual(RELEASED_BETA_RUN_CONFIG.initial_mu, 0.75)
        self.assertEqual(beta.initial_mu, 0.75)
        self.assertEqual(beta.shape_scale_lambda, 20.0)
        self.assertEqual(beta.required_relative_decrease, 0.03)
        self.assertEqual(beta.initial_step, 1.0)
        self.assertEqual(beta.reduction_ratio, 0.8)
        self.assertEqual(beta.minimum_step, 0.4)
        self.assertEqual(beta.mu_decrement, 0.05)
        self.assertTrue(beta.freeze_targets_after_initial_refresh)
        self.assertTrue(beta.released_beta_semantics)
        self.assertTrue(beta.released_beta_outer_loop)
        self.assertEqual(beta.maximum_carving_blocks, 7)
        self.assertEqual(RELEASED_BETA_RUN_CONFIG.maximum_carving_blocks, 20)

        self.assertEqual(modern.initial_mu, 0.2)
        self.assertEqual(modern.shape_scale_lambda, 1.0)
        self.assertEqual(modern.required_relative_decrease, 0.01)
        self.assertFalse(modern.freeze_targets_after_initial_refresh)
        self.assertFalse(modern.released_beta_semantics)
        self.assertFalse(modern.released_beta_outer_loop)

    def test_longest_extent_and_center_match_beta_convention(self) -> None:
        vertices = np.asarray(
            [
                [10.0, 20.0, 30.0],
                [14.0, 22.0, 31.0],
                [12.0, 20.5, 30.25],
            ],
            dtype=np.float64,
        )
        transform = PaperNormalization.from_vertices(vertices)
        normalized = transform.to_normalized(vertices)
        extent = normalized.max(axis=0) - normalized.min(axis=0)
        center = 0.5 * (normalized.max(axis=0) + normalized.min(axis=0))
        self.assertAlmostEqual(float(extent.max()), 0.95, places=14)
        np.testing.assert_allclose(center, (0.5, 0.5, 0.5), atol=1.0e-15)
        np.testing.assert_allclose(
            transform.to_original(normalized),
            vertices,
            atol=1.0e-14,
        )

    def test_support_band_is_really_projected_to_one_plane(self) -> None:
        vertices, _ = _sloped_box()
        transform = PaperNormalization.from_vertices(vertices)
        normalized = transform.to_normalized(vertices)
        flattening = SupportFlattening.from_vertices(
            normalized,
            (0.0, -2.0, 0.0),
            0.06,
            1.0e-4,
        )
        self.assertEqual(flattening.vertex_indices, (0, 1, 4, 5))
        self.assertEqual(flattening.minimum_stability_inset, 1.0e-4)
        projected = flattening.project(normalized)
        heights = -projected[list(flattening.vertex_indices)] @ np.asarray(
            flattening.gravity
        )
        np.testing.assert_allclose(
            heights,
            flattening.support_height,
            atol=1.0e-15,
        )


@unittest.skipUnless(_native_available(), "optional native extension is not built")
class PaperCompatibilityPipelineTests(unittest.TestCase):
    def test_native_material_shell_audit_measures_clearance_and_containment(
        self,
    ) -> None:
        import _make_it_stand as native

        outer_vertices, faces = _sloped_box()
        center = outer_vertices.mean(axis=0)
        inner_vertices = center + 0.5 * (outer_vertices - center)
        audit = native.audit_material_shells(
            outer_vertices,
            faces,
            inner_vertices,
            faces[:, [0, 2, 1]],
            0.01,
        )
        self.assertTrue(audit["wall_clearance"]["completed"])
        self.assertTrue(audit["wall_clearance"]["intersection_free"])
        self.assertTrue(audit["wall_clearance"]["meets_requested_thickness"])
        self.assertGreater(
            audit["wall_clearance"]["measured_minimum_clearance"],
            0.01,
        )
        self.assertTrue(audit["containment"]["completed"])
        self.assertTrue(audit["containment"]["all_cavities_inside_outer"])

    def test_normalized_optimization_restores_and_flattens_output(self) -> None:
        vertices, faces = _sloped_box()
        model = PaperCompatibilityModel(
            vertices,
            faces,
            config=ModelConfig(path="<memory>", voxel_resolution=8),
        )
        normalized = model.normalized_vertices
        extent = normalized.max(axis=0) - normalized.min(axis=0)
        center = 0.5 * (normalized.max(axis=0) + normalized.min(axis=0))
        self.assertAlmostEqual(float(extent.max()), 0.95, places=14)
        np.testing.assert_allclose(center, (0.5, 0.5, 0.5), atol=1.0e-15)
        self.assertTrue(model.native.config.fixed_unit_cube_voxel_grid)

        model.set_gravity((0.0, -1.0, 0.0))
        model.set_support_from_mesh(
            0.06,
            0.0,
            beta_minimum_inset=True,
        )
        model.add_handle((1.0, 11.0, 0.5))
        result = model.optimize_released_beta(
            maximum_carving_blocks=1,
            maximum_deformation_iterations_per_block=2,
            enable_carving=False,
            enable_deformation=True,
            released_beta_outer_loop=False,
        )

        self.assertTrue(result.flatten_applied)
        self.assertEqual(result.flattenings[0].minimum_stability_inset, 1.0e-4)
        support = result.flattenings[0]
        expected_original_height = (
            support.support_height - 0.5
        ) / model.normalization.scale + model.normalization.source_center[1]
        restored = result.outer_vertices
        indices = list(support.vertex_indices)
        np.testing.assert_allclose(
            restored[indices, 1],
            expected_original_height,
            atol=1.0e-12,
        )
        self.assertTrue(
            result.post_flatten_outer_boundary_audit["strict_boundary_valid"]
        )
        self.assertFalse(result.manufacturing_ready)
        self.assertIn(
            "post_export_flattening_requires_manufacturing_reaudit",
            result.manufacturing_failures,
        )
        self.assertNotIn(
            "source_hypothetical_flattening_not_applied",
            result.manufacturing_failures,
        )
        native_configuration = result.native_result.diagnostics["configuration"]
        native_metrics = result.native_result.diagnostics["r10_metrics"]
        self.assertEqual(native_configuration["initial_mu"], 0.75)
        self.assertEqual(native_configuration["shape_scale_lambda"], 20.0)
        self.assertTrue(native_configuration["fixed_unit_cube_voxel_grid"])
        self.assertEqual(native_metrics["grid_dimensions"], [8, 8, 8])
        np.testing.assert_allclose(
            native_metrics["voxel_size"],
            [0.125, 0.125, 0.125],
            atol=0.0,
        )
        self.assertTrue(
            native_configuration["freeze_targets_after_initial_refresh"]
        )
        self.assertFalse(native_configuration["released_beta_outer_loop"])
        self.assertTrue(
            native_configuration["line_search"]["released_beta_semantics"]
        )

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            obj = root / "restored.obj"
            stl = root / "restored.stl"
            diagnostics = root / "paper_compat.json"
            result.save_outer_mesh(obj)
            result.save_outer_mesh(stl)
            result.save_diagnostics(diagnostics)
            self.assertTrue(obj.read_text(encoding="utf-8").startswith("v "))
            self.assertTrue(
                stl.read_text(encoding="utf-8").startswith(
                    "solid make_it_stand_paper_compat"
                )
            )
            reloaded = PaperCompatibilityModel(
                stl,
                config=ModelConfig(path=str(stl), voxel_resolution=8),
            )
            self.assertEqual(reloaded.source_vertices.shape[1], 3)
            self.assertEqual(reloaded.faces.shape, result.outer_faces.shape)
            reloaded_extent = (
                reloaded.normalized_vertices.max(axis=0)
                - reloaded.normalized_vertices.min(axis=0)
            )
            self.assertAlmostEqual(float(reloaded_extent.max()), 0.95, places=14)
            document = json.loads(diagnostics.read_text(encoding="utf-8"))
            self.assertEqual(
                document["coordinate_contract"]["returned_meshes"],
                "source_coordinates",
            )
            self.assertTrue(
                document["support_flattenings"][0]["projection_applied"]
            )

    def test_sealed_artifact_recomputes_moments_after_support_flattening(self) -> None:
        import _make_it_stand as native

        vertices, faces = _sloped_box()
        normalization = PaperNormalization.from_vertices(vertices)
        normalized = normalization.to_normalized(vertices)
        flattening = SupportFlattening.from_vertices(
            normalized,
            (0.0, -1.0, 0.0),
            0.06,
            1.0e-4,
        )
        sealed = SealedBalanceArtifact.create(
            outer_vertices=normalized,
            outer_faces=faces,
            inner_vertices=np.empty((0, 3), dtype=np.float64),
            inner_faces=np.empty((0, 3), dtype=np.int32),
            flattenings=[flattening],
            occupancy_hash="alpha-fixture",
            source_geometry_hash="geometry-before-flattening",
            gravity=(0.0, -1.0, 0.0),
            support_polygon=normalized[list(flattening.vertex_indices)],
            evidence_class="EXACT_CACHE",
            native_module=native,
            safe_region_gamma_radians=0.0,
        )
        self.assertTrue(sealed.ready)
        projected = sealed.outer_vertices[list(flattening.vertex_indices)]
        heights = -projected @ np.asarray(flattening.gravity)
        np.testing.assert_allclose(
            heights,
            flattening.support_height,
            atol=1.0e-15,
        )
        self.assertNotEqual(
            sealed.geometry_sha256,
            "geometry-before-flattening",
        )
        self.assertGreater(sealed.material_mass, 0.0)
        support_heights = -sealed.support_polygon @ np.asarray(sealed.gravity)
        np.testing.assert_allclose(
            support_heights,
            flattening.support_height,
            atol=1.0e-15,
        )
        self.assertEqual(sealed.safe_region_state, "SAFE")
        self.assertEqual(sealed.safe_region_residual, 0.0)
        self.assertEqual(sealed.mass_properties_sha256, sealed.pybullet_entry()[
            "mass_properties_sha256"
        ])

        with tempfile.TemporaryDirectory() as directory:
            manifest = sealed.write_directory(Path(directory))
            payload = json.loads(manifest.read_text(encoding="utf-8"))
            self.assertEqual(payload["schema_version"], 1)
            self.assertEqual(payload["status"], "SEALED")
            self.assertTrue((Path(directory) / "outer-flattened.obj").is_file())
            entry = load_sealed_pybullet_entry(manifest)
            self.assertEqual(entry["occupancy_hash"], "alpha-fixture")
            self.assertEqual(
                Path(entry["collision_mesh"]),
                Path(directory) / "outer-flattened.obj",
            )
            (Path(directory) / "outer-flattened.obj").write_text(
                "tampered\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "collision mesh hash"):
                load_sealed_pybullet_entry(manifest)

    def test_sealed_artifact_remeshes_a_folded_flat_support_patch(self) -> None:
        import _make_it_stand as native

        vertices = np.asarray(
            [
                [-1.0, -1.0, 0.0], [1.0, -1.0, 0.0],
                [1.0, 1.0, 0.0], [-1.0, 1.0, 0.0],
                [-0.5, -0.5, 0.1], [0.5, -0.5, 0.1],
                [0.5, 0.5, 0.1], [-0.5, 0.5, 0.1],
                [-0.5, -0.5, 1.0], [0.5, -0.5, 1.0],
                [0.5, 0.5, 1.0], [-0.5, 0.5, 1.0],
            ],
            dtype=np.float64,
        )
        faces: list[list[int]] = [[0, 2, 1], [0, 3, 2]]
        for lower, upper in ((0, 4), (4, 8)):
            for index in range(4):
                following = (index + 1) % 4
                faces.extend(
                    [
                        [lower + index, lower + following, upper + following],
                        [lower + index, upper + following, upper + index],
                    ]
                )
        faces.extend([[8, 9, 10], [8, 10, 11]])
        face_array = np.asarray(faces, dtype=np.int32)
        flattening = SupportFlattening.from_vertices(
            vertices,
            (0.0, 0.0, -1.0),
            0.2,
        )
        folded = flattening.project(vertices)
        self.assertFalse(
            native.audit_triangle_boundary(folded, face_array)[
                "strict_boundary_valid"
            ]
        )

        sealed = SealedBalanceArtifact.create(
            outer_vertices=vertices,
            outer_faces=face_array,
            inner_vertices=np.empty((0, 3), dtype=np.float64),
            inner_faces=np.empty((0, 3), dtype=np.int32),
            flattenings=[flattening],
            occupancy_hash="alpha-folded-support",
            source_geometry_hash="geometry-folded-support",
            gravity=(0.0, 0.0, -1.0),
            support_polygon=vertices[4:8],
            evidence_class="RECONSTRUCTED_PARAMETER",
            native_module=native,
            safe_region_gamma_radians=0.0,
        )

        self.assertTrue(sealed.ready)
        self.assertTrue(sealed.outer_boundary_audit["strict_boundary_valid"])
        self.assertEqual(
            sealed.outer_flattening_repair_audit["status"],
            "APPLIED",
        )
        self.assertEqual(
            sealed.outer_flattening_repair_audit["patches"][0][
                "boundary_loop_lengths"
            ],
            [4],
        )
        self.assertEqual(sealed.outer_faces.shape, (12, 3))

    def test_invalid_flattened_boundary_cannot_feed_pybullet(self) -> None:
        import _make_it_stand as native

        vertices, faces = _sloped_box()
        invalid_faces = faces[:-1]
        sealed = SealedBalanceArtifact.create(
            outer_vertices=vertices,
            outer_faces=invalid_faces,
            inner_vertices=np.empty((0, 3), dtype=np.float64),
            inner_faces=np.empty((0, 3), dtype=np.int32),
            flattenings=[],
            occupancy_hash="alpha-invalid",
            source_geometry_hash="geometry-invalid",
            gravity=(0.0, -1.0, 0.0),
            support_polygon=vertices[[0, 1, 4, 5]],
            evidence_class="RECONSTRUCTED_ASSET",
            native_module=native,
        )
        self.assertFalse(sealed.ready)
        self.assertIn("outer_boundary_invalid", sealed.failures)
        with self.assertRaisesRegex(RuntimeError, "not sealed"):
            sealed.pybullet_entry()

    def test_sealed_artifact_compacts_unreferenced_cavity_grid_vertices(self) -> None:
        import _make_it_stand as native

        vertices, faces = _sloped_box()
        center = vertices.mean(axis=0)
        inner = center + 0.25 * (vertices - center)
        inner = np.vstack((inner, np.asarray([[99.0, 99.0, 99.0]])))
        sealed = SealedBalanceArtifact.create(
            outer_vertices=vertices,
            outer_faces=faces,
            inner_vertices=inner,
            inner_faces=faces[:, [0, 2, 1]],
            flattenings=[],
            occupancy_hash="alpha-with-grid-corners",
            source_geometry_hash="geometry-with-grid-corners",
            gravity=(0.0, -1.0, 0.0),
            support_polygon=vertices[[0, 1, 4, 5]],
            evidence_class="EXACT_CACHE",
            native_module=native,
        )
        self.assertTrue(sealed.ready)
        self.assertEqual(sealed.inner_vertices.shape, (8, 3))
        self.assertEqual(
            sealed.inner_boundary_audit["unreferenced_vertex_count"],
            0,
        )

    def test_sealed_artifact_repairs_checkerboard_cavity_by_refill(self) -> None:
        import _make_it_stand as native

        dimensions = (2, 2, 1)
        corner_dimensions = tuple(value + 1 for value in dimensions)
        grid_corners = np.asarray(
            [
                (float(x), float(y), float(z))
                for z in range(corner_dimensions[2])
                for y in range(corner_dimensions[1])
                for x in range(corner_dimensions[0])
            ],
            dtype=np.float64,
        )
        checkerboard_void = {(0, 0, 0), (1, 1, 0)}
        raw_inner_faces = _structured_cavity_boundary(
            dimensions,
            checkerboard_void,
        )
        raw_void_audit = native.audit_triangle_boundary(
            grid_corners,
            raw_inner_faces[:, [0, 2, 1]],
        )
        self.assertFalse(raw_void_audit["strict_boundary_valid"])
        self.assertEqual(raw_void_audit["non_manifold_edge_count"], 1)

        outer_vertices = np.asarray(
            [
                [-1.0, -1.0, -1.0], [3.0, -1.0, -1.0],
                [3.0, 3.0, -1.0], [-1.0, 3.0, -1.0],
                [-1.0, -1.0, 2.0], [3.0, -1.0, 2.0],
                [3.0, 3.0, 2.0], [-1.0, 3.0, 2.0],
            ],
            dtype=np.float64,
        )
        outer_faces = _sloped_box()[1]
        source_hash = _structured_occupancy_hash(dimensions, checkerboard_void)
        sealed = SealedBalanceArtifact.create(
            outer_vertices=outer_vertices,
            outer_faces=outer_faces,
            inner_vertices=grid_corners,
            inner_faces=raw_inner_faces,
            flattenings=[],
            occupancy_hash=source_hash,
            source_geometry_hash="checkerboard-cavity",
            gravity=(0.0, -1.0, 0.0),
            support_polygon=outer_vertices[[0, 1, 4, 5]],
            evidence_class="RECONSTRUCTED_PARAMETER",
            native_module=native,
            structured_grid_dimensions=dimensions,
            repair_nonmanifold_cavity=True,
        )

        self.assertTrue(sealed.ready)
        self.assertNotEqual(sealed.occupancy_hash, source_hash)
        self.assertEqual(sealed.source_occupancy_hash, source_hash)
        self.assertEqual(
            sealed.occupancy_repair_audit["algorithm"],
            "deterministic_nonincreasing_refill_v2",
        )
        self.assertEqual(
            len(sealed.occupancy_repair_audit["refilled_cells_xyz"]),
            1,
        )
        self.assertTrue(sealed.inner_boundary_audit["strict_boundary_valid"])


if __name__ == "__main__":
    unittest.main()
