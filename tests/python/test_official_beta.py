"""Released-beta archive parsing, coordinate mapping, and opt-in gate tests."""

from __future__ import annotations

import importlib.util
import os
from pathlib import Path
import tempfile
import unittest

import numpy as np

from make_it_stand.official_beta import (
    ReleasedBetaArchive,
    normalize_and_map_beta_bbw,
    parse_released_beta_mis,
)
from make_it_stand.model import RunConfig
from make_it_stand.paper_compat import paper_run_config


def _discover_official_root() -> Path | None:
    configured = os.environ.get("MAKE_IT_STAND_BETA_ROOT")
    candidates = []
    if configured:
        candidates.append(Path(configured).expanduser())
    repository = Path(__file__).resolve().parents[2]
    candidates.append(
        repository.parent
        / "reference-only"
        / "make-it-stand"
        / "beta-extracted"
        / "make-it-stand_src_beta"
    )
    for candidate in candidates:
        if candidate.is_dir():
            return candidate.resolve()
    return None


OFFICIAL_ROOT = _discover_official_root()


def _native_available() -> bool:
    for name in ("_make_it_stand", "make_it_stand._make_it_stand"):
        try:
            if importlib.util.find_spec(name) is not None:
                return True
        except (ImportError, ModuleNotFoundError, ValueError):
            continue
    return False


class ReleasedBetaParserTests(unittest.TestCase):
    def test_mis_is_strict_and_preserves_handle_order(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "fixture.mis"
            path.write_text(
                "mesh fixture.obj\n"
                "mesh_rotation 1 0 0 0\n"
                "flatten_threshold 0.005\n"
                "vertex_handle 0.2 0.3 0.4\n"
                "vertex_handle 0.7 0.6 0.5\n",
                encoding="utf-8",
            )
            result = parse_released_beta_mis(path)
            self.assertEqual(result.mesh, "fixture.obj")
            self.assertEqual(result.mesh_rotation_wxyz, (1.0, 0.0, 0.0, 0.0))
            self.assertEqual(result.flatten_threshold, 0.005)
            self.assertEqual(
                result.vertex_handles,
                ((0.2, 0.3, 0.4), (0.7, 0.6, 0.5)),
            )
            path.write_text("mesh fixture.obj\nunknown 1\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "unsupported compatibility tag"):
                parse_released_beta_mis(path)

    def test_bbw_rows_are_normalized_then_mapped_by_coordinate(self) -> None:
        resolution = 3
        ids = np.full(resolution ** 3, -1, dtype=np.int32)
        ids[(0 * resolution + 1) * resolution + 2] = 0
        occupied = ids.reshape((resolution,) * 3) >= 0
        active = np.zeros((resolution + 1,) * 3, dtype=np.bool_)
        for dx in (0, 1):
            for dy in (0, 1):
                for dz in (0, 1):
                    active[
                        dx:dx + resolution,
                        dy:dy + resolution,
                        dz:dz + resolution,
                    ] |= occupied
        beta_coordinates = np.argwhere(active)
        raw = np.asarray(
            [[float(index + 1), 1.0] for index in range(len(beta_coordinates))],
            dtype=np.float64,
        )
        normalized, mapped = normalize_and_map_beta_bbw(ids, raw, resolution)
        np.testing.assert_allclose(normalized.sum(axis=1), 1.0, atol=1.0e-15)
        np.testing.assert_allclose(mapped.sum(axis=1), 1.0, atol=1.0e-15)

        side = resolution + 1
        keys = (
            beta_coordinates[:, 0]
            + side * (beta_coordinates[:, 1] + side * beta_coordinates[:, 2])
        )
        expected_order = np.argsort(keys, kind="stable")
        np.testing.assert_array_equal(mapped, normalized[expected_order])
        self.assertFalse(np.array_equal(expected_order, np.arange(8)))

        bad = raw.copy()
        bad[0] = 0.0
        with self.assertRaisesRegex(ValueError, "positive finite sum"):
            normalize_and_map_beta_bbw(ids, bad, resolution)
        bad = raw.copy()
        bad[0, 0] = -1.0
        with self.assertRaisesRegex(ValueError, "finite, nonnegative"):
            normalize_and_map_beta_bbw(ids, bad, resolution)


@unittest.skipUnless(
    _native_available() and OFFICIAL_ROOT is not None,
    "native extension or reference-only released beta is unavailable",
)
class OfficialArchiveGateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        import _make_it_stand as native

        cls.native = native

    def test_spheres_cache_contract_and_explicit_gate(self) -> None:
        assert OFFICIAL_ROOT is not None
        archive = ReleasedBetaArchive.load(
            OFFICIAL_ROOT / "data/configs/spheres.mis",
            beta_root=OFFICIAL_ROOT,
        )
        self.assertEqual(archive.voxel_resolution, 80)
        self.assertEqual(archive.occupied_cell_count, 83623)
        self.assertEqual(archive.active_node_count, 92809)
        self.assertEqual(archive.beta_weights_z_fast.shape, (92809, 4))
        np.testing.assert_allclose(
            archive.beta_normalized_weights.sum(axis=1),
            1.0,
            atol=1.0e-12,
        )

        boundary = self.native.audit_triangle_boundary(
            archive.normalized_vertices,
            archive.faces,
        )
        self.assertTrue(boundary["strict_boundary_valid"])
        self.assertEqual(boundary["self_intersection_pair_count"], 0)

        model = archive.build_model(native_module=self.native)
        self.assertEqual(model.native.config.voxel_resolution, 80)
        self.assertTrue(model.native.config.fixed_unit_cube_voxel_grid)
        paper = paper_run_config()
        self.assertFalse(paper.freeze_targets_after_initial_refresh)
        self.assertFalse(paper.released_beta_semantics)
        self.assertFalse(paper.released_beta_outer_loop)

    def test_horse_true_crossing_remains_fail_closed(self) -> None:
        assert OFFICIAL_ROOT is not None
        archive = ReleasedBetaArchive.load(
            OFFICIAL_ROOT / "data/configs/horse.mis",
            beta_root=OFFICIAL_ROOT,
        )
        with self.assertRaisesRegex(ValueError, "compatibility boundary gate"):
            archive.build_model(native_module=self.native)

    def test_horse_occupancy_reconstruction_is_explicit_and_buildable(self) -> None:
        assert OFFICIAL_ROOT is not None
        archive = ReleasedBetaArchive.load(
            OFFICIAL_ROOT / "data/configs/horse.mis",
            beta_root=OFFICIAL_ROOT,
        )
        replacement = (
            Path(__file__).resolve().parents[2]
            / "experiments/paper_sim_v1/reconstructed_assets"
            / "horse_100_occupancy_outer.obj"
        )
        self.assertTrue(replacement.is_file())
        reconstructed = archive.with_reconstructed_normalized_surface(replacement)
        self.assertIsNot(reconstructed, archive)
        np.testing.assert_array_equal(
            reconstructed.source_vertices,
            archive.source_vertices,
        )
        self.assertNotEqual(reconstructed.faces.shape, archive.faces.shape)
        model = reconstructed.build_model(native_module=self.native)
        result = model.optimize(
            paper_run_config(
                maximum_carving_blocks=1,
                maximum_deformation_iterations_per_block=1,
            )
        )
        self.assertEqual(
            result.diagnostics["configuration"]["voxel_resolution"],
            100,
        )

    def test_official_complete_archive_build_matrix(self) -> None:
        assert OFFICIAL_ROOT is not None
        expected = {
            "armadillo1": (80, 21636, 27607, 8, True),
            "armadillo2": (80, 20340, 25927, 4, True),
            "dinosaur": (100, 37549, 46095, 5, True),
            "horse": (100, 26597, 33646, 4, False),
            "spheres": (80, 83623, 92809, 4, True),
        }
        for name, (resolution, cells, nodes, columns, accepted) in expected.items():
            with self.subTest(name=name):
                archive = ReleasedBetaArchive.load(
                    OFFICIAL_ROOT / f"data/configs/{name}.mis",
                    beta_root=OFFICIAL_ROOT,
                )
                self.assertEqual(archive.voxel_resolution, resolution)
                self.assertEqual(archive.occupied_cell_count, cells)
                self.assertEqual(archive.active_node_count, nodes)
                self.assertEqual(archive.beta_weights_z_fast.shape, (nodes, columns))
                if accepted:
                    archive.build_model(native_module=self.native)
                else:
                    with self.assertRaisesRegex(
                        ValueError,
                        "compatibility boundary gate",
                    ):
                        archive.build_model(native_module=self.native)

    @unittest.skipUnless(
        os.environ.get("MAKE_IT_STAND_RUN_OFFICIAL_BETA_SMOKE") == "1",
        "set MAKE_IT_STAND_RUN_OFFICIAL_BETA_SMOKE=1 for the N80 optimizer smoke",
    )
    def test_spheres_enters_one_bounded_optimization(self) -> None:
        assert OFFICIAL_ROOT is not None
        archive = ReleasedBetaArchive.load(
            OFFICIAL_ROOT / "data/configs/spheres.mis",
            beta_root=OFFICIAL_ROOT,
        )
        model = archive.build_model(native_module=self.native)
        result = model.optimize_released_beta(
            maximum_carving_blocks=1,
            maximum_deformation_iterations_per_block=1,
            enable_carving=False,
            enable_deformation=True,
            enable_laplacian=False,
            released_beta_outer_loop=False,
        )
        configuration = result.diagnostics["configuration"]
        metrics = result.diagnostics["r10_metrics"]
        self.assertTrue(configuration["released_beta_archive_import"])
        self.assertEqual(metrics["grid_dimensions"], [80, 80, 80])
        self.assertGreaterEqual(len(result.diagnostics["r10_events"]), 1)

    @unittest.skipUnless(
        os.environ.get("MAKE_IT_STAND_RUN_OFFICIAL_BETA_FULL_SMOKE") == "1",
        "set MAKE_IT_STAND_RUN_OFFICIAL_BETA_FULL_SMOKE=1 for the full N80 smoke",
    )
    def test_spheres_enters_bounded_full_method(self) -> None:
        assert OFFICIAL_ROOT is not None
        archive = ReleasedBetaArchive.load(
            OFFICIAL_ROOT / "data/configs/spheres.mis",
            beta_root=OFFICIAL_ROOT,
        )
        model = archive.build_model(native_module=self.native)
        result = model.optimize_released_beta(
            maximum_carving_blocks=1,
            maximum_deformation_iterations_per_block=1,
            enable_carving=True,
            enable_deformation=True,
            enable_laplacian=True,
            released_beta_outer_loop=True,
        )
        diagnostics = result.diagnostics
        configuration = diagnostics["configuration"]
        self.assertTrue(configuration["released_beta_archive_import"])
        self.assertTrue(configuration["enable_carving"])
        self.assertTrue(configuration["enable_deformation"])
        self.assertTrue(configuration["enable_laplacian"])
        self.assertTrue(configuration["released_beta_outer_loop"])
        self.assertEqual(configuration["effective_shape_scale_lambda"], 20)
        self.assertEqual(
            diagnostics["r10_metrics"]["grid_dimensions"],
            [80, 80, 80],
        )
        self.assertGreaterEqual(len(diagnostics["carving_blocks"]), 1)
        self.assertGreaterEqual(len(diagnostics["carving_iterations"]), 1)
        self.assertEqual(len(diagnostics["deformation_iterations"]), 1)
        deformation = diagnostics["deformation_iterations"][0]
        self.assertGreaterEqual(deformation["carving_calls"], 1)
        self.assertGreaterEqual(deformation["trial_balance_calls"], 1)
        event_kinds = [event["kind"] for event in diagnostics["r10_events"]]
        self.assertIn("carving_iteration", event_kinds)
        self.assertIn("accepted_deformation", event_kinds)
        self.assertEqual(diagnostics["status"], "maximum_deformation_iterations")
        self.assertFalse(diagnostics["succeeded"])


if __name__ == "__main__":
    unittest.main()
