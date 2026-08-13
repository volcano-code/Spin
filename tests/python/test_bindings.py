"""Native pybind11 API and owning-lifetime regression tests."""

from __future__ import annotations

import importlib
import importlib.util
import math
import gc
import json
import threading
import tempfile
import time
import unittest

from pathlib import Path


def _find_native_module() -> str | None:
    """Return the importable extension name without importing it."""

    candidates = ("_make_it_stand", "make_it_stand._make_it_stand")
    for candidate in candidates:
        try:
            if importlib.util.find_spec(candidate) is not None:
                return candidate
        except (ImportError, ModuleNotFoundError, ValueError):
            # A package-local extension cannot be discovered before the package
            # exists.  The top-level candidate remains the canonical build-tree
            # layout used by these tests.
            continue
    return None


_NATIVE_MODULE_NAME = _find_native_module()


@unittest.skipIf(
    _NATIVE_MODULE_NAME is None,
    "optional _make_it_stand extension is not built",
)
class NativeBindingSmokeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        # If a module was discoverable but cannot be loaded (missing shared
        # library, ABI mismatch, unresolved symbol, ...), this import must fail
        # the test rather than silently turning a broken build into a skip.
        cls.native = importlib.import_module(_NATIVE_MODULE_NAME)  # type: ignore[arg-type]
        # Load the C++ extension first.  Some Conda NumPy builds preload an
        # older libstdc++ whose symbol set cannot satisfy the extension.
        # NumPy remains a required runtime dependency when the extension exists:
        # its ImportError intentionally fails this class rather than skipping it.
        import numpy

        cls.np = numpy

    def test_required_types_are_exported(self) -> None:
        for name in (
            "BalanceModel",
            "BalanceModelConfig",
            "BalanceRunConfig",
            "BalanceResult",
        ):
            with self.subTest(name=name):
                self.assertTrue(hasattr(self.native, name), name)

    def test_model_surface_is_present(self) -> None:
        required_methods = (
            "set_gravity",
            "set_support_polygon",
            "set_support_from_mesh",
            "add_standing_constraint_polygon",
            "add_standing_constraint_from_mesh",
            "add_handle",
            "add_suspension_constraint",
            "clear_handles",
            "clear_constraints",
            "optimize",
        )
        for name in required_methods:
            with self.subTest(name=name):
                self.assertTrue(hasattr(self.native.BalanceModel, name), name)

    def test_result_surface_is_present(self) -> None:
        required_properties = (
            "status",
            "succeeded",
            "outer_vertices",
            "outer_faces",
            "inner_vertices",
            "inner_faces",
            "manufacturing_ready",
            "manufacturing_failures",
            "manufacturing_report_json",
            "diagnostics_json",
            "performance_report_json",
        )
        required_methods = (
            "save_outer_mesh",
            "save_inner_mesh",
            "save_printable_mesh",
            "save_diagnostics",
        )
        for name in required_properties + required_methods:
            with self.subTest(name=name):
                self.assertTrue(hasattr(self.native.BalanceResult, name), name)

    def test_model_config_defaults_and_fields(self) -> None:
        config = self.native.BalanceModelConfig()
        expected = {
            "voxel_resolution": 50,
            "minimum_wall_thickness": 0.0,
            "bbw_max_iterations": 200,
            "bbw_tolerance": 1.0e-8,
        }
        for name, value in expected.items():
            with self.subTest(name=name):
                self.assertTrue(hasattr(config, name), name)
                self.assertEqual(getattr(config, name), value)

    def test_run_config_fields_have_finite_defaults(self) -> None:
        config = self.native.BalanceRunConfig()
        fields = (
            "initial_mu",
            "maximum_carving_blocks",
            "maximum_deformation_iterations_per_block",
            "gradient_tolerance",
            "shape_scale_lambda",
            "enable_carving",
            "enable_deformation",
            "enable_laplacian",
            "density",
            "minimum_positive_mass",
            "minimum_positive_cell_mass",
            "minimum_carving_removals_to_continue",
            "initial_step",
            "minimum_step",
            "reduction_ratio",
            "required_relative_decrease",
            "mu_decrement",
        )
        for name in fields:
            with self.subTest(name=name):
                self.assertTrue(hasattr(config, name), name)
                value = getattr(config, name)
                if name.startswith("enable_"):
                    self.assertIsInstance(value, bool)
                else:
                    self.assertIsInstance(value, (int, float))
                    self.assertNotIsInstance(value, bool)
                    self.assertTrue(math.isfinite(value), name)

    @classmethod
    def _box_mesh(cls):
        vertices = cls.np.asarray(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [1.0, 1.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
                [1.0, 0.0, 1.0],
                [1.0, 1.0, 1.0],
                [0.0, 1.0, 1.0],
            ],
            dtype=cls.np.float64,
        )
        faces = cls.np.asarray(
            [
                [0, 2, 1], [0, 3, 2], [4, 5, 6], [4, 6, 7],
                [0, 1, 5], [0, 5, 4], [3, 7, 6], [3, 6, 2],
                [0, 4, 7], [0, 7, 3], [1, 2, 6], [1, 6, 5],
            ],
            dtype=cls.np.int32,
        )
        return vertices, faces

    def _stable_box_result(self):
        vertices, faces = self._box_mesh()
        model_config = self.native.BalanceModelConfig()
        model_config.voxel_resolution = 8
        model = self.native.BalanceModel(vertices, faces, model_config)
        model.set_gravity(self.np.asarray([0.0, -1.0, 0.0]))
        model.set_support_polygon(
            self.np.asarray(
                [
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [1.0, 0.0, 1.0],
                    [0.0, 0.0, 1.0],
                ],
                dtype=self.np.float64,
            ),
            0.0,
        )
        run_config = self.native.BalanceRunConfig()
        run_config.maximum_carving_blocks = 2
        run_config.maximum_deformation_iterations_per_block = 4
        return model, model.optimize(run_config)

    def test_invalid_numpy_inputs_fail_before_optimization(self) -> None:
        vertices, faces = self._box_mesh()
        with self.assertRaises((ValueError, RuntimeError)):
            self.native.BalanceModel(vertices[:, :2], faces)
        damaged = vertices.copy()
        damaged[0, 0] = math.nan
        with self.assertRaises((ValueError, RuntimeError)):
            self.native.BalanceModel(damaged, faces)
        bad_faces = faces.astype(self.np.int64)
        bad_faces[0, 0] = len(vertices)
        with self.assertRaises((ValueError, RuntimeError)):
            self.native.BalanceModel(vertices, bad_faces)

        with self.assertRaises((TypeError, ValueError, RuntimeError)):
            self.native.BalanceModel(
                vertices.astype(self.np.complex128),
                faces,
            )
        with self.assertRaises((TypeError, ValueError, RuntimeError)):
            self.native.BalanceModel(vertices, faces.astype(self.np.float64))

        wrapped_unsigned = faces.astype(self.np.uint64)
        wrapped_unsigned[0, 0] = self.np.iinfo(self.np.uint64).max
        with self.assertRaises((TypeError, ValueError, RuntimeError)):
            self.native.BalanceModel(vertices, wrapped_unsigned)

        oversized_signed = faces.astype(self.np.int64)
        oversized_signed[0, 0] = self.np.iinfo(self.np.int32).max + 1
        with self.assertRaises((TypeError, ValueError, RuntimeError)):
            self.native.BalanceModel(vertices, oversized_signed)

        # Wider integer dtypes remain valid when every value is representable.
        accepted = self.native.BalanceModel(vertices, faces.astype(self.np.uint64))
        self.assertTrue(self.np.array_equal(accepted.faces, faces))

    def test_optimize_snapshots_python_config_before_releasing_gil(self) -> None:
        vertices, faces = self._box_mesh()
        model_config = self.native.BalanceModelConfig()
        model_config.voxel_resolution = 8
        model = self.native.BalanceModel(vertices, faces, model_config)
        model.set_support_from_mesh(0.0, 1.0)
        model.add_handle(self.np.asarray([0.5, 0.5, 0.5]))

        run_config = self.native.BalanceRunConfig()
        run_config.maximum_carving_blocks = 5
        run_config.maximum_deformation_iterations_per_block = 100
        mutator_started = threading.Event()
        mutation_times: list[float] = []

        def mutate_after_native_entry() -> None:
            mutator_started.set()
            time.sleep(0.005)
            # Invalid values make any post-release read of the Python-owned
            # object observable.  The C++ run must only see the prior snapshot.
            run_config.maximum_carving_blocks = 0
            run_config.maximum_deformation_iterations_per_block = 0
            mutation_times.append(time.monotonic())

        mutator = threading.Thread(target=mutate_after_native_entry)
        mutator.start()
        self.assertTrue(mutator_started.wait(timeout=1.0))
        result = model.optimize(run_config)
        returned_at = time.monotonic()
        mutator.join(timeout=1.0)
        self.assertFalse(mutator.is_alive())
        self.assertTrue(mutation_times)
        self.assertLess(mutation_times[0], returned_at)
        configuration = json.loads(result.diagnostics_json)["configuration"]
        self.assertEqual(configuration["maximum_carving_blocks"], 5)
        self.assertEqual(
            configuration["maximum_deformation_iterations_per_block"],
            100,
        )

    def test_real_cpp_run_is_stable_and_deterministic(self) -> None:
        first_model, first = self._stable_box_result()
        second_model, second = self._stable_box_result()
        self.assertTrue(first.succeeded)
        self.assertEqual(first.status, "stable_after_carving")
        self.assertTrue(first.manufacturing_ready)
        self.assertEqual(first.diagnostics_json, second.diagnostics_json)
        diagnostic = json.loads(first.diagnostics_json)
        self.assertEqual(
            diagnostic["schema"],
            "make-it-stand-r9-balance-result-v1",
        )
        self.assertEqual(diagnostic["schema_version"], 1)
        self.assertEqual(diagnostic["status"], first.status)
        self.assertTrue(diagnostic["succeeded"])
        self.assertGreater(diagnostic["mass_properties"]["mass"], 0.0)
        self.assertEqual(len(diagnostic["center_of_mass"]), 3)
        self.assertIn("r10_metrics", diagnostic)
        self.assertIn("r10_events", diagnostic)
        self.assertEqual(diagnostic["r10_events"][0]["kind"], "initial")
        self.assertEqual(diagnostic["r10_metrics"]["first_safe_event"], 0)
        self.assertGreaterEqual(
            diagnostic["r10_metrics"]["final_cavity_volume_fraction"],
            0.0,
        )
        self.assertGreaterEqual(
            diagnostic["r10_metrics"]["maximum_user_handle_scale_deviation"],
            0.0,
        )
        performance = json.loads(first.performance_report_json)
        self.assertEqual(
            performance["schema"],
            "make-it-stand-r10-performance-report-v1",
        )
        self.assertEqual(performance["clock"], "steady_clock")
        self.assertEqual(
            performance["scope"],
            "BalanceModel.optimize_after_config_validation_"
            "through_diagnostics_serialization",
        )
        self.assertGreaterEqual(performance["total_seconds"], 0.0)
        for value in performance["phases"].values():
            self.assertTrue(math.isfinite(value))
            self.assertGreaterEqual(value, 0.0)
        del first_model, second_model

    def test_result_owns_state_and_arrays_are_copies(self) -> None:
        model, result = self._stable_box_result()
        del model
        gc.collect()
        baseline = result.outer_vertices
        changed = result.outer_vertices
        changed[0, 0] += 1000.0
        self.assertTrue(self.np.array_equal(result.outer_vertices, baseline))
        self.assertTrue(result.manufacturing_ready)

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            result.save_outer_mesh(str(root / "outer.obj"))
            result.save_inner_mesh(str(root / "inner.obj"))
            result.save_printable_mesh(str(root / "printable.stl"))
            result.save_diagnostics(str(root / "diagnostics.json"))
            self.assertTrue((root / "outer.obj").read_bytes())
            self.assertTrue((root / "printable.stl").read_bytes())
            self.assertEqual(
                json.loads((root / "diagnostics.json").read_text()),
                json.loads(result.diagnostics_json),
            )

    def test_configuration_and_support_contracts_fail_closed(self) -> None:
        vertices, faces = self._box_mesh()
        model = self.native.BalanceModel(vertices, faces)
        with self.assertRaises((ValueError, RuntimeError)):
            model.set_gravity(self.np.asarray([0.0, 0.0, 0.0]))
        model.set_gravity(self.np.asarray([0.0, -1.0, 0.0]))
        floating = self.np.asarray(
            [[0.0, 0.2, 0.0], [1.0, 0.2, 0.0], [1.0, 0.2, 1.0]],
            dtype=self.np.float64,
        )
        with self.assertRaises((ValueError, RuntimeError)):
            model.set_support_polygon(floating, 0.0)
        model.set_support_from_mesh(0.0, 0.0)
        with self.assertRaises((ValueError, RuntimeError)):
            model.set_gravity(self.np.asarray([1.0, 0.0, 0.0]))

    def test_multiple_standing_poses_are_exposed_and_auditable(self) -> None:
        vertices, faces = self._box_mesh()
        model_config = self.native.BalanceModelConfig()
        model_config.voxel_resolution = 8
        model = self.native.BalanceModel(vertices, faces, model_config)
        model.add_standing_constraint_from_mesh(
            self.np.asarray([0.0, -1.0, 0.0]),
            0.0,
            0.0,
        )
        model.add_standing_constraint_from_mesh(
            self.np.asarray([-1.0, 0.0, 0.0]),
            0.0,
            0.0,
        )
        run_config = self.native.BalanceRunConfig()
        run_config.maximum_carving_blocks = 2
        run_config.maximum_deformation_iterations_per_block = 2
        result = model.optimize(run_config)
        diagnostic = json.loads(result.diagnostics_json)
        self.assertEqual(len(diagnostic["problem"]["standing_constraints"]), 2)
        self.assertIsNone(diagnostic["problem"]["standing_constraint"])
        self.assertEqual(len(diagnostic["targets"]["constraints"]), 2)
        metrics = diagnostic["r10_metrics"]
        self.assertGreater(metrics["active_deformation_nodes"], 0)
        self.assertGreater(metrics["fixed_deformation_nodes"], 0)
        self.assertGreater(metrics["fixed_deformation_node_fraction"], 0.0)


if __name__ == "__main__":
    unittest.main()
