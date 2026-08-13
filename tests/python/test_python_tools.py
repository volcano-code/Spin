"""Unit tests for the Python-only R9 experiment layer.

The native optimizer is replaced by a small fake in workflow tests.  This
keeps schema, orchestration, reporting, and filesystem contracts testable on a
machine that has not built the optional pybind11 extension.
"""

from __future__ import annotations

import csv
import copy
import dataclasses
import importlib
import json
import math
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock


def _valid_config_dict() -> dict[str, object]:
    """A deliberately complete schema-v1 document.

    Keeping every optional numeric field explicit lets the rejection tests
    exercise the exact path they intend instead of accidentally testing a
    default or an earlier missing-field error.
    """

    return {
        "schema": "make-it-stand-balance-config-v1",
        "schema_version": 1,
        "model": {
            "path": "fixtures/model.obj",
            "voxel_resolution": 24,
            "minimum_wall_thickness": 0.01,
            "bbw_max_iterations": 40,
            "bbw_tolerance": 1.0e-7,
        },
        "gravity": [0.0, -9.81, 0.0],
        "support": {
            "kind": "polygon",
            "gamma_radians": 0.1,
            "polygon": [
                [-1.0, 0.0, -1.0],
                [1.0, 0.0, -1.0],
                [0.0, 0.0, 1.0],
            ],
        },
        "handles": [
            {
                "rest_pivot": [0.0, 1.0, 0.0],
                "fixed_rotation": [
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 1.0],
                ],
            }
        ],
        "suspensions": [
            {
                "attachment_point": [0.0, 2.0, 0.0],
                "allowed_deviation_radians": 0.2,
            }
        ],
        "run": {
            "initial_mu": 0.2,
            "maximum_carving_blocks": 3,
            "maximum_deformation_iterations_per_block": 4,
            "gradient_tolerance": 1.0e-8,
            "shape_scale_lambda": 1.0,
            "enable_carving": True,
            "enable_deformation": True,
            "enable_laplacian": True,
            "density": 1.0,
            "minimum_positive_mass": 0.0,
            "minimum_positive_cell_mass": 0.0,
            "minimum_carving_removals_to_continue": 1,
            "initial_step": 1.0,
            "minimum_step": 0.4,
            "reduction_ratio": 0.8,
            "required_relative_decrease": 0.01,
            "mu_decrement": 0.05,
        },
        "outputs": {
            "diagnostics": "diagnostics.json",
            "outer_mesh": "outer.obj",
            "inner_mesh": "inner.obj",
            "printable_mesh": "printable.stl",
        },
    }


def _assign_path(document: object, path: tuple[object, ...], value: object) -> None:
    target = document
    for component in path[:-1]:
        target = target[component]  # type: ignore[index]
    target[path[-1]] = value  # type: ignore[index]


def _diagnostics(**overrides: object) -> dict[str, object]:
    value: dict[str, object] = {
        "schema": "make-it-stand-r9-balance-result-v1",
        "schema_version": 1,
        "status": "stable_after_carving",
        "succeeded": True,
        "manufacturing": {"ready": False},
        "mass_properties": {"mass": 1.0},
        "carving_iterations": [],
        "deformation_iterations": [],
    }
    value.update(overrides)
    return value


class ConfigSchemaTests(unittest.TestCase):
    def setUp(self) -> None:
        self.api = importlib.import_module("make_it_stand")

    def test_complete_schema_v1_normalizes_to_deeply_immutable_dataclasses(self) -> None:
        config = self.api.validate_config(_valid_config_dict())

        self.assertIsInstance(config, self.api.BalanceConfig)
        self.assertTrue(dataclasses.is_dataclass(config))
        self.assertEqual(config.schema_version, 1)
        self.assertEqual(config.model.path, "fixtures/model.obj")
        self.assertEqual(config.gravity, (0.0, -9.81, 0.0))
        self.assertEqual(config.support.kind, "polygon")
        self.assertIsInstance(config.handles, tuple)
        self.assertIsInstance(config.support.polygon, tuple)
        self.assertIsInstance(config.handles[0].fixed_rotation, tuple)

        with self.assertRaises(dataclasses.FrozenInstanceError):
            config.schema_version = 2
        with self.assertRaises(dataclasses.FrozenInstanceError):
            config.model.path = "changed.obj"

    def test_defaults_are_versioned_and_not_shared_mutable_state(self) -> None:
        raw = {
            "schema": "make-it-stand-balance-config-v1",
            "schema_version": 1,
            "model": {"path": "model.obj"},
            "gravity": [0, -1, 0],
            "support": {
                "kind": "mesh",
                "gamma_radians": 0,
            },
            "outputs": {},
        }
        first = self.api.validate_config(raw)
        second = self.api.validate_config(copy.deepcopy(raw))

        self.assertEqual(first, second)
        self.assertEqual(first.model.voxel_resolution, 50)
        self.assertEqual(first.handles, ())
        self.assertEqual(first.suspensions, ())
        self.assertEqual(first.outputs.diagnostics, "diagnostics.json")
        self.assertEqual(first.support.contact_threshold, 0.0)

    def test_r10_standing_field_preserves_legacy_positional_constructor(self) -> None:
        normalized = self.api.validate_config(_valid_config_dict())
        legacy = self.api.BalanceConfig(
            normalized.model,
            normalized.gravity,
            normalized.support,
            normalized.handles,
            normalized.suspensions,
            normalized.run,
            normalized.outputs,
            normalized.schema,
            normalized.schema_version,
        )
        self.assertEqual(legacy.run, normalized.run)
        self.assertEqual(legacy.outputs, normalized.outputs)
        self.assertEqual(legacy.standing_constraints, ())

        legacy_run = self.api.RunConfig(
            0.3,
            4,
            5,
            1.0e-7,
            2.0,
            1.2,
            0.01,
            0.02,
            3,
            0.8,
            0.2,
            0.7,
            0.03,
            0.04,
        )
        self.assertEqual(legacy_run.density, 1.2)
        self.assertEqual(legacy_run.minimum_positive_mass, 0.01)
        self.assertEqual(legacy_run.mu_decrement, 0.04)
        self.assertTrue(legacy_run.enable_carving)

    def test_multiple_standings_validate_and_preserve_declared_order(self) -> None:
        raw = _valid_config_dict()
        raw["support"] = None
        raw["suspensions"] = []
        raw["standing_constraints"] = [
            {
                "kind": "polygon",
                "gravity": [0.0, -1.0, 0.0],
                "gamma_radians": 0.1,
                "polygon": [
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 0.0, 1.0],
                ],
            },
            {
                "kind": "mesh",
                "gravity": [-1.0, 0.0, 0.0],
                "gamma_radians": 0.2,
                "contact_threshold": 0.03,
            },
        ]
        config = self.api.validate_config(raw)
        self.assertEqual(
            tuple(item.gravity for item in config.standing_constraints),
            ((0.0, -1.0, 0.0), (-1.0, 0.0, 0.0)),
        )
        native = _fake_native_module()
        self.api.build_model(config, native_module=native)
        calls = native.state.models[0].calls
        self.assertEqual(calls[0][0], "set_gravity")
        self.assertEqual(calls[1][0], "add_standing_constraint_polygon")
        self.assertEqual(calls[1][2], (0.0, -1.0, 0.0))
        self.assertEqual(calls[2], (
            "add_standing_constraint_from_mesh",
            (-1.0, 0.0, 0.0),
            0.03,
            0.2,
        ))

    def test_unknown_keys_are_rejected_at_every_schema_object(self) -> None:
        paths = (
            (),
            ("model",),
            ("support",),
            ("handles", 0),
            ("suspensions", 0),
            ("run",),
            ("outputs",),
        )
        for path in paths:
            raw = _valid_config_dict()
            target = raw
            for component in path:
                target = target[component]  # type: ignore[index]
            target["unexpected"] = 1  # type: ignore[index]
            with self.subTest(path=path), self.assertRaisesRegex(
                self.api.ConfigurationError,
                r"unknown field.*unexpected",
            ):
                self.api.validate_config(raw)

    def test_schema_marker_version_required_fields_and_root_type_are_strict(self) -> None:
        cases = (
            ([], r"config: expected an object"),
            ({}, r"missing field"),
            (
                {key: value for key, value in _valid_config_dict().items() if key != "schema"},
                r"missing field.*schema",
            ),
            ({**_valid_config_dict(), "schema": "v2"}, r"config\.schema"),
            ({**_valid_config_dict(), "schema_version": 2}, r"schema_version"),
            ({**_valid_config_dict(), "schema_version": True}, r"expected an integer"),
            ({**_valid_config_dict(), "schema_version": 1.0}, r"expected an integer"),
        )
        for raw, message in cases:
            with self.subTest(raw=raw), self.assertRaisesRegex(
                self.api.ConfigurationError,
                message,
            ):
                self.api.validate_config(raw)

    def test_wrong_container_scalar_and_cross_field_types_are_rejected(self) -> None:
        changes = (
            (("model",), []),
            (("model", "path"), 7),
            (("model", "voxel_resolution"), 3.0),
            (("model", "voxel_resolution"), True),
            (("gravity",), [0, -1]),
            (("gravity", 1), False),
            (("support", "kind"), "unknown"),
            (("support", "polygon"), [[0, 0, 0], [1, 0, 0]]),
            (("handles",), {}),
            (("handles", 0), "handle"),
            (("handles", 0, "rest_pivot"), [0, 1]),
            (("handles", 0, "fixed_rotation"), [[1, 0], [0, 1]]),
            (("suspensions",), {}),
            (("suspensions", 0, "allowed_deviation_radians"), "0.2"),
            (("run", "maximum_carving_blocks"), 2.0),
            (("run", "density"), True),
            (("run", "enable_carving"), 1),
            (("outputs",), []),
            (("outputs", "diagnostics"), 4),
        )
        for path, value in changes:
            raw = _valid_config_dict()
            _assign_path(raw, path, value)
            with self.subTest(path=path, value=value), self.assertRaises(
                self.api.ConfigurationError
            ):
                self.api.validate_config(raw)

        mesh_with_polygon = _valid_config_dict()
        mesh_with_polygon["support"] = {
            "kind": "mesh",
            "gamma_radians": 0.0,
            "polygon": [[0, 0, 0], [1, 0, 0], [0, 0, 1]],
        }
        with self.assertRaises(self.api.ConfigurationError):
            self.api.validate_config(mesh_with_polygon)

    def test_all_numeric_sections_reject_non_finite_values(self) -> None:
        numeric_paths = (
            ("model", "minimum_wall_thickness"),
            ("model", "bbw_tolerance"),
            ("gravity", 0),
            ("support", "gamma_radians"),
            ("support", "polygon", 0, 0),
            ("handles", 0, "rest_pivot", 0),
            ("handles", 0, "fixed_rotation", 0, 0),
            ("suspensions", 0, "attachment_point", 0),
            ("suspensions", 0, "allowed_deviation_radians"),
            ("run", "initial_mu"),
            ("run", "gradient_tolerance"),
            ("run", "shape_scale_lambda"),
            ("run", "density"),
            ("run", "minimum_positive_mass"),
            ("run", "minimum_positive_cell_mass"),
            ("run", "initial_step"),
            ("run", "minimum_step"),
            ("run", "reduction_ratio"),
            ("run", "required_relative_decrease"),
            ("run", "mu_decrement"),
        )
        for path in numeric_paths:
            for value in (math.nan, math.inf, -math.inf):
                raw = _valid_config_dict()
                _assign_path(raw, path, value)
                with self.subTest(path=path, value=value), self.assertRaisesRegex(
                    self.api.ConfigurationError,
                    r"finite number",
                ):
                    self.api.validate_config(raw)

    def test_numeric_domains_and_output_path_containment_are_enforced(self) -> None:
        changes = (
            (("model", "voxel_resolution"), 0),
            (("model", "minimum_wall_thickness"), -0.1),
            (("model", "bbw_max_iterations"), 0),
            (("model", "bbw_tolerance"), 0.0),
            (("gravity",), [0.0, 0.0, 0.0]),
            (("support", "gamma_radians"), math.pi),
            (("suspensions", 0, "allowed_deviation_radians"), math.pi + 0.1),
            (("run", "maximum_carving_blocks"), 0),
            (("run", "maximum_deformation_iterations_per_block"), 0),
            (("run", "shape_scale_lambda"), -1.0),
            (("run", "minimum_positive_mass"), -1.0),
            (("run", "reduction_ratio"), 1.0),
            (("run", "required_relative_decrease"), -0.1),
            (("outputs", "diagnostics"), "../escape.json"),
            (("outputs", "outer_mesh"), "/absolute.obj"),
            (("outputs", "inner_mesh"), ""),
            (("outputs", "diagnostics"), None),
        )
        for path, value in changes:
            raw = _valid_config_dict()
            _assign_path(raw, path, value)
            with self.subTest(path=path, value=value), self.assertRaises(
                self.api.ConfigurationError
            ):
                self.api.validate_config(raw)

        duplicate = _valid_config_dict()
        duplicate["outputs"]["inner_mesh"] = duplicate["outputs"]["outer_mesh"]
        with self.assertRaisesRegex(self.api.ConfigurationError, "distinct"):
            self.api.validate_config(duplicate)

        nested_collision = _valid_config_dict()
        nested_collision["outputs"]["outer_mesh"] = "diagnostics.json/outer.obj"
        with self.assertRaisesRegex(self.api.ConfigurationError, "contain"):
            self.api.validate_config(nested_collision)

    def test_public_dataclass_instances_are_revalidated_before_native_use(self) -> None:
        config = self.api.validate_config(_valid_config_dict())
        native = _fake_native_module()
        bad_outputs = dataclasses.replace(
            config.outputs,
            diagnostics="../escape.json",
        )
        with self.assertRaisesRegex(self.api.ConfigurationError, "parent traversal"):
            self.api.run_config(
                dataclasses.replace(config, outputs=bad_outputs),
                native_module=native,
            )
        duplicate_outputs = dataclasses.replace(
            config.outputs,
            inner_mesh=config.outputs.outer_mesh,
        )
        with self.assertRaisesRegex(self.api.ConfigurationError, "distinct"):
            self.api.run_config(
                dataclasses.replace(config, outputs=duplicate_outputs),
                native_module=native,
            )
        invalid_run = dataclasses.replace(config.run, density=math.nan)
        with self.assertRaisesRegex(self.api.ConfigurationError, "finite"):
            self.api.run_config(
                dataclasses.replace(config, run=invalid_run),
                native_module=native,
            )
        for invalid_standings in (None, ("bad",)):
            with self.subTest(invalid_standings=invalid_standings), self.assertRaises(
                self.api.ConfigurationError
            ):
                self.api.build_model(
                    dataclasses.replace(
                        config,
                        standing_constraints=invalid_standings,
                    ),
                    native_module=native,
                )
        with self.assertRaisesRegex(self.api.ConfigurationError, "at least 3"):
            self.api.BalanceModel(
                "model.obj",
                config=dataclasses.replace(config.model, voxel_resolution=2),
                native_module=native,
            )
        model = self.api.BalanceModel("model.obj", native_module=native)
        with self.assertRaisesRegex(self.api.ConfigurationError, "greater than zero"):
            model.optimize(dataclasses.replace(config.run, density=0.0))
        self.assertEqual(native.state.models[-1].calls, [])

    def test_at_least_one_standing_or_suspension_target_is_required(self) -> None:
        raw = _valid_config_dict()
        raw["support"] = None
        raw["suspensions"] = []
        with self.assertRaisesRegex(
            self.api.ConfigurationError,
            "support, standing, or suspension",
        ):
            self.api.validate_config(raw)

    def test_load_config_uses_strict_json_and_rejects_duplicate_or_nonfinite_tokens(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "config.json"
            path.write_text(json.dumps(_valid_config_dict()), encoding="utf-8")
            loaded = self.api.load_config(path)
            self.assertEqual(
                loaded.model.path,
                str((path.parent / "fixtures/model.obj").resolve()),
            )
            self.assertEqual(
                dataclasses.replace(loaded.model, path="fixtures/model.obj"),
                self.api.validate_config(_valid_config_dict()).model,
            )

            path.write_text(
                '{"schema_version": 1, "schema_version": 1}',
                encoding="utf-8",
            )
            with self.assertRaisesRegex(self.api.DataFormatError, "duplicate JSON key"):
                self.api.load_config(path)

            path.write_text('{"schema_version": NaN}', encoding="utf-8")
            with self.assertRaisesRegex(self.api.DataFormatError, "non-finite"):
                self.api.load_config(path)

    def test_public_vector_builder_accepts_indexable_numeric_sequences(self) -> None:
        class NumericVector:
            def __init__(self, values: tuple[float, ...]):
                self.values = values

            def __len__(self) -> int:
                return len(self.values)

            def __getitem__(self, index: int) -> float:
                return self.values[index]

        native = _fake_native_module()
        model = self.api.BalanceModel("model.obj", native_module=native)
        model.set_gravity(NumericVector((0.0, -1.0, 0.0)))
        self.assertEqual(
            native.state.models[0].calls,
            [("set_gravity", (0.0, -1.0, 0.0))],
        )

        try:
            import numpy as np
        except ImportError:
            return
        model.set_gravity(np.asarray([0.0, -1.0, 0.0], dtype=np.float32))
        model.set_support_polygon(
            np.asarray(
                [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
                dtype=np.float64,
            ),
            np.float64(0.0),
        )
        self.assertEqual(native.state.models[0].calls[-1][0], "set_support_polygon")


class AtomicOutputTests(unittest.TestCase):
    def setUp(self) -> None:
        self.api = importlib.import_module("make_it_stand")
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary_directory.cleanup)
        self.root = Path(self.temporary_directory.name)

    def test_text_writer_creates_parents_and_replaces_existing_file(self) -> None:
        destination = self.root / "nested" / "report.svg"
        self.api.write_text_atomic(destination, "old")
        self.api.write_text_atomic(destination, "<svg>new</svg>\n")

        self.assertEqual(destination.read_text(encoding="utf-8"), "<svg>new</svg>\n")
        self.assertEqual(list(destination.parent.iterdir()), [destination])

    def test_json_writer_is_utf8_deterministic_and_strict(self) -> None:
        destination = self.root / "summary.json"
        value = {"z": [1, 2], "a": "雪"}

        self.api.write_json_atomic(destination, value)
        first = destination.read_bytes()
        self.api.write_json_atomic(destination, value)

        self.assertEqual(destination.read_bytes(), first)
        self.assertEqual(json.loads(first.decode("utf-8")), value)
        self.assertNotIn(b"NaN", first)

        with self.assertRaises((TypeError, ValueError)):
            self.api.write_json_atomic(destination, {"bad": math.nan})
        self.assertEqual(destination.read_bytes(), first)
        self.assertEqual(list(self.root.iterdir()), [destination])

    def test_json_serialization_failure_preserves_existing_destination(self) -> None:
        destination = self.root / "summary.json"
        destination.write_text("sentinel\n", encoding="utf-8")

        with self.assertRaises((TypeError, ValueError)):
            self.api.write_json_atomic(destination, {"not_json": object()})

        self.assertEqual(destination.read_text(encoding="utf-8"), "sentinel\n")
        self.assertEqual(list(self.root.iterdir()), [destination])

    def test_replace_failure_preserves_destination_and_removes_temporary_file(self) -> None:
        destination = self.root / "report.txt"
        destination.write_text("sentinel\n", encoding="utf-8")

        with mock.patch(
            "make_it_stand.io.os.replace",
            side_effect=OSError("synthetic replace failure"),
        ), self.assertRaisesRegex(OSError, "synthetic replace failure"):
            self.api.write_text_atomic(destination, "replacement\n")

        self.assertEqual(destination.read_text(encoding="utf-8"), "sentinel\n")
        self.assertEqual(list(self.root.iterdir()), [destination])

    def test_csv_writer_uses_declared_column_order_and_csv_escaping(self) -> None:
        destination = self.root / "records.csv"
        rows = iter(
            (
                {"name": "alpha,one", "score": 1.25},
                {"name": 'beta "two"', "score": 2.5},
            )
        )

        self.api.write_csv_atomic(destination, ("name", "score"), rows)

        with destination.open("r", encoding="utf-8", newline="") as stream:
            parsed = list(csv.DictReader(stream))
        self.assertEqual(
            parsed,
            [
                {"name": "alpha,one", "score": "1.25"},
                {"name": 'beta "two"', "score": "2.5"},
            ],
        )
        self.assertEqual(list(self.root.iterdir()), [destination])

    def test_csv_validation_failure_preserves_existing_destination(self) -> None:
        destination = self.root / "records.csv"
        destination.write_text("sentinel\n", encoding="utf-8")

        with self.assertRaises(ValueError):
            self.api.write_csv_atomic(
                destination,
                ("known",),
                ({"known": 1, "unexpected": 2},),
            )

        self.assertEqual(destination.read_text(encoding="utf-8"), "sentinel\n")
        self.assertEqual(list(self.root.iterdir()), [destination])

    def test_svg_report_is_deterministic_and_escapes_user_text(self) -> None:
        destination = self.root / "plot.svg"
        diagnostics = {
            "status": "control <A>",
            "carving_iterations": [
                {"iteration": 0, "energy_after": 1.0},
                {"iteration": 1, "energy_after": 0.5},
                {"iteration": 2, "energy_after": 0.25},
            ],
        }

        self.api.write_svg_report(destination, diagnostics, title="A & B")
        first = destination.read_text(encoding="utf-8")
        self.api.write_svg_report(destination, diagnostics, title="A & B")

        self.assertEqual(destination.read_text(encoding="utf-8"), first)
        self.assertIn("<svg", first)
        self.assertIn("</svg>", first)
        self.assertIn("A &amp; B", first)
        self.assertIn("control &lt;A&gt;", first)
        self.assertIn("metric: energy_after", first)
        self.assertIn("<polyline", first)
        self.assertNotIn("A & B", first)

        self.api.write_svg_report(
            destination,
            {"status": "empty", "carving_iterations": [{"energy_after": math.inf}]},
        )
        filtered = destination.read_text(encoding="utf-8").lower()
        self.assertNotIn("<polyline", filtered)
        self.assertNotIn("inf", filtered)
        self.assertNotIn("nan", filtered)


def _fake_native_module() -> types.SimpleNamespace:
    state = types.SimpleNamespace(
        models=[],
        saved=[],
        fail_on_model_indices=set(),
        fail_printable=False,
    )

    class FakeModelConfig:
        pass

    class FakeRunConfig:
        pass

    class FakeResult:
        status = "stable_after_carving"
        succeeded = True
        manufacturing_ready = False
        manufacturing_failures = ("synthetic fixture is not printable",)
        manufacturing_report_json = '{"ready":false}'
        diagnostics_json = json.dumps(
            {
                "schema": "make-it-stand-r9-balance-result-v1",
                "schema_version": 1,
                "status": "stable_after_carving",
                "succeeded": True,
                "manufacturing_ready": False,
                "metrics": {
                    "projected_distance": 0.125,
                    "projected_energy": 0.0078125,
                },
            },
            allow_nan=False,
            sort_keys=True,
        )

        def _save(self, kind: str, path: str) -> None:
            state.saved.append((kind, path))

        def save_outer_mesh(self, path: str) -> None:
            self._save("outer", path)

        def save_inner_mesh(self, path: str) -> None:
            self._save("inner", path)

        def save_printable_mesh(self, path: str) -> None:
            self._save("printable", path)
            if state.fail_printable:
                raise RuntimeError("synthetic guarded printable rejection")

        def save_diagnostics(self, path: str) -> None:
            self._save("diagnostics", path)
            Path(path).write_text(self.diagnostics_json, encoding="utf-8")

    class FakeModel:
        def __init__(self, *arguments: object):
            self.arguments = arguments
            self.calls: list[tuple[object, ...]] = []
            self.index = len(state.models)
            state.models.append(self)

        def set_gravity(self, gravity: object) -> None:
            self.calls.append(("set_gravity", gravity))

        def set_support_polygon(self, polygon: object, gamma: object) -> None:
            self.calls.append(("set_support_polygon", polygon, gamma))

        def set_support_from_mesh(self, threshold: object, gamma: object) -> None:
            self.calls.append(("set_support_from_mesh", threshold, gamma))

        def add_standing_constraint_polygon(self, *arguments: object) -> None:
            self.calls.append(("add_standing_constraint_polygon", *arguments))

        def add_standing_constraint_from_mesh(self, *arguments: object) -> None:
            self.calls.append(("add_standing_constraint_from_mesh", *arguments))

        def add_handle(self, *arguments: object) -> None:
            self.calls.append(("add_handle", *arguments))

        def add_suspension_constraint(self, *arguments: object) -> None:
            self.calls.append(("add_suspension_constraint", *arguments))

        def optimize(self, run_config: object) -> FakeResult:
            self.calls.append(("optimize", run_config))
            if self.index in state.fail_on_model_indices:
                raise RuntimeError(f"synthetic failure for model {self.index}")
            return FakeResult()

    return types.SimpleNamespace(
        BalanceModel=FakeModel,
        BalanceModelConfig=FakeModelConfig,
        BalanceRunConfig=FakeRunConfig,
        BalanceResult=FakeResult,
        state=state,
    )


class MockedNativeWorkflowTests(unittest.TestCase):
    def setUp(self) -> None:
        self.api = importlib.import_module("make_it_stand")
        self.native = _fake_native_module()
        self.config = self.api.validate_config(_valid_config_dict())

    def test_single_run_builds_exactly_one_native_model_and_writes_named_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output_root = Path(directory) / "case-a"
            result = self.api.run_config(
                self.config,
                output_directory=output_root,
                native_module=self.native,
            )

        self.assertEqual(len(self.native.state.models), 1)
        model = self.native.state.models[0]
        self.assertEqual(model.arguments[0], "fixtures/model.obj")
        native_model_config = model.arguments[1]
        for name in (
            "voxel_resolution",
            "minimum_wall_thickness",
            "bbw_max_iterations",
            "bbw_tolerance",
        ):
            self.assertEqual(
                getattr(native_model_config, name),
                getattr(self.config.model, name),
            )

        self.assertEqual(model.calls[0], ("set_gravity", (0.0, -9.81, 0.0)))
        self.assertEqual(model.calls[1][0], "set_support_polygon")
        self.assertEqual(model.calls[2][0], "add_suspension_constraint")
        self.assertEqual(model.calls[3][0], "add_handle")
        self.assertEqual(model.calls[4][0], "optimize")
        native_run_config = model.calls[4][1]
        for field in dataclasses.fields(self.config.run):
            self.assertEqual(
                getattr(native_run_config, field.name),
                getattr(self.config.run, field.name),
            )

        expected_root = output_root
        self.assertEqual(
            self.native.state.saved,
            [
                ("diagnostics", str(expected_root / "diagnostics.json")),
                ("outer", str(expected_root / "outer.obj")),
                ("inner", str(expected_root / "inner.obj")),
                ("printable", str(expected_root / "printable.stl")),
            ],
        )
        self.assertIsInstance(result, self.api.BalanceResult)
        self.assertEqual(result.status, "stable_after_carving")
        self.assertTrue(result.succeeded)
        self.assertFalse(result.manufacturing_ready)
        self.assertEqual(
            result.manufacturing_failures,
            ("synthetic fixture is not printable",),
        )
        self.assertEqual(result.diagnostics["metrics"]["projected_distance"], 0.125)

    def test_mesh_support_uses_only_the_native_mesh_support_operation(self) -> None:
        raw = _valid_config_dict()
        raw["support"] = {
            "kind": "mesh",
            "contact_threshold": 0.02,
            "gamma_radians": 0.15,
        }
        config = self.api.validate_config(raw)
        self.api.build_model(config, native_module=self.native)

        calls = self.native.state.models[0].calls
        support_calls = [call for call in calls if str(call[0]).startswith("set_support")]
        self.assertEqual(
            support_calls,
            [("set_support_from_mesh", 0.02, 0.15)],
        )

    def test_in_memory_model_constructor_is_forwarded_once(self) -> None:
        vertices = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
        faces = [[0, 1, 2]]
        model = self.api.BalanceModel(
            vertices,
            faces,
            native_module=self.native,
        )

        self.assertIsNotNone(model)
        self.assertEqual(len(self.native.state.models), 1)
        arguments = self.native.state.models[0].arguments
        self.assertIs(arguments[0], vertices)
        self.assertIs(arguments[1], faces)

    def test_wrapper_rejects_bad_vectors_before_calling_native_code(self) -> None:
        model = self.api.BalanceModel("model.obj", native_module=self.native)
        with self.assertRaises(self.api.ConfigurationError):
            model.set_gravity([math.nan, -1.0, 0.0])
        with self.assertRaises(self.api.ConfigurationError):
            model.add_handle([0.0, 1.0])
        with self.assertRaises(self.api.ConfigurationError):
            model.add_suspension_constraint([0.0, 1.0, 0.0], math.inf)

        self.assertEqual(self.native.state.models[0].calls, [])

    def test_result_diagnostics_reject_duplicate_keys_and_non_object_json(self) -> None:
        native_result = self.native.BalanceResult()
        wrapped = self.api.BalanceResult(native_result)

        native_result.diagnostics_json = '{"status":"a","status":"b"}'
        with self.assertRaises(self.api.DataFormatError):
            _ = wrapped.diagnostics

        native_result.diagnostics_json = "[]"
        with self.assertRaises(self.api.DataFormatError):
            _ = wrapped.diagnostics

        native_result.diagnostics_json = json.dumps(
            {"status": "stable", "succeeded": True}
        )
        with self.assertRaisesRegex(self.api.DataFormatError, "schema"):
            _ = wrapped.diagnostics

    def test_nested_outputs_create_parents_and_diagnostics_is_attempted_first(self) -> None:
        raw = _valid_config_dict()
        raw["outputs"] = {
            "diagnostics": "records/diagnostics.json",
            "outer_mesh": "meshes/outer.obj",
            "inner_mesh": "meshes/cavity/inner.obj",
            "printable_mesh": "print/ready.stl",
        }
        config = self.api.validate_config(raw)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.api.run_config(config, output_directory=root, native_module=self.native)
            for parent in ("records", "meshes", "meshes/cavity", "print"):
                self.assertTrue((root / parent).is_dir())
        self.assertEqual(self.native.state.saved[0][0], "diagnostics")

    def test_guarded_printable_failure_preserves_published_diagnostics(self) -> None:
        raw = _valid_config_dict()
        raw["outputs"] = {
            "diagnostics": "records/diagnostics.json",
            "printable_mesh": "print/ready.stl",
        }
        config = self.api.validate_config(raw)
        self.native.state.fail_printable = True
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with self.assertRaisesRegex(RuntimeError, "guarded printable"):
                self.api.run_config(
                    config,
                    output_directory=root,
                    native_module=self.native,
                )
            diagnostic_path = root / "records/diagnostics.json"
            self.assertTrue(diagnostic_path.is_file())
            self.assertEqual(
                json.loads(diagnostic_path.read_text(encoding="utf-8"))["status"],
                "stable_after_carving",
            )
        self.assertEqual(
            [kind for kind, _ in self.native.state.saved],
            ["diagnostics", "printable"],
        )

    def test_nested_output_parent_and_final_symlinks_are_rejected(self) -> None:
        raw = _valid_config_dict()
        raw["outputs"] = {"diagnostics": "records/diagnostics.json"}
        config = self.api.validate_config(raw)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "selected-root"
            root.mkdir()
            outside = Path(directory) / "outside"
            outside.mkdir()
            (root / "records").symlink_to(outside, target_is_directory=True)
            with self.assertRaisesRegex(self.api.UnsafePathError, "symbolic-link"):
                self.api.run_config(
                    config,
                    output_directory=root,
                    native_module=self.native,
                )
            self.assertEqual(self.native.state.models, [])

            (root / "records").unlink()
            (root / "records").mkdir()
            (root / "records" / "diagnostics.json").symlink_to(
                outside / "escaped.json"
            )
            with self.assertRaisesRegex(self.api.UnsafePathError, "symbolic-link"):
                self.api.run_config(
                    config,
                    output_directory=root,
                    native_module=self.native,
                )
            self.assertFalse((outside / "escaped.json").exists())

    def test_caller_selected_output_root_may_itself_be_a_symlink(self) -> None:
        raw = _valid_config_dict()
        raw["outputs"] = {"diagnostics": "records/diagnostics.json"}
        config = self.api.validate_config(raw)
        with tempfile.TemporaryDirectory() as directory:
            actual = Path(directory) / "actual-root"
            actual.mkdir()
            selected = Path(directory) / "selected-root"
            selected.symlink_to(actual, target_is_directory=True)
            self.api.run_config(
                config,
                output_directory=selected,
                native_module=self.native,
            )
            self.assertTrue((actual / "records" / "diagnostics.json").is_file())


class BatchWorkflowTests(unittest.TestCase):
    def setUp(self) -> None:
        self.api = importlib.import_module("make_it_stand")
        self.config = self.api.validate_config(_valid_config_dict())
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary_directory.cleanup)
        self.root = Path(self.temporary_directory.name)

    def test_batch_runs_cases_in_order_in_deterministic_subdirectories(self) -> None:
        native = _fake_native_module()
        summary = self.api.run_batch(
            [self.config, self.config],
            output_directory=self.root,
            native_module=native,
        )

        self.assertEqual(summary["schema"], "make-it-stand-r9-batch-summary-v1")
        self.assertEqual(summary["schema_version"], 1)
        self.assertEqual(summary["case_count"], 2)
        self.assertEqual(summary["completed_count"], 2)
        self.assertEqual([case["index"] for case in summary["cases"]], [0, 1])
        self.assertTrue(all(case["completed"] for case in summary["cases"]))
        self.assertEqual(len(native.state.models), 2)

        diagnostics_paths = [
            Path(path)
            for kind, path in native.state.saved
            if kind == "diagnostics"
        ]
        self.assertEqual(
            diagnostics_paths,
            [
                self.root / "case-0000" / "diagnostics.json",
                self.root / "case-0001" / "diagnostics.json",
            ],
        )
        persisted = json.loads(
            (self.root / "batch_summary.json").read_text(encoding="utf-8")
        )
        self.assertEqual(persisted, summary)
        serialized = json.dumps(summary, sort_keys=True)
        self.assertNotIn(str(self.root), serialized)

    def test_keep_going_records_one_failure_without_losing_other_cases(self) -> None:
        native = _fake_native_module()
        native.state.fail_on_model_indices.add(1)

        summary = self.api.run_batch(
            [self.config, self.config, self.config],
            output_directory=self.root,
            native_module=native,
            keep_going=True,
        )

        self.assertEqual(summary["case_count"], 3)
        self.assertEqual(summary["completed_count"], 2)
        self.assertEqual(len(native.state.models), 3)
        failed = summary["cases"][1]
        self.assertEqual(
            failed,
            {
                "index": 1,
                "completed": False,
                "error_type": "RuntimeError",
                "error": "synthetic failure for model 1",
            },
        )
        self.assertTrue(summary["cases"][2]["completed"])
        self.assertEqual(
            json.loads((self.root / "batch_summary.json").read_text()),
            summary,
        )

    def test_fail_fast_reraises_and_does_not_publish_a_partial_summary(self) -> None:
        native = _fake_native_module()
        native.state.fail_on_model_indices.add(1)

        with self.assertRaisesRegex(RuntimeError, "synthetic failure"):
            self.api.run_batch(
                [self.config, self.config, self.config],
                output_directory=self.root,
                native_module=native,
            )

        self.assertEqual(len(native.state.models), 2)
        self.assertFalse((self.root / "batch_summary.json").exists())

    def test_batch_refuses_existing_summary_or_case_path(self) -> None:
        (self.root / "case-0000").mkdir()
        (self.root / "case-0000" / "stale.json").write_text("stale\n")
        native = _fake_native_module()
        with self.assertRaisesRegex(FileExistsError, "existing batch case"):
            self.api.run_batch(
                [self.config],
                output_directory=self.root,
                native_module=native,
            )
        self.assertEqual(native.state.models, [])

        (self.root / "case-0000" / "stale.json").unlink()
        (self.root / "case-0000").rmdir()
        (self.root / "batch_summary.json").write_text("stale\n")
        with self.assertRaisesRegex(FileExistsError, "batch artifact"):
            self.api.run_batch(
                [self.config],
                output_directory=self.root,
                native_module=native,
            )
        self.assertEqual(native.state.models, [])

    def test_batch_refuses_existing_exclusive_claim(self) -> None:
        claim = self.root / ".r9-batch.claim"
        claim.write_text("another process\n", encoding="utf-8")
        native = _fake_native_module()
        with self.assertRaisesRegex(FileExistsError, "batch output claim"):
            self.api.run_batch(
                [self.config],
                output_directory=self.root,
                native_module=native,
            )
        self.assertEqual(native.state.models, [])
        self.assertEqual(claim.read_text(encoding="utf-8"), "another process\n")

    def test_batch_revalidates_programmatic_configs_before_claiming_output(self) -> None:
        invalid = dataclasses.replace(
            self.config,
            outputs=dataclasses.replace(
                self.config.outputs,
                diagnostics="../escape.json",
            ),
        )
        native = _fake_native_module()
        with self.assertRaisesRegex(self.api.ConfigurationError, "parent traversal"):
            self.api.run_batch(
                [invalid],
                output_directory=self.root,
                native_module=native,
            )
        self.assertEqual(native.state.models, [])
        self.assertEqual(list(self.root.iterdir()), [])


class DiagnosticsToolsTests(unittest.TestCase):
    def setUp(self) -> None:
        self.api = importlib.import_module("make_it_stand")
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary_directory.cleanup)
        self.root = Path(self.temporary_directory.name)
        self.records = [
            _diagnostics(**{
                "status": "stable",
                "manufacturing_ready": True,
                "manufacturing": {"ready": True},
                "carving_iterations": [{"iteration": 0}, {"iteration": 1}],
                "deformation_iterations": [{}, {}, {}],
                "mass_properties": {"mass": 2.0},
            }),
            _diagnostics(**{
                "status": "maximum_iterations",
                "succeeded": False,
                "manufacturing": {"ready": True},
                "carving_iterations": [{}, {}, {}, {}],
                "deformation_iterations": [{"iteration": 0}],
                "mass": 4.0,
            }),
            _diagnostics(**{
                "status": "stable",
                "succeeded": False,
                "manufacturing_ready": False,
                "carving_iterations": [],
                "deformation_iterations": [],
            }),
        ]

    def test_load_diagnostics_parses_strict_json_and_validates_core_fields(self) -> None:
        path = self.root / "diagnostics.json"
        self.api.write_json_atomic(path, self.records[0])
        self.assertEqual(self.api.load_diagnostics(path), self.records[0])

        invalid_values = (
            [],
            {},
            {
                key: value
                for key, value in _diagnostics().items()
                if key != "schema"
            },
            _diagnostics(schema="wrong"),
            _diagnostics(schema_version=2),
            _diagnostics(schema_version=True),
            _diagnostics(schema_version=1.0),
            _diagnostics(status=""),
            _diagnostics(status=1),
            _diagnostics(succeeded=1),
            _diagnostics(manufacturing_ready=1),
            {key: item for key, item in _diagnostics().items() if key != "manufacturing"},
            _diagnostics(manufacturing=[]),
            _diagnostics(manufacturing={"ready": 1}),
            _diagnostics(
                manufacturing_ready=True,
                manufacturing={"ready": False},
            ),
            _diagnostics(carving_iterations=3),
            _diagnostics(deformation_iterations={}),
            _diagnostics(mass="not-a-number"),
            {key: item for key, item in _diagnostics().items() if key != "mass_properties"},
            _diagnostics(mass_properties={}),
            _diagnostics(mass_properties={"mass": True}),
        )
        for index, value in enumerate(invalid_values):
            path = self.root / f"invalid-{index}.json"
            self.api.write_json_atomic(path, value)
            with self.subTest(value=value), self.assertRaises(self.api.DiagnosticsError):
                self.api.load_diagnostics(path)

        path.write_text('{"status":"a","status":"b"}', encoding="utf-8")
        with self.assertRaises(self.api.DataFormatError):
            self.api.load_diagnostics(path)

    def test_summary_counts_statuses_traces_successes_and_only_finite_masses(self) -> None:
        summary = self.api.summarize_diagnostics(iter(self.records))

        self.assertEqual(summary["schema"], "make-it-stand-r9-diagnostics-summary-v1")
        self.assertEqual(summary["schema_version"], 1)
        self.assertEqual(summary["record_count"], 3)
        self.assertEqual(summary["succeeded_count"], 1)
        self.assertEqual(summary["manufacturing_ready_count"], 2)
        self.assertEqual(
            summary["status_counts"],
            {"maximum_iterations": 1, "stable": 2},
        )
        self.assertEqual(summary["carving_iteration_count"], 6)
        self.assertEqual(summary["deformation_iteration_count"], 4)
        self.assertEqual(summary["mean_mass"], 7.0 / 3.0)

    def test_summary_rejects_non_objects_and_missing_or_empty_status(self) -> None:
        for records in (
            ["bad"],
            [{}],
            [_diagnostics(status="")],
            [_diagnostics(mass=math.nan)],
        ):
            with self.subTest(records=records), self.assertRaises(
                self.api.DiagnosticsError
            ):
                self.api.summarize_diagnostics(records)

    def test_compare_is_sorted_path_free_and_explicitly_not_an_r10_ablation(self) -> None:
        comparison = self.api.compare_diagnostics(
            {"z-control": self.records[:1], "a-treatment": self.records[1:]}
        )

        self.assertEqual(
            comparison["schema"],
            "make-it-stand-r9-generic-diagnostics-comparison-v1",
        )
        self.assertEqual(list(comparison["groups"]), ["a-treatment", "z-control"])
        self.assertIn("no R10 paper ablation", comparison["claim"])
        self.assertEqual(comparison["groups"]["z-control"]["record_count"], 1)

        for groups in ({}, {"": self.records}, {"bad": "records"}):
            with self.subTest(groups=groups), self.assertRaises(
                self.api.DiagnosticsError
            ):
                self.api.compare_diagnostics(groups)

    def test_trace_flattening_csv_and_svg_ignore_boolean_and_nonfinite_metrics(self) -> None:
        diagnostics = _diagnostics(**{
            "status": "stable",
            "carving_iterations": [
                {
                    "block": 1,
                    "iteration": 2,
                    "energy_after": 0.25,
                    "accepted": True,
                    "bad": math.inf,
                }
            ],
            "deformation_iterations": [
                {"block": 1, "iteration": 3, "gradient_norm": 0.125}
            ],
        })
        rows = self.api.diagnostics_trace_rows(diagnostics)
        metrics = {row["metric"] for row in rows}
        self.assertIn("energy_after", metrics)
        self.assertIn("gradient_norm", metrics)
        self.assertNotIn("accepted", metrics)
        self.assertNotIn("bad", metrics)

        csv_path = self.root / "trace.csv"
        svg_path = self.root / "trace.svg"
        self.api.write_trace_csv(csv_path, diagnostics)
        self.api.write_svg_report(svg_path, diagnostics)
        with csv_path.open("r", encoding="utf-8", newline="") as stream:
            csv_rows = list(csv.DictReader(stream))
        self.assertEqual(
            tuple(csv_rows[0]),
            ("trace", "block", "iteration", "metric", "value"),
        )
        self.assertIn("<polyline", svg_path.read_text(encoding="utf-8"))

    def test_svg_prefers_the_actual_cpp_carving_energy_field(self) -> None:
        diagnostics = _diagnostics(**{
            "carving_iterations": [
                {"iteration": 0, "carved_voxels": 10, "refreshed_energy_after": 1.0},
                {"iteration": 1, "carved_voxels": 20, "refreshed_energy_after": 0.5},
            ],
        })
        svg_path = self.root / "cpp-trace.svg"
        self.api.write_svg_report(svg_path, diagnostics)
        self.assertIn(
            "metric: refreshed_energy_after",
            svg_path.read_text(encoding="utf-8"),
        )


class CliContractTests(unittest.TestCase):
    def setUp(self) -> None:
        self.cli = importlib.import_module("make_it_stand.cli")

    def test_run_cli_reports_readiness_and_requires_both_success_conditions(self) -> None:
        for succeeded, ready, expected_code in (
            (True, True, 0),
            (True, False, 3),
            (False, True, 3),
        ):
            result = types.SimpleNamespace(
                status="fixture",
                succeeded=succeeded,
                manufacturing_ready=ready,
            )
            with self.subTest(succeeded=succeeded, ready=ready), mock.patch(
                "make_it_stand.cli.load_config", return_value=object()
            ), mock.patch(
                "make_it_stand.cli.run_config", return_value=result
            ), mock.patch("builtins.print") as printed:
                code = self.cli.run_balance_main(["config.json"])
            self.assertEqual(code, expected_code)
            payload = json.loads(printed.call_args.args[0])
            self.assertIs(payload["manufacturing_ready"], ready)
            self.assertIs(payload["succeeded"], succeeded)

    def test_batch_cli_requires_every_completed_case_to_be_manufacturing_ready(self) -> None:
        summary = {
            "case_count": 1,
            "completed_count": 1,
            "cases": [
                {
                    "completed": True,
                    "succeeded": True,
                    "manufacturing_ready": False,
                }
            ],
        }
        with mock.patch(
            "make_it_stand.cli._load_batch_manifest", return_value=[]
        ), mock.patch(
            "make_it_stand.cli.run_batch", return_value=summary
        ), mock.patch("builtins.print"):
            code = self.cli.batch_experiment_main(
                ["manifest.json", "--output-directory", "results"]
            )
        self.assertEqual(code, 3)

    def test_batch_and_comparison_schema_versions_reject_bool_and_float(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            batch = root / "batch.json"
            comparison = root / "comparison.json"
            for version in (True, 1.0):
                batch.write_text(
                    json.dumps(
                        {
                            "schema": "make-it-stand-r9-batch-config-v1",
                            "schema_version": version,
                            "configs": ["case.json"],
                        }
                    ),
                    encoding="utf-8",
                )
                comparison.write_text(
                    json.dumps(
                        {
                            "schema": "make-it-stand-r9-comparison-config-v1",
                            "schema_version": version,
                            "groups": {"a": ["diagnostics.json"]},
                        }
                    ),
                    encoding="utf-8",
                )
                with self.subTest(kind="batch", version=version), self.assertRaisesRegex(
                    ValueError, "integer 1"
                ):
                    self.cli._load_batch_manifest(batch)
                with self.subTest(
                    kind="comparison", version=version
                ), self.assertRaisesRegex(ValueError, "integer 1"):
                    self.cli._load_comparison_manifest(comparison)


if __name__ == "__main__":
    unittest.main()
