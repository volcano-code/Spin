from __future__ import annotations

import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from make_it_stand.pybullet_runner import (
    SimulationConfig,
    _alignment_rotation,
    _contact_records,
    _load_mesh_properties,
    _mass_properties_fingerprint,
    _tilt_from_native_gravity_deg,
    load_manifest,
)

import numpy as np


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _tetra(path: Path, scale: float, inward: bool = False, open_mesh: bool = False) -> None:
    faces = [(1, 3, 2), (1, 2, 4), (1, 4, 3), (2, 3, 4)]
    if inward:
        faces = [(a, c, b) for a, b, c in faces]
    if open_mesh:
        faces.pop()
    lines = [
        "v 0 0 0",
        f"v {scale} 0 0",
        f"v 0 {scale} 0",
        f"v 0 0 {scale}",
    ] + [f"f {a} {b} {c}" for a, b, c in faces]
    path.write_text("\n".join(lines) + "\n", encoding="ascii")


class PyBulletMaterialContractTests(unittest.TestCase):
    def test_occupancy_bound_analytic_moments_do_not_use_open_inner_mesh(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            collision = root / "open-collision.obj"
            _tetra(collision, 1.0, open_mesh=True)
            entry = {
                "name": "analytic",
                "collision_mesh": collision.name,
                "visual_mesh": collision.name,
                "occupancy_hash": "alpha-analytic",
                "geometry_hash": "geometry-analytic",
                "collision_mesh_sha256": _sha256(collision),
                "material_mass_native": 2.0,
                "center_of_mass_native": [0.25, 0.25, 0.25],
                "inertia_tensor_about_center_native": [
                    [0.15, 0.01, 0.0],
                    [0.01, 0.16, 0.0],
                    [0.0, 0.0, 0.17],
                ],
                "gravity_direction_native": [0.0, -1.0, 0.0],
            }
            entry["mass_properties_sha256"] = _mass_properties_fingerprint(entry)
            manifest = root / "manifest.json"
            manifest.write_text(json.dumps({"cases": [entry]}), encoding="utf-8")

            case = load_manifest(manifest)[0]
            properties = _load_mesh_properties(case, 0.1, 0.2)
            self.assertEqual(
                properties["mass_properties_source"],
                "occupancy_hash_bound_analytic_volume_moments",
            )
            self.assertEqual(properties["occupancy_hash"], "alpha-analytic")
            self.assertAlmostEqual(properties["material_volume_native"], 2.0)
            self.assertEqual(case.gravity_direction_native, (0.0, -1.0, 0.0))

            entry["material_mass_native"] = 3.0
            manifest.write_text(json.dumps({"cases": [entry]}), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "mass-properties hash"):
                load_manifest(manifest)

    def test_occupancy_bound_closed_material_model_is_accepted(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            outer = root / "outer.obj"
            inner = root / "inner.obj"
            _tetra(outer, 1.0)
            _tetra(inner, 0.25, inward=True)
            manifest = root / "manifest.json"
            manifest.write_text(
                json.dumps({"cases": [{
                    "name": "closed",
                    "collision_mesh": outer.name,
                    "visual_mesh": outer.name,
                    "inner_void_mesh": inner.name,
                    "occupancy_hash": "alpha-123",
                    "collision_mesh_sha256": _sha256(outer),
                    "inner_void_mesh_sha256": _sha256(inner),
                }]}),
                encoding="utf-8",
            )
            case = load_manifest(manifest)[0]
            properties = _load_mesh_properties(case, 0.1, 0.2)
            self.assertEqual(properties["occupancy_hash"], "alpha-123")
            self.assertEqual(
                properties["mass_properties_source"],
                "occupancy_hash_bound_closed_outer_minus_inner",
            )

    def test_occupancy_bound_model_rejects_open_inner_boundary(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            outer = root / "outer.obj"
            inner = root / "inner.obj"
            _tetra(outer, 1.0)
            _tetra(inner, 0.25, inward=True, open_mesh=True)
            manifest = root / "manifest.json"
            manifest.write_text(
                json.dumps({"cases": [{
                    "name": "open",
                    "collision_mesh": outer.name,
                    "inner_void_mesh": inner.name,
                    "occupancy_hash": "alpha-456",
                    "collision_mesh_sha256": _sha256(outer),
                    "inner_void_mesh_sha256": _sha256(inner),
                }]}),
                encoding="utf-8",
            )
            case = load_manifest(manifest)[0]
            with self.assertRaisesRegex(ValueError, "closed material boundary"):
                _load_mesh_properties(case, 0.1, 0.2)

    def test_perturbation_configuration_is_explicit(self) -> None:
        config = SimulationConfig(
            initial_tilt_deg=2.0,
            initial_tilt_azimuth_deg=315.0,
        )
        self.assertEqual(config.initial_tilt_deg, 2.0)
        self.assertEqual(config.initial_tilt_azimuth_deg, 315.0)
        self.assertEqual(
            config.perturbation_protocol,
            "reseeded_initial_pose_drop",
        )

    def test_contact_records_preserve_points_normals_and_forces(self) -> None:
        contact = (
            0, 1, 2, -1, -1,
            (0.1, 0.2, 0.3),
            (0.1, 0.2, 0.0),
            (0.0, 0.0, 1.0),
            -0.001,
            3.5,
            0.4,
            (1.0, 0.0, 0.0),
            0.2,
            (0.0, 1.0, 0.0),
        )
        records = _contact_records([contact])
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0]["position_on_object_world_m"], [0.1, 0.2, 0.3])
        self.assertEqual(records[0]["position_on_support_world_m"], [0.1, 0.2, 0.0])
        self.assertEqual(records[0]["support_normal_world"], [0.0, 0.0, 1.0])
        self.assertEqual(records[0]["normal_force_n"], 3.5)

    def test_native_gravity_is_aligned_with_world_gravity(self) -> None:
        gravity = np.asarray([0.0, -1.0, 0.0])
        alignment = _alignment_rotation(gravity)
        np.testing.assert_allclose(
            alignment @ gravity,
            np.asarray([0.0, 0.0, -1.0]),
            atol=1.0e-12,
        )
        np.testing.assert_allclose(
            alignment.T @ alignment,
            np.eye(3),
            atol=1.0e-12,
        )
        self.assertAlmostEqual(float(np.linalg.det(alignment)), 1.0)
        self.assertAlmostEqual(
            _tilt_from_native_gravity_deg(alignment, gravity),
            0.0,
            places=10,
        )

        angle = np.deg2rad(2.0)
        perturbation = np.asarray([
            [1.0, 0.0, 0.0],
            [0.0, np.cos(angle), -np.sin(angle)],
            [0.0, np.sin(angle), np.cos(angle)],
        ])
        self.assertAlmostEqual(
            _tilt_from_native_gravity_deg(perturbation @ alignment, gravity),
            2.0,
            places=10,
        )

    def test_manifest_rejects_invalid_native_gravity(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            collision = root / "tetra.obj"
            _tetra(collision, 1.0)
            manifest = root / "manifest.json"
            manifest.write_text(json.dumps({"cases": [{
                "name": "bad-gravity",
                "collision_mesh": collision.name,
                "gravity_direction_native": [0.0, 0.0, 0.0],
            }]}), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "gravity_direction_native"):
                load_manifest(manifest)


if __name__ == "__main__":
    unittest.main()
