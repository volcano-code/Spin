import hashlib
import json
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiments/teapot_spin2014"))

from render_reconstructed_g7_video import (  # noqa: E402
    interpolate_recorded_poses,
    load_sealed_mesh,
    model_origins_from_com_poses,
    select_contact_run,
    validate_seal_report_binding,
)


class ReconstructedG7VideoEvidenceTests(unittest.TestCase):
    def _report(self):
        return {
            "schema": "spin2014-reconstructed-g7-free-contact-v1",
            "scope": "RECONSTRUCTED_G7_SAME_SEALED_MATERIAL_CONTACT_OBSERVATION",
            "evidence_level": "RECONSTRUCTED",
            "original_paper_case": False,
            "artifact_digest": "a" * 64,
            "protocol": {
                "type": "PYBULLET_FREE_TABLE_CONTACT_NO_FIXED_CONSTRAINT",
                "duration_s": 8.0,
                "dt_s": 0.0005,
                "spin_hz": 35.0,
                "initial_tilt_deg": 2.0,
                "directions": 8,
            },
            "runs": [
                {
                    "profile": "default",
                    "direction": 0,
                    "protocol": "PYBULLET_FREE_TABLE_CONTACT_NO_FIXED_CONSTRAINT",
                    "artifact_digest": "a" * 64,
                    "samples": [
                        {"time": 0.0, "position": [0.0, 0.0, 0.0],
                         "model_quaternion": [0.0, 0.0, 0.0, 1.0]},
                        {"time": 1.0, "position": [1.0, 2.0, 3.0],
                         "model_quaternion": [0.0, 0.0, 1.0, 0.0]},
                    ],
                }
            ],
        }

    def test_selects_only_the_requested_run_bound_to_the_sealed_artifact(self):
        report = self._report()
        selected = select_contact_run(report, "a" * 64, profile="default", direction=0)
        self.assertEqual((selected["profile"], selected["direction"]), ("default", 0))
        with self.assertRaisesRegex(ValueError, "ARTIFACT_DIGEST_MISMATCH"):
            select_contact_run(report, "b" * 64, profile="default", direction=0)
        report["runs"].append(dict(report["runs"][0]))
        with self.assertRaisesRegex(ValueError, "CONTACT_RUN_NOT_UNIQUE"):
            select_contact_run(report, "a" * 64, profile="default", direction=0)

    def test_interpolates_position_and_quaternion_without_extrapolation(self):
        record = self._report()["runs"][0]
        poses = interpolate_recorded_poses(record, np.array([0.0, 0.5, 1.0]))
        np.testing.assert_allclose(poses["positions"][1], [0.5, 1.0, 1.5], atol=1e-14)
        # Halfway between identity and a 180-degree z rotation is 90 degrees.
        np.testing.assert_allclose(poses["quaternions"][1], [0.0, 0.0, 2**-0.5, 2**-0.5], atol=1e-12)
        with self.assertRaisesRegex(ValueError, "POSE_TIME_OUTSIDE_RECORDED_INTERVAL"):
            interpolate_recorded_poses(record, np.array([-0.001]))
        record["samples"][1]["time"] = 0.0
        with self.assertRaisesRegex(ValueError, "INVALID_RECORDED_TIME_SEQUENCE"):
            interpolate_recorded_poses(record, np.array([0.0]))

    def test_converts_pybullet_com_pose_to_sealed_mesh_origin(self):
        positions = np.array([[1.0, 2.0, 3.0]])
        # 90 degrees around z: local center [2, 0, 0] becomes world offset [0, 2, 0].
        quaternions = np.array([[0.0, 0.0, 2**-0.5, 2**-0.5]])
        origins = model_origins_from_com_poses(positions, quaternions, [2.0, 0.0, 0.0])
        np.testing.assert_allclose(origins, [[1.0, 0.0, 3.0]], atol=1e-12)

    def test_reads_geometry_digest_from_independent_seal_section(self):
        report = {
            "schema": "spin2014-reconstructed-g7-sealed-geometry-v1",
            "source_binding": "VERIFIED",
            "problem_digest": "p" * 64,
            "independent_seal": {
                "status": "MESH_MATERIAL_VALIDATED",
                "scientific_digest": "a" * 64,
            },
        }
        self.assertTrue(validate_seal_report_binding(report, "a" * 64, "p" * 64))
        self.assertFalse(validate_seal_report_binding(report, "b" * 64, "p" * 64))
        self.assertFalse(validate_seal_report_binding(report, "a" * 64, "q" * 64))

    def test_sealed_mesh_loader_checks_hash_shape_and_topology(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            vertices = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]])
            faces = np.array([[0, 1, 2]], dtype=np.int32)
            np.save(root / "v.npy", vertices, allow_pickle=False)
            np.save(root / "f.npy", faces, allow_pickle=False)
            def descriptor(name, array):
                path = root / name
                return {"name": name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                        "size": path.stat().st_size, "shape": list(array.shape),
                        "dtype": array.dtype.str}
            manifest = {"final": {"outer": {
                "vertices": descriptor("v.npy", vertices),
                "faces": descriptor("f.npy", faces),
            }}}
            actual_v, actual_f = load_sealed_mesh(root, manifest, "outer")
            np.testing.assert_array_equal(actual_v, vertices)
            np.testing.assert_array_equal(actual_f, faces)
            np.save(root / "f.npy", np.array([[0, 1, 9]], dtype=np.int32), allow_pickle=False)
            manifest["final"]["outer"]["faces"] = descriptor("f.npy", np.array([[0, 1, 9]], dtype=np.int32))
            with self.assertRaisesRegex(ValueError, "SEALED_MESH_TOPOLOGY_INVALID"):
                load_sealed_mesh(root, manifest, "outer")


if __name__ == "__main__":
    unittest.main()
