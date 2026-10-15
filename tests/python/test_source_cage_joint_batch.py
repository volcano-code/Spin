"""Exercise the native experiment protocol, not a manufactured success flag."""
import json
import importlib.util
import os
from pathlib import Path
import subprocess
import tarfile
import tempfile
import unittest


class JointBatchProtocolTests(unittest.TestCase):
    def test_full_attempt_retains_terminal_reasons_and_pass_metadata(self):
        root = Path(__file__).resolve().parents[2]
        binary = Path(os.environ.get("SPIN_JOINT_BATCH_BINARY", str(
            root / "build/strict-algorithms-exact-current-20260922/"
            "spin-it-2014-source-cage-asymmetric-joint-regression")))
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "run"
            process = subprocess.run([str(binary), "--batch-run", str(output)],
                                     capture_output=True, text=True, timeout=180)
            self.assertEqual(process.returncode, 0, process.stdout + process.stderr)
            result = json.loads((output / "result.json").read_text())
            self.assertIn("adaptive_reason", result)
            self.assertIn("binary_reason", result)
            self.assertTrue((output / "pass_1_metadata.json").is_file())
            metadata = json.loads((output / "pass_1_metadata.json").read_text())
            self.assertIn("rebind", metadata)
            if result["binary_status"] == "NOT_ADMITTED":
                self.assertTrue(result["binary_reason"])
                self.assertFalse((output / "binary_material.obj").exists())

    def test_first_primitive_subset_rejection_exports_a_replayable_witness(self):
        """A geometry failure must identify the same state and mesh inputs."""
        root = Path(__file__).resolve().parents[2]
        binary = Path(os.environ.get("SPIN_JOINT_BATCH_BINARY", str(
            root / "build/strict-algorithms-exact-current-20260922/"
            "spin-it-2014-source-cage-asymmetric-joint-regression")))
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "run"
            process = subprocess.run([str(binary), "--batch-run", str(output)],
                                     capture_output=True, text=True, timeout=180)
            self.assertEqual(process.returncode, 0, process.stdout + process.stderr)
            witnesses = output / "primitive_subset_witnesses.jsonl"
            self.assertTrue(witnesses.is_file(),
                            "rejected primitive containment must preserve a replay witness")
            rows = [json.loads(line) for line in witnesses.read_text().splitlines()]
            self.assertTrue(rows)
            first = rows[0]
            self.assertEqual(first["primitive"], 0)
            self.assertGreaterEqual(first["callback"], 0)
            self.assertEqual(first["replayed_evaluation_output_sha256"],
                             first["geometry_output_sha256"],
                             "the witness must replay the exact evaluation consumed by the guard")
            self.assertEqual(len(first["state"]), first["state_dimension"])
            self.assertTrue((output / first["outer_obj"]).is_file())
            self.assertTrue((output / first["primitive_obj"]).is_file())
            self.assertEqual(len(first["primitive_source_indices"]),
                             first["primitive_vertex_count"])
            self.assertEqual(first["weight_rows"], first["primitive_vertex_count"])
            self.assertEqual(first["weight_columns"], 98)

    def test_failure_archive_preserves_empty_and_truncated_outputs(self):
        root = Path(__file__).resolve().parents[2]
        path = root / "experiments/spin_it_reproduction/strict/run_source_cage_joint_batch.py"
        spec = importlib.util.spec_from_file_location("joint_batch_runner", path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        self.assertTrue(hasattr(module, "pack_raw_outputs"), "partial output must remain bindable evidence")
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            (output / "native").mkdir()
            with (output / "native.stdout").open("w"):
                pass
            with (output / "native/result.json").open("w") as stream:
                stream.write('{"status":')
            archive = module.pack_raw_outputs(output)
            with tarfile.open(archive) as stream:
                self.assertEqual(stream.extractfile("native.stdout").read(), b"")
                self.assertEqual(stream.extractfile("native/result.json").read(), b'{"status":')

    def test_scientific_digest_excludes_diagnostic_witness_attachments(self):
        root = Path(__file__).resolve().parents[2]
        path = root / "experiments/spin_it_reproduction/strict/run_source_cage_joint_batch.py"
        spec = importlib.util.spec_from_file_location("joint_batch_runner", path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        self.assertTrue(hasattr(module, "scientific_file_digests"),
                        "the evidence runner must distinguish solver output from replay attachments")
        with tempfile.TemporaryDirectory() as temporary:
            native = Path(temporary)
            (native / "result.json").write_text('{"status":"NOT_ADMITTED"}\n')
            (native / "pass_1_callbacks.jsonl").write_text('{"callback":1}\n')
            witness = native / "primitive_subset_witnesses.jsonl"
            witness.write_text('{"weights":[1.0]}\n')
            (native / "primitive_subset_callback_1_outer.obj").write_text("v 0 0 0\n")
            science_a, diagnostic_a = module.scientific_file_digests(native)
            witness.write_text('{"weights":[1.0000000001]}\n')
            science_b, diagnostic_b = module.scientific_file_digests(native)
            self.assertEqual(science_a, science_b)
            self.assertNotEqual(diagnostic_a, diagnostic_b)
            self.assertIn("primitive_subset_witnesses.jsonl", diagnostic_a)

    def test_geometry_checker_recomputes_axis_constraints(self):
        root = Path(__file__).resolve().parents[2]
        path = root / "experiments/spin_it_reproduction/strict/run_source_cage_joint_batch.py"
        spec = importlib.util.spec_from_file_location("joint_batch_runner", path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        self.assertTrue(hasattr(module, "independent_top_constraints"),
                        "constraints must be recomputed from exported geometry")
        from make_it_stand.verification.mass_moments import TriangleMesh, compute_mass_properties
        vertices = [(float(i & 1) + .1, float((i >> 1) & 1), float((i >> 2) & 1)) for i in range(8)]
        faces = [(0,2,1),(1,2,3),(4,5,6),(5,7,6),(0,1,4),(1,5,4),
                 (2,6,3),(3,6,7),(0,4,2),(2,4,6),(1,3,5),(3,7,5)]
        properties = compute_mass_properties(TriangleMesh(vertices, faces), 1.3)
        constraints = module.independent_top_constraints(properties, 0.)
        self.assertAlmostEqual(constraints[0], .1, places=13)
        self.assertTrue(all(abs(v) < 1e-13 for v in constraints[1:]))

    def test_checker_rejects_a_claim_without_admitted_material(self):
        root = Path(__file__).resolve().parents[2]
        path = root / "experiments/spin_it_reproduction/strict/run_source_cage_joint_batch.py"
        self.assertTrue(path.is_file(), "batch checker must exist before it can grant validation")
        spec = importlib.util.spec_from_file_location("joint_batch_runner", path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            with (output / "result.json").open("w") as stream:
                json.dump({"binary_status": "BINARY_MATERIAL_VALIDATED_PRE_SEAL",
                           "binary_candidate_available": False, "passed": True}, stream)
            with self.assertRaises(ValueError):
                module.check_material_output(output)

    def test_preflight_freezes_nontrivial_input_without_optimizing(self):
        root = Path(__file__).resolve().parents[2]
        binary = Path(os.environ.get("SPIN_JOINT_BATCH_BINARY", str(
            root / "build/strict-algorithms-exact-current-20260922/"
            "spin-it-2014-source-cage-asymmetric-joint-regression")))
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "preflight"
            completed = subprocess.run(
                [str(binary), "--batch-preflight", str(output)],
                capture_output=True, text=True, timeout=180)
            self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
            self.assertTrue((output / "input.json").is_file(),
                            "native runner must publish pre-solve input, not run the short regression")
            data = json.loads((output / "input.json").read_text())
            self.assertEqual(data["classification"], "SYNTHETIC_REGRESSION")
            self.assertGreater(data["initial_constraint_inf"], 1e-5)
            self.assertEqual(data["handle_count"], 98)
            self.assertGreater(data["beta_count"], 0)
            self.assertEqual(len(data["initial_state"]), data["beta_count"] + 1 + 4 * 98)
            self.assertFalse((output / "result.json").exists())
            original = (output / "input.json").read_bytes()
            retry = subprocess.run([str(binary), "--batch-preflight", str(output)],
                                   capture_output=True, text=True, timeout=180)
            self.assertNotEqual(retry.returncode, 0, "existing run must not be overwritten")
            self.assertEqual((output / "input.json").read_bytes(), original)


if __name__ == "__main__":
    unittest.main()
