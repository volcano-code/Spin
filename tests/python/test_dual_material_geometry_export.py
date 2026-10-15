import json
import importlib.util
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

from make_it_stand.verification.mass_moments import TriangleMesh, compute_mass_properties


def read_obj(path):
    vertices, faces = [], []
    for line in Path(path).read_text().splitlines():
        values = line.split()
        if values and values[0] == "v":
            vertices.append(tuple(map(float, values[1:])))
        elif values and values[0] == "f":
            faces.append(tuple(int(x) - 1 for x in values[1:]))
    return TriangleMesh(vertices, faces)


class DualMaterialGeometryTests(unittest.TestCase):
    def test_separate_exported_materials_reproduce_solver_moments(self):
        root = Path(__file__).resolve().parents[2]
        binary = Path(os.environ.get("SPIN_DUAL_EXPORT_BINARY", str(root /
            "build/strict-algorithms-exact-current-20260922/spin-it-2014-dual-material-regression")))
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "exports"
            process = subprocess.run([str(binary), "--export-directory", str(output)],
                                     capture_output=True, text=True, timeout=180, cwd=root)
            self.assertEqual(process.returncode, 0, process.stdout + process.stderr)
            self.assertTrue((output / "dense_regression/manifest.json").is_file(),
                            "the two actual material boundaries must be emitted")
            path = root / "experiments/spin_it_reproduction/strict/run_spin2014_dual_material.py"
            spec = importlib.util.spec_from_file_location("dual_runner", path)
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            self.assertTrue(hasattr(module, "check_exported_geometry"), "strict runner needs independent mesh checks")
            check = module.check_exported_geometry(output)
            self.assertEqual(check["sumone_reconstructed"]["status"], "NO_RETAINED_DENSE_INSERT")
            for case in ("dense_regression", "sumone_reconstructed"):
                data = json.loads((output / case / "manifest.json").read_text())
                self.assertEqual(len(data["input_material"]), len(data["labels"]))
                self.assertEqual(len(data["input_carveable"]), len(data["labels"]))
                self.assertEqual(check[case]["partition"]["label_partition_status"],
                                 "DISJOINT_VOXEL_PARTITION_VALIDATED")
                self.assertLessEqual(check[case]["constraints"]["recomputation_error"], 1e-10)
                combined = [0.] * 10
                for entry in data["materials"]:
                    if entry["cell_count"] == 0:
                        self.assertIsNone(entry["mesh"])
                        # Native subtraction may leave roundoff for a truly
                        # empty label. Preserve it and apply the fixed gate.
                        self.assertLessEqual(max(map(abs, entry["moments"])), 1e-10)
                        continue
                    region = next(item for item in check[case]["regions"]
                                  if item["label"] == entry["label"])
                    self.assertEqual(region["label_surface_status"],
                                     "EXACT_TRIANGLE_MULTISET_MATCH")
                    props = compute_mass_properties(read_obj(output / case / entry["mesh"]), entry["density"])
                    q = props.second_moment
                    moments = [props.mass, *props.first_moment, q[0][0], q[1][1], q[2][2], q[1][2], q[0][2], q[0][1]]
                    self.assertLessEqual(max(abs(a-b) for a,b in zip(moments,entry["moments"])), 1e-10)
                    combined = [a+b for a,b in zip(combined,moments)]
                self.assertLessEqual(max(abs(a-b) for a,b in zip(combined,data["total_moments"])), 1e-10)
                if case == "sumone_reconstructed":
                    self.assertEqual(data["materials"][1]["cell_count"], 0)
                    self.assertEqual(data["strict_result"], "NO_RETAINED_DENSE_INSERT")
            self.assertGreater(
                check["dense_regression"]["partition"]["shared_material_interface_area"], 0.0,
                "the disjoint material regions should meet at a nonzero-area interface",
            )
            self.assertEqual(
                check["sumone_reconstructed"]["partition"]["shared_material_interface_area"], 0.0,
                "an empty dense insert has no shared material interface",
            )
            mesh_path = output / "dense_regression/low_density_shell.obj"
            text = mesh_path.read_text()
            with mesh_path.open("w") as stream:
                stream.write(text.replace("v -1 ", "v -1.1 ", 1))
            with self.assertRaises(ValueError):
                module.check_exported_geometry(output)

    def test_exported_region_mesh_must_match_its_actual_label_cells(self):
        root = Path(__file__).resolve().parents[2]
        binary = Path(os.environ.get("SPIN_DUAL_EXPORT_BINARY", str(root /
            "build/strict-algorithms-exact-current-20260922/spin-it-2014-dual-material-regression")))
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "exports"
            process = subprocess.run([str(binary), "--export-directory", str(output)],
                                     capture_output=True, text=True, timeout=180, cwd=root)
            self.assertEqual(process.returncode, 0, process.stdout + process.stderr)
            source = output / "dense_regression/low_density_shell.obj"
            target = output / "dense_regression/high_density_insert.obj"
            target.write_bytes(source.read_bytes())
            path = root / "experiments/spin_it_reproduction/strict/run_spin2014_dual_material.py"
            spec = importlib.util.spec_from_file_location("dual_runner_mesh_label_mismatch", path)
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            with self.assertRaisesRegex(ValueError, "does not match label-derived voxel surface"):
                module.check_exported_geometry(output)


if __name__ == "__main__":
    unittest.main()
