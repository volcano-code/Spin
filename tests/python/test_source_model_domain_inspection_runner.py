from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

from experiments.spin_it_reproduction.strict.run_source_model_domain_inspection import run


class SourceModelDomainInspectionRunnerTests(unittest.TestCase):
    def test_preflight_failure_is_preserved_as_a_failure_manifest(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            mesh = root / "invalid.mesh"
            mesh.write_text("not a Medit mesh\n", encoding="ascii")
            output = root / "run"

            with self.assertRaisesRegex(ValueError, "SHA-256 does not match"):
                run(output, mesh_path=mesh, expected_sha256="0" * 64,
                    build_directory=root, initial_level=1, maximum_depth=1,
                    maximum_leaf_count=100, dependency_roots={})

            failure = json.loads((output / "failure.json").read_text(encoding="ascii"))
            self.assertEqual(failure["status"], "FAILED")
            self.assertEqual(failure["stage"], "RUNNER")
            self.assertIn("SHA-256 does not match", failure["error"])


if __name__ == "__main__":
    unittest.main()
