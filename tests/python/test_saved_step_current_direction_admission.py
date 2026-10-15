"""The supplementary mode must enforce admission before opening a fixture."""
import subprocess
import sys
import unittest


class AdmissionTest(unittest.TestCase):
    def test_current_direction_mode_rejects_unbounded_execution_before_input(self):
        result = subprocess.run(
            [sys.argv[1], "/nonexistent-frozen-fixture.json", "--actual-current-direction"],
            capture_output=True, text=True, timeout=10,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertNotIn("usage:", result.stderr)
        self.assertNotIn("missing frozen fixture", result.stderr)
        self.assertTrue(any(message in result.stderr for message in (
            "finite memory limit required", "requires memory <= 24 GiB and swap = 0",
            "cannot read cgroup resource limit", "unambiguous cgroup v2 path required",
        )), result.stderr)


if __name__ == "__main__":
    unittest.main(argv=[sys.argv[0]])
