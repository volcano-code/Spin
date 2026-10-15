"""Negative safety regressions for the retained Cherchi probe wrapper."""

from __future__ import annotations

from pathlib import Path
import subprocess
import sys
import unittest

from test_probe import box

EXECUTABLE = Path(sys.argv[1]).resolve()
BASE = Path(sys.argv[2]).resolve()


class SafetyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.executable = EXECUTABLE
        cls.base = BASE
        cls.base.mkdir(parents=True, exist_ok=False)

    def test_existing_workdir_is_not_reused(self) -> None:
        work = self.base / "existing-workdir"
        work.mkdir()
        sentinel = work / "sentinel.txt"
        sentinel.write_text("retain this evidence", encoding="ascii")
        result = subprocess.run((sys.executable, str(Path(__file__).with_name("test_probe.py")),
                                 str(self.executable), str(work)), capture_output=True, text=True, timeout=120)
        self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(sentinel.read_text(encoding="ascii"), "retain this evidence")
        self.assertEqual(sorted(path.name for path in work.iterdir()), ["sentinel.txt"])

    def test_existing_output_is_not_overwritten(self) -> None:
        work = self.base / "existing-output"
        work.mkdir()
        left, right, output = work / "left.obj", work / "right.obj", work / "result.obj"
        box(left, (0, 0, 0), (3, 2, 1.5))
        box(right, (1, .5, .25), (4, 2.5, 1.25))
        output.write_text("retain this output", encoding="ascii")
        result = subprocess.run((str(self.executable), "intersection", str(left), str(right), str(output)),
                                capture_output=True, text=True, timeout=120)
        self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(output.read_text(encoding="ascii"), "retain this output")

    def test_unsupported_operation_is_rejected_without_output(self) -> None:
        work = self.base / "unsupported"
        work.mkdir()
        output = work / "result.obj"
        result = subprocess.run((str(self.executable), "union", "not-needed.obj", "also-not-needed.obj", str(output)),
                                capture_output=True, text=True, timeout=120)
        self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
        self.assertFalse(output.exists())


if __name__ == "__main__":
    sys.argv = [sys.argv[0]]
    unittest.main(verbosity=2)
