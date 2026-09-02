from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY_ROOT / "scripts"))

import create_pre_spin_checkpoint as checkpoint


class PreSpinCheckpointTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name) / "repository"
        self.root.mkdir()
        self._git("init", "-b", "main")
        self._git("config", "user.name", "Fixture User")
        self._git("config", "user.email", "fixture@example.com")
        (self.root / ".gitignore").write_text(
            "/experiments/paper_sim_v1/results*/\n",
            encoding="utf-8",
        )
        (self.root / "tracked.txt").write_text("before\n", encoding="utf-8")
        self._git("add", ".gitignore", "tracked.txt")
        self._git("commit", "-m", "fixture baseline")
        (self.root / "tracked.txt").write_text("after\n", encoding="utf-8")
        (self.root / "untracked.txt").write_text("new\n", encoding="utf-8")
        ignored_result = (
            self.root
            / "experiments"
            / "paper_sim_v1"
            / "results_current"
            / "result.bin"
        )
        ignored_result.parent.mkdir(parents=True)
        ignored_result.write_bytes(b"preserve me\n")
        self.output = self.root / "evidence" / "pre-spin"

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def _git(self, *arguments: str) -> str:
        return subprocess.run(
            ["git", *arguments],
            cwd=self.root,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        ).stdout.strip()

    def test_checkpoint_records_real_diff_and_sorted_untracked_hashes(self) -> None:
        manifest = checkpoint.create_checkpoint(self.root, self.output)

        self.assertEqual(manifest["schema"], "spin-it-pre-checkpoint-v1")
        self.assertEqual(manifest["git"]["branch"], "main")
        self.assertEqual(manifest["git"]["user_name"], "Fixture User")
        self.assertEqual(manifest["git"]["user_email"], "fixture@example.com")

        patch_bytes = (self.output / "tracked.patch").read_bytes()
        self.assertIn(b"-before", patch_bytes)
        self.assertIn(b"+after", patch_bytes)
        self.assertEqual(
            manifest["tracked_patch_sha256"],
            hashlib.sha256(patch_bytes).hexdigest(),
        )

        untracked = manifest["untracked_files"]
        self.assertEqual(
            untracked,
            sorted(untracked, key=lambda item: item["path"].encode("utf-8")),
        )
        self.assertEqual(
            untracked,
            [
                {
                    "path": "experiments/paper_sim_v1/results_current/result.bin",
                    "size": 12,
                    "sha256": hashlib.sha256(b"preserve me\n").hexdigest(),
                },
                {
                    "path": "untracked.txt",
                    "size": 4,
                    "sha256": hashlib.sha256(b"new\n").hexdigest(),
                },
            ],
        )
        self.assertNotIn("evidence/pre-spin/manifest.json", {
            item["path"] for item in untracked
        })

        persisted = json.loads(
            (self.output / "manifest.json").read_text(encoding="utf-8")
        )
        self.assertEqual(persisted, manifest)
        sums = (self.output / "SHA256SUMS").read_text(encoding="utf-8")
        self.assertIn("  tracked.patch\n", sums)
        self.assertIn("  manifest.json\n", sums)

    def test_checkpoint_rejects_output_outside_repository(self) -> None:
        outside = Path(self.temporary.name) / "outside"

        with self.assertRaises(ValueError):
            checkpoint.create_checkpoint(self.root, outside)

    def test_repository_ignores_legacy_raw_results_but_not_reports(self) -> None:
        ignored = subprocess.run(
            [
                "git",
                "check-ignore",
                "--no-index",
                "--quiet",
                "experiments/paper_sim_v1/results_current/run/result.json",
            ],
            cwd=REPOSITORY_ROOT,
            check=False,
        )
        report = subprocess.run(
            [
                "git",
                "check-ignore",
                "--no-index",
                "--quiet",
                "experiments/paper_sim_v1/REPORT.md",
            ],
            cwd=REPOSITORY_ROOT,
            check=False,
        )
        new_evidence = subprocess.run(
            [
                "git",
                "check-ignore",
                "--no-index",
                "--quiet",
                "experiments/spin_it_reproduction/reachability/result.json",
            ],
            cwd=REPOSITORY_ROOT,
            check=False,
        )

        self.assertEqual(ignored.returncode, 0)
        self.assertEqual(report.returncode, 1)
        self.assertEqual(new_evidence.returncode, 1)


if __name__ == "__main__":
    unittest.main()
