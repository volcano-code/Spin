"""Behavioral contracts for observed provenance, using real temporary Git trees."""

from __future__ import annotations

import dataclasses
import hashlib
import importlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch


class RunEvidenceV2Tests(unittest.TestCase):
    def setUp(self):
        name = "make_it_stand.verification.run_evidence"
        self.assertIsNotNone(importlib.util.find_spec(name), "observed V2 evidence module is missing")
        self.api = importlib.import_module(name)
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name) / "repo"
        self.root.mkdir()
        for name, content in {
            "src/solver.cpp": "int value = 1;\n",
            "include/model.h": "#pragma once\n",
            "CMakeLists.txt": "project(fixture)\n",
            "python/pkg/__init__.py": "VALUE = 1\n",
            "experiments/runner.py": "print('run')\n",
            "experiments/checker.py": "print('check')\n",
            "data/asset.bin": "asset bytes\n",
            "config/build.json": '{"optimization":"O2"}\n',
            "reports/final.json": '{"mass":1.0}\n',
        }.items():
            self.write(name, content)
        self.dep = Path(self.temp.name) / "dependency"
        self.dep.mkdir()
        (self.dep / "VERSION").write_text("1.2.3\n")
        (self.dep / "source.h").write_text("int dependency;\n")
        self.git("init", "-q")
        self.git("add", ".")
        self.commit("fixture")
        self.params = {"density": 1.0, "resolution": [2, 3]}
        self.seed = {"seed": 7, "protocol": "single deterministic fixture seed"}
        self.events = [
            {"iteration": 0, "status": "failed", "data": {"reason": "trial rejected"}},
            {"iteration": 1, "status": "ok", "data": {"energy": 0.25}},
        ]

    def write(self, name, content):
        path = self.root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content)
        return path

    def git(self, *args):
        return subprocess.check_output(["git", "-C", str(self.root), *args], stderr=subprocess.PIPE).decode().strip()

    def commit(self, message):
        self.git("-c", "user.name=Evidence Fixture", "-c", "user.email=fixture@example.invalid", "commit", "-qm", message)

    def inputs(self, **changes):
        values = dict(
            repository_root=self.root,
            runner_files=[self.root / "experiments/runner.py"],
            checker_files=[self.root / "experiments/checker.py"],
            asset_files=[self.root / "data/asset.bin"],
            build_configuration_files=[self.root / "config/build.json"],
            parameters=self.params, seed_protocol=self.seed,
            dependency_roots={"fixture-dependency": self.dep},
            native_modules=[], pure_python=True,
        )
        values.update(changes)
        return values

    def begin(self, **changes):
        return self.api.begin_run(**self.inputs(**changes))

    def finish(self, begin=None, **changes):
        values = dict(events=self.events, output_files=[self.root / "reports/final.json"],
                      measurements={"timestamp": "2030-01-01T00:00:00Z", "elapsed_seconds": 1.0})
        values.update(changes)
        return self.api.finish_run(begin or self.begin(), **values)

    def test_valid_capture_is_stable_and_round_trips_without_attestation(self):
        first, second = self.finish(), self.finish()
        self.assertEqual(first.to_json(), second.to_json())
        restored = self.api.PaperRunEvidenceV2.from_json(first.to_json())
        audit = self.api.check_freshness(restored)
        self.assertEqual(audit.status, "CURRENT")
        self.assertEqual(audit.build_source_attestation, "NOT_VERIFIED")
        self.assertEqual(audit.physical_validation, "NOT_IN_SCOPE")
        self.assertEqual(restored.as_dict()["build_source_attestation"], "NOT_VERIFIED")
        observation = first.as_dict()["begin"]["observation"]
        self.assertEqual(observation["source"]["git"]["head"], self.git("rev-parse", "HEAD"))
        self.assertEqual(observation["runtime"]["executable"]["path"], str(Path(sys.executable).resolve()))
        self.assertIn("disk", observation["runtime"]["measurement_scope"])
        self.assertEqual(observation["dependencies"]["coverage"], "DECLARED_ROOTS_ONLY")
        self.assertIn("1.2.3", str(observation["dependencies"]))

    def test_source_edit_add_delete_and_untracked_changes_are_stale(self):
        for action in ("edit", "add", "delete", "untracked_edit", "untracked_delete"):
            with self.subTest(action=action):
                extra = self.write("src/new.cpp", "old\n") if action.startswith("untracked") else None
                evidence = self.finish()
                path = self.root / "src/solver.cpp"
                original = path.read_text()
                if action == "edit":
                    path.write_text("changed\n")
                elif action == "add":
                    extra = self.write("src/new.cpp", "new\n")
                elif action == "delete":
                    path.unlink()
                elif action == "untracked_edit":
                    extra.write_text("changed\n")
                else:
                    extra.unlink()
                result = self.api.check_freshness(evidence)
                self.assertEqual(result.status, "STALE_RESULT")
                self.assertTrue(any("src/" in item for item in result.mismatches), result)
                path.write_text(original)
                if extra and extra.exists():
                    extra.unlink()

    def test_staged_and_unstaged_binary_diffs_are_observed(self):
        self.write("src/solver.cpp", "staged\n")
        self.git("add", "src/solver.cpp")
        self.write("src/solver.cpp", "unstaged\n")
        git = self.begin().as_dict()["observation"]["source"]["git"]
        import base64
        for key, args in (("head_to_index", ["--cached", "HEAD"]),
                          ("index_to_worktree", []), ("head_to_worktree", ["HEAD"])):
            observed = base64.b64decode(git[key + "_diff_base64"])
            expected = subprocess.check_output(["git", "-C", str(self.root), "diff", "--binary", "--no-ext-diff", "--no-textconv", *args, "--", "src/solver.cpp"])
            self.assertEqual(observed, expected)
            self.assertTrue(observed)

    def test_staged_addition_deleted_from_worktree_remains_in_git_diff(self):
        path = self.write("src/staged-only.cpp", "int staged_only;\n")
        self.git("add", "src/staged-only.cpp")
        path.unlink()
        context = self.begin().as_dict()["observation"]["source"]["git"]
        import base64
        self.assertIn("src/staged-only.cpp", context["diff_paths"])
        self.assertIn(b"int staged_only;", base64.b64decode(context["head_to_index_diff_base64"]))
        self.assertIn(b"int staged_only;", base64.b64decode(context["index_to_worktree_diff_base64"]))

    def test_source_policy_covers_effective_code_build_and_test_additions(self):
        for name in ("bindings/new.cpp", "cmake/options.cmake", "python/pkg/new.py",
                     "scripts/check.sh", "tests/fixture.json", "experiments/new/run.py",
                     "CMakeUserPresets.json"):
            with self.subTest(name=name):
                evidence = self.finish()
                path = self.write(name, "new effective input\n")
                self.assertEqual(self.api.check_freshness(evidence).status, "STALE_RESULT")
                path.unlink()

    def test_report_cache_and_build_additions_do_not_change_source_digest(self):
        before = self.finish()
        for name in ("reports/new.md", "python/pkg/__pycache__/temp.pyc", "build/temp.cpp", "experiments/run/reports/generated.py"):
            self.write(name, "ignored generated data\n")
        after = self.finish()
        self.assertEqual(before.as_dict()["begin"]["observation"]["source"]["content_sha256"],
                         after.as_dict()["begin"]["observation"]["source"]["content_sha256"])
        self.assertEqual(self.api.check_freshness(before).status, "CURRENT")

    def test_report_only_commit_is_contextual_revision_change(self):
        evidence = self.finish()
        self.write("reports/note.md", "new report\n")
        self.git("add", "reports/note.md")
        self.commit("report only")
        result = self.api.check_freshness(evidence)
        self.assertEqual(result.status, "CURRENT")
        self.assertTrue(result.contextual_changes)
        self.assertIn("head", " ".join(result.contextual_changes).lower())

    def test_input_file_edits_and_disappearance_are_stale(self):
        paths = [self.dep / "source.h", self.dep / "VERSION", self.root / "data/asset.bin",
                 self.root / "experiments/runner.py", self.root / "experiments/checker.py",
                 self.root / "config/build.json", self.root / "reports/final.json"]
        for path in paths:
            for change in ("edit", "delete"):
                with self.subTest(path=path, change=change):
                    evidence = self.finish()
                    original = path.read_bytes()
                    path.write_bytes(b"changed\n") if change == "edit" else path.unlink()
                    result = self.api.check_freshness(evidence)
                    self.assertEqual(result.status, "STALE_RESULT")
                    self.assertIn(path.name, " ".join(result.mismatches))
                    path.write_bytes(original)

    def test_dependency_additions_are_detected_and_siblings_are_out_of_scope(self):
        evidence = self.finish()
        (self.dep.parent / "outside.h").write_text("outside declared root")
        self.assertEqual(self.api.check_freshness(evidence).status, "CURRENT")
        (self.dep / "added.h").write_text("new dependency content")
        self.assertEqual(self.api.check_freshness(evidence).status, "STALE_RESULT")

    def test_end_capture_keeps_midrun_mismatch_even_after_file_is_restored(self):
        begin = self.begin()
        path = self.root / "data/asset.bin"
        original = path.read_text()
        path.write_text("changed during run")
        evidence = self.finish(begin)
        path.write_text(original)
        self.assertEqual(self.api.check_freshness(evidence).status, "STALE_RESULT")

    def test_caller_and_returned_data_mutation_cannot_change_evidence(self):
        begin = self.begin()
        self.params["resolution"].append(99)
        self.seed["seed"] = 99
        evidence = self.finish(begin)
        saved = evidence.to_json()
        self.events[0]["data"]["reason"] = "mutated"
        data = evidence.as_dict()
        data["begin"]["request"]["parameters"]["resolution"].append(100)
        self.assertEqual(evidence.to_json(), saved)
        self.assertEqual(data["begin"]["request"]["seed_protocol"]["seed"], 7)
        self.assertEqual(evidence.as_dict()["begin"]["request"]["parameters"]["resolution"], [2, 3])
        with self.assertRaises((dataclasses.FrozenInstanceError, AttributeError)):
            evidence.payload_json = "{}"

    def test_missing_legacy_tampered_or_empty_evidence_cannot_be_current(self):
        valid = self.finish().as_dict()
        tampered = json.loads(json.dumps(valid))
        tampered["scientific"]["events"][0]["data"]["reason"] = "tampered"
        for bad in (None, {}, "{}", {"schema": "paper-evidence-v1"}, {"passed": True},
                    {**valid, "scientific": {}}, tampered, {**valid, "verified": True}):
            with self.subTest(bad=str(bad)[:100]):
                self.assertNotEqual(self.api.check_freshness(bad).status, "CURRENT")
        for raw in ('{"schema":1,"schema":2}', '{"x":NaN}'):
            with self.assertRaises(ValueError):
                self.api.PaperRunEvidenceV2.from_json(raw)

    def test_required_inputs_cannot_be_empty_or_missing(self):
        for key, bad in (("runner_files", []), ("checker_files", []), ("asset_files", []),
                         ("build_configuration_files", []), ("parameters", {}),
                         ("seed_protocol", {}), ("dependency_roots", {}),
                         ("runner_files", [self.root / "missing.py"]),
                         ("dependency_roots", {"missing": self.dep / "missing"})):
            with self.subTest(key=key, bad=bad), self.assertRaises(ValueError):
                self.begin(**{key: bad})

    def test_nonfinite_non_json_and_boolean_seeds_are_rejected(self):
        for bad in (float("nan"), float("inf"), -float("inf"), object(), {1, 2}):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                self.begin(parameters={"bad": bad})
        for bad in (True, False, 1.5, "1"):
            with self.subTest(seed=bad), self.assertRaises(ValueError):
                self.begin(seed_protocol={"seed": bad, "protocol": "fixture"})

    def test_duplicates_symlinks_and_overlapping_dependency_roots_are_rejected(self):
        asset = self.root / "data/asset.bin"
        for paths in ([asset, asset], [asset, asset.parent / "." / asset.name]):
            with self.assertRaises(ValueError):
                self.begin(asset_files=paths)
        alias = self.root / "data/link.bin"
        alias.symlink_to(asset)
        with self.assertRaises(ValueError):
            self.begin(asset_files=[alias])
        source_link = self.root / "src/linked.cpp"
        source_link.symlink_to(self.dep / "source.h")
        with self.assertRaises(ValueError):
            self.begin()
        source_link.unlink()
        with self.assertRaises(ValueError):
            self.begin(dependency_roots={"a": self.dep, "b": self.dep})

    def test_git_metadata_is_required(self):
        shutil.rmtree(self.root / ".git")
        with self.assertRaises(ValueError):
            self.begin()

    def test_parent_segments_cannot_hide_a_symlink_input(self):
        alias = self.root / "alias"
        alias.symlink_to(self.dep, target_is_directory=True)
        # abspath would erase alias/.. before checking for the symlink.
        with self.assertRaises(ValueError):
            self.begin(asset_files=[alias / ".." / "data/asset.bin"])

    def test_native_modules_must_be_loaded_native_objects(self):
        import _json
        evidence = self.finish(self.begin(native_modules=[_json], pure_python=False))
        self.assertEqual(self.api.check_freshness(evidence).status, "CURRENT")
        for modules in (["not_loaded_module"], ["/tmp/fake.so"], [json], [object()], [_json, _json]):
            with self.subTest(modules=modules), self.assertRaises(ValueError):
                self.begin(native_modules=modules, pure_python=False)
        with self.assertRaises(ValueError):
            self.begin(native_modules=[], pure_python=False)
        with self.assertRaises(ValueError):
            self.begin(native_modules=[_json], pure_python=True)
        with patch.object(_json, "__file__", str(self.root / "missing.so")):
            with self.assertRaises(ValueError):
                self.begin(native_modules=[_json], pure_python=False)
        with patch.object(sys, "executable", ""):
            with self.assertRaises(ValueError):
                self.begin()

    def test_runtime_file_changes_are_stale_in_actual_copied_interpreter(self):
        # Exercise actual launched/loaded binaries. Replace paths after loading;
        # never edit the host interpreter or a mapped shared object's bytes.
        import _json
        executable = Path(self.temp.name) / "python-copy"
        module = Path(self.temp.name) / Path(_json.__file__).name
        shutil.copy2(sys.executable, executable)
        shutil.copy2(_json.__file__, module)
        script = """
import importlib.util, json, pathlib, shutil, sys
from make_it_stand.verification.run_evidence import begin_run, finish_run, check_freshness
module_path, request_json = sys.argv[1:]
spec = importlib.util.spec_from_file_location('_json', module_path)
native = importlib.util.module_from_spec(spec)
sys.modules['_json'] = native
spec.loader.exec_module(native)
request = json.loads(request_json)
request['native_modules'] = [native]
begin = begin_run(**request)
evidence = finish_run(begin, events=[{'iteration':0,'status':'ok','data':{'mass':1}}],
                      output_files=[str(pathlib.Path(request['repository_root'])/'reports/final.json')],
                      measurements={'elapsed_seconds':0.1})
for path in (pathlib.Path(module_path), pathlib.Path(sys.executable)):
    old = path.with_name(path.name + '.old')
    path.rename(old)
    shutil.copy2(old, path)
    with path.open('ab') as stream:
        stream.write(b'changed on disk')
    result = check_freshness(evidence)
    assert result.status == 'STALE_RESULT', result
    assert path.name in ' '.join(result.mismatches), result
    path.unlink()
    old.rename(path)
print('runtime replacements detected')
"""
        inputs = self.inputs(pure_python=False)
        for key in ("repository_root",):
            inputs[key] = str(inputs[key])
        for key in ("runner_files", "checker_files", "asset_files", "build_configuration_files"):
            inputs[key] = list(map(str, inputs[key]))
        inputs["dependency_roots"] = {k: str(v) for k, v in inputs["dependency_roots"].items()}
        completed = subprocess.run([str(executable), "-c", script, str(module), json.dumps(inputs)],
                                   text=True, capture_output=True)
        self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
        self.assertIn("runtime replacements detected", completed.stdout)

    def test_event_order_failure_and_output_content_affect_scientific_digest(self):
        original = self.finish()
        reordered = self.finish(events=list(reversed(self.events)))
        no_failure = self.finish(events=self.events[1:])
        altered = json.loads(json.dumps(self.events))
        altered[0]["status"] = "ok"
        changed_failure = self.finish(events=altered)
        self.write("reports/final.json", '{"mass":2.0}\n')
        changed_output = self.finish()
        digests = [item.as_dict()["scientific"]["determinism_sha256"]
                   for item in (original, reordered, no_failure, changed_failure, changed_output)]
        self.assertEqual(len(set(digests)), 5)

    def test_timestamps_and_resource_measurements_are_preserved_outside_digest(self):
        first = self.finish()
        measurements = {"timestamp": "2031-02-03T04:05:06Z", "elapsed_seconds": 17.0, "peak_rss_bytes": 300}
        second = self.finish(measurements=measurements)
        self.assertEqual(first.as_dict()["scientific"]["determinism_sha256"],
                         second.as_dict()["scientific"]["determinism_sha256"])
        self.assertEqual(second.as_dict()["measurements"], measurements)

    def test_empty_malformed_nonfinite_scientific_results_fail_closed(self):
        for changes in ({"events": []}, {"events": [{}]}, {"output_files": []},
                        {"events": [{"iteration": True, "status": "ok", "data": {"x": 1}}]},
                        {"events": [{"iteration": 0, "status": "ok", "data": {"x": float("nan")}}]},
                        {"events": [{"iteration": 0, "status": "ok", "data": {}}]},
                        {"measurements": {"elapsed_seconds": float("inf")}}):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                self.finish(**changes)
        self.write("reports/final.json", "")
        with self.assertRaises(ValueError):
            self.finish()

    def test_json_scientific_output_must_be_finite_nonempty_and_well_formed(self):
        for content in ('{"mass":NaN}', '{"mass":Infinity}', '{"mass":1,"mass":2}',
                        'invalid json', '{}', '[]', 'null'):
            with self.subTest(content=content):
                self.write("reports/final.json", content)
                with self.assertRaises(ValueError):
                    self.finish()

    def test_nested_schema_cannot_hide_boolean_counters_or_fake_git_observations(self):
        def seal(value):
            value.pop("integrity_sha256", None)
            payload = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False)
            value["integrity_sha256"] = hashlib.sha256(payload.encode()).hexdigest()
        for kind in ("runtime_boolean", "dependency_revision", "duplicate_diff_path", "version_status"):
            with self.subTest(kind=kind):
                value = self.finish().as_dict()
                for stage in ("begin", "end"):
                    obs = value[stage]["observation"]
                    if kind == "runtime_boolean":
                        obs["runtime"]["pure_python"] = 1
                    elif kind == "dependency_revision":
                        obs["dependencies"]["roots"]["fixture-dependency"]["git_head"] = True
                    elif kind == "duplicate_diff_path":
                        obs["source"]["git"]["diff_paths"].append("src/solver.cpp")
                    else:
                        obs["dependencies"]["roots"]["fixture-dependency"]["version_observation_status"] = "UNAVAILABLE"
                    seal(value[stage])
                seal(value)
                self.assertEqual(self.api.check_freshness(value).status, "INVALID_EVIDENCE")


if __name__ == "__main__":
    unittest.main()
