"""Real GNU compiler fixtures; none of these fixtures attest linking or science."""

from concurrent.futures import ThreadPoolExecutor
from contextlib import redirect_stderr, redirect_stdout
import hashlib
import importlib
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "python/make_it_stand/verification/compiler_observer.py"


class CompilerObserverTests(unittest.TestCase):
    def setUp(self):
        self.assertTrue(SCRIPT.is_file(), "compiler observer is missing")
        self.api = importlib.import_module("make_it_stand.verification.compiler_observer")
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.cc = shutil.which("cc")
        self.assertIsNotNone(self.cc, "real GNU cc is required; do not skip fixtures")
        self.write("main.c", '#include <stddef.h>\n#include "outer.h"\nint main(void) { return VALUE != 42; }\n')
        self.write("outer.h", '#include "generated.h"\n')
        self.write("generated.h", '#define VALUE 42\n')
        self.evidence = self.root / "evidence"

    def write(self, name, text):
        path = self.root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
        return path

    def argv(self, output="main.o", extra=()):
        return [self.cc, *extra, "-c", "main.c", "-o", output]

    def launch(self, args=None, timeout="10", extra=(), env=None):
        return subprocess.run([sys.executable, str(SCRIPT), "--evidence-dir", str(self.evidence),
                               "--timeout", timeout, *extra, "--", *(args or self.argv())],
                              cwd=self.root, env=env, capture_output=True, text=True, timeout=20)

    def receipts(self):
        return [json.loads(p.read_text()) for p in sorted(self.evidence.glob("*/receipt.json"))]

    def test_real_object_and_runnable_link_include_transitive_system_generated_inputs(self):
        result = self.launch()
        self.assertEqual(result.returncode, 0, result.stderr)
        receipt, = self.receipts()
        self.assertEqual(receipt["status"], "OBSERVED_SUCCESS")
        self.assertEqual(receipt["build_source_attestation"], "NOT_VERIFIED")
        self.assertEqual(receipt["compiler_exit_status"], 0)
        self.assertEqual(receipt["argv"], self.argv())
        self.assertEqual(receipt["cwd"], str(self.root))
        paths = [p["logical_path"] for p in receipt["inputs_before"]]
        for name in ("main.c", "outer.h", "generated.h"):
            self.assertIn(name, paths)
        self.assertTrue(any(p.endswith("stddef.h") for p in paths))
        self.assertEqual(receipt["object"]["sha256"], hashlib.sha256((self.root / "main.o").read_bytes()).hexdigest())
        self.assertIn("Free Software Foundation", receipt["compiler_before"]["version"])
        subprocess.run([self.cc, "main.o", "-o", "tiny"], cwd=self.root, check=True, capture_output=True)
        self.assertEqual(subprocess.run([str(self.root / "tiny")]).returncode, 0)

    def test_real_cmake_lto_flags_produce_observed_object_without_claiming_link_attestation(self):
        cxx = shutil.which("c++")
        self.assertIsNotNone(cxx, "real GNU c++ is required")
        self.write("main.cpp", '#include "outer.h"\nint main() { return VALUE != 42; }\n')
        args = [cxx, "-O3", "-flto=auto", "-fno-fat-lto-objects", "-std=gnu++17",
                "-c", "main.cpp", "-o", "main-lto.o"]
        result = self.launch(args)
        self.assertEqual(result.returncode, 0, result.stderr)
        receipt, = self.receipts()
        self.assertEqual(receipt["status"], "OBSERVED_SUCCESS")
        self.assertEqual(receipt["argv"], args)
        self.assertIn("generated.h", [row["logical_path"] for row in receipt["inputs_before"]])
        self.assertEqual(receipt["object"]["sha256"],
                         hashlib.sha256((self.root / "main-lto.o").read_bytes()).hexdigest())
        self.assertEqual(receipt["link_coverage"], "NOT_VERIFIED")
        self.assertEqual(receipt["build_source_attestation"], "NOT_VERIFIED")
        # The recorded object must actually link and execute; this separate
        # test operation does not upgrade the compiler-only receipt's scope.
        subprocess.run([cxx, "-flto=auto", "main-lto.o", "-o", "lto-program"],
                       cwd=self.root, check=True, capture_output=True)
        self.assertEqual(subprocess.run([str(self.root / "lto-program")]).returncode, 0)

    def test_configured_occt_strip_option_keeps_real_compile_observation(self):
        args = self.argv(output="occt.o", extra=["-O3", "-DNDEBUG", "-s"])
        result = self.launch(args)
        self.assertEqual(result.returncode, 0, result.stderr)
        receipt, = self.receipts()
        self.assertEqual(receipt["argv"], args)
        self.assertEqual(receipt["status"], "OBSERVED_SUCCESS")
        self.assertIn("generated.h", [r["logical_path"] for r in receipt["inputs_before"]])
        subprocess.run([self.cc, "occt.o", "-o", "occt-control"],
                       cwd=self.root, check=True, capture_output=True)
        self.assertEqual(subprocess.run([str(self.root / "occt-control")]).returncode, 0)
        from make_it_stand.verification.compiler_dependencies import parse_command, ObservationError
        with self.assertRaises(ObservationError):
            parse_command(["clang", "--no-default-config", "-x", "c", "-s",
                           "-c", "main.c", "-o", "out.o"], profile="clang-c-v1")

    def test_configured_tbb_gnu_flags_keep_feature_selected_dependency_closure(self):
        cxx = shutil.which("c++")
        self.assertIsNotNone(cxx)
        self.write("rtm-selected.h", "#define RTM_VALUE 17\n")
        self.write("waitpkg-selected.h", "#define WAIT_VALUE 25\n")
        self.write("tbb-flags.cpp", '#if !defined(__RTM__) || !defined(__WAITPKG__)\n'
                   '#error feature flags lost\n#endif\n#include "rtm-selected.h"\n'
                   '#include "waitpkg-selected.h"\nint main(){return RTM_VALUE+WAIT_VALUE!=42;}\n')
        args = [cxx, "-O3", "-flifetime-dse=1", "-mrtm", "-mwaitpkg",
                "-flto=auto", "-fno-fat-lto-objects", "-std=c++20",
                "-c", "tbb-flags.cpp", "-o", "tbb-flags.o"]
        result = self.launch(args)
        self.assertEqual(result.returncode, 0, result.stderr)
        receipt, = self.receipts()
        self.assertEqual(receipt["argv"], args)
        self.assertEqual(receipt["status"], "OBSERVED_SUCCESS")
        for collection in ("inputs_before", "inputs_after"):
            paths = [row["logical_path"] for row in receipt[collection]]
            self.assertIn("rtm-selected.h", paths)
            self.assertIn("waitpkg-selected.h", paths)
        subprocess.run([cxx, "-flto=1", "tbb-flags.o", "-o", "tbb-flags"],
                       cwd=self.root, check=True, capture_output=True)
        self.assertEqual(subprocess.run([str(self.root / "tbb-flags")]).returncode, 0)
        from make_it_stand.verification.compiler_dependencies import parse_command, ObservationError
        for flag in ("-flifetime-dse=2", "-mrtm-unknown", "-mwaitpkg=1"):
            with self.subTest(flag=flag), self.assertRaises(ObservationError):
                parse_command([cxx, flag, "-c", "main.c", "-o", "out.o"])
        for flag in ("-flifetime-dse=1", "-mrtm", "-mwaitpkg"):
            with self.subTest(clang_flag=flag), self.assertRaises(ObservationError):
                parse_command(["clang", "--no-default-config", "-x", "c", flag,
                               "-c", "main.c", "-o", "out.o"], profile="clang-c-v1")

    def test_include_defines_language_and_dependency_flags_preserved_without_probe_overwrites(self):
        self.write("chosen/selected.h", "#define SELECTED 42\n")
        self.write("main.c", '#include HEADER\nint main(void) { return SELECTED != 42; }\n')
        args = self.argv(extra=["-I", "chosen", '-DHEADER="selected.h"', "-x", "c", "-MD", "-MF", "original.d", "-MT", "main.o"])
        result = self.launch(args)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue((self.root / "original.d").read_text().startswith("main.o:"))
        self.assertIn("chosen/selected.h", [p["logical_path"] for p in self.receipts()[0]["inputs_before"]])

    def test_gnu_escaped_dependency_paths(self):
        for dirname in ("space dir", "hash#dir", "dollar$dir", "back\\slash", "mixed\\ space#$dir"):
            with self.subTest(dirname=dirname):
                self.write(dirname + "/outer.h", '#include "generated.h"\n')
                self.write(dirname + "/generated.h", '#define VALUE 42\n')
                self.write("main.c", '#include <outer.h>\nint main(void) { return VALUE != 42; }\n')
                result = self.launch(self.argv(output=dirname + "/main.o", extra=["-I", dirname]))
                self.assertEqual(result.returncode, 0, result.stderr)
                receipt = self.receipts()[-1] if len(self.receipts()) == 1 else next(r for r in self.receipts() if r["argv"][-1] == dirname + "/main.o")
                self.assertIn(dirname + "/generated.h", [p["logical_path"] for p in receipt["inputs_before"]])

    def test_unsupported_arguments_fail_closed(self):
        for extra in (["@args"], ["-fmodules-ts"], ["-include-pch", "x"], ["-fplugin=plugin.so"],
                      ["-Xpreprocessor", "-include"], ["-Wp,-MMD,unsafe.d"], ["-save-temps"], ["other.c"], ["-specs=x"], ["-MG"]):
            with self.subTest(extra=extra):
                result = self.launch(self.argv(extra=extra))
                self.assertNotEqual(result.returncode, 0)
                self.assertFalse((self.root / "main.o").exists())
        self.assertTrue(all(r["status"] == "FAILED" for r in self.receipts()))

    def test_compile_and_probe_failures_preserve_diagnostics(self):
        self.write("main.c", "syntax error here\n")
        result = self.launch()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("error:", result.stderr)
        receipt, = self.receipts()
        self.assertEqual(receipt["compiler_exit_status"], result.returncode)
        self.assertEqual(receipt["status"], "FAILED")
        self.write("main.c", '#include "absent.h"\n')
        result = self.launch(self.argv(output="missing-header.o"))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("absent.h", result.stderr)
        self.assertTrue(any(r["compiler_exit_status"] is None for r in self.receipts()))
        self.assertEqual(len(list(self.evidence.glob("*/in-progress.json"))), 2)

    def test_preexisting_output_is_preserved_and_never_success(self):
        target = self.write("main.o", "unrelated previous bytes")
        result = self.launch()
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(target.read_text(), "unrelated previous bytes")
        self.assertEqual(self.receipts()[0]["status"], "FAILED")

    def test_repeated_and_concurrent_invocations_have_distinct_receipts(self):
        with ThreadPoolExecutor(max_workers=4) as pool:
            results = list(pool.map(lambda n: self.launch(self.argv(output=f"{n}.o")), range(4)))
        self.assertTrue(all(r.returncode == 0 for r in results), [r.stderr for r in results])
        self.assertEqual(len(self.receipts()), 4)
        self.assertEqual(len({r["invocation_id"] for r in self.receipts()}), 4)

    def test_invalid_timeouts_and_caller_success_claim_are_rejected(self):
        for value in ("0", "-1", "nan", "inf", "garbage"):
            self.assertNotEqual(self.launch(timeout=value).returncode, 0)
        self.assertNotEqual(self.launch(extra=["--passed"]).returncode, 0)
        self.assertFalse((self.root / "main.o").exists())

    def run_with_compile_hook(self, hook, after=False, profile=None, args=None):
        original = self.api._run
        def controlled(argv, **kwargs):
            if kwargs["phase"] == "compile" and not after:
                hook()
            result = original(argv, **kwargs)
            if kwargs["phase"] == "compile" and after:
                hook()
            return result
        previous = Path.cwd()
        try:
            os.chdir(self.root)
            with patch.object(self.api, "_run", side_effect=controlled), redirect_stderr(io.StringIO()), redirect_stdout(io.StringIO()):
                options = {} if profile is None else {'profile': profile}
                return self.api.observe(args or self.argv(), evidence_dir=self.evidence, timeout=10, **options)
        finally:
            os.chdir(previous)

    def test_mutated_dependency_is_failed_despite_real_successful_compile(self):
        result, path = self.run_with_compile_hook(lambda: self.write("generated.h", "#define VALUE 43\n"))
        self.assertNotEqual(result, 0)
        receipt = json.loads(path.read_text())
        self.assertEqual(receipt["compiler_exit_status"], 0)
        self.assertIn("changed", receipt["error"])
        self.assertTrue((self.root / "main.o").stat().st_size)

    def test_new_closure_dependency_is_detected_even_if_original_files_are_identical(self):
        self.write("main.c", '#if __has_include("appeared.h")\n#include "appeared.h"\n#endif\nint main(void) { return 0; }\n')
        result, path = self.run_with_compile_hook(lambda: self.write("appeared.h", "#define NEW 1\n"))
        receipt = json.loads(path.read_text())
        self.assertNotEqual(result, 0)
        self.assertEqual(receipt["compiler_exit_status"], 0)
        self.assertIn("inputs_rehashed_after", receipt)
        self.assertEqual(receipt["inputs_before"], receipt["inputs_rehashed_after"])
        self.assertIn("appeared.h", [r["logical_path"] for r in receipt["inputs_after"]])

    def test_removed_dependency_is_rehashed_and_reported_when_post_probe_fails(self):
        result, path = self.run_with_compile_hook(lambda: (self.root / "generated.h").unlink(), after=True)
        receipt = json.loads(path.read_text())
        self.assertNotEqual(result, 0)
        self.assertEqual(receipt["compiler_exit_status"], 0)
        self.assertIn("inputs_rehashed_after", receipt)
        self.assertTrue(any(r["logical_path"] == "generated.h" and "error" in r for r in receipt["inputs_rehashed_after"]))

    def test_inline_asm_hidden_file_is_explicitly_outside_observed_closure(self):
        self.write("payload.bin", "abc")
        self.write("main.c", '__asm__(".pushsection .rodata\\n.incbin \\\"payload.bin\\\"\\n.popsection");\nint main(void) { return 0; }\n')
        result = self.launch()
        self.assertEqual(result.returncode, 0, result.stderr)
        receipt, = self.receipts()
        self.assertEqual(receipt.get("assembler_input_coverage"), "NOT_VERIFIED")
        self.assertEqual(receipt["build_source_attestation"], "NOT_VERIFIED")
        self.assertNotIn("payload.bin", [r["logical_path"] for r in receipt["inputs_before"]])

    def test_object_is_exclusively_reserved_before_compiler_launch(self):
        def verify_reservation():
            with self.assertRaises(FileExistsError):
                with (self.root / "main.o").open("xb"):
                    pass
        result, path = self.run_with_compile_hook(verify_reservation)
        self.assertEqual(result, 0, json.loads(path.read_text()))

    def test_symlink_retarget_with_identical_bytes_is_detected(self):
        self.write("version-a.h", "#define VALUE 42\n")
        self.write("version-b.h", "#define VALUE 42\n")
        (self.root / "generated.h").unlink()
        (self.root / "generated.h").symlink_to("version-a.h")
        def retarget():
            (self.root / "generated.h").unlink()
            (self.root / "generated.h").symlink_to("version-b.h")
        result, path = self.run_with_compile_hook(retarget)
        self.assertNotEqual(result, 0)
        receipt = json.loads(path.read_text())
        before = next(r for r in receipt["inputs_before"] if r["logical_path"] == "generated.h")
        after = next(r for r in receipt["inputs_after"] if r["logical_path"] == "generated.h")
        self.assertEqual(before["sha256"], after["sha256"])
        self.assertNotEqual(before["resolved_path"], after["resolved_path"])

    def test_symlink_retarget_to_same_final_file_is_detected(self):
        self.write("version.h", "#define VALUE 42\n")
        (self.root / "alias.h").symlink_to("version.h")
        (self.root / "generated.h").unlink()
        (self.root / "generated.h").symlink_to("version.h")
        def retarget():
            (self.root / "generated.h").unlink()
            (self.root / "generated.h").symlink_to("alias.h")
        result, path = self.run_with_compile_hook(retarget)
        self.assertNotEqual(result, 0)
        self.assertIn("changed", json.loads(path.read_text())["error"])

    def test_missing_and_empty_outputs_never_succeed(self):
        for remove in (True, False):
            with self.subTest(remove=remove):
                def invalidate():
                    if remove:
                        (self.root / "main.o").unlink()
                    else:
                        (self.root / "main.o").write_bytes(b"")
                result, path = self.run_with_compile_hook(invalidate, after=True)
                self.assertNotEqual(result, 0)
                self.assertEqual(json.loads(path.read_text())["compiler_exit_status"], 0)

    def test_implicit_precompiled_header_is_rejected(self):
        subprocess.run([self.cc, "-x", "c-header", "outer.h", "-o", "outer.h.gch"], cwd=self.root, check=True, capture_output=True)
        result = self.launch()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("precompiled", result.stderr)
        self.assertFalse((self.root / "main.o").exists())

    def test_outputs_cannot_overwrite_inputs_or_receipts(self):
        for depfile in ("generated.h", "main.o", str(self.evidence / "shared.d")):
            with self.subTest(depfile=depfile):
                result = self.launch(self.argv(extra=["-MD", "-MF", depfile]))
                self.assertNotEqual(result.returncode, 0)
                self.assertEqual((self.root / "generated.h").read_text(), "#define VALUE 42\n")
                self.assertFalse((self.root / "main.o").exists())

    def test_compiler_symlink_retarget_is_detected(self):
        compiler = self.root / "compiler"
        compiler.symlink_to(self.cc)
        self.cc = str(compiler)
        alternate = shutil.which("c++")
        self.assertIsNotNone(alternate)
        def retarget():
            compiler.unlink()
            compiler.symlink_to(alternate)
        result, path = self.run_with_compile_hook(retarget, after=True)
        self.assertNotEqual(result, 0)
        receipt = json.loads(path.read_text())
        self.assertIn("compiler identity", receipt["error"])

    def test_compiler_path_search_shadowing_is_detected(self):
        (self.root / "first").mkdir()
        (self.root / "second").mkdir()
        (self.root / "second/compiler").symlink_to(self.cc)
        self.cc = "compiler"
        with patch.dict(os.environ, {"PATH": str(self.root / "first") + ":" + str(self.root / "second") + ":" + os.environ["PATH"]}):
            result, path = self.run_with_compile_hook(lambda: (self.root / "first/compiler").symlink_to(shutil.which("c++")))
        self.assertNotEqual(result, 0)
        self.assertIn("compiler", json.loads(path.read_text())["error"])

    def test_pch_in_earlier_search_directory_is_rejected(self):
        self.write("early/outer.h", '#include "../generated.h"\n')
        subprocess.run([self.cc, "-x", "c-header", "early/outer.h", "-o", "early/outer.h.gch"], cwd=self.root, check=True, capture_output=True)
        (self.root / "early/outer.h").unlink()
        self.write("main.c", '#include <outer.h>\nint main(void) { return VALUE != 42; }\n')
        result = self.launch(self.argv(extra=["-Iearly", "-I."]))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("precompiled", result.stderr)

    def test_timeout_retains_receipt_and_terminates_child_group(self):
        # A controlled subprocess hook exercises real process-group cleanup;
        # this timeout fixture is NOT represented as an actual compilation.
        original = self.api._run
        code = 'import subprocess,sys,time; p=subprocess.Popen([sys.executable,"-c","import time; time.sleep(30)"]); open("child.pid","w").write(str(p.pid)); time.sleep(30)'
        def timeout_child(argv, **kwargs):
            if kwargs["phase"] == "compile":
                kwargs["deadline"] = time.monotonic() + 0.3
                argv = [sys.executable, "-c", code]
            return original(argv, **kwargs)
        previous = Path.cwd()
        try:
            os.chdir(self.root)
            with patch.object(self.api, "_run", side_effect=timeout_child), redirect_stderr(io.StringIO()), redirect_stdout(io.StringIO()):
                result, path = self.api.observe(self.argv(), evidence_dir=self.evidence, timeout=10)
        finally:
            os.chdir(previous)
        self.assertEqual(result, 124)
        receipt = json.loads(path.read_text())
        self.assertEqual(receipt["status"], "FAILED")
        self.assertTrue(any(op.get("incomplete") for op in receipt["operations"]))
        pid = int((self.root / "child.pid").read_text())
        state = Path(f"/proc/{pid}/stat")
        self.assertTrue(not state.exists() or state.read_text().split()[2] == "Z")

    def test_environment_include_selection_is_recorded_without_unrelated_secrets(self):
        self.write("env/selected.h", "#define VALUE 42\n")
        self.write("main.c", '#include <selected.h>\nint main(void) { return VALUE != 42; }\n')
        env = dict(os.environ, CPATH=str(self.root / "env"), OBSERVER_TEST_SECRET="never-record-this-token")
        result = self.launch(env=env)
        self.assertEqual(result.returncode, 0, result.stderr)
        receipt, = self.receipts()
        self.assertEqual(receipt["environment"]["CPATH"], str(self.root / "env"))
        self.assertNotIn("never-record-this-token", json.dumps(receipt))
        self.assertTrue(any(r["logical_path"].endswith("env/selected.h") for r in receipt["inputs_before"]))

    def test_compiler_wrappers_cannot_masquerade_as_observed_binary(self):
        wrapper = self.write("wrapper", '#!/bin/sh\nexec "' + self.cc + '" "$@"\n')
        wrapper.chmod(0o700)
        result = self.launch([str(wrapper), "-c", "main.c", "-o", "main.o"])
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("ELF", result.stderr)

    def test_real_cxx_current_build_flags_and_forced_headers(self):
        self.write("nested/selected.h", "#define SELECTED 42\n")
        self.write("forced.h", "#define FORCED 42\n")
        self.write("main.cpp", '#include <selected.h>\nint main() { return SELECTED != FORCED; }\n')
        args = [shutil.which("c++"), "-isystem", "nested", "-include", "forced.h", "-fPIC", "-O2", "-g", "-DNDEBUG", "-std=c++17", "-frounding-math", "-Wno-invalid-offsetof", "-pthread", "-MMD", "-MFmain.d", "-MQmain.o", "-c", "main.cpp", "-omain.o"]
        result = self.launch(args)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("forced.h", [r["logical_path"] for r in self.receipts()[0]["inputs_before"]])

    def test_real_cxx_popcnt_flag_produces_observed_object_and_runnable_executable(self):
        self.write("main.cpp", "int main() { return __builtin_popcount(0xf0f0u) != 8; }\n")
        self.cc = shutil.which("c++")
        args = self.argv(extra=["-mpopcnt", "-MMD", "-MF", "main.d"])
        result = self.launch(args)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue((self.root / "main.o").is_file())
        self.assertTrue((self.root / "main.d").is_file())
        receipt, = self.receipts()
        self.assertEqual(receipt["status"], "OBSERVED_SUCCESS")
        self.assertEqual(receipt["object"]["sha256"], hashlib.sha256((self.root / "main.o").read_bytes()).hexdigest())
        subprocess.run([self.cc, "main.o", "-o", "popcnt"], cwd=self.root, check=True, capture_output=True)
        self.assertEqual(subprocess.run([str(self.root / "popcnt")]).returncode, 0)

    def test_existing_depfile_is_untouched_by_both_probes(self):
        self.write("original.d", "untouched before compilation\n")
        original = self.api._run
        def checked(argv, **kwargs):
            if kwargs["phase"] == "compile":
                self.assertEqual((self.root / "original.d").read_text(), "untouched before compilation\n")
            result = original(argv, **kwargs)
            if kwargs["phase"] == "compile":
                self.compiled_depfile = (self.root / "original.d").read_bytes()
            elif kwargs["phase"] == "dependencies-after":
                self.assertEqual((self.root / "original.d").read_bytes(), self.compiled_depfile)
            return result
        previous = Path.cwd()
        try:
            os.chdir(self.root)
            with patch.object(self.api, "_run", side_effect=checked):
                result, _ = self.api.observe(self.argv(extra=["-MD", "-MF", "original.d"]), evidence_dir=self.evidence, timeout=10)
        finally:
            os.chdir(previous)
        self.assertEqual(result, 0)

    def test_real_probe_timeout_on_fifo_header_retains_failure(self):
        os.mkfifo(self.root / "fifo.h")
        self.write("main.c", '#include "fifo.h"\nint main(void) { return 0; }\n')
        result = self.launch(timeout="0.2")
        self.assertEqual(result.returncode, 124, result.stderr)
        receipt, = self.receipts()
        self.assertEqual(receipt["status"], "FAILED")
        self.assertIsNone(receipt["compiler_exit_status"])

    def test_real_probe_interruption_retains_failure_and_stops_group(self):
        os.mkfifo(self.root / "fifo.h")
        self.write("main.c", '#include "fifo.h"\n')
        child = subprocess.Popen([sys.executable, str(SCRIPT), "--evidence-dir", str(self.evidence), "--timeout", "10", "--", *self.argv()], cwd=self.root, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        try:
            deadline = time.monotonic() + 5
            while not list(self.evidence.glob("*/dependencies-before.stdout")) and child.poll() is None and time.monotonic() < deadline:
                time.sleep(0.01)
            self.assertTrue(list(self.evidence.glob("*/dependencies-before.stdout")))
            child.terminate()
            _, err = child.communicate(timeout=5)
            self.assertEqual(child.returncode, 130, err)
            self.assertEqual(self.receipts()[0]["status"], "FAILED")
        finally:
            if child.poll() is None:
                child.kill()
                child.communicate()

    def clang_argv(self, output='clang.o', extra=()):
        self.assertTrue(Path('/usr/bin/clang').is_file(), 'real Clang C required; do not skip')
        return ['/usr/bin/clang', '--no-default-config', '-x', 'c', *extra,
                '-c', 'main.c', '-o', output]

    def test_clang_profile_observes_actual_c_headers_and_runnable_object(self):
        args = self.clang_argv(extra=['-O2', '-fPIC', '-MD', '-MF', 'clang.d'])
        result = self.launch(args, extra=['--profile', 'clang-c-v1'])
        self.assertEqual(result.returncode, 0, result.stderr)
        receipt, = self.receipts()
        self.assertEqual(receipt['settings']['profile'], 'clang-c-v1')
        self.assertEqual(receipt['argv'], args)
        self.assertEqual(receipt['compiler_before']['family'], 'clang')
        self.assertEqual(receipt['inputs_before'], receipt['inputs_after'])
        paths = [r['logical_path'] for r in receipt['inputs_before']]
        for name in ('main.c', 'outer.h', 'generated.h'):
            self.assertIn(name, paths)
        self.assertTrue(any(p.endswith('stddef.h') for p in paths))
        self.assertEqual(receipt['object']['sha256'], hashlib.sha256((self.root/'clang.o').read_bytes()).hexdigest())
        self.assertEqual(receipt['build_source_attestation'], 'NOT_VERIFIED')
        self.assertTrue((self.root/'clang.d').read_text().startswith('clang.o:'))
        subprocess.run(['/usr/bin/clang', '--no-default-config', 'clang.o', '-o', 'clang-tiny'], cwd=self.root, check=True, capture_output=True)
        self.assertEqual(subprocess.run([str(self.root/'clang-tiny')]).returncode, 0)

    def test_clang_mutated_header_and_new_closure_cannot_be_accepted(self):
        for new in (False, True):
            self.write('main.c', '#include "outer.h"\n#if __has_include("new.h")\n#include "new.h"\n#endif\nint answer(void){return VALUE;}\n')
            def mutate():
                self.write('new.h' if new else 'generated.h', '#define NEW 1\n' if new else '#define VALUE 43\n')
            result, path = self.run_with_compile_hook(mutate, profile='clang-c-v1', args=self.clang_argv(output=f'mutate-{new}.o'))
            receipt = json.loads(path.read_text())
            self.assertNotEqual(result, 0)
            self.assertEqual(receipt['compiler_exit_status'], 0)
            self.assertIn('changed', receipt['error'])

    def test_clang_profile_requires_explicit_c_and_no_default_config(self):
        base = self.clang_argv()
        for args in ([a for a in base if a != '--no-default-config'],
                     ['/usr/bin/clang', '--no-default-config', '-c', 'main.c', '-o', 'clang.o'],
                     [a if a != 'c' else 'c++' for a in base]):
            result = self.launch(args, extra=['--profile', 'clang-c-v1'])
            self.assertNotEqual(result.returncode, 0)
            self.assertFalse((self.root/'clang.o').exists())
        self.assertTrue(all(r['compiler_exit_status'] is None and r['settings']['profile'] == 'clang-c-v1'
                            for r in self.receipts()))
        self.assertEqual(len(self.receipts()), 3)

    def test_clang_profile_rejects_unobserved_input_mechanisms(self):
        for flags in (['-include', 'outer.h'], ['-imacros', 'outer.h'], ['-include-pch', 'outer.pch'],
                      ['-fmodules'], ['-fmodule-file=outer.pcm'], ['--config=extra.cfg'],
                      ['-Xclang', '-include-pch'], ['@args']):
            with self.subTest(flags=flags):
                result = self.launch(self.clang_argv(extra=flags), extra=['--profile', 'clang-c-v1'])
                self.assertNotEqual(result.returncode, 0)
                self.assertFalse((self.root/'clang.o').exists())
        self.assertEqual(len(self.receipts()), 8)
        self.assertTrue(all(r['compiler_exit_status'] is None for r in self.receipts()))

    def test_clang_hidden_argument_edits_are_rejected(self):
        result = self.launch(self.clang_argv(), extra=['--profile', 'clang-c-v1'],
                             env=dict(os.environ, CCC_OVERRIDE_OPTIONS='+-DVALUE=7'))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('CCC_OVERRIDE_OPTIONS', result.stderr)
        self.assertFalse((self.root/'clang.o').exists())

    def test_gnu_default_does_not_silently_admit_clang(self):
        result = self.launch(['/usr/bin/clang', '-c', 'main.c', '-o', 'clang.o'])
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('GNU compiler required', result.stderr)
        self.assertFalse((self.root/'clang.o').exists())

    def test_clang_compile_failure_retains_real_diagnostics(self):
        self.write('main.c', 'syntax error here\n')
        result = self.launch(self.clang_argv(), extra=['--profile', 'clang-c-v1'])
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('error:', result.stderr)
        receipt, = self.receipts()
        self.assertEqual(receipt['compiler_exit_status'], result.returncode)
        self.assertEqual(receipt['status'], 'FAILED')

    def test_clang_late_language_cannot_masquerade_as_observed_c(self):
        self.assertTrue(Path('/usr/bin/clang++').is_file(), 'real Clang++ required')
        self.write('main.c', '#ifndef __cplusplus\n#error requires C++\n#endif\nint answer(){return 42;}\n')
        args = ['/usr/bin/clang++', '--no-default-config', '-c', 'main.c', '-x', 'c', '-o', 'late.o']
        result = self.launch(args, extra=['--profile', 'clang-c-v1'])
        self.assertNotEqual(result.returncode, 0, 'C++ source was incorrectly accepted under C-only observation')
        self.assertFalse((self.root/'late.o').exists())
        receipt, = self.receipts()
        self.assertIsNone(receipt['compiler_exit_status'])

    def test_clang_explicit_early_c_overrides_cxx_driver(self):
        self.write('main.c', '#ifdef __cplusplus\n#error requires C\n#endif\nint answer(void){return 42;}\n')
        args = ['/usr/bin/clang++', '--no-default-config', '-x', 'c', '-c', 'main.c', '-o', 'early.o']
        result = self.launch(args, extra=['--profile', 'clang-c-v1'])
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.receipts()[0]['status'], 'OBSERVED_SUCCESS')


class CompilerObserverCMakeTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.assertIsNotNone(shutil.which("cmake"), "cmake required; do not skip integration")
        (self.root / "main.c").write_text("int answer(void) { return 42; }\n")
        (self.root / "main.cpp").write_text('extern "C" int answer(void); int main() { return answer()!=42; }\n')
        (self.root / "CMakeLists.txt").write_text('cmake_minimum_required(VERSION 3.16)\nproject(observer_fixture LANGUAGES C CXX)\ninclude("' + str(ROOT / "cmake/CompilerObservation.cmake") + '")\nadd_executable(tiny main.c main.cpp)\n')

    def configure(self, *flags):
        return subprocess.run(["cmake", "-S", str(self.root), "-B", str(self.root / "build"), *flags], text=True, capture_output=True, timeout=30)

    def test_real_cmake_launcher_observes_c_and_cxx_compilation(self):
        evidence = self.root / "evidence"
        result = self.configure("-DMAKE_IT_STAND_OBSERVE_COMPILERS=ON", "-DMAKE_IT_STAND_COMPILER_EVIDENCE_DIR=" + str(evidence), "-DCMAKE_BUILD_TYPE=RelWithDebInfo")
        self.assertEqual(result.returncode, 0, result.stderr)
        built = subprocess.run(["cmake", "--build", str(self.root / "build"), "-j", "2"], capture_output=True, text=True, timeout=30)
        self.assertEqual(built.returncode, 0, built.stdout + built.stderr)
        self.assertEqual(subprocess.run([str(self.root / "build/tiny")]).returncode, 0)
        receipts = [json.loads(p.read_text()) for p in evidence.glob("*/receipt.json")]
        self.assertEqual(len(receipts), 2)
        self.assertTrue(all(r["status"] == "OBSERVED_SUCCESS" for r in receipts))
        self.assertTrue((self.root / "build/compile_commands.json").is_file())

    def test_default_is_off_and_legacy_build_succeeds_without_evidence(self):
        result = self.configure()
        self.assertEqual(result.returncode, 0, result.stderr)
        result = subprocess.run(["cmake", "--build", str(self.root / "build"), "-j", "2"], capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertFalse((self.root / "evidence").exists())

    def test_configuration_rejects_missing_root_bad_timeout_and_existing_launcher(self):
        for flags in ([], ["-DMAKE_IT_STAND_COMPILER_EVIDENCE_DIR=" + str(self.root / "evidence"), "-DMAKE_IT_STAND_COMPILER_TIMEOUT=nan"],
                      ["-DMAKE_IT_STAND_COMPILER_TIMEOUT=10", "-DCMAKE_CXX_COMPILER_LAUNCHER=ccache"]):
            with self.subTest(flags=flags):
                result = self.configure("-DMAKE_IT_STAND_OBSERVE_COMPILERS=ON", *flags)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("Compiler observation", result.stderr)

    def test_configuration_rejects_false_like_existing_launcher(self):
        result = self.configure("-DMAKE_IT_STAND_OBSERVE_COMPILERS=ON", "-DMAKE_IT_STAND_COMPILER_EVIDENCE_DIR=" + str(self.root / "evidence"), "-DCMAKE_CXX_COMPILER_LAUNCHER=OFF")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("existing", result.stderr)

    def test_configuration_rejects_non_linux_target(self):
        result = self.configure("-DMAKE_IT_STAND_OBSERVE_COMPILERS=ON", "-DMAKE_IT_STAND_COMPILER_EVIDENCE_DIR=" + str(self.root / "evidence"), "-DCMAKE_SYSTEM_NAME=Generic")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("GNU/Linux", result.stderr)

    def test_configuration_rejects_nonfinite_overflowing_timeout(self):
        result = self.configure("-DMAKE_IT_STAND_OBSERVE_COMPILERS=ON", "-DMAKE_IT_STAND_COMPILER_EVIDENCE_DIR=" + str(self.root / "evidence"), "-DMAKE_IT_STAND_COMPILER_TIMEOUT=" + "9" * 320)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("timeout", result.stderr.lower())


if __name__ == "__main__":
    unittest.main()
