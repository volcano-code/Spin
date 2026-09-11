# Task5n: retained isolated Spin2014 dependency-test workspaces

Approved continuation of Task8 test correctness and Task5m followup. Work only /data/xvxinyao/make-it-balance-paper-sim, existing spin-it-strict-computational-v2 linked worktree. BASE99cfbf0. User approved starting the next work. This is a bounded safety/verification fix, not paper-algorithm completion.

## Root cause and required change

Existing tests/cmake/assert_spin_it_2014_dependency.cmake unconditionally REMOVE_RECURSEs fixed SOURCE_DIR/build/spin-2014-{disabled,enabled}-contract. Prior run interrupted after30 passing CTests; preserve it and both old directories. Do NOT run the old script against the real checkout for RED. No deletion or cleanup, including temporary historical results, in this task.

Replace only this2014 dependency script's fixed scratch with unique per-invocation retained directories under the invoking build directory. CMake test registration must explicitly pass the current binary root; direct invocation missing/invalid root fails before writes. Avoid new Python requirement for headless/native2014. A small CMake helper may allocate a unique run directory with CMake random token and short-lived allocation file lock in a function (release before expensive build); reject collisions, never reuse/overwrite run directories. Validate root as absolute existing directory, not source root, and ensure created run resolves under it (reject scratch symlink escape). A persistent allocator lock is allowed, existing experiment data is not modified. Do not build a generalized test framework.

Each actual invocation preserves disabled/enabled build trees and configure/build/runtime stdout/stderr and exit metadata, including failures. Print the retained directory before launching work so CTest interruption still leaves discoverable evidence. Keep all current scientific/dependency assertions:2014OFF excludes implementation;2014ON/FasterOFF/reachabilityOFF builds actual audited target and runs it; no hidden CGAL dependency. No --passed flags, skip shortcuts or assertion removal. Do not change algorithm tolerances or claim old interrupted run became successful. Existing timeout240 may be increased only if actual clean compilation shows resource time exhaustion, with controller ruling; do not preemptively modify it.

Files: existing script; CMakeLists.txt scoped registration; at most one small reusable allocator helper and one focused behavior test. Documentation docs/spin2014_dependency_test_isolation.md. Do not touch other dependency scripts/Faster algorithms or implement source-adaptive work yet.

## TDD and validation

Behavioral regression must exercise allocation/script behavior, not grep source. Safely reproduce old fixed-directory deletion in a newly generated minimal fake repository layout/copy of script, or an equally safe controlled fixture; never target real historical cache dirs. Real CMake configure on a tiny fixture can supply required build graph strings and a tiny audited executable; distinguish harness fixtures from the mandatory actual project build. No fake success counted as actual dependency validation.

Verify: existing sentinel files remain unchanged; two actual helper/script invocations get distinct dirs and do not alter first output; missing/invalid/symlink-escaping destination rejected before mutation; concurrent allocations do not share a directory; explicit failing child leaves logs/nonzero status and no pass. Tests need actual behavior and TDD RED/GREEN logs. Avoid arbitrary large test-only injection API.

Then run the actual registered dependency contract in build/task-5c-no-faster-20260909. Permit its new nested build dirs ONLY within that build root. Capture actual configure/build/linked audited result and logs. Run scoped 2014 group once including source-spin-material-python and new isolation test; keep all logs and count actual tests, no NotRun/skips. You own builds until report, -j3 max, ulimit-v25165824, timeout3600, sequential heavy work. Do not invoke unrelated destructive contracts. Preserve preexisting dependency warnings accurately.

## Protected work and reporting

Unrelated dirty files: include/spin/QuadricMomentEvaluator.h; src/spin/QuadricHollowingSolver.cpp; src/spin/QuadricMomentEvaluator.cpp; tests/quadric_moment_evaluator_regression.cpp; untracked src/spin/QuadricMaterialExportGuard.h and tests/quadric_certified_empty_cases.h. Do not stage/edit these, or protected /data/xvxinyao/make-it-balance. Do not rebuild build/strict-v2-spin. No merge/push/tag/release/deletion. Use apply_patch. User identity volcano <xuxinyao307@gmail.com> authorized for scoped commits only.

No subagents: implement and self-review yourself; controller performs independent checks and fresh review. Commit explicit owned files. Full report with commands/RED/GREEN/current evidence/limitations/hashes at .superpowers/sdd/2026-09-05-spin-it-strict-computational-v2/task-5n-report.md; return under15 lines status,commit,tests,concerns,path. OverallPARTIAL_BLOCKED. Ask controller for unclear contracts; don't silently weaken requirements.
