# Task5n report: retained isolated Spin2014 dependency-test workspaces

## Status and scope

**Overall status: PARTIAL_BLOCKED.** The bounded dependency-test isolation fix is
implemented and verified. This does not complete the paper algorithm,
source-adaptive binding, release evidence, or the previously interrupted Task5m
suite. No numerical solver, tolerance, Faster algorithm, paper asset, protected
main checkout, merge, push, tag, release, deletion, or cleanup was performed.

Worktree: `/data/xvxinyao/make-it-balance-paper-sim`, branch
`spin-it-strict-computational-v2`, base and starting HEAD
`99cfbf07b58967a8243108fee2f7590b2a021bd8`. The worktree is linked
(`GIT_DIR` differs from `GIT_COMMON`).

## Implemented contract

- `spin-2014-dependency-contract` now receives the invoking
  `CMAKE_CURRENT_BINARY_DIR` explicitly. Missing, relative, nonexistent, and
  source-root destinations fail before allocation.
- `make_it_stand_allocate_retained_test_directory(root, out)` validates the
  existing absolute root, rejects a symlink at the reserved run namespace,
  takes a short function-scoped file lock, rejects collisions, and creates a
  random never-reused run beneath
  `spin-2014-dependency-contract-runs/`. The lock is released on function return
  before configure/build/runtime work.
- The contract prints `SPIN_IT_2014_DEPENDENCY_RETAINED_DIR=<path>` before child
  work. Each configure/build/audited phase durably records `.command`, direct
  `.stdout.log`, direct `.stderr.log`, and `.status`; status is `RUNNING` before
  launch and then the actual exit result. `contract.status` records the durable
  overall state. An interrupted CTest therefore leaves the run location, live
  partial output, and active phase discoverable.
- Disabled and enabled builds are retained under the unique run. The original
  assertions remain: 2014-OFF excludes implementation; 2014-ON with Faster-OFF
  and reachability-OFF includes `AdaptiveVoxelSpinSolver.cpp`, excludes hidden
  CGAL/GMP/MPFR graph dependencies, builds the real `spin-it-2014-audited`
  target at `-j3`, and runs it.
- The focused CMake behavior test creates only unique tiny fake repositories
  beneath the owned build root. It checks historical sentinels, two distinct
  script runs and an unchanged first run, rejection-before-mutation, namespace
  symlink escape, concurrent allocation, numeric child failure evidence, and
  live stdout/stderr plus `RUNNING` status before a sleeping child exits. These
  fixtures are harness checks and are not counted as actual project dependency
  validation.

Owned source files:

- `CMakeLists.txt`
- `tests/cmake/assert_spin_it_2014_dependency.cmake`
- `tests/cmake/allocate_retained_test_directory.cmake`
- `tests/cmake/assert_spin_it_2014_dependency_isolation.cmake`
- `docs/spin2014_dependency_test_isolation.md`
- this report

## Behavioral TDD

All commands ran with `ulimit -v 25165824`; generated data was retained under
`build/task-5c-no-faster-20260909/`.

1. Initial RED, against an unchanged copy of the old real script in a generated
   fake repository:

   ```sh
   cmake -DTEST_BINARY_ROOT=$PWD/build/task-5c-no-faster-20260909 \
     -DDEPENDENCY_SCRIPT=$PWD/tests/cmake/assert_spin_it_2014_dependency.cmake \
     -P tests/cmake/assert_spin_it_2014_dependency_isolation.cmake
   ```

   Exit 1. The old script deleted the controlled
   `success-source/build/spin-2014-disabled-contract/sentinel.txt`, exactly
   reproducing the root cause without touching the real historical directories.
   Raw log: `task-5n-evidence/tdd-red-sentinel.log`, SHA-256
   `4c8bf391e1f3efb07a5789471c91880cfb32b22326a0d07aef48e61abfe76d80`.

2. Initial GREEN after unique allocation and retained result files: exit 0.
   Raw log `task-5n-evidence/tdd-green-01.log`, SHA-256
   `e4ad69ec48fa225eece8303861964ece4c2ee56ff51adf0e094ab56e59be1cc0`.

3. Live-output RED after adding the controlled sleeping-child observation:
   exit 1 with `child logs and RUNNING status were not visible before child
   exit`. Raw log `task-5n-evidence/live-log-red.log`, SHA-256
   `13d7c266fa0a5564b7b64ba9fb02e00290171a4ef40ae4453d48e45e37b94a59`.

4. Live-output GREEN after direct `OUTPUT_FILE`/`ERROR_FILE` capture and
   prelaunch phase status: exit 0. Raw log
   `task-5n-evidence/live-log-green.log`, SHA-256
   `67ee3c6e5e5b87a61f790139a411b193505bdc4e85ec40ff776e2034159d0791`.

5. Exact-final registered behavior test:

   ```sh
   ctest --test-dir build/task-5c-no-faster-20260909 \
     -R '^spin-2014-dependency-isolation$' --output-on-failure
   ```

   1/1 passed in 6.32 s. Raw log
   `task-5n-evidence/final2-registered-isolation.log`, SHA-256
   `4576496272ca586b7c93d59fed90f381cadcd72815ad9a8d13c9f4e93f0b2302`.

## Actual project dependency validation

The exact-final production contract was run through its real CTest registration,
not a fixture:

```sh
ctest --test-dir build/task-5c-no-faster-20260909 \
  -R '^spin-2014-dependency-contract$' -V
```

Result: 1/1 passed in 89.25 s with CTest timeout still 240 s. Retained run:
`build/task-5c-no-faster-20260909/spin-2014-dependency-contract-runs/run-c94c38ea111af14c/`.
All four child statuses are 0 and `contract.status` is `PASSED`. Audited output
reported `status=CONVERGED x=0.75 residual=0`, rejected the subnormal-curvature
update while retaining inverse `0.125`, exercised restoration/invalid-trial and
LP-cap outcomes, reported Givens KKT `6.3137399822221774e-11`, and printed the
actual audited PASS line. CTest log:
`task-5n-evidence/final2-actual-dependency-contract.log`, SHA-256
`e555afe2b0389b97302c3721a1a16ddfbcb3d05c95fed6fdd93d23ed486410b9`.

Key exact-final retained hashes:

- `audited.stdout.log`:
  `01d25b5c038ccda81ff3d24055fce47769ba5a480d4d6e5d98946c1f5499994b`
- `build-enabled.stdout.log`:
  `217810504e46f38f08f0c20440d7b2e9b3d16fb94640e9d31b59a386a29ec833`
- `configure-disabled.stdout.log`:
  `ba3a1605c3091aa8d9d2fd7a722b3102cf078ccbcb2c9b0f03e5e1dac9322280`
- `configure-enabled.stdout.log`:
  `fd3f7b4b8bbce574174445c5681d7dc6892ee9fe06504b5baa4f966365392493`
- both configure stderr logs:
  `036e239575841fdf40ef29449f6c438983422df4989d27df0b975639e51d2f5e`
- empty build/runtime stderr:
  `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`

The configure stderr retains the two pre-existing CMake minimum-version
deprecation warnings. Dependencies and warnings were not suppressed or changed.

## Scoped 2014 tests

After registration and before the later live-output retention correction, the
requested sequential label run was executed once:

```sh
timeout 3600 ctest --test-dir build/task-5c-no-faster-20260909 \
  -L '^spin-it-2014$' -j1 --output-on-failure
```

Exactly 29/29 unique tests passed in 217.40 s, with zero failure, NotRun, or skip
records. It includes `source-spin-material-python` (23.35 s), the dependency
contract (90.86 s), and the registered isolation test (1.96 s). This retained
label result is accurately a pre-live-output-fix snapshot; it is not represented
as final production-contract evidence. Raw log
`task-5n-evidence/scoped-spin-2014.log`, SHA-256
`853558d6e29d298592eae98278047ede1a1ad3eb17efafaf3416e24258833c85`.

The label intentionally omits three existing 2014 tests. They were run once on
the final source snapshot:

```sh
ctest --test-dir build/task-5c-no-faster-20260909 \
  -R '^(spin-it-2014-surface-value|spin-it-2014-surface-value-allocation|spin-it-2014-audited)$' \
  -j1 --output-on-failure
```

3/3 passed in 0.02 s. Raw log
`task-5n-evidence/final-unlabelled-2014.log`, SHA-256
`afc55063181605d987e382f4cb29123879475d7c13b4a74615a208c34894af87`.
No broad all-tests suite was run, and no unbuilt unrelated target is claimed.

## Preservation, hashes, and self-review

The old real `build/spin-2014-disabled-contract` and
`build/spin-2014-enabled-contract` directories remain present. The destructive
old script was never invoked against the real checkout. This task does not claim
to recover cache contents already removed/recreated during Task5m, and does not
reinterpret that interrupted run as passing.

The six externally dirty Faster paths retained their preflight SHA-256 values
through final verification:

- `include/spin/QuadricMomentEvaluator.h`:
  `85776f2175a2f9de528cf3aa28eb3baf8cef08bfe7c626b2f2523a491d733bd7`
- `src/spin/QuadricHollowingSolver.cpp`:
  `47d3f0f371d5766f8f1a2c6d0e1dc74ba39f94d5c90bdea018ec02aec5234490`
- `src/spin/QuadricMomentEvaluator.cpp`:
  `44f72591b0e5eab1be36a277aa67148e40a750a207b2ded163d7107754212ec7`
- `tests/quadric_moment_evaluator_regression.cpp`:
  `38130fc195cb3d1555b96c518d68b967d932a77bd1de691de5ad4946642dd640`
- `src/spin/QuadricMaterialExportGuard.h`:
  `4accc55b5157df0a79d385e1ab3062ddc3a7330750224c57c118f545a2927cf4`
- `tests/quadric_certified_empty_cases.h`:
  `9ac869c24b14f1c07e3b3c7ecd076cef1ef3fef16e28f91a67a4a9bdd856beb7`

Final owned source hashes before commit:

- `CMakeLists.txt`:
  `93ba14bcb83c0140e45c2670f0a9089a9b7e2dcbefbc890cbc97046e7029ee8a`
- dependency script:
  `ad69496196ba5f2c160865678d85c2ef24bcd0c9eb9c304e7985d19b4e1a0baf`
- allocator helper:
  `13845b7aed6cab20f4294fadd7a9c3b6ae12d2a754ebec21c1cabc6d4d204f93`
- behavior test:
  `9ee9ace4c13bcc9657b32e43468099ec4dd07a44a5241ffa8e1df5b50db8c95f`
- documentation:
  `9a9262a0322175d4bebc1cac9840b80b0fd64d67f8b2c5b7d8a850196192a045`

`git diff --check` passed. Self-review confirmed no `REMOVE_RECURSE`, no cleanup
path, no Python dependency, no solver/tolerance/assertion change, no allocator
lock held during child work, live direct logs with durable active status, actual
failure propagation, and explicit CTest binary-root registration. Generated
fixture and real run directories are intentionally retained and ignored. The
commit uses authorized identity `volcano <xuxinyao307@gmail.com>`; its hash is
reported in the handoff because a commit cannot embed its own hash.

## Limitations

This is dependency-test workspace safety and verification only. It is not a
generalized test framework, cleanup policy, proof against hostile concurrent
filesystem replacement after validation, full-suite result, source-adaptive
binding implementation, scientific result, paper-model reproduction, or release
approval. Overall status therefore remains **PARTIAL_BLOCKED**.
