# Task5y implementer report

Scope: actual deformed binary numerical/geometric **pre-seal** validation.
No SealedSpinArtifact, exact final region partition, continuous BBW Jacobian,
full cavity-web thickness, feasible witness, spin/physics or paper grant.

Status: DONE. Scoped commit `5acbfa8` —
`feat(spin2014): validate actual deformed binary material endpoints`, author
and committer volcano <xuxinyao307@gmail.com>. The commit contains exactly the
five production/test/doc files listed below. Post-commit status retains only
the pre-existing foreign Faster dirt/untracked experiments; no scoped file
remains dirty. Build ownership has been released to the controller.

## Implementation

- `validateSourceCageBinaryMaterial` accepts only the private immutable adaptive
  audit, requires complete/last-solved/no-pending and terminal identity agreement.
- Private-constructor shared immutable result and material; copied facts cannot
  publish, failed results have null `validatedMaterial()`, and `candidate()`
  retains actual emitted mesh/state evidence after numerical/geometry failure.
- Only beta is rounded (`>=.5` hollowed), preserving every angle/handle bit.
  Actual endpoint dimensions are used even when adaptive merging changes them.
- Actual `evaluateValues(binary)` deformed outer and conforming points feed
  emission. Summed oriented primitive references cancel internal faces;
  compaction follows ascending referenced point IDs. The final mesh is outer
  plus reversed carved boundary. Empty cavity is explicit; the nonempty outer
  material still receives exact audit.
- Independent signed tetrahedra integrate all ten emitted material moments.
  Fixed input density and root-side scaling enforces the 1e-10 moment gate.
  Independent fixed-frame transformation recomputes Givens constraints and
  central inertia. Raw and normalized constraints are retained. Normalized
  constraints use rho L^4 / rho L^5 and tolerance 1e-6 by default.
- Central inertia eigenvalues are checked against `-128*epsilon*rho*L^5`,
  matching the prior source binary moment PSD allowance. Raw eigenvalues,
  fixed normalized tolerance and PSD result are retained and hashed.
- Fresh actual source geometry checker plus exact final material boundary
  audit AnyNonzero run even on constraint failure. Established numerical
  failures remain BINARY_VALIDATION_FAILED if a later checker cannot complete.
- Checked pre-allocation count/index/byte/diagnostic bounds and full immutable
  input/config/functional/partition/state/output/scientific hashes.

Final API additions after controller gate review: `normalizedConstraints`
(VectorXd), `centralInertiaEigenvalues` (Vector3d), `centralInertiaPsd` (bool),
`normalizedInertiaTolerance` (double, 128*binary64 epsilon). Existing raw
moment fields are fixed Eigen ten-vectors; centralInertia is in the fixed
paper coordinate frame. Header and documentation describe candidate versus
validated ownership and all pre-seal limitations.

## TDD and verification evidence

Every command ran in the allowed working checkout with:

```sh
ulimit -v 25165824
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
set -o pipefail
```

Build command (all iterations, maximum -j3):

```sh
timeout 3600 cmake --build build/task-5c-no-faster-20260909 \
  --target spin-it-2014-source-cage-binary-material-regression -j3
```

Direct focused test command (all iterations):

```sh
timeout 3600 build/task-5c-no-faster-20260909/spin-it-2014-source-cage-binary-material-regression
```

Both outputs were piped through `tee` into the task-5y log files in this SDD
directory, preserving stderr and actual failures. The final wrapper explicitly
exits if the build pipeline fails, before attempting CTest.

- **Initial RED:** `task-5y-red-build.log` and `task-5y-red-test.log`.
  Tests were written before the production header/source. The test executable
  compiled and exited 1: `FAIL: actual deformed binary material admission and
  validation layer is absent`. The absence fallback makes the missing feature
  an executable expected failure instead of an unresolved-header error.
- Interim raw logs: `task-5y-build-1.log` through `task-5y-build-5.log`,
  `task-5y-test-2.log` through `task-5y-test-4.log`. Build 2 exposed missing
  Eigen Geometry definitions for long-double cross product at link time;
  adding its standard header fixed that. Tests 2 and 3 exposed an incorrect
  fixture assumption that zero/one endpoints retain eight beta variables:
  real adaptation merges them to one. The test now uses the actual endpoint
  layout. No production admission gate was weakened. Build 1 reused the
  pre-header absence-fallback object until the test source changed; no pass
  claim was based on that cached executable.
- **First GREEN:** `task-5y-test-4.log`, exit 0, includes zero/.25/.5/.75/one
  stationary validations and the full adversarial fixture checks.
- **Normalization RED:** `task-5y-normalization-red-build.log` and
  `task-5y-normalization-red-test.log`, exit 1 at `fixed input density
  normalization must reject despite small raw residual`. An actual mixed
  endpoint at rho=.001 has raw constraints below1e-6 but fixed-input normalized
  constraints above1e-6. This demonstrates the raw-only gate bug directly.
- **Normalization GREEN:** `task-5y-normalization-green-build.log` and
  `task-5y-normalization-green-test.log`, exit 0. Includes the low-density
  rejection, fixed-frame inertia spectrum and normalized-vector checks, and
  retention of a known constraint failure when the geometry cap is exhausted.
- **Final verification:** `task-5y-final-build.log` and
  `task-5y-final-ctest.log`, command:

```sh
timeout 3600 ctest --test-dir build/task-5c-no-faster-20260909 \
  -R '^spin-it-2014-source-cage-binary-material$' --output-on-failure -j1
```

Final result: final build exit 0; CTest **1/1 PASS**, zero failures, 85.48
seconds test / 85.49 seconds total. The assigned build was explicitly released
to the controller immediately after this run. The controller owns final
all-target build, broader related suite and independent rational/BBW replay;
the implementer did not duplicate that work.

Coverage: actual zero/one/fractional/mixed beta and exact-half tie; nonidentity
all-handle deformation; nonzero phi and nonidentity fixed principal frame;
known stationary success; genuine rounding-induced constraint failure;
low-density normalization distinction; pending, nonconverged and resource
endpoint non-admission; actual geometry rejection; final nonmanifold union
rejection even though primitive geometry checker passed; nonfinite/invalid
tolerances and invalid/overflowing budgets; deterministic repeat; immutable
copied-facts/move/state-delta/owned functional lifetime. Stationary fixtures
are KNOWN_SOLUTION_VALIDATION only, with zero objective weights and an explicit
fixed principal frame. They do not assert an optimization achievement.

## Files and self-review

Only these production/test/doc files are to be committed:

- include/spin2014/SourceCageBinaryMaterial.h
- src/spin2014/SourceCageBinaryMaterial.cpp
- tests/spin_it_2014_source_cage_binary_material_regression.cpp
- docs/spin2014_source_cage_binary_material.md
- CMakeLists.txt (one scoped source/header and test registration)

Self-review read the complete new implementation, header, regression and
documentation. Checked actual deformed-only emission; deterministic cancellation
and indexing; density/root-derived ten-moment and constraint scales; fixed-frame
Givens formulas; PSD and positive mass; failure precedence; candidate/private
ownership; pre-allocation arithmetic and diagnostic bounds; all digests; and
the bounded pre-seal limitation wording. No upstream solver or rest geometry
implementation was changed. All foreign dirty Faster files and unrelated
untracked experiments remain unstaged and untouched. No branch/worktree,
dependency, download, merge, push, tag or protected-main operation occurred.

Known bounds: exact final region partition, full cavity/web thickness and
Python sealing remain subsequent gates. Native CGAL memory is opaque and
subject to the external process cap; the owned/diagnostic budgets explicitly
do not claim an RSS guarantee. Exact predicates can return unknown/resource
failure without publishing a product. No unresolved implementation concern
was identified in self-review.
