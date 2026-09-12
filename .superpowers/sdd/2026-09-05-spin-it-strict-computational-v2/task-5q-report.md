# Task 5q — actual source/full-BBW material moment derivatives

Base: `6f61946cbb2b7d3329321f6741c4667ad616042f`, branch
`spin-it-strict-computational-v2`. Global status remains `PARTIAL_BLOCKED`.
This report covers the completed analytic material moment layer, not a joint
solver or geometry/paper completion. Implementation commit:
`e6fc2d07549a9b96ca699a845670b14381b0dc2a` (nine explicit owned paths).

## Implementation

Added immutable `SourceCageMomentFunctional(boundary, density, config)`, with
`SourceCageMomentFunctionalConfig` and `SourceCageMomentEvaluation`. Parameters
are exactly `[scale, tx, ty, tz]` per actual original cage handle, all rotations
explicitly identity, initial scale one and translation zero. Each beta is one
actual adaptive leaf's hollowed fraction. The result includes final native raw
10 moments, all beta and handle columns, moved original-plus-material rows,
evaluation-mode flag and a versioned scientific identity.

The retained successful BBW final field is interpolated at cached **rest**
eight-node embeddings, including boundary face centers and hanging points.
Only one full handle-weight row is allocated as scratch. No global P-by-H or
P-by-4H matrix is formed, and no weights are normalized, clamped or selected.
The private supplement Algorithms 1/2 triangle kernel is shared with both
existing whole-surface APIs; their closed-oriented-chain validation remains.
The value-only path uses identical polynomial values and geometry, with absent
Jacobians represented by `hasJacobians=false` and two `0 x 0` matrices.

Sign convention was clarified with the controller before implementation:
canonical triangles point along the positive axis. The coefficient **added to
outer** is `beta[positiveOwner] - beta[negativeOwner]`. This is exactly outer
minus each actual primitive's outward boundary, not leaf bounding boxes.
Density multiplies all accumulated moments and derivatives once.

Counts, derived byte sizes, Eigen indices and upstream applicable budgets are
checked. The configured Jacobian cap covers `10B+40H` combined; value-only does
not require that capacity. Every triangle is checked even for zero coefficient.
State/density and collapsed geometry errors reject; nonfinite arithmetic and
nonpositive final mass reject. Resource failures use `RESOURCE_LIMIT`
`length_error`. Hashes bind method/version/layout, identity rotations, density,
budgets, counts, upstream actual boundary and BBW, input state, evaluation mode,
and all produced scientific arrays. Caller-visible output is separate from
immutable retained input.

## Actual TDD and focused checks

The new native target was registered and built **before** its production
interface or kernel edits. Command (all build commands use `ulimit -v 25165824`
and at most `-j3`):

```
cmake --build build/task-5c-no-faster-20260909 \
  --target spin-it-2014-source-cage-moments-regression -j3
```

`task-5q-red.log`: exit 1, test translation unit line 1 could not include the
not-yet-created `spin2014/SourceCageMomentFunctional.h`. This is the requested
missing-interface RED, not a typo or an unrelated dependency failure.

Focused native GREEN is retained in `task-5q-final-focused.log` (exit 0).
Tests use the actual baseline outer cube `(0.3,0.7)`, design `(0.34,0.66)`,
source depth 5, adaptive initial 2 / max 5, cage 2 / grid 3 and successful
coupled BBW: **26 handles, 8 beta leaves, 496 real primitives**. Beta zero is
analytic outer material; beta one subtracts exactly `[11/32,21/32]^3`.
Nonuniform beta and all-handle deformation agree with an independent native
signed-tetrahedron integration of the actual primitive boundaries. A second
translated/scaled nonsymmetric-design fixture repeats the derivative checks.

All 112 columns on both fixtures use centered differences at the three
preselected native-unit steps `3e-4, 1e-4, 3e-5`. The fixed entry rule, chosen
before observations, is relative error `<=1e-6` for `abs(J)>1e-7`, otherwise
absolute error `<=1e-11`. No threshold or denominator was changed after results.

| Fixture / columns | h=3e-4 worst relative | h=1e-4 | h=3e-5 |
|---|---:|---:|---:|
| Baseline beta | 1.67491e-10 | 3.66044e-10 | 1.97233e-9 |
| Baseline handles | 7.69325e-8 | 1.25577e-7 | 4.15404e-7 |
| Transformed beta | 4.13260e-10 | 9.45012e-10 | 4.41479e-9 |
| Transformed handles | 1.43882e-7 | 3.56380e-7 | 5.56981e-7 |

The baseline handle worst absolute error improves from `1.26546e-10` to
`1.45642e-11` (8.69x for a 3x step reduction, consistent with h-squared
truncation); the finest step reaches `1.63991e-11`. Beta is analytically linear
and its errors are roundoff-dominated. Other entries, particularly the
translated fixture at its finest step, are also noise-dominated. There is no
claim that every entry improves monotonically.

Additional native checks cover all output dimensions; exact value-only versus
derivative values/points; repeat identity; density scaling; uniform affine
covariance; full one-row weight budget and insufficient Jacobian output budget;
malformed/overflow budgets; beta dimensions/range/NaN; handle
dimensions/scales/NaN/Inf; density; finite arithmetic overflow; collapsed output
with zero beta; actual BBW orientation reversal producing nonpositive final
mass; empty carveable material; and state/config-sensitive identities.

## Independent controller evidence

The controller compiled its separate exporter against stable core archive
`52754814994bfa24dfd84d335e9cfad1c64b9e330cb94510d94254a36a98eb94`.
The retained directory is `build/task-5q-controller-DhEfuX`.
It independently reconstructs eight-node interpolation from actual exported
grid/node data and computes signed-tetrahedron moments and complex-step
derivatives without the production polynomial kernel. Both actual 26-handle,
8-beta fixtures were exported twice, including a nonsymmetric outer/design
case distinct from the native test's transformed case.

`check1.log` and `check2.log` report normalized moment error at most
`5.583e-16`, all-column complex-step Jacobian error at most `3.717e-16`, and
exactly zero double-precision independently reconstructed motion error.
All mixed second moments are checked. Both exports are byte-identical SHA-256
`3c4610f96af31682ea7073d1d269f13292b8a7ce335ee51d191c9cc6933e4617`.
The controller also reports all 12 deliberate evidence corruptions rejected
(`negative.log`). This evidence establishes the derivative layer only.

## External relocation and isolated build recovery

During the task, an external actor moved the worktree to
`/data/xvxinyao/moment-guided-spin-v1/code/make-it-balance-paper-sim` and left
`/data/xvxinyao/make-it-balance-paper-sim` as a symlink. Neither alias nor Git
metadata was changed by this agent. Automatic regeneration in the **assigned**
build changed libigl dependency stamp paths and re-cloned CGAL before the
controller's stop instruction arrived. The generated script removed/replaced
the prior dependency directory; we made no backup and do not claim otherwise.
Eigen's same-time updated stamps suggest it may also have been re-cloned.
CGAL was clean at its exact pinned commit
`f7c3c8212b56c0d6dae63787efc99093f4383415`; Eigen's pinned identity was unchanged.
The interrupted build and partial files were preserved in `task-5q-build.log`.

After read-only controller audit, the same assigned build was configured with
`-DLIBIGL_SKIP_DOWNLOAD=ON`, retaining 2014 ON, exact ON, Faster OFF. The
successful resume logs are `task-5q-resume-configure.log` and
`task-5q-resume-build.log`. The full assigned build then relinked all affected
executables and Python module (`task-5q-safe-build.log`, 109/109 steps).

The controller authorized a narrow dependency-test infrastructure fix: both
fresh configure child invocations now pass `-DLIBIGL_SKIP_DOWNLOAD=ON`.
Before changing those calls, the actual generated isolation fixture was made
to reject a missing flag. Its RED was observed: exit 8 and the retained child
error `fixture dependency downloads were not disabled`.
`task-5q-download-guard-red.log` and its retained fixture directory preserve
that evidence. Adding the two flags produced isolation GREEN, 1/1 in 7.32 s
(`task-5q-download-guard-green.log`). Existing feature, scratch-isolation,
concurrency, child-status and runtime checks are unchanged.

The first safe suite then exposed the external alias in existing Python
evidence capture, whose no-symlink validator correctly rejected it. This
partial, interrupted run is preserved in `task-5q-safe-scoped-suite.log` and is
not a passing suite. No evidence-validator code or test threshold was weakened.
At the controller's direction, CMake was rerun with canonical `-S` and `-B`
paths in the same build. It accepted the existing cache without deletion or
manual cache editing. CMake's internal historical `CMAKE_HOME_DIRECTORY` still
names the alias, while generated CTest Python paths and Ninja's regeneration
command use canonical paths. `task-5q-canonical-configure.log` records success.
The full canonical rebuild completed 362/362 steps, including HiGHS, all
affected executables and the Python module (`task-5q-canonical-build.log`).
The previously failing `source-spin-material-python` then passed 1/1 in
23.85 s (`task-5q-canonical-python.log`), with no source changes to validation.
The core archive remained byte-identical to the independently checked archive.
The controller additionally relinked its exporter with the canonical build,
repeated both exports and independent checks, and obtained the same exporter
binary and scientific export bytes. Its original and post-move evidence thus
bind the same final production implementation and actual dependency build.

## Final verification, identity and self-review

The final command ran from the canonical worktree, with the same memory and
time bounds, against only the assigned build:

```
timeout 3600 ctest --test-dir build/task-5c-no-faster-20260909 \
  -R '^(spin-it-2014-|spin-2014-|source-spin-material-python$)' \
  --output-on-failure -j1
```

**35/35 PASS, exit 0, 237.47 s**, retained in
`task-5q-canonical-safe-suite.log`. This includes the new source cage moment
regression (1.58 s), existing surface value/differential compatibility and
value-only allocation guard, source/binding/BBW regressions, rebuilt Python
bridge, and both dependency checks. All executables were relinked against the
current core before this suite; no old static linkage is counted as evidence.

The actual retained dependency run is
`build/task-5c-no-faster-20260909/spin-2014-dependency-contract-runs/run-74d7a050dfe48d0e`.
Its contract status is `PASSED`, audited child status `0`, both fresh child
caches have `LIBIGL_SKIP_DOWNLOAD=ON`, and both have Faster/exact OFF. The
enabled child has 2014 ON; its actual generated Ninja graph contains neither
`SourceCageMomentFunctional` nor CGAL. Graph SHA-256:
`e9abb084824a4edc82a41a655111fdaca359d539dee0faac77e043f5118071fa`.

| Retained evidence | SHA-256 |
|---|---|
| `task-5q-red.log` | `409bee7404affd863fd51bc10099886f6f03836324274578c18491440ae0363f` |
| `task-5q-final-focused.log` | `286b9ad1b0f032ac74f9979d29f6898d0685767fca3960a4ec9bdb328c226bb5` |
| `task-5q-download-guard-red.log` | `a918cec58ab68a69b88ab7f8ad28586b2a5e1604e7eb70a6e2e9252a399ca9b0` |
| `task-5q-download-guard-green.log` | `7fd664f1a001a3010c737275a1baa9c0164596cc17ff813d4c0138feaac883aa` |
| `task-5q-canonical-python.log` | `b8b91820c413ae3012911e10d5ed7cbff7920f04c4de49b1224b488c3d37d6d3` |
| `task-5q-canonical-safe-suite.log` | `279d8f67221c9bcc9d6c5534cdda5c94dc3c8f85c2d88dbc1424929237b04f13` |
| Final `libmake_it_stand_core.a` | `52754814994bfa24dfd84d335e9cfad1c64b9e330cb94510d94254a36a98eb94` |

Self-review covered canonical ownership signs, native row ordering, coordinate
chain rule and density application, fixed rest embeddings/all handle order,
zero-coefficient triangle checks, immutable snapshots, every retained/dynamic
allocation's size checks, value-only absence, exception categories, hash inputs,
the CMake feature conjunction, and the complete owned diff. `git diff --check`
and staged diff checks passed. The new functional implementation is 238 lines;
constructor/embedding/digest helpers remain separated and readable. No Task 5p
readability refactor or upstream source/cage/binding changes were included.
No unresolved algorithmic concern was found. Resource caps bound value/count
allocations; they do not claim an exact process RSS bound.

No main strict build, external Faster/Plane/Polyhedral/Quadric files, external
experiments, solver implementation or upstream source/cage/binding code was
edited or staged. The build is owned only at `build/task-5c-no-faster-20260909`.
No subagents, push, merge, tag, or manual cleanup was used.

Positive mass, finite triangle geometry and an analytic derivative do not
establish deformed inversion/intersection/overlap legality, minimum shell
thickness, exact final partition/seal, physical feasibility, regularizer/SLQP
completeness or an original-paper reproduction. Every fixture here remains
`SYNTHETIC_REGRESSION`. The remaining strict scientific gates are unchanged.
