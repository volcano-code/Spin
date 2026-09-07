# Task 6c report — reject non-finite Faster numerical states

## Scope and result

This bounded change starts from
`03dccd0cd9098e3d6ec699e0f8d754e1adbe34ea` and fixes the identified
finite-input false convergence in the Faster spherical Newton solver. It does
not claim completion of Faster, the original paper method, or a release.

The solver now fails closed when finite evaluator inputs overflow while
constructing or consuming a numerical state. Invalid initial evaluations are
retained as failed events and permit the next deterministic restart. Invalid
line-search evaluations consume the existing budget, remain failed events,
and permit the same line search to try a smaller step. No invalid state enters
the accepted trace or produces `CONVERGED`.

## Root cause

`evaluateState` checked only the raw moment vector, raw moment Jacobian, final
normalized/physical residual scalars, and multiplier. The derived
`SpinMomentEvaluation` was not validated. In the existing Fidget Spinner
fixture, finite `Qxx=Qyy=1e308` produced objective
`-(Qxx + Qyy) = -infinity`, while the objective gradient, equality residual,
and positive multiplier remained finite. The convergence predicate therefore
accepted the state. Task 6a only made the unavailable objective safe in the
event representation; it intentionally did not change this solver behavior.

Related construction gaps allowed finite entries to overflow aggregate norms
and allowed a finite objective Hessian and moment Jacobian to overflow their
matrix product before LU factorization.

## Files changed

- `src/spin/QuadricHollowingSolver.cpp`
- `tests/quadric_hollowing_solver_regression.cpp`
- this report and the separately named Task 6c raw logs

No public header, status enum, tolerance, model input, restart limit, iteration
limit, line-search reduction, geodesic-step cap, Python API, historical output,
or main checkout was changed.

## RED

The regression was changed before production code. It uses only the existing
public `SpinMomentEvaluator` seam.

```text
cmake --build build/strict-v2-spin --target quadric-hollowing-solver-regression -j4
ctest --test-dir build/strict-v2-spin -R '^spin-faster-newton$' --output-on-failure
```

The first run failed for the actual legacy defect:

```text
[FAIL] scenario=quadric-hollowing-solver overflowed initial objective was accepted
0% tests passed, 1 tests failed out of 1
```

After aggregating the boundary fixtures, the second diagnostic RED showed:

```text
[RED] line-search numerical state outcome=accepted reason=residual_decreased objective_present=0
[RED] invalid residual normalization state was not rejected at evaluation
[RED] invalid moment-Jacobian norm state was not rejected at evaluation
[RED] invalid derived derivative state was not rejected at evaluation
[RED] overflowed Newton system was not rejected before factorization
```

Retained RED logs:

- `task-6c-red-build.log`, `task-6c-red-focused.log`
- `task-6c-red-line-build.log`, `task-6c-red-line-focused.log`
- `task-6c-red-boundaries-build.log`, `task-6c-red-boundaries-focused.log`

The first implementation exposed two test-fixture defects without weakening
the production gates: the nominal Fidget restart evaluator returned zero mass,
and the first Newton-product fixture overflowed the Jacobian norm before
reaching matrix assembly. Those failures are retained in
`task-6c-green-focused.log` and its build log. The corrected fixtures use an
explicit finite-mass Fidget state and individually finite/norm-finite values
whose Hessian-Jacobian product overflows. The corrected focused run is retained
as `task-6c-green-focused-2.log` and its build log.

## GREEN and final verification

The final focused build/test, after fixture refinement and source cleanup:

```text
cmake --build build/strict-v2-spin --target quadric-hollowing-solver-regression -j4
ctest --test-dir build/strict-v2-spin -R '^spin-faster-newton$' --output-on-failure
```

Observed `1/1` passed, zero failures. Logs:
`task-6c-green-focused-final-build.log` and
`task-6c-green-focused-final.log`.

The registered Faster suite was enumerated and then run using its actual label:

```text
ctest --test-dir build/strict-v2-spin -N -L spin-it-faster
ctest --test-dir build/strict-v2-spin -L spin-it-faster --output-on-failure
```

Observed `14/14` passed, zero failures in 18.83 seconds, including native,
Python API, dependency, and miniature paper-runner tests. Logs:
`task-6c-faster-suite-list.log` and `task-6c-faster-suite.log`.

With source frozen, the required final build and full suite were run once:

```text
cmake --build build/strict-v2-spin -j4
ctest --test-dir build/strict-v2-spin --output-on-failure
```

The all-target build exited 0. The full registered suite reported `182/182`
passed, zero failures in 187.84 seconds. A scan of the raw log found no
`Not Run`, skipped, or `SKIP` marker. Logs:
`task-6c-final-all-target-build.log` and
`task-6c-final-full-ctest.log`. `git diff --check` produced no output.

The verified cache has Faster, Spin-It 2014, Python bindings, and exact
geometry checks enabled, with the viewer disabled, as required. Its optional
PyBullet tests are disabled; therefore `182/182` means the complete registered
suite for this configured build, not a simulation-complete or release claim.

## Scientific behavior changes

- All scalar/matrix fields returned by the moment problem and consumed by
  objective reporting, multiplier/residual construction, Newton assembly, and
  convergence are required to be finite.
- Aggregate raw-moment, moment-Jacobian, residual, and normalization quantities
  must remain representable; finite entries that overflow these calculations
  are rejected rather than normalized into a plausible zero.
- The assembled Lagrangian Hessian, bordered Newton matrix, right-hand side,
  and solved direction are checked before downstream acceptance/use.
- Evaluation events for numerical rejection use
  `reason=invalid_numerical_state`; invalid initial attempts terminate with
  `initial_numerical_state_invalid`. Failed event metrics remain absent, while
  finite coefficients and material hashes remain available.
- Evaluation and per-restart counts are unchanged in meaning. Rejected calls
  consume budget. Finite restart and smaller-step line-search recovery are
  covered, as is deterministic hashing of invalid-state events.
- Numerical failures return existing non-infeasibility statuses
  (`GEOMETRY_EVALUATION_FAILED` for all-invalid initial evaluations or
  `OPTIMIZATION_DID_NOT_CONVERGE` for an invalid Newton system). No
  `INFEASIBLE` certificate is manufactured.

## Hashes

Final source and primary build artifacts used by the recorded suite:

```text
8983b40c6006ae4d3d8217b81055d25377c193c7e19fc05130eefaded94f92c1  src/spin/QuadricHollowingSolver.cpp
91d19f0bda4928902ae425d25f9f97ac26c06c765ec387d7fcd085987b9083f4  tests/quadric_hollowing_solver_regression.cpp
c6b25d0ebc4f7c32c9ffd1d0a6266c367933cda7f93ee7637091ea932d1fd96f  build/strict-v2-spin/quadric-hollowing-solver-regression
e4dfb67f535764947ec3b67cb55485328f58fbff271dd281e54827c6c358bd1f  build/strict-v2-spin/libmake_it_stand_core.a
2c3f5a72c82d0ba5d8b5d0cd8869d13d89ea8a21e0daa405e88242e1778e517c  build/strict-v2-spin/CMakeCache.txt
```

The uncommitted two-file production/test diff from the stated base had SHA-256
`d7893d8a86dff4e0832ebc7cce3eb8be4d3b149d77b9e03b34bbc203b84f73ea`
at final verification.

## Self-review and numerical caveats

The final two-file source diff was reviewed for every construction and
acceptance path: raw map, derived problem fields, multiplier solve, residual
and scale construction, initial acceptance, line-search acceptance, Newton
matrix construction, direction solve, event optional fields, counters,
restart recovery, and terminal status/reasons. Mutation checks are represented
by distinct fixtures for derived objective/derivative overflow, normalization
overflow, Jacobian-norm overflow, and a later matrix-product overflow.

This is a fail-closed representability policy, not a proof of Eq. 14, a proof
of convergence, a feasibility certificate, or validation of the paper's full
discrete derivative. It does not fix the separate buoyancy adapter or
inverse-Gauss/discrete-derivative limitations. Very large finite physical
states whose required double-precision norms or Newton algebra overflow are
reported as numerical failure; values are not clamped or rescaled into a
plausible result. There is no separate best-candidate store in this solver;
only accepted finite states enter its legacy trace/final selection.

---

# Task 6c independent-review fix round 1

## Finding and resolution

The independent review correctly found that the original cases labelled
"residual normalization" and "derived objective/derivative" were
mutation-blind. The former contained two `1e308` raw moments and stopped at the
earlier raw-moment norm guard. The latter made `r3/r0` itself infinite and
stopped at the objective guard before independently protecting derivative
checks.

This fix round changes no production code. It adds three discriminating
public-API fixtures to `tests/quadric_hollowing_solver_regression.cpp`, adds
explicit precondition checks proving which earlier quantities remain finite,
and renames the former normalization helper/case to its accurate role:
`overflowedRawMomentNorm` / `raw-moment norm`.

## Discriminating fixtures and preconditions

All fixtures enter through the existing public `SpinMomentEvaluator` and
`SpinMomentProblem` APIs; no test-only production seam was added.

1. **Derived objective Hessian only.** `LowerCenterOfMass` receives
   `r0=1e-110`, `r3=0`, identity moment Jacobian, and initial quadric `-e3`.
   The test directly evaluates the public problem and asserts finite raw-map
   entries/norms, finite objective, finite objective-gradient entries and norm,
   finite constraint fields/Hessians, and a non-finite objective Hessian. The
   automatic-differentiation reciprocal forms `2/r0^3` before the zero
   numerator is applied, so this is a real formula-path Hessian overflow. With
   only the objective-Hessian guard removed, the zero residual and positive
   multiplier allow the invalid state to cross the evaluation boundary.
2. **Objective-gradient norm / equation scale.** `SpinningTop(1e155)` receives
   mass 1, zero first moments, diagonal `Q=(0.5,0.5,0.5)`, identity moment
   Jacobian, and `q=-e3`. The public problem produces a finite objective,
   finite individual gradient and Hessian entries, finite constraint fields,
   and finite raw/Jacobian norms, while the objective-gradient Euclidean norm
   overflows. This reaches the later scale block without an earlier rejection.
3. **Constraint residual norm.** `BuoyantStability(targetMass=1e155,
   horizontalTarget=(0,0))` receives unit mass, otherwise zero raw moments,
   identity moment Jacobian, and `q=-e3`. Objective, gradient norm, Hessian,
   constraint entries/Jacobian/Hessians, and raw/Jacobian norms are finite;
   only `constraints.squaredNorm()` overflows. This independently protects the
   residual-normalization boundary before invalid state acceptance.

An isolated non-finite *gradient entry* with a finite objective Hessian was not
manufactured. For the controllable public quotient family `r3/r0`, gradient
overflow from `1/r0` or `r3/r0^2` necessarily reaches the same or a higher
inverse power in the Hessian (`1/r0^2` or `2r3/r0^3`) first; the linear Fidget
and MaximumMass gradients are fixed and finite. There is no public arbitrary
`SpinMomentEvaluation` injection seam. The test therefore retains the honest
combined objective/derivative case and adds the discriminating finite-entry
gradient-*norm* case above rather than inventing fault injection.

## Mutation RED proof

The tests first passed against the production guards. Three temporary
production mutations were then applied one at a time with `apply_patch`, each
followed by a target rebuild and the real registered native test. Each mutant
was an expected RED, and the production line(s) were restored with
`apply_patch` immediately afterward.

### Objective-Hessian guard mutant

Removed only `!state.problem.objectiveHessian.allFinite()`:

```text
cmake --build build/strict-v2-spin --target quadric-hollowing-solver-regression -j4
ctest --test-dir build/strict-v2-spin -R '^spin-faster-newton$' --output-on-failure
[DETAIL] invalid derived objective Hessian state was not rejected at evaluation
0% tests passed, 1 tests failed out of 1
```

Logs: `task-6c-fix1-red-hessian-mutant-build.log` and
`task-6c-fix1-red-hessian-mutant.log`.

### Gradient-norm/equation-scale guard mutant

Removed only the `objectiveGradientNorm` and `equationScale` non-finite
clauses. These are paired because `equationScale` is constructed directly from
the gradient norm:

```text
[DETAIL] invalid objective-gradient norm state was not rejected at evaluation
0% tests passed, 1 tests failed out of 1
```

Logs: `task-6c-fix1-red-gradient-norm-mutant-build.log` and
`task-6c-fix1-red-gradient-norm-mutant.log`.

### Constraint-residual validation mutant

The constraint path has intentional redundant checks: non-finite
`constraintSquared` also makes `combinedSquared`, normalized residual, and
physical residual non-finite. The targeted mutation therefore removed that
four-clause residual-validation cluster together. A one-clause mutation cannot
bypass the later guards and would not test acceptance behavior.

```text
[DETAIL] invalid constraint residual norm state was not rejected at evaluation
0% tests passed, 1 tests failed out of 1
```

Logs: `task-6c-fix1-red-constraint-residual-mutant-build.log` and
`task-6c-fix1-red-constraint-residual-mutant.log`.

After every mutant, `src/spin/QuadricHollowingSolver.cpp` returned to SHA-256
`8983b40c6006ae4d3d8217b81055d25377c193c7e19fc05130eefaded94f92c1`.
`git diff -- src/spin/QuadricHollowingSolver.cpp` was empty. The final restored
target rebuild prevents any mutant object/library from leaking into GREEN
verification.

## Fix-round GREEN verification

Restored and final focused runs:

```text
cmake --build build/strict-v2-spin --target quadric-hollowing-solver-regression -j4
ctest --test-dir build/strict-v2-spin -R '^spin-faster-newton$' --output-on-failure
```

Observed `1/1` passed, zero failures. Final logs:
`task-6c-fix1-final-focused-build.log` and
`task-6c-fix1-final-focused.log`. The earlier restored check is retained as
`task-6c-fix1-green-restored-build.log` and
`task-6c-fix1-green-restored-focused.log`.

Relevant native/Python Faster suite:

```text
ctest --test-dir build/strict-v2-spin -L spin-it-faster --output-on-failure
```

Observed `14/14` passed, zero failures in 18.65 seconds. Log:
`task-6c-fix1-final-faster.log`.

Fresh final registered suite after production restoration and the final test
build:

```text
ctest --test-dir build/strict-v2-spin --output-on-failure
```

Observed `182/182` passed, zero failures in 183.96 seconds. The raw log contains
no `Not Run`, skipped, or `SKIP` marker. Log:
`task-6c-fix1-final-full-ctest.log`.

The original all-target build and full-suite evidence at commit `f352332`
remains valid because production bytes did not change. No new all-target build
was required for this test-only review correction. PyBullet remains disabled
in this configured cache, so neither the original nor fix-round registered
suite is a simulation-complete or release claim.

## Fix-round hashes and self-review

```text
8983b40c6006ae4d3d8217b81055d25377c193c7e19fc05130eefaded94f92c1  src/spin/QuadricHollowingSolver.cpp
b7249cc37db0a5d26679a8c460e9a8b6be641fdc97fc0c460c2070a4aed16a25  tests/quadric_hollowing_solver_regression.cpp
e92ad7054d9dd325faa4c9a80bc07c2d4235a0c3f942b3782749f0ca56b66bab  build/strict-v2-spin/quadric-hollowing-solver-regression
e4dfb67f535764947ec3b67cb55485328f58fbff271dd281e54827c6c358bd1f  build/strict-v2-spin/libmake_it_stand_core.a
```

The final test-only diff from `f352332` has SHA-256
`ae6defdafdd0bcad29a6a27555781f0f73af18569bb44e96556f27db2799f261`.
Self-review confirmed that each new case asserts its gate-order preconditions,
uses a supported problem formula and the real solver, and requires the same
failed-event/empty-trace/terminal-reason behavior as the original numerical
acceptance regressions. Existing cases and RED evidence were retained. No
threshold, model input, solver status, production API, or production byte was
changed in this fix round.
