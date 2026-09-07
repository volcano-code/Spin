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
