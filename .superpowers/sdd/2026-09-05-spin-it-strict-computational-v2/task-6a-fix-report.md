# Task 6a review fix report

## Scope and result

This is the scoped follow-up to `task-6a-independent-review.md` for commit
`063575e94cae58a55c2c9805678006475f34aae2`. It fixes only the two confirmed
Important representation/integrity defects:

1. Python event coefficients are converted to a fresh, validated, finite
   10-tuple before constructing the frozen `SpinSolverEvent`. Evaluation-event
   coefficients can no longer be mutated or alias the cached diagnostics tree;
   terminal-event absence remains `None`.
2. Optional scientific event metrics are engaged only for finite values. The
   reviewer overflow fixture now retains the pre-existing legacy
   `CONVERGED` result while representing its `-inf` objective as unavailable,
   so native optional availability, event hashing, JSON `null`, and Python
   `None` are consistent.

No Eq. 14 equation, tolerance, step, restart logic, accepted legacy trace,
known-solution selection, or critical-parameter semantic was changed. The
legacy overflow false-convergence remains a separate solver correctness
blocker and is not claimed fixed.

## Files changed

- `src/spin/QuadricHollowingSolver.cpp`
- `python/make_it_stand/spin.py`
- `tests/quadric_hollowing_solver_regression.cpp`
- `tests/python/test_spin_it_faster.py`
- this report and eight separately named fix-round logs

The original Task 6a report and logs were not edited in this fix round.

## RED

Native build and regression:

```text
cmake --build build/strict-v2-spin --target quadric-hollowing-solver-regression -j2
ctest --test-dir build/strict-v2-spin -R '^spin-faster-newton$' --output-on-failure
```

Observed expected failure:

```text
[FAIL] scenario=quadric-hollowing-solver scientific event contains a non-finite metric
0% tests passed, 1 tests failed out of 1
```

Python regression:

```text
ctest --test-dir build/strict-v2-spin -R '^spin-faster-python-api$' --output-on-failure
```

Observed expected failure in the real native-backed artifact test:

```text
AssertionError: [1, 0, 0, 0, 0, 0, 0, 0, 0, 0] is not an instance of <class 'tuple'>
FAILED (failures=1)
```

Retained logs: `task-6a-fix-red-build.log`,
`task-6a-fix-red-native.log`, and `task-6a-fix-red-python.log`.

## GREEN

Focused build and test:

```text
cmake --build build/strict-v2-spin --target quadric-hollowing-solver-regression _make_it_stand -j2
ctest --test-dir build/strict-v2-spin -R '^(spin-faster-newton|spin-faster-python-api)$' --output-on-failure
```

Observed `2/2` passed, `100% tests passed, 0 tests failed`, 0.16 seconds.
Retained logs: `task-6a-fix-green-build.log` and
`task-6a-fix-green-focused.log`.

The native test uses the reviewer's actual finite-input Fidget Spinner map
(`mass=1`, `r4=r5=1e308`, identity Jacobian, normalized `q4=q5`) and verifies:

- all present event metrics are finite;
- the overflowed objective optional is absent;
- the evaluator event is still present; and
- the pre-existing `CONVERGED` status is unchanged, preventing an accidental
  algorithmic fix in this instrumentation-only patch.

The Python test verifies the real failed-evaluation artifact provides a tuple
of exactly 10 finite coefficients, element assignment raises `TypeError`, a
mutated diagnostics copy cannot change either the event view or later
diagnostics snapshots, event and quadric hashes remain stable, and terminal
coefficient absence remains `None`.

## Final-source verification

All-target build after the focused GREEN:

```text
cmake --build build/strict-v2-spin -j2
```

Observed exit 0 and all affected dependents relinked. Retained log:
`task-6a-fix-full-build.log`.

Relevant Faster suite:

```text
ctest --test-dir build/strict-v2-spin -L spin-it-faster --output-on-failure
```

Observed `13/13` passed, zero failures, 20.89 seconds. Retained log:
`task-6a-fix-faster-suite.log`.

Full registered suite, after the all-target build and with sources frozen:

```text
ctest --test-dir build/strict-v2-spin --output-on-failure
```

Observed `181/181` passed, `100% tests passed, 0 tests failed out of 181`,
181.47 seconds, with no skipped or Not Run tests. Retained log:
`task-6a-fix-full-suite.log`.

`git diff --check` produced no output.

## Self-review and remaining blockers

The four-file production/test diff was reviewed for tuple detachment,
coefficient length/finite validation, terminal `None` preservation, finite-only
optional metrics, native/hash/JSON/Python availability consistency, and scope
containment. No production or test source changed after the final all-target
build and test suites.

Remaining independent blockers are unchanged: the legacy overflow fixture can
still converge numerically with an infinite objective; the whole-map
inverse-Gauss/discrete-derivative diagnostic remains approximately
`0.00658869`; the known one-candidate buoyancy adapter issue remains open; and
the `newton_direction_vanished` event branch still lacks a deterministic seam
fixture. None is claimed fixed by this representation patch.
