# Task 6a report: complete native Faster attempt/evaluation evidence

## Status and scope

Implemented the bounded Task 6a instrumentation/accounting change on base
`7057c8120997c7e0d70b8d1a11833822be3c9ed9`. The spherical Newton equations,
step/reduction constants, tolerances, initial-guess selection, known-solution
selection, and critical-parameter semantics were not changed. Existing
scientific files and historical artifacts were not rewritten.

This work does **not** claim full Eq. 14 verification, original-paper case
reproduction, a PaperRunEvidenceV2 grant, or a full-build attestation.

## Files changed

- `include/spin/SpinDesignRequest.h`
- `include/spin/QuadricHollowingSolver.h`
- `src/spin/QuadricHollowingSolver.cpp`
- `include/api/BalanceModel.h`
- `src/api/BalanceModel.cpp`
- `bindings/pybind_module.cpp`
- `python/make_it_stand/spin.py`
- `tests/quadric_hollowing_solver_regression.cpp`
- `tests/python/test_spin_it_faster.py`
- this report and the retained Task 6a test/build logs

## New public contract

`SpinNewtonResult` now exposes event stream version 1, the full ordered
`SpinSolverEvent` stream, a separate SHA-256 `eventHash`, total actual Eq. 14
evaluator calls, per-restart evaluator calls, evaluator wall seconds, and a
per-call `SpinEvaluationTiming` sidecar. `QuadricHollowingResult` propagates
those fields and separately reports post-solve material-export evaluation
count and timing.

Each restart has explicit `attempt_start` and `attempt_terminal` events. Each
actual evaluator entry has one `evaluation` event with global and per-restart
indices, restart/iteration, optional trial step, attempted normalized quadric
coefficients and quadric hash, outcome, stable reason, and diagnostic detail.
Successful evaluations retain optional scientific metrics and actual
geometry/topology hashes even when rejected. Failed evaluations use explicit
absence for unavailable scientific values; no NaN/Infinity sentinel is used.
Standard evaluator exceptions become failed events; non-standard exceptions
are rethrown.

Terminal reasons implemented are `initial_evaluation_failed`,
`initial_sensitivity_vanished`, `converged`,
`bordered_system_rank_deficient`, `newton_direction_vanished`,
`maximum_iterations_reached`, `maximum_evaluations_reached`, and
`line_search_exhausted`. A later successful restart clears earlier diagnostic
text.

The existing `trace` and `traceHash` keep their accepted-state-only v1
semantics. Serialized artifacts explicitly publish
`trace_hash_semantics: legacy_accepted_states_v1`; no old trace is upgraded.
The new output hash schema is `QuadricHollowingResult-v2` and binds the separate
event hash.

Native artifact getters and Python `SpinItFasterResult` expose the event hash,
events, solve-only evaluator count/wall time, evaluation timing sidecars, and
post-solve material-export count/wall time. Native result fields additionally
separate evaluator-reported component totals for the solve and post-solve
material export. JSON event strings are escaped, absent values serialize as
`null`, and Python tests parse both success and failed-evaluation artifacts
with `allow_nan=False`.

## Deterministic hash coverage

`SpinSolverEvents-v1` hashes the complete ordered scientific event stream:
event type, both evaluation indices, restart/iteration, optional step, outcome,
stable reason, diagnostic detail, optional attempted coefficients, attempted
quadric hash, optional metrics, and optional geometry/topology hashes. Every
variable-length string is length-prefixed. Timing sidecars are deliberately
excluded. A controlled evaluator repeat with different evaluator-reported
timings produces byte-equal scientific events and the same event hash.

## TDD evidence

### Focused baseline

Command:

```text
cmake --build build/spin-faster-python --target quadric-hollowing-solver-regression -j2 && ctest --test-dir build/spin-faster-python -R '^spin-faster-newton$' --output-on-failure
```

Observed before test edits: `1/1` passed, `100% tests passed, 0 tests failed`.
The configure/build was exact-geometry-OFF and used the cached historical
`/home/lin/xvxinyao/...` source spelling, so it is recorded only as the focused
baseline, not the final all-feature attestation.

### RED

Command:

```text
cmake --build build/spin-faster-python --target quadric-hollowing-solver-regression -j2 2>&1 | tee /tmp/task6a-red.log
```

Observed expected compile failure in the real native regression: missing
`spin::SpinSolverEvent` and missing `eventStreamVersion`, `events`, `eventHash`,
`evaluationTimings`, and `restartFunctionEvaluations` fields on
`SpinNewtonResult`. Full retained output:
`task-6a-red.log`.

### GREEN and intermediate fixture correction

Initial native GREEN command:

```text
ctest --test-dir build/spin-faster-python -R '^spin-faster-newton$' --output-on-failure
```

Observed `1/1` passed. After adding a rank-terminal regression, the combined
focused run correctly failed `spin-faster-newton`; retained as
`task-6a-rank-fixture-failure.log`. Investigation showed the proposed rank
fixture supplied zero raw mass and therefore correctly terminated earlier as
`initial_evaluation_failed` with `raw moments require positive mass/volume`.
The fixture—not production behavior—was corrected to keep positive raw mass
while constructing the dependent bordered system. The direct native binary
then passed and reached `bordered_system_rank_deficient`.

Final focused command on the final source:

```text
cmake --build build/strict-v2-spin --target quadric-hollowing-solver-regression _make_it_stand -j2 && ctest --test-dir build/strict-v2-spin -R '^(spin-faster-newton|spin-faster-python-api)$' --output-on-failure
```

Observed `2/2` passed, `100% tests passed, 0 tests failed`. A verbose Python
rerun observed all 6 unittest methods `ok`, `Ran 6 tests`, `OK`, with no skip.
Full focused log from the earlier equivalent build is retained as
`task-6a-green-focused.log`; final-source results are also present in the full
suite log.

Native cases independently count evaluator invocations and cover initial throw
then successful restart, all initial evaluations failing, finite rejected
trial, exception trial with later success, zero sensitivity, evaluation budget,
iteration budget, deterministic rank deficiency, line-search exhaustion,
non-standard exception propagation, deterministic event equality/hash despite
timing changes, and preservation of the analytic solution and legacy accepted
trace behavior.

## Final verification on final source

All-feature configure: `build/strict-v2-spin`, `BUILD_TESTING=ON`, Python ON,
Spin-It Faster ON, exact geometry checks ON. Inventory: 181 registered tests.

Build command:

```text
cmake --build build/strict-v2-spin -j2
```

Observed exit 0. Full output: `task-6a-full-build.log`.

Relevant Faster suite command:

```text
ctest --test-dir build/strict-v2-spin -L spin-it-faster --output-on-failure
```

Observed `13/13` passed, `100% tests passed, 0 tests failed`, 19.09 seconds.
Full output: `task-6a-faster-suite.log`.

Full registered suite command:

```text
ctest --test-dir build/strict-v2-spin --output-on-failure
```

Observed `181/181` passed, `100% tests passed, 0 tests failed out of 181`,
181.06 seconds, with no skipped or Not Run tests. This run followed the final
all-target build, which refreshed every executable affected by the public
header changes. An earlier 182.39-second run occurred before that final
all-target relink and is not used as final-source evidence. Full output:
`task-6a-full-suite.log`.

`git diff --check` produced no output.

## Numerical behavior and remaining gaps

Only accounting, diagnostic, public serialization, timing separation, and hash
identity were intended to change. The existing analytic identity-map solution,
accepted residual monotonicity, accepted trace determinism, and legacy trace
hash repeatability remain covered. Corrected accounting now counts initial
throws and every attempted line-search evaluation against the restart budget.

No deterministic native-seam fixture was found for the
`newton_direction_vanished` branch; its terminal instrumentation is present but
that branch is not claimed as regression-covered. The known one-candidate
buoyancy-to-`Infeasible` artifact-adapter issue is unchanged and remains outside
Task 6a. The separate controller derivative check reported historical error
`0.00658869` against threshold `0.08`; inverse-Gauss/discrete-derivative
mathematical verification remains a later Task 6 concern and was not changed
here.

## Self-review

Reviewed the final diff for evaluator-entry count placement, budget boundaries,
one timing row per returned result, optional failure metrics, exception
propagation, stale diagnostics, terminal-event completeness, string/hash
ambiguity, timing exclusion from event/output scientific identity, JSON
escaping/null behavior, legacy trace preservation, Python immutability, and
the two post-solve export calls. No production/test source changed after the
final Faster and 181-test suite runs.
