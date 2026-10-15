# Reconstructed Heart G8 attempt — 2026-09-23

## Scope

This is a reconstructed method-validation attempt, not an author Heart model or
paper-case reproduction. The G4 parent geometry, objective, target, density,
wall setting, initial fill/angle, solver tolerances, iteration limit and leaf
cap were held fixed; source-geometry and optimizer maximum levels were raised
to 8 as preregistered. No lower-resolution result is upgraded to G8.

## Launcher repairs and retained attempts

The scientific configuration was not changed across launcher repairs. Three
fresh directories contain pre-optimization failures and are retained:

1. `build/reconstructed-heart-g8-20260923-run1/launcher-failure.json` —
   Python package import failed before runner initialization.
2. `build/reconstructed-heart-g8-20260923-run2/failure.json` — `begin_run`
   rejected an empty `dependency_roots` map.
3. `build/reconstructed-heart-g8-20260923-run3/failure.json` — evidence
   validation rejected tuple-valued identity fields before native solving.

The runner now records the actual pybind11 include root and normalizes evidence
parameters through the project's strict JSON parser. Tests first reproduced
the missing helper/normalization behaviors, then passed after the repairs.
The native optimization was launched only in `run4`.

## Scientific result

Run directory: `build/reconstructed-heart-g8-20260923-run4/`.

- `source.stats.maximum_visited_level = 8`.
- The six optimizer passes reached maximum active/output levels 4, 5, 6, 7,
  8, and 8 respectively.
- Passes 0–4 reported `CONVERGED` under their local first-order KKT test.
- Pass 5 reached its configured 160-iteration limit:
  `SLQP reached its iteration limit`.
- The last pass was operating at level 8 with 3,520 level-8 leaves, of which
  2,776 remained fractional. Its independently recorded KKT witness has
  equality residual `1.48e-14`, bound residual `0`, and stationarity/full
  residual `1.05e-5`—above the frozen `1e-7` projected-gradient tolerance.
  The trace contains 160 optimizer iterations plus a `final_state` KKT snapshot
  at index 160. The snapshot's `accepted_step=0` and zero line-search count are
  bookkeeping fields, not a rejected line search. The actual final optimizer
  iteration (159) accepted a full step after one line-search evaluation.
  Across the final 40 actual iterations, the active-membership hash was stable,
  every step was accepted (minimum step 0.5), and no restoration or first-order
  fallback was used; KKT infinity norm ranged from `3.81e-6` to `7.00e-5`.
  The run's direct failure cause is therefore the fixed 160-iteration cap being
  reached before the first-order criterion, not a terminal line-search failure.
  These observations do not establish a specific SLQP implementation defect.
- Overall status: `RELAXATION_UNKNOWN`; terminal relaxation KKT does not apply.
- No binary check or binary material was produced; no material artifact was
  sealed. The run is not a successful G8 optimization or a manufacturing result.
- The independent repeat was not launched because the preregistered repeat
  gate only applies to a valid first result. No solver setting was changed to
  force convergence.

The original run summary's maximum-level parser accepted only JSON objects,
while the native solver emits per-level counts as arrays. A new failing test
reproduced that reporting bug; the parser now handles both representations.
Independent inspection of the untouched native output confirms level 8 was
actually reached before the later pass failed. Because the runner itself is a
bound evidence input, correcting its parser makes the old run's strict
freshness check report `STALE_RESULT`; the historical `evidence.json` and all
large trace/output files remain unchanged. The run is retained as diagnostic
evidence only.

## Verification and limits

Focused runner tests: **4/4 passed** after the fixes. The G8 attempt used a
24 GiB cgroup memory ceiling, zero swap, and a 60-minute timeout; measured
peak RSS was about 0.66 GiB, with zero swaps. This removes resource exhaustion
as the cause of the iteration-limit result, but does not establish why the
last SLQP pass failed to meet its termination conditions. Active membership
changed six times (iterations 4, 5, 8, 9, 104, and 105); it was stable for the
last 55 actual iterations. The available trace is insufficient to distinguish
slow convergence from a model/derivative or numerical conditioning limitation.

The configured 282-test CTest suite was run with `-j4`: 281 passed and
`spin-faster-dependency-contract` failed in that parallel invocation. The
same registered test passed in two subsequent isolated invocations. Since the
parallel failure output was superseded by the later CTest log and its exact
cause was not recovered, the full suite is recorded as **not fully green**;
the isolated passes do not retroactively change the 281/282 result.

Next valid action is a targeted trace/partition diagnosis, not an automatic
same-configuration rerun. The outcome remains `PARTIAL_BLOCKED`; no author-case
credit is granted.
