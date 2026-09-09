# Task5b implementation report

Base: `080b1f1d1f275cc873b1ef9317af0e0153e85d50`, isolated branch
`spin-it-strict-computational-v2` in `/data/xvxinyao/make-it-balance-paper-sim`.
This is the approved bounded Task5b increment, not complete Task5/2014
reproduction. Implementation and builds had one writer; the controller owns
independent review. No model assets, historical experiments, reachability
algorithm/pin, dirty main source, release, push, merge, or tag were changed.
TDD and verification-before-completion skills were read and applied.

## Implemented architecture and evidence

One shared SLQP outer loop consumes `EqualityBoxProblem`; both the real
AdaptiveFillFunctional adapter and tiny differentiable analytic fixtures use
it. The default remains LegacyV3. AuditedKktRestorationV1 is an explicit opt-in,
accepted only with PaperGivens at spin/adaptive interfaces. New controls and
fields are appended, preserving positional meanings. Standalone KKT residuals,
sparse pinned-HiGHS witness/restoration LPs, actual L-BFGS history, audit hashes,
and JSON serialization live in focused separate sources. Python and native
JSON expose selected controls and actual state/gradient/Jacobian/multiplier
vectors, full membership and history/trial evidence. Invalid diagnostic
quantities serialize as null with invalid/UNKNOWN status, never zero success.

The witness LP minimizes a stationarity/complementarity max residual using
free lambda; signs match L=f+lambdaᵀc−alphaᵀx+etaᵀ(x−1). A separate evaluator
checks all five components, dimensions and derived overflow. Restoration
minimizes the actual linearized equality L1 residual in beta/trust-step bounds
without the optimization active mask, then requires actual nonlinear L1
Armijo decrease. Its target is unchanged. Invalid initial callbacks throw;
nonfinite or dimension-changing trial callbacks are retained as rejected
trials and cannot produce a successful step. Equality dimensions are stable
throughout a solve. No restoration result is called an infeasibility proof.

Membership changes reset real history before metric use, including changes
with equal active counts. Restoration and new solves reset history. Curvature
pairs use the current lambda at both endpoints and reject nonpositive/invalid
curvature. Every audited return binds the actual final state, all residuals,
witness vectors and controls, including bounded-prefix identity. Elapsed time
is excluded. Adaptive profile/angle carry, binary validation and unsupported
cage/dual consumers retain fail-closed behavior.

The public usage and complete numerical choices are documented in
`docs/spin_it_2014_audited_slqp.md`. New defaults are Delta=1, LP iterations
10000, LP time cap 10 seconds, LP primal/dual tolerance 1e-9; HiGHS simplex,
presolve on, threads=1, parallel=off, seed=0, output=off. All are
IMPLEMENTATION_DEFAULT, not recovered author values. Existing controls remain
120 iterations, 24 trials, memory 8, equality 1e-8, projected 1e-7, active
1e-10, merit 100, reduction .5, Armijo 1e-4. The document also records inherited
numerical cutoffs and allowed ranges. LP assembly is O(n*m), no n-by-n dense
matrix. The existing pinned HiGHS revision was not changed.

## RED/GREEN record (all files are siblings of this report)

1. `task-5b-red-availability.log` was a build setup error (new target not yet
   configured), not behavioral evidence. After configuration,
   `task-5b-red-availability-2.log` proves EqualityBoxProblem was unavailable.
   A thin adapter exposed the existing real engine without audited changes.
   `task-5b-red-behavior.log` and its console capture the fixed actual failure:
   min f=x, c=x−.75, x0=0 returned LINE_SEARCH_EXHAUSTED, x=0, residual=.75.
2. `task-5b-red-audit-api.log` records the absent KKT interface. After the
   modules and policy branch, `task-5b-green-core.log` passed audited release
   and the existing literal legacy SLQP hash, 2/2.
3. Expanded fixed analytic and real Givens tests reached an actual missing
   adaptive guard: `task-5b-red-integration.log`, invalid input accepted.
   Python opt-in availability failed in `task-5b-red-python.log` before facade
   changes. `task-5b-green-integration.log` and `task-5b-green-python.log` pass.
   `task-5b-build-integration.log` additionally records a mistaken `_native`
   target name after the native executable built; `_make_it_stand` is the
   correct target, used successfully in `task-5b-build-python.log`.
4. `task-5b-red-history-api.log` records the missing shared history interface.
   The actual history now rejects zero/negative curvature and clears for
   [lower,free] -> [free,lower], even though both masks have one active member.
   `task-5b-red-boundaries.log` then caught a callback changing its equality
   dimension; stable-dimension validation made the 4-test focused group pass
   in `task-5b-green-boundaries.log`.
5. `task-5b-red-identity.log` is intentionally retained even though that first
   supplementary probe passed: x=0 already changed membership, so it could
   not expose the missing prefix identity. The independent policy-only probe
   at the same interior state (.5,.5) failed in `task-5b-red-identity-2.log`.
   Explicit bounded-prefix hashing fixed it. None of the brief's fixed
   optimization inputs or acceptance tolerances was retuned.
6. `task-5b-red-trial-reason.log` records the absent per-trial failure reason.
   Reasons are now exported and hashed. Final focused run
   `task-5b-final-focused.log` / `task-5b-final-focused-console.log`: 4/4 pass,
   including nine Python unittest methods and all fixed/supplemental cases.
   All intermediate raw build/console logs remain retained; no failing
   scenario was unregistered or silently skipped.
7. The first candidate full suite passed 190/190, but self-review found that
   an appended maximum-iteration snapshot left merit/count/history-size fields
   at defaults while recording correct state/KKT. `task-5b-red-final-snapshot.log`
   fails on the literal fixed x=.75 merit expectation; the unchanged box
   quadratic with maximumIterations=1 fails its final active-count/history
   expectation in `task-5b-red-final-membership.log`. The snapshot now evaluates
   actual merit, counts the recorded membership, and reports the reset metric's
   size. MaximumIterations remains a nonconverged termination status.
   `task-5b-final-2-focused.log` passes 4/4 after this last production fix.
   Earlier logs retain their existing names and are prior-candidate evidence.

The final analytic suite independently recomputes KKT components with literal
C++ loops over the actual exported vectors; Python independently recomputes
each JSON trace row without calling the native checker. Fixed outcomes:
x=.75 and lambda=−1 for required-bound release; duplicate equalities converge
with lambda1+2*lambda2=−1; inconsistent equalities retain min L1=.5 and remain
nonconverged; x²−.25 at zero returns RESTORATION_STALLED; the box quadratic
reaches (.8,0) with alpha_y=.2; mixed prefix/angle quadratic reaches (.4,3).
The real Givens fixture reaches phi=1.5707963268031664 and full KKT
6.3137399822221774e-11, with deterministic repeated states/hashes. Its adaptive
carry and deliberate 1e-14 binary failure remain tested. A one-iteration exit
audits the actual accepted x=.75 but keeps MAXIMUM_ITERATIONS status.

Supplemental adverse cases retain 24 invalid trials, and 24 finite trials
with no actual constraint decrease. The latter uses the explicitly ill-scaled
analytic c=.25−x+1e30*x² at x=0 to exercise bounded backtracking exhaustion,
without changing the mandated fixtures. A fixed 12-variable/5-equality
integer-Jacobian witness LP with iteration limit 1 actually returns
`WITNESS_LP_UNKNOWN: Iteration limit reached`. Supplied-witness tests cover
primal/dual/complementarity violations, angle duals, missing dimensions,
nonfinite objective and derived overflow. Positive unconstrained/lower/upper/
interior-equality witnesses also pass. Stationarity zero with c=.25 reports
full residual >=.25.

## Compatibility, build isolation, resources and final validation

The literal legacy SLQP hash stays
`ba7371831d27b970aeab54b8c1c14f8b411023e16bd06c7fc0c8a9d27f9d9aae`.
For stronger default-Givens evidence, BASE solver source was read directly
from git and compiled with renamed function symbols against the same current
unchanged functional and appended-compatible types. The retained driver
`task-5b-baseline-identity.cpp` compares BASE/current values and hashes using
the same -O2, -frounding-math compiler options. Both hashes are
`13fcb7597f2bc6a0275921674dd5f109f63abeb0013487e5216bfba14aa99eb5`, with exactly
equal beta/phi (`task-5b-default-givens-baseline.log`, repeated after relinking
against the final core in `task-5b-default-givens-baseline-2.log`). This compares the solver
behavior, not a separately rebuilt entire historical application. Default
native JSON writes remain byte-for-byte conditional on the old profile;
no new audit fields are emitted for default legacy or default Givens runs.

Main configuration: 2014/Faster/exact/Python ON; viewer/PyBullet OFF.
`task-5b-final-build.log`: prior-candidate all-target build passes, -j3, 58.14 seconds,
maximum RSS 1,756,220 KiB. Heavy tests are serial. A separate dependency
contract now builds and runs the actual audited executable with 2014 ON,
reachability/Faster/CGAL OFF; `task-5b-isolation.log` and timed console pass
in 91.26 seconds, maximum RSS 787,632 KiB. It exercises both pinned LPs,
so a static archive alone cannot hide a missing HiGHS link.
After the final-snapshot fix, `task-5b-final-2-build.log` records the final
all-target -j3 build: 34.46 seconds, maximum RSS 1,758,032 KiB, exit zero.

The prior-candidate full configured suite (`task-5b-final-full.log` and its
timed console) passed 190/190 in 312.80 seconds, peak RSS 4,230,800 KiB. It
predates the final-snapshot metadata fix and is not final-code verification.
Final-code full configured suite: **190/190 passed**, zero failed or skipped
tests, 311.74 seconds, maximum RSS 4,253,072 KiB. No failing scenario was
removed from test registration.
Its unique raw/timed-console paths are `task-5b-final-2-full.log` and
`task-5b-final-2-full-console.log`; default CTest LastTest.log is not relied on.
Final source/binary/cache identities are recorded in `task-5b-final-2-sha256.log`.
All 26 entries were verified after the full suite; raw verification is
`task-5b-final-2-sha256-check.log`. The complete earlier manifest is
retained as `task-5b-final-sha256-complete.log` and refers to the prior candidate.
The first partial manifest (`task-5b-final-sha256.log`) used an incorrect
Python 3.11 extension filename; the observed native extension is CPython 3.13,
and the complete manifest includes its actual file. This was a manifest
collection error, not a solver/build/test failure.

Final audited executable SHA256:
`4561b5d547086cbfae95875d285ad50d09708ec9694c37a39c5ed020edf1273f`.
Final native extension SHA256:
`5998f82afe48a005ae5ad065f7a7744cc80a6433d3c66b432156bffb0aaf3c35`.
Final core archive SHA256:
`7ae5d1a22477d49302565d5581a1a216046e93981623ea84eec38340a41113ae`.
Main CMake cache SHA256:
`2c3f5a72c82d0ba5d8b5d0cd8869d13d89ea8a21e0daa405e88242e1778e517c`.

Commands: `cmake -S . -B build/strict-v2-spin`; scoped
`cmake --build build/strict-v2-spin --target <targets> -j3`; focused
`ctest --test-dir build/strict-v2-spin -R <exact-regex> -V --output-log <unique>`;
`/usr/bin/time -v cmake --build build/strict-v2-spin -j3`;
`/usr/bin/time -v ctest --test-dir build/strict-v2-spin -j1 --output-on-failure
--output-log <unique-final-log>`; `git diff --check`; `sha256sum -c <manifest>`.
All implementation edits and configured build directories stayed within the
isolated tree. Shared system toolchain/dependency paths were only consumed.
The controller retains before/after evidence that dirty main was preserved.
Captured SDD failure logs and historical experiments were not deleted.
CMake's old-libigl compatibility warning
remains unchanged. No new dependency upgrade or warning suppression occurred.

## Limits and handoff

Self-review checked the entire scoped diff and new modules: multiplier signs,
near-face complementarity, independent residual validation, sparse LP row
assembly, actual boundary release, callback dimensions, current-lambda pairs,
membership reset before metric use, restoration evidence, final-state hashing,
default identity, Python/native guards, and independent recomputation tests.
Fresh independent review remains the controller's responsibility after handoff.

This is first-order evidence only. It does not prove a local minimum,
second-order conditions, global SLQP globalization/subproblem correctness,
nonlinear infeasibility, a high-resolution original-surface result, or a
physical certificate. Full sparse source geometry/cage work, original assets,
author trajectories and the remainder of Task5 remain separate obligations.
Generic trace hashes do not attest arbitrary callbacks. Numerical caps and
ill-scaled or inconsistent problems can return UNKNOWN or stalled. Existing
Faster result bundles become stale under their CMake input binding; they are
preserved, and the controller owns any unchanged-spec refresh after review.

## Fix round 1: reject reciprocal-curvature overflow before mutation

Fix base: `c624c845578ddb6c738caf1c9a4fca988ae0f843`. The independent review's
important finding was verified against `SlqpHistory.cpp` and reproduced in
the real history implementation. With fixed one-dimensional s=y=1e-160,
curvature is positive and finite, approximately 1e-320, while its reciprocal
overflows. The prior code accepted the pair and could evict valid history
before storing infinity. This fix computes and checks the reciprocal before
any deque eviction or insertion, then stores the already-validated value.
No threshold, objective, formulation, API or normal finite-curvature behavior
is changed. The deferred minor findings remain outside this fix's scope.

The regression prepopulates capacity-one history with s=2, y=4, reciprocal
.125. It then attempts the fixed tiny pair and checks rejection, all three
deque sizes and original values, plus the unchanged finite inverse action
H*2=1. An explicit fixture precondition requires finite positive curvature
and a nonfinite reciprocal, so flush-to-zero/underflow cannot make the test
pass through the older zero-curvature guard. The precondition was added
following the controller's clarification, after the initial RED observation.

`task-5b-fix1-red.log` and its console show the actual failure before the fix:
`subnormal_curvature=9.9998886718268301e-321 accepted=1 retained_inverse=inf
retained_metric=-nan`. The test fails with “nonfinite reciprocal curvature
accepted”; its RED build is `task-5b-fix1-red-build.log`.
After the fix, `task-5b-fix1-green.log` passes the four focused tests. The final
test including the explicit precondition was rebuilt in
`task-5b-fix1-final-test-build.log`; `task-5b-fix1-final-focused.log` and its
console again pass 4/4 (including nine Python unittest methods). Actual GREEN:
`subnormal_curvature=9.9998886718268301e-321 accepted=0 retained_inverse=0.125
retained_metric=1`.

The all-target `task-5b-fix1-build.log` records -j3, exit zero, 26.24 seconds,
peak RSS 1,757,380 KiB. The last edit after this build only adds the test's
arithmetic precondition; that executable was rebuilt before the final focused
and full runs. Normal legacy SLQP still has hash
`ba7371831d27b970aeab54b8c1c14f8b411023e16bd06c7fc0c8a9d27f9d9aae`.
The retained BASE comparison driver was relinked against the fixed core;
`task-5b-fix1-default-givens.log` confirms equal BASE/current default Givens
hashes `13fcb7597f2bc6a0275921674dd5f109f63abeb0013487e5216bfba14aa99eb5` and
exactly equal final beta/phi. Original fixed fixtures and tolerances remain
unchanged.

Final fix1 full suite passed 190/190 tests on unchanged fixed binaries, with
190 actual Passed records, zero failures, exit zero, 314.01 seconds elapsed,
and peak RSS 4,253,644 KiB. Unique logs: `task-5b-fix1-full.log` and
`task-5b-fix1-full-console.log`. The registered 2014-ON/reachability-OFF
dependency contract built and executed the actual audited solver, including
this boundary regression, successfully as part of the full run.
Previous full runs and raw failure logs were retained under their old names.
The refreshed 26-entry source/binary/cache manifest is
`task-5b-fix1-sha256.log`; all 26 entries were rechecked after completion and
match, as recorded in `task-5b-fix1-sha256-check.log`.

`task-5b-fix1-loaded-dso.log` records an actual import using the configured
CTest Python/environment, the module's resolved path, SHA256, and its live
`/proc/<python-pid>/maps` entries. Interpreter:
`/home/lin/anaconda3/bin/python3`. Actual DSO:
`build/strict-v2-spin/python/_make_it_stand.cpython-313-x86_64-linux-gnu.so`,
SHA256 `54d47b3f8ee834d4ceaac145d7491a692f07ae9038687a38103166a0d980519c`.
The import used the source/build/test PYTHONPATH and configured GCC 13
libstdc++ LD_PRELOAD from CTestTestfile.cmake, not a guessed installed extension.

Commands retain the original isolated build configuration and resource limits:
`cmake --build build/strict-v2-spin --target spin-it-2014-audited -j3`;
`/usr/bin/time -v cmake --build build/strict-v2-spin -j3`;
`ctest --test-dir build/strict-v2-spin -R '^spin-it-2014-(audited|slqp|givens|python)$'
-V --output-log <unique-fix1-log>`;
`/usr/bin/time -v ctest --test-dir build/strict-v2-spin -j1 --output-on-failure
--output-log <unique-fix1-full-log>`; `sha256sum -c <fix1-manifest>`.
The loaded-DSO record imports `_native_module()`, prints `__file__` resolved
through realpath, and invokes `sha256sum` and `rg` on the live process maps.
All code edits are limited to the reciprocal guard and its actual-history
regression. No dirty-main, historical-evidence, dependency, asset or unrelated
trace cleanup was performed. First-order-only and incomplete whole-method
limitations above remain; the controller owns the scoped re-review.
