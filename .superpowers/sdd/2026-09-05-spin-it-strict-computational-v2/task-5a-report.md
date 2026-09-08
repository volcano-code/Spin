# Task5a — explicit paper Givens moments and joint fill/angle integration

BASE `3114cda0e36feef538732c55d8ef32620813a861`, branch
`spin-it-strict-computational-v2`, workspace
`/data/xvxinyao/make-it-balance-paper-sim`. This is a bounded Task5 increment.
The full reproduction remains **PARTIAL_BLOCKED**. No original model,
physical result, full paper SLQP restoration, cage, surface-source octree,
level-9 case, or dual digital-cut completion is asserted. No subagents,
push, merge, tag, release, dirty-main edits, or historical experiment edits.

## Implemented contract

`Spin2014GivensProblem` implements the explicit equations from authenticated
`build/strict-paper-sources/spin2014.txt` lines 208–275. The project order is
`(M,Lx,Ly,Lz,Qxx,Qyy,Qzz,Qyz,Qxz,Qxy,phi)`. A focused second-order scalar
uses analytic sum/product/unary chain rules for all objective gradients and
Hessians and ordered constraint Jacobians/Hessians. It neither invokes
finite differences nor copies/refactors the Faster derivative evaluator.

Yo-yo uses the axis-origin inertia; all top variants subtract `Lz²/M` only
from transverse diagonal entries. The raw constraints are ordered
`Lx,Ly,Lz,Qxz,Qyz,g` for yo-yo, `Lx,Ly,Qxz,Qyz,g` for tops, and `Lx,Ly` for
static balance. With `R=[[c,-s],[s,c]]`, the actual rotated offdiagonal is
`-g`. Independent inertia/center weights and all five objectives follow the
brief. Nonpositive mass, nonpositive used axial inertia, nonfinite state,
or nonfinite analytic outputs reject. Repeated transverse eigenvalues
remain finite and angle invariant, with no unique-axis assertion.

`AdaptiveFillFunctional::evaluatePaper(beta,phi)` shares geometry/moment
assembly and the existing Laplacian with the legacy evaluator. Its chain
uses `-momentColumns` and a separate phi column. Only beta receives the
Laplacian term. `betaCount()` and `paperHasAngle()` distinguish fill leaves
from the joint coordinate. Direct boundary/interior density weighting
continues to work.

The existing SLQP loop is now one internal callback implementation with a
beta bound count. `solveSlqpLbfgs` is the legacy compatibility wrapper;
`solveSlqpLbfgsPaper` solves `[beta,phi]`, returning beta without the angle.
All active-bound logic and bound step/clipping apply only to beta. Phi is
finite, unbounded, never clamped or wrapped, and participates directly in
accepted steps and L-BFGS history. There is no diagonalization postprocess.
Static solves have beta alone and require initial phi zero. The standalone
unified static moment evaluator permits finite phi with zero derivatives;
all static/legacy request interfaces require zero to avoid irrelevant
angle trajectories. New C++ aggregate fields and Python dataclass fields
are appended to preserve existing positional meanings; ABI compatibility
of extended C++ structs is not asserted.

`AdaptiveVoxelSpinRequest` defaults to `LegacyAngleEliminated`; opt-in
`PaperGivens` and `initialPhi` select the new path. Optimized phi is carried
between adaptive passes and retained on relaxation/resource failures.
Each explicit trace row carries the current angle; a row describes the
state before its recorded accepted step. Final phi is separately retained
and included in the scientific trace identity, including max-iteration
exits immediately after an accepted step. Request, adaptive trace, SLQP
trace, and output hashes use distinct explicit profile identities and
bind the actual angle. Legacy hash preimages remain unchanged.

Final binary constraints use the paper evaluator at final phi. The old
`beta >= 0.5` rule and `BINARY_VALIDATION_FAILED` are retained, without beta
repair or a changed objective. Nonconvergence is never infeasibility.

The native Python binding accepts `formulation` and `initial_phi`; the
Python facade exposes `formulation`, `initial_phi`, `final_phi`, and the
angle-bearing diagnostic trace. Native JSON adds fields only for the
explicit profile. `original_phi_trajectory` remains false: this is an
existing-engine integration, not recovered author iteration semantics.
The native and Python cage consumer rejects explicit requests/artifacts.
The dedicated dual-material product/export wrapper also rejects the new
profile because that wrapper is not validated in this increment; this
does not remove direct adaptive density weighting.

## Actual RED, numerical refinement, and GREEN evidence

Raw logs originated under `build/strict-v2-spin/task5a-*-20260908.log` and
are preserved. The selected tracked evidence copies beside this report
retain the same names and bytes. No failed attempt was overwritten.

1. Before production implementation, registered the new C++ test and wrote
   fixed five-objective rational expectations plus all-coordinate numerical
   checks. An initial build before reconfigure had no target (`red-build`),
   and the reconfigured build correctly lacked the new header (`red-compile`).
   A temporary feature-presence branch made the executable return the actual
   runtime RED `[FAIL] explicit Givens evaluator is unavailable`
   (`red-runtime-build`, `red-runtime`). The branch was removed once the
   evaluator existed. `math-first` preserves a stale build caused by that
   temporary `__has_include` dependency; the test was then rebuilt normally.
2. Before functional/solver/API integration edits, the Python behavior test
   failed its explicit-path availability assertion (`red-python`). The same
   test subsequently traverses the actual native solve, checks angle change
   and accepted steps, repeated output identity, and cage rejections.
3. `math-second` records the initial three-point centered-difference gate
   failure at h=1e-3. No production formula or fixed moment fixture changed
   to address it. After observing this truncation effect, the checker was
   explicitly changed to fourth-order centered differences at the same
   h={1e-3,3e-4,1e-4}, retaining the 1e-6 acceptance threshold. Controller
   accepted this as numerical-check refinement, not tolerance relaxation.
   All eleven objective gradient/Hessian coordinates and all constraint
   Jacobian/Hessian entries are checked. Nonzero derivatives use relative
   error; entries with magnitude <=1e-7 use absolute error. Maximum accepted
   error is about 1.213e-10. The fourth-order values are at roundoff level,
   so small nonmonotonicity below a documented 1e-8 trend floor is allowed;
   the 1e-6 accuracy gate is unchanged. Original second-order objective
   Hessian diagnostics are retained for all three steps: approximately
   2.27567e-6, 2.04811e-7, 2.27564e-8 for the rotational cases, showing the
   expected second-order decrease. Static values are approximately
   3.75e-7, 3.375e-8, 3.752e-9.
4. The rational fixture is independently supplied by the controller:
   M=2, L=(1/3,-1/4,1/2), Qdiag=(3,2,4),
   (Qyz,Qxz,Qxy)=(1/7,-1/5,3/4), phi=atan2(4,3), gammaI=2, gammaC=3.
   Objective literals are 107474/15625, 1655209/250000, 425521/62500,
   1842709/250000, and 3/16. The test also checks independent angle
   derivatives, raw constraint order/sign, density covariance, periodicity,
   negative and >1-rad angles, repeated eigenvalues, invalid state, legacy
   off-manifold difference, and separate feasible diagonalization agreement.
5. `first-joint-tests` records a fixed synthetic 3x3x3 grid whose low-level
   joint solve converged to residual 1.69808e-14. Its later adaptive export
   rejected the nondivisible source dimensions. The corrected positive
   integration fixture is 4x4x4, keeping bounds (-1.5,-3,-1.5) to
   (1.5,3,1.5), with the eight center cells carveable and all 64 material
   cells present. This correction satisfies the export precondition; it did
   not improve optimizer feasibility. The original 3x3x3 dimensions and
   center-cell carveability remain an expected export-rejection test.
6. The valid joint fixture uses top inertia only, initial beta=.2,
   initial phi=1.7, Laplacian weight .03, and equality tolerance 1e-9.
   It accepts four unit steps and reaches phi approximately pi/2, with
   equality residual 1.42131e-11. The actual rows and beta values are printed,
   not inferred from populated result fields. A separate nonconstant beta
   field checks beta/phi chain rules and nonzero beta-only regularization.
   A two-pass adaptive variant checks exact angle carry. Repeated solve
   hashes/results match; formulation/initial angle alter input identities.
   A 1e-14 binary tolerance deliberately returns binary validation failure;
   max-iteration and resource exits retain angle provenance and no success
   certification. Synthetic fixtures are integration tests only.
7. During self-review, `dual-red-tests` shows the dedicated dual wrapper
   accepted the explicit profile before its guard. The guarded wrapper
   passes in `dual-green-tests`; legacy wrapper tests remain registered.

## Compatibility and scope checks

Before production edits, `baseline-valid` ran the actual BASE SLQP and
end-to-end executables. Their scientific hashes are now literal regression
expectations: SLQP `ba7371831d27b970aeab54b8c1c14f8b411023e16bd06c7fc0c8a9d27f9d9aae`
and end-to-end output
`d1721e8f1564bedf596f9eb88469a1b28bf9040930de3dbfbff9a6ba56275ad2`.
The earlier empty regex run is also retained as `baseline`; it is not counted
as validation. Existing default Python diagnostics retain their schema and
semantics. No legacy JSON fields are added; source review verifies that the
equivalent split string writes preserve existing legacy bytes. No independent
whole-JSON BASE fingerprint was captured in this increment, so that stronger
byte-comparison claim is not made. `git diff BASE -- experiments` is empty.

## Verification and artifact identities

Commands use the primary `build/strict-v2-spin` build. CMake feature flags
are 2014/Faster/exact/Python ON, viewer/PyBullet OFF. Build maximum was -j3;
heavy tests run serially. The all-target build exited zero in 33.96 seconds,
maximum RSS 1,746,984 KiB (`final-all-build`). The full configured suite uses
an explicit unique `ctest --output-log`; `LastTest.log` is not cited as full
raw evidence. Full final run: **189/189 passed**, zero failures or skipped
tests, 260.62 seconds, maximum RSS 4,253,708 KiB. The final raw console is
`task5a-final-full-raw-20260908.log`; the accompanying timed console is
`task5a-final-full-console-20260908.log`. This was the one complete configured
suite run after the final implementation. Its 2014 label contains 18 tests;
the Python API CTest executes 8 unittest methods. No tests were unregistered.
All 23 final source/binary/cache manifest entries were rechecked with
`sha256sum -c` and matched after the suite. Tracked evidence copies are
byte-identical to the original build logs (`cmp`, exit zero).

```sh
cmake -S . -B build/strict-v2-spin
cmake --build build/strict-v2-spin --target spin-it-2014-givens-regression -j3
ctest --test-dir build/strict-v2-spin -R '^spin-it-2014-givens$' -V --output-log <unique-log>
ctest --test-dir build/strict-v2-spin -R 'spin-it-2014-python' -V --output-log <unique-log>
/usr/bin/time -v cmake --build build/strict-v2-spin -j3
/usr/bin/time -v ctest --test-dir build/strict-v2-spin -j1 --output-on-failure --output-log <unique-full-log>
git diff --check
```

`final-sha256` lists every changed production/test source and final binaries.
Givens executable SHA-256:
`e7d4af01c365dbaeb09b6ebc90f48d9c41752774e43eda47586fe0d216ffe632`.
Native extension:
`4052d643674d5538851f4aa2f711b2b960ecdf1a26c5b6242265aa80aea0317a`.
Core archive:
`b65a3357d4306ec96f597b9128578ea539ba7b98c316e5d60905a83fce509b5f`.
CMake cache:
`2c3f5a72c82d0ba5d8b5d0cd8869d13d89ea8a21e0daa405e88242e1778e517c`.

Existing dependency warning: configure prints libigl's old-CMake
compatibility deprecation. No dependency upgrade or warning suppression
was made. Source/test whitespace checks are clean; raw console whitespace
is retained exactly in evidence.

Self-review covered the new evaluator and complete scoped diff: raw versus
full-COM formulas, rotation sign, all derivatives, beta bound membership in
both normal and fallback line searches, absence of angle wrapping,
append-only public fields, final accepted-state hash, adaptive angle carry,
binary selected evaluator, failure provenance, and unsupported consumers.
The controller owns independent review after this handoff. This does not
close future full KKT/restoration or active-set history obligations.
