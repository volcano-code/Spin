# Task5w — atomic actual-source adaptive partition rebind

Base: 60b31e0f3d8568b045a0428b94c54e973c3228a2. Continuation of approved Task5,
not a new paper profile. Read this brief first; parent-task extract is background.
Work only in /data/xvxinyao/moment-guided-spin-v1/code/make-it-balance-paper-sim.

## Bounded delivery

Add a native immutable rebind transaction between SourceAdaptiveDomain and the
existing actual full-BBW guarded source solver. This is NOT the complete adaptive
optimizer, binary material validation, a new sealed artifact, or paper acceptance.

Inputs: an actual SourceCagePaperFunctional, its full state, fractional epsilon
(default .05, explicit IMPLEMENTATION_DEFAULT), and actual geometry checker settings.
Use the existing splitFractionalLeaves followed by mergeBinarySiblings on a LOCAL
copy of its binding domain and beta. Preserve the full phi and every handle's four
parameters exactly, taking offsets from old/new functionals (including static mode).
Do not accept an externally supplied partition, replacement BBW, passing bool or
stale audit as authorization. Re-check the input's actual represented geometry.

Reconstruct SourceCageBinding -> SourceConformingBoundary -> SourceCageMomentFunctional
-> SourceCagePaperFunctional with the actual retained source/BBW and unchanged binding,
boundary, moment, paper configs, density, frame. All handles remain present. Recompute
new points, embeddings, primitive ownership, moment/constraint derivatives and fill
Laplacian as needed, not a resized copy of old derived data. Freshly check the resulting
actual geometry, including no-change transactions; check actual material values and
retain their hashes and ten moments. Publish the new immutable functional and remapped
state ONLY after all applicable checks succeed. Old source/state remain bit-identical
on success or any failure. No caller-accessible mutable fields that can turn a failed
transaction into an authorized checkpoint. Prefer const shared immutable storage to
avoid moved-from stale passing references. Numerical/history data are never inherited.

Failure dispositions: resource, rejected geometry, unknown; invalid API settings may
throw. Do not manufacture a new usable checkpoint after failure. Resource exceptions
may propagate if diagnostic allocation itself fails. Preserve enough staged diagnostics
to distinguish input check, split, merge, rebuild, output check/material evaluation.
Use upstream caps (including SourceAdaptiveDomain leaf budget) plus checked bounded
owned transaction evidence/state storage; document caps are not native RSS guarantee.
Caller must continue using 24 GiB / 60 minute process guards. No silent budget raises.

Bind deterministic scientific hashes to old/new source/partition/state/config,
transition maps and checker/material outputs. Timing is not in scientific hash.
Expose old/new beta/leaf identity and split/merge facts needed for independent checking.
No previous numerical convergence/KKT/geometry conclusion carries into a new solve.
Demonstrate the resulting checkpoint feeding solveSourceCageGuardedSlqp, whose initial
history is empty and actual geometry is checked again. No need to add an automatic
multi-pass loop this increment.

Important mathematical scope: SourceAdaptiveDomain's removedBefore/After are REST
geometric moments. In general non-affine BBW, refinement can change the represented
polyhedron; report before/after deformed moments and their difference, do NOT require
them equal or label rest conservation as deformed conservation. Preserve exact old
BBW field; rebuild samples at new points. A changed discretization requires new
geometry verification and numerical solve, even when inherited beta is identical.

## Tests (TDD, actual source fixtures, no mocked acceptance)

Write an executable red test that fails because the rebind feature is absent before
production implementation. Then cover:
- Fractional split with independently checked ownership/inherited beta, full phi and
  all handle params preserved, changed partition/material bindings; old objects intact.
- Equal binary sibling merge and unequal near-bound values not merged; no-change path.
- Static no-angle and non-static angle layouts; immutable config/frame/density retention.
- Rest or globally affine ten-moment conservation checked independently; non-affine
  output freshly reintegrated rather than asserting arbitrary conservation.
- Actual fresh input/output geometry rejection and resource exhaustion cannot publish.
  Use attainable source fixtures; if a required negative cannot be constructed, report
  it explicitly rather than forge flags. Invalid dimensions, NaN and epsilon rejected.
- Actual new checkpoint -> guarded solve, new state dimension, history reset, hashes
  and final recheck; numerical timeout/stall is not convergence/feasibility.
- Identical runs have identical scientific hashes; changed input/config changes binding.
- Test output exposing genuine before/after geometric data for controller verification
  is welcome, but avoid building a general exporter or coupling production to tests.

Allowed production files: new include/spin2014/SourceCageAdaptiveRebind.h,
src/spin2014/SourceCageAdaptiveRebind.cpp, scoped CMake registration, new tests and
docs/spin2014_source_cage_adaptive_rebind.md. If existing production must change, ask
controller first. No Python binding required this increment. Do not touch foreign
include/spin or src/spin or their tests, nor experiments/m1_* work.

## Engineering and report

No subagents; independent review is controller-owned. apply_patch edits only.
Build only build/task-5c-no-faster-20260909 (2014/exact ON, Faster OFF, skip downloads).
At most 3 compile jobs, all heavy runs sequential, single BLAS/OpenMP thread,
ulimit -v 25165824 and timeout 3600. Do not modify dependencies or main worktree.
Do not delete old logs. Commit only explicit task-owned files as authorized identity
volcano <xuxinyao307@gmail.com>; do not stage generic progress.md or foreign paths.

Report to task-5w-report.md in this directory: implementation, exact API, tests and
commands/outputs, RED/GREEN evidence, file list, self-review, remaining limitations.
Return only status, commits, short test summary, concerns, report path.
Controller concurrently owns Build Graph evidence and independent validation artifacts.

Old deferred findings: SourceCageSlqp.cpp:143 / SourceCageGeometry.cpp:172 readability;
two actual-source negative geometry fixtures (exact contact and primitive positive
overlap) remain open. Do not change their status without actual evidence.
