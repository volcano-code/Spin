# Task5v — actual-source guarded SLQP integration

Status: DONE_WITH_CONCERNS. The bounded guarded integration and scoped worker
verification are complete. Native allocator internals remain externally guarded
and the explicitly open Task5u negative-fixture gaps are not claimed closed.
This is synthetic engineering integration, feasibility UNKNOWN, with no paper,
physics or seal grant. Final controller comparisons/review are separate work.

## Scope and implementation

The controller approved the typed API and simultaneous resource proposal before
production. `SourceCageGuardedSlqpAudit` is privately constructed, exposes const
facts and returns by move. The entry point consumes the immutable actual paper
functional, full state, original solver settings and actual geometry settings.
The original source adapter was factored into one file-private common runner,
with an internal guard pointer. The unguarded route passes null. There is no
duplicated generic engine, cloned source adapter, or user callback admission.

Every initial, objective/fallback, restoration and final source evaluation is
gated by the actual `SourceCageGeometryChecker`. A completed passing checker
audit is required before derivatives run. Rejection/UNKNOWN throws through the
existing search machinery, which shrinks or exhausts while retaining the last
accepted state, moments, objective and curvature history. Initial failure
retains the actual attempted state and stages without manufacturing a numerical
result, KKT witness or final state. Final acceptance requires a new actual
checker invocation and recomputed actual material/derivative output.

Resource failure is permanently latched in the common source outcome. Later
generic callback attempts retain a disposition but invoke neither checker nor
source. State-hash and record failures receive the same treatment; incomplete
checker publication is never a passing audit. Failure outside recoverable
publication/storage propagates without a returned pass. A numerical termination
is preserved separately, including convergence before a failed final operation.
The guarded root's resource/final status is authoritative for geometry execution.
The nested legacy geometry status stays NOT_EVALUATED; feasibility stays UNKNOWN.

The flat stage tape retains the original checker records, exact rational
clearance strings and witnesses, primitive/beta IDs, complete failure/evidence
hashes, and their order. Each fixed-size attempt indexes the existing raw source
callback-state tape and its stage span. The checker scientific/input/source
output hashes and full derivative output hashes are retained separately. The
scientific digest binds source/config/frame, all settings/budgets, raw states,
guard records, solver trial/acceptance/history evidence, failures and final
material output, without timestamps.

Counts distinguish actual checker invocations, geometry source-value calls,
full derivative source calls, and actual final checker calls. Preflight rejection
counts a checker call but zero source-value evaluations. A checker exception
before stage publication marks the source-value count incomplete; no total
source count is claimed in that case. A skipped final attempt after a resource
latch does not increment finalRecheckCount.

Owned implementation paths changed:

- `src/spin2014/SourceCageSlqp.cpp`
- `include/spin2014/SourceCageGuardedSlqp.h` (new)
- `tests/spin_it_2014_source_cage_guarded_slqp_regression.cpp` (new)
- `docs/spin2014_source_cage_guarded_slqp.md` (new)
- Eight gated registration lines in `CMakeLists.txt`

`include/spin2014/SourceCageSlqp.h` remains byte-identical. Its public struct
layout, legacy accounting identity and old expected hashes were not altered.
No generic solver, exact checker, functional, geometry backend, foreign Faster/
Canonical/Plane/Quadric source, dependency, Python API or viewer was edited.
No extra implementation file, worktree, installation, push, merge or tag.

## Resource envelope

Finite IMPLEMENTATION_DEFAULT limits are 1024 actual checker calls, 65536
retained stage records and 2 GiB guarded known buffers. The original adapter
planner is unchanged. Guarded preflight validates actual checker settings and
checked count/products before capacity allocation. Actual primitive count P
gives at most `3+3P+P(P-1)/2` stages per guard. The complete declared guard-call
budget times that count must fit the retained-stage cap; no budget is silently
lowered. Runtime guard exhaustion remains an attempted resource outcome.

The added envelope includes preallocated fixed attempts and flat stage records;
twice bounded exact-string payload per primitive per declared guard call;
twice the checker buffer cap plus diagnostic/metadata reserve for one current
checker and its returned audit; the unchanged adapter/source/LP/history/hash
workspace simultaneously; and fixed guarded result/settings/hash overhead.
There are no per-trial retained geometry matrices or Nef objects. Return by
move prevents deep return-copy growth. Caller copies of extracted facts are
outside the solve's own allocation envelope.

Opaque CGAL/GMP/HiGHS allocations and allocator metadata still require external
process guards; this is not an exact RSS guarantee. Every heavy run used the
assigned `build/task-5c-no-faster-20260909`, 2014ON/exactON/FasterOFF,
LIBIGL_SKIP_DOWNLOAD=ON, OMP/OPENBLAS/MKL1, `ulimit -v 25165824`, timeout3600,
at most3 compilation jobs and sequential tests.

## TDD and focused evidence

The implementer template, TDD skill and writing-good-tests instructions were
read before production. The registered test preceded the new header/production
implementation and used the existing missing-component behavioral pattern.

Initial RED, `task-5v-red-test.log`, command:

```
ctest --test-dir build/task-5c-no-faster-20260909 \
  -R '^spin-it-2014-source-cage-guarded-slqp$' --output-on-failure
FAIL: actual-source geometry gate cannot reject an invalid initial state before SLQP evaluation
0% tests passed, 1 tests failed
```

After implementing the gate, `task-5v-green-test-2.log` is 1/1 PASS in .21s.
The initial file named `task-5v-green-test.log` is intentionally preserved as a
failure: the missing-header `__has_include` branch was still in the existing
test object because adding a previously absent header did not invalidate its
dependency list. Recompiling the changed actual test body exercised the new API
and produced the passing result; this stale build observation is not hidden.

The real one-shot allocator test overrides only test-binary operator new,
injecting a failure at each of48 declared 65-byte SHA-allocation ordinals.
Its RED `task-5v-allocation-red.log` caught:

```
FAIL: every attempted callback has a guard disposition, including bookkeeping resource failure
```

The guard disposition is now preallocated/started before state hashing. The
GREEN sweep `task-5v-allocation-green.log` reports48 injected failures,47
returned resource audits and1 propagated allocation failure; no returned
success. The returned failures pass the raw state/guard links and resource-latch
assertions. A separate moved-from audit RED and GREEN are preserved in
`task-5v-move-red.log` and `task-5v-move-green.log`.

Comprehensive focused test `task-5v-comprehensive-test.log`:1/1 PASS,66.19s.
Its real-source assertions cover:

- All five objectives and all112/113 beta/phi/full-handle variables, each run
  twice with identical complete scientific digest and final state.
- Every attempted state/source hash, guard stage span and ordered solver trial;
  all accepted trace states match a passing guard/full evaluation.
- Actual solver-generated geometry rejection followed by reduced accepted
  objective step, checking the entire raw state against direction times step.
- Actual restoration trial states against recorded direction and step, with
  rejected candidates followed by a guarded accepted restoration step.
- Real exhaustion preserving last accepted state, moments, objective and
  empty curvature history; no manufactured convergence or infeasibility.
- Initial rejected geometry, initial source-value UNKNOWN, checker preflight
  and native Nef resource failure, derivative resource failure after a pass.
- Permanent guard/full-source budget latches, later generic attempts without
  further geometry/source work, and truthful separate source-work counts.
- Exact/one-below/zero/overflow guard-call, retained-stage and total-byte caps;
  failed final guard or final material recomputation never yields final pass.
- Const/private audit construction constraints and moved-from audit access.

## Preserved synthetic calibration

Both calibration protocols were declared and approved before running. The
fixture uses actual full26-handle BBW, source depth3/adaptive levels2..3,
outer cube[.3125,.6875]^3 and design cube[.34,.66]^3. Identity handles,
beta[.17,.24,.31,.38,.45,.52,.59,.66], optional phi.23, frame origin(.48,.51,.3),
weights inertia/center/fill/shape/displacement .8/.6/.03/.04/.05 are fixed.
There is no mock predicate, forged passing audit or injected solver direction.

`task-5v-calibration-25.log` retains all25 cases: five objectives times clearance
{0,.005,.02,.031,.04}, one iteration and8 backtracks. Selected mechanism cases:

- TopMassAndLowering, clearance0: steps1/.5 UNKNOWN due to actual source mass;
  steps.25/.125 are geometry rejection at PrimitiveSubset, then.0625 passes
  and is accepted. UNKNOWN is never counted as a geometric rejection.
- YoYo, clearance.04: actual rejected restoration trials shrink to accepted
  step.015625. Full raw restoration direction/step linkage is checked.

`task-5v-calibration-exhaustion-4.log` retains all4 separately approved follow-up
cases: TopMassAndLowering, clearance.04, maximum backtracks{1,2,3,4}. All4
reported RESTORATION_STALLED and zero final-state displacement. The selected
four-backtrack case includes actual geometry-rejected objective and restoration
trials, retains last valid state/moments/history, and freshly checks that final
accepted state. No open-ended search, tolerance relaxation or geometry repair.

Task5u's actual-source exact-contact and positive-volume pair-overlap negative
coverage gaps remain open; this integration does not claim to close them.
Conservative domain/outer clearance applies to all cells regardless of beta,
not a continuous BBW Jacobian or cavity-web thickness. None is an original-paper
experiment, feasible witness, seal or physical result.

## Final verification, hashes and commits

The complete related registered group passed **43/43 in347.63s** after scoped
relink of all existing2014 executable targets and `_make_it_stand` at-j3:

```
ctest --test-dir build/task-5c-no-faster-20260909 --output-on-failure -j1 \
  -R '^(spin-it-2014-|spin-2014-|source-spin-material-python$)'
```

`task-5v-final-group.log` and preserved full `task-5v-final-raw.log` are the
results. `task-5v-final-relink.log` records33 scoped relinks. No skipped,
disabled, failed or NotRun case was counted as a pass. Source manifests
`task-5v-source-begin.sha256` / `task-5v-source-end.sha256` matched byte-for-byte,
with manifest SHA159ceac6df27b853365a685f30ee6e0d2603248836dc2865f80a7f3fe8850d9e.

After that frozen run, the controller explicitly approved one narrow header
hardening with a focused rerun instead of repeating all43 tests. A facts
reference obtained before moving the audit could still observe stale scalar
availability beside transferred vectors. New behavioral RED is preserved in
`task-5v-retained-reference-red.log`:

```
FAIL: pre-move const reference cannot retain stale passing scalar or material availability
```

The final move constructor now resets the moved-from root final status to
Unavailable and both numerical/final material availability flags to false,
in addition to rejecting future moved-from accessor calls. The new assertion
checks that retained pre-move reference. No numerical algorithm, checker,
resource formula or source evaluation behavior changed. The final header and
test were frozen, every scoped2014/Python target relinked, and the covering
registered group passed **5/5 in95.00s**:

```
ctest --test-dir build/task-5c-no-faster-20260909 --output-on-failure -j1 \
  -R '^(spin-it-2014-source-cage-guarded-slqp.*|source-spin-material-python|spin-it-2014-python)$'
```

The final five cases are guarded64.59s, allocation2.21s, move3.96s,
source-material Python24.12s, and2014 Python.11s. They are preserved in
`task-5v-hardened-tests.log` and `task-5v-hardened-raw.log`, with relink evidence
in `task-5v-hardened-relink.log`. The full43 group is explicitly pre-hardening
evidence; there was no second43-test run. The final core and module bytes are
identical to the broad-run bytes, while the guarded test binary changed.

Current-source Python evidence is from the final narrow run, not reused from
before the header change. `task-5v-hardened-raw.log` lines122/124 report10tests
OK in24.036s for source-spin-material-python; lines151/153 report9tests OK in
.029s for2014 Python. Final source manifests
`task-5v-hardened-source-begin.sha256` / `task-5v-hardened-source-end.sha256`
match byte-for-byte, with SHA:

```
288e2e899f723442059d29e3847ed9f778501ac79e120f328ca8a555f24354d7
```

The actual fresh exact-OFF child created by the broad registered dependency
test is:

```
build/task-5c-no-faster-20260909/spin-2014-dependency-contract-runs/run-ef5629b3f01b7249/enabled
```

The retained parent is named at line1130 of the broad raw log. Its cache
explicitly records2014ON/exactOFF/FasterOFF/LIBIGL_SKIP_DOWNLOAD=ON. Actual
graph inspection in `task-5v-exact-off-graph.log` found no guarded/source-SLQP/
geometry checker source, header or target name. Its graph SHA is
b3cb021ac901a4d87be256fe6aa1c14a0e91153e36d3c0110e30410b0687c42f.
The assigned active graph also excludes foreign Canonical/FullRank/Plane/Quadric
implementation paths. This graph inspection precedes the final move-only
header-body change; no new child was generated afterward, and no CMake/source
graph membership changed. Both registered dependency tests passed.

Final core SHA:
`b436244dc23e8fde6da4fd1c29b9af33dcaaf3b525005cce122c31b9eb14dd74`.
Final module SHA:
`5b1092d3baf5dafca9cfe8bd27df1469e583f1b28bef750aaa6302b4c0b34453`.
Final guarded native binary SHA:
`45e7566664123915585bc104e2ea1ac2e740167c2ee67ec2a1f6b7f920d917f8`.
The full source/checker/core/module/binary/log/probe manifest is
`task-5v-final-artifacts.sha256`, manifest SHA
`d286edf9ae2f7f1351cc7ee1139f9063daf4fc7ec05e19de4d877960849656dd`.

Implementation commit: `c52fdbd feat(spin2014): gate actual source SLQP trials with geometry`.
This report is a subsequent separate scoped commit. Controller helpers and
`task-5v-controller-*` artifacts are not staged or owned by this implementation.

## Self-review and independent work boundary

Read the full scoped diff and checked whitespace/staged paths. Verified the
legacy header remains unchanged, source adapter bookkeeping and generic engine
are factored/reused, original geometry records are retained, all later work
stops under the resource latch, final source evaluation is fresh, and every
accepted native state is linked to a passing actual checker. Fixed both the
allocation-before-guard-record gap and moved-from/retained-reference stale
availability before handoff. All touched implementation paths are in scope;
foreign dirty files remain untouched.

The controller independently owns raw exports, scalar KKT/link verification,
rational clearance/material reintegration and corrupted/dropped/mismatched
evidence checks. Its preliminary8-case result was communicated as passing,
but its final10-case two-run export and Task5t archived legacy byte comparisons
are pending in separate controller evidence at this report's commit. They are
not presented as worker-run or already-final evidence here. Existing registered
legacy source-SLQP tests did pass in the43-case group. The controller explicitly
requested scoped commits without waiting for its final exporters so a fresh
reviewer can examine the frozen implementation in parallel with those checks.
