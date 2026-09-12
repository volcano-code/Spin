# Task5v — actual-source SLQP geometric trial admission

Continue approved Task5/M2–M3 integration after checker Task5u. BASE c8813a5.
This is the already approved reject/shrink/rollback/final-recheck step; not a
new solver, a paper experiment, or full binary/sealing completion.

## Scope and protection

Canonical worktree /data/xvxinyao/moment-guided-spin-v1/code/make-it-balance-paper-sim,
branch spin-it-strict-computational-v2. Protected main is read-only. Foreign
Faster/Canonical/Plane/Quadric pending sources and experiment directories are
excluded. No dependency installation/download/change, merge/push/tag/deletion.
Use apply_patch. Identity already authorized volcano <xuxinyao307@gmail.com>.

Owned paths: include/spin2014/SourceCageSlqp.h,
src/spin2014/SourceCageSlqp.cpp, new include/spin2014/SourceCageGuardedSlqp.h,
new tests/spin_it_2014_source_cage_guarded_slqp_regression.cpp,
new docs/spin2014_source_cage_guarded_slqp.md, minimal CMakeLists.txt registration.
If a small private implementation header/source is necessary, propose it first.
Do not clone the generic SLQP engine or duplicate the entire existing adapter.
Do not modify existing exact checker, paper functional, geometry backend or
generic solver without a concrete controller-approved reason. No Python API or
viewer/new artifact work here. Report task-5v-report.md beside this brief,
separate scoped commit from implementation. Controller owns independent helpers.

## Required contract

Add a clearly named native guarded entry point consuming the immutable actual
SourceCagePaperFunctional, full initial state, solver settings and actual
SourceCageGeometrySettings. Reuse the current audited KKT/restoration/L-BFGS
engine and existing source adapter bookkeeping through factored common code.
Keep unguarded API behavior/scientific identity compatible and explicitly
NOT_EVALUATED; never silently promote old results. No supplied user callback,
mesh, evaluation or passed boolean can authorize geometry in the guarded route.

Propose the minimal typed API and simultaneous resource accounting BEFORE
production. Prefer a privately constructed guarded audit with const facts,
retaining numerical diagnostics separately from its trustworthy guard execution
history. Existing mutable numerical results alone must not forge admission.

1. Check every actual initial/trial state with SourceCageGeometryChecker before
   it may yield an accepted objective/gradient evaluation. Only CHECKS_PASSED
   permits the source's full derivative evaluation for that state. No reuse of
   pass facts for a changed state. All beta, phi and full handles remain.
2. Geometry rejection or incomplete UNKNOWN trial is invalid for line search,
   preserving exact cause and ordered stage evidence. Existing backtracking and
   restoration paths must shrink/reject it; failed candidates never replace
   current state, moments, objective, or L-BFGS history. Exhaustion leaves last
   valid accepted state, with actual failure termination, not convergence.
3. Reject invalid initial geometry before an accepted numerical iteration;
   retain the actual attempted state/audit without manufactured KKT/final state.
4. Resource failure latches permanently for this solve. Subsequent generic
   callback attempts may be recorded but must not invoke geometry or source
   again. Generic catch-all must not turn resource exhaustion into convergence
   or infeasibility. Allocation failures during hashes/records/source mapping
   receive the same treatment. No incomplete record becomes a successful audit.
5. Recheck the final accepted state with the actual checker and recompute its
   material evaluation from that state; explicitly distinguish final recheck
   success, unavailable/incomplete and resource failure. Do not return a final
   guarded pass merely because a prior candidate passed. Numerical convergence
   remains a numerical statement only; feasibility UNKNOWN, no seal/physics/
   stable-spin/original-paper grant even if geometry checks pass.
6. Bind raw state, actual source/config/frame, all settings/budgets, actual
   source-output hashes, guard attempts/stages, solver trial/acceptance history,
   failures and final material output in scientific digests. No timestamps.
   Preserve enough typed raw evidence to independently link every accepted
   solver state to its matching passing check and every failed/shrunk trial to
   its actual guard result. Counts distinguish geometry value-only evaluations
   from full derivative evaluations and final rechecks. Do not claim a total
   source-evaluation count which omits geometry evaluation work.

## Resources / compatibility

Retaining full audits per trial is potentially quadratic. Preflight all counts,
checked products, index narrowing, owned state/audit/record/string capacities,
return-copy growth and existing adapter/checker coexistence before allocation.
Finite explicit caps for guard calls and total retained guard evidence are
required. Prefer sequential one-current-audit scratch plus bounded compact
retained evidence; if full audits are kept, account every dynamic payload and
copy. Existing checker known-buffer cap cannot substitute for total solve budget.
No dense state-by-state all-geometry copies or silent cap lowering. Opaque
CGAL/GMP/HiGHS allocation internals still need external process guards, not an
exact RSS assertion. API preflight may throw; attempted failures should retain
complete evidence where memory permits, otherwise propagate without pass.

Legacy scientific outputs must remain identical when guarded route is unused.
If factoring changes resource-accounting identity, explain and ask first; do not
silently rewrite old expected output hashes. Current SourceCageSlqp initial
exception handling and generic catch-all deserve particular attention.

## Actual tests and independent verification

TDD: add registered failing behavioral test before production. Use actual
source octree/conforming/full26-handle BBW fixtures (Task5u source depth3 is an
explicit synthetic engineering fixture). No fake accepted geometry, mocked
predicate or injected known solver direction may count as integrated coverage.

- All five objectives remain supported; actual runs retain all113/112 variables.
- Guarded initial geometry rejection/source failure/resource failure accurately
  distinguished; no fabricated accepted result or KKT.
- At least one actual solver-generated trial rejected by geometry, followed by
  observed reduced step and a checked accepted state; validate actual step/state
  relationship and retained cause. Another trial-exhaustion case must show last
  valid state/moments/history unchanged. Restoration trials cannot bypass gate.
- Exact/one-below/zero/overflow count and memory caps; transient allocation
  failure if current test utilities support it. Resource latch stops later work.
- Final check + recomputed moments/state hashes; two identical full runs.
- Unguarded existing tests and archived Task5t scientific output remain unchanged.
- Independent controller checks link raw solver/guard histories and reintegrate
  exported actual geometry; corrupted/dropped/mismatched evidence must fail.

Synthetic fixture calibration is allowed only to exercise a named mechanism,
never to manufacture an original-paper result. Predeclare a bounded sweep if
needed, preserve all failed probes and chosen rationale; no open-ended hunt,
relaxed tolerances or geometry repair. Tell controller if an actual integrated
negative/rollback case cannot be achieved. Task5u exact-contact and positive-volume
pair-overlap end-to-end negatives remain open; do not claim current integration
closes them unless it actually does. Domain/outer clearance is conservative
for all cells regardless beta; not continuous BBW Jacobian or cavity-web thickness.

Only build/task-5c-no-faster-20260909, 2014ON/exactON/FasterOFF,
LIBIGL_SKIP_DOWNLOAD=ON. Baseline geometry+sourceSLQP tests run by controller.
Compile at most3 jobs (controller may use1), all heavy tests sequential,
OMP/OPENBLAS/MKL1, ulimit-v25165824 and timeout3600. Focused tests while
iterating; after freeze relink scoped2014/Python targets and run registered
related group once. Inspect actual fresh exactOFF child graph exclusion/current
sourcePython evidence; preserve logs, source/core/module/binary hashes.
No subagents or worker-owned reviewer. Read implementer prompt/TDD instructions.

## Approved proposal before production

Privately constructed SourceCageGuardedSlqpAudit with const facts wrapping
legacy numerical diagnostics and fixed-size guard attempt rows, plus a flat
ordered vector of original SourceCageGeometryRecord stages. Attempts index
existing raw callback-state tape; no per-trial geometry matrices. Common private
runner in SourceCageSlqp.cpp takes internal GuardContext pointer; legacy null
preserves resource/config/scientific identity. No extra implementation file.
Defaults IMPLEMENTATION_DEFAULT:1024maximumGuardCalls,65536retained stage
records,2GiB maximumGuardedBytes. Check full declared budget before allocation,
include nativechecker/adapter coexistence and record/string/copy capacities.
Count checker invocations separately from actual geometry source evaluations:
preflight failure can precede SourceEvaluation. Also count full derivative and
final-recheck work. Existing deliberate membership/history resets are distinct
from a failed candidate illegitimately updating accepted state/curvature history.
