# Actual-source guarded SLQP

`solveSourceCageGuardedSlqp` runs the immutable `SourceCagePaperFunctional`
through `SourceCageGeometryChecker` before each full derivative evaluation.
It uses the same private source adapter runner and audited KKT/restoration/
L-BFGS engine as `solveSourceCageSlqp`. All beta variables, the optional phi,
and all original `(scale, tx, ty, tz)` handle coordinates participate.
There is no user geometry callback, supplied mesh, supplied evaluation, or
acceptance boolean. Construct a new source functional to rebind upstream data.

The returned `SourceCageGuardedSlqpAudit` has a private constructor, exposes
only const facts, and is move-only. Moved-from audit access throws: a stale
passing scalar cannot be read beside transferred evidence. An extracted copy
of numerical diagnostics or facts cannot construct a guarded audit.

## Admission and termination

Only a completed checker audit with `ChecksPassed` admits the matching raw
state to the full derivative evaluation. Rejected geometry and `Unknown`
both invalidate that solver trial. The existing merit search, first-order
fallback, and restoration backtrack on the failure. Unaccepted trials never
replace the accepted objective, moments, state, or curvature pairs. The
engine's deliberate active-membership history resets remain its normal
iteration bookkeeping, distinct from failed-candidate updates.

An initial failure retains the attempted full state, stage evidence, and
source failure cause without manufacturing a numerical iteration, final state,
termination, or KKT result. Source value failure is `Unknown` geometry, not
`RejectedGeometry`. A checker resource failure or adapter resource failure
sets the solve's permanent resource latch. Later generic callback attempts
are recorded but execute neither checker nor full source. The resource
outcome overrides any underlying numerical termination, including a numerical
convergence reached before a failed final operation.

The last accepted state receives a new actual checker invocation and then a
new full material/derivative evaluation. `finalStatus` distinguishes a fresh
pass from unavailable, rejected, incomplete, and resource outcomes. A prior
trial pass cannot supply the final pass. The legacy nested `source` result
keeps its `NotEvaluated` geometry identity; only the private guarded root
owns the actual geometry history. Feasibility remains `Unknown`. Numerical
convergence means first-order tolerances only. These results grant no seal,
stable spin, physical behavior, original-paper case, or original reproduction.

## Raw evidence and counts

Each `attempts[i].callbackRecord` identifies the actual state column and
source evaluation record in `source.callbackStates` and `source.records`.
The attempt's `[firstStage, firstStage + stageCount)` interval addresses
the original ordered checker stage records in `stages`, including primitive
IDs/beta ownership, outcome, exact clearance string/witness, and evidence
and full-reason hashes. The complete original checker scientific hash and
its actual source-value output hash are retained. `checkCompleted` means
the checker audit was completely published; its outcome still must pass.

The first solver callback is the initial evaluation. Thereafter concatenate
each numerical trace row's `objectiveTrials`, then `restoration.trials`, in
row order to link callback states to attempted solver steps. The final
accepted-state callback has its own purpose and is not a solver trial.
For restoration, raw state equals the row state plus recorded step times
recorded restoration direction. Objective backtracks within one direction
have the same relation, with the beta prefix clamped to `[0,1]`; the fallback
may start a new direction. Accepted trace rows must match an earlier callback
whose geometry passed and whose full evaluation completed.

`checkerInvocationCount` counts actual checker calls, including preflight
failures. `geometryValueEvaluationCount` counts returned checker stage tapes
that reached `SourceEvaluation`, including value failures. It excludes checks
that failed preflight. `source.sourceEvaluationCount` counts full derivative
calls, including failures and final material recomputation.
`finalRecheckCount` counts only actual final checker invocations; a final
attempt skipped after a resource latch contributes zero. If a checker throws
before publishing its stage tape, `geometryValueEvaluationCountComplete` is
false: the known value count is a lower bound and no total source count is
claimed. Otherwise value plus derivative calls gives the full actual source
evaluation count; it is never the derivative count alone.

Scientific digests bind the raw initial/callback/accepted states, actual
source/configuration/frame, all solver/geometry/budget settings, source-output
hashes, ordered stage and trial histories, failure causes, and final material
output. No timestamps participate. The unguarded source adapter's public
structs, accounting formula/version, config hash, and scientific serialization
are unchanged when the private gate pointer is null.

## Simultaneous resources

All defaults and conservative factors are `IMPLEMENTATION_DEFAULT`.
The guarded planner first applies the unchanged legacy planner and actual
checker settings validation. Checked arithmetic bounds the actual primitive
pair count, `3 + 3P + P(P-1)/2` stages per checker call, the declared guard-call
budget times that stage count, stage and attempt vector capacities, and
primitive count times the exact-string bound. No cap is silently lowered.
The complete declared stage capacity must fit `maximumRetainedStageRecords`.
Guard calls may exhaust at runtime, preserving their actual last valid state.

The known-buffer envelope adds:

- The complete legacy adapter envelope, including source outputs, source
  workspace, raw state tape, generic trace/hash/LP workspace, and history.
- Preallocated flat stage and attempt capacities, two times the bounded exact
  clearance string payload per primitive per declared guard call, and fixed
  metadata reserve. There are no retained per-trial point matrices or Nef
  geometries. The raw callback states reuse the legacy tape.
- Twice the declared checker buffer cap plus metadata reserve for the one
  current checker, its returned audit, record capacity, and diagnostic
  snapshots. This conservatively allows checker/audit copies to coexist
  with retained stages and all adapter buffers.
- Guarded settings/results/hash metadata and transfer overhead. The public
  audit forbids deep copies and returns by move. Caller-created copies of
  extracted facts are caller allocations outside the solve's envelope.

Defaults are 1,024 checker calls, 65,536 retained stages, and 2 GiB total
guarded known buffers. The configured checker cap does not replace the
whole-solve cap. Opaque CGAL/GMP/HiGHS allocation internals and allocator
metadata are not an exact RSS assertion; use the external 24 GiB address-space,
one-thread, 3,600-second process guards. Invalid API/preflight caps may throw.
Attempted allocation failures retain resource evidence where recovery succeeds;
failures outside recoverable storage/hash publication propagate without any
returned guarded pass.

## Synthetic integration coverage

The registered test uses the actual source octree, conforming boundary, and
full 26-handle BBW pipeline on Task5u's explicitly synthetic depth-3 boxes.
It tests all five objectives (112/113 variables), two identical complete
runs, raw callback/step links, initial and final gates, resource latching,
exact/one-below/zero/overflow caps, and real allocator failures. Separate
registered modes test 48 actual SHA-allocation failure ordinals and moved-from
audit rejection. They use no mock checker and no injected solver direction.

The predeclared 25-case calibration used five objectives, clearances
`{0, .005, .02, .031, .04}`, one iteration, eight backtracks, beta
`[.17,.24,.31,.38,.45,.52,.59,.66]`, identity handles, phi `.23` where present,
frame origin `(.48,.51,.3)`, and weights `.8/.6/.03/.04/.05` for
inertia/center/fill/shape/displacement. All probe output is preserved in
the Task5v evidence directory; these are mechanism tests, not paper results.

The selected TopMassAndLowering zero-clearance case has source-value
`Unknown` trials at steps `1` and `.5`, actual primitive-subset geometry
rejections at `.25` and `.125`, then a checked accepted `.0625` step.
The separate approved four-case exhaustion sweep used that objective,
clearance `.04`, and backtrack counts `{1,2,3,4}`. The selected four-backtrack
case rejects objective and restoration geometry trials, reports
`RestorationStalled`, and preserves the initial state, moments, objective,
and empty curvature history. The YoYo `.04` case exercises actual restoration
rejection followed by accepted step `.015625`.

This integration does not close Task5u's open exact-contact or positive-volume
primitive-pair-overlap negative cases. Outer/domain clearance is conservative
for every cell regardless of beta. It does not establish a continuous BBW
Jacobian condition or cavity-web thickness.
