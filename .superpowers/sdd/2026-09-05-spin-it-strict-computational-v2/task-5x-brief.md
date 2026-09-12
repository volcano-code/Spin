# Task5x — actual-source guarded adaptive multi-pass orchestration

Approved continuation of Task5 (parent extraction alongside this file), not full
Spin-It reproduction. BASE d81bc808806fbff24fd40fad4c6d5e16e5dc3dc0.
Canonical worktree /data/xvxinyao/moment-guided-spin-v1/code/make-it-balance-paper-sim.
Main sibling checkout is read-only. Existing foreign Faster tracked/untracked
work is excluded. No new worktree, dependencies, downloads, merge, push or tags.

## Scope and interfaces

Add SourceCageAdaptiveSolver.h/.cpp and one native regression plus scoped CMake
registration and docs/spin2014_source_cage_adaptive_solver.md. Use existing actual
SourceCagePaperFunctional, solveSourceCageGuardedSlqp and rebindSourceCageAdaptive.
Do not replace legacy fill-only SourceAdaptiveSpinSolver. No binary/seal/Python
facade or paper asset work in this increment. No modifications to Faster files.

Every pass starts a fresh guarded SLQP solve with full beta/phi/all handles. Only
numerical Converged plus final freshly passed geometry AND material evaluation
can enter rebind. Preserve raw numerical termination, callbacks, KKT, trial and
geometry records for every attempted pass. Incomplete/limits/stalls/geometry
failures are never infeasibility. No user callbacks/booleans authorize progress.

Successful rebind changing the actual partition feeds its private immutable
checkpoint into the next fresh solve, resetting history. Preserve phi and every
handle by actual offsets. No copied old KKT/geometry grant for changed partition.
If no partition changes, require input/output identity and no refinable fractional
leaf before reporting relaxed adaptive completion; maximum-depth fractional
leaves may remain, and this is NOT binary validation, feasible witness, spin
stability, geometry seal, or paper-case completion. State this in the API/doc.

Keep last successfully solved source/state separate from a pending changed
rebound checkpoint not yet solved. At pass limit or later failure preserve both
honestly; failed rebind candidates must never become published checkpoints.
Immutable owned audit/endpoint access; a copied facts struct grants nothing.
Record each pass number, old/new partition/state/functional hashes, guarded solve
and rebind scientific hashes, starting dimensions and fresh history, cumulative
resources and precise terminal reason. Scientific digest excludes time.

Bound aggregate retained solve/rebind evidence and checkpoint payload across
passes BEFORE admitting work, not only per-pass caps. Use checked arithmetic and
explicit ResourceLimit (preserving prior evidence). State exact coverage of
budget accounting; upstream immutable allocations/checker scratch still have
their existing independent caps and external process memory limit. Avoid
pretending a known-buffer bound is an RSS guarantee. If needed escalate a clear
bounded resource design before coding instead of hiding unaccounted retention.

## Tests, honesty, and execution

TDD: capture real expected RED before production code. Cover actual multiple
passes, fresh histories and identities; unchanged convergence; max passes with
pending state; actual nonconvergence prevents adaptation; initial geometry reject;
rebind failure after successful solve; aggregate resource exhaustion and overflow;
immutable lifetime/move behavior. Real callbacks/checkers, no injected pass flags.
Use full BBW and actual source fixtures. If positive multi-pass test requires zero
weights or known stationary input, label KNOWN_SOLUTION_VALIDATION; it only tests
orchestration, not nontrivial optimization or original-paper convergence. Keep
nontrivial failure cases, no tolerance relaxation to manufacture convergence.

Use build/task-5c-no-faster-20260909 only (2014/exact ON, Faster OFF, skip download).
Compile at most -j3; native runs sequential, coordinate ownership with controller.
ulimit -v 25165824; OMP/OPENBLAS/MKL_NUM_THREADS=1; timeout3600. No full related
suite needed on worker: controller owns final all-target build and related suite.
Controller owns task-5x-controller-* artifacts, public remaining-gap doc, ledger
and independent evidence checks. Do not touch those or old evidence bytes.

Read relevant skills and implementer prompt. Do all work yourself, no subagents
or reviews. Ask for uncertain design decisions. Commit only explicitly owned
production/test/doc paths, user identity volcano <xuxinyao307@gmail.com>. Do not
stage progress.md, generic task-5-brief.md, other untracked work or failed logs.
Write full report with RED/GREEN commands, results, limits, files, self-review to
task-5x-report.md in this workspace. Controller dispatches fresh review. Return
short status, SHA, test count, concerns and report path (under15lines).
