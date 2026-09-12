# Task5t: actual-source audited joint SLQP integration

## Agreed API / resource supplement (before implementation)

SourceCageSlqpSettings contains audited-only solver settings and budget caps:
maximumStateValues100000, maximumSourceEvaluations10000, maximumTraceRows1000,
maximumCallbackRecords100000, maximumAdapterBytes512MiB. A public resource plan
computes the known-array envelope; positive caps and all checked sizes enforced.
Typed result keeps generic numerical termination separate from adapter
Completed/ResourceLimit, geometryNotEvaluated and feasibilityUnknown.

Capacity A=2+3*I*L includes explicit final evaluation, all attempted callback
states and latched skips, not only successful calls. Preallocate bounded tape
and state storage; no post-latch source work/unbounded allocation. Bound failure
strings with deterministic categories and a streaming full-message digest.
Check generic int and HiGHS n+2K capacity, in addition to size_t/Eigen limits.

Known-envelope components: source output world/point/J and paper outputs;
generic trace growing vector + finish return copies (3T+2 row envelopes);
L-BFGS 2MN+M; conservative32N+12KN+16K²+8M solver double workspace;
explicit LP assembly/model copies; whole generic hexfloat serialization with
string growth/copies (8 times upper length at32bytes per numeric field plus
bounded metadata); adapter fixed-endian hashes streamed. No exact claim for
opaque HiGHS internals; process24GiB ceiling and LP caps remain mandatory.
Upstream immutable geometry and previously separately bounded source/shape
scratch remain explicit separate budgets. Resource equation detailed in doc
and tested exact boundary/one-below, not merely a hand-set maxN.

Final source evaluation consumes source budget and has a separate-purpose tape
record. Unavailable final evaluation on resource failure must be explicit,
never fill a success-shaped default. Scale positivity remains source domain
validation, not a new bound/clamp in SLQP.

Initial-failure clarification: preflight/state validation BEFORE any source call
may throw. Once an actual initial source callback has executed, preserve its
typed record even when generic solver propagates that exception. Expose whether
a generic numerical result actually exists (false until solver returned);
on first callback failure return initial state + recorded reason/digest/tape,
no invented generic iteration/KKT/termination. Distinguish initial numerical
failure from RESOURCE_LIMIT, keep finalEvaluationAvailable=false. This is not
a generic engine change or a new feasibility classification. Include real
first-call derivative-cap and numerical-domain failure regressions.

Read this first. Binding scope is the user's approved strict M0–M8 reproduction
specification, specifically M3 joint beta/phi/full-cage solve. This is the next
bounded integration after accepted SourceCagePaperFunctional + shape terms.
It is NOT full original-paper reproduction or legal deformation admission.

## Workspace / ownership

Work only /data/xvxinyao/moment-guided-spin-v1/code/make-it-balance-paper-sim,
branch spin-it-strict-computational-v2, BASE fb6aace (resolve full SHA).
Protected main /data/xvxinyao/moment-guided-spin-v1/code/make-it-balance is read-only.
Foreign include/spin, src/spin, tests/quadric*, tests/*canonical*, *full_rank*,
*plane*, *polyhedral*, *truncation*, experiments/m1_* changes are excluded.
Do not modify, stage, build, clean, merge, or promote them. No new dependencies,
downloads, tags, pushes, worktrees or broad staging. No subagents.

Owned implementation paths:
- include/spin2014/SourceCageSlqp.h (new typed adapter/result/resource contract)
- src/spin2014/SourceCageSlqp.cpp (new implementation)
- tests/spin_it_2014_source_cage_slqp_regression.cpp (new actual-source tests)
- docs/spin2014_source_cage_slqp.md (new numerical status/limits/source choices)
- CMakeLists.txt (only new source/header/test under 2014 AND exact gate)

Existing generic SLQP, LP, history, source geometry/paper/shape/legacy APIs remain
unchanged unless you identify a concrete correctness blocker and first discuss
the smallest scoped fix with controller. Do not copy the entire SLQP engine.

## Required behavior

Expose a typed source-bound solve using immutable SourceCagePaperFunctional and
its full state. Default and require AuditedKktRestorationV1, never silently use
LegacyV3. Reuse existing EqualityBoxProblem + solveSlqpLbfgs outer algorithm.
Bounded prefix is beta ONLY. State is [beta, optional phi, all original handles
(scale,tx,ty,tz)]. Static has no phi: generic solver's first-unbounded-as-phi
diagnostic must not mislabel a handle. Result/trace phi corrected explicitly and
scientific hash recomputed. Preserve full raw state and all handles, no sample,
clamp of unbounded parameters, phi wrap, supplied weights/moments or toy callback.

Each callback evaluates the actual source paper functional, maps totalObjective,
paperObjective and sum(fill/shape/displacement penalties), jointGradient,
constraints/J and correctly framed mass moments into AdaptiveFillEvaluation.
Never discard regularizers or recompute target origin from moving CoM.
Final source evaluation comes from the final accepted state, not a trial/cache
selected for a nicer output. Output actual world/spin moments and geometry in
the existing source evaluation; do not call them a sealed material object.

Retain generic first-order KKT witness components, active membership, L-BFGS
history changes, rejected objective/restoration trials, accepted states, final
state and underlying numerical termination. Bind every callback invocation
(including invalid attempts) to its state + source output digest or failure
reason, with deterministic sequence and call count. This provenance is internal
scientific identity, not a source/build executable attestation or manual pass.
Bind config/upstream/source/BBW/frame/initial state/solver settings/budgets,
result state, final geometry+moments digest and generic trace into source result
identity. No timestamps/wall elapsed in deterministic scientific hash. All
settings that affect execution are explicit IMPLEMENTATION_DEFAULT and hashed.

Separate numerical termination and geometry-validation status. Every result
from this task has geometry NOT_EVALUATED and feasibility UNKNOWN. No result
can be passed to a physical/sealed artifact release gate on numeric status.
Converged means the existing FIRST_ORDER_KKT_TOLERANCES_MET only: no local-minimum,
second-order, continuous feasibility, author trajectory or original case grant.
Partition/BBW changes require new immutable functional and new solve/history;
adaptive transition orchestration, binary material validation and actual
precommit intersection/inversion/wall guards remain explicit downstream work.

## Resources / exceptions — agree API before implementation

Send controller proposed API and concrete simultaneous-buffer accounting first.
Use checked arithmetic and explicit finite positive caps for state dimensions,
callback evaluations, stored trace and adapter-owned buffers before allocation.
Count N=(B+angle+4H), K=2/5/6, generic (iterations+1) trace and L-BFGS memory,
callback call-tape arrays/digests, final output duplication and hash workspace.
Avoid dense NxN, full Hessian and duplicate world/Jacobian buffers; stream hashes.
Existing upstream geometry/derivative budgets remain separate. HiGHS opaque
internal allocations need not be falsely claimed exact; state the limitation,
keep process 24GiB ceiling and LP time/iteration controls. Accounting must not
ignore known Eigen packing or full trace/serialization copies.

Resource failures never become ordinary invalid geometry or convergence. The
generic engine catches callback exceptions during trial evaluation: adapter
must retain a resource-limit latch/evidence and must not report success after
one, even if the engine later terminates. Subsequent callbacks after exhaustion
must not evaluate source geometry or grow an unbounded tape. Preserve last
available numerical trace/state; if preflight cannot start, throw clear
RESOURCE_LIMIT length_error. Ordinary invalid numerical trial failures are
recorded and can be backtracked, without certifying infeasibility. Do not add
a user callback that can hand out an acceptance boolean or substitute geometry.

## Tests (TDD)

First execute a real registered RED for missing typed adapter. Then real native
actual-source fixtures with complete cage/FEM/BBW (reuse test mesh utility only).
Small synthetic boxes are independently specified validation inputs, NEVER paper
cases. Do not tune starting geometry, target, weights or tolerances to force a
successful original experiment. Declare fixed test settings before observations.

Cover all five objectives with nonstationary full states and nonzero fill,
shape/displacement weights; dimensions 112 static/113 spinning for 8beta/26handles
fixture. Verify actual nonzero handle and beta step participation where mathematically
expected, unbounded scale remains distinct from angle, bound residuals, final
evaluation equality against direct source at actual final state, complete trace
and KKT residual from independent direct algebra (not computeKktWitness itself).
Do not assert every mode must converge in a small test budget. Test honest
iteration-limit/nonprogress statuses; at least demonstrate accepted objective
or restoration steps rather than only exact known stationary initializations.

Two repeat runs: equal scientific hashes/states/events, excluding times.
Source/frame/weight/density/budget changes bind identity. Check all actual
handle positions preserved, no omitted columns. Static trace phi exactly0;
rotational phi tracks state at B. Fresh solve history starts empty, never reused
after source rebinding. Bad state/NaN/invalid algorithm/settings, zero/overflow/
one-below/exact caps and callback exhaustion remain honest. Raw trace corruption
or false convergence must fail independent checker tests. Geometric validity
must remain NOT_EVALUATED even on numeric convergence, no sealed/witness grant.

## Build / evidence / report

Assigned sole shared build: build/task-5c-no-faster-20260909, canonical -S/-B.
MAKE_IT_STAND_BUILD_SPIN_IT_2014=ON, EXACT_GEOMETRY_CHECKS=ON,
MAKE_IT_STAND_BUILD_SPIN_IT_FASTER=OFF, LIBIGL_SKIP_DOWNLOAD=ON.
Use correct option names from actual cache; no dependency repinning/rebuild.
Compile <=3 jobs, experiments/tests serial, OMP/OPENBLAS/MKL=1, ulimit-v25165824,
timeout3600. Do not run foreign tests. Iterate focused new test; once stable
notify controller of source/archive hashes for independent check. Freeze own
sources while controller checks. Relink selected relevant existing native and
Python targets once and run registered 2014/source-material group once:
ctest --test-dir build/task-5c-no-faster-20260909 --output-on-failure -j1
-R '^(spin-it-2014-|spin-2014-|source-spin-material-python$)'.
Inspect actual exact-OFF child graph excludes new adapter/test and retain current
Python success logs (not old file). No skip/NotRun disguised as pass. Record
begin/end scoped source hashes, core/test binary hashes and original log paths.
Foreign concurrent source changes outside actual graphs are not your edits.

Use apply_patch for source/file edits. Scoped implementation commit only owned
files (volcano <xuxinyao307@gmail.com>, authorized). Write detailed report to
.superpowers/sdd/2026-09-05-spin-it-strict-computational-v2/task-5t-report.md,
commit report separately. Controller owns separate independent checker/review,
ledger/public status docs. Return only status, commits, tests and concerns.
