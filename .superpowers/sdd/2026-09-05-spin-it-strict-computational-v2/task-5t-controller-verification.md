# Task5t controller verification — scoped acceptance

BASE fb6aace8884cfc98dde4ceded9f51a6f1afac201. Canonical existing linked
worktree /data/xvxinyao/moment-guided-spin-v1/code/make-it-balance-paper-sim,
branch spin-it-strict-computational-v2. Gitdir differs from common, no
superproject, no relevant AGENTS. Protected main tracked diff SHA
93e32af47ab1957a1efcb6716437656be129053f4269d2169378fcd95b05a35e.
/data1.4TiB free, foreign concurrent Faster/Canonical/etc changes excluded.

This continues the explicitly approved strict M0–M8 specification and the
previously named next joint solve; no new original-case or publication grant.
Skills: existing-worktree verification, SDD implementation/review, TDD and
independent verification. Existing dependencies already configured; no package
install or download. Worker spin2014_source_cage_slqp owns five scoped code/
test/doc/CMake paths, controller owns independent helpers/evidence/status.

## Source and integration audit

Local Spin-It2014 text build/strict-paper-sources/spin2014.txt Section5.4 and6.1
read: box fill variables, CoM/inertia equalities, active set/SLQP and L-BFGS,
joint fill/handle gradients and the two shape regularizers. This is formula/
numerical reconstruction with IMPLEMENTATION_DEFAULT controls, not recovered
author subproblems/globalization/trajectories. Original asset gaps unchanged.

Inspected complete generic SlqpLbfgsSolver, SlqpHistory and SlqpAuditLp and
source paper integration. Generic first-unbounded phi diagnostic needs static
correction; callback catch-all needs a persistent resource failure latch.
Existing source handles require strictly positive scales; no new clamp/bound
is imposed on generic unbounded coordinates. Numerical output remains geometric
NOT_EVALUATED / feasibility UNKNOWN until actual geometry admission exists.

Baseline commands:
- ctest --test-dir build/task-5c-no-faster-20260909 --output-on-failure -j1
  -R '^(spin-it-2014-source-cage-shape|spin-it-2014-audited-slqp)$': only existing
  shape1/1PASS1.08s; second regex name was nonexistent, not a two-test batch.
- Corrected -R '^spin-it-2014-audited$': actual1/1PASS0.01s.
- Core6450f4e751ca0b3a61c5c4a7ab9c4f22cc9193715f2fd96a6cb320a101ded854.

## Agreed resource design

Typed actual-source adapter reuses audited generic engine; no copied optimizer.
API/resource proposal and approval are captured in task-5t-brief.md supplement.
Preallocated callback tape capacity2+3IL includes final evaluation. Caps include
known generic trace/return copies, source outputs, history/solver scratch,
explicit LP arrays and full generic hash serialization. Opaque HiGHS storage
is not falsely called exactly bounded; process24GiB/LP limits remain required.
Callback failure strings must have deterministic categories/bounded excerpts
and full-message digest. No source evaluation after resource latch.

## Independent checks prepared before actual runs

Controller planned two existing independently specified source fixtures from
Task5s (8beta, all26handles), all five modes, six iteration limit fixed before
observations. Source/shape weights unchanged1.7/.6/.23/.37/.19, the same fixed
frame and nonstationary initial beta/scales/translations. No required all-mode
convergence; require actual accepted progress and retained nonconvergence.

task-5t-controller-check.py recomputes KKT by direct multiplier algebra, beta-only
bounds, phi semantics, actual accepted-state merit/restoration decrease,
membership hash and history reset constraints, final-state consistency and
first-order-only interpretation. It does not import the witness LP/solver.
task-5t-controller-negative.py mutates46 actual observations across both end
fixtures, including forged geometry/witness status and false residuals.
Python syntax checks passed; these are prepared tools, not scientific passes.
Final actual material endpoints will separately use the existing independent
Task5s point-frame/tetrahedron/raw-coordinate QR/complex-step oracle.
The exporter additionally evaluates every recorded accepted state and successful
callback directly through the source functional to check adapter mapping; this
is direct integration consistency, distinct from the independent endpoint oracle.
Call-state digests are independently recomputed from explicit IEEE big-endian
serialization; counters, order, final-purpose record and trace-state membership
are checked. Added six call-tape corruptions before any execution, no tolerance
or fixture changes. Native exporter now bound to actual header names, awaiting
the worker's stable core (no shared build by controller).

## First focused execution and initial-failure boundary

Actual expected missing-adapter RED at build/task-5t-red.log read by controller.
The header-present behavioral test branch was written before production; no
placeholder zero-return branch was run as GREEN. First real native GREEN1/1
passed5.43s, then strengthened converged/final-cap checks passed5.77s. Initial
core24d51bbd0be8bd3edd9c15f1d2f60fa5a49ff22b05355cd1365b3d22f3da67ac.

Controller identified a boundary mismatch: docs proposed throwing all initial
source failures, losing already-generated call tape, whereas the brief requires
all actual attempts retained. Clarification recorded before final verification:
preflight/bad-state-before-call may throw; actual first callback failure returns
its typed record with numericalResultAvailable=false and no invented engine
result. Worker adds real resource/numerical first-evaluation failure tests first.
Initial freeze cancelled before independent link or broad suite. No old result
promoted to amended source; previous logs retained. Worker hit one transient
model-capacity failure and same worker was resumed, not duplicated.

## First independent actual execution (pre allocator-hardening)

Core c8b535e1778919357acfbd2f9308564660771a2fac663c26254a4c195775e920
read before linking /usr/bin/c++17 O2/g/DNDEBUG/frounding-math exporter to
actual core/HiGHS/Boost/GMP/MPFR/zlib. New mktemp directory
build/task-5t-controller-7GEDzt. ulimit-v25165824 and OMP/BLAS/MKL1.
Native run1 exited0, two actual first-call resource/numerical-failure checks
passed with no invented numerical result. Generated10 endpoint cases and10
trace JSON lines. All source-successful callback hashes and every accepted
state's objective/gradient/J matched direct source reevaluation.

task-5t-controller-check.py trace1.jsonl:10/10PASS, each7 snapshots/6 accepted
steps, all MAXIMUM_ITERATIONS. Final first-order residual range roughly
0.00473954–1054.21: NOT convergence or stable physical geometry. Trace checks
validate recorded arithmetic/acceptance rules, not that six steps suffice.
task-5s-controller-check.py run1.txt:10/10PASS, same predeclared independent
tetrahedron/raw-coordinateQR/complex-step gates, max nonzero derivative relative
5.836e-9; shape relative2.604e-11, near-zero joint absolute7.614e-18.
task-5t-controller-negative.py trace1.jsonl:46/46 actual corruptions rejected.
No tolerance/input/weights change was made after observing these results.
Original logs are native1.log,trace-check1.log,material-check1.log,
trace-negative1.log in that directory. No determinism claim yet (one export).

Worker self-review identified a real narrow edge: failed small allocation in
state/hash serialization can leave an attempted record Unused after generic
catch-all. Controller inspected allocation/status ordering and requires
no-allocation failure tagging plus digest recovery or explicit resource error,
not a returned incomplete tape. Test-only one-shot allocation fault to be used
if feasible; no production injection hooks. Ongoing39suite remains historical
pre-hardening evidence; final source must be rebound and checked afresh.

## Hardened final-source verification

Real test-only allocator fault RED failed at ordinal8/event1 with an Unused
attempted record (5.68s). Hardened GREEN1/1 passed6.34s; source first-failure
RED5.72s and GREEN5.69s likewise retained. Controller read all these raw logs.
Final core3d174485c2a0aad05e6cdcebdc39c080ce4249dd102860f467d4680e89f68105.
Fresh independent export compiled into build/task-5t-controller-xXmsd4, using
the same compiler, libraries, caps, fixtures and controls as the earlier run.
Two complete executions exited0, both initial actual source failure checks
passed, and both full independent trace and endpoint checks passed10/10.
Endpoint and trace files respectively cmp-identical across the two runs.
46/46 trace/call/status corruptions and36/36 material/derivative corruptions
were rejected. No tolerance, objective or input changes were made.

Both independent runs again had7 snapshots/6 accepted steps per case and all
MAXIMUM_ITERATIONS. Largest nonzero derivative relative error5.836e-9,
shape2.604e-11, near-zero joint absolute7.614e-18. The reused Task5s checker
checks the endpoint layer, not optimizer convergence; the new Task5t checker
separately checks actual joint optimization traces. Neither is geometry admission.

Controller read final hardened covering log:39/39 actual PASS253.73s, no skipped
or NotRun records; this is the related2014/FasterOFF suite, not the foreign
Faster graph or entire project. Initial pre-hardening39/39PASS255.27s remains
historical, not substituted. Raw TDD, final covering, native and independent
check logs copied using apply_patch and cmp-verified byte-for-byte into this
SDD directory. Full numerical output remains under the named build directory.

Final independent scientific/build hashes:
- exporter21553761c407e479eb4ab44a9893f81611b5c098a8086908443be48bf555e543
- run1=run2:2f1d95f023f7068b120b10fc079286533aacc3daaac083c6a8b49cad5b90d092
- trace1=trace2:d18a8cbdc6a9ef2ac4a0221df402beaf60557a13585870956fa487bf59511abd
- native1=native2:435bb1ddf11a17372c949693a67b659f8a516c68196804a01ed74d342b80aeec
- trace-check1=2:1b001cfd012e5a3c85e03f07ddf9c0fea4cdc8c0cc8726e1212b295bc36cb190
- material-check1=2:b76046c16b8872eb2d8aafd8475972e25e18a14407913a6c96101fd75ad84d76
- trace-negative:ea9b37f6b66cbeb9390cb9315d3fda5e7bedb229526c4996eb3df67e35e292b5
- material-negative:0bda1819aedbc16f60e6bd7fb37c2ec956f84306460a37d8980dca3e0c36644d

Independent code review remains pending; no full task acceptance yet.

Implementation commit8806dd6 and report8caf013 inspected: only five owned paths
and report, no foreign paths. Reviewer independently reads the frozen package
review-fb6aace..8caf013.diff; math/resource/API risks require the SDD high reasoning
review seat. Task5t-report.md read fully, with current source Python10tests OK
and dependency-child raw evidence checked at their stated paths. Begin/end
manifests cmp-identical and sha256sum-c all five current files OK. Actual exactOFF
child build graph excludes this adapter/test and named foreign modules. Current
main tracked diff remains93e32af47ab1957a1efcb6716437656be129053f4269d2169378fcd95b05a35e.

task-5t-build-identity.log captures GCC13.3/CMake3.28.3/Python3.13.5/Ninja,
actual cache/rules/build graph hashes, linked native libraries and scoped source
diff hash03731b961af152adbde2a9a41894ce4db1cb5d77628eda32d8bb2e69f98db566.
This is scoped local verification evidence, not a clean release manifest or
automatic complete-reproduction admission. The predeclared synthetic inputs
are embedded in the exporter and bound by its SHA:
aafee522036ee3956a034de3a32c2983d57373aa0ff6ca579f51f8fdffa186c6.
Independent trace checker6d634ed3d904b5db4bc607103fb33d88ebd17d1f7d2c9502ab8880079c1679c7;
negative checker99d5c28440d1e395455b6ac1a4fb28d92d62b6b8e2a0d1c89e093dcaf3009c7e.

## Independent review and scoped integer-cap fix

Fresh reviewer review_spin2014_source_cage_slqp reviewed BASE..8caf013. One
Important finding: native generic SlqpAuditTrace.cpp:15 uses int for K*N matrix
coefficients but adapter only guarded N<=INT_MAX/4. Controller read the actual
upstream loop and agreed the K5/6 edge can exceed its range if resource caps are
raised; default/small-case results unaffected. Full review retained in
task-5t-review.md. Receiving-code-review skill used to verify the finding before
dispatching the original worker; no new implementation agent or engine rewrite.
Required small settings-only boundary RED/GREEN, conservative max state cap
INT_MAX/6 plus explicit checked serialized-product guard. No huge allocation.
39/39 and xXmsd4 results now explicitly pre-review-fix evidence; affected native/
audited tests and independent final-core rebind required, then scoped re-review.

Review fix cccc1b9: real settings-boundary RED0/1 failed5.37s, GREEN2/2passed
6.35s (actual source6.34s and audited0.01s), raw logs read and preserved. Five
owned begin/end manifests cmp-identical and sha256sum-c OK. New core
c3d33e9fb849ce959e54a9eeb24a51c24b4cdfbddd7c97151740b6d4a0a332be;
native76e1cc783ea28b1e469fc745520a654cf149c334ee0ae6f27d07da3955d82b4b.
Fresh compiler link into build/task-5t-controller-reviewfix-J6wJC4, exporter
7e09a092b1682963006f22649e2b7ddd14e7b7f3f99779c88c34336a9010de98.
Two real executions and all independent checks completed exit0:10trace/10endpoint
per run,46trace+36material corruptions rejected, endpoint/trace exactcmp. Native
first-source failure records also passed both runs. Full output/log SHAs equal
the previous xXmsd4 values above; only binary changed. All cases still six-step
MaximumIterations. No loss of failed trials, no relaxation and no convergence
or geometry admission. Original39suite expressly not rerun for this guard-only
fix; affected2test coverage and fresh independent verification are separated.

Fix report fdd85c2 read, package8caf013..fdd85c2 reviewed by an independent
scoped re-review seat. Verdict: finding ADDRESSED, no new breakage or out-of-scope
observations (task-5t-rereview.md). Controller accepts this actual-source joint
adapter layer only. Full original-paper SLQP, geometry commit gate, adaptive
joint orchestration, binary verification and original cases remain open.
All ten independent six-step cases remain MAXIMUM_ITERATIONS. Overall strict
reproduction PARTIAL_BLOCKED, no physical/sealed/release grant. No OpenGL result
added in this task. TDD and independent-review skills led to real first-failure,
allocator and native-index fixes before this scoped acceptance.
Preserved final review-fix raw log initially missed a trailing blank line during
copy; corrected with apply_patch and cmp-verified against the original. All
retained independent logs cmp-identical with the final review-fix run. Source/
test/checker assertions were not changed to alter scientific outcomes.
