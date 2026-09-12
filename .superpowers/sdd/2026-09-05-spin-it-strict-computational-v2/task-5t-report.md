# Task5t implementation report

## Scope and implementation

Base: `fb6aace8884cfc98dde4ceded9f51a6f1afac201`, branch
`spin-it-strict-computational-v2`, canonical worktree
`/data/xvxinyao/moment-guided-spin-v1/code/make-it-balance-paper-sim`.
Only the five owned implementation paths are changed:

- `include/spin2014/SourceCageSlqp.h`
- `src/spin2014/SourceCageSlqp.cpp`
- `tests/spin_it_2014_source_cage_slqp_regression.cpp`
- `docs/spin2014_source_cage_slqp.md`
- `CMakeLists.txt` (six additions, all under 2014 AND exact gates)

The report is separately committed. Protected main and concurrent foreign Faster
work are excluded. No generic SLQP, LP, KKT/history, source geometry, paper,
shape, legacy API, dependency, or existing test source was modified.

Implementation commit: `8806dd6f8380ac12df717194e7525071bc583f59`,
`feat(spin2014): bind audited joint SLQP to actual source cage`, author
`volcano <xuxinyao307@gmail.com>`. The staged implementation diff contained
exactly the five owned files (775 inserted lines).

The typed adapter requires `AuditedKktRestorationV1` and invokes the existing
`EqualityBoxProblem` / `solveSlqpLbfgs` algorithm with the actual immutable source
paper functional. Full states preserve beta, optional raw phi, and all original
handles. Only beta is bounded; scale-domain violations remain invalid source
trials without clamping. Generic objective and gradient/J map the full paper,
fill, shape and displacement terms; generic moments use the fixed spin frame.
Final source geometry and world/spin moments come from explicit evaluation at
the accepted final state. Static phi diagnostics are corrected to zero and the
generic scientific trace digest is recomputed.

Every returned result has geometry `NotEvaluated` and feasibility `Unknown`.
There is no sealed/material-release admission, minimum claim, second-order
certificate, continuous geometry legality, author trajectory, original paper
case, or full original-paper SLQP claim. Partition/BBW rebinding requires fresh
immutable source and solve/history. Adaptive transition orchestration and
precommit geometry guards remain downstream work.

The preallocated tape records full state columns, deterministic sequence and
purpose, actual source-call flag, state/output or full-reason digest, bounded
display reason and category. Callback counts include invalid/resource attempts;
actual source-call counts include final accepted-state evaluation. Upstream
identity/frame/weights/density, initial state, solver/LP settings, all budgets,
trace and final source geometry/moments bind scientific identity. Elapsed time
and timestamps do not. Identity is scientific provenance, not executable/source
attestation or manual approval.

Resource planning checks N=B+a+4H, K=2/5/6, T=I+1 and A=2+3IL before allocation.
The public breakdown includes full tape, source outputs and shape/split scratch,
generic trace growth and whole-result copy, L-BFGS history, solver/Eigen packing,
explicit LP arrays/model copies, full generic hash serialization copies, and
result storage. Adapter hashing streams; it does not introduce N×N arrays or
duplicate mapped world/J matrices. Upstream immutable storage/construction caps
remain separate. HiGHS opaque internals and allocator bookkeeping are not
claimed exact; all runs use the mandated 24 GiB process ceiling and LP controls.

`ResourceLimit` latches through generic catch-all trial recovery, prevents further
source calls, and keeps returned numerical state/trace. A final resource failure
cannot promote an underlying `Converged` status. If the very first actual source
call fails, the typed result retains its attempted-state tape and initial state,
with `numericalResultAvailable=false`; no termination, final state, trace, or
KKT witness is invented. Ordinary first-call failure has
`InitialEvaluationFailure`; resource failure has `ResourceLimit`. Only preflight
resource/state/settings errors before a source call throw without a tape.

## TDD and focused evidence

The complete task brief, implementer prompt, TDD skill and writing-good-tests
reference were read. API and simultaneous-buffer accounting were agreed with
the controller before production implementation.

Registered missing-adapter RED (before production header/source):

```sh
cmake -S /data/xvxinyao/moment-guided-spin-v1/code/make-it-balance-paper-sim \
  -B /data/xvxinyao/moment-guided-spin-v1/code/make-it-balance-paper-sim/build/task-5c-no-faster-20260909 \
  -DMAKE_IT_STAND_BUILD_SPIN_IT_2014=ON \
  -DMAKE_IT_STAND_BUILD_EXACT_GEOMETRY_CHECKS=ON \
  -DMAKE_IT_STAND_BUILD_SPIN_IT_FASTER=OFF -DLIBIGL_SKIP_DOWNLOAD=ON
cmake --build build/task-5c-no-faster-20260909 \
  --target spin-it-2014-source-cage-slqp-regression -j3
ctest --test-dir build/task-5c-no-faster-20260909 --output-on-failure -j1 \
  -R '^spin-it-2014-source-cage-slqp$'
```

`build/task-5t-red.log`: **0/1 passed**, expected runtime assertion
`FAIL: typed actual-source audited joint SLQP adapter is not implemented`.
The header-present placeholder was never run as GREEN. Actual-source behavioral
tests were written before production code. `build/task-5t-configure.log` and
`build/task-5t-red-build.log` retain configuration and build evidence.

First substantive GREEN: `build/task-5t-green.log`, **1/1 PASS, 5.43s**.
Expanded initial/final-resource and convergence boundary GREEN:
`build/task-5t-green-final.log`, **1/1 PASS, 5.77s**.

The controller identified that a first source callback exception lost the local
tape when the generic engine threw before producing a result. The native test
was changed first to require returned evidence for an actual source paper
`maximumJointDerivativeValues=1` failure. Same registered command then produced
`build/task-5t-initial-failure-red.log`: **0/1 passed, 5.72s**, with
`RESOURCE_LIMIT: source cage paper combined joint outputs and fill scratch budget`
escaping instead of the required typed result. The small adapter-only fix added
`numericalResultAvailable`, `InitialEvaluationFailure`, and return of the already
recorded first attempt. No engine change was needed.

Targeted final GREEN: `build/task-5t-initial-failure-green.log`, **1/1 PASS,
5.69s**. Its build is `build/task-5t-initial-failure-green-build.log`.
The test also exercises ordinary first-call numerical overflow using finite
positive handle scales of 1e200 and verifies failure provenance, absence of
manufactured generic output, and deterministic repeated initial resource failure.

All test executions/builds use `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
MKL_NUM_THREADS=1`, `ulimit -v 25165824`, `timeout 3600`, at most three compilation
jobs and serial experiments/tests.

## Actual-source regression coverage

The fixtures use original outer cube [.3,.7] and design cube [.34,.66], plus an
independently specified rectangular design [.34,.36,.38]–[.66,.64,.62]. They build
the real source octree, adaptive domain, full cage/grid/FEM/BBW and conforming
boundary. They are synthetic validation inputs, never original paper cases.
Eight betas and 26 original handles produce 112 static or 113 rotational state
coordinates. All five objectives start with nonuniform beta, positive nonunit
scales, distinct small translations, and nonzero fill/shape/displacement weights.
Fixed settings declared before observations are six iterations, twelve search
steps and L-BFGS memory four, with standard solver tolerances.

Observed five-mode final statuses and tape counts:

| Objective | Numerical status | Records | Invalid source trials |
| --- | --- | ---: | ---: |
| StaticLowering | MaximumIterations | 12 | 0 |
| TopMassAndLowering | MaximumIterations | 26 | 16 |
| TopInertiaOnly | MaximumIterations | 26 | 16 |
| TopLoweringOnly | MaximumIterations | 20 | 12 |
| YoYo | MaximumIterations | 11 | 3 |

Each mode has an accepted objective/restoration step, nonzero beta participation
and nonzero handle participation. Tests compare every available trace state
against direct source objective/gradient/J evaluation and independently derive
KKT equality, bound, stationarity, dual-sign and complementarity residuals by
scalar algebra using the reported multipliers. They check full membership and
raw phi, all final geometry rows/columns, moments, summed regularizers and actual
final-state equality, repeated full-state/event/scientific identities, raw trace
corruption and false-convergence rejection. Every successful callback output is
re-evaluated directly at its recorded raw state.

Additional cases cover one-iteration nonprogress, fixed symmetric first-order
convergence with unknown geometry, callback exhaustion, final-source exact and
one-below caps, underlying convergence plus final resource failure, all zero and
overflow caps, exact/one-below state/trace/tape/byte caps, bad dimensions/beta/
NaN/nonpositive scale, forbidden algorithms and invalid settings, changed
source/frame/weights/density/budget identity, and fresh history after rebinding.

## Final covering verification

The first covering run completed **39/39 PASS, 255.27s**, with the source files
matching `build/task-5t-final-source-begin.sha256` and
`build/task-5t-final-source-end.sha256` exactly.
`build/task-5t-covering.log` SHA256 is
`80960c564963e257a7aeeb1b87db1dec2591df3b69648b5e875cd663baac9385`;
the fresh full raw log `build/task-5t-covering-raw.log` SHA256 is
`adb590fd2dc49e026982b0fa62441c5930fb2e5657d4e84333dcfb11f07c06c5`.
Those are pre-allocator-hardening results, not claimed as final source evidence.
The pre-hardening manifest SHA256 is
`ccd6fea7abfbc4647985a51ccc7c928c59b8a54a633ff89a4b23251aa1237d73`.
Core archive SHA256 at that freeze was
`c8b535e1778919357acfbd2f9308564660771a2fac663c26254a4c195775e920`;
the native binary SHA256 was
`0f92f70f82461b050969ad955c03eec5a70b92b33c43aadc1133a1fc6ed2a9d5`.

## Allocator provenance hardening: additional real RED/GREEN

Self-review identified a narrow actual allocator failure case in the adapter's
callback state hash: its small digest string allocation could fail after the
record count increment but before the event status/output had been filled.
The generic catch-all line search could continue and eventually return a
ResourceLimit result whose attempted record was still `Unused`. This was
reported to the controller before altering the frozen source. The controller
required a targeted hardening and a real test-only allocator fault regression.

Only the existing owned test file was extended, with a global new/delete wrapper
and thread-local one-shot fault injection for real 65-byte allocations. The
bounded sweep calibrates actual digest allocation ordinals; no production hook,
replacement geometry or user callback exists. Against pre-hardening production,
the registered test failed at **allocation ordinal 8, event 1 remained Unused**:
`build/task-5t-allocator-red.log`, **0/1 passed, 5.68s**. Its compiler evidence is
`build/task-5t-allocator-red-build.log`.

The adapter now marks a fixed resource status/reason without allocation in the
callback recovery path. Finalization fills missing state and reason digests if
memory is available, also hashes the raw state column, and rejects incomplete
attempt records with an explicit RESOURCE_LIMIT exception if proof finalization
cannot finish. No incomplete returned proof or successful outcome is allowed.

Focused hardened GREEN: `build/task-5t-allocator-green.log`, **1/1 PASS, 6.34s**;
build log `build/task-5t-allocator-green-build.log`. The real fault sweep exercised
**30 one-shot allocation failures**, with **17 returned complete ResourceLimit
results** and the remaining failures reported by explicit resource exceptions.
All five actual-source modes and earlier tests continue to pass.

The hardened source freeze manifest is
`build/task-5t-hardened-source-begin.sha256`, SHA256
`ccd212fe118535a4f91be2754c0c385a38a3bd94d705b4826087e413b929493a`.
Core SHA256 is
`3d174485c2a0aad05e6cdcebdc39c080ce4249dd102860f467d4680e89f68105`;
new native binary SHA256 is
`dfb05755e5c810942f41d85e9f206ef56b5239e92cc7a5992dac45640a0acfe1`.
The required scoped relink and covering group completed once against this
hardened source: **39/39 PASS, 253.73s**. All 31 selected link operations were
existing 2014 native targets or the Python module; no unrelated target or
foreign source was compiled. The new native test passed in 6.34s. No skipped,
NotRun, failed or disabled test is counted as a pass.

```sh
mapfile -t task5t_targets < <(ninja -C build/task-5c-no-faster-20260909 -t targets all \
  | sed -n 's/^\(spin-it-2014[^:]*\): CXX_EXECUTABLE.*/\1/p')
cmake --build /data/xvxinyao/moment-guided-spin-v1/code/make-it-balance-paper-sim/build/task-5c-no-faster-20260909 \
  --target "${task5t_targets[@]}" _make_it_stand -j3
ctest --test-dir build/task-5c-no-faster-20260909 --output-on-failure -j1 \
  -R '^(spin-it-2014-|spin-2014-|source-spin-material-python$)'
```

The target list was obtained with `ninja -C build/task-5c-no-faster-20260909 -t
targets all`, selecting only `spin-it-2014` CXX_EXECUTABLE targets; `_make_it_stand`
was explicitly added. Complete link output is
`build/task-5t-hardened-covering-relink.log`. Test summary
`build/task-5t-hardened-covering.log` SHA256 is
`dd7c26a37584cf26269b3d166d950c55f5a2ae7827869ea11c00c42f4419f306`.
Fresh full raw CTest output was preserved in
`build/task-5t-hardened-covering-raw.log`, SHA256
`09cc4dcc4c3834816105218dc2170d0d797a1e6e782fef0961043fb35d6f4d31`.

Current Python success evidence in that final raw log: lines 459–461 report
**10 tests OK**, 23.169s, for `source-spin-material-python`; lines 975–977 report
**9 tests OK**, 0.032s, for `spin-it-2014-python`; lines 996–998 report the
existing paper-miniature **1 test OK**, 0.073s. These are current-run logs, not
prior-task files. The existing miniature test does not confer a physical or
sealed grant on this new adapter. The current Python module SHA256 is
`5b1092d3baf5dafca9cfe8bd27df1469e583f1b28bef750aaa6302b4c0b34453`.

`build/task-5t-hardened-source-end.sha256` matched the begin manifest exactly;
both have SHA256
`ccd212fe118535a4f91be2754c0c385a38a3bd94d705b4826087e413b929493a`.
Core and new native binary still match the frozen SHA256 values above. This
asserts scoped source quiescence, not stability of concurrent foreign files.

The actual exact-OFF child from the final registered dependency run is
`build/task-5c-no-faster-20260909/spin-2014-dependency-contract-runs/run-c2585319cac27b2a/enabled`.
Its cache records 2014 ON, exact OFF, Faster OFF and downloads skipped. Its
`build.ninja` SHA256 is
`b0e633950dc36b6a72935caf259d8c89ac62ae55601f2791460a4c84b358149f`.
Actual graph inspection confirms the new adapter/header/test names are absent.
Concurrent foreign Quadric/Canonical/FullRank/Plane/Polyhedral/Truncation source
and corresponding test paths are absent from that child and the assigned graph.
Evidence: `build/task-5t-hardened-exact-off-inspection.log`, SHA256
`418b49b65f3bddaa3d28ace0f1984dbf52eb2d759393adb4d802dcd0d39cd7eb`.

## Independent controller evidence

The controller independently linked its exporter against final core
`3d174485c2a0aad05e6cdcebdc39c080ce4249dd102860f467d4680e89f68105`.
The retained directory is `build/task-5t-controller-xXmsd4`. Two native exports
passed independently authored trace checks (10/10 each) and final endpoint
raw-coordinate QR/tetra-moment/complex-step checks (10/10 each). Exact `cmp`
passed for `run1.txt` versus `run2.txt` and `trace1.jsonl` versus `trace2.jsonl`.
All six-step cases honestly remain MaximumIterations, with no geometry grant.
The checker rejected **46 trace corruptions and 36 material corruptions**.
Evidence is in `native{1,2}.log`, `trace-check{1,2}.log`,
`material-check{1,2}.log`, `trace-negative1.log` and `material-negative1.log`.
The implementer inspected the resulting logs and exact-repeat comparisons;
the controller owns those checkers and any further review/status documentation.

## Self-review / remaining boundaries

The CMake diff contains only the two gated core entries and four gated target/
test/property lines. Generic engine APIs and algorithms are unchanged. Actual
allocation failure inside the generic engine can still throw without exposing
its in-progress trace, because the engine has no such recovery API; it is
converted to a RESOURCE_LIMIT exception and never success. Opaque HiGHS memory
and time-dependent LP termination remain explicit limitations.

The identified actual-OOM record bookkeeping edge was fixed and its real
allocator-fault RED/GREEN evidence is recorded above. Known configured-budget
exhaustion and source resource failures remain separately covered natively.

## Independent review fix: native trace coefficient indexing

Review identified that bounding N by `INT_MAX/4` did not establish K×N ≤ INT_MAX
for rotational objectives with K=5 or K=6. The unchanged generic
`SlqpAuditTrace.cpp` serializer uses `for(int i=0; i<a.size(); ++i)` for the full
constraint Jacobian. Defaults and the existing small fixtures were unaffected,
but the adapter's accepted-cap contract was too broad. The finding was verified
against that actual upstream loop before implementing the correction.

Fix commit: `cccc1b9ee1a548e2830cdad9115f459258dd1e41`,
`fix(spin2014): bound source SLQP native trace coefficient indices`, author
`volcano <xuxinyao307@gmail.com>`. Exactly four owned files changed: source,
header, regression and numerical documentation (19 inserted / 3 removed lines).
No generic engine, geometry, objective, target, tolerance, fixture or foreign
source was changed; CMake was unchanged.

The public `maximumStateValues` cap must now be ≤ `INT_MAX/6`, conservatively
covering all five objectives. Preflight also separately checks actual K×N with
checked multiplication against INT_MAX, in addition to its existing state-loop,
LP indexing, Eigen and byte limits. The header and numerical documentation state
this explicit public-cap constraint.

TDD used a small settings-only boundary case on the existing actual-source
112-state fixture: `maximumStateValues=INT_MAX/6` must permit resource planning;
one above must throw RESOURCE_LIMIT before allocation. No huge source mesh or
matrix was constructed to test this integer boundary. This test was added before
the production guard. Registered RED produced **0/1 passed, 5.37s**, failing
`one-above native serialization state cap rejected before allocation` because
the old production accepted that cap. Evidence:
`build/task-5t-review-int-red-build.log` and
`build/task-5t-review-int-red.log` (SHA256
`8cbceeb9ddca18b656db67d4b67c7b40a23e50641fbc1c7b462f85eeb75d7c00`).

```sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
ulimit -v 25165824
timeout 3600 cmake --build build/task-5c-no-faster-20260909 \
  --target spin-it-2014-source-cage-slqp-regression -j3
timeout 3600 ctest --test-dir build/task-5c-no-faster-20260909 \
  --output-on-failure -j1 -R '^spin-it-2014-source-cage-slqp$'
```

After the fix, the required affected native regression and relevant audited
solver test were rebuilt and passed **2/2, 6.35s** (source-cage 6.34s, audited
0.01s). The source-cage test still covers all five nonstationary objectives,
initial failures, configured caps and the real one-shot allocator sweep.

```sh
timeout 3600 cmake --build build/task-5c-no-faster-20260909 \
  --target spin-it-2014-source-cage-slqp-regression spin-it-2014-audited -j3
timeout 3600 ctest --test-dir build/task-5c-no-faster-20260909 \
  --output-on-failure -j1 \
  -R '^(spin-it-2014-source-cage-slqp|spin-it-2014-audited)$'
```

GREEN build: `build/task-5t-review-int-green-build.log`.
GREEN result: `build/task-5t-review-int-green.log`, SHA256
`8b33667eda88c8fc5850459a77025bf7580fa771cc89607dc2208b4e7e943732`.
Fresh raw result: `build/task-5t-review-int-green-raw.log`, SHA256
`e5563e93057030ef9e68e3067ee40b33e55b032d0a8c50f1cb3a5746e15c5a84`.

Scoped `build/task-5t-review-int-source-begin.sha256` and
`build/task-5t-review-int-source-end.sha256` compare exactly, with SHA256
`996d30fba57f9bf31dbec5fdfe4d9111fecc33f582f0e980d1e890a7806a1824`.
Final review-fix core SHA256:
`c3d33e9fb849ce959e54a9eeb24a51c24b4cdfbddd7c97151740b6d4a0a332be`.
Source-cage regression SHA256:
`76e1cc783ea28b1e469fc745520a654cf149c334ee0ae6f27d07da3955d82b4b`.
Audited regression SHA256:
`a6bafe9c26bfd93cc35c3f8bd27bedee128f8b65772a3c34457ef2a1082f3131`.

Per the controller's narrowed review-fix validation scope, the 39-test group was
not rerun for this guard. The 39/39, Python and exact-OFF results above remain
historical **pre-review-fix** evidence; the fresh post-fix native evidence is
the 2/2 run recorded here. The controller independently relinked and ran its
actual-source exporter twice against the final review-fix core in
`build/task-5t-controller-reviewfix-J6wJC4`. Both 10-trace and 10-endpoint checks
passed; 46 trace and 36 material corruptions were rejected. Exact `cmp` passed
for both repeated endpoint and raw trace files. The implementer inspected the
four check logs, both negative-check logs and the exact-repeat comparisons.
The exporter binary SHA256 reported by the controller is
`7e09a092b1682963006f22649e2b7ddd14e7b7f3f99779c88c34336a9010de98`.
The controller also confirmed unchanged scientific hashes relative to the
pre-guard cases. Scientific and geometry-admission limits are unchanged.
