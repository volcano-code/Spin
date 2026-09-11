# Task5k — actual source-cage BBW solves and independent FEM KKT

Implementation baseline: `29faa6423d96bccaf7b80a16a502c4686715884a`, isolated
`/data/xvxinyao/make-it-balance-paper-sim`, branch
`spin-it-strict-computational-v2`. This is the bounded Task5k increment; the
full Task5 joint material/shape optimization and original-case reproduction
remain open. No sourceFEM/grid replacement, legacy BBW changes, Python API,
geometry/material certificate, dependency update, helper agent, release,
merge, push, tag, historical-output deletion, or main-checkout write occurred.

## Interfaces and actual optimization

`include/spin2014/SourceCageBbw.h` and `src/spin2014/SourceCageBbw.cpp` add an
immutable actual-run result constructed only from immutable `SourceCageFem`
and configuration. Its const accessors expose raw/final matrices, separate
assessments, individual backend attempts, exact input energy scale, resource
counts, backend settings, config/output digests and terminal reason. There is
no constructor accepting candidate weights or successful report flags.
`assessSourceCageBbw` separately evaluates arbitrary candidates without granting
run authority. Invalid shape/nonfinite data yields `available=false` and false
acceptance flags. Failed runs expose no final validated field; incomplete raw
columns are explicitly identified by observed attempts, with raw full-field
assessment unavailable. Every failure has UNKNOWN feasibility semantics.

The native module and one new regression target are gated by 2014+exact,
independently of Faster/reachability. Existing pinned HiGHS
`7df0786de3088c832297e5ed821db236d8fab281` is used without modification.
All actual free nodes and complete cage handles participate. Sparse triangular
Q is accumulated from row outer products of K weighted by inverse FEM mass;
all prescribed boundary variables are eliminated before optimization. The
reduced Hessian is Q_FF and linear coefficient is Q_FB B_h. Boundary values are
scattered unchanged as part of elimination, never repaired after solving.

`COUPLED_PARTITION` solves all free-node/handle variables in one QP with a
partition equality at each free node. It returns identical raw and final
weights, without clamp or normalization. `INDEPENDENT_NORMALIZED` actually
solves every handle's box QP and independently checks each observed column.
Only all-column raw KKT success permits recorded row normalization. Denominators
must be finite and positive. Raw/final energy, partition/boundary/bounds
residuals, and maximum weight change remain distinct. Final coupled KKT is
always measured; its failure does not falsely invalidate a correctly computed
published independent-normalized variant or inherit raw optimality.

The input-only positive scale is `max diag(Q)` over the original complete FEM.
HiGHS coefficients and KKT use this same scale; geometry and weights are not
scaled. Original FEM `energy()` supplies unscaled scientific energies and
`applyQ()` independently recomputes every assessment gradient, without trusting
the backend Hessian, duals, or `kOptimal`. Coupled multipliers are derived from
interior gradients or valid bound intervals. Boundary nodes are fixed and are
not incorrectly required to have zero gradient. Stationarity, bound-dual sign,
complementarity, primal constraints and projected-gradient mapping are
separately recorded. The unscaled projected mapping is scale times the mapping
with step `1/scale`, documented in the header.

Default acceptance remains primal/boundary/partition `1e-9` and normalized
stationarity/dual/complementarity/projected gradient `1e-8`. Backend primal,
dual, and optimality tolerances are `1e-10`; none was relaxed. Backend settings
are recorded and config-hashed: `qpasm`, one thread, parallel off, seed0,
presolve off, zero QP regularization, small-matrix threshold1e-12, explicit
iteration/nullspace/per-call time limits, and bounded timeless diagnostics.
The time limit is per actual backend call, not a claim of a total wall-time
budget. No global scheduler reset or settings mutation is used.

## Diagnosed backend limitation and approved exact adapter representation

The initial real adapter duplicated HiGHS's preinitialized zero in both sparse
start arrays. Raw diagnostics showed 7 instead of8 constraint entries and
misplaced Hessian diagonals. Clearing those arrays before writing explicit
starts corrected this interface error; its failing logs remain intact.

The corrected zero-general-row independent QP exposed a separate pinned-backend
defect: medium handle4 returned `Optimal` with recomputed normalized KKT
`0.00513806`, so our assessment correctly rejected it. Pinned
`highs/qpsolver/a_quass.cpp` calls `computeStartingPointBounded` when there are
zero general rows and at most15000 variables. The implementation in
`feasibility_bounded.hpp` allocates N² storage, computes only sparse-listed
Cholesky entries without fill, suppresses values below1e-4, and can return its
interior point as optimal without QP iterations. This is incompatible with
the strict contract.

After source inspection, the controller explicitly approved a single redundant
`0=0` row in EVERY independent QP from the outset. It preserves the exact
objective and feasible set while using the general feasible-start/active-set
path. This representation is documented, config-hashed, and included in sparse
resource accounting; it is never an implicit retry/fallback. HiGHS and its
dependency revision were not patched or upgraded. The independent KKT gate
and original FEM equations are unchanged. Both tiny and medium regressions
now require actual backend calls; every tiny call has positive actual
`qp_iteration_count`, with no known optimized weights supplied as initialization.

## Scientific observations on fixed synthetic fixtures

Both mass rules and both protocols solve the supplied `.3–.4` cube with
27 nodes,48 tetrahedra,26 boundary rows,one free node,eight handles. Every nodal
boundary row is compared with independent literal cube-triangle formulas.
The actual optimized free row has `5/16` at000 and111 and `1/16` at the other
six canonical handles. These values are test expectations, never solver input.

The fixed `.3–.7` cube has64 cells,125 nodes,98 boundary rows,27 free nodes and
26 handles: one actual702-variable coupled solve or all26 actual independent
column solves. Observed unscaled energies:

| Mass rule | Protocol | Raw energy | Final energy |
|---|---|---:|---:|
| Restricted-grid Voronoi | Coupled | 986.3219182459500 | 986.3219182459500 |
| Restricted-grid Voronoi | Independent normalized | 986.0962610813783 | 987.1999527812048 |
| Barycentric tet | Coupled | 1003.3219182459500 | 1003.3219182459500 |
| Barycentric tet | Independent normalized | 1003.0962610813783 | 1004.1999527812048 |

For both mass rules, coupled projected-KKT residual is about3.33e-16,
partition residual1.33e-15, and352 free-node weights are bound-active within
1e-9. Independent raw projected-KKT residual is about1.11e-16,310 weights are
bound-active, raw partition residual0.0986206844013331, and maximum explicit
normalization change0.0323777096292423. Final normalized coupled projected KKT
is0.0146489541988485 and `optimalityPassed=false`, while final primal validity
passes. Tests explicitly reject transferring raw optimality to this field.

Both protocols/mass rules also solve two disconnected actual components with
54 total nodes,16 complete handles,two free nodes and zero cross-component
influence within acceptance. Source arrays and FEM coefficients remain
unchanged. Two actual medium runs have identical weights and scientific
digests. Exact dyadic translation/scale preserves weights within1e-9 and
unscaled energy follows inverse uniform scale. Actual solved final fields
transfer by the existing trilinear API to original surface vertices and a new
dyadic vertex with finite bounds/row sums. This remains nodal triangular
boundary data plus trilinear transfer, not continuously identical off-node
triangle functions.

## Resource, trace, and failure contract

Checked multiplication/addition guards cover Q pair contributions, unique
stored sparse entries (including structural zeros), complete coupled variable
indices, dense nodal fields, sparse block replication, starts/constraint
arrays and byte estimates. Q pair counts are checked before allocating Q;
each new unique map node is checked before allocation. No dense N×N Q or dense
root grid is created. Dense retained/work arrays have a conservative16*N*H
value budget, and FEM evaluation's own2*N*H budget is checked. Sparse model
and HiGHS model-copy arrays are included; map-node and diagnostic reserves are
conservative platform estimates. `checkedStorageBytes()` does NOT claim exact
HiGHS internal allocator accounting. All heavy commands run with the specified
24GiB process cap, timeout3600 and sequential `-j3` builds. No runtime
coarsening, handle reduction, changed grid or dropped boundary is used.

Each attempt preserves handle scope, pass-model return, actual run invocation,
run return/model status, actual QP iteration count, candidate availability,
terminal reason, recomputed column/coupled assessment, Unix start/end and
monotonic elapsed seconds. Diagnostic text is capped at65536 bytes with an
explicit truncation flag. Trace scope is backend calls/statuses/counts and
available diagnostic output, not every internal optimizer iterate. Timing and
diagnostic text are excluded from scientific digests; actual scientific event
sequence, controls, statuses/counts, assessments, upstream identities and
raw/final fields are included. No optimizer trajectory was fabricated.

Regression failure paths include actual zero iteration and tiny positive time
limits on the medium problem under both protocols, plus Q-pair/Q-entry,
variable,dense-field,byte and size_t/index-overflow preflight failures. They
never return success or a final normalized field. Arbitrary candidate checks
reject feasible nonoptimal interior fields, wrong pure one-hot active fields,
wrong boundary/partition,bounds violations,NaN,Inf and wrong shape. Independent
raw partition violations can legitimately pass independent KKT and correctly
fail coupled primal assessment. There is no status flag or candidate-loading
API through which an incomplete or invalid normalization can be promoted.

## TDD and verification evidence

All new evidence is retained in fresh unique
`build/task-5k-evidence-CX9SJksa/`. From the worktree root, under
`ulimit -v 25165824`, these commands were used:

```sh
timeout 3600 cmake -S . -B build/strict-v2-spin
timeout 3600 cmake --build build/strict-v2-spin --target spin-it-2014-source-cage-bbw-regression -j3
timeout 3600 build/strict-v2-spin/spin-it-2014-source-cage-bbw-regression
timeout 3600 cmake -S . -B build/task-5c-no-faster-20260909
timeout 3600 cmake --build build/task-5c-no-faster-20260909 --target spin-it-2014-source-cage-bbw-regression spin-it-2014-source-cage-fem-regression spin-it-2014-source-cage-grid-regression spin-it-2014-source-surface-cage-regression -j3
timeout 3600 ctest --test-dir build/task-5c-no-faster-20260909 -R '^spin-it-2014-source-(cage-bbw|cage-fem|cage-grid|surface-cage)$' -V
timeout 3600 cmake --build build/strict-v2-spin -j3
timeout 3600 ctest --test-dir build/strict-v2-spin --output-on-failure
```

1. `red-configure.log`, `red-build.log`, `red-behavior.log`: header/scaffold
   compiled, then all18 behavioral groups failed because actual solving was
   absent. This was the genuine RED before implementation, not just a missing
   include/link error.
2. `green-build-01.log`, `green-behavior-01.log`: first real implementation
   compiled,18 failed groups. `diagnostic-build-02.log` and
   `diagnostic-behavior-02.log` retain the actual malformed sparse-handoff
   diagnostics used to identify the duplicate initial starts.
3. `green-build-03.log`, `green-behavior-03.log`: coupled and tiny/disconnected
   paths passed; six groups exposed the independent bounded-start defect or
   the upstream nondyadic FEM issue. Both failures are preserved.
4. `green-build-04.log`, `green-behavior-04.log`: all18 groups passed after the
   approved exact redundant-row representation and scoped dyadic covariance
   fixture. No acceptance threshold changed.
5. `trace-red-build-05.log`, `trace-red-behavior-05.log`: new actual-run trace
   assertions failed17 groups before invocation metadata was populated.
   The same log also preserves an explicit rerun of the known upstream
   nondyadic failure. That known-bug assertion was not retained as a permanent
   regression requiring a bug to persist.
6. `trace-green-build-06.log`, `trace-green-behavior-06.log`: final trace/resource
   implementation, all18 groups passed; no compiler warnings in new code.
7. `faster-off-configure.log`, `faster-off-build.log`, `faster-off-ctest.log`:
   FasterOFF compatibility passed4/4 tests,6.74s CTest elapsed. Both cached
   feature flags2014/exact are ON and Faster is OFF in that build.
   The all-target main build also completed successfully in54.412s before the
   sole full CTest invocation.
8. `all-target-build.log` and `full-ctest.log`: both commands exited0. All203
   actual registered tests passed,0 failed,499.12s CTest elapsed. The raw log
   contains203 separate passed result rows and no skipped/disabled tests.

Controller-owned independent final evidence is under
`build/task-5k-controller-Nmogqv/`: all8 tiny/medium mass/protocol configurations
were run twice from a final statically linked native probe, compared against
an independently assembled geometry/literal-stencil reference and polished
original KKT equations. All passed; maximum weight difference1.42e-15 and
independent raw stationarity2.98e-16. Seven deliberately corrupted candidates
were rejected. `probe-final-bindings.sha256` binds the final header/source/test,
core/probe/reference/checker artifacts. These are controller results, not a
replacement for the native regressions or an original-paper claim.

The controller also ran an actual partial-column failure probe under BOTH
mass rules with fixed medium independent iterationLimit30. Column0 passed at
26 actual iterations; column1 returned IterationLimit at31 actual iterations.
The pinned backend can overshoot that requested limit by one, and the trace
retains31 rather than falsifying the count. Both runs expose exactly2 attempts,
unavailable complete raw assessment, zero final values, and UNKNOWN
feasibility. Evidence: `partial-failure-final.log` and `partial-failure-probe`
in the controller directory; source
`.superpowers/sdd/2026-09-05-spin-it-strict-computational-v2/task-5k-partial-failure.cpp`.
This additional check did not alter the frozen production or registered tests.

Existing configure warnings are preserved: libigl's old CMake minimum-version
deprecation warning in the main configure, and existing dependency CMake
minimum-version warnings in FasterOFF. No broad CMake cleanup was attempted.
The new regression adds one CTest registration to baseline202; no existing
tests were removed or disabled.

## Self-review and remaining concern

Self-review traced the exact original-FEM objective through boundary elimination,
triangular Hessian indexing, replicated constraints, backend status/candidate
checks, per-column completion, row normalization, final coupled KKT and immutable
digest construction. False-acceptance tests exercise both all-interior and
all-at-bound candidates. Actual backend `Optimal` with incorrect KKT was
observed and rejected, demonstrating why status-only acceptance is insufficient.
The independent controller owns geometry-derived reference checking and later
fresh review; its code/docs are not included in this implementation commit.

**Open upstream FEM arithmetic limitation:** constructing the valid medium
cube after scale1.7 and origin(-.3,-.2,-.1) yields a valid125-node source grid
but `SourceCageFem` throws `source cage FEM incomplete boundary constraints`.
This is not mislabeled invalid geometry. Reproducer: use source cube(.3,.7),
multiply its vertices by1.7 and add that origin, configure cageLevel2/rootSide1.7
and matching rootOrigin, fineLevel3, then construct `SourceCageFem`. The failed
original transform is preserved in logs03 and05; no upstream code was changed.
The controller explicitly accepted an exact dyadic positive covariance test
(scale2, origin(-.5,-.25,-.125)) in this bounded increment and required this
limitation to remain open for a separate FEM audit. This prevents a blanket
claim of arbitrary floating-coordinate covariance even though the BBW adapter
passes all supported fixed fixtures.

## Completion addendum

Status: **DONE_WITH_CONCERNS** for the bounded Task5k increment. The open
nondyadic upstream FEM limitation and unrecovered original-author optimizer/
internal trajectory remain explicit. Full Task5/M0–M8 is not declared complete.

Scoped implementation commit:
`6ba2aa1f8e693694b98fd87c7d3f9de1a0917a07`
(`feat(spin2014): solve source-cage BBW QPs with independent FEM KKT`). Its only
four paths are the new header/source/native regression and six CMake additions.
This report is a separate scoped documentation commit. The independent
controller verification document is intentionally left for its owner.

Full final raw logs are bound in
`build/task-5k-evidence-CX9SJksa/final-logs.sha256`.
`final-source-binary-dependency.sha256` binds all11 source/binary/dependency/
cache artifacts and was rechecked after the full suite; every entry matched.
Key SHA-256 values:

```text
5feb7e19e4c497a594f9609709cd974b5f119e35447484a380570e2ad2010814  include/spin2014/SourceCageBbw.h
fa8f7cf7d5ee6905af2fb605b9eef9e03157d19c91d952d3d9fbdb31dfbaeed5  src/spin2014/SourceCageBbw.cpp
73857bdde619beffbe448cf87dbee9f65ba160773e82ea07dcbf1ad3d7d5cbfc  tests/spin_it_2014_source_cage_bbw_regression.cpp
c2fce7457ba0e060848efeeb8bdcb7e5b3402a53c0efa64c14fc8840052389e9  CMakeLists.txt
72cb65c0105338299f6058dad00b39dea4b8debb42239d9841a455657955006e  build/strict-v2-spin/libmake_it_stand_core.a
98126ba3cc7c777e011bac7309ecd351007cb91453d0e695eebde48d1afcfc3a  build/strict-v2-spin/spin-it-2014-source-cage-bbw-regression
718f57b3abf97e5f5b1165198c7c01c904cf24b87ca5ca35a016526d3ac35cc5  build/task-5c-no-faster-20260909/spin-it-2014-source-cage-bbw-regression
def6e678af29a2dae4c980c334f2b53f77664814c52f60fac8622b01381ec110  build/strict-v2-spin/lib/libhighs.a
```

Final digest examples for the regression's medium fixture are retained in
`trace-green-behavior-06.log`: coupled Voronoi
`9319cffd4ba3224b49b311c635496c1c63f2e1f099aa0d52bf970e753a4d541e`,
independent Voronoi
`a73cac5126cd74696b2159546c82e1a48bb056f57cadcb453300dcf5b89831d9`,
coupled barycentric
`50517eae88e8078f29f2bdce72fc00ffd96047ffb57460050565e33fe941be3d`,
independent barycentric
`981c0905b2c78a70c84a4d7533494bfa552e4e99e591ad60c605c7caa2ee3bc2`.
These are config/upstream-bound computational digests, not geometry seals.

Unrelated worktree edits appeared DURING the full CTest run at18:12–18:14:
`include/spin/QuadricMomentEvaluator.h`, `src/spin/QuadricMomentEvaluator.cpp`,
`tests/quadric_moment_evaluator_regression.cpp`, and new
`tests/quadric_certified_empty_cases.h`. Neither this implementer nor the
controller made those edits. They are preserved and excluded from both scoped
commits. The203/203 result applies to the exact built binary/source checkpoint
and hashes above, **not** those newly edited/unbuilt Faster files or a claim
that the entire current worktree is clean/tested. No rebuild incorporating
that unrelated work was attempted. The controller explicitly confirmed this
scoping before the implementation commit. Fresh independent review follows
this report; this implementer did not spawn a helper or reviewer.
