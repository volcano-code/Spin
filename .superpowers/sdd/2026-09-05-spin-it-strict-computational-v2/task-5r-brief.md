# Task5r: actual-source joint paper objective and constraints

Approved strict plan Task5/M3, following completed Task5q. BASE86234ff.
Canonical worktree `/data/xvxinyao/moment-guided-spin-v1/code/make-it-balance-paper-sim`,
branch `spin-it-strict-computational-v2`. Main checkout is read-only.

This bounded integration connects the actual source/full-BBW moment functional
to five existing paper Givens objectives, a fixed axis/support coordinate frame,
all joint first derivatives, and existing sparse source-fill regularization.
It is NOT the complete shape regularizer, joint SLQP, or geometry-commit gate.
No reconstructed asset or this synthetic regression grants original-paper success.

## Owned files and interfaces

Create `include/spin2014/SourceCagePaperFunctional.h`,
`src/spin2014/SourceCagePaperFunctional.cpp`,
`tests/spin_it_2014_source_cage_paper_regression.cpp`,
`docs/spin2014_source_cage_paper_functional.md`, scoped CMake source/header/test lines.
Do not edit existing source moments, Givens formulas, legacy deformation,
solvers, source partition/BBW or validators. Propose exact public API promptly.
The new immutable functional consumes an actual SourceCageMomentFunctional
snapshot, config with objective and independent inertia/center/fill weights,
and explicit fixed reference frame. No supplied replacement mass or Jacobian.

Existing building blocks:
- SourceCageMomentFunctional: evaluate/evaluateValues(beta, all handle parameters),
  actual world-coordinate RawMassMoments10, beta/handle Jacobians, moved points,
  input/config/output SHA identities. betaCount B, handleCount H, handle order
  [scale,tx,ty,tz], all fixed rotations identity, strictly positive scales.
- SourceConformingBoundary.binding().domain(): immutable actual partition and
  sparse uniformLaplacian()/edges(). Own a snapshot; no external mutations.
- Spin2014GivensProblem: five objectives, eleven input coordinates (10 moments
  plus phi), analytic objective/constraint gradients. Reuse this mathematics,
  not the legacy angle-eliminated evaluator. No new paper-formula copies.

Joint state layout is [beta B, phi only for non-static objective, handle4H].
Static has NO redundant angle variable. Split-state evaluate must reject nonzero
static phi and nonfinite phi. Expose counts/offsets, settings, identity and an
initial-state builder from explicit beta/phi plus actual initial handles.
State validation rejects wrong lengths/nonfinite/beta outside[0,1]/nonpositive
scale; no silently clamping or dropping handles. Static count B+4H, others B+1+4H.

Output retains actual world moment evaluation (including moved points),
spin-frame moments, unregularized paper objective, fill penalty, total objective,
constraints, joint gradient, joint constraint Jacobian, and spin-frame raw
moment Jacobian with exact zero phi column. A value-only path returns identical
moments/geometry/objective/constraints but no global derivative arrays, and can
work when a derivative cap prevents the derivative path. No full joint Hessian
claim: moment second derivatives are not available; later solver uses L-BFGS.
Output identity binds actual upstream evaluator output, full state, objective,
all weights, frame, budgets, mode, results; no timestamp or pass boolean.

## Fixed frame and exact transformation

Given origin o in source/world coordinates and a proper orthonormal 3x3 matrix R,
spin coordinates are q=R*(p-o). +z is the chosen axis; origin is the selected
support/axis point. Explicit defaults o=0,R=I are implementation defaults, not
recovered paper poses. Validate finite entries, R*R^T=I and detR=+1 within1e-10;
reject scaling/shear/reflection. No recentering by the current center of mass.

For world density-weighted (M,L,Q):
M'=M, L'=R*(L-M*o),
Q'=R*(Q-o*L^T-L*o^T+M*o*o^T)*R^T.
Build/use the exact linear10x10 map in native order
[M,LX,LY,LZ,QXX,QYY,QZZ,QYZ,QXZ,QXY]. No density multiplier here.
Apply same constant map to actual source Jacobians. No finite differences in
production. Never confuse tensor off-diagonals and doubled polynomial terms.
Phi is the existing transverse Givens angle, independent of fixed R.

## Paper objective, constraints and sparse fill term

Use Spin2014GivensProblem with transformed moments and actual phi. Preserve
all existing five objective distinctions and mode constraints:
- static: LX,LY (2);
- top: LX,LY,QXZ,QYZ,transverse diagonalization (5);
- yo-yo additionally LZ (6, existing ordering).
Don't normalize constraints with result-dependent scales or replace equalities
with success flags. Chain all moment and phi derivatives to the joint state.

Fill penalty uses actual fixed-rest source partition L:
Efill=0.5*gammaL*beta^T*L*beta,
gradient=gammaL*L*beta, handle and phi columns zero.
Independent edge form is0.5*gammaL*sum_edges(beta_i-beta_j)^2.
Keep sparse L, no dense BxB. Explicitly label this source unit-adjacency
discretization IMPLEMENTATION_DEFAULT, not the full recovered paper Laplacian.
No shape/displacement energy field implicitly claiming those terms complete.
When partition/BBW changes a NEW functional must be constructed. No adaptive
split/rebind integrated in this task; document requirement clearly.

## Resource and arithmetic contracts

Finite nonnegative weights (existing objective validation reused); finite proper
frame. Positive explicit capacity for all new joint derivative allocations,
proposed maximumJointDerivativeValues=10000000. Check index/size products and
all retained/workspace derivative array capacities before allocation; keep
source moments' own budgets. No P*H/P*4H or dense NxN arrays. O((10+K)*N)
joint Jacobians are permitted, plus bounded source scratch. Value-only must not
allocate joint matrices or require their capacity. Invalid config/state=>
invalid_argument; bad capacities/allocation=>length_error RESOURCE_LIMIT;
nonfinite arithmetic=>fail closed, no NaN result. Do not alter legacy exception
contracts. No inferred feasibility from finite objective/positive mass.

## Test-first acceptance

1. Register real missing-interface test, observe expected RED before production.
2. Actual two source/BBW pipeline fixtures from Task5q native regression:
   baseline and transformed nonsymmetric, all26handles,8beta; sourceD5/cage2/grid3.
   Keep nonuniform beta and distinct all-handle scale/translations. All cases
   SYNTHETIC_REGRESSION only. Reuse fixture utilities, no fake supplied weights.
3. Five objectives, nonzero phi for nonstatic, nonzero origin and nontrivial
   proper rotation. Independently transform moved POINTS, integrate signed tetra
   moments, and compute objectives/constraints without production Givens calls.
   Test both world and spin moments, mixed entries and inertia conventions.
4. Check objective gradient, constraint Jacobian and transformed raw Jacobian
   for ALL113 nonstatic/112static columns; central FD at three fixed steps
   e.g3e-4,1e-4,3e-5 with fixed input-derived scales. Nonzero derivative relative
   target1e-6; predeclare near-zero absolute rule, don't tune after failure.
   Supplement with independent complex-step/direct derivatives where necessary
   to resolve near-zero entries; no changing fixture or tolerances to pass.
5. Fill energy/gradient from actual edge differences, not L under test. Constant
   beta cancels; weights independent; density covariance: inertia and pure height
   penalties invariant, mass-height scales density squared, fill unchanged.
6. Rigid covariance: rotating/translating world geometry and corresponding frame
   preserves spin moments and objective. At minimum directly test exact frame
   map against independent transformed geometry; no result recentering. Static
   no angle column; other raw moment phi columns exactly zero. Phi affects actual
   constraints. Value-only arrays absent and same results; derivative cap guard.
7. Wrong states, invalid objectives/weights, nonfinite/reflective/nonorthogonal
   frames, overflowing origin, budget and invalid arithmetic reject. Frame,
   state/config/weights/mode identity sensitivity plus identical-run identity.
   Don't require geometric feasibility for a value/derivative fixture.
8. Existing source-moment/Givens/fill/legacy tests remain valid. Full assigned
   build relink then safe scoped suite once; actual exact-OFF graph excludes
   new functional. No full global-suite run across unrelated dependency scripts.
9. Self-review, scoped commit, detailed report including raw RED/GREEN/logs and
   hashes. Fresh task reviewer is controller's responsibility.

## Execution safeguards

Only existing canonical build `build/task-5c-no-faster-20260909`, downloads
disabled (`LIBIGL_SKIP_DOWNLOAD=ON`), FasterOFF, exactON,2014ON. max-j3,
24GiB address-space limit, serialheavy tests,60minute batch. No manual cache
deletion/repair, no downloads/pin changes. ExactOFF test children already disable
downloads. CMake hunks must be solely new source/test registration.
Safe suite regex `^(spin-it-2014-|spin-2014-|source-spin-material-python$)`.
Foreign Faster/Plane/Polyhedral/FullRank/Truncation files and experiments actively
change: never edit/stage them. Main tracked-diff SHA93e32af47ab1957a1efcb6716437656be129053f4269d2169378fcd95b05a35e.
If source freshness test detects external edits, preserve failure and inspect
begin/end evidence; at most one separately recorded quiet retry, never disable
validator or rewrite a failed suite as all-green.
No subagents, merge/push/tag/release/deletion. Commit explicit ownedpaths using
authorized volcano <xuxinyao307@gmail.com>. Controller independently verifies
exported data; tell it exact API and when source/core stable for linking.
Report `.superpowers/sdd/2026-09-05-spin-it-strict-computational-v2/task-5r-report.md`,
return <=15lines status/commit/tests/concerns. Global PARTIAL_BLOCKED unchanged.

## Numerical verification addendum

On the unchanged original fixtures and three preselected step sizes, ordinary
double-precision two-point differences showed truncation and subtractive-noise
limitations, including finest-step relative errors above1e-6 for very small
derivatives. Their original failed logs and all sample counts remain evidence.
Both native and controller independent complex-step point/tetra/inertia paths
agree with production analytic derivatives at relative errors below1.3e-10.
No production formula, fixture, weight, step or tolerance is changed to fit a
desired result. Ordinary production-double differences are retained diagnostics,
not relabelled as passing1e-6 checks.

Ruling: a genuine long-double independent point/tetra/inertia reference and
five-point centered stencil at the SAME three fixed steps are permitted to
resolve this verification problem. Precision must persist throughout the
independent evaluation, not be cast onto double results. Require actual>=64bit
mantissa support and record it. Production analytic derivatives still must
meet1e-6 relative and the predeclared native near-zero absolute2e-9 for
|derivative|<=1e-7, with no enlarged cutoff; controller separately checks the
stronger near-zero1e-11. All three high-precision stencils, independent
complex-step and original double FD diagnostics must be recorded. If the
high-precision check fails, keep it failed and investigate; no success grant
from this addendum itself. This improves checker accuracy rather than loosens
the derivative requirement. Cost if wrong: checker/algorithm correction, not a
paper-experiment exception.
# Fix-round1 resource-regression addendum

Independent review found that the two fixed10x10-to-dynamic Jacobian products
invoke Eigen GEMM packing outside the declared derivative budget. The original
worker owns this repair in the same production/test/doc files and may extend
only Task5r's CMake registration with one dedicated allocation-probe test.
Use the actual production column-transform helper with preallocated10x8 and
10x104 inputs/outputs. A target-local Eigen runtime no-malloc guard, with stack
packing forced through its checked allocator, must reject the old GEMM route
and accept the bounded fixed-ten-vector route. Keep probe macro settings local;
do not alter core/global Eigen flags or rely on incompatible runtime state
across translation units. Normal functional tests continue to link actual core.
Register the probe under the same2014+exact gate and include it in isolation
verification. This expands only resource-test coverage, not scientific scope,
formula, tolerance, input or the declared cap. Preserve all pre-repair logs and
manifest; bind repair evidence to new hashes and run fresh scoped tests.
