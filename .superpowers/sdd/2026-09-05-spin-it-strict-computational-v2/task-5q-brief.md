# Task5q — analytic source-boundary / full-BBW material moments

Approved strict plan Task5/M3 derivative integration following completed Task5p.
BASE6f61946. Worktree `/data/xvxinyao/make-it-balance-paper-sim`, branch
spin-it-strict-computational-v2. Global PARTIAL_BLOCKED. This is the source
material derivative layer required before the full joint solver, not an attempt
to relabel the legacy dense void-cell adapter as source integration.

## Interface and owned files

Create `include/spin2014/SourceCageMomentFunctional.h`,
`src/spin2014/SourceCageMomentFunctional.cpp`,
`tests/spin_it_2014_source_cage_moments_regression.cpp`, and
`docs/spin2014_source_cage_moments.md`; scoped CMake source/header/test lines.
To reuse supplement Algorithms1/2 without duplicate formulas or dense global
vertex Jacobians, refactor `src/spin2014/SurfaceRawMomentDifferential.cpp` into
a private reusable per-triangle polynomial kernel in
`src/spin2014/SurfaceTriangleMoments.h` (and .cpp if needed). Keep public existing
surface APIs and behavior compatible. No source-tree/cage/binding/solver changes.
Propose exact API early; controller will compile an independent exporter.

Immutable SourceCageMomentFunctional consumes actual SourceConformingBoundary
and positive finite density, with explicit checked config. No supplied-weight,
mesh, precomputed fixed-shell moment or passed-boolean constructor. Retain
boundary/config/counts/upstream SHA. Evaluation accepts beta (one per actual
adaptive leaf, hollowed fraction in [0,1]) and all H handle parameters. Define
and document one deterministic layout [scale,tx,ty,tz] for each actual cage
handle in original order. Initial scales1/translations0. Scales must be finite
and strictly positive; translations finite. All handles are retained. Fixed
rotations are identity for this source interface (explicit assumption, not
an unimplemented author rotation claim).

Output: density-weighted RawMassMoments10 finalMoments, 10xB betaJacobian,
10x(4H) handleJacobian, deformed original+material point coordinates, and
versioned output identity binding input state/config/upstream and arrays.
Provide a value-only evaluation path that uses the same polynomial values but
does not allocate Jacobians; useful for solver line search and independent FD.
It must return the same values/geometry as derivative evaluation and distinguish
its absent Jacobians explicitly. Exact signature decision may be proposed.

## Mathematics and sparse evaluation

For rest point p and fixed actual handle position ch:
  p' = sum_h w_h(p) [ ch + scale_h*(p-ch) + translation_h ].
Weights come from the ACTUAL successful final BBW field via rest-space
SourceCageGrid embeddings, including new face centers/hanging points. Do not
re-embed moved points, renormalize weights, subsample handles, or allocate global
point-by-handle / point-by-parameter matrices. Cache bounded eight-node
embeddings if helpful; use bounded chunks or local point weight rows.
Point derivatives: dp'/dscale_h=w_h*(p-ch), dp'/dtranslation_h=w_h*I.

Actual material moments:
  s = density * (s_outer(p') - sum_primitive beta[owner] s_primitive(p')).
All10 components use native order M,LX,LY,LZ,QXX,QYY,QZZ,QYZ,QXZ,QXY.
Outer surface is the binding source's original faces at moved outer rows.
Primitive surfaces are the actual conforming boundary, not leaf bounding boxes.
For a canonical shared triangle, signed cavity coefficient is
  (positiveOwner present ? beta[pos.betaIndex] : 0)
  - (negativeOwner present ? beta[neg.betaIndex] : 0).
One triangle evaluation supports both owners with opposite signs. Beta
derivatives use the corresponding +/- triangle contribution, handle derivatives
contract local triangle coordinate derivatives and full fixed-rest weights.
Accumulate via reusable supplement surface-polynomial kernel, not finite
differences or the independent tetrahedral checker as production implementation.
Original whole-surface integrators still validate closed chains; the private
triangle kernel makes no closure claim. Every deformed triangle is checked for
finite/nondegenerate geometry even if beta coefficient is zero. Invalid arithmetic
or nonpositive final mass rejects, but positive mass is NOT geometric legality.
Self-intersection, inversion relative to rest, overlap, shell thickness and
final exact partition are separate precommit gates, explicitly not granted here.

## Resources / identity / compatibility

Suggested config maximumWeightValues=1000000, maximumJacobianValues=10000000;
positive finite density. Bound all weights, Jacobians, embeddings, points,
indices and derived allocation sizes before allocating; respect upstream
applicable budgets. No silent reduction/partial output. Invalid state/density:
invalid_argument. Nonfinite arithmetic: overflow_error. Invalid/overflowing
budgets or allocation failure: length_error with RESOURCE_LIMIT. Large values
must not wrap Eigen or size_t arithmetic. Empty carveable material is valid:
no primitive contributions; any actual beta leaf count remains domain-defined.
Value-only path must not require Jacobian capacity or allocate dense derivatives.
Hash version/method/layout/density/budgets, upstream actual boundary/BBW, beta,
handle parameters and produced scientific arrays deterministically; timing is
excluded. State hashes differ for distinct states even when geometric motion
cancels. Caller cannot mutate retained inputs or obtain an upgraded certificate.
2014 AND exact required, Faster independent. Existing legacy APIs remain.

## Test-first execution and fixed checks

- [ ] Missing-interface RED from newly registered actual native target before
  creating production interface. Preserve log; not a typo failure.
- [ ] Actual baseline cube outer(.3,.7), design(.34,.66), sourceD5, adaptive
  initial2/max5, cage2/grid3/actual coupledBBW, all26handles,8beta,496primitives.
  Rest beta0 => analytic outer moments. Beta1 => outer minus accepted interior
  [11/32,21/32]^3. Nonuniform beta uses independently integrated real primitive
  boundaries and ownership; no bounding-cell approximation.
- [ ] Full-BBW distinct handle scales/translations; transformed nonsymmetric
  fixture. Independent signed-tetrahedral final10moments match1e-10 normalized.
  Finite difference ALL beta/handle columns at THREE preselected steps, e.g.
  3e-4,1e-4,3e-5 times appropriate fixed input-derived scales. For nonzero
  entries check relative1e-6 with explicit near-zero absolute rule, not a
  denominator chosen after seeing failures. Truncation-dominated errors must
  show reasonable h^2 improvement; analytically linear beta/translation entries
  may already be at roundoff and must be recorded rather than forced monotone.
- [ ] Independently differentiate signed-tetra formulas on a fixture or use
  complex-step test oracle to verify near-zero derivative entries without
  sharing the production kernel. Check mixed moments/non-diagonal inertia.
- [ ] Density scaling, joint uniform translation/scale covariance, beta signs,
  constant beta shared-face cancellation, full26handles and rest fixed weights.
  Value-only matches evaluate and works with a Jacobian budget too small for
  derivative output. No global P*H storage; resource caps/exhaustion tested on
  actual pipeline, not grep-only tests.
- [ ] Invalid beta size/bounds/NaN; invalid handle size/scale/NaN/Inf; density
  invalid; malformed/overflow budgets; finite input that creates nonfinite or
  collapsed output rejects. Empty carveable binding if actual upstream permits.
  Deterministic two identical evaluations and state/config-sensitive identities.
- [ ] Existing surface integration/differential regressions pass after local
  kernel factoring, including value-only allocation guard. Focused GREEN then
  safe suite once using assigned Faster-OFF build; actual retained exactOFF
  dependency graph excludes new functional.
- [ ] Self-review, scoped commit, report with commands/RED/GREEN/SHA/limitations.

Do NOT weaken fixtures or thresholds. If the intended FD steps leave noise,
retain results and ask controller with concrete evidence, do not claim failure
as infeasibility. Sparse sourceD9 geometry is optional this task: no fulllevel9
optimization claim from any fixture. No original-paper case counted.

## Ownership and reporting

Sole production/build writer; max-j3, ulimit-v25165824, heavy tests serial,
timeout3600. Only build/task-5c-no-faster-20260909. Do not configure/build the
main strict-v2-spin build across external Faster edits. Safe final command:
ctest --test-dir build/task-5c-no-faster-20260909 -R
'^(spin-it-2014-|spin-2014-|source-spin-material-python$)' --output-on-failure -j1.
Other dependency test scripts are not audited for scratch safety: no broad suite.
Tell controller when stable core archive is ready for independent linking.
Preserve all external Faster/PlaneMoment files, tests and experiments; no staging
their files or CMake hunks. Check CMake diff ownership before build/commit.
No subagents. Commit explicit owned paths only, authorized volcano
<xuxinyao307@gmail.com>. No push/merge/tag/delete. Full report path:
.superpowers/sdd/2026-09-05-spin-it-strict-computational-v2/task-5q-report.md.
Return under15lines status/commits/testsummary/concerns/buildownership.

Ruling: fixed-rest full-BBW affine motion and analytic surface moment chain
are a separate necessary layer; not full paper regularizer/SLQP/geometry seal.
Cost if wrong: source functional/derivative rework, not relaxed reproduction
criteria. The complete approved user spec still governs later gates.

## Evidence-driven infrastructure addendum (external relocation)

The external move invalidated download metadata and owned build regeneration
re-created pinned CGAL (and may have re-created Eigen); no previous dependency
directory backup is claimed. Parent build is now verified SKIP_DOWNLOAD=ON.
Fresh dependency-contract children still default OFF. Before final suite, also
own ONLY the following focused behavior change:

- `tests/cmake/assert_spin_it_2014_dependency_isolation.cmake`: in generated
  SUCCESS_SOURCE/CMakeLists.txt, add
  `option(LIBIGL_SKIP_DOWNLOAD "" OFF)` and fail configuration if it is OFF.
  Run existing isolation test with the unmodified dependency script; expect
  RED because actual child invocation has not disabled download.
- `tests/cmake/assert_spin_it_2014_dependency.cmake`: add
  `-DLIBIGL_SKIP_DOWNLOAD=ON` to both configure child invocations. All feature,
  retained-scratch, build and runtime assertions stay unchanged. Run isolation
  GREEN, then actual final safe suite; inspect both fresh child caches ON.

This is test harness isolation, not dependency availability being skipped: a
required unavailable library must fail, never download/replace a source tree
inside a registered regression. Keep all RED/GREEN/previous interrupted logs.
No other dependency scripts or global defaults are changed.
