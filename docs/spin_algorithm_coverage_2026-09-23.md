# Spin-It / Spin-It Faster algorithm coverage baseline (2026-09-23)

This is an evidence index, not a completion percentage. A row is marked
`implemented_and_regressed` only when the current registered suite exercises
the path; it is not a claim that the corresponding author experiment has been
reproduced.

| Area | Current status | Evidence | Remaining strict gap |
|---|---|---|---|
| Spin-It mass moments and non-diagonal inertia | `implemented_and_regressed` | Spin-It 2014 label and current full suite: 278/278 | Author mesh/units/density and paper numerical comparison are missing for most cases |
| Spin-It five objective/constraint modes | `implemented_and_regressed` | `spin-it-2014-objectives`, `spin-it-2014-modes`, `spin-it-2014-givens` | Need author parameter table and full case matrix |
| Spin-It SLQP, active set, restoration and KKT | `implemented_and_regressed` | source-cage guarded SLQP and Spin-It 2014 tests | Need paper-scale runs and published iteration/time comparison |
| Spin-It adaptive source-surface octree | `implemented_and_regressed` | source-adaptive, geometry9 and surface-octree-level9 tests | Unit-level level 9 is not a complete paper-scale level-9 experiment |
| Spin-It cage/BBW and joint deformation | `implemented_and_regressed` (geometry guard only; joint terminal remains open) | source-cage BBW/FEM/grid and deformation groups; exact primitive-pair overlap proofs are reused only for identical region digests (2026-09-26) | Need nontrivial adaptive-to-binary joint output and independent final-material validation; no certified paper asset/case output |
| Spin-It dual material | `implemented_and_regressed` | source-cage-binary-material and dual-material tests | Original Sumone inputs/densities and paper output comparison are not available |
| Spin-It reachability classification | `implemented_and_regressed` | source affine certificate plus reconstructed Heart/Cuboid probes; non-principal Cuboid positive case has two identical sealed runs | Certificate is for a discrete affine partition only; continuous infeasibility remains `UNKNOWN`; positive case is reconstructed and level 4, not a paper level-9 case |
| Faster quadric coefficients and inverse Gauss map | `implemented_and_regressed` | Faster label, 39/39; full suite | Need paper input domains and numerical table/figure reproduction |
| Faster Boolean and discrete moment differential | `implemented_and_regressed` | registered polyhedral derivative, cross-topology Newton guard, the 2026-09-23 40-test Faster suite, and the 2026-09-24 reduced-stratum evaluator regressions | Full-rank ellipsoid/hyperboloid paths use all ten coefficients; plane, simple elliptic/hyperbolic/parabolic cylinders, and simple-spectrum paraboloid return true-dimension fixed-stratum Jacobians through the native evaluator. Topology transitions fail closed; rank-deficient optimization, other degenerate families, and original paper inputs remain open |
| Faster Newton, dual coordinates, restarts | `implemented_and_regressed` | Faster canonical/Newton/end-to-end tests; explicit nontrivial `DiscreteFullRankQuadric` Boolean-to-Newton solve converges and repeats identically on a reconstructed hyperbolic yo-yo fixture (2026-09-24); strict regressions reject nonpositive initial multipliers and restart at a synthetic zero-residual/nonpositive-multiplier stationary candidate (2026-09-26); per-restart event ledger and budget validator | Need original-input 50-seed runs, broad start/topology coverage, and independent repeated full failure/restart statistics |
| Faster 2D polygon--disk path | `implemented_and_regressed` (partial application chain) | Dogleg returns exact oriented line/circular-arc cycles; exact extrusion result packages the design sweep and per-density fixed polygonal prisms with independent ten-moment records; a real layered dogleg result now passes through pinned OCCT to multi-solid B-rep/STEP and reloads with all moments within `5e-8`; wrong density metadata is rejected (2026-09-26) | Synthetic method fixture only; no author fidget input/parameter protocol; arbitrary layer arrangements and material provenance in downstream CAD remain limited |
| Faster static-stability plane objectives | `implemented_and_regressed` (method-level paths) | `spin-faster-static-plane-balance`; analytic fixtures exercise maximum mass, lower center of mass, Pareto `kappa=0` and `kappa=1` endpoints, and an interior Pareto KKT solve with independent balance and ten-moment checks (refreshed 2026-09-26) | Dedicated §4.2 plane-family solvers only; they are not the general 10-D Eq. 14 solver, and there is no author vase input/parameter protocol |
| Faster hydrostatics / floating applications | `implemented_and_regressed` (synthetic two-stage method chain) | Measured A-side prescribed-waterline solve and material seal; B consumes A's frozen geometry/moments/digest, runs its own plane solve, and seals the result; pose-stability adapter recomputes hydrostatic geometry after pose changes; bridge suite passed 24/24 on 2026-09-26 | Author Floating Island geometry, pose/protocol and paper-specific parameters remain absent; synthetic A→B is not the paper experiment or a broad perturbation matrix |
| Faster CAD output | `implemented_and_regressed` (partial analytic-quadric and layered-2D method paths) | pinned OCCT 7.9.3 round-trips coefficient-driven spheres, rotated ellipsoids, rotated elliptic cylinders/paraboloids, one-/two-sheet hyperboloids, hyperbolic paraboloids, exact layered 2D solver output, and (2026-09-27) a sealed reduced-stratum plane candidate through exact half-space Boolean, multi-material B-rep/STEP, and independent ten-moment validation; plane fixture contains one editable plus six fixed solids at distinct densities and passes the `1e-6` material-moment gate; full Faster material bridge passes 26/26 | STEP does not yet embed native material labels/densities; Cherchi face/source provenance, general multi-shell author material mapping, broader solver-result coverage, remaining degenerate quadric families, the full-rank triangulated-to-analytic CAD mismatch, and author CAD/design-domain comparison remain open |
| Spin-It original case matrix | `blocked_missing_asset` | model acquisition audit and Armadillo domain diagnostic; reconstructed Heart and non-principal Cuboid now provide method-validation cases | Heart/Elephant/Woven Ring/Armadillo/T-Rex author assets, poses and critical parameters |
| Faster Table 1 and Figures 10--13 | `not_started_strict` | only miniature paper tests are registered | Need all original inputs, 800 top-level runs, Pareto scans and degeneracy protocols |
| Physical printing and measured experiments | `not_in_scope` | project scope | Must remain separate from computational reproduction |

## Current evidence identity

- Build: `build/strict-algorithms-exact-current-20260922`.
- Full registered run: `build/registered-ctest-evidence-20260923-cuboid-target-axis-topology-contract/ctest.log`.
- Full run result: **278/278 passed**, exit code 0, wall time `869.73 s`; log SHA-256 is `26c3c7df9ad820d1ed833b455fa440836e7fd7deca5c4ba813636a60bc3b6d10`.
- This run includes the source-adaptive target-frame runner contract and both Faster topology-strata outcomes: decreasing transitions are relinearized on their new derivative map, while non-decreasing transitions are rejected. It is a test-suite receipt, not a replacement for the per-case two-process scientific-output checks.
- Project-level status: `PARTIAL_BLOCKED`; no `author_bitwise_exact` or
  `VERIFIED_COMPUTATIONAL_REPRODUCTION` claim is made.

## Fresh 99-row strict baseline — 2026-09-23

The fixed catalogue currently contains **99 requirements** (40 original
experiments, 28 methods, 11 formula rows, 6 algorithm rows, 4 theorem rows, 7
evidence-chain rows, 1 independent-check row, and 2 physical rows explicitly
out of scope). A new fail-closed audit was generated at
`build/strict-current-baseline-audit-20260923/audit.json` from the current
catalogue. It reports `PARTIAL_BLOCKED`, `required=99`,
`blocked_or_not_implemented=97`, `physical_measurements_out_of_scope=2`, and
`verified=0` because no stale checker result was admitted.

The previous checker bundle was deliberately rejected with
`STALE_RESULT: cmake_configuration content binding changed`; its historical
49/98 figure remains historical evidence only. This is a freshness correction,
not a regression claim. The next evidence task is to rerun the affected
registered checkers against the current binary and rebuild their bindings.

## Current Faster checker refresh — 2026-09-23

Four additional checker families were rerun twice against the current binary:
discrete moment/inverse-Gauss fixtures, Equation-14 Newton/restart fixtures,
the published Section 6.6 degenerate examples, and the pinned Cherchi/CGAL
Boolean cross-check. All evidence has freshness status `CURRENT`.

The resulting audit is
`build/strict-current-baseline-audit-20260923/audit-after-faster-four-20260923.json`
(SHA-256
`de2faf90ccca02851815fbf2e12f1367132bae2dfc7fed58639e7f4c360a7bff`). It
contains **18 independently verified rows**, **79 not verified rows**, and 2
physical rows explicitly out of scope. No original paper-case or Table 1
credit is inferred.

## Faster application-contract refresh — 2026-09-23

Static stability, buoyant stability, yo-yo, and spinning-top synthetic
application contracts were each run twice against the current build. All
eight evidence files are `CURRENT` and retain full executable/source
bindings. The resulting audit is
`build/strict-current-baseline-audit-20260923/audit-after-faster-applications-20260923.json`
(SHA-256
`2afda8a4d83d9aa80b5a2275a1a0f8db37819a288e435c36a5dcae38afaa7033`). It
contains **22 independently verified rows**, **75 not verified rows**, and 2
physical rows out of scope. These are synthetic application contracts, not
original Faster Table 1 experiments.

The pinned OCCT 7.9.3 synthetic Box-minus-Sphere B-rep/STEP transport was
also run twice with current evidence and independent analytic moment checks.
The resulting audit is
`build/strict-current-baseline-audit-20260923/audit-after-faster-cad-20260923.json`
(SHA-256
`492d4b21187377870f92444baff88de1edeca55237013ef65fbe4a3ffb66fa61`):
**23 independently verified rows**, **74 not verified rows**, and 2 physical
rows out of scope. This is a synthetic CAD-method validation; it does not
provide the author CAD assets needed for Faster's original cases.

The isolated quadric-stratum matrix was also rerun twice, covering all 14
predeclared supported and unsupported surface categories. The resulting
audit is
`build/strict-current-baseline-audit-20260923/audit-after-faster-stratum-20260923.json`
(SHA-256
`d54826c1b7f63afd2ca7ffc91e78cce215b3c9b9d03f5ccf148b7cbd6d7ffc3`):
**26 independently verified rows**, **71 not verified rows**, and 2 physical
rows out of scope. This raises method-level geometry coverage only; it does
not provide original Faster inputs or prove derivatives across topology
changes.

The post-case Faster rerun is separately recorded in
`experiments/spin_it_reproduction/multi_case_2026-09-22/evidence/ctest-faster-after-cuboid-20260923.log`:
**39/39 passed**, zero failures, wall time `127.29 s`, log SHA-256
`8441f326a5d1bdd406cb1c39b72743e4607e431ba450ccf99340a9fc3a5c09f1`.

## Non-principal Cuboid yo-yo method-validation evidence — 2026-09-23

The fixed reconstructed Cuboid case in
`experiments/spin_it_reproduction/multi_case_2026-09-22/cuboid-yoyo-offaxis-pairs-20260923-manifest.json`
is the first current-version positive yo-yo optimization outside the trivial
initially-satisfied fixtures. It uses a five-degree non-principal rotation,
four adaptive passes, two bounded void components, and two independent
processes with identical problem/scientific digests. The binary normalized
constraint residual is `9.647203315061265e-7` (gate `1e-6`) and the
independent moment error is `4.607859233063393e-16`.

This closes the immediate method-validation gate only. It is a reconstructed
level-4 asset, not an author input or a level-9 paper experiment; wall-thickness
and physical validation remain unknown/out of scope.

## Next execution gate

The next algorithmic gate is the Armadillo internal design-domain audit. If
its tetrahedral path cannot produce a closed legal domain without arbitrary
input changes, preserve that blocker and move to the Faster discrete-derivative
and application integration track. Original-paper credit remains blocked until
assets and critical parameters are certified.

## Spin-It method refresh — 2026-09-23

Current-extension evidence was refreshed twice for source cage/BBW,
joint-deformation guards, source material/SLQP, dual-density material, and
adaptive split/transition inspection. The audit is
`build/strict-current-baseline-audit-20260923/audit-after-spin2014-methods-20260923.json`
(SHA-256
`185da222ebae6af48daefbfa84ff420d95c73a08dd756d0cfa56210a1187f278`):
**40 independently verified rows**, **57 not verified rows**, and 2 physical
rows out of scope. The adaptive run's `BINARY_VALIDATION_FAILED` result is
retained as an expected post-transition diagnostic, not counted as a solved
paper example.

The supplement-material runner remains blocked because the available local
PDF digest (`3add70de...`) does not match the frozen exact-material digest
(`f0acc779...`). No approximate PDF was admitted.

After adding current Givens-equation evidence, the batch audit is
`build/strict-current-baseline-audit-20260923/audit-after-spin2014-givens-20260923.json`
(SHA-256
`a54e6b7425fa4ebe77291b2a9304ea375893b749d7099257b26af7e1c1994c33`). It
contains **49 independently verified rows**, **48 not verified rows**, and 2
physical rows out of scope. This is method-level coverage only; no original
case-matrix credit is inferred.

The public official Spin-It paper and supplement are now stored under
`reference-only/spin-it-2014/` with verified digests
`SpinItPaper.pdf=3add70de...`, `SpinItPaper.txt=5cf5ab...`, and
`SpinItSupMat.pdf=f0acc779...`. Two current runs of the supplement runner
recomputed Algorithms 1 and 2 with identical scientific digest
`fe8a23d7afcb692003d63653bd2874738f435292712e117ae1731f242d3551a7`.
The project-bound audit is
`build/strict-current-baseline-audit-20260923/audit-after-spin-supplement-project-20260923.json`
(SHA-256
`015a8778b902d4887c495d120b7426c347aa51965a17049cd7474f5e61546c8f`):
**52 independently verified rows**, **45 not verified rows**, and 2 physical
rows out of scope. This verifies the exact mass-property supplement methods,
not the original object case matrix.

## Current non-trivial Heart validation

Two current-source runs of the reconstructed Heart completed with identical
scientific digest, transition trace, final partition, and exported OBJ hashes.
Both are `BINARY_MOMENT_VALIDATED`; the binary normalized constraint residual
is `3.931478211419742e-16` and the independent moment error is
`1.9189039931792828e-17`. The material volume changes from the outer volume
`0.9687810692756251` to `0.7908025536506251`, so this is an actual hollowing
run rather than an initially satisfied fixture. The machine-readable case
manifest is
`experiments/spin_it_reproduction/multi_case_2026-09-22/reconstructed-heart-current-20260923-manifest.json`.

This remains `RECONSTRUCTED_METHOD_VALIDATION` (self-built Heart, level 4),
not an author-case result. The automatically generated task split is
`build/strict-current-baseline-audit-20260923/algorithm-gap-table-20260923.json`:
52 verified current-evidence rows, 10 implemented/covered-but-unverified
rows, 35 asset/original-case blockers, and 2 explicitly out-of-scope physical
rows. These categories are planning evidence, not a completion percentage.

## Rank-deficient Faster strata: explicit fail-closed boundary

The source-motion regression now exercises elliptic-cylinder and paraboloid
inputs and requires the explicit
`RANK_DEFICIENT_STRATUM_DERIVATIVE_UNDEFINED` classification. A rank-2
quadric is a lower-dimensional stratum in the ten-dimensional normalized
coefficient sphere; a generic coefficient perturbation leaves that stratum
(for example, adding an axial linear term changes a cylinder into a
paraboloid). Returning a nominal 10-by-10 Jacobian for such a point would
silently feed a derivative from the wrong topology into Newton. The native
source-motion and Faster end-to-end regressions both pass after this guard.

This closes an evidence/contract gap, not the full rank-deficient Faster
solver. A stratum-aware parameterization and its own optimizer remain required
before the paraboloid/cylinder cases can enter the original paper matrix.

The latest guard was captured twice as strict evidence in
`build/faster-source-stratum-guard-evidence-20260923-current-run-3/` and
`build/faster-source-stratum-guard-evidence-20260923-current-run-4/`; both
reports are `CURRENT` and have identical scientific stdout hash
`4d96d74371bf78b76133c77281b1d18e926759a4ef79169f3d61cf23637c5602`. After
the change, the exact-current native/Python preset completed **275/275
tests**. This is a regression and provenance gate, not evidence that the
rank-deficient Faster optimizer or the original Faster case matrix is done.

The fixed-topology discrete Boolean differential path was also refreshed from
the current build in two independent runs:
`build/faster-discrete-moment-evidence-20260923-current-run-3/` and
`build/faster-discrete-moment-evidence-20260923-current-run-4/`. Both are
`CURRENT`, with identical transcript hash
`4f4ce25501ba93086dbac6dfc8491c707ec2f45fbdd47fed4ff48e7f3c3b02cf`.
The contact fixture validates only the registered coefficient-0 normal slice
on either side of the contact; the non-contact fixtures retain the full
ten-dimensional fixed-topology checks. This strengthens method evidence only;
topology-changing derivatives and original Faster cases remain open.

The Armadillo tetrahedral design-domain audit now has a pre-registered layer
sweep (dual boundary-distance thresholds 1 through 10). It is deterministic
across two runs, with identical result hash
`324975fe3b515f3696e2d3a65142be0b133692bf2ad79250ef059655a2faec48`, and
returns `BLOCKED_NO_MANIFOLD_CORE_IN_PRE_REGISTERED_RANGE`. No optimizer was
invoked and no disconnected component was deleted. The case remains
`SOURCE_MODEL_ONLY`; this diagnostic does not add original-paper credit.

The frozen three-layer candidate was independently replayed again on
2026-09-23. The fresh machine-readable decision is
`experiments/spin_it_reproduction/multi_case_2026-09-22/armadillo-design-domain-refresh-20260923.json`:
the candidate remains rejected with `82` non-manifold boundary edges and `165`
bad vertex links. The largest selected component itself contains non-manifold
witnesses, so deleting small islands would not be a justified repair. No
optimizer run was started. The project therefore moves to the Faster
integration track while retaining Armadillo as `SOURCE_MODEL_ONLY` and the
original paper rows as asset-blocked.

## Faster full-rank discrete Boolean-to-moment differential — 2026-09-23

The discrete material derivative now supports both ellipsoids and full-rank
hyperboloids. The hyperboloid path reuses the native source-motion differential,
including moving truncation roots and cap centroids, and verifies that the
reconstructed source triangulation exactly matches the mesh used by the
Boolean evaluator before differentiating the Boolean vertices. Rank-deficient
quadric strata remain explicitly rejected and are classified separately from
Boolean topology transitions. The previous public enum spelling
`DiscreteEllipsoid` remains as a source-compatible alias for the generalized
`DiscreteFullRankQuadric` method.

The registered hyperboloid integration fixture uses an actual capped
hyperboloid/design-domain Boolean, compares all ten columns of the complete
mass-moment map against central differences at three step sizes, checks
unchanged provenance topology, checks moments with an independent signed-
tetrahedron integrator, and repeats the evaluation for determinism. Relative
full-map errors were `1.6083125573298586e-8`, `1.8001270686077591e-9`, and
`2.5455573060025292e-10`. This validates only the fixed-topology, fixed-radius
selection branch; it does not claim a derivative across Boolean topology or
radius-selection changes.

The full registered Faster CTest label completed **40/40 tests** in `114.59 s`
(exit code 0). Its log is
`build/faster-discrete-quadric-ctest-20260923.log` (SHA-256
`efb016c6f9ba1a1929e89bbb191b231e1ecb33b6939e3508f266d35173f7cada`). The
suite covers the existing Faster integration group, not the full project test
matrix. This closes a discrete derivative integration gap, not the Faster
paper's original case matrix, 2D/floating/CAD requirements, or strict final
reproduction status.

## Targeted Spin-It algorithm refresh — 2026-09-24

Against the rebuilt `build/strict-algorithms-exact-current-20260922` tree, the
following current-source regressions were built and executed after the
source-cage SLQP regression had been extended with an explicit check that its
unsupported disk-trace mode rejects before any solve or post-solve state
mutation:

- `spin-it-2014-source-cage-slqp`: **1/1 passed**. It exercised all five
  objective modes, the actual-source callback tape, the complete 26-handle
  state, independently recomputed KKT quantities, resource-latch behavior,
  fail-closed geometry/feasibility status, and same-process deterministic
  repeats. This is a bounded synthetic method regression, not a paper case.
- `spin-it-2014-source-cage-binary-material`,
  `spin-it-2014-source-cage-adaptive-solver`,
  `spin-it-2014-source-cage-adaptive-rebind`, and
  `spin-it-2014-dual-material-geometry-export`: **4/4 passed**. These exercise
  adaptive pass/rebind behavior, binary endpoint admission and exact boundary
  auditing, the rounding-failure and nonmanifold negative controls, and the
  independent dual-material geometry exporter.

The binary-material positive controls use a known stationary configuration;
they establish the implemented adaptive-to-binary geometry/export path, but
they do not alone establish a nontrivial successful joint-deformation
optimization. A second fixed synthetic method fixture now applies a
shape-regularized cage objective from a uniformly scaled outer cage, requires
an accepted SLQP objective trial and a changed handle state, then runs adaptive
completion, beta rounding, binary mesh emission, independent mass-moment
integration, exact boundary audit, and the `1e-6` binary constraint gate. The
registered `spin-it-2014-source-cage-binary-material` test passed with this
additional chain (**129.13 s**). Its evidence is a single deterministic test
execution and is not an author-paper case or a two-process repeat.

The previously recorded 98-handle joint case remains
`RESTORATION_STALLED` at a verified containment boundary and was not rerun.
Neither that failure nor this synthetic positive fixture is a continuous
infeasibility certificate. The new chain grants no sealed artifact,
manufacturing clearance, physical validation, or original-paper case credit.
Overall paper-reproduction status remains `PARTIAL_BLOCKED`.

The full exact-current CTest preset was rerun after the topology-side
validation change: **275/275 tests passed** in 596.84 seconds. The strict
audit remains `PARTIAL_BLOCKED`; no test count is being used as a completion
percentage.

## Fixed-stratum tangent and current test gate

`QuadricStratumTangent` now records the local tangent constraints for the
normalized coefficient strata used by degenerate quadrics. The registered
regression checks dimensions 3 (plane), 7 (elliptic cylinder), 8
(paraboloid), and 9 (full-rank ellipsoid), together with finite rows and a
translation pullback. `QuadricSourceMotion` consumes this descriptor for a
diagnostic, but remains fail-closed on rank-deficient source motion because a
stratum tangent is not yet a complete rank-deficient Newton solver.

The intermediate `run-5` and `run-6` batches are retained but became stale
when the test-budget experiment changed the source binding. The final
registered runner uses CTest's explicit `--timeout 420` only for the evidence
command while retaining the source CMake timeout and sequential `-j1`
execution. The two current batches in
`build/registered-ctest-evidence-20260923-stratum-run-7/` and
`build/registered-ctest-evidence-20260923-stratum-run-8/` each completed
276/276 with no failed or `Not Run` tests. Their transcript hashes are
`be79ccc15942f28663ec2cada5d5d527eca02816fc2f21040e9f6eabcecee9e5` and
`018af8b2c0f3d8b9b4d458854cf5b317893735ed9609a5554dc601a594ed42ee`.
Both are `CURRENT` under the registered-test checker. The timeout batch is
retained as failure evidence, not relabeled.

The intermediate audit attempt was intentionally retained as stale evidence
after the source-binding change. It is superseded by the current audit refresh
below: all mapped checker batches were rerun against the current source before
the audit was accepted. The current strict status is `PARTIAL_BLOCKED`, not
`INVALID_BLOCKED`.

These are implementation and regression gates, not original-paper case
results; the strict status remains `PARTIAL_BLOCKED`.

## Current audit after independent refresh

All 21 registered independent checker entries were regenerated twice against
the current source tree. The valid audit is
`build/strict-current-baseline-audit-20260923/audit-after-full-independent-refresh-20260923.json`
with SHA-256
`6eb2cef9beb4bd6e87066d7021d77d519069fb3095470a4b414a154a1d26ea9f`.
It reports `PARTIAL_BLOCKED`: 53 verified requirements, 44
blocked/not-implemented requirements, and 2 physical measurements outside
scope, out of 99 total. This is evidence coverage, not a completion
percentage.

The remaining blockers are requirement-specific: missing original Spin-It
assets and parameters, Faster Table 1/Figures 10–13 and theorem gates,
shared sealed-geometry/moment/simulation/provenance checks, and the remaining
manufacturing/case-matrix requirements. Synthetic regressions remain clearly
separated from original-paper credit.

The exact-rational Faster theorem-condition runner was refreshed twice and is
attached as supporting evidence in
`build/strict-current-baseline-audit-20260923/audit-after-full-refresh-with-theorem-support-20260923.json`
(SHA-256
`77b007f65ab6cc80d48b1686ca5fc3479f218a4c45a71b384bab0511aa0d8621`). It
checks the published numerical hypotheses and canonical replay, but its audit
effect remains supporting-only: `theorem_coverage=NOT_GRANTED_BY_THIS_RUNNER`.

## Independent mass-moment gate and final current audit — 2026-09-23

The shared independent-moment checker validates closed-form cube and
tetrahedron moments plus a translated two-cube composite with non-diagonal
central inertia. Both current runs are `CURRENT`, with maximum absolute error
below `5e-13`, and remain synthetic method evidence:

`build/shared-independent-moments-evidence-20260923-current-final-a/`  
`build/shared-independent-moments-evidence-20260923-current-final-b/`

The accepted current audit is
`build/strict-current-baseline-audit-20260923/audit-after-independent-moments-final-20260923.json`
with SHA-256
`71226549f91c9891d860ad5af09e8db0b7243bb7a272486cb0b06285c9e50034`.
It reports **54 verified**, **43 blocked/not implemented**, and **2 physical
measurements out of scope** out of 99 requirements. The status remains
`PARTIAL_BLOCKED`; this count is evidence coverage, not a reproduction
percentage. Original case matrices, CAD/asset gates, general Faster theorem
gates, sealed geometry, simulations and visualization remain open.

## Algorithm-first baseline and test refresh — 2026-09-23

The current audit was converted into a machine-generated planning table (no
hand-entered pass flags):

`build/strict-current-baseline-audit-20260923/algorithm-gap-table-final-20260923.json`

Its source audit is bound by SHA-256
`71226549f91c9891d860ad5af09e8db0b7243bb7a272486cb0b06285c9e50034` and its
four buckets are:

| Bucket | Count | Meaning |
|---|---:|---|
| `verified_current_evidence` | 54 | Current independent evidence, including method-level or synthetic fixtures |
| `implemented_or_covered_but_unverified` | 8 | Code or declared scope exists, but an independent gate is still missing |
| `asset_or_original_case_blocked` | 35 | Original asset/parameter, case matrix, or downstream evidence is unavailable |
| `postponed_out_of_scope` | 2 | Physical measurements explicitly excluded from this computational phase |

The two native standalone regressions that previously failed before reaching
the checker were rerun after rebuilding the historical compatibility build and
preloading its matching GCC runtime. Both passed: the fixed-topology Faster
discrete-moment test and the pinned reference-Boolean test. The complete
`test_spin_it_strict_audit.py` discovery then ran **43/43** tests in
`285.966 s` with exit code 0. This repairs a stale-build/runtime evidence
problem; it does not change the strict audit or promote any original-paper
case.

## Readiness-bound baseline refresh — 2026-09-23

The case readiness inventory was re-bound to the current source-audit bytes;
the previous readiness hash was stale and the audit correctly rejected it.
After updating only that binding (case statuses and missing parameters were
unchanged), the fresh accepted audit is

`build/strict-current-baseline-audit-20260923/audit-with-readiness-20260923.json`

with SHA-256
`a55d1b73d71a25d5e1e8240be8f025d658d73ea9a33faa7979a9bcf6a1198cdd`.
It remains `PARTIAL_BLOCKED` with the same `54/43/2` requirement split, but now
records `case_asset_readiness.status=CURRENT` instead of `NOT_SUPPLIED`. The
corresponding machine-generated table is

`build/strict-current-baseline-audit-20260923/algorithm-gap-table-readiness-final-20260923.json`

and is bound to the fresh audit rather than the superseded no-readiness
invocation.

## Faster full-rank integration evidence — 2026-09-23

The native full-rank chain now has a separate evidence runner:
`experiments/spin_it_reproduction/strict/run_faster_fullrank_end_to_end.py`.
Two independent processes executed the eight registered modes (adapter,
Boolean/moment chain, solver guard, independent integration, topology-change
negative control, theorem-condition audit, replay, and negative-path guard).
Both are `CURRENT` and share normalized scientific digest
`3df418a61872b7125599935dbcb03bb4aaf84d94b6885ab1ddf3a464106b81dc`:

- `build/faster-fullrank-end-to-end-evidence-20260923-c/`
- `build/faster-fullrank-end-to-end-evidence-20260923-d/`

The older off-target yo-yo record was produced before residual-decreasing
topology transitions could be relinearized. Its `LINE_SEARCH_EXHAUSTED` state
is preserved as historical evidence, but it is not the current solver contract
or an infeasibility certificate. The current regression verifies both the
accepted relinearization and the non-decreasing rejection branches. This
establishes a synthetic full-rank integration chain, not a Faster author case,
global optimality theorem, CAD result, or Table-1 reproduction.

## Reference Boolean geometry gate — 2026-09-23 update

The earlier line above describing a Cherchi/CGAL cross-check as validated was
too broad: the prior runner established deterministic output and matching raw
moments, but did not establish equality of the represented solids. The
updated runner and strict audit now compare the actual output geometry using
the exact represented-set relation, and retain triangle-boundary signatures
as diagnostics rather than incorrectly requiring identical triangulation.

Two new independent batches are in
`build/faster-cherchi-geometry-final-20260923-a/` and
`build/faster-cherchi-geometry-final-20260923-b/`. Both are `CURRENT`; every
Cherchi fixture is byte-identical across within-batch repeats, and all ten
raw moments agree within `3.55e-15`. Exact represented-set comparison agrees
for the ellipsoid, but reports a mismatch for the one-sheet hyperboloid,
two-sheet hyperboloid, plane, and paraboloid. Those four outputs also have
different oriented triangle sets, despite equal vertex/face counts.

The independent audit deliberately refuses to grant the reference-Boolean
requirement for these results. This is a synthetic backend discrepancy, not
an author-case result and not yet a root-cause attribution. The required next
step is a bounded geometric comparison for the alternate triangulations; until
that comparison is independently justified, Boolean equivalence remains
blocked. See
`docs/faster_boolean_geometry_crosscheck_2026-09-23.md` for the per-fixture
results and exact scope.

The strict-audit Python file was freshly rerun after the evidence-binding
changes: **45/45 tests passed** in `284.128 s`. This includes a regression
that verifies matching moments do not override a geometric mismatch, and the
two-process Spin-It synthetic material audit. It is not a replacement for the
full native CTest suite or any original-paper experiment.

## Final registered-test evidence refresh — 2026-09-23

Adding the public-asset and Armadillo refresh records changed the effective
source digest, so earlier registered-test receipts were intentionally rejected
as stale. Two complete sequential reruns were then performed against the same
strict build and configuration:

- `build/registered-ctest-evidence-20260923-current-refresh-1a/`
- `build/registered-ctest-evidence-20260923-current-refresh-1b/`

Each executed all `276/276` registered tests with no failures and no `Not Run`
entries. Elapsed times were `1006.44 s` and `916.70 s`. The current strict
audit is
`build/strict-current-baseline-audit-20260923/audit-final-after-ctest-refresh-20260923.json`
(SHA-256 `c446f13c847aaf952b652b6e50b351d407fbb5c54b25348917fa86cffdb91108`),
with readiness `CURRENT` and status `PARTIAL_BLOCKED`: 54 verified, 43
blocked/not implemented, and 2 physical requirements out of scope. The fresh
gap table is
`build/strict-current-baseline-audit-20260923/algorithm-gap-table-final-ctest-refresh-20260923.json`
(SHA-256 `26a66a522e384f43d207548d3051a6c07e4cafab34b700d139dfd6f0a6740a74`).
This restores a valid current baseline; it does not grant original-paper case
credit.

## Multi-case Faster algorithm evidence — 2026-09-23

As the next algorithm-first step, the current-source reconstructed Faster
runner completed the eight-case matrix with 50 deterministic random starts per
case and up to eight restarts. The inputs are explicitly
`RECONSTRUCTED_ASSET`; they are not author models or a paper-parameter
recovery. Result SHA-256:
`58c030452ae14855c2ea41898828e01663fbcd2985531d1443aa4b9679b34ed3`.

The nominal success counts were 49/50 (static tilted vase), 48/50 (floating
tea infuser), and 50/50 for the two floating-island orientations, heart,
noodles, jellyfish, and flower spinner. The two incomplete cases contain
only retained restart failures (`LOCAL_MAXIMUM`); the flower run also retains
five `INVALID_MATERIAL_ARTIFACT` and one `RESOURCE_LIMIT_REACHED` attempt.
No failure was removed from the denominator. This closes a meaningful
multi-model algorithm regression and restart audit, but it is not evidence of
the Faster paper's original Table 1 and does not change the strict audit
status (`54` verified, `43` blocked/not implemented, `2` physical out of
scope).

The reconstructed runner also received a fail-closed `passed` correction:
partial convergence is no longer encoded as a pass. The corrected current
source matrix is stored at
`build/faster-reconstructed-table1-20260923-strict-pass-semantics/`; its
output SHA-256 is
`9f622cf510bfe334330e5b61ceaa5f7fe781341329e9b0b1c9f05e7a9df3be3d`.

## Registered-test evidence after runner correction — 2026-09-23

Because the runner change affected the declared source observation, the prior
registered-test receipts were not reused. Two new complete sequential runs
were executed:

- `build/registered-ctest-evidence-20260923-after-runner-fix-1a/`
- `build/registered-ctest-evidence-20260923-after-runner-fix-1b/`

Both executed all `276/276` tests, with no failures and no `Not Run` entries;
they are both `CURRENT`. Their elapsed times were `893.42 s` and `889.74 s`.
The normalized executed-test-ID digest matched exactly in both runs. The
comparison is recorded in
`build/registered-ctest-evidence-20260923-after-runner-fix/repeat-comparison.json`.

The refreshed strict audit remains `PARTIAL_BLOCKED`, with `54` verified,
`43` blocked/not implemented, and `2` physical requirements out of scope:

- [audit-after-runner-fix-20260923.json](/data/xvxinyao/moment-guided-spin-v1/code/make-it-balance-paper-sim/build/strict-current-baseline-audit-20260923/audit-after-runner-fix-20260923.json)
- [algorithm-gap-table-after-runner-fix-20260923.json](/data/xvxinyao/moment-guided-spin-v1/code/make-it-balance-paper-sim/build/strict-current-baseline-audit-20260923/algorithm-gap-table-after-runner-fix-20260923.json)

## Current Spin-It reconstructed fixtures and Armadillo gate — 2026-09-23

The Heart-top and non-principal Cuboid yo-yo fixtures were rerun from the
current source after the evidence refresh. Each was executed twice at the
registered level-4 cap; both repeats converged, sealed, and produced identical
scientific/output hashes. The result is
`build/algorithm-first-current-refresh-20260923-after-runner-fix/heart-cuboid-level4.json`
(SHA-256
`19bbefca9e43706e57b8b270d04a7ad01a803a6a8b6d4c52086d9657234715d2`). These
remain `RECONSTRUCTED_ASSET` method regressions, not original Heart/Cuboid
paper cases or level-9 evidence.

## Restart event ledger — 2026-09-24

`run_paper_table1.py` now stores and checks the complete deterministic solver
event stream for each generated start, including failed evaluations, in the
scientific result payload. The exact change and test receipts are recorded in
[`faster_restart_event_ledger_2026-09-24.md`](faster_restart_event_ledger_2026-09-24.md).
This closes an experiment-runner evidence gap only; it adds no original Table 1
case or reproduction credit.

Independent-process replay found and fixed a second issue: aggregate
`*_seconds` fields had leaked into the scientific hash. The regression was
observed before the fix and then passed after it. Two retained mini-run payloads
now have identical scientific hash `2502022eb005299c2e28e3beb9e753da1427574e2cd893ebc58c6392540cf402`,
while their full JSON file hashes differ only because timing remains in the
payload. The replay is synthetic and does not add author-case credit.

The Armadillo three-layer candidate was also recomputed from the unchanged FAST
volume mesh. The result is
`build/armadillo-design-domain-refresh-20260923-current/contact-topology.json`
(SHA-256
`5b59c4c8f4cc647de194385cee127c2ab13814aa69424f58aa2ba37320c87a5f`). It
still has 53 face-connected design components, 82 non-manifold boundary edges,
165 bad vertex links, and 52 cross-component edge contacts. No optimizer was
invoked; the candidate remains rejected as a legal design domain.

## Joint deformation and dual-material fallback path — 2026-09-23

Because the Armadillo design-domain gate remains blocked, the next algorithm
path was exercised on registered synthetic fixtures. The current build ran
the asymmetric complete-cage joint deformation, cage binding, BBW, and
dual-material regressions twice; both runs passed `4/4` with no skipped tests.
The bound transcripts and executable hashes are recorded in
`build/spin2014-joint-dual-current-20260923/RUN_MANIFEST.md`.

This verifies execution and repeatability of the joint deformation and
dual-material implementation paths. It does not upgrade the missing Armadillo
or Sumone original cases, and the joint fixture remains synthetic method
evidence.

## Faster discrete Boolean differential focus — 2026-09-23

The next algorithm-first gate was run twice from the current strict build,
covering the polyhedral differential, source-feature provenance, canonical
ordering/adapter/chain, negative-path classification, solver integration,
Boolean/triangulation/inverse-Gauss/source-motion, discrete moments, stratum
tangents, and quadric-moment evaluation tests (18 registered tests per run).
Both runs passed `18/18`, with no skipped or `Not Run` tests. The raw receipts
are:

- `build/faster-discrete-derivative-focused-20260923-a/ctest.log`
  (SHA-256 `ae775c861baa714c9f71f0765f39dabf48644a7d769b8d5128a2c4b8961035ea`)
- `build/faster-discrete-derivative-focused-20260923-b/ctest.log`
  (SHA-256 `11c1b9632b279fb992b0337d3fc3567a9811d369fd4a53a29daf4d51c832e781`)

The focused evidence confirms the currently implemented scope: a fixed
finite-ellipsoid Boolean topology has a finite-difference-checked discrete
Jacobian and deterministic source-feature identities; rank-deficient and
unsupported surface types and input contacts fail closed. A changed Boolean
topology is accepted only after an independent residual decrease and then uses
a fresh derivative map; otherwise it fails closed and is not passed to Newton.
This is a real Faster algorithm
verification increment, but it is not yet the paper's full discrete Boolean
derivative implementation: paraboloid/cylinder/hyperboloid variants,
topology-changing derivatives, the Cherchi reference backend, and the full
Table 1 inputs remain open or blocked.

The strict audit was regenerated after this focus run at
`build/strict-current-baseline-audit-20260923/audit-after-discrete-focus-20260923.json`
(SHA-256
`53ed341c9d0cefea6b5e221b3a8ce75c7f0f5361873f10f2aed723a314298c30`). It
remains `PARTIAL_BLOCKED` with 54 verified requirements, 43 not implemented or
blocked requirements, and 2 physical requirements out of scope; the focused
test receipts do not change the original-case admission status.

## Non-principal Cuboid yo-yo reachability diagnostic — 2026-09-23

The previously reported reconstructed Cuboid yo-yo was not a sufficient
non-trivial optimizer demonstration: its registered target frame is the
axis-aligned frame and its final constraints are already zero. A new,
pre-registered synthetic diagnostic used a cuboid rotated by 31 degrees about
the normalized `(1,2,-1)` axis while retaining the target z axis. Its geometry,
weights, tolerance and initial fill were frozen before execution.

The source-grid/maximum-level validation correctly rejected the accidental
`24^3`/level-4 request. Two valid `32^3`/level-5 runs then tested the same
geometry with initial levels 2 and 4. They stopped at the first adaptive pass
with `LINE_SEARCH_EXHAUSTED` and residuals about `3.07e-2` and `2.42e-2`,
respectively. The relevant result files are:

- `build/rotated-cuboid-yoyo-g5-20260923.json`
  (SHA-256 `6f90340757afa347ceee94c120f95cab8030833b73cc94efda5f73a96d95e05b`)
- `build/rotated-cuboid-yoyo-g5-initial4-20260923.json`
  (SHA-256 `9b763c9bab0b03ebc2e582b22618ed51811d8592e52dc008d21a38c7b9558f6a`)

This is evidence that the present coarse-to-fine route does not yet supply a
feasible direction for this non-principal fixture; it is **not** an
infeasibility certificate and is not counted as a successful method case. The
next implementation target is therefore a frame-aware moment functional (or
an equivalent verified coordinate transform) plus an explicit linear
reachability preflight, rather than another unrecorded parameter retry.

## Frame-aware adaptive Spin-It 2014 functional — 2026-09-23

The adaptive voxel request now carries an explicit right-handed orthonormal
target frame. The material grid, source geometry, and exported artifact remain
in model coordinates; only the raw ten-moment vector, objective, equality
constraints, and their derivatives are evaluated in target coordinates. The
identity frame keeps the legacy input hash and behavior unchanged. A
non-identity frame is part of the request hash and artifact diagnostics, so it
cannot be silently mixed with a model-axis result.

The native regression checks transformed objective/constraints and a
finite-difference Jacobian, plus rejection of non-orthonormal frames. The
Python facade independently validates and forwards the frame. The focused
native/objective/Givens/end-to-end/Python group passed `6/6` from
`build/strict-algorithms-exact-current-20260922` after the change.

This closes an interface-semantic gap: a fixed nonprincipal yo-yo target can
now be expressed without pre-rotating a mesh before voxelization. It does not
turn either preserved `LINE_SEARCH_EXHAUSTED` Cuboid run into a convergence or
feasibility claim; the next required evidence is a fixed-domain reachability
preflight and a new, independently repeated solve.

## Frame-aware source-surface adaptive path — 2026-09-23

The same target-frame contract now applies to the source-surface adaptive
solver used by the reconstructed Cuboid workflow. Source octree coordinates,
emitted cuboids, and independent signed-tetrahedron material integration stay
in model coordinates. The fill functional transforms only the ten raw moments
for objective/constraint evaluation; the independent binary check transforms
its separately integrated moments before evaluating the same paper Givens
constraints. This prevents an optimizer and its final binary check from
silently using different axes.

The request validates a finite, right-handed orthonormal frame. Identity keeps
the legacy request hash, while a nonidentity frame is bound to the hash. The
Python source-material bridge requires the declared `SpinArtifactProblem`
origin and axis to agree exactly with the frame origin and its third column;
the native diagnostics retain both frame arrays. This is algorithm-path
coverage, not a claim that the reconstructed Cuboid is an author asset or that
the preserved 31-degree nonconvergence has become feasible.

## Faster Boolean topology strata relinearization — 2026-09-23

The full-rank Faster Newton path previously rejected every candidate whose
Boolean output had a different topology hash, even when that candidate had a
freshly evaluated discrete mass-moment map and a lower independent residual.
The result was an artificial `RESOURCE_LIMIT_REACHED` after consuming the
paper-limit function evaluations on a perturbed spinning-top fixture.

The solver now accepts a topology-changing candidate only if its independently
evaluated stopping residual decreases. The event is explicitly recorded as
`topology_transition_relinearized`; the following Newton iteration starts from
the candidate's new Boolean geometry and derivative map, rather than reusing a
derivative from the old stratum. A non-decreasing transition remains rejected.
This preserves the required fail-closed behavior for invalid or unavailable
derivatives while allowing the solver to traverse valid, separately evaluated
topology strata.

The Faster API regression, Table-1 miniature runner, Figures miniature runner,
and actual material bridge all pass after this change. This validates the
algorithmic transition handling on reconstructed fixtures only; it does not
replace the missing author CAD/mesh inputs or the required full paper matrix.

## Target-axis Cuboid runner and binary-gate result — 2026-09-23

`experiments/reconstructed_cuboid_yoyo/run.py` now accepts a frozen physical
`--target-axis` and constructs a deterministic right-handed target frame whose
third column is that axis. The geometry, material export and independent
integration remain in model coordinates; only the moment objective and Givens
constraints are transformed into that target frame. The registered
`reconstructed-cuboid-target-frame-python` test checks normalization,
orthogonality, orientation and this column convention.

The pre-registered reconstructed offaxis-pairs Cuboid was then run twice with
a 5-degree model rotation and the 1-degree target axis
`(0.01745240643728351, 0, 0.9998476951563913)`. Both output directories,
`build/reconstructed-cuboid-yoyo-target-axis-1deg-preflight-20260923/` and
`build/reconstructed-cuboid-yoyo-target-axis-1deg-preflight-20260923-repeat/`,
have identical problem digest
`d29713c2c20d910b8807970a12adbb575db0c57b883a4f64e0937143cbaa88a0`,
scientific digest
`d9e07705488de9db29e5f237c96db2516e0f66d6c20356b0010091610a5c2e24`,
trace hash
`b163b53c47553e288707ff765ba01b57d195c579124327e8fd1f3e4cb562d625`,
and output hash
`6cdf7e71925b10c06b63e63ea9b84f13dfb7860628f5da9a7c80bc26a06e436c`.

This is not a successful method case: every continuous adaptive pass reports
the first-order KKT witness, but the prescribed `beta >= 0.5` binary rounding
fails the independent paper constraint gate. The raw constraint residual is
`8.585414966673097e-3` (normalized `2.55865304668936e-5`), above the fixed
`1e-6` binary tolerance. Its independently reintegrated normalized moment
error is `4.0657581468206406e-16`, which isolates the observed failure to the
rounded constraint rather than material-moment disagreement. It remains a
`RECONSTRUCTED` level-4 diagnostic, not an author case, an infeasibility
certificate or a sealed artifact.

The same frozen target/domain case was evaluated once at G5 because its G4
trace showed every fractional leaf already at the source-depth cap. The G5
result completed six converged passes but retained 264 fractional level-5
leaves. Its binary normalized constraint residual decreased to
`1.447240363817546e-5` yet remains above the frozen `1e-6` gate; independent
moment error is `9.215718466126786e-16`. This localizes the failure to the
binary material for this configuration and does not justify changing the
threshold or marking the artifact sealed.

## Faster full-rank relinearization replay — 2026-09-23

The eight-mode Faster evidence runner's prior status field expected an
off-target yo-yo line-search failure. That expectation became stale after the
solver began accepting independently verified residual-decreasing topology
transitions and relinearizing from the new Boolean geometry. The runner now
requires the resulting `CONVERGED` status and independently verified final
geometry. Two current runs—`build/faster-fullrank-end-to-end-evidence-20260923-e/`
and `...-20260923-f/`—completed all eight modes and share scientific digest
`18290eb2bec815fab58a197d484ed91a5910881423750972a770b20f02cbaafd`. The
transition result has residual `3.2089567488379244e-15`, and matching output
and trace hashes. The enclosing receipt hashes differ because its output
records include absolute transcript paths; the path-independent scientific
digest and all native transcript content hashes match. This is synthetic
integration evidence only; author cases and the complete Faster paper matrix
remain blocked.

## Faster discrete full-rank quadric derivative refresh — 2026-09-23

The opt-in discrete Boolean moment differential now covers finite ellipsoids
and full-rank hyperboloids. For a hyperboloid, it reuses the native source
motion map for fixed-normal surface vertices, truncation-root vertices, and
cap centroids. Before using that derivative, it rebuilds the source
triangulation and requires exact agreement with the mesh and provenance used
by the Boolean operation. Paraboloid, elliptic-cylinder, and plane strata
remain fail-closed and are reported as rank deficient rather than topology
transitions. The old `DiscreteEllipsoid` C++/Python spelling remains an alias;
the canonical name is now `DiscreteFullRankQuadric` and its method identity is
versioned as `discrete-fullrank-quadric-boolean-transport-v2`.

A capped hyperboloid integration fixture exercises all ten coefficient
directions against the complete Boolean-output mass-moment map, with three
central-difference steps, fixed source-feature topology, an independent
signed-tetrahedron moment oracle, and deterministic replay. Relative full-map
errors are `1.6083125573298586e-8`, `1.8001270686077591e-9`, and
`2.5455573060025292e-10`. These results apply only within one fixed Boolean
topology and truncation-radius selection branch; derivatives across topology
or radius-branch transitions remain unsupported.

Verification: the full Faster CTest label passed **40/40** in `114.59 s`; the
latest complete build/test run passed **278/278** in `886.30 s`, exit code 0.
The complete log is
`build/full-ctest-evidence-20260923-discrete-hyperboloid-final.log`, SHA-256
`42fb0d2dce7fd7ca5ef508abd9b81975598592043376d80d73b6d8b90a0dea18`.

This suite result is not paper-case completion. The strict project status
remains `PARTIAL_BLOCKED`: original inputs/parameters, rank-deficient solvers,
full 2014 and Faster experiment matrices, CAD/floating application integration,
and final material/physical gates remain outstanding.

Two independent strict method runs are retained under
`build/faster-discrete-fullrank-quadric-evidence-20260923-a/` and
`build/faster-discrete-fullrank-quadric-evidence-20260923-b/`. Both are
`CURRENT`, both report
`SYNTHETIC_DISCRETE_FULLRANK_QUADRIC_METHOD_VALIDATED`, and both native
transcripts have SHA-256
`96605ada7257b1dbc8b68682ae6c5d54331c34983313b314a8aa05fbf135d51c`.
The strict audit integration test admits exactly
`faster.discrete_moment_map.differentiable` and
`faster.geometry.inverse_gauss_map`; it does not grant original-case or
independent-Boolean verification.

## Full joint attempt, dual-region exports, and rounding diagnosis — 2026-09-23

The follow-on batch is documented in
[spin_joint_batch_2026-09-23.md](spin_joint_batch_2026-09-23.md). This is an
append-only scope update, not a replacement for historical results above.

- Two fresh complete-attempt runs of the frozen asymmetric 98-handle/actual-BBW
  joint fixture have identical native scientific files. Both stop at
  `RESTORATION_STALLED`, principally after internal primitive/outer-region
  containment rejections. The objective and handles do change, but the binary
  material remains `NOT_ADMITTED`. No complete successful joint case is granted.
- The existing dual-material regression now exports the actual low-density shell
  and high-density insert separately. Two independent processes reproduce the
  same meshes, and a separate Python tetrahedral integrator reproduces their
  combined moments to about 1.8e-15. The reconstructed Sumone's dense insert is
  still empty; it remains a negative case with its legacy 0.005 binary tolerance
  explicitly disclosed. This is not complete dual-material optimization or
  manufacturing verification.
- A bounded historical-receipt diagnostic explains the saved Cuboid G4/G5 binary
  constraint errors by per-cell, per-level original-rule rounding moment deltas.
  It does not reintegrate historical geometry or certify infeasibility. G8/G9
  experiments have not been run by this batch.

No original-paper case count, strict catalogue completion grant, or final
project status is upgraded. The current project remains `PARTIAL_BLOCKED`.

After the new protocol/export/rounding registrations, the current exact-geometry
build passed **281/281** CTest entries in **903.24 s**, exit 0, with no observed
Not Run or skipped entries. Log: `build/spin-joint-dual-batch-ctest-20260923-b.log`,
SHA-256 `fd15ded230a2964f8daf75018bce4c642d9db9c59252c12540a1037b3f185cf3`.
The two joint runs and two dual-export runs remain `CURRENT` under a fresh
check with the declared native extension loaded. These test/freshness results
do not change the failed scientific outcomes recorded above.

## Faster default-start follow-up — 2026-09-23

An additional frozen-input synthetic spinning-top campaign uses the strict
full-rank Boolean-to-moment evaluator and supplies no initial quadric. Its new
seed set 20–24 was preregistered before these ten fresh-process runs; all
results and restart events are retained under
`build/faster-default-start-seeds-20-24-20260923/`. Seeds 21–24 converged and
passed final geometry checks, with exact within-seed repeat agreement across
the event stream, trace, output, restart evaluations, and status. Seed 20
reached its evaluation budget in both repeats after the truncation cap could
not be certified outside the design domain. That is a solver/geometry
evaluation blocker for this initialization, not an infeasibility certificate.

The pre-registration and values are detailed in
[faster_default_start_batch_2026-09-23.md](faster_default_start_batch_2026-09-23.md).
This adds bounded evidence that the integrated discrete-map Newton path can
solve this synthetic problem from several answer-independent deterministic
starts; it does not establish robustness over the quadric space, validate
rank-deficient solvers, or reproduce any Faster author case. Project status
remains `PARTIAL_BLOCKED`.

Targeted post-batch regression: **19/19** exact-region, source/dual-material,
canonical Faster, discrete moment, stratum-tangent, plane-family and full-rank
end-to-end tests passed in **102.76 s** under
`build/strict-algorithms-exact-current-20260922`.

## Faster degeneracy scope audit against the paper — 2026-09-24

The primary paper's Eq. 14 optimizes a normalized coefficient vector in
`S^(n-1)` (here `n=10`), but Section 5.2 explicitly introduces the
inverse-Gauss triangulation for **non-degenerate** quadrics. The supplementary
discretization derives its point map using `A^-1`, which is not defined for a
singular quadratic block. At the same time, Sections 4.2–4.3 describe static
and buoyant solutions with planar interfaces. Therefore the implementation
must distinguish (a) the full ten-coefficient Newton path for supported
non-degenerate states, (b) any paper application whose optimum is a plane, and
(c) other classified degenerate quadrics. A four-parameter plane-family
Jacobian is useful evidence for (b)'s moment map, but is not an Eq. 14
ten-dimensional solver. No general rank-deficient solver is currently
claimed. [Primary paper](../reference-only/spin-it-faster-2024/sif-final.pdf),
Section 4.2–4.3 and 5.1–5.2; [supplement](../reference-only/spin-it-faster-2024/sif-supp-final.pdf),
Section 1.

This source audit narrows the next implementation question without reducing
the requested reproduction scope: the planar static/buoyant paper path needs a
faithful solver/terminal treatment or evidence that the authors' full-space
iteration reaches and verifies that planar solution; cylinders and
paraboloids remain separate classified strata. The existing plane-family
value/10x4 derivative test does not close that gap.

## Algorithm-first integration refresh — 2026-09-24

The native unknown-answer probe now exposes both deterministic starts, per-start
evaluation totals, and attempt-start/terminal plus rejected/failed solver
events. Its strict Python parser checks the frozen two-start development budget
and compares complete start/event evidence between independent process
transcripts. The parser's seven focused unit tests pass; the native diagnostic
target builds and a single seed-2 run was parsed from its retained transcript.
That run fails the first start's Boolean-boundary audit and exhausts the second
start's 60-evaluation budget (`[1,60]`); no final geometry is produced. This is
not the five-seed repeated protocol and provides no paper-case credit.

Four Spin-It source-cage integration regressions passed on the current build:
component clearance, exact material regions, binary material, and the adaptive
solver. Together they exercise adaptive-solver output through binary material
validation and exact region/clearance checks on synthetic and known-solution
fixtures. They do not validate the separate nontrivial asymmetric joint
deformation batch: its retained result remains `RESTORATION_STALLED` before
binary admission.

The focused diagnostic source and transcript are documented in
[`faster_two_restart_diagnostic_2026-09-24.md`](faster_two_restart_diagnostic_2026-09-24.md).
The overall requirement catalogue still has unmet algorithm integrations,
author-case assets, scale runs, and final geometry/simulation gates; overall
status remains `PARTIAL_BLOCKED`.

## Hyperbolic-cylinder source and Boolean derivative — 2026-09-24

The native finite triangulator now constructs the two simple-spectrum
hyperbolic-cylinder branches, closes their finite grid boundaries beyond the
design region, and records interface/cap topology. The reduced source-motion
map uses the seven-dimensional fixed-stratum chart; its source vertices match
central finite differences with relative error `4.42e-10` at step `1e-5`.
The actual CGAL design/quadric Boolean material map now uses the same reduced
chart and passes independent ten-moment Jacobian finite differences at three
step sizes; the `1e-5` relative error is `2.71e-10` with 32 constructed
intersection vertices.

This uncovered and fixed a shared Boolean-transport classification defect:
dimension 3 is the missing-feature sentinel, not a constructed feature. A
single-source Boolean vertex now inherits the moving source-face velocity by
barycentric interpolation (or zero velocity for a design-only feature), with
its provenance still emitted. This does not extend the full-rank Eq. 14 Newton
solver to a rank-deficient stratum, certify across stratum/topology changes, or
reproduce a Faster author case. The surface is the project's finite grid-and-cap
construction, not an exact unbounded cylinder representation.

## Parabolic-cylinder source and Boolean derivative — 2026-09-24

The rank-one real case `x^2 + 2y = 0` is now classified separately from
rank-one quadrics with no transverse linear term (parallel/repeated/empty
plane families, which remain fail-closed). The supported fixed stratum has a
six-dimensional normalized coefficient tangent. A finite parabolic graph
patch is closed by translation across its transverse linear axis, avoiding
degenerate side faces from extruding along the cylinder axis.

The analytic source-motion Jacobian agrees with fixed-stratum retraction
finite differences at step sizes `1e-4`, `3e-5`, and `1e-5`; the last relative
error is `2.34e-10`. The actual CGAL Boolean-to-ten-moment Jacobian has 30
constructed intersection vertices and reaches `6.52e-10` relative error at
`1e-5`. The same reduced `10 x 6` result is carried by the public unified
moment evaluator without populating its ambient `10 x 10` derivative.

This is a finite synthetic method-validation path, not an unbounded exact
surface or paper case. It does not supply a rank-deficient optimizer, handle
stratum/topology transitions, cover other rank-one degeneracies, or reproduce
an author result.

## Exact truncation-root certificate test registration — 2026-09-24

The existing `tests/truncation_root_certificate_regression.cpp` exercised the
exact single-root certificate, rejection of multiple admissible roots,
endpoint roots, and tangent/repeated roots, but was not registered in CMake and
therefore did not run in the Faster suite. It is now a dedicated
`spin-faster-truncation-root-certificate` CTest. The test passed, and the full
Faster label excluding only the long unknown-seed sweep passed **45/45**.
This closes a test-registration gap for the current root-certificate algorithm;
it does not change solver convergence or reproduction status.

## Prescribed-waterline buoyant plane solve — 2026-09-24

Added a dedicated three-constraint/three-parameter planar solver for Faster
§4.3. It uses the same native plane Boolean/moment map and its reduced
coefficient derivative, enforces displaced-liquid total mass plus horizontal
center-of-buoyancy moments, and independently compares center height against
the supplied-waterline metacenter. A non-solution asymmetric-box start
recovers the analytic half-box plane; a separate waterline input yields a
negative margin without conflating that with optimizer failure. Iteration
budget exhaustion remains non-converged. The focused neighboring plane,
moment-problem, and hydrostatic tests passed **5/5**.

This fills a specialized application-solver path, not the general ten-variable
Eq. 14 solver, not free-waterline hydrostatics, and not a Faster author case.
The output verifies a nonsingular reduced plane-family equilibrium only; it
does not provide a continuous global infeasibility certificate. Overall
reproduction status remains `PARTIAL_BLOCKED`.

## Explicit digital cutting plane — 2026-09-24

The Spin-It digital-cut path now accepts a caller-supplied plane point, normal,
grid origin, anisotropic voxel spacing, and positive kerf width. It normalizes
the plane normal and selects cells by closed axis-aligned cell-box intersection
with the infinite kerf slab (rather than center-only sampling). The generated
cut is passed through the existing face-connected void/cavity audit; fixed
material is immutable, and the result reports the remaining discrete material,
material components, and which pre-existing cavities became exterior-connected.

The four new plane-rasterization/integration cases and five prior explicit-cell
cut-plan cases pass in the registered `spin-it-2014-discrete-cut-plan-python`
CTest. This is a discrete occupancy algorithm check only: no finite saw/tool
geometry, continuous mesh Boolean, automatic cut selection, wall-thickness
claim, or author-selected paper cut position is implemented by this addition.
It does not alter any historical optimization result or paper-case readiness.
