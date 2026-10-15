# Faster degenerate-quadric Newton audit — 2026-09-23

## 2026-10-03 finite-vertex coverage obstruction

The eight frozen states from the local-cap diagnostic were independently
rechecked against the convex hull of ALL their existing finite vertices, using
exact represented-coordinate projection and clipping. At historical journal
index 57, 0.7948714226% of the projected design area lies outside that hull;
37.4407354062% lies outside the actual nonconvex patch. A rational separating
edge/domain-vertex witness is saved. Consequently, neither retriangulation of
the same vertices nor longer exact axial extrusion can repair this state.
Indices 18, 58 and 59 also have genuinely incomplete patches, but their vertex
hulls contain the design; that necessary condition does not validate a new
triangulation. These are historical frozen states, not the post-cap-fix seed
trajectory's indices.

No production algorithm was changed or rerun. The 5.44-second audit and three
diagnostic tests passed; they do not grant optimization or paper completion.
Evidence and scope: `.superpowers/sdd/2026-09-29-seven-day-faster-sprint/finite-patch-hull-audit/REPORT.md`.
The next source construction requires a distinct design decision, not disabling
coverage or repeating radius-only trials. Rank-two boundary mapping and its
material derivatives remain unverified.

## 2026-10-03 local far-cap certificate

Targeted replay selected seed 1 evaluations 0, 1, 2, 18, 56, 57, 58 and 59
(zero-based journal indices): the initial/first accepted state, first rejected
trial, first coverage rejection, final accepted state and its remaining trials.
An instrumented build preserved the original predicates and reproduced the
saved successful triangulation hashes and failure messages exactly. Independent
rational clipping found genuine nonzero uncovered area at 18 and 57--59.
Those trials remain invalid; increasing the radius cannot fill missing finite
source faces.

Evaluation 2 was different: all projected design area was covered, and all
existing cap-separation predicates passed at every allowed radius. The only
failure was the global far-vertex side condition. At radius 64, a remote far
vertex has axial height `-58.65581741284322`, but the actual represented far
graph clipped above the design projection has minimum height
`255.3884521723049`, versus design maximum `0.26985517352434785`.
Thus this is a sufficient-predicate false rejection, not proof of cap intrusion.

The bounded profile retains its existing fast global test and adds an exact
local sufficient certificate. A simple enclosing rim, consistently oriented
nondegenerate projected faces, cancelling interior edges and exact rim-edge
matching establish degree-one coverage. Each far triangle is clipped to the
convex projected design hull, interpolating its actual represented axial
height with rational coordinates. Positive gaps must satisfy
`gap^2 > margin^2 * axis.squaredNorm()` at every clipped vertex. All subsequent
cap-separation checks, normal sampling and radius budgets remain unchanged.
This proof concerns the represented finite graph, not complete ideal points or
ambient coefficient derivatives. Previously successful states can select an
earlier allowed radius and acquire different cap/source hashes; old scientific
results are not relabeled as current.

The frozen failed trial first failed its new native regression against the old
implementation, then passed after the change. Analytic tests reject wrong-side,
crossing, sub-margin, missing-face, folded, duplicated and invalid-index graphs;
oblique translated geometry and nonunit axes retain the clearance contract.
The independent read-only review found no actionable issue and reran the native
regression successfully. Existing distant-closure and eight-direction source
derivative tests remain in that regression.

An isolated 187-test Python run passed 186 checks. Its one error was a stale
test expecting the same seed-0 start to fail the former cap predicate: it now
returns a genuine material map. The updated test verifies that map and still
rejects solved/final-material status at its one-evaluation budget. A genuine
uncovered saved trial was added as a negative. The eight-test affected group
then passed, leaving 188 distinct current checks across the two jobs; this is
not a new whole-project green run. The intermediate assertion mistake
(`None` versus the existing empty-string success error field) is retained in
`python-02`, not removed from history.

Evidence: `.superpowers/sdd/2026-09-29-seven-day-faster-sprint/finite-cap-diagnostic/`.
No shared runtime or qualified material has been promoted.

### Unchanged-budget method outcomes after local certification

Both original seeds ran once with the corrected code, unchanged starts,
normal level, four doublings, 60 evaluations per attempt and eight attempts.
All evaluations, failures and terminal reasons are retained. Seed 0 makes
eight evaluations with two valid material maps but no accepted update:
nonpositive multiplier and zero sensitivity are rejected correctly, alongside
four coverage failures, one empty intersection and one unavailable native-side
witness. It remains `OPTIMIZATION_DID_NOT_CONVERGE`.

Seed 1 makes 67 evaluations, including 14 valid maps and 10 accepted updates.
The first attempt exhausts its 60-evaluation budget at raw residual
`0.5210635572406017`, not the required `1e-12`. Its entire ledger contains
50 coverage failures and three empty intersections; all later failed starts
remain counted. No cap-certification failure remains in either run. The
existing map auditor binds every valid evaluation to its actual Cherchi
geometry, moments and five active derivative directions. This binding is not
an assertion of generic ten-coefficient support or final-material validity.

The native regression also passes a separately recorded independent process
with byte-identical output. `finite-cap-diagnostic/collect.py` rechecks receipts,
source identities (explicitly accounting for the corrected test-only input),
fixed seed/configuration, per-attempt budgets and all successful mapping
bindings. Its `summary.json` preserves both failed solver results and the
targeted diagnostic's original identities. No further unchanged full solve is
justified: the next algorithm work concerns finite-normal coverage/ideal-map
limitations, not another arbitrary radius or restart increase. Project status
remains `PARTIAL_BLOCKED`.

## 2026-10-03 correction: finite coverage versus polygon kernel

The winding batch below recorded rejection by the then-current coverage
predicate. That rejection was **not proof of missing geometric coverage**.
Independent rational clipping of its saved source now establishes that all
134 projected source triangles have disjoint interiors, and their intersection
with the projected frozen design hull has exactly the hull's area. The old
edge-halfspace test required inclusion in the nonconvex rim's visibility
kernel, which is strictly stronger than actual inclusion.

The explicit bounded profile now checks a simple projected rim containing the
entire convex design hull, not just its vertices: all hull edges must remain
disjoint from the rim and exceed the unchanged clearance margin. Projection,
incidence and squared plane-metric distances use rationals of the represented
coordinates. A cap-separation fallback likewise tests the actual represented
cap triangle's projected convex hull, not its infinite supporting plane.
The far-cap side condition, radius budget and fixed normal selection remain
unchanged. Point/line design bounds retain compatibility; crossing a concavity,
touching, sub-margin clearance and self-crossing rims are rejected.

The saved-coefficient native regression with radius 16 / six allowed doublings
now constructs an exactly audited source over the original design vertices.
Separately, replaying **all 16 saved states with their actual original settings**
(radius 4 / four doublings, normal level 2) yields four defined source maps,
five coverage-gate rejections and seven cap-gate rejections. This replay performs
no Boolean, Newton or material approval. In particular, seed 0's first state
still fails, now at the cap guard; its failure-ledger test was updated to name
the actual failure rather than preserve the incorrect historical reason.

Evidence and isolated builds are under
`.superpowers/sdd/2026-09-29-seven-day-faster-sprint/finite-patch-coverage/`.
This correction changes admission of covered finite patches, not the source
discretization or its derivative formula. It does not supply complete ideal
points, prove convergence or grant paper reproduction.

### Verification and unchanged seeded runs

The isolated current native build passes all **187 affected Python tests**
(193.483 seconds, no skips). The native source regression passes in a second
independent process with byte-identical output. An independent rational
area-certificate job confirms the 134-triangle coverage calculation above.
`finite-patch-coverage/collect.py` rehashes current job inputs and completion
receipts, audits all attempts/evaluations and saves `summary.json`. Historical
RED runs, the original failure-message regression and its correction remain
preserved; this is not a full-project test or optimizer determinism claim.

Both frozen seeds ran with the original radius 4, four allowed doublings,
normal level 2, 60 evaluations per attempt and eight attempts:

| Seed | Evaluations | Actual outcome |
| --- | ---: | --- |
| 0 | 8 | Eight initial failures: four coverage, three cap certification, one empty Cherchi intersection; no valid material map. |
| 1 | 67 | First attempt accepts 15 updates then reaches its 60-evaluation budget; six later initial geometry failures and one nonpositive-multiplier rejection remain in the ledger. |

Seed 1 has 18 valid material evaluations, each bound to its actual Cherchi
geometry, moments and `10 x 5` active-stratum Jacobian by the existing map
auditor. Of 67 total evaluations, 29 fail coverage, 19 fail cap certification
and one produces an empty intersection. The accepted first-attempt raw
residual decreases from `0.6014658629333803` to `0.5367984577799955`;
this is **not** convergence to the unchanged `1e-12` threshold. The final
accepted step is `0.001953125`. The terminal status is
`RESOURCE_LIMIT_REACHED` (the evaluation budget, not wall time or RAM);
the last diagnostic describes a later rejected initial multiplier and must
not replace the full attempt ledger. Neither seed has verified final material.

The next bounded diagnostic should replay the first rejected trial and the
last accepted/trial states, distinguishing genuinely insufficient finite
coverage from conservative cap certification before another code change.
Do not increase radius budgets, change normal sampling or rerun these same
failed configurations without a confirmed correction. Complete ideal-point
construction, ambient coefficient directions and paper experiments remain
separate open requirements. Scientific status remains `PARTIAL_BLOCKED`.

## 2026-10-03 bounded finite-patch winding correction

The first saved seed-0 boundary-loop failure (evaluation 3) was replayed
without another optimization run. Two coarse saddle triangles disagreed with
the analytic gradient at their centroids. Individually flipping these chords
by that gradient broke the otherwise coherent finite graph boundary. The
fixed-normal paraboloid path now orients its lower graph boundary against the
positive null axis, using exact rational projected-area signs of represented
vertices. Other surface paths, normal selection, coverage and cap guards are
unchanged; a zero projected area is rejected.

The frozen optimization domain still fails finite-patch coverage. A separate
local geometry control, fixed at the uniquely derived transverse stationary
point with half-width 0.02, passes exact closed-boundary auditing. It is not a
replacement optimizer input. Source-motion aggregate Jacobian errors at
steps `3e-6, 1e-6, 3e-7` are respectively `6.34793e-7, 7.05446e-8,
6.31792e-9`. The coarser `1e-5` error `7.05328e-6` and both initial failed
test runs are retained; refinement exhibits central-difference truncation
behavior, not a relaxed scientific threshold. These aggregate norms do not
establish a per-entry bound for this new control.

The isolated native build passes 187 affected Python tests. The source
regression also passes in a second independent process with byte-identical
output. This is source-regression repeatability, not optimizer determinism.
Receipts, RED/GREEN logs, source/binary identities and a recomputable summary
are under `.superpowers/sdd/2026-09-29-seven-day-faster-sprint/finite-patch-winding/`.
Read-only review found no blocking issue in the two-file implementation/test
delta. No shared runtime, final material, full-project test status or paper
completion was promoted. Status remains `PARTIAL_BLOCKED`; finite-patch
coverage and complete singular/ideal mapping remain open.

## 2026-10-03 scope clarification: singular-A fixed-normal integration

The supplement's Section 1 explicitly uses `A^-1`, and its truncation recipe
describes hyperboloids with two remote loops. That formula does not directly
define the inverse map or ideal boundary of a rank-two paraboloid. The main
paper lists paraboloids among possible optimal surfaces; this is not evidence
that the singular-A implementation details have been recovered.

The existing finite-normal paraboloid is therefore retained as the opt-in
`BOUNDED_FINITE_NORMAL_PARABOLOID_V1` method profile. Its reduced derivative
has eight true directions; the pure-height application can transport those
into five active tangent directions only when the existing span checks pass.
Integrating this profile into a nonlinear runner does not grant complete
ideal-point construction, rank-transition support or ambient Newton coverage.

The current source integration exposed overly restrictive AABB rejection of
oblique lateral cap triangles. An alternative exact-rational projection
interval test can establish that the represented cap triangle and the convex
hull of all design vertices are disjoint. This changes neither the finite
patch nor its normal discretization. It does not bypass the patch coverage or
oriented far-cap conditions. Initial seed-0 and seed-1 evaluations still fail
finite-patch coverage and must remain failures, not infeasibility certificates.
Batch receipts and subsequent verification are under
`.superpowers/sdd/2026-09-29-seven-day-faster-sprint/fixed-normal-runner/`.

## Question

Can the existing strict full-rank Boolean/moment/Newton path directly
reproduce the paper's planar optimum for a static-balance case, and does the
paper specify a lower-dimensional Newton method that would justify replacing
Eq. (14) by a plane-only solve?

## Source finding

The official paper's Section 4.2 derives planar optimal interfaces for both
static-balance objectives. Section 5.1, however, states Eq. (14) for the
coefficient vector `p in R^n` on `S^(n-1)` and forms the Newton system using
`Jac_p r`. Section 5.2 describes its inverse-Gauss triangulation for
**non-degenerate** quadrics, and the supplement's formula uses `A^-1`; neither
source provides a separate four-parameter plane Newton system or a full
ten-coefficient derivative at the degenerate plane stratum. The paper does
not make a plane-only stationary point equivalent to the ambient Eq. (14)
system by merely restricting `p` to the plane family.

Consequently, the existing plane evaluator's `10 x 4` Jacobian is a valid
restricted-family derivative but is not by itself enough to claim that the
paper's ten-dimensional Eq. (14) has been solved. Padding it to `10 x 10`
would fabricate sensitivities and is prohibited. A plane-only Newton method
would need to be separately identified as a restricted-stratum method unless
its ambient stationarity is also independently established.

## Fixed reconstructed static-case diagnostic

One answer-independent, fixed-start diagnostic was run on the reconstructed
`static-tilted-vase` input using the strict full-rank solver and its discrete
Boolean/moment Jacobian. The supplied initial coefficients were the fixed
full-rank fixture
`(-0.31, 0.07, -0.09, 0.11, 0.83, 1.17, -0.72, 0.06, -0.04, 0.03)`;
the solver normalized them internally. The request used `DUAL_MONOMIAL_RAW_V1`,
60 iterations, 1 restart, at most 60 evaluations, Eq. (14) tolerance `1e-12`,
physical tolerance `1e-6`, icosphere subdivision 1, and truncation radius 8.
This is `RECONSTRUCTED_ASSET` method evidence, not an author case.

| Observation | Result |
|---|---:|
| Status | `OPTIMIZATION_DID_NOT_CONVERGE` |
| Terminal diagnostic | `bordered Newton system is rank deficient` |
| Function evaluations / accepted trace states | 27 / 11 |
| Raw residual, initial → last accepted | `1.1179871222` → `0.1984387135` |
| Physical residual at termination | `1.4193358578e-2` |
| Final geometry verified | no; no final artifact returned |
| Symmetric quadratic-block Frobenius norm at last accepted state | `1.1779133441e-2` |

The residual decreased while the accepted coefficient state approached a
near-planar quadric; the bordered Newton system then lost rank. This is
consistent with a degeneracy-conditioning issue, but does not prove that
degeneracy alone caused the rank loss. No plane approximation, regularization,
solver tolerance change, or four-dimensional fallback was applied.

## Reproducibility bindings

- Paper asset manifest SHA-256: `6f041346f9254fcd6141a2349039ab48a7420abaffe7c7757e1da9be2c232142`
- Table-1 request builder SHA-256: `ae11c2822b4eda9670dc061a92587514e440ae694bad60e23b998b4091931596`
- Reconstructed outer OBJ SHA-256: `8171810d7f142c158b9c9658a45ef14f02681f081ca0d534ee140a63a0dd4690`
- Reconstructed design OBJ SHA-256: `60ce2c441719197d743dc106b808260d792f1619adab0ec9a3dbe6563a2979d5`
- Native extension SHA-256: `2ec0efd0c508b2cc15af17b04547185e11dee20bc2fc0408f26ee51eead21949`
- Git commit at run: `3f930fa379bc06dcd747f99999459ea12705cff2`

## Next diagnostic before any solver change

The fixed-start run was repeated only to recover its accepted-state coefficient
vector, then the same state was reevaluated without another optimizer step. This
is a diagnostic replay, not a new success trial. At the last accepted state,
the actual discrete Boolean/material moment Jacobian and the `STATIC_BALANCE /
MAXIMUM_MASS` Eq. (14) bordered matrix were reconstructed using the same
dual-coordinate map as the solver.

The accepted native coefficient vector was
`(-0.38816061514940425, -0.9189865991007677, -0.06823414037580754,
0.0008053601163443, 0.006076811131052338, 0.00856610725702556,
-0.005271450619708028, 0.0004392875516423374, -0.0002928583677615571,
0.00021964377582116827)`. Re-evaluation returned the same geometry hash
`3de83e42c8a5aea8549aed722997213a6aacb06acbdb7fa7d0dc7ab1f267214d` and a
defined fixed-topology discrete derivative. The values below are diagnostic
linear algebra reconstructed from that derivative, not fields emitted by the
optimizer.

| Matrix diagnostic | Result |
| --- | ---: |
| KKT multipliers `(lambda_1, lambda_2, mu)` | `(2.36534345, 0.17562517, -2.57386192)` |
| Singular values of `G J_s` | `54.1816, 3.78033` |
| Singular values of the 13x13 bordered matrix, largest → smallest | `54.2427, 4.57350, 3.44490, 2.76072, 2.60075, 2.57386 (×5), 0.786173, 0.0471718, 1.36e-16` |
| Smallest / largest bordered singular value | `2.51e-18` |
| Final full-pivot elimination pivot | `1.91e-17` |
| Default Eigen `FullPivLU` relative rank threshold | `13 * epsilon` (about `2.89e-15`) |
| Effective rank at that threshold | `12 / 13` |

The row for the second first-moment constraint in `G J_s` is numerically zero
(entries at approximately `1e-16`); its other singular direction remains
well-conditioned. Thus the rank rejection is reproduced by the actual
discrete-map Jacobian and is not explained by the Boolean adapter returning a
different coordinate map. The state is close to a vertical plane. A separate
input audit found that both source meshes' vertex sets are exactly invariant
under `y -> -y` (98/98 vertices), but their triangulated face sets are not
(64/192 reflected faces match). Therefore symmetry is a plausible contributor,
not an established cause for the vanishing sensitivity. This is evidence for a
constrained degeneracy at this particular reconstructed state, not a general
failure of Eq. (14), not an infeasibility certificate, and not proof that all
planar optima have this rank defect.

The corresponding `10 x 10` moment Jacobian has only three singular values
above numerical noise at this near-plane state. This is consistent with the
loss of independent moment sensitivities as the quadric tends to a plane, but
does not alone identify which part is intrinsic to the limiting geometry and
which part is specific to this reconstructed input.

No solver behavior was changed. The next valid implementation task is to add
an explicitly typed, lower-dimensional planar stationary-system diagnostic
only if it can be tied to the paper's Eq. (14) conditions; it must not be
represented as the full ten-parameter Newton solve. Otherwise this case stays
`OPTIMIZATION_DID_NOT_CONVERGE`, and work proceeds on cylinder/paraboloid
strata and the reference-Boolean comparison. Do not reinterpret this failed
start as infeasibility or as evidence that the paper's whole algorithm fails.

## Current strict-path coverage boundary

The current source and focused CTest inventory separate implementation from
integration as follows:

| Quadric family | Current value/derivative path | Eq. (14) solver status |
|---|---|---|
| Ellipsoid, one-/two-sheet hyperboloid | Full-rank triangulation, source motion, discrete Boolean moment Jacobian | Integrated into strict full-rank Newton; synthetic end-to-end evidence exists |
| Plane | Plane triangulation plus dedicated `10 x 4` plane-family moment derivative (`PlaneMomentEvaluator`) | No restricted stationary solver; deliberately not padded into the full-rank solver |
| Elliptic cylinder, paraboloid | Surface classification and fixed-stratum tangent dimensions exist | `QuadricSourceMotion` and the strict full-rank differential fail closed; no Newton integration |

The `spin-faster-plane-family` test passes and validates the local plane-family
derivative, not optimization. The 17 focused tests run for this audit also pass;
they do not close the missing cylinder/paraboloid paths, plane stationarity, or
the reference-Boolean discrepancy. The pinned Cherchi executable remains a
separate same-operand probe: four non-ellipsoid outputs are still not exact
represented-set matches to CGAL, so that probe is not yet a drop-in paper
solver backend. These boundaries are the immediate implementation queue, not
completed paper requirements.

### Primary sources

- Hafner, Ly, and Wojtan, *Spin-It Faster*, official paper, Sections 4.2,
  5.1–5.2: <https://research-explorer.ista.ac.at/download/17203/17204/sif-final.pdf>
- Supplement, Section 1: <https://research-explorer.ista.ac.at/download/17203/17205/sif-supp-final.pdf>

## 2026-09-24 implementation update

The coverage table above records the state at the time of the 2026-09-23
audit and is superseded for local geometry derivatives by the follow-up work.
As of 2026-09-24:

| Quadric family | Current value/derivative evidence | Solver boundary |
|---|---|---|
| Plane | Native plane-family values and a restricted `10 x 4` derivative, checked by multi-step finite differences. | No full ten-dimensional Eq. (14) claim; contact and basis-branch events fail closed. |
| Simple elliptic cylinder | Source-motion and clipped Boolean-to-moment derivatives in a reduced `10 x 7` chart, checked by multi-step finite differences. | Not promoted to the full-rank Eq. (14) solver. Repeated-eigenvalue charts remain unsupported. |
| Simple-spectrum paraboloid | Source-motion and clipped Boolean-to-moment derivatives in a reduced `10 x 8` chart, including translated closure motion. | Not promoted to the full-rank Eq. (14) solver. |
| Ellipsoid / one- and two-sheet hyperboloids | Full-rank discrete Boolean/moment Jacobian and strict Newton integration. | The current frozen answer-independent five-seed run is reproducibly unsuccessful: four hit the 60-evaluation cap and one failed the initial geometry audit. |

The affected current-source CTest selection
(`spin-it-2014-dual-material*`, `spin-faster-source-motion`,
`spin-faster-discrete-moments`, `spin-faster-stratum-tangent`,
`spin-faster-plane-family`, `spin-faster-canonical-integration`, and
`spin-faster-end-to-end`) passed **8/8** on 2026-09-24. This verifies the
registered local regressions only. The current strict full-rank runner and its
five-seed outcomes are recorded in
[`experiments/spin_it_reproduction/faster/REPORT.md`](../experiments/spin_it_reproduction/faster/REPORT.md).

These updates close local derivative implementation gaps for simple cylinder
and paraboloid strata; they do not close rank-deficient optimization, prove
cross-stratum differentiability, resolve the Cherchi/CGAL represented-set
discrepancy, or add author-original case evidence. In particular, the reduced
Jacobians retain their actual dimensions and must not be padded into the
ten-dimensional Eq. (14) system. Overall status remains `PARTIAL_BLOCKED`.

## 2026-10-03: bounded fixed-normal paraboloid source increment

An explicit C++ `fixedNormalParaboloid` profile now constructs a finite
icosphere-normal patch with conservative domain coverage and oriented exterior
closure checks. Its reduced source motion retains a genuine 10×8 basis.
Four curvature/polarity controls pass exact boundary auditing and three-step
finite differences; independent source-process outputs repeat exactly.
The parameter-grid default remains unchanged. This is not complete ideal-point
construction or actual Cherchi material-derivative admission, and it is not
promoted to generic Newton. Python/request integration remains pending.

The broader projective-source regression exposed 15 multiple-crossing derivative
discrepancies; rebuilding the before-images reproduces all diagnostic lines.
This existing failure remains unresolved. See the bounded implementation and
verification [report](../.superpowers/sdd/2026-09-29-seven-day-faster-sprint/fixed-normal-paraboloid/REPORT.md).

## 2026-10-03: actual Cherchi fixed-normal material derivative

The preceding Python/material admission gaps are closed for one frozen synthetic
elliptic-paraboloid fixture. The reduced interface explicitly exposes
`fixed_normal_paraboloid`; legacy defaults remain unchanged. Source normals,
vertex identities and connectivity bind the protocol and derivative branch.
Invalid profiles, protocol downgrade, contradictory surface types, missing
directions and padded Jacobians are rejected.

Two independent processes each performed 49 actual Cherchi Boolean evaluations:
one base and eight reduced directions at three central-difference step sizes.
The material Jacobian remains genuinely 10×8. Maximum derivative relative error
was `1.3826896968215547e-8`, below `1e-6`. Scientific records, events and all 246
geometry files per process match. The separate saved-output checker reloaded all
98 rational materials, independently integrated moments and rehashed 694 output
files. Double-converted moment differences were zero on this fixture, not a
general exact-arithmetic theorem.

Final affected Python regression: 24 tests, no skips, all passed; 79.37 seconds,
226,947,072 bytes peak memory, zero swap. CMake registration was added; actual
tests used the isolated extension directly, not a fresh whole-project CTest
build. The shared runtime was not replaced. Fresh scoped review identified a
surface-type consistency defect, reproduced before repair. No scoped Important
issue remained after repair; legacy generic checker hardening remains pending.
Historical numerical runs were rechecked, not rewritten as new runs.

See the [recomputable summary](../.superpowers/sdd/2026-09-29-seven-day-faster-sprint/fixed-normal-cherchi/summary.json)
and adjacent `collect.py`, `build-01`, `regression-02`. Complete ideal-point
coverage, cross-type/full-direction Newton, the existing 15 multiple-crossing
derivative discrepancies, near-cone CAD topology failure and original paper
cases remain open. Neither manufacturing acceptance nor paper reproduction is
granted; overall status remains `PARTIAL_BLOCKED`.

## 2026-10-03: multiple-crossing failure localization and actual material check

The 15 failing native source differences are no longer unexplained. A frozen
diagnostic uses the original polynomial, normal resolution 2, effective radius
16384 and steps `1e-8, 3e-9, 1e-9`. It reproduces all 15 failures, splitting
finite-point, projective-ray and cap-centroid contributions. No production
source, original test, threshold or input was changed.

An independent arbitrary-precision implementation reconstructs finite inverse
Gauss points, ideal edge roots, straight-ray sphere cuts and cap centroids.
At both 80 and 120 decimal digits, the reported results agree:

- Native binary64 finite differences: 15 failures.
- High-precision geometry using the saved binary64 normalized coefficients:
  6 failures remain.
- Retain binary64 coefficient addition but normalize/evaluate in high
  precision: 0 failures. Entirely high-precision perturbation also has 0.
- Tiny-step reference derivatives versus the native analytic Jacobian:
  maximum aggregate relative discrepancy `1.7104e-14`; all 2,840 individual
  vertex/direction comparisons pass (`3.7238e-11` maximum nonzero relative
  discrepancy; `2.7285e-12` maximum near-zero absolute discrepancy).
- Across all 60 perturbed states, independently reconstructed geometry differs
  from native coordinates by at most `1.0942e-11` absolute and `6.6781e-16`
  scaled by each point's norm or one, whichever is larger.

This is strong numerical evidence of coefficient-normalization and geometric
evaluation/representation roundoff, not evidence for replacing the analytic
Jacobian. It is not an interval-certified error bound. The original source
finite-difference gate remains failing; this diagnostic does not silently
convert that test into a pass or relax its `1e-6` threshold.

The corresponding preregistered clipped-material diagnostic uses the same
polynomial and design cube with the existing full tangent-basis checker:
9 tangent directions plus the radial null direction, 3 steps, 61 actual Cherchi
evaluations per process. Two independent processes pass; maximum nonzero
relative error `4.1728e-7`, radial normalized absolute error `1.9931e-17`.
Scientific records/events and 306 geometry files match. A separate checker
reloads all 122 rational materials and rehashes 862 scientific output files;
independent moments agree within the unchanged `1e-10` check. This verifies one
fixed full-rank near-transition material map, not a transition/Newton solve.

The [summary and receipts](../.superpowers/sdd/2026-09-29-seven-day-faster-sprint/multiple-crossing-diagnostic/summary.json)
retain all output, decomposition and identities. Scoped independent review found
no blocking issue in the diagnostic conclusion; generalizing the parser and
near-zero global comparison requires additional checks. Exact ideal-region
coverage, cross-type optimization, near-cone CAD, original cases and final
acceptance remain outstanding. Whole-project green status is not granted.
