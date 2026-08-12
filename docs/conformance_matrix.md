# Paper conformance matrix

Last audited: 2026-08-12 (Asia/Shanghai)

Normative source: Prévost et al., *Make It Stand: Balancing Shapes for
3D Fabrication*, ACM TOG 32(4), 2013. The official beta archive is used only
as a read-only behavioral oracle; it is not part of this repository.

## Status vocabulary

- `NOT_STARTED`: no local implementation.
- `PARTIAL`: some required behavior exists, but the requirement is incomplete.
- `IMPLEMENTED_UNVERIFIED`: code exists without the required verification.
- `VERIFIED_UNIT`: focused analytic/invariant tests pass.
- `VERIFIED_REGRESSION`: deterministic end-to-end regression evidence passes.
- `VERIFIED_PAPER`: paper cases, metrics and expected trends have been matched.
- `NONCONFORMANT`: current behavior is known to contradict the paper.

No row may advance to `VERIFIED_PAPER` without an identified fixture, metric,
tolerance and reproducible result artifact.

## Input, support and mass properties

| ID | Paper anchor and requirement | Beta oracle | Local symbol | Required verification | Status / next phase |
|---|---|---|---|---|---|
| INPUT-01 | §3.1: input is a closed outer triangle mesh; outputs are deformed outer surface `M_O` and void boundary `M_I`. | `Mesh`, `InnerMesh` | Input mesh, `VoidMeshExtractor`, `DeformationSystem::evaluate` | Synchronized rest-binding/connectivity tests at R5; watertight/orientation fixtures and signed material volume between both final surfaces remain required. | `PARTIAL` / R8; R5 emits synchronized `M_O`/`M_I` snapshots but does not certify their final manufacturing geometry. |
| SUPPORT-01 | §3: support polygon is the convex hull of the lowest/contact vertices relative to gravity. | `Support` | `SupportRegion::fromMesh`; viewer contact-threshold input | Analytic contact/hull fixtures including inclusive threshold, duplicates, ordering, rigid transforms, gravity scaling, large translations and invalid/degenerate inputs. | `VERIFIED_UNIT` at R3; the geometry is fixed during one solver construction. |
| SUPPORT-02 | §3 and Fig. 2: safe region is defined by toppling angle and current CoM height; target is the CoM projection into it. | `Support` | `makeStandingTarget`, `StandingTarget`, viewer standing-target construction | Square/triangle/rectangle offsets; outside/feasible/safe states; closed boundaries; nearest edge/vertex; `gamma=0`; maximum and unattainable insets. | `VERIFIED_UNIT` at R3; target rescheduling during alternation remains R7. |
| SUPPORT-03 | §3: selected support remains fixed during optimization; support voxels are fixed handles. | `Support`, `Handles` | `makeFixedSupportHandle`, `HandleRole::Support`, `HandleSet`, `DeformationSystem` | Contact-cell corners have exact one-hot support weights; support pose is locked identity and remains fixed under user-handle updates. | `VERIFIED_UNIT` at R5 for the representation; optimizer-wide invariance during R7 alternation remains unverified. |
| SUSP-01 | §3: suspended CoM projects along gravity to a selected suspension point. | suspension config paths | `makeSuspensionObjective` | 3-4-5 analytic fixture, gravity-line zero energy, gravity scaling and invalid inputs; cord validity remains user responsibility. | `VERIFIED_UNIT` at R3; no cord-clearance claim. |
| MULTI-01 | §3, Eq. 5: sum projected CoM objectives for any number of bases/gravity vectors. | beta `multipleObj` is limited to two | `BalanceObjectives`, aggregate energy/distance, `CarvingSolver::balanceObjectives` | `K=1..4`, independent gravity scaling, permutation/duplicate behavior, empty-set rejection and production multi-target carving. | `VERIFIED_UNIT` at R3 for the CoM term and `VERIFIED_REGRESSION` at R4 for production carving. |
| MASS-01 | §3 volume-integral definition and Appendix: uniform-density signed mass is evaluated over a closed oriented boundary; outward `M_O` is positive and inward `M_I` is negative. The printed coefficient in Eq. 3 is treated as a documented discrepancy (item 4 below). | `mass_functions.h`, `Optimizer::compute_COM` | `MassProperties::mass`, `computeSurfaceMassProperties`, `computeAxisAlignedBoxMassProperties` | Analytic cuboid/tetrahedron, density, winding, inward-cavity and translation/rotation/scale fixtures; `1e-12` absolute plus `1e-12` relative tolerance. | `VERIFIED_UNIT` at R2; final clean gate passed. |
| MOMENT-01 | §3, Eq. 4 and Appendix: accumulate the signed first moment `Q = integral rho*x dV` independently of the derived CoM, with the same outer/inner orientation convention. | `mass_functions.h`, `Optimizer::compute_COM` | `MassProperties::firstMoment`, signed-tetrahedron surface integration, axis-aligned-box contribution | The MASS-01 fixtures additionally check analytic `Q`, `Q'=Q+mt`, `Q'=RQ`, `Q'=s^4Q`, inward-cavity subtraction and deletion prefixes. | `VERIFIED_UNIT` at R2; final clean gate passed. |
| COM-01 | §3, Eq. 4: derive `c=Q/m`; carving must update mass and first moment as a pair rather than directly adding/subtracting CoM values. | `mass_functions.h`, `Optimizer::balanceByPlaneCarving`, `Optimizer::compute_COM` | `MassProperties::tryCenterOfMass`, `PlaneCarvingPrefixResult::bestMassProperties`, `CarvingSolver::recomputeCarvingState` | Analytic fixtures plus every-step comparison of incremental `m`, `Q`, and `c` against an independent long-double voxel-state oracle during one-step/full Bunny regressions. | `VERIFIED_REGRESSION` at R2 after the per-step oracle; final clean gate passed. |

## Energy and inner carving

| ID | Paper anchor and requirement | Beta oracle | Local symbol | Required verification | Status / next phase |
|---|---|---|---|---|---|
| ENERGY-01 | §3.1, Eq. 2: `0.5 * ||(c-c*)^(perp g)||²`; gravity magnitude must not affect the perpendicular projection. | `Optimizer::compute_EC` | `normalizedGravityDirection`, `projectPerpendicularToGravity`, `projectedCenterOfMassDistance`, `projectedCenterOfMassEnergy`, `CarvingSolver::getCoMEnergy` | Axis-aligned, oblique, gravity-parallel, tiny and near-`DBL_MAX` scaled-gravity fixtures; reject zero/non-finite gravity and verify `E=0.5*d²`. | `VERIFIED_UNIT` at R2; final clean gate passed. |
| WALL-01 | §4: never empty voxels within user-set distance `t_min` from fixed `M_O`. | `BoxGrid`, hull voxel thickness | configurable `VoxelGrid::build`; conservative center-SDF/AABB bound | Same model-space `t_min` at resolutions 24/48; every eligible voxel is checked by an independent unit-box clearance oracle. | `VERIFIED_UNIT` at R4 for fixed-geometry voxel eligibility; physical units require calibrated input coordinates and post-deformation thickness remains R7/R8. |
| CARVE-01 | §4, Eq. 8: cutting plane through target projection, perpendicular to projected current-target direction. | `Optimizer::balanceByPlaneCarving` | `makeCuttingPlanes`, `aggregatePositivePlaneScore` | Analytic plane normal, gravity-line, half-space and satisfied-target fixtures. | `VERIFIED_UNIT` at R4. |
| CARVE-02 | §4: ignore negative `d_i`, sort positive distances in decreasing order. | beta normalizes the residual | `sortPlaneCarvingCandidates` | Analytic raw Eq. 8 score and deterministic `(score,x,y,z)` ordering. | `VERIFIED_UNIT` at R4. |
| CARVE-03 | §4: simulate every ordered deletion, record true energy, select the best prefix under the documented floating-point comparison; do not stop at transient increase. | beta scans all prefixes | `selectBestCarvingPrefix`, transactional `CarvingSolver::carveStep` | Prefix-zero, shortest tie, coordinate-scale invariance, temporary-increase, middle optimum and 1001-candidate fixtures; production energy non-increase and exception rollback. | `VERIFIED_UNIT` and `VERIFIED_REGRESSION` at R4. |
| CARVE-04 | §3.2/§4: reset/recompute binary fills during alternating iterations; carved voxels need not stay empty. | `Optimizer` | None | Two-alternation reset fixture. | `NOT_STARTED` / R7 |
| CARVE-05 | §4: repeat with recomputed CoM/cutting plane; stop when too few voxels are removed. | beta lacks the paper's repeated-stage threshold | `CarvingSolver::stopReason`, `carveAll`, configurable minimum removals | No-candidate, zero-best-prefix, iterated Bunny planning, terminal idempotence, `p+1` commit-then-stop and strict `p` threshold boundary. | `VERIFIED_REGRESSION` at R4; default threshold 1 is a documented local operationalization. |
| CARVE-06 | §4: for multiple bases sort by the sum of positive per-target distances. | beta requires both two-target distances positive | `aggregatePositivePlaneScore`, production arbitrary-`K` solver | 3-target one-positive/two-negative analytic fixture and production three-target run. | `VERIFIED_UNIT` and `VERIFIED_REGRESSION` at R4. |
| CARVE-07 | §4 prose: ordered plane carving avoids floating components. | no explicit general proof | `VoxelTopologyReport`, per-step component diagnostics | Six-neighbor analytic report and Bunny step audit. | `VERIFIED_UNIT` only for component reporting at R4; no universal no-floating guarantee is claimed. |

## Deformation and alternating optimization

| ID | Paper anchor and requirement | Beta oracle | Local symbol | Required verification | Status / next phase |
|---|---|---|---|---|---|
| BBW-01 | §5, Eq. 9: bounded biharmonic weights on initial voxel grid; partition of unity and handle interpolation. | `BoxGrid::computeBBW`, MOSEK QP | `GridDeformationDomain`, `solveBoundedBiharmonicWeights` using libigl open `active_set` | Analytic path and active-bound solutions; bounds, row sums, one-hot constraints, determinism, handle permutation, component anchoring and trilinear transfer. External official `.bbw` alignment remains a separate paper-equivalence gate. | `VERIFIED_UNIT` at R5. The official beta fixture has not been aligned, so this row is explicitly below `VERIFIED_PAPER`. |
| LBS-01 | §5, Eq. 9–10: deform outer mesh, inner mesh and grid by LBS with fixed connectivity. | `Deformable`, `Handles` | `deformLinearBlendSkinning`, `ReducedLbsModel`, `DeformationSystem::evaluate` | Identity, translation, uniform scale, fixed rotation, direct Eq. 9/reduced Eq. 10 equivalence, stale-basis rejection and synchronized `M_O`/`M_I`/grid snapshots. | `VERIFIED_UNIT` at R5; no optimization or paper-case output match is claimed. |
| HANDLE-01 | §5: optimizer uses 3 translations + 1 uniform scale per handle; rotation is user-controlled; support is fixed. | `Handles`, `Optimizer` | `Handle`, `HandleSet`, `makeFixedSupportHandle` | Exact free-parameter count/order, transactional parameter update, fixed support/contact-cell interpolation and invalid-pose rejection. | `VERIFIED_UNIT` at R5 for the parameterization; choosing parameters remains R6. |
| LAPLACE-01 | §5.1, Eq. 12–13: outer-surface Laplacian editing PSD quadratic with implicit rotations. | `Deformable`, `Optimizer` | None | PSD/symmetry and translation/uniform-scale invariance tests. | `NOT_STARTED` / R6 |
| GRAD-01 | §5.2, Eq. 14 and Appendix: analytic mass, CoM and deformation gradients. | `Optimizer`, `mass_functions.h` | None | Finite-difference or complex-step relative error ≤ `1e-5`. | `NOT_STARTED` / R6 |
| STEP-01 | §5.2–5.3: rollback insufficient decrease; step `1.0`, ratio `0.8`, minimum `0.4`, reduce `mu` by `0.05`. | `Optimizer`; beta decrease threshold differs from paper | None | Deterministic line-search state-machine test; paper's 1% threshold is normative. | `NOT_STARTED` / R6 |
| ALT-01 | §3.2/Fig. 3: block-coordinate alternation between carving and deformation; accepted steps decrease the objective. | `Optimizer` | `make_it_stand()` placeholder | End-to-end energy trace with fresh CoM/target/grid state each block. | `NOT_STARTED` / R7 |

## Diagnostics, outputs and paper experiments

| ID | Requirement | Local symbol | Required verification | Status / next phase |
|---|---|---|---|---|
| DIAG-01 | Per-step auditable CoM, first moment, projected distance/energy, prefix selection, topology, mass state, counts and stop reason. | `CarvingDiagnostics`, 44-column CSV writer | Exact schema and row count; per-column round trip; finite values; energy/prefix/topology and state-continuity checks. | `VERIFIED_REGRESSION` through R4. Diagnostic scope only. |
| TEST-01 | Zero/one-step state and diagnostic contracts. | `carving-initial-state`, `carving-one-step` | Clean configure/build/CTest; zero-step header only, one-step continuity/count checks and R2 mass-property oracle checks. | `VERIFIED_REGRESSION` at R1 and extended at R2; final clean gate passed. |
| TEST-02 | Historical deterministic full-run characterization. | `tests/baselines/bunny_r0_known_failure.csv`; retired `carving-full-known-baseline` scenario | Preserve the R0/R1 defect artifact and hash, but do not use its divergent floating-point values as an executable target. | `VERIFIED_REGRESSION` at R1 only; retired from the active R2 suite and retained as historical evidence. |
| TEST-03 | Initial empty-candidate safety contract. | `carving-no-candidates` | One search, zero removals, unchanged mass/first moment/CoM/energy, idempotent terminal `carveAll`, and `no_candidates` reason. | `VERIFIED_REGRESSION` at R1 and extended at R2; final clean gate passed. |
| TEST-04 | Analytic R2 mass-property and Eq. 2 contracts. | `mass-surface-analytic`, `mass-voxel-analytic`, `mass-deletion-exhaustive`, `balance-energy-analytic` | Four focused tests covering signed closed surfaces, axis-aligned boxes, all `8!` deletion orders, and normalized-gravity energy. | `VERIFIED_UNIT` at R2; final clean gate passed. |
| TEST-05 | Full-run mass-property and carving consistency. | `carving-full-mass-consistency` plus the other three carving regressions | At every solver step, compare `(m,Q,c)` against a long-double voxel oracle; validate non-increasing accepted energy, topology counts, determinism and the 44-column CSV. | `VERIFIED_REGRESSION` through R4. |
| TEST-06 | R3 support, target and Eq. 5 objective contracts. | Five `stability-regression` scenarios plus `multi-target-solver-boundary` | Analytic geometry/energy values, invariances, boundary/degeneracy behavior, arbitrary `K`, and R4 production carving. | `VERIFIED_UNIT` at R3 and production integration extended at R4. |
| TEST-07 | R4 fixed-geometry carving contracts. | Eleven focused R4 CTests plus Bunny single/multi-target integration | Eq. 8 planes/scores/order, all-prefix selection, shortest-prefix ties, coordinate-scale invariance, transient increase, no 1000 cap, exception rollback, model-space wall eligibility and topology reporting. | `VERIFIED_UNIT`/`VERIFIED_REGRESSION` at R4; local final gates passed. |
| TEST-08 | R5 Section 5 deformation-representation contracts. | Eight `deformation-regression` CTests | Analytic BBW path/active-bound oracles, feasibility and permutation contracts, handle parameters, Eq. 9/10 agreement, trilinear transfer, synchronized geometry and invalid-input rejection. | `VERIFIED_UNIT` at R5; Release 33/33 passed locally, while Debug/viewer/clean-clone gates are pending in the provisional R5 acceptance record. |
| OUTPUT-01 | §3.1/§4: emit coherent `M_O` and `M_I`; printing preparation remains a separate engineering layer. | Input `M_O`, `VoidMeshExtractor`, `DeformationSystem::evaluate` | R5 fixed-connectivity/synchronized-snapshot tests; signed-volume, watertightness, orientation, deformed-wall checks and slicer validation remain required. | `PARTIAL` / R8 |
| EXP-01 | §6/Table 1: reproduce optimization performance/convergence cases on modern hardware. | None | Recorded configurations, hardware, seeds, metrics and discrepancy report. | `NOT_STARTED` / R9 |
| EXP-02 | §6/Fig. 14: carving/deformation/Laplacian ablations. | None | Same inputs and comparable outcome trends. | `NOT_STARTED` / R9 |
| EXP-03 | §6: standing, suspension, multiple bases, handle noise and resolution studies. | None | Reproducible case matrix and quantitative summaries. | `NOT_STARTED` / R9 |

## Known paper/beta discrepancies

1. The paper permits any number of bases; the beta implementation exposes a
   two-object/multiple-object path. The local abstraction will follow the paper.
2. The paper's multi-base carving score sums positive `d_i` contributions. The
   beta two-target path applies a stricter joint-positive condition. The local
   implementation will use the paper formula and retain a beta-compatibility
   comparison only as an optional diagnostic.
3. Printed Equation (8) uses the raw projected residual, while the prose calls
   the value a signed distance and the beta normalizes that residual. R4 follows
   the printed equation; this affects relative multi-target ordering but not a
   one-target order.
4. The paper states that a deformation step must decrease energy by at least
   1%; the inspected beta path uses a different threshold. The paper value is
   normative and the choice must remain visible in configuration and tests.
5. In the printed Eq. 3, applying the coefficient `rho/6` to the symmetric raw
   cross-product sum `((vj-vi)x(vk-vi)) dot (vi+vj+vk)` yields three times the
   physical signed volume. The volume-integral definition, the Appendix, and
   the beta's one-coordinate divergence form agree with the standard signed
   tetrahedron contribution `rho * vi dot (vj x vk) / 6` (equivalently a
   `rho/18` coefficient for that symmetric raw-cross expression). R2 follows
   this physically consistent interpretation and records the printed-factor
   discrepancy rather than reproducing it.
6. The paper states that BBW is computed on the initial voxel grid but does not
   fully prescribe the active graph outside the object. R5 uses corners incident
   to object-intersecting cells under a signed-distance half-diagonal test. The
   object-domain and cavity weights are solved/transferred on that graph. Only
   complete-grid output corners outside it use a nearest-active-node extension;
   that compatibility convention is not treated as a paper result.
