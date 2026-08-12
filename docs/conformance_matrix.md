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
| INPUT-01 | §3.1: input is a closed outer triangle mesh; outputs are deformed outer surface `M_O` and void boundary `M_I`. | `Mesh`, `InnerMesh` | Input mesh plus `VoidMeshExtractor` | Watertight/orientation fixtures; verify signed material volume between both surfaces. | `PARTIAL` / R8; only the input outer mesh and cavity display exist. |
| SUPPORT-01 | §3: support polygon is the convex hull of the lowest/contact vertices relative to gravity. | `Support` | Single `balancePoint` in viewer/solver | Analytic plane/contact/hull fixtures. | `NONCONFORMANT` / R3 |
| SUPPORT-02 | §3 and Fig. 2: safe region is defined by toppling angle and current CoM height; target is the CoM projection into it. | `Support` | None | Polygon-offset fixtures and three-state outside/feasible/safe tests. | `NOT_STARTED` / R3 |
| SUPPORT-03 | §3: selected support remains fixed during optimization; support voxels are fixed handles. | `Support`, `Handles` | None | Fixed-base deformation invariant. | `NOT_STARTED` / R5–R7 |
| SUSP-01 | §3: suspended CoM projects along gravity to a selected suspension point. | suspension config paths | None | Synthetic suspension fixture; cord validity remains user responsibility. | `NOT_STARTED` / R3 |
| MULTI-01 | §3, Eq. 5: sum projected CoM objectives for any number of bases/gravity vectors. | beta `multipleObj` is limited to two | Single target/gravity | Vector-of-targets unit tests for 1, 2 and 3 targets. | `NONCONFORMANT` / R3 |
| MASS-01 | §3 volume-integral definition and Appendix: uniform-density signed mass is evaluated over a closed oriented boundary; outward `M_O` is positive and inward `M_I` is negative. The printed coefficient in Eq. 3 is treated as a documented discrepancy (item 4 below). | `mass_functions.h`, `Optimizer::compute_COM` | `MassProperties::mass`, `computeSurfaceMassProperties`, `computeAxisAlignedBoxMassProperties` | Analytic cuboid/tetrahedron, density, winding, inward-cavity and translation/rotation/scale fixtures; `1e-12` absolute plus `1e-12` relative tolerance. | `VERIFIED_UNIT` at R2; final clean gate passed. |
| MOMENT-01 | §3, Eq. 4 and Appendix: accumulate the signed first moment `Q = integral rho*x dV` independently of the derived CoM, with the same outer/inner orientation convention. | `mass_functions.h`, `Optimizer::compute_COM` | `MassProperties::firstMoment`, signed-tetrahedron surface integration, axis-aligned-box contribution | The MASS-01 fixtures additionally check analytic `Q`, `Q'=Q+mt`, `Q'=RQ`, `Q'=s^4Q`, inward-cavity subtraction and deletion prefixes. | `VERIFIED_UNIT` at R2; final clean gate passed. |
| COM-01 | §3, Eq. 4: derive `c=Q/m`; carving must update mass and first moment as a pair rather than directly adding/subtracting CoM values. | `mass_functions.h`, `Optimizer::balanceByPlaneCarving`, `Optimizer::compute_COM` | `MassProperties::tryCenterOfMass`, `CarvingSolver::subtractVoxelMassProperties`, `CarvingSolver::recomputeCarvingState` | Analytic fixtures plus every-step comparison of incremental `m`, `Q`, and `c` against an independent long-double voxel-state oracle during one-step/full Bunny regressions. | `VERIFIED_REGRESSION` at R2 after the per-step oracle; final clean gate passed. |

## Energy and inner carving

| ID | Paper anchor and requirement | Beta oracle | Local symbol | Required verification | Status / next phase |
|---|---|---|---|---|---|
| ENERGY-01 | §3.1, Eq. 2: `0.5 * ||(c-c*)^(perp g)||²`; gravity magnitude must not affect the perpendicular projection. | `Optimizer::compute_EC` | `normalizedGravityDirection`, `projectPerpendicularToGravity`, `projectedCenterOfMassDistance`, `projectedCenterOfMassEnergy`, `CarvingSolver::getCoMEnergy` | Axis-aligned, oblique, gravity-parallel, tiny and near-`DBL_MAX` scaled-gravity fixtures; reject zero/non-finite gravity and verify `E=0.5*d²`. | `VERIFIED_UNIT` at R2; final clean gate passed. |
| WALL-01 | §4: never empty voxels within user-set distance `t_min` from `M_O`. | `BoxGrid`, hull voxel thickness | `VoxelGrid::build` | Resolution-independent distance/shell fixtures. | `PARTIAL` / R4; current one-voxel shell is only a proxy. |
| CARVE-01 | §4, Eq. 8: cutting plane through target projection, perpendicular to projected current-target direction. | `InnerMesh`, `Optimizer` | `CarvingSolver::carveStep` | Synthetic plane normal/orientation/distance tests. | `PARTIAL` / R4 |
| CARVE-02 | §4: ignore negative `d_i`, sort positive distances in decreasing order. | `InnerMesh` | `CompareVoxelDistance`, candidate construction | Deterministic ordered-candidate fixture. | `IMPLEMENTED_UNVERIFIED` / R4; pointer-address sort was corrected. |
| CARVE-03 | §4: simulate every ordered deletion, record true energy, select the globally best prefix; do not stop at transient increase. | `InnerMesh`, `Optimizer` | Fixed batch of up to 1000 voxels | Hand-computable prefix and temporary-increase fixtures; selected energy ≤ initial. | `NONCONFORMANT` / R4 |
| CARVE-04 | §3.2/§4: reset/recompute binary fills during alternating iterations; carved voxels need not stay empty. | `Optimizer` | None | Two-alternation reset fixture. | `NOT_STARTED` / R4–R7 |
| CARVE-05 | §4: repeat with recomputed CoM/cutting plane; stop when too few voxels are removed. | `InnerMesh`, `Optimizer` | `CarvingSolver::stopReason`, `carveAll` | No-candidate, too-few and converged-energy tests. | `PARTIAL` / R1 verifies immediate no-candidate termination; too-few/best-prefix behavior remains R4. |
| CARVE-06 | §4: for multiple bases sort by the sum of positive per-target distances. | beta behavior differs for its two-target path | None | 2/3-target ordering fixtures governed by paper formula. | `NOT_STARTED` / R4 |
| CARVE-07 | §4: ordered plane carving avoids floating components. | `InnerMesh` | No explicit connectivity assertion | Connected-material/cavity topology test after each accepted prefix. | `IMPLEMENTED_UNVERIFIED` / R4 |

## Deformation and alternating optimization

| ID | Paper anchor and requirement | Beta oracle | Local symbol | Required verification | Status / next phase |
|---|---|---|---|---|---|
| BBW-01 | §5, Eq. 9: bounded biharmonic weights on initial voxel grid; partition of unity and handle interpolation. | `BoxGrid::computeBBW`, MOSEK QP | None | Compare open-solver weights with external official `.bbw` fixture; bounds/sum/interpolation tests. | `NOT_STARTED` / R5 |
| LBS-01 | §5, Eq. 9–10: deform outer mesh, inner mesh and grid by LBS with fixed connectivity. | `Deformable`, `Handles` | None | Identity, translation, scale and interpolation tests. | `NOT_STARTED` / R5 |
| HANDLE-01 | §5: optimizer uses 3 translations + 1 uniform scale per handle; rotation is user-controlled; support is fixed. | `Handles`, `Optimizer` | None | Parameter-count and fixed-support tests. | `NOT_STARTED` / R5 |
| LAPLACE-01 | §5.1, Eq. 12–13: outer-surface Laplacian editing PSD quadratic with implicit rotations. | `Deformable`, `Optimizer` | None | PSD/symmetry and translation/uniform-scale invariance tests. | `NOT_STARTED` / R6 |
| GRAD-01 | §5.2, Eq. 14 and Appendix: analytic mass, CoM and deformation gradients. | `Optimizer`, `mass_functions.h` | None | Finite-difference or complex-step relative error ≤ `1e-5`. | `NOT_STARTED` / R6 |
| STEP-01 | §5.2–5.3: rollback insufficient decrease; step `1.0`, ratio `0.8`, minimum `0.4`, reduce `mu` by `0.05`. | `Optimizer`; beta decrease threshold differs from paper | None | Deterministic line-search state-machine test; paper's 1% threshold is normative. | `NOT_STARTED` / R6 |
| ALT-01 | §3.2/Fig. 3: block-coordinate alternation between carving and deformation; accepted steps decrease the objective. | `Optimizer` | `make_it_stand()` placeholder | End-to-end energy trace with fresh CoM/target/grid state each block. | `NOT_STARTED` / R7 |

## Diagnostics, outputs and paper experiments

| ID | Requirement | Local symbol | Required verification | Status / next phase |
|---|---|---|---|---|
| DIAG-01 | Per-step auditable CoM, first moment, projected distance/energy, incremental/recomputed discrepancies, mass state, counts and stop reason. | `CarvingDiagnostics`, 39-column CSV writer | Exact schema and row count; per-column round trip; finite values; `E=0.5*d²`; stop-token and state-continuity checks. | `VERIFIED_REGRESSION` at R2; final clean gate passed. Diagnostic scope only. |
| TEST-01 | Zero/one-step state and diagnostic contracts. | `carving-initial-state`, `carving-one-step` | Clean configure/build/CTest; zero-step header only, one-step continuity/count checks and R2 mass-property oracle checks. | `VERIFIED_REGRESSION` at R1 and extended at R2; final clean gate passed. |
| TEST-02 | Historical deterministic full-run characterization. | `tests/baselines/bunny_r0_known_failure.csv`; retired `carving-full-known-baseline` scenario | Preserve the R0/R1 defect artifact and hash, but do not use its divergent floating-point values as an executable target. | `VERIFIED_REGRESSION` at R1 only; retired from the active R2 suite and retained as historical evidence. |
| TEST-03 | Initial empty-candidate safety contract. | `carving-no-candidates` | One search, zero removals, unchanged mass/first moment/CoM/energy, idempotent terminal `carveAll`, and `no_candidates` reason. | `VERIFIED_REGRESSION` at R1 and extended at R2; final clean gate passed. |
| TEST-04 | Analytic R2 mass-property and Eq. 2 contracts. | `mass-surface-analytic`, `mass-voxel-analytic`, `mass-deletion-exhaustive`, `balance-energy-analytic` | Four focused tests covering signed closed surfaces, axis-aligned boxes, all `8!` deletion orders, and normalized-gravity energy. | `VERIFIED_UNIT` at R2; final clean gate passed. |
| TEST-05 | R2 full-run mass-property consistency. | `carving-full-mass-consistency` plus the other three carving regressions | At every solver step, compare incremental mass, first moment and derived CoM against a separately accumulated long-double voxel-state oracle; validate the 39-column CSV. | `VERIFIED_REGRESSION` at R2; final clean gate passed. This is not paper-carving verification. |
| OUTPUT-01 | §3.1/§4: emit coherent `M_O` and `M_I`; printing preparation remains a separate engineering layer. | Input `M_O`, `VoidMeshExtractor` | Signed-volume, watertightness, orientation, wall checks and slicer validation. | `PARTIAL` / R8 |
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
3. The paper states that a deformation step must decrease energy by at least
   1%; the inspected beta path uses a different threshold. The paper value is
   normative and the choice must remain visible in configuration and tests.
4. In the printed Eq. 3, applying the coefficient `rho/6` to the symmetric raw
   cross-product sum `((vj-vi)x(vk-vi)) dot (vi+vj+vk)` yields three times the
   physical signed volume. The volume-integral definition, the Appendix, and
   the beta's one-coordinate divergence form agree with the standard signed
   tetrahedron contribution `rho * vi dot (vj x vk) / 6` (equivalently a
   `rho/18` coefficient for that symmetric raw-cross expression). R2 follows
   this physically consistent interpretation and records the printed-factor
   discrepancy rather than reproducing it.
