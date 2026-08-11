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
| MASS-01 | §3, Eq. 3: signed surface integral for mass over combined outer/inner faces. | `mass_functions.h` | `center_of_mass.cpp` | Cube, cuboid, tetrahedron, hollow shell; translation/rotation/scale invariants, relative error ≤ `1e-8`. | `IMPLEMENTED_UNVERIFIED` / R2 |
| COM-01 | §3, Eq. 4: signed first moment and CoM use outward `M_O`, inward `M_I`. | `mass_functions.h` | `center_of_mass.cpp`, legacy incremental updates in `CarvingSolver` | Analytic fixtures plus incremental/full equality for randomized deletion sequences. | `NONCONFORMANT` / R2; current incremental drift is large. |

## Energy and inner carving

| ID | Paper anchor and requirement | Beta oracle | Local symbol | Required verification | Status / next phase |
|---|---|---|---|---|---|
| ENERGY-01 | §3.1, Eq. 2: `0.5 * ||(c-c*)^(perp g)||²`, with normalized gravity. | `Optimizer` | `CarvingSolver::getCoMDistance`, diagnostic geometric metric | Analytic projection tests for horizontal/oblique gravity. | `NONCONFORMANT` / R2; legacy expression repeats a square. |
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
| DIAG-01 | Per-step auditable CoM, legacy/geometric projected-error metrics, counts, mass proxy and stop reason. | `CarvingDiagnostics`, CSV writer | Clean headless CTest, exact 26-column schema, row count, per-column round-trip, finite fields and stop-token checks. | `VERIFIED_REGRESSION` at R1; diagnostic scope only. |
| TEST-01 | Zero/one-step state and diagnostic contracts. | `carving-initial-state`, `carving-one-step` | Clean configure/build/CTest; zero-step header only and one-step continuity/count checks. | `VERIFIED_REGRESSION` at R1. |
| TEST-02 | Deterministic full-run characterization. | `carving-full-known-baseline` | Stops at first empty candidate set; frozen integer outcomes and `1e-6` numerical snapshot; repeat ten times. | `VERIFIED_REGRESSION` at R1; explicitly preserves a known-divergent baseline, not correctness. |
| TEST-03 | Initial empty-candidate safety contract. | `carving-no-candidates` | One search, zero removals, unchanged CoM/errors and `no_candidates` terminal reason. | `VERIFIED_REGRESSION` at R1. |
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
