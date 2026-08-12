# Reproduction change manifest

This manifest describes the reproduction work relative to upstream commit
`71e08361ca255a1cf6b89de6e4c92ef33024f432` on branch `reproduction-v1`.
It does not claim that the paper has been fully reproduced.

## Modified upstream files

| File | Change | Status |
|---|---|---|
| `.gitignore` | Ignores isolated CMake/Ninja build products and user presets. | R0 implemented |
| `CMakeLists.txt` | Adds explicit core/viewer/test targets, pins the local libigl submodule, C++17 and target-local dependencies. | R0 verified by clean headless and viewer builds |
| `README.md` | Replaces the stale monolithic description with the audited scope, build instructions and known non-conformances. | R0 implemented |
| `main.cpp` | Includes `<cstdint>` before legacy libigl headers for GCC 13 compatibility. | R0 compatibility change |
| `include/inner_void_mesh.h` | Converts `InnerVoidMesh` into the GUI-compatible facade over `VoxelGrid` and `CarvingSolver`. | Phase 1.5 implemented |
| `src/inner_void_mesh.cpp` | Delegates carving and cavity extraction to the modular components. | Phase 1.5 implemented |

The temporary executable-bit changes on `src/center_of_mass.cpp`,
`src/grid_util.cpp`, and `src/make_it_stand.cpp` were removed; their contents
remain identical to the upstream baseline.

## New implementation files

| Files | Responsibility | Current evidence |
|---|---|---|
| `include/voxel/Voxel.h`, `src/voxel/Voxel.cpp` | Voxel state extracted from the original facade. | Compiles in core target |
| `include/voxel/VoxelGrid.h`, `src/voxel/VoxelGrid.cpp` | Regular-grid construction, SDF classification and corner associations. | Compiles; paper fidelity not yet verified |
| `include/voxel/VoxelSurface.h`, `src/voxel/VoxelSurface.cpp` | Shared-face triangulation and cavity orientation. | Compiles; regression evidence pending |
| `include/extraction/VoidMeshExtractor.h`, `src/extraction/VoidMeshExtractor.cpp` | Carved-cavity boundary extraction. | Compiles; printable-output claim not made |
| `include/carving/CarvingSolver.h`, `src/carving/CarvingSolver.cpp` | Carving state, R3 objective collection and inherited one-target heuristic. | R3 objective integration verified; carving remains algorithmically non-conformant until R4. |
| `include/carving/CarvingDiagnostics.h`, `src/carving/CarvingDiagnostics.cpp` | Per-step diagnostics, stop reasons and CSV serialization. | Verified by clean headless regression |
| `tests/carving_regression.cpp` | Headless Bunny diagnostic smoke test. | CTest passes; does not establish paper correctness |

## New R0 files

| File | Purpose |
|---|---|
| `CMakePresets.json` | Separates reproducible headless and viewer build trees. |
| `docs/conformance_matrix.md` | Maps paper requirements to local symbols, tests, tolerances and status. |
| `docs/reproduction_environment.md` | Locks the verified toolchain and resolved dependency commits. |
| `docs/provenance.md` | Records source hierarchy, hashes and license boundaries. |
| `docs/r0_baseline.md` | Records the pre-fix numerical baseline and R0 acceptance evidence. |

## Deliberately unchanged behavior

R0 freezes and documents the existing Phase 1 implementation. It does not fix
the known center-of-mass, energy, stopping or best-prefix carving defects. Those
changes belong to the following numerical and algorithm-conformance phases.

## Submodule state

The repository gitlink remains pinned to libigl
`f06fb13014aadff7eb1d7423ac018620441a3657`. No local libigl patch is part of
R0; compatibility includes live in this project's own translation units.

## R1 diagnostic and termination changes

| File | R1 change | Evidence |
|---|---|---|
| `include/carving/CarvingDiagnostics.h`, `src/carving/CarvingDiagnostics.cpp` | Adds the explicit `NoCandidatesAvailable` / `no_candidates` terminal state. | CSV and solver stop tokens agree in regression tests. |
| `include/carving/CarvingSolver.h`, `src/carving/CarvingSolver.cpp` | Persists whether candidate search ran and stops immediately when the resulting set is empty. | Full Bunny run records the terminal zero-candidate probe at iteration 17 and eliminates no-op iterations 18–21. |
| `tests/carving_regression.cpp` | Adds initial, one-step, full-known-baseline and synthetic no-candidate scenarios; validates state continuity, finite values and CSV schema. | Four CTest cases pass and survive ten repeated runs. |
| `CMakeLists.txt` | Registers four independently named R1 tests with separate CSV outputs and timeouts. | `ctest --preset headless` reports `4/4 Passed`. |
| `tests/baselines/bunny_r0_known_failure.csv` | Preserves the exact pre-R1 21-step divergent trace. | SHA-256 `3f9108649fc329a136a5371597335a99275021b6ddd9a3984198ce2b8e8c1856`. |
| `tests/baselines/README.md` | Defines the snapshot as historical failure evidence, not a target result. | Documentation audit. |
| `docs/r1_acceptance.md` | Records R1 scope, commands, outcomes and remaining non-conformances. | R1 acceptance record. |

R1 does not modify the legacy projected-energy formula, mass/first-moment
updates, fixed batch size or best-prefix policy. Passing these control-flow
tests therefore does not establish paper equivalence.

## R2 mass-property and projected-energy changes

| Files | R2 change | Evidence / status |
|---|---|---|
| `include/physics/MassProperties.h`, `src/physics/MassProperties.cpp` | Adds the paired signed state `(m,Q)`, derives CoM only as `Q/m`, evaluates a consistently oriented closed triangle boundary by signed tetrahedra, validates closed edge orientation, and provides exact uniform axis-aligned-box contributions. | Four analytic physics tests pass in the final clean verification. |
| `include/physics/BalanceEnergy.h`, `src/physics/BalanceEnergy.cpp` | Implements the normalized-gravity perpendicular projection, projected distance, and paper Eq. 2 energy `0.5*d^2`; rejects zero/non-finite gravity and normalizes finite non-zero directions without scale overflow/underflow. | `balance-energy-analytic` passes, including tiny and near-`DBL_MAX` direction scaling. |
| `include/center_of_mass.h`, `src/center_of_mass.cpp` | Retains the inherited API as a compatibility wrapper over `MassProperties` instead of maintaining a second mass/CoM implementation. | Covered by the signed-surface implementation and final clean suite. |
| `include/carving/CarvingSolver.h`, `src/carving/CarvingSolver.cpp` | Replaces direct CoM adjustment with exact subtraction of each carved axis-aligned voxel's mass and first moment; requires an initially filled grid and a closed, consistently edge-oriented, positive-mass outer mesh; derives all subsequent CoM states from `(m,Q)`. | Every-step comparison with a separately accumulated long-double voxel-state oracle passes. |
| `include/carving/CarvingDiagnostics.h`, `src/carving/CarvingDiagnostics.cpp` | Replaces the legacy/repeated-square fields with mass, first moment, recomputed state, normalized-gravity distance/energy, discrepancies, counts and terminal state in a 39-column CSV schema. | Per-column round-trip and invariant checks pass in all four carving regressions. |
| `tests/mass_properties_regression.cpp` | Adds `mass-surface-analytic`, `mass-voxel-analytic`, `mass-deletion-exhaustive`, and `balance-energy-analytic`. The deletion fixture enumerates all `8!` orders. | `VERIFIED_UNIT` is the maximum R2 claim for mass, first moment and Eq. 2. |
| `tests/carving_regression.cpp` | Retains initial, one-step and no-candidate contracts, replaces the hard-coded divergent full-run target with `carving-full-mass-consistency`, and checks incremental `(m,Q,c)` against the long-double oracle after every step. | `VERIFIED_REGRESSION` is the R2 CoM claim after this oracle. |
| `CMakeLists.txt` | Registers four analytic physics tests and four carving regressions as separate CTest cases. | Final clean headless suite: eight of eight tests passed. |
| `tests/baselines/README.md` | Keeps the R0/R1 divergent CSV immutable as historical defect evidence while removing its values from active assertions. | Historical artifact only; it is not an R2 target. |
| `docs/r2_acceptance.md` | Records R2 equations, preconditions, test matrix, tolerances, status ceiling and remaining non-conformances. | R2 accepted after headless, repeated, Debug, viewer and clean-state gates. |

R2 uses density one inside `CarvingSolver`, so its reported mass is signed
volume in model units rather than calibrated kilograms. The surface helper can
scale a uniform positive density, but neither heterogeneous density nor physical
units are claimed. The solver precondition is a closed, consistently
edge-oriented, outward/positive-mass outer surface; R2 does not claim to repair
open or misoriented input meshes.

R2 does not implement support polygons or safe regions, multiple targets,
paper best-prefix carving, the paper's too-few-removals stopping rule, wall-
thickness guarantees, deformation, analytic gradients, or alternating
optimization. Consequently, its unit/regression evidence is not a claim of
paper-equivalent carving, objective monotonicity, physical stability,
printability, or `VERIFIED_PAPER` status.

## R3 support and objective changes

| Files | R3 change | Evidence / status |
|---|---|---|
| `include/stability/SupportRegion.h`, `src/stability/SupportRegion.cpp` | Extracts the inclusive contact hull in a gravity-orthogonal local frame, constructs `h*tan(gamma)` safe offsets, classifies outside/feasible/safe states, projects to the closed safe set, and explicitly reports unattainable margins. | `support-polygon-analytic`, `safe-region-analytic`, and `standing-target-analytic` pass. |
| `include/physics/BalanceObjective.h`, `src/physics/BalanceObjective.cpp` | Represents normalized target/gravity pairs, suspension targets and the unweighted arbitrary-`K` projected CoM sum from paper Eq. 5. | Suspension and `K=1..4` analytic fixtures pass. |
| `include/carving/CarvingSolver.h`, `src/carving/CarvingSolver.cpp`, `include/inner_void_mesh.h`, `src/inner_void_mesh.cpp` | Makes a non-empty objective vector the single target state while preserving the legacy one-target constructor; aggregate queries consume all objectives. Incomplete multi-target carving is rejected before mutation and remains R4. | `multi-target-solver-boundary` and all four R1/R2 carving regressions pass. |
| `main.cpp` | Treats the clicked vertex as a lateral posture anchor, places the actual lowest geometry on the ground, exposes contact/toppling inputs and constructs the carving target from the standing support geometry. | Viewer target compiles in the final viewer and clean-clone gates. |
| `tests/stability_regression.cpp`, `tests/objective_integration_regression.cpp` | Adds five analytic geometry/objective scenarios and one production integration/R4-boundary scenario. | Six R3 CTests pass; active suite is `14/14`. |
| `docs/r3_acceptance.md`, `docs/conformance_matrix.md`, `README.md`, `CMakeLists.txt` | Records paper/beta differences, status ceilings, test targets and remaining phase boundaries. | Final R3 verification record. |

R3 is limited to support/target construction and the CoM term of Eq. 5. It
does not claim multi-target carving, target updates during alternation,
deformation, cord clearance, printability, physical stability or complete
paper reproduction. The current exact inset implementation also has a
documented worst-case `O(H^4)` performance ceiling in the number of support
hull edges; normal and analytic fixtures are verified, not unbounded hull
scalability.
