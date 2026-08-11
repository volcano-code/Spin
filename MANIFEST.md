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
| `include/carving/CarvingSolver.h`, `src/carving/CarvingSolver.cpp` | Single-target carving state and inherited heuristic. | Implemented but algorithmically non-conformant |
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
