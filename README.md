# Make It Stand — Auditable Reproduction

This repository is an independent, staged reproduction of the method in
Prévost et al., *Make It Stand: Balancing Shapes for 3D Fabrication*
(ACM TOG / SIGGRAPH 2013). It starts from the simplified
[`helloxinw/make-it-balance`](https://github.com/helloxinw/make-it-balance)
carving prototype. It is not yet a complete reproduction of the paper.

## Source hierarchy

Implementation decisions follow this order:

1. the [published paper](https://igl.ethz.ch/projects/make-it-stand/make-it-stand-siggraph-2013-prevost-et-al.pdf)
   is the normative algorithm specification;
2. the authors' beta source is a read-only behavioral oracle for details that
   the paper does not fully specify;
3. this repository is an independent implementation and the only code that is
   built or distributed here.

Paper/beta discrepancies and every local implementation claim are tracked in
[`docs/conformance_matrix.md`](docs/conformance_matrix.md). Provenance and the
license boundary are recorded in [`docs/provenance.md`](docs/provenance.md).

## Current status

The current code modularizes and instruments the original carving prototype:

- `VoxelGrid` owns voxelization and volumetric state;
- `VoxelSurface` triangulates exposed voxel faces;
- `VoidMeshExtractor` extracts the carved cavity boundary;
- `CarvingSolver` owns the current single-target carving state;
- `CarvingDiagnostics` records per-step metrics and CSV traces;
- `InnerVoidMesh` remains a compatibility facade for the existing viewer;
- `carving-regression` provides deterministic initial, one-step, full-run and
  empty-candidate headless scenarios.

The R1 regression suite verifies execution, finite diagnostics, count/state
consistency and CSV round-tripping. It does **not** verify paper-equivalent
optimization. In particular, the current solver still has known
non-conformances:

- the legacy projected center-of-mass expression contains a repeated square;
- mass and first moment are not updated as a consistent pair;
- carving removes a fixed-size batch instead of selecting the best energy
  prefix over the ordered candidates;
- the old stopping ratio is based on the entire regular grid;
- the target is a single point rather than a support polygon/safe region;
- suspension, multiple bases, BBW/LBS deformation, Laplacian preservation,
  analytic gradients and deformation/carving alternation are not implemented;
- `make_it_stand()` remains a placeholder API.

These limitations are intentional baseline facts, not accepted final behavior.

## Reproducible environment

The R0 reference environment is:

- Ubuntu 24.04.3 LTS, x86-64;
- `/usr/bin/cmake` 3.28.3;
- GCC/G++ 13.3.0;
- Ninja 1.13.0;
- C++17, `RelWithDebInfo`;
- libigl `f06fb13014aadff7eb1d7423ac018620441a3657` from the
  repository submodule.

CMake 4.x is deliberately rejected because the pinned 2020 libigl download
machinery uses compatibility behavior removed in CMake 4. Exact dependency
commits and the first-config network behavior are documented in
[`docs/reproduction_environment.md`](docs/reproduction_environment.md).

## Clone, build and test

Clone with the pinned submodule:

```bash
git clone --recursive <repository-url>
cd make-it-balance
git submodule update --init --recursive
```

Configure, build and test the headless target:

```bash
/usr/bin/cmake --preset headless
/usr/bin/cmake --build --preset headless --parallel
/usr/bin/ctest --preset headless
```

The first configuration may download the dependency revisions selected by the
pinned libigl commit. Subsequent builds use those cached source trees.

Build the interactive viewer separately:

```bash
/usr/bin/cmake --preset viewer
/usr/bin/cmake --build --preset viewer --parallel
./build/viewer/make-it-stand data/bunny.off
```

The viewer build is a compile/regression checkpoint for the inherited GUI, not
evidence that the full paper method has been reproduced.

R1 registers four CTest cases:

```text
carving-initial-state
carving-one-step
carving-full-known-baseline
carving-no-candidates
```

They verify state/count continuity, finite diagnostics, the exact 26-column
CSV schema and safe termination when a completed candidate search is empty.
The full scenario deliberately characterizes the known divergent R0 numerical
state; it does not assert that the balancing objective improves.

To generate the full diagnostic trace manually after building the headless
target:

```bash
./build/headless/carving-regression \
  full-known-baseline \
  data/bunny.off \
  build/headless/carving_full.csv
```

The resulting R1 trace terminates at the first empty candidate set. The
pre-R1 21-step trace is retained in
`tests/baselines/bunny_r0_known_failure.csv` as historical evidence and still
exposes the center-of-mass drift that R2 must replace with correctness tests.

## Repository structure

```text
include/, src/
  carving/       single-target solver and diagnostics
  extraction/    carved-cavity surface extraction
  voxel/         voxel data, grid construction and exposed faces
  center_of_mass legacy signed-surface mass helper
  grid_util      alignment, transforms and regular-grid construction
  inner_void_mesh compatibility facade used by the viewer
  make_it_stand  placeholder public API
tests/
  carving_regression.cpp
  baselines/       frozen known-failure traces; never target outputs
data/
  inherited sample meshes; Bunny is used by the R1 regressions
docs/
  conformance, environment, provenance and R0/R1 acceptance records
```

`list3d.h` is header-only; there is no `list3d.cpp`.

## Verification policy

The project uses the following status levels:

`NOT_STARTED`, `PARTIAL`, `IMPLEMENTED_UNVERIFIED`, `VERIFIED_UNIT`,
`VERIFIED_REGRESSION`, `VERIFIED_PAPER`, and `NONCONFORMANT`.

No feature is described as paper-reproduced until its governing equation,
implementation symbol, fixture, metric and tolerance are all present in the
conformance matrix. A passing diagnostic regression is never promoted to a
correctness or paper-equivalence claim.

## License

This repository retains the upstream Mozilla Public License 2.0. The official
2013 beta source has a separate research-only, non-commercial license and is
not included, linked or redistributed here. Sample-mesh provenance must be
rechecked before redistributing a derived release.

## Reference

Romain Prévost, Emily Whiting, Sylvain Lefebvre, and Olga Sorkine-Hornung.
2013. *Make It Stand: Balancing Shapes for 3D Fabrication*. ACM Transactions on
Graphics 32(4), Article 81. https://doi.org/10.1145/2461912.2461957
