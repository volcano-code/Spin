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

- `VoxelGrid` owns configurable-resolution voxelization, volumetric state and
  conservative fixed-surface `t_min` eligibility;
- `VoxelSurface` triangulates exposed voxel faces;
- `VoidMeshExtractor` extracts the carved cavity boundary;
- `PlaneCarvingPolicy` implements paper Eq. (8), deterministic multi-target
  ordering and complete true-energy prefix selection;
- `CarvingSolver` owns transactional carving state and consumes arbitrary
  objective collections;
- `CarvingDiagnostics` records per-step metrics and CSV traces;
- `MassProperties` keeps signed mass and first moment as one state;
- `BalanceEnergy` implements the normalized-gravity projection in Eq. (2);
- `SupportRegion` extracts contact hulls and constructs toppling-angle-safe
  standing targets;
- `BalanceObjective` represents standing or suspension targets and sums the
  projected CoM term of Eq. (5) for any finite number of objectives;
- `InnerVoidMesh` remains a compatibility facade for the existing viewer;
- `VoxelTopology` reports six-neighbor material and carved-void components;
- `GridDeformationDomain` builds the deterministic active voxel-corner graph
  and trilinear transfer used by the Section 5 representation;
- `BoundedBiharmonicWeights` solves graph BBW with libigl's open active-set
  solver and audits bounds, partition of unity and handle interpolation;
- `HandleSet` separates locked support handles from user handles and exposes
  only three translations plus one uniform scale per free user handle;
- `LinearBlendSkinning` evaluates paper Eq. (9) directly and through the
  reduced Eq. (10) matrix;
- `DeformationSystem` evaluates `M_O`, `M_I` and the voxel lattice from the
  same immutable rest binding without mutating the R4 carving state;
- the headless regressions cover analytic policy, manufacturing constraints,
  deterministic initial/step/full behavior and empty candidates.

R2 adds analytic mass/first-moment/energy tests and verifies the Bunny solver's
incremental state against an independent full recompute after every carving
step. R3 adds paper-governed support/safe-region geometry, standing and
suspension target construction, and arbitrary-`K` Eq. (5) CoM evaluation.
R4 replaces the inherited finite-plane SDF and fixed 1000-voxel batch with
the paper's analytic cutting planes, `sum(max(0,d_ik))` multi-target score,
a complete prefix scan and the lowest computed true Eq. (5) energy under the
documented floating-point comparison. Prefix
zero remains eligible, and temporary energy increases do not end the scan.
The paper leaves “too few voxels” unspecified; the local configurable
threshold defaults to one. R4 also adds a model-space wall-thickness input and
component-count diagnostics.

R5 adds the paper Section 5 deformation representation. It binds fixed support
contacts and user handles to an active voxel-corner graph, computes bounded
biharmonic weights with an open solver, transfers those weights to arbitrary
outer-surface points by trilinear interpolation, and evaluates synchronized
outer, cavity and grid snapshots by linear blend skinning. The direct Eq. (9)
evaluation and reduced Eq. (10) matrix are checked against one another. This
stage only supplies the deformation state and parameterization: it does not
yet choose an energy-minimizing handle update.

The active graph contains corners incident to object-intersecting cells under
the documented signed-distance/half-diagonal construction. Outer-mesh points
and carved-cavity corners inside that domain use the solved BBW field exactly.
For a complete regular-grid snapshot, global corners outside the active domain
receive the weights of the nearest active node; this is a deterministic
compatibility extension for non-object grid storage, not an additional BBW
claim and not a field used to bind `M_O` or the cavity.

The viewer now constructs a standing target from the actual lowest contacts
when its flatten-band input is zero, rather than treating the clicked posture
anchor as the target. A positive band previews the support that a later
flattening stage would create; R3 does not mutate the mesh. These stages still
do not establish the complete paper optimizer. Remaining known
non-conformances include:

- target refresh and voxel refill across carving/deformation alternations are
  not scheduled yet;
- the Laplacian preservation energy, analytic mass/CoM/deformation gradients,
  deformation optimizer and rollback/step-size policy remain R6;
- carving/deformation alternation and target/grid refresh remain R7;
- final deformed-wall, connectivity and fabrication checks remain R8;
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

The active headless suite registers 33 tests through R5:

```text
mass-surface-analytic
mass-voxel-analytic
mass-deletion-exhaustive
balance-energy-analytic
support-polygon-analytic
safe-region-analytic
standing-target-analytic
suspension-target-analytic
multi-target-energy-analytic
multi-target-solver-boundary
cutting-plane-analytic
multi-plane-score-analytic
carving-order-deterministic
global-prefix-transient
global-prefix-middle
global-prefix-zero
global-prefix-shortest-tie
global-prefix-scale-invariant
global-prefix-over-1000
minimum-wall-thickness-analytic
voxel-topology-report
carving-initial-state
carving-one-step
carving-full-mass-consistency
carving-no-candidates
bbw-path-analytic
bbw-active-bound-analytic
bbw-contract-analytic
handle-parameter-contract
lbs-equation-9-10
grid-weight-trilinear
deformation-geometry-sync
deformation-invalid-inputs
```

They cover closed oriented surfaces, analytic boxes, all `8!` orders of an
eight-voxel deletion fixture, oblique/scaled gravity, state/count continuity,
the exact 44-column CSV schema, per-step incremental/full consistency and safe
termination on an empty candidate search. The R5 additions cover analytic
bounded-BBW solutions and feasibility, handle contracts, Eq. (9)/(10)
equivalence, trilinear weight transfer, synchronized geometry and invalid
state rejection. Density is fixed to one, so reported mass is a signed-volume
proxy in model units, not a calibrated value in kg.

To generate the full diagnostic trace manually after building the headless
target:

```bash
./build/headless/carving-regression \
  full-mass-consistency \
  data/bunny.off \
  build/headless/carving_full.csv
```

The historical pre-R1 trace remains in
`tests/baselines/bunny_r0_known_failure.csv`; R2 preserves it as defect evidence
but no longer treats its divergent numbers as an executable target.

## Repository structure

```text
include/, src/
  carving/       R4 plane policy, solver state and diagnostics
  deformation/   R5 graph BBW, handles, Eq. (9)/(10) LBS and synced snapshots
  extraction/    carved-cavity surface extraction
  physics/       signed mass properties and projected CoM objectives
  stability/     standing support, safe-region and target construction
  voxel/         grid, wall eligibility, topology reports and exposed faces
  center_of_mass compatibility wrapper over signed mass properties
  grid_util      alignment, transforms and regular-grid construction
  inner_void_mesh compatibility facade used by the viewer
  make_it_stand  placeholder public API
tests/
  carving_regression.cpp
  deformation_regression.cpp
  mass_properties_regression.cpp
  objective_integration_regression.cpp
  plane_carving_policy_regression.cpp
  stability_regression.cpp
  voxel_constraints_regression.cpp
  baselines/       frozen known-failure traces; never target outputs
data/
  inherited sample meshes; Bunny is used by the carving regressions
docs/
  conformance, environment, provenance and staged acceptance records
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
