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
- `CarvingSolver` owns transactional rest/deformed-lattice carving state and
  accepts target-context updates between committed steps;
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
- `LaplacianEditingEnergy` builds the Section 5 implicit-similarity residual
  operator and its positive-semidefinite quadratic form;
- `SurfaceMassDifferential` analytically differentiates signed mass, first
  moment and the outer-minus-cavity center of mass, including closed oriented
  voxel triangle chains whose components touch along an edge;
- `DeformationObjective` composes the multi-target CoM and shape terms over
  the free handle parameters while keeping the R6 voxel fill and targets fixed;
- `DeformationOptimizer` implements transactional gradient trials with the
  paper's 1% decrease rule, step schedule and `mu` schedule;
- `DeformedVoxelMassField` integrates exact mass, first moment and volume
  centroid for every active synchronously deformed hexahedral cell;
- `BalanceTargetSchedule` preserves fixed support/attachment specifications,
  refreshes standing targets from committed CoM states and applies all-of
  terminal semantics;
- `AlternatingBalanceOptimizer` executes all-filled carving, fixed-`alpha`
  deformation and recarving blocks with transactional traces and explicit
  terminal states;
- `ManufacturingGeometry` fail-closes on final outer/cavity boundary topology,
  orientation, self-intersection, containment, wall clearance, deformed-cell
  validity, material connectivity and signed-mass consistency before export;
- `SlicerValidation` provides an opt-in, hash/version/profile-pinned local
  Prusa-style slicer smoke-test interface with bounded process execution and
  auditable G-code hashes;
- `ExternalPartitionManifest` validates a user-authored split declaration,
  contained artifact paths, hashes and rigid assembly transforms without
  pretending to design the split;
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
stage supplies the deformation state and parameterization consumed by R6.

R6 adds a fixed-`alpha` deformation block for paper Section 5. It constructs
the implicit-similarity Laplacian editing energy, differentiates signed surface
mass, first moment and CoM analytically, maps those derivatives through the
free translation/scale LBS parameters, and evaluates the weighted Eq. (11)
objective for arbitrary finite target collections. The transactional optimizer
uses the paper's 1% acceptance boundary, the `1.0, 0.8, ..., 0.4` step policy
and `mu -= 0.05`; infeasible physical trials roll back while internal errors
propagate. The cavity is compacted to its referenced boundary corners for the
objective, with its mass and CoM checked against the complete synchronized R5
snapshot. This fixed-fill machinery is the foundation consumed by R7.

R7 adds the Section 3.2/Figure 3 block-coordinate controller. Each carving
block starts from an all-filled grid at the current deformation, refreshes the
standing target from that solid CoM, and refreshes it again after every
committed carving iteration. Each following deformation block freezes the
complete voxel mask: line-search trials neither refill nor recarve. A new
carving block begins only after the fixed-mask gradient has converged, and may
therefore refill previously carved voxels before selecting a different
pattern. Deformed carving uses exact signed mass and first moment for each
piecewise-linear hexahedral cell rather than rest-grid box approximations.
The controller distinguishes stable, stalled, `mu`-exhausted and local
iteration-limit terminals; no limit is reported as stability.

The production integration fixture completes Carving--Deformation--Carving,
refilling all 304 voxels carved by the first block before producing a second
pattern at Hamming distance 48. Its trace contains 27 deformation trials,
including accepted commits and complete rejected-trial rollback. Repeated
runs are deterministic, the support stays fixed, and final signed mass,
first moment and CoM agree with an independent outer-plus-cavity surface
oracle. These are local `VERIFIED_REGRESSION` results, not official asset or
`VERIFIED_PAPER` alignment.

The active graph contains corners incident to object-intersecting cells under
the documented signed-distance/half-diagonal construction. Outer-mesh points
and carved-cavity corners inside that domain use the solved BBW field exactly.
For a complete regular-grid snapshot, global corners outside the active domain
receive the weights of the nearest active node; this is a deterministic
compatibility extension for non-object grid storage, not an additional BBW
claim and not a field used to bind `M_O` or the cavity.

R8 adds a fail-closed final geometry gate. It accepts one outward outer shell,
zero or more inward nonnested cavity shells, a strict coherent triangle
boundary, exact closest-feature wall clearance, cavity containment, a single
six-neighbor material component and positive consistent signed material mass.
The R7-bound overload additionally checks supplied-result terminal, target,
trace, epoch, occupancy, parameter, geometry and mass consistency. This is an
internal-consistency contract, not proof that a publicly constructible result
originated from the optimizer. Only a sealed ready result can be atomically
exported as deterministic OBJ, ASCII STL or binary STL.

The optional slicer adapter executes an absolute pinned binary with a pinned
profile and exact version, never a shell or PATH lookup. Its fake-slicer tests
verify adapter behavior only; no real PrusaSlicer/profile artifact is bundled,
so real-slicer compatibility remains unverified. The fork/exec call requires a
single-threaded caller with exclusive, conventional `SIGCHLD` ownership. The
external partition interface records hashes and rigid transforms but does not generate cuts,
connectors, supports, drainage, coverage or assembly feasibility. Detailed
scope and evidence are in [`docs/r8_acceptance.md`](docs/r8_acceptance.md).

The viewer now constructs a standing target from the actual lowest contacts
when its flatten-band input is zero, rather than treating the clicked posture
anchor as the target. A positive band previews the support that a later
flattening stage would create; R3 does not mutate the mesh. These stages still
do not establish paper-asset equivalence or physical manufacturing readiness.
Remaining known non-conformances include:

- real process-specific slicing, printing, assembly and physical standing
  evidence are outside the local R8 unit-contract ceiling;
- official paper cases, metrics, ablations and performance alignment remain
  R9;
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

The active headless suite registers 100 tests through R8. The 38 pre-R7 tests
are:

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
laplacian-editing-analytic
surface-mass-gradient-analytic
deformation-objective-gradient-analytic
paper-line-search-state-machine
deformation-optimizer-end-to-end
```

R7 adds one closed-triangle-chain extension to the R6 differential tests and
16 tests carrying the `r7` label:

```text
surface-mass-triangle-chain-analytic
deformed-voxel-mass-identity-similarity
deformed-voxel-mass-cavity-sum
deformed-voxel-mass-invalid
deformed-carving-identity
deformed-carving-mass-field
deformed-carving-transaction
deformed-carving-objective-update
alternating-target-refresh
alternating-stability-terminal
alternating-multi-target-terminal
alternating-refill-reset
alternating-fixed-alpha-transaction
alternating-paper-beta-profile
alternating-controller-refill-block
alternating-controller-transaction
alternating-controller-end-to-end
```

The triangle-chain extension is labelled `r6`; the remaining 16 comprise the
R7 label. R8 adds 45 labelled tests: 26 geometry/export contracts, 18 pinned
slicer-adapter contracts, and one external-partition-interface contract. The
complete total is `38 + 1 + 16 + 45 = 100`.

They cover closed oriented surfaces, analytic boxes, all `8!` orders of an
eight-voxel deletion fixture, oblique/scaled gravity, state/count continuity,
the exact 44-column CSV schema, per-step incremental/full consistency and safe
termination on an empty candidate search. The R5 additions cover analytic
bounded-BBW solutions and feasibility, handle contracts, Eq. (9)/(10)
equivalence, trilinear weight transfer, synchronized geometry and invalid
state rejection. Density is fixed to one, so reported mass is a signed-volume
proxy in model units, not a calibrated value in kg.

The five R6 additions check the implicit-similarity residual and PSD quadratic,
analytic mass/first-moment/CoM derivatives and the complete fixed-`alpha`
objective gradient against central differences at relative tolerance `1e-5`.
They also check the exact 1% acceptance boundary, rollback of insufficient and
physically infeasible trials, the clamped step/`mu` state machine, synchronized
cavity mass/CoM, fixed support and immutable voxel fill in an end-to-end
optimizer fixture. A change of `mu` reweights the objective, so energies across
that state transition are deliberately not compared as a descent pair.

The R7 additions check exact deformed-cell mass integration, transactional
target updates, fixed support/specifications, all-of stability, complete voxel
refill, paper fixed-`alpha` versus beta recarving semantics, per-step target
refresh, same-context energy descent, line-search rollback, repeated
recarving, deterministic reruns and the final independent surface-mass oracle.
The R8 additions check closed/manifold/oriented/intersection-free boundaries,
outer/cavity coherence, exact wall clearance, material connectivity,
deformation inversion/collapse rejection, supplied-result consistency,
transactional deterministic export, binary round trips, executable/profile/
version pins, timeout and failure classification, auto-repair diagnostics,
minimal layer-and-XY-extrusion G-code semantics, artifact hashes and contained
external partition manifests on a trusted, non-mutating local filesystem.

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
  deformation/   R5 representation, R6 objective and R7 cell-mass sidecar
  extraction/    carved-cavity surface extraction
  optimization/  R7 block-coordinate carving/deformation controller
  physics/       signed mass/CoM values, differentials and projected objectives
  stability/     standing support, safe-region and R7 target schedule
  voxel/         grid, wall eligibility, topology reports and exposed faces
  center_of_mass compatibility wrapper over signed mass properties
  grid_util      alignment, transforms and regular-grid construction
  inner_void_mesh compatibility facade used by the viewer
  manufacturing/ final geometry gate, deterministic export, slicer adapter
                 and external partition declaration audit
  make_it_stand  placeholder public API
tests/
  carving_regression.cpp
  alternating_balance_optimizer_regression.cpp
  alternating_block_semantics_regression.cpp
  balance_target_schedule_regression.cpp
  deformed_carving_solver_regression.cpp
  deformed_voxel_mass_field_regression.cpp
  deformation_optimization_regression.cpp
  deformation_regression.cpp
  mass_properties_regression.cpp
  objective_integration_regression.cpp
  manufacturing_geometry_regression.cpp
  slicer_validation_regression.cpp
  partition_manifest_regression.cpp
  fake_slicer.cpp
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
