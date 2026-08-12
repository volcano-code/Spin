# R5 deformation-representation acceptance record

Specification audit date: 2026-08-12 (Asia/Shanghai)

Acceptance status: **accepted at the documented R5 status ceiling**

The RelWithDebInfo and independent Debug headless suites pass 33/33 tests,
including all eight R5-focused tests. The viewer and committed-tree recursive
clean-clone gates also pass. This accepts the implementation and focused
verification at the unit-status ceiling described below; it is not a claim of
full-paper or fabrication reproduction.

## Authority and phase boundary

The normative source is Prévost et al., *Make It Stand: Balancing Shapes for
3D Fabrication*, especially Section 5 and Equations (9)--(10). The official
beta archive remains a read-only behavioral oracle and is not copied or linked
into this repository. Where the beta implementation and the paper differ, the
paper governs the local clean-room implementation.

R5 implements the deformation **representation**:

- a fixed-connectivity graph on the initial voxel-corner lattice;
- locked support handles and movable user handles;
- graph bounded biharmonic weights with an open active-set solver;
- trilinear transfer of grid-node weights to `M_O`;
- direct paper Eq. (9) linear blend skinning;
- the fixed reduced matrix in paper Eq. (10); and
- coherent value snapshots of deformed `M_O`, `M_I` and voxel-grid geometry.

R5 deliberately does not implement:

- the Section 5.1 Laplacian preservation energy;
- analytic mass, CoM or deformation gradients;
- optimization of the free handle parameters, line search or rollback;
- carving/deformation alternation, voxel refill or target refresh; or
- final deformed-wall, watertightness, connectivity, slicing or fabrication
  validation.

Those responsibilities remain R6, R7 and R8 respectively. Consequently R5 is
not evidence that the full optimizer converges, that an output balances, or
that the result is printable.

## Initial deformation graph

`GridDeformationDomain::fromVoxelGrid` constructs the graph once from the
initial grid. A cell is active when it is inside the input or can intersect the
surface according to the signed-distance half-diagonal test. The graph nodes
are the union of the corners incident to those active cells. Undirected edges
connect adjacent lattice corners, with deterministic local indices ordered by
their global regular-grid corner indices. Connected-component metadata is
stored and every component must be anchored by a handle constraint before BBW
is solved.

This active-domain construction is a documented local operationalization: the
paper requires weights on the initial voxel grid but does not completely
prescribe the graph outside the object. It must not be presented as an exact
recovery of an unspecified beta voxel-domain convention.

For arbitrary points contained in an active cell, the domain stores the eight
node indices and standard trilinear coefficients. `M_O` weights are transferred
with this interpolation. Carved-cavity vertices are regular-grid corners and
therefore use the corresponding graph-node weights exactly.

`DeformationSystem` also returns the complete global regular-grid corner array
for compatibility with `VoxelGrid`. A global corner outside the active object
domain receives the weight vector of its nearest active node, with lowest-index
tie breaking. This extension is used only to produce a total output snapshot;
it is not part of the solved BBW field, and `M_O` or the object-domain cavity
does not depend on it.

## Handles and fixed support

A `Handle` stores a rest pivot, a current translation, a fixed user-selected
rotation, a uniform scale, a role and constrained graph nodes. `HandleSet`
enforces the R5 parameter contract:

```text
free user handle q_j = [t_x, t_y, t_z, s]
```

Only unlocked user handles appear in the optimization vector. Rotation is
fixed in the linear basis rather than optimized. A support handle must be
locked and have identity translation, rotation and scale.

`makeFixedSupportHandle` embeds every supplied contact into an active cell and
constrains all eight corners of every such cell to the single support handle.
Because those corners receive exact one-hot Dirichlet weights, every point in
the contact cells remains fixed under trilinear interpolation while other user
handles move. Duplicate constraints are removed deterministically. Overlapping
constraints between distinct handles are rejected.

R5 proves this representation-level fixed-base invariant. It does not yet
prove that a later R7 optimization schedule never rebuilds or reassigns the
support.

## Open graph-BBW solve

For the unweighted six-neighbor graph, R5 builds the combinatorial Laplacian

```text
L_ii = valence(i),  L_ij = -1 for graph edge (i,j)
```

and minimizes each handle field with quadratic matrix `Q=L^T L`, subject to

```text
0 <= w_ij <= 1,
w_ij = 1 at handle-j constraint nodes,
w_ik = 0 there for k != j.
```

The implementation uses libigl's open `active_set` path, not proprietary
MOSEK. After all columns are solved it clamps only tolerance-scale boundary
noise, normalizes each row to partition of unity, rewrites the exact Dirichlet
rows and audits maximum bound, row-sum and constraint errors. Solver terminal
states distinguish convergence, maximum iterations and errors; a maximum-
iteration result is not silently reported as converged.

The focused analytic tests establish the local graph formulation and its
contracts. They do **not** establish numerical equivalence to the official
beta's MOSEK output. No external official `.bbw` fixture has yet been mapped
to an identical graph, constraint order and tolerance. `BBW-01` therefore has
a hard ceiling of `VERIFIED_UNIT` at R5, not `VERIFIED_PAPER`.

## Equations (9) and (10)

For rest position `x`, rest pivot `p_j`, translation `t_j`, fixed rotation
`R_j`, scale `s_j` and weight `w_j`, direct LBS evaluates

```text
x' = sum_j w_j [p_j + t_j + s_j R_j (x - p_j)].
```

`buildReducedLbsMatrix` bakes fixed rest positions, weights, pivots and
rotations into a `3n x 4N` matrix. Multiplication by the stacked per-handle
parameters `[p_j+t_j, s_j]` is the local Equation (10) form. The reduced model
rejects a handle set whose rest binding or fixed rotations differ from the
basis used to construct that matrix.

Identity, translation, uniform scale, fixed rotation and blended poses are
tested. Direct Eq. (9), raw Eq. (10), and the reduction from only free user
parameters must agree within the scale-aware test tolerance. All deformation
is reevaluated from immutable rest geometry, so repeated calls do not
accumulate incremental skinning drift.

## Synchronized geometry snapshot

`DeformationSystem` precomputes the domain, node weights, complete-corner
weights and outer-surface weights from one immutable rest binding. Its
`evaluate` method applies the current handle pose to:

- the outer rest vertices, while preserving outer face connectivity;
- the complete regular-grid corner array;
- each voxel center as the average of its eight deformed corners; and
- the current carved-cavity face connectivity over the same deformed corners.

Only voxel fill flags may differ between the current grid and the rest grid.
Lattice geometry, voxel metadata and the handle binding must match; otherwise
evaluation fails rather than combining stale states. Evaluation returns a
value snapshot and never writes deformed coordinates into `VoxelGrid`. This is
intentional at R5 because R4's mass, SDF and `t_min` data describe the original
axis-aligned grid. Safe writeback and recomputation belong to R7.

## Registered R5 test matrix

| CTest | Contract | Current result |
|---|---|---|
| `bbw-path-analytic` | Closed-form four-node path weights, exact endpoint interpolation and solver diagnostics. | Passed |
| `bbw-active-bound-analytic` | A graph whose unconstrained solution violates a bound; verifies the active bound and analytic bounded result. | Passed |
| `bbw-contract-analytic` | Bounds, row sums, one-hot constraints, determinism, constraint-order invariance and handle-column permutation. | Passed |
| `handle-parameter-contract` | Locked support identity, free `[tx,ty,tz,s]` order/count, transactionality and invalid-pose rejection. | Passed |
| `lbs-equation-9-10` | Identity, translation, scale, fixed rotation, mixed weights, global translation equivariance and Eq. (9)/(10) agreement. | Passed |
| `grid-weight-trilinear` | Deterministic active graph, exact corner/contact interpolation, trilinear coefficients and out-of-domain rejection. | Passed |
| `deformation-geometry-sync` | One rest binding drives `M_O`, `M_I`, grid corners and centers; connectivity/fill state is preserved and support stays fixed. | Passed |
| `deformation-invalid-inputs` | Empty/unanchored/overlapping constraints, invalid weights/poses/bases/grids and stale-state rejection. | Passed |

The complete RelWithDebInfo suite also re-runs all R1--R4 tests, giving 33/33
passing tests in the current working-tree verification.

## Numerical and claim boundaries

The BBW fixtures use analytic, scale-aware tolerances and explicitly audit the
maximum feasibility errors recorded in `WeightSolveDiagnostics`. The geometry
fixtures compare direct and matrix LBS at tighter floating-point tolerances.
These are deterministic unit and integration invariants on local fixtures;
they are not cross-platform bitwise baselines and do not reproduce a paper
experiment.

The following claims are permitted after the acceptance gates:

- the R5 local graph-BBW solution is bounded, partition-of-unity and
  handle-interpolating on the tested domains;
- the fixed support representation is invariant under tested user-handle
  transforms;
- the direct and reduced local LBS formulations agree; and
- `M_O`, `M_I` and the grid are evaluated coherently from one rest binding.

The following claims are not permitted:

- equality to the official beta or its proprietary-solver `.bbw` values;
- paper-level deformation quality or convergence;
- successful balancing after deformation;
- preservation of `t_min` after deformation; or
- final watertightness, structural strength, printability or fabrication
  success.

## Final verification gate

Observed on 2026-08-12:

```text
RelWithDebInfo headless configure/build                passed
RelWithDebInfo complete CTest                          33/33 passed
R5 tests, each repeated until-fail:10                  8/8 passed
Independent Debug headless CTest                       33/33 passed
RelWithDebInfo viewer configure/build                  passed
git diff --check                                       passed
Pinned libigl status                                   clean
Historical R0 trace SHA-256                            unchanged
```

The committed tree at `de3d3170f88e86b94590f54e703659ad8df800ca` was then
checked through a recursive clean clone with hard links disabled. Its headless
suite passed 33/33, its viewer target built, the parent and pinned submodule
remained clean, and the historical trace stayed byte-identical at SHA-256
`3f9108649fc329a136a5371597335a99275021b6ddd9a3984198ce2b8e8c1856`.

An external official-beta `.bbw` comparison is optional for accepting the R5
unit-status ceiling, but mandatory before considering any promotion of
`BBW-01` to `VERIFIED_PAPER`.
