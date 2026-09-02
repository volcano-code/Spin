# Spin-It Reproduction Architecture

Date: 2026-09-02  
Repository: `/home/lin/xvxinyao/make-it-balance-paper-sim`  
Target branch: `spin-it-reproduction`

## 1. Objective

Extend the existing Make-It-Balance project with three independently auditable capabilities:

1. a carving-only reachability analyzer for the current static-balancing pipeline;
2. a reproduction of *Spin-It Faster* based on quadric material-air interfaces;
3. a reproduction of *Spin-It* 2014, including adaptive-voxel hollowing, cage deformation, and dual-material optimization.

The work remains in the Make-It-Balance repository. It does not change the status of the unfinished Make-It-Stand paper-simulation release gate. Static-balancing results, Spin-It results, and local extensions must remain separately labeled.

## 2. Source Scope and Claim Levels

The implementation is based on:

- Bächer et al., *Spin-It: Optimizing Moment of Inertia for Spinnable Objects*, ACM TOG 2014, DOI `10.1145/2601097.2601157`;
- Hafner et al., *Spin-It Faster: Quadrics Solve All Topology Optimization Problems That Depend Only on Mass Moments*, ACM TOG 2024, DOI `10.1145/3658194`;
- the existing Make-It-Balance mass integration, voxel, BBW, deformation, sealed-artifact, and OpenGL infrastructure.

Every experiment receives one claim label:

- `PAPER_FORMULA_REPRODUCTION`: paper equations and stated algorithm are implemented, but author source code is unavailable;
- `RECONSTRUCTED_ASSET`: one or more paper inputs were rebuilt or substituted;
- `LOCAL_EXTENSION`: behavior or validation not specified by the paper;
- `BLOCKED_MISSING_ASSET`: an experiment cannot be represented honestly with available data.

The project must not claim bitwise or author-code equivalence without the original code, configuration, and assets.

## 3. Repository and Provenance Strategy

Before implementation:

1. record the current branch, tracked diff, untracked-path list, and SHA-256 manifest;
2. commit the existing paper-simulation state as a named checkpoint after verifying the manifest;
3. create `spin-it-reproduction` in the same repository;
4. keep all generated results under `experiments/spin_it_reproduction/`;
5. bind every run to source commit, dirty diff hash, executable hash, configuration hash, dependency versions, and input hashes.

The local Git identity is `volcano <xuxinyao307@gmail.com>`. No global Git configuration is changed.

## 4. Selected Architecture

The project uses a native C++ layered design. A shared mass-moment kernel feeds three independent solver families instead of extending the current plane-carving heuristic with mode-specific conditions.

```text
OBJ + shell/domain + target frame
                |
                v
       RawMassMoments10
                |
       +--------+------------------+
       |                           |
       v                           v
Reachability LP             SpinDesignRequest
                                   |
                    +--------------+-------------+
                    |                            |
                    v                            v
          QuadricHollowingSolver       AdaptiveVoxelSpinSolver
              (Faster 2024)                (Spin-It 2014)
                                                   |
                                      +------------+----------+
                                      |                       |
                                      v                       v
                           CageSpinDeformationSolver  DualMaterialSpinSolver
```

This isolates solver semantics while reusing exact mass integration, input normalization, material geometry auditing, experiment provenance, and visualization.

## 5. Shared Data Model

### 5.1 RawMassMoments10

`RawMassMoments10` stores the ten independent raw moments in the paper ordering:

```text
(v, lx, ly, lz, Qxx, Qyy, Qzz, Qyz, Qxz, Qxy)
```

It provides:

- addition and subtraction for fixed-shell and variable-domain composition;
- uniform-density scaling;
- transformations under translation, rotation, and uniform scale;
- conversion to mass, center of mass, inertia about the origin, and inertia about the center of mass;
- finite-value and physical-validity checks;
- deterministic serialization and hashing.

The existing `MassProperties` representation remains supported. Conversion between the two forms is explicit and covered by regression tests.

### 5.2 SpinDesignRequest

`SpinDesignRequest` contains:

- immutable outer shell;
- interior design domain after applying minimum shell thickness;
- base and optional interior material densities;
- normalized target frame;
- contact point or axle point;
- requested mode;
- objective weights stated in physical units;
- numerical tolerances, maximum iterations, restart seed sequence, and resource limits.

Supported modes are:

- `STATIC_BALANCE`;
- `SPINNING_TOP`;
- `YO_YO`;
- `FIDGET_SPINNER`.

Requests are rejected before optimization if the mesh, target frame, density, unit scale, design domain, or shell relationship is invalid.

### 5.3 Shared Result Contract

All solvers return a result containing:

- explicit terminal status;
- final raw moments, center of mass, inertia tensor, eigenvalues, and principal axes;
- objective and constraint residuals evaluated independently from the optimizer;
- outer shell, material boundary, and cavity boundary;
- watertightness, winding, self-intersection, non-manifold, minimum-thickness, and mass-consistency audits;
- iteration or restart trace;
- input, configuration, geometry, and output SHA-256 hashes;
- source and claim labels.

## 6. Carving Reachability Analyzer

### 6.1 Mathematical Model

The analyzer asks whether a fixed shell plus continuously retained interior cell material can place the projected center of mass in the safe support region.

For interior cell `i`, let `alpha_i` be its retained fraction, `m_i` its mass, and `q_i` its first moment. With fixed-shell mass `m_f` and first moment `q_f`, the material moments are

```text
m = m_f + sum(alpha_i m_i)
q = q_f + sum(alpha_i q_i),  0 <= alpha_i <= 1.
```

For every safe-region halfspace `a_j^T c_xy <= b_j`, positivity of `m` yields the linear inequality

```text
a_j^T (q_f,xy + sum(alpha_i q_i,xy))
    <= b_j (m_f + sum(alpha_i m_i)).
```

The complete system is a linear feasibility problem. A fixed HiGHS version is used so infeasibility can be accompanied by a solver certificate instead of inferred from heuristic failure.

### 6.2 Status Semantics

- `ALREADY_SAFE`: the uncarved model satisfies the requested safe region;
- `PROVEN_UNREACHABLE`: the continuous relaxation is infeasible and the certificate passes independent checks;
- `RELAXATION_REACHABLE`: the relaxation is feasible; this does not prove a discrete, connected, manufacturable cavity exists;
- `INDETERMINATE`: numerical failure, invalid certificate, or resource exhaustion.

If the existing plane-carving algorithm fails after `RELAXATION_REACHABLE`, the report records `SOLVER_GAP` and proceeds. No further heuristic tuning is required before Spin-It development.

## 7. Spin-It Faster Solver

### 7.1 Representation

The interior material distribution is represented by

```text
beta(x) = a + b^T x + x^T A x,
```

where `A` is symmetric and the ten coefficients are collected in `p`. Material occupies the portion of the design domain where `beta(x) > 0`. Because positive scaling does not change the material set, `p` is constrained to the unit 9-sphere.

Fixed shell moments are computed once and added as a constant offset to all variable-domain moments.

### 7.2 Optimization

`QuadricHollowingSolver` implements the paper's optimality system:

- compute objective and constraint derivatives with respect to raw moments;
- recover the Lagrange multipliers and positive scale by orthogonal projection;
- solve the Newton system for a tangent direction on the coefficient sphere;
- update along a geodesic arc;
- halve the step until the full optimality residual decreases;
- restart from a deterministic application-aware sample if line search or local optimality checks fail.

The target equation residual uses `1e-12` on normalized unit-scale inputs. Independently recomputed physical constraints must meet a relative tolerance of at most `1e-6`.

### 7.3 Geometry and Moment Evaluation

The quadric is triangulated from a fixed icosphere through the inverse Gauss-map construction described by the paper and supplement. The triangulated quadric is intersected with the design-domain boundary using fixed CGAL/libigl Boolean dependencies. Raw moments are integrated from the sealed result using exact-degree triangle quadrature after divergence-theorem reduction.

The optimizer must not silently replace Boolean integration with voxel sampling. Geometry-evaluation failure is a first-class terminal status.

### 7.4 Application Objectives

- Static balance: prescribed horizontal center of mass with maximum-mass or minimum-height objective; expected interface is planar.
- Fidget spinner: prescribed planar center of mass with maximum polar moment; expected void interface is circular before extrusion.
- Yo-yo: center of mass at the axle, target axis a principal axis, and dominant axial inertia objective.
- Spinning top: horizontal center of mass on the gravity/spin axis, target axis a principal axis, dominant axial inertia objective, and center-of-mass lowering term.

## 8. Spin-It 2014 Solver

### 8.1 Adaptive Hollowing

`AdaptiveVoxelSpinSolver` represents the design domain with an octree. Each active cell has a continuous fill value in `[0,1]`. Boundary-overlapping cells use finest-level subcell integration against the input surface.

Optimization uses:

- the 2014 top and yo-yo objectives and equality constraints;
- a symmetric cell-neighbor Laplacian regularizer;
- an active-set sequential linear-quadratic method;
- an L-BFGS approximation instead of a dense Hessian;
- refinement of fractional cells until all fractional variables occur at maximum resolution;
- deterministic binary rounding followed by independent constraint and geometry validation.

Continuous fill output is diagnostic only and cannot be published as final printable geometry.

### 8.2 Cage Deformation

`CageSpinDeformationSolver` reuses the existing bounded biharmonic weights, linear blend skinning, and handle infrastructure. It jointly optimizes fill values and handle scale/translation parameters. The objective adds the paper's bi-Laplacian smoothness and displacement penalties.

Every candidate deformation is rejected before commit if it introduces inverted or degenerate faces, self-intersection, an invalid design domain, or a shell-thickness violation. Hollowing-only and deformation results remain separate products so silhouette changes are explicit.

### 8.3 Dual Material

`DualMaterialSpinSolver` accepts a specified boundary density and a specified denser interior material. It optimizes the interior material/void distribution with correct density-weighted raw moments. It does not introduce an unsupported continuously variable density model.

The one-material limit must reduce to the hollowing formulation. Output geometry includes material labels and separate mass-property contributions.

## 9. Failure Handling

Numerical and geometric failures are never converted into physical infeasibility. Terminal statuses include:

- `OPTIMIZATION_DID_NOT_CONVERGE`;
- `LINE_SEARCH_EXHAUSTED`;
- `GEOMETRY_EVALUATION_FAILED`;
- `INVALID_MATERIAL_ARTIFACT`;
- `RESOURCE_LIMIT_REACHED`;
- `INDETERMINATE`.

All input geometry is normalized to approximately unit scale for optimization and restored to the original frame at export. Restart seeds are explicit. Repeated runs under the same environment must produce identical event sequences and hashes.

Spin-It Faster and 2014 hollowing preserve the exact outer vertex and triangle arrays. Only the deformation solver may alter the outer geometry, and its result must state that the silhouette may change.

## 10. OpenGL Integration

The existing result viewer is extended rather than replaced. It can independently toggle:

- original outer surface;
- optimized outer surface;
- fixed shell and design domain;
- material and cavity boundaries;
- quadric surface and sign regions;
- center of mass;
- target spin axis;
- three principal inertia axes;
- support region and target point for static cases.

The viewer reads only sealed result artifacts. A real-desktop interactive smoke test is required for every solver family after its non-GUI tests pass.

## 11. Verification Strategy

### 11.1 Raw Moments

- analytic sphere, cuboid, and tetrahedron results;
- translation, rotation, scale, sum, and difference properties;
- parallel-axis theorem;
- mesh integration versus high-resolution voxel integration;
- round-trip conversion with existing `MassProperties`.

### 11.2 Reachability

- analytically reachable and unreachable fixtures;
- comparison with exhaustive enumeration on small binary domains;
- certificate validation;
- deterministic classification of the existing failed orientation experiments.

### 11.3 Spin-It Faster

- planar static interface and circular fidget-spinner interface;
- synthetic top and yo-yo with known symmetry;
- gradient and Jacobian checks;
- Boolean moment integration cross-checked against refinement studies;
- equation residual, independent physical residual, local-minimum sign, and repeatability checks.

### 11.4 Spin-It 2014

- finite-difference checks of objectives, constraints, and deformation derivatives;
- octree partition and refinement invariants;
- box-constraint and active-set behavior;
- binary output verification;
- convergence toward the corresponding Faster solution on representative mass-moment-only cases as resolution increases.

### 11.5 Extensions and End-to-End Tests

- deformation topology and surface-quality checks;
- correct two-density mass integration and one-material reduction;
- at least one case where each extension expands the hollowing-only feasible set;
- at least one analytic fixture and one complex OBJ for every supported mode;
- deterministic artifact and trace hashes;
- real-desktop OpenGL interaction.

## 12. Milestones and Gates

The work is decomposed into sequential implementation milestones:

1. preserve the current state; add shared raw moments and reachability;
2. implement Spin-It Faster for static balance, fidget spinners, yo-yos, and tops;
3. implement Spin-It 2014 adaptive hollowing;
4. implement cage deformation and dual-material optimization;
5. finish OpenGL integration, experiment matrices, and audit reports.

Every milestone must have all newly registered tests running and passing before the next begins. A physically unreachable input is a valid analyzed outcome; a missing test, numerical ambiguity, invalid geometry, or unrecorded dependency is not.

The prior Make-It-Stand release blockers remain reported but do not block these milestones. Completion of this design means the implementations and simulation evidence exist at their declared claim levels; it does not imply reproduction of undocumented author parameters or physical fabrication experiments.

