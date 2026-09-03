# Spin-It and Spin-It Faster Public-Complete Reproduction Design

Date: 2026-09-04  
Repository: `/home/lin/xvxinyao/make-it-balance-paper-sim`  
Branch: `spin-it-reproduction`

## 1. Approved Objective

Complete the hybrid reproduction approved by the user:

1. reproduce every algorithmic and numerical result that can be verified from
   the public Spin-It (2014) and Spin-It Faster (2024) papers, supplements and
   videos;
2. use deterministic reconstructed assets where author CAD or interaction
   histories are unavailable;
3. preserve an author-asset compatibility layer so unpublished inputs can be
   substituted later without changing solver semantics;
4. never call reconstructed results author-exact or physical validation.

The release target is `PUBLIC_COMPLETE_REPRODUCTION`.  The independent flags
`author_exact` and `physical_validation` remain false until the corresponding
external evidence exists.

## 2. Evidence and Claim Model

Every case and aggregate gate carries all of these independent fields:

- `formula_coverage`: whether the published equations and algorithm are
  implemented and tested;
- `experiment_coverage`: whether the paper's numerical protocol has been run;
- `asset_class`: `AUTHOR_ASSET`, `RECONSTRUCTED_ASSET`, or
  `BLOCKED_MISSING_ASSET`;
- `author_exact`: true only with author code, assets, parameters and matching
  output evidence;
- `physical_validation`: true only for an executed fabrication or physical
  experiment, never for simulation or OpenGL visualization.

The existing `passed=true` milestones remain historical reconstructed-core
gates.  New paper-complete gates must not infer completeness from them.

## 3. Spin-It Faster Completion

### 3.1 Problem coverage

Extend the shared request and moment-problem layer with `BUOYANT_STABILITY`.
The request stores liquid density, displaced-liquid volume and center of
buoyancy, waterline area moments and shell volume.  Validation rejects an
inconsistent prescribed waterline or non-positive density.

Retain both published static objectives as distinct cases:

- maximum material volume/mass;
- minimum vertical center of mass.

The existing fidget-spinner, yo-yo and spinning-top objectives remain separate
and are checked against independently evaluated moment constraints.

### 3.2 Numerical method fidelity

The default paper profile uses the published Eq. 14 spherical Newton method,
geodesic updates, `alpha = 1/2`, `t_max = 1/2`, residual threshold `1e-12`,
application-aware random initial guesses, 60 function evaluations per attempt,
and restart on non-convergence or non-minimum `mu <= 0`.

Production results must use exact mesh Boolean moment integration and inverse
Gauss-map quadric triangulation.  Analytic oracles and voxel estimates may be
used for tests only and receive separate labels.

### 3.3 Published experiment matrix

Reconstruct and run the eight Table 1 examples:

1. static tilted vase;
2. floating tea infuser;
3. floating island, cat side;
4. floating island, dog and human side;
5. heart spinning top;
6. noodles spinning top;
7. jellyfish yo-yo;
8. flower fidget spinner.

Each baseline is run 50 times with recorded initial points, restarts, Newton
iterations, function evaluations, residual time and derivative time.  Geometry
counts and timing are reported as local reconstructed results, not copied paper
numbers.

Also reproduce:

- quadratic-convergence plots grouped by application;
- the static mass-versus-center-of-mass Pareto front;
- the spinning-top inertia-versus-center-of-mass Pareto front, including
  quadric type transitions;
- the nine near-infeasible tilted-cube cases;
- the two one-dimensional degenerate examples;
- analytic interface-type checks for planes, circles/elliptic cylinders,
  ellipsoids and one-/two-sheeted hyperboloids.

## 4. Spin-It 2014 Completion

### 4.1 Solver robustness

Preserve the published adaptive-octree, continuous fill, active-set SLQP,
L-BFGS, Laplacian, cage deformation and dual-material formulation.  Replace
the current small-fixture milestone with a resolution ladder that reaches the
paper's maximum refinement level nine where memory permits.  Every refinement
level records relaxation residual, binary projection residual, topology,
runtime and peak memory.

Failed high-resolution runs remain explicit failures.  The gate cannot be
passed by reducing the resolution, loosening the equality tolerance or
substituting a Faster quadric result.

### 4.2 Published experiment matrix

Build deterministic reconstructed counterparts for the public paper examples:

- asymmetric ellipsoid;
- heart;
- Armadillo spinning on hand and shell;
- teapot inertia-only, lowering-only and mass-reduction variants;
- cuboid and woven-ring yo-yos;
- dancing couple with hollowing-only and cage deformation;
- Sumone dual-density example;
- T-Rex static balance comparison;
- elephant and the paper's void/cut accounting set.

Reproduce the reported qualitative and quantitative checks that are possible
digitally:

- principal-axis alignment and center-of-mass constraints;
- the teapot angular-velocity reduction ratios;
- T-Rex normalized center-of-mass height and simulated incline threshold;
- number of connected cavities and accessibility cuts;
- hollowing-only versus cage-deformed feasibility;
- one-material reduction and dual-density benefit;
- octree-resolution convergence and runtime statistics.

Printer-specific photographs, manual spin duration and original cut placement
are recorded as `PHYSICAL_VALIDATION_MISSING` unless actually performed.

## 5. Asset Strategy

All reconstructed inputs are generated from versioned parametric scripts or
legally redistributable meshes.  Each asset manifest stores license/source,
units, coordinate frame, shell thickness, target axis, density, construction
parameters and SHA-256.

An `author_assets/` lookup path is supported but ignored by Git.  If matching
author assets later appear, the same experiment schema can run them under
`AUTHOR_ASSET` without changing objectives or tolerances.

## 6. Experiment and Artifact Contract

Every run emits:

- source commit and dirty-diff hash;
- executable, dependency, input and configuration hashes;
- full restart/iteration trace;
- normalized and physical residuals;
- sealed outer, material and cavity meshes;
- independent mass-moment and topology audit;
- deterministic output hash excluding wall-clock measurements;
- evidence and validation flags.

OpenGL is a visualization check only.  It must display outer shell, design
domain, retained material, cavity, target axis, center of mass and principal
axes on the real desktop, while preserving the saved result.

## 7. Verification and Release Gates

### Spin-It Faster gate

- all five problem families pass derivative and analytic-oracle tests;
- all eight reconstructed baselines converge in two deterministic repetitions;
- the 50-start statistics, convergence plots, Pareto fronts, feasibility and
  degenerate experiments are generated from machine-readable traces;
- Eq. 14 residual is at most `1e-12` and independent physical relative
  residual is at most `1e-6` for successful normalized cases;
- sealed geometry and moment audits pass.

### Spin-It 2014 gate

- all paper modules pass analytic and finite-difference tests;
- all reconstructed example families have a result or an explicit blocker;
- required high-resolution cases converge without tolerance weakening;
- binary material constraints and sealed geometry audits pass;
- the qualitative ordering of hollowing, deformation and dual-material
  variants matches the paper;
- repeated deterministic runs have identical event and output hashes.

### Aggregate release

`PUBLIC_COMPLETE_REPRODUCTION` is allowed only when both computational gates
pass and the full registered test suite actually runs with no unexpected skip
or `Not Run`.  The final report lists paper claims and local measurements side
by side.  `author_exact=false` and `physical_validation=false` remain visible
at the top level unless external evidence changes them.

