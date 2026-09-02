# Spin-It Faster Formula-Reproduction Implementation Plan

> Execute every task with test-driven development.  A numerical or geometry
> failure is never reclassified as physical infeasibility.

**Goal:** Reproduce the 2024 *Spin-It Faster* method for static balance,
fidget spinners, yo-yos, and spinning tops while preserving the exact input
outer surface and returning independently audited material geometry and mass
moments.

**Primary sources:** Hafner, Ly, and Wojtan, ACM TOG 43(4), Article 78,
DOI `10.1145/3658194`, plus the two-page supplemental material.  The authors'
Matlab application and CAD assets are not publicly released, so results are
limited to `PAPER_FORMULA_REPRODUCTION` or `RECONSTRUCTED_ASSET`.

**Architecture:** A strict `SpinDesignRequest` supplies a fixed shell, a
watertight interior design domain, a normalized target frame, and a mode.
`QuadricCoefficients10` represents Eq. 5 on the unit 9-sphere.  A differentiable
quadric triangulator and fixed mesh Boolean produce the material domain
`Omega intersection {beta > 0}`.  Exact surface mass integration provides raw
moments; interface quadrature provides their shape Jacobian.  A spherical
Newton solver implements Eq. 14 with deterministic restarts and fail-closed
geometry handling.

## Source-locked constants and equations

- Raw moment ordering (Eq. 3):
  `(v,l1,l2,l3,Q11,Q22,Q33,Q23,Q13,Q12)`.
- Quadric ordering (Eq. 5):
  `(a,b1,b2,b3,A11,A22,A33,A23,A13,A12)`.
- Material sign convention (Eq. 6): material iff `beta(x) > 0`.
- Optimality residual (Eq. 14):
  `grad_r f + mu*p - Jac_r(g)^T*lambda = 0`, `g=0`, `||p||=1`.
- Multipliers are recomputed by orthogonal least squares each iteration.
- Geodesic update:
  `p_next=cos(t)p+sin(t)Delta/||Delta||`.
- `alpha=1/2`, `t_max=1/2`, normalized residual threshold `1e-12`.
- Restart after at most 60 function evaluations, require `mu>0`, and reject an
  initial point whose moment Jacobian vanishes.
- The paper's moment evaluation uses a triangulated quadric, a mesh Boolean,
  one divergence-theorem reduction, and exact-degree triangle quadrature.

## Global gates

- New code is behind `MAKE_IT_STAND_BUILD_SPIN_IT_FASTER` until complete.
- Fixed shell and original outer vertex/face arrays are immutable.
- The optimizer may not silently fall back to voxel sampling.
- Every returned material product must pass watertight, winding,
  self-intersection, non-manifold, containment, minimum-wall, and
  independently recomputed moment checks.
- Normalized Eq. 14 residual is at most `1e-12`; physical constraints are at
  most `1e-6` relative to their declared scale.
- Identical runs produce identical traces and output hashes.

## Task 1: Lock primary-source provenance

**Create:**

- `experiments/spin_it_reproduction/sources/manifest.json`
- `docs/spin_it_faster_equations.md`
- `tests/python/test_spin_it_sources.py`

Write a failing test for exact DOI, URLs, SHA-256 values, equation ordering,
and claim boundary.  Record the open-access main-paper and supplement hashes;
do not commit copyrighted PDF bytes.  Commit as
`docs: lock Spin-It Faster source contract`.

## Task 2: Add strict request and result contracts

**Create:**

- `include/spin/SpinDesignRequest.h`
- `src/spin/SpinDesignRequest.cpp`
- `tests/spin_design_request_regression.cpp`

Define `SpinMode {StaticBalance, FidgetSpinner, YoYo, SpinningTop}`,
`SpinTargetFrame`, `SpinDesignRequest`, `SpinSolverSettings`, terminal status,
and `QuadricHollowingResult`.  Reject non-orthonormal/right-handed frames,
non-positive densities/scales, non-finite weights, invalid meshes, an empty
design domain, shell/domain overlap errors, and unsupported mode fields.
Hash the complete normalized request deterministically.  Commit as
`feat: add Spin-It Faster request contract`.

## Task 3: Implement four paper objective systems

**Create:**

- `include/spin/SpinMomentProblem.h`
- `src/spin/SpinMomentProblem.cpp`
- `tests/spin_moment_problem_regression.cpp`

Expose objective value, `grad_r f`, `Hess_r f`, constraint vector `g`, and
`Jac_r g` in the canonical raw-moment coordinates.  Cover:

- static maximum mass: `min -v`, `l1=l2=0`;
- optional static lower-CoM variant: `min c3`, `l1=l2=0`;
- fidget spinner: `min -(Q11+Q22)`, `l1=l2=0` in the 2D/extruded contract;
- yo-yo: paper Eq. 4.5 with `l=0`, `Ic13=Ic23=0`;
- spinning top: yo-yo objective plus `wc*c3`, with `l1=l2=0` and
  `Ic13=Ic23=0`.

Use analytic cuboids and finite-difference derivative checks.  Commit as
`feat: implement Spin-It Faster moment objectives`.

## Task 4: Implement quadric coefficient algebra

**Create:**

- `include/spin/QuadricCoefficients.h`
- `src/spin/QuadricCoefficients.cpp`
- `tests/quadric_coefficients_regression.cpp`

Implement Eq. 4/5 evaluation, gradient, homogeneous matrix conversion,
normalization on `S^9`, tangent projection, geodesic updates, rigid/similarity
coordinate transforms, sign inversion, and deterministic hashing.  Analytic
fixtures cover planes, spheres, ellipsoids, one/two-sheet hyperboloids,
paraboloids, and degenerate rejection.  Commit as
`feat: add normalized quadric representation`.

## Task 5: Pin the Boolean geometry stack

**Create/modify:**

- `cmake/SpinItFasterDependencies.cmake`
- `tests/cmake/assert_spin_it_faster_dependency.cmake`
- `CMakeLists.txt`
- `docs/reproduction_environment.md`

Enable libigl's CGAL module only behind the Faster option and lock the CGAL
commit already named by pinned libigl,
`f7c3c8212b56c0d6dae63787efc99093f4383415`.  Record GMP/MPFR versions and
prove that the disabled build neither downloads nor links them.  Add Boolean
regressions for disjoint, contained, tangent, and intersecting watertight
meshes.  Commit as `build: pin Spin-It Faster Boolean dependencies`.

## Task 6: Reproduce quadric discretization

**Create:**

- `include/spin/QuadricTriangulator.h`
- `src/spin/QuadricTriangulator.cpp`
- `tests/quadric_triangulator_regression.cpp`

Use a fixed deterministic icosphere.  Implement the supplement's inverse
Gauss map for nonsingular quadrics, ideal-normal edge intersections, radius
doubling from beyond the design-domain diameter, two-loop verification, and
truncation caps.  Add explicit plane/paraboloid handling without mislabeling
them as the nonsingular construction.  Validate normals, topology, and
convergence under icosphere refinement.  Commit as
`feat: reproduce quadric triangulation`.

## Task 7: Build differentiable material moment evaluation

**Create:**

- `include/spin/QuadricMomentEvaluator.h`
- `src/spin/QuadricMomentEvaluator.cpp`
- `tests/quadric_moment_evaluator_regression.cpp`

Boolean the design domain with the positive quadric region.  Compute raw
moments from the closed boundary using the existing exact surface integrator.
Compute `Jac_p r` from the quadric-interface shape derivative and fixed
high-order triangle quadrature; reject topology events for which the local
derivative is undefined.  Cross-check derivatives against symmetric finite
differences only as a test oracle.  Return topology/geometry hashes and never
use voxel moments.  Commit as `feat: evaluate quadric material moments`.

## Task 8: Implement spherical Newton and deterministic restarts

**Create:**

- `include/spin/QuadricHollowingSolver.h`
- `src/spin/QuadricHollowingSolver.cpp`
- `tests/quadric_hollowing_solver_regression.cpp`

Implement the least-squares multiplier recovery, Eq. 14 residual, bordered
Newton system, tangent/geodesic update, halving line search, 60-evaluation
restart boundary, deterministic mode-aware initial samples, `mu>0` local
minimum check, resource limits, and full trace.  Test quadratic convergence
and every terminal status.  Commit as
`feat: implement Spin-It Faster Newton solver`.

## Task 9: Add analytic four-mode end-to-end fixtures

**Create:**

- `tests/spin_it_faster_end_to_end.cpp`
- `experiments/spin_it_reproduction/faster/analytic_cases.json`
- `experiments/spin_it_reproduction/faster/run_analytic.py`

Require the planar static interface, circular fidget void, symmetric yo-yo,
and top interface, with physical residual `<=1e-6`, Eq. 14 residual `<=1e-12`,
immutable outer arrays, valid material geometry, and byte-identical repeats.
Commit as `test: verify Spin-It Faster analytic modes`.

## Task 10: Expose native/Python artifact API

**Create/modify:**

- `include/api/BalanceModel.h`
- `src/api/BalanceModel.cpp`
- `bindings/pybind_module.cpp`
- `python/make_it_stand/spin.py`
- `python/make_it_stand/__init__.py`
- `tests/python/test_spin_it_faster.py`

Expose strict immutable schemas for requests, traces, results, and sealed
geometry.  Results are explicitly labeled `PAPER_FORMULA_REPRODUCTION` or
`RECONSTRUCTED_ASSET`.  Commit as `feat: expose Spin-It Faster API`.

## Task 11: Run complex OBJ reproduction matrix

**Create:**

- `experiments/spin_it_reproduction/faster/cases.json`
- `experiments/spin_it_reproduction/faster/REPORT.md`
- `experiments/spin_it_reproduction/faster/release-gate.json`

Use available license-compatible or reconstructed objects for all four modes,
run twice, and report source/config/executable/output hashes, residuals,
restarts, timings, interface type, and geometry audits.  Missing author CAD
assets stay `RECONSTRUCTED_ASSET`; no author-code equivalence claim is allowed.
Commit as `exp: verify Spin-It Faster reproduction`.

## Task 12: Integrate OpenGL and close the milestone

Extend `ResultView`/viewer toggles for fixed shell, design domain, material,
cavity, signed quadric, CoM, target axis, and principal axes.  Add non-GUI
view-state tests, then run a real-desktop interactive smoke for every mode.
Run Faster ON/OFF builds, all CTest/Python suites, deterministic evidence,
source hashes, and clean-tree checks.  Commit as
`feat: visualize Spin-It Faster artifacts`, then set the Faster gate true only
if every requirement above passes.
