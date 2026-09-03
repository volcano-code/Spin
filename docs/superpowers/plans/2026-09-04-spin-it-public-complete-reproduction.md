# Spin-It Public-Complete Reproduction Implementation Plan

> **Execution rule:** apply test-driven development to each behavior change and
> verify each commit before moving to the next task.

**Goal:** Close every publicly verifiable computational experiment in Spin-It
(2014) and Spin-It Faster (2024), while keeping author-exact and physical claims
explicitly false when source assets or fabrication evidence are unavailable.

**Architecture:** Extend the existing shared ten-moment kernel and immutable
artifact APIs.  Add the missing Faster buoyancy problem without coupling it to
the 2014 solver.  Build paper experiment runners on top of public APIs, then
strengthen the 2014 adaptive solver using resolution-ladder regression tests.
All generated evidence flows through a new aggregate gate with multidimensional
claim fields.

**Technology:** C++17, Eigen, libigl/CGAL mesh Booleans, pybind11, Python 3,
NumPy, Matplotlib, CTest, existing OpenGL viewer.

---

## Task 1: Freeze source and experiment inventory

**Files:**

- Modify: `experiments/spin_it_reproduction/sources/manifest.json`
- Create: `experiments/spin_it_reproduction/paper_inventory.json`
- Create: `tests/python/test_spin_it_paper_inventory.py`
- Modify: `CMakeLists.txt`

1. Write a failing schema test requiring both papers, supplement/video URLs,
   checksums or publisher identifiers, all paper examples, and independent
   claim fields.
2. Run the inventory test and confirm failure.
3. Add the versioned inventory and author-asset lookup policy.
4. Register the test and verify it passes.
5. Commit as `test: inventory complete Spin-It reproduction scope`.

## Task 2: Add Faster buoyant-stability request contract

**Files:**

- Modify: `include/spin/SpinDesignRequest.h`
- Modify: `src/spin/SpinDesignRequest.cpp`
- Modify: `tests/spin_design_request_regression.cpp`
- Modify: `python/make_it_stand/spin.py`
- Modify: `tests/python/test_spin_it_faster.py`
- Modify: native pybind source located from the existing `SpinMode` binding

1. Add failing C++ and Python tests for `BUOYANT_STABILITY`, liquid density,
   displaced volume/center, level-area moments and shell volume validation.
2. Verify failures are due to absent API.
3. Add the enum, token, immutable fields, request hashing and mode-specific
   validation.
4. Bind and expose the new fields in Python.
5. Run focused C++ and Python tests.
6. Commit as `feat: add Faster buoyancy request contract`.

## Task 3: Implement Faster buoyant moment problem

**Files:**

- Modify: `include/spin/SpinMomentProblem.h`
- Modify: `src/spin/SpinMomentProblem.cpp`
- Modify: `src/spin/QuadricHollowingSolver.cpp`
- Modify: `tests/spin_moment_problem_regression.cpp`
- Modify: `tests/quadric_hollowing_solver_regression.cpp`

1. Add analytic failing tests for equilibrium constraints: horizontal center
   of mass and prescribed total material volume.
2. Add independent metacenter feasibility tests.
3. Implement the published minimum-center objective and three constraints.
4. Route buoyancy through Eq. 14 and add application-aware initial quadrics.
5. Check objective gradient, Hessian and constraint Jacobian by centered finite
   differences.
6. Commit as `feat: implement Faster buoyant stability`.

## Task 4: Reproduce all Faster interface families

**Files:**

- Modify: `tests/spin_it_faster_end_to_end.cpp`
- Create: `tests/spin_it_faster_interface_types.cpp`
- Modify: `CMakeLists.txt`

1. Add failing oracles for both static objectives and buoyancy planes.
2. Add fidget circular-cylinder, yo-yo ellipsoid/two-sheet hyperboloid, and top
   one-sheet/two-sheet hyperboloid fixtures.
3. Fix classification, initialization or Boolean differentiation defects found
   by these tests without adding case-specific result constants.
4. Require Eq. 14 `<=1e-12`, physical residual `<=1e-6`, sealed geometry and
   identical repeat hashes.
5. Commit as `test: cover Faster paper interface families`.

## Task 5: Build eight Faster reconstructed assets

**Files:**

- Create: `experiments/spin_it_reproduction/faster/paper_cases.json`
- Create: `experiments/spin_it_reproduction/faster/generate_paper_assets.py`
- Create: generated asset manifests under
  `experiments/spin_it_reproduction/faster/paper-assets/`
- Create: `tests/python/test_spin_it_faster_paper_assets.py`

1. Test deterministic asset hashes, units, closed meshes, shell/domain nesting
   and required attribution fields.
2. Construct parameterized counterparts for the eight Table 1 examples.
3. Independently integrate fixed-shell and domain moments.
4. Record target frames, waterlines, densities and objective weights.
5. Commit as `test: add reconstructed Faster paper assets`.

## Task 6: Implement Faster 50-start Table 1 protocol

**Files:**

- Create: `experiments/spin_it_reproduction/faster/run_paper_table1.py`
- Create: `experiments/spin_it_reproduction/faster/table1_schema.json`
- Create: `tests/python/test_spin_it_faster_table1_runner.py`

1. Test a two-start miniature run for restart accounting and timing separation.
2. Implement 50 stored random initial points per case and the paper's 60
   function-evaluation cutoff.
3. Record success, local maximum, non-convergence, iterations, evaluations,
   residual and derivative timings.
4. Exclude wall time from deterministic hashes.
5. Commit as `test: reproduce Faster Table 1 protocol`.

## Task 7: Reproduce Faster Figures 10--13

**Files:**

- Create: `experiments/spin_it_reproduction/faster/run_paper_figures.py`
- Create: `experiments/spin_it_reproduction/faster/figures/README.md`
- Create: `tests/python/test_spin_it_faster_paper_figures.py`

1. Test plot inputs and monotone/finite curve invariants.
2. Generate application-grouped quadratic convergence data.
3. Sweep both static and top Pareto weights with warm starts and record quadric
   type transitions.
4. Binary-search the tilted-cube feasibility limit and run distances
   `d_crit - 10^-m`, `m=1..9`.
5. Reproduce both degenerate one-dimensional examples analytically.
6. Commit as `test: reproduce Faster paper figures`.

## Task 8: Strengthen Spin-It 2014 adaptive convergence

**Files:**

- Modify: `include/spin2014/SlqpLbfgsSolver.h`
- Modify: `src/spin2014/SlqpLbfgsSolver.cpp`
- Modify: `include/spin2014/AdaptiveVoxelSpinSolver.h`
- Modify: `src/spin2014/AdaptiveVoxelSpinSolver.cpp`
- Create: `tests/spin_it_2014_resolution_ladder.cpp`
- Modify: `CMakeLists.txt`

1. Add failing 4/8/16-resolution top and yo-yo convergence tests.
2. Trace KKT residual, active-set transitions, line search, refinement and
   binary projection to identify the first divergence.
3. Correct scaling, projected step, merit-function or warm-start defects rather
   than relaxing tolerances.
4. Add bounded level-9 smoke coverage and resource-limit reporting.
5. Commit as `fix: stabilize Spin-It 2014 adaptive convergence`.

## Task 9: Build Spin-It 2014 reconstructed paper cases

**Files:**

- Create: `experiments/spin_it_reproduction/spin2014/paper_cases.json`
- Create: `experiments/spin_it_reproduction/spin2014/generate_paper_assets.py`
- Create: `tests/python/test_spin_it_2014_paper_assets.py`

1. Define the complete public example inventory and expected experiment role.
2. Generate or import attributable deterministic meshes.
3. Store target axes, densities, cages, weights, cut accessibility metadata and
   asset evidence level.
4. Verify manifoldness, units, repeat hashes and outer/domain nesting.
5. Commit as `test: add reconstructed Spin-It paper assets`.

## Task 10: Reproduce Spin-It 2014 paper comparisons

**Files:**

- Create: `experiments/spin_it_reproduction/spin2014/run_paper_cases.py`
- Create: `tests/python/test_spin_it_2014_paper_runner.py`
- Modify solver files only when a tested general defect is exposed

1. Run ellipsoid, heart, two Armadillos, three teapot variants, two yo-yos,
   dancing couple variants, dual-density Sumone, T-Rex and elephant families.
2. Compute alignment, objective, mass, center height, cavity components and
   accessibility-cut metrics.
3. Add quasi-static incline simulation for reconstructed T-Rex.
4. Require two deterministic repetitions and sealed geometry.
5. Commit as `test: reproduce Spin-It paper comparisons`.

## Task 11: Aggregate evidence and honest release gates

**Files:**

- Create: `experiments/spin_it_reproduction/build_release_gate.py`
- Create: `experiments/spin_it_reproduction/REPRODUCTION_REPORT.md`
- Modify: both solver-specific `REPORT.md` and `release-gate.json`
- Create: `tests/python/test_spin_it_complete_release_gate.py`

1. Add a failing test that rejects a gate lacking formula, experiment, asset,
   author-exact or physical flags.
2. Aggregate machine-readable results without treating blocked author assets as
   computational failures.
3. Require all public computational rows to pass before setting
   `public_complete_reproduction=true`.
4. Keep `author_exact=false` and `physical_validation=false` prominent.
5. Commit as `test: add complete Spin-It reproduction gate`.

## Task 12: OpenGL, full verification and freeze

**Files:**

- Modify viewer code only for missing buoyancy/waterline controls
- Create: `experiments/spin_it_reproduction/viewer/README.md`
- Create: real-desktop screenshots and machine-readable smoke reports

1. Add a viewer regression for buoyancy waterline and displaced-volume overlay.
2. Run representative Faster and 2014 results on host `DISPLAY=:1`.
3. Run all registered tests; require zero `Not Run` and unexpected skips.
4. Regenerate reports from the clean source commit and verify SHA-256 manifests.
5. Request code review, address findings, and rerun verification.
6. Commit verified evidence and freeze a reproducible release commit; create a
   tag only after every public-complete gate passes.

