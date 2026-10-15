# Spin-It Fixed-Shape Quadric Reimplementation Plan

> **For agentic workers:** Use superpowers:executing-plans inline, with test-driven-development. Preserve the dirty worktree and all historical evidence. No automatic commits, merges, or deletion of scratch evidence.

**Goal:** Implement the five Spin-It 2014 objectives and their physical constraints using Faster quadric materials, with independently checked top, yo-yo and static results.

**Architecture:** An owned callable mass-moment problem feeds the existing spherical Newton loop. An explicit `SPIN2014_OBJECTIVES_QUADRIC_V1` adapter reuses `SpinIt2014MomentProblem` and applies a declared-unit chain rule. The static plane path obtains candidates but freshly verifies squared-height optimality. Existing material and evidence systems remain authoritative.

**Tech Stack:** C++17, Eigen, CGAL (default), optional explicit Cherchi, pybind11, Python, CMake/Ninja and existing source-job resource limits.

**Spec:** The user's latest full text, “新目标实施规划：用二次曲面求解实现 Spin-It 的固定外形功能”, supplied and approved in this chat. It supersedes the LEGACY predictive transaction task, not the historical files.

## Global Constraints

- Work only in `/data/xvxinyao/moment-guided-spin-v1/code/make-it-balance-paper-sim`; preserve existing changes. The linked worktree on `spin-it-strict-computational-v2` is reused.
- Five objectives: YOYO and TOP_INERTIA_ONLY = gamma_I R; TOP_LOWERING_ONLY = gamma_I R + gamma_c h^2; TOP_MASS_AND_LOWERING = gamma_I R + gamma_c (M h)^2; STATIC_LOWERING = gamma_c h^2.
- Coordinates are about the target origin/axis. Mass moments already contain density; the objective adapter must never multiply density again. Default declared mass/length scale factors are one. A normalized caller must explicitly supply its inverse unit map.
- Public input fixed material must be measured from geometry, never trusted hand-filled moments. Density applies separately once to shell and design material. Outer vertices and faces are unchanged.
- Tops/yo-yo use full-rank discrete quadrics and the strict raw/dual Newton path; CGAL default. No silent backend fallback. Static uses its plane-family solver with fresh squared-height verification for positive mass and h >= 0; h=0 is a distinct degenerate case.
- Per-run normal subdivision is fixed at level 2. Topology/cap/parameter branch changes require fresh derivative evaluation. Unsupported ranks reject/restart; low-dimensional Jacobians are not padded.
- Analytic moments <=1e-10; derivatives <=1e-6 at three steps; final independently integrated moments/constraints <=1e-6. Newton raw residual <=1e-12, alpha=t_max=0.5, <=60 evaluations per attempt, <=8 attempts.
- Compile <=4 workers; sequential large jobs, numerical libraries single-thread; job 24 GiB, swap zero, 60 min. Start/ongoing free-space gates remain those of the registered batch, ongoing >=20 GiB.
- No voxel/SLQP/octree/G8/G9/cage development, new videos, manufacturing or complete-paper experiments. Keep all old artifacts. Printing thickness remains unverified, not silently passed.
- Success name is fixed-shape quadric reimplementation, not perfect reproduction of either paper. All five strategies and all three frozen functions need real updated material, independent checks and process repeats.

## Review Focus

1. Mass-weighted height is L_z^2, not h or h^2; test unequal masses and declared units.
2. Unit conversion must transform gradient, Hessian and constraints, not only the objective value.
3. Historical Newton wrappers must still use the same kernel and retain their scientific traces.
4. A plane height candidate cannot borrow its old success flag for the new squared objective; check h=0 separately.
5. A rejected floating mesh must not become a qualified artifact just because exact material succeeded.

## Frozen first batch (freeze before any optimization)

| Input | Source/configuration | Objectives |
|---|---|---|
| Heart | `experiments/reconstructed_heart_top/heart_top.py::heart_mesh(48,12)`, existing design scale 0.72 and original pivot/axis | Three top strategies |
| Cuboid | `experiments/reconstructed_cuboid_yoyo/run.py`, Y rotation 5 degrees, `offaxis-pairs`, original target axis | YOYO |
| Static | Existing faceted-shell outer/domain; target `(0.05,0,-1.2)`, real support hull | STATIC_LOWERING |

Weights gamma_I=gamma_c=1 in declared units; inactive terms absent. Seeds 0 and 1 for 3-D, request-derived deterministic static initial guess. No old endpoint, known answer or checkpoint may initialize optimization. Failed seeds stay in the ledger; no replacement or input tuning. Every successful configuration repeats in another process, excluding timing from scientific hashes.

### Task 1: Shared moment problem and declared-unit 2014 adapter

**Files:** `include/spin/MassMomentProblem.h`, `include/spin2014/QuadricMomentProblem.h`, `src/spin2014/QuadricMomentProblem.cpp`, existing solver header/source, CMake and `tests/spin2014_quadric_moment_adapter_regression.cpp`.

**Interfaces:** `MassMomentProblem` owns `std::function<SpinMomentEvaluation(const Vector10&)>` and can copy an existing `SpinMomentProblem`. `makeSpin2014QuadricMomentProblem(settings)` returns it; `spin2014QuadricObjectiveHash(settings)` binds objective, weights, unit labels and mass/length scale. `solveMassMomentSystem(...)` calls the SAME Newton implementation as the historical entry.

- [ ] Write tests for five actual formulas, distinct constraints, non-diagonal inertia, zero/negative/nonfinite scales, transformed gradient/Hessian and three-step differences. Build RED: missing new interface, then implement reuse of the existing 2014 evaluator, not new formulas.
- [ ] Unit map D is diagonal: mass factor m, first factors m*l, second factors m*l*l. Evaluate at D*s; return D*g, D*H*D and corresponding constraint derivative transforms. Density is already in s.
- [ ] Compare legacy and callable wrapper trace/event/status/evaluation counts on the same existing native regression map; capture actual 2014 objective values through the new Newton entry. No synthetic success is counted as a geometry case.
- [ ] Build and run affected 2014/Newton/subspace tests under the source-job envelope; save exact command and result. Do not call the overall goal complete.

### Task 2: Request-backed top and yo-yo adapter

**Files:** existing solver and request types, bindings, `python/make_it_stand/spin.py`, existing artifact problem/evidence adapters; a focused 2014 quadric adapter module if needed.

**Interfaces:** explicit objective settings plus a `SpinDesignRequest`; original Faster call signature remains unchanged. Bind the objective settings hash to evaluator request identity and all result attachments.

- [ ] Write mode-mismatch, hand-filled shell rejection, unit metadata and backend identity tests before adding request-backed integration.
- [ ] Reuse the existing full-rank CGAL evaluator, material export, complement and event ledger. Fixed shell measured from actual geometry and transformed once. Explicit discrete derivative + strict dual/raw settings for this profile.
- [ ] Native/Python adapter returns original outer, fixed shell, retained design material, cavity, coefficients, moments, all trials/restarts and classified failures.
- [ ] Run local real-Boolean derivative/moment checks; unsupported strata remain failures, never disguised as infeasible.

### Task 3: Static squared-height candidate and fresh verification

**Files:** existing plane solver and 2014 adapter; existing support checker; focused tests.

- [ ] Tests reject mass <=0 and material below support; distinguish h=0, h>0 and gamma_c=0.
- [ ] For h>0 obtain the existing height-minimizing plane candidate, freshly recompute gamma_c*h^2 and its analytic derivatives, equality residuals and scaled KKT multipliers/positive optimality condition; no inherited pass boolean.
- [ ] Check actual support hull and target projection independently. Keep static solver tolerances separate from Newton's raw Eq.14 gate.

### Task 4: Freeze inputs, run five strategies and reload independent materials

**Files:** a new versioned function batch under `experiments/spin_it_reproduction`, existing `SpinArtifactProblem`, `SealedSpinArtifact` and independent verifier.

- [ ] Persist all three geometries/configs and hashes before solving. Register seed 0/1, eight attempts and level 2. Seal density, units and target axis/origin.
- [ ] Run five strategies from independent initialization; require accepted objective updates and actual material change, not only initial feasibility. Preserve each failed seed.
- [ ] Reload saved geometry and independently integrate full moments, density, inertia and constraints. Audit exact partition, containment, overlap, winding, degeneracy and unchanged outer surface. Float-derived outputs get separate pass/failure witness; thickness NOT_VERIFIED.
- [ ] Repeat each successful configuration in a distinct process and compare events/material/scientific hashes. A failed case remains open with its actual failure stage; do not tune input/weights/tolerances.

### Task 5: Whole-goal audit and delivery

- [ ] Run all required affected/core algorithm suites at the milestone, report any historical failures by name, perform final fresh-context review, fix confirmed important issues with RED→GREEN tests.
- [ ] Produce the five-objective/three-function acceptance table, implementation/support limits, evidence hashes, all failures and reproduction commands.
- [ ] Do not declare complete until all five actual strategies, independent materials and repeated processes meet the user's gates. Do not claim manufacture, simulation, arbitrary-model success or perfect paper reproduction.

## Approved 2026-10-08 precision completion sequence

The user's latest implementation request approves the following remaining work on this SAME fixed-shape goal. Preserve all earlier evidence and legacy defaults.

1. Explicit rational source profile: finite and truncation coordinates with absolute enclosures <=2^-384; the same rational source, Boolean, material, moments and derivative feed Newton and saved outputs. Start with ellipsoids, then the existing hyperboloid normal-edge clipping/caps; preserve level 2 and radius budgets. Rational points are approximations of analytic roots, not exact analytic quadric points. Rational identity/ordering must not weld distinct vertices whose display coordinates coincide.
2. Preserve all original three-step central-difference failures. A separately registered five-point protocol keeps each original h and uses outer scales h,h/2,h/4 with nodes {-H,-H/2,0,H/2,H}. Actual stored raw abscissae define rational Lagrange derivative weights; apply the same weights to actual normalized q and independently integrated moments and compare to the production base Jacobian. All ten columns and three scales retain 1e-6. The full production batch needs 81 evaluations; independent 61-point precision probes are not production evidence. No post-hoc shrinking or subtraction of estimated truncation errors to grant a pass.
3. Only after those gates, freeze and re-run affected Heart/Cuboid configurations with original seeds 0/1 and independent successful repeats; retain all failures. Re-audit unchanged static saved evidence rather than rewriting its plane path. Final core tests, independent review and five-strategy material/trace audit remain necessary.
4. Authority and floating derivatives retain separate dispositions; rejected float outputs are prohibited from qualified downstream use, but repairing every float OBJ is not an added algorithm-goal prerequisite. No manufacturing/simulation/voxel/rank-two generalization is added.
