# Spin-It strict computational reproduction implementation ledger and plan

The binding specification is the user's approved “Spin-It 与 Spin-It Faster 严格完整计算复现方案” of 2026-09-05. This execution breakdown does not supersede or reduce its M0–M8 scope.

## Global constraints

- Work only in the existing linked worktree, branch `spin-it-strict-computational-v2`, based on `e1190569294a8853bc12a86c6fd9d31b2428de05`.
- Preserve all historical scientific outputs; old synthetic runs do not prove paper reproduction. Do not modify the dirty main checkout.
- No simplified or reconstructed asset may satisfy an original-paper experiment. Unknown critical parameters block that experiment.
- Complete computational reproduction is distinct from author bitwise equivalence; physical validation is `NOT_IN_SCOPE`.
- No manual booleans grant verification. No missing requirement may disappear. No failed solve is an infeasibility certificate.
- Native/Python compatibility remains; legacy behavior must be labelled and cannot be promoted into strict results.
- Small analytic fixtures are tests only. Never modify a test model or tolerance to manufacture a paper success.
- Compile with at most 4 jobs; heavy experiments sequentially, 24 GiB / 60 minutes per case batch. Full resolution batches require 100 GiB free, stop below 20 GiB. Never delete old files to make space.
- Use test-first development and fresh verification. Commit only explicitly scoped changes, using the already configured user identity.

## Task 1: Historical claims correction and fail-closed strict requirements audit (M0 / M1 foundation)

Implement in `experiments/spin_it_reproduction/strict/` a versioned, fixed requirements catalogue and audit CLI. Keep the catalogue separate from its audit logic. Cover individually:

- 2014 Eq.1–8, static objective, continuous fill and binary semantics, boundary subtrees, split/merge, Laplacian, SLQP/active set/L-BFGS, complete octree cage and BBW interpolation, joint deformation, dual density, supplement Algorithms 1 and 2.
- 2014 original cases: ellipsoid, heart, armadillo hand/shell, teapot three distinct objectives, cuboid, woven ring, dancing couple hollowing and joint deformation, Sumone dual density, T-Rex matched Make-It-Stand comparison, elephant, and each Table 1 cavity/cut family.
- Faster raw moments/coefficient duality, theorems and degenerate conditions, Eq.14 Newton/geodesic/multiplier/restart, differentiable discrete moment map, inverse Gauss map/ideal points/caps/degenerate types, independent Boolean audit, separate 2D dogleg/polygon-circle/extrusion, actual hydrostatics and chained floating-island assets, CAD roundtrip.
- All eight Faster Table 1 cases, Figures 10,11,12,13, both section 6.6 examples. Shared sealed geometry/independent moments/provenance/determinism/registered tests/independent simulations/real desktop visualization.
- Physical measurements retained explicitly as out of scope.

Each catalogue row contains requirement_id, paper_locator, requirement_kind, expected_outcome, dependencies, and a required evidence description. Current evidence must state NOT_VERIFIED, NOT_IMPLEMENTED or BLOCKED honestly; do not infer verification from existing tests. Differentiate asset class and parameter evidence. Exact required IDs come from the internal catalogue, not caller-selected rows. Reject missing, duplicate and unknown IDs, nonfinite data and malformed evidence.

Provide asset preflight that allows only AUTHOR_ORIGINAL or PUBLISHED_EXACT_SPECIFICATION with checked file SHA-256 and source reference and all declared required critical parameters with PUBLISHED/AUTHOR_PROVIDED/EXACTLY_DERIVED evidence. No UNKNOWN_CRITICAL, reconstructed or source-only inputs. Readiness is not verification. Until independent requirement checkers are implemented, the audit must remain PARTIAL_BLOCKED even if input says passed=true or status=VERIFIED. This phase must not implement a generic trust-the-JSON positive publication path.

Retire the old `build_release_gate.py` grant mechanism: keep import compatibility if needed but legacy aggregate can never grant public_complete_reproduction, formula_coverage or experiment_coverage from caller booleans. Remove --tests-passed CLI, explain supersession and legacy scope, preserve old result JSON bytes. Update its tests to exercise rejection and compatibility, not source-text matching. Add new strict audit/preflight behavioral tests, register a lightweight CTest within existing Python paper-test configuration, without requiring unavailable native extensions.

Tests first: existing synthetic all-true payload cannot publish; incomplete/duplicate catalogue, invented pass flags, unknown parameter, bad digest, source-only/reconstructed asset, NaN and empty evidence rejected; genuine exact-specification fixture can be READY but cannot publish; physical entries cannot be dropped; CLI exits nonzero with detailed blockers. Existing raw scientific result files must be unchanged.

The controller handles human correction/provenance documents separately. Do not edit those documents or solver code. Verify focused Python tests, commit scoped changes, report red/green evidence and limitations. Never call this task full M1 or full evidence-chain completion.

## Task 2: Independent signed-tetrahedron mass-moment checker (M2 subset)

Add a standard-library Python independent checker under `python/make_it_stand/verification/`, not importing the native solver integrator. Explicit geometric versus density-weighted moments; finite and oriented closed manifold validation, signed cavity components; mass, first and symmetric second moments and central inertia from signed tetrahedron integration. Validate closure/winding and reject invalid/zero-mass input, do not claim full self-intersection or wall-thickness audit. State limits explicitly.

Tests use hand-derived translated tetrahedron, rotated/non-diagonal inertia, box, hollow box with reversed inner component, multi-density aggregation, translation/scale/density covariance, invalid faces/nonfinite vertices/open boundary/zero mass. Error target 1e-10 normalized for analytic nondegenerate polyhedra. No solver reuse or test-derived expected values. Keep this checker out of strict release success until complete sealed geometry checks exist.

## Task 3: Distinct 2014 objectives (M3 subset)

Introduce Spin2014Objective with YOYO, TOP_INERTIA_ONLY, TOP_LOWERING_ONLY, TOP_MASS_AND_LOWERING, STATIC_LOWERING. Expose independent inertia and center weights. Implement fyo=gammaI*((Ia/Ic)^2+(Ib/Ic)^2), inertia-only=fyo, lowering-only=fyo+gammaC*l^2, full=fyo+gammaC*(l*M)^2, static=gammaC*l^2 with analytical gradients/Hessians. Preserve existing factory API defaults and explicitly identify old/angle-eliminated semantics. Thread objective selection through adaptive request hash, functional, Python facade and paper runner; teapot strategies must call different objectives. Tests compare independent formula values, multiple weights, finite-difference derivatives and density scaling (height penalty density-independent, full penalty density squared). Do not claim Givens-phi, octree or cage completion from this task.

## Task 4: Remaining common verified artifact chain (M2)

Implement typed geometric moments; frozen source/build/runner/checker provenance; SealedSpinArtifact for raw/final outer, shell/design/material/cavities/densities/targets/hydrostatics and independent reintegration. Add watertight, winding, nonmanifold vertices/edges, degeneracy, self intersection, wall thickness, containment, cavity components and mesh/CAD audits. No silent repair or stale mass properties. Hollow-only outer must remain byte-identical. Feasibility must distinguish verified witness, applicable infeasibility certificate and UNKNOWN. Discrete certificates cannot imply continuous infeasibility.

## Task 5: Remaining complete 2014 method (M3)

Restore explicit Givens phi and equality constraints; source-surface native sparse octree with boundary subtrees, initial midlevel, genuine level9 refinement, continuous solve, split/merge, Laplacian and binary validation. Implement full SLQP/active-set/feasibility-restoration/L-BFGS/merit/KKT and history changes. Extract complete cage from intersecting cells (no handle subsampling), finer-grid BBW, trilinear interpolation, joint beta/h optimization and published regularizers with precommit geometry checks. Full dual-density material products, measured cavities and original manual-cut digital replay. All original cases need original asset/pose/thickness/weights; missing data blocks, expected hollow-only insufficiency remains a valid comparison outcome. T-Rex baseline must use identical inputs and discretization.

## Task 6: Complete Faster method (M4)

Audit coefficient/moment metric, Eq14 on sphere, per-step multipliers, 0.5 geodesic cap/backtrack, 1e-12 residual, 60 function evaluations/attempt, mode-aware positive-mu/nonzero-Jacobian initials and recorded restarts. Triangulate via fixed sphere normals, inverse Gauss, ideal points/caps outside domain and explicit degenerate types. Primary Cherchi public commit 7bd6c2697695f555f9a15d7e154e084d841316f0, existing pinned CGAL independent secondary; no fallback. Differentiate actual discrete moment map, including intersection vertices under fixed topology, verify 1e-6 derivatives and refined geometry error independently. Actual 2D six moments, polygon/disk and dogleg; extrusion equivalence. Actual geometric hydrostatics, floodable parts, chained A-then-B island. OCCT 7.9.3 commit a016080bf6738d6aeae020badee4e888ad1540a5 for original CAD quadric/B-rep/STEP roundtrip and independent masses. No original CAD means CAD experiment blocked.

## Task 7: Original experiments, independent simulation, desktop viewer (M5/M6)

Eight Faster Table1 cases, each 50 top-level random runs, entire sequence repeated twice. Retain all filtered starts/restarts/failures, paper definition first success=1 attempt, times/geometry/residuals and all outputs. Figures10–13: actual quadratic convergence estimates above noise; original eggs/heart Pareto at 41 initial weights plus refinement and exact endpoint objectives; retain dominated rows; true geometry subtype changes; certified near-feasibility bracket and nine offsets using common original initialization with sufficiently small boundary uncertainty. Both 1D p=0 cases with correct outcome. Original 2014 matrix from Task1 including computed angular-velocity comparisons, matching T-Rex and true cavity/cut counts. Missing critical data cannot be replaced by chosen values.

Separate ideal torque-free and fixed-pivot dynamics, support-convex-hull tilt, re-evaluated hydrostatic perturbations, eight directions, energy/momentum/time-step convergence. Collision-oracle mismatch remains blocker; no friction tuning to pass. Real desktop OpenGL loads only sealed products, original pose, shell/material/cavity layers, cutaway, axes/CoM/support/waterline, trace replay and persistent evidence status. CAD, display and physics must bind the same artifact.

## Task 8: Strict complete verification and integration (M7/M8)

All fixed requirements verified with independent checkers and no missing original assets or critical parameters; complete original experiment outcomes, two same-environment scientific hashes, all registered required tests actually run without unexpected skip/Not Run. Mutation tests must reject missing/duplicate/unknown cases, proxy source/upsampling, identical teapot objectives, known-solution promotion, copied cavities, manual true, NaN/empty statistics, dropped failures, geometry/mass mismatch, stale build/input, wrong expected-outcome logic and single-run determinism. Report paper/local/difference/precision/distributions; distinguish computations from physical references and implementation performance from original timings. Release only from clean audited scoped commits, preserve main checkout, no publish/tag/merge with unresolved blockers. Final state PARTIAL_BLOCKED otherwise; no complete claim merely for code completion.
