# Task 6d implementation report — finite discrete ellipsoid/Boolean moments

Bounded implementation on strict branch, starting at `ea478dd4b68e0a990c3ad8a0340076f1b68524ca`. Worktree `/data/xvxinyao/make-it-balance-paper-sim`; sole implementation/configure/build writer. Controller independently reviews this increment. No original-paper, whole-method, M4, release, merge, tag, or push claim follows. No original experiment or main-checkout geometry/output was edited. No subagents were used.

## Implemented contract and mathematical scope

Native `QuadricMomentEvaluatorSettings::derivativeMethod` selects `LegacyContinuousImplicitShape` (unchanged default) or `DiscreteEllipsoid`. Selected method, derivative applicability reason, source-feature topology, and a separately versioned scientific identity are returned. `solveQuadricHollowing` already accepts these evaluator settings; its returned native product now retains selected method and last returned/final-solution derivative diagnostic/scientific identity. Undefined derivative is rejected with its reason, including the unsupported-plane regression; it is never silently replaced by legacy quadrature or accepted as a zero derivative. The generic callback's historical empty-diagnostic error text is preserved.

`QuadricTriangulator` optionally retains each actual ellipsoid vertex's fixed icosphere normal. Metadata does not modify geometry or historical triangulation hash. Existing Task6b `QuadricInverseGaussMap::coefficientDifferential` supplies each moving source vertex's normalized coefficient differential. Its symmetric offdiagonal factor two and coefficient tangent map are used exactly once; no extra tangent projection is applied in the new path. Fixed design vertices have zero velocity.

The focused new helper builds constraints from actual CGAL birth facets and source vertex/edge/face membership. For an input facet, `N·dx = N·da - dN·(x-a)`. A source edge uses `(I-eeᵀ/|e|²)dx = (I-eeᵀ/|e|²)((1-s)da+s db)`; these constraints handle coplanar domain diagonals. Original source vertices supply their three velocity components. SVD must have three sufficiently independent directions; all ten right-hand sides must have a consistent finite solution. Rank-deficient or inconsistent systems, ambiguous features, missing construction provenance, nonfinite data, and input-vertex/other-operand facet contacts explicitly fail. The contact check scans opposite facets even when only one Boolean birth side survives. New conservative applicability certificates are scaled/barycentric membership `2e-9`, minimum singular ratio `1e-10`, and relative constraint residual `2e-9`; no old point cutoff, Boolean geometry setting, normalization, equation, or convergence tolerance was changed.

On the actual oriented output triangles, Reynolds transport integrates each of the ten raw monomials times interpolated vertex velocity dotted with the oriented area vector. Degree-two monomials times affine velocity have degree three. The four-point cubic-exact rule uses centroid weight `-27/48` and permutations of `(3/5,1/5,1/5)` with weight `25/48`. Density multiplies once. Production does not use finite differences or the 2014 moment derivative.

Source identities use fixed design/icosphere-derived vertex, edge, and facet indices, not coordinate-sorted Boolean vertex indices. Returned vertex provenance plus the sorted, orientation-preserving source-feature face connectivity hash expose combinatorial changes. Every evaluation rebuilds provenance; there is no stale derivative cache. The scientific v2 hash includes method, applicability, request identity, quadric coefficients, triangulation identity, actual material geometry, and provenance topology. Legacy geometry/topology/output hash serializers retain their formats. This is evaluator scientific metadata, not full build/run attestation; exact source/binary/cache hashes are recorded separately below.

Unsupported: ideal inverse-Gauss points; moving truncation/radius intersections and caps; plane, paraboloid, cylinder and hyperboloid moment differentials; contact/coincidence or rank/consistency failures; constructions without sufficient source feature provenance; arbitrary topology transitions. Python/UI bindings remain default-only and do not expose the new method selection/diagnostics in this bounded native increment. Cherchi primary adapter, arbitrary-quadric support, coefficient/raw-moment duality audit, geometry-refinement studies, and original experiment gates remain outside this increment.

## RED, implementation progression, and preserved failures

The TDD skill was read with its test-writing reference. The existing sphere fixture and old default `0.08` gate are unchanged. An opt-in `--strict-discrete-negative-control` asserts `1e-6` on that same legacy full-map test. Actual RED is retained in `task-6d-red-legacy.log`: exit 1, historical mismatch printed `0.006589` (legacy normal PASS prints about `0.00658869`), using the unchanged old `h=2e-5`. This failure preceded new production implementation. The dedicated desired-behavior native regression was written before the new helper/mode implementation. It independently calls the public Boolean/material integrator for all coefficient finite differences and uses a separate signed-tetrahedron moment oracle, not the production moment kernel.

The first attempt to build the new target before regenerating CMake failed with “No rule to make target”; retained `task-6d-first-build.log`. Configuration then succeeded. The first compiled native run passed all fifteen prescribed full-map FD checks and then failed because the test had set `maximumRestarts=0`, which existing request validation disallows. `task-6d-first-native.log` is preserved; only that test setting changed to valid `1`, not solver validation/constants. `task-6d-second-native.log` passed, including unsupported plane and actual sphere/cube tangency rejection. No clipped fixture or resolution was replaced after results.

Before the first native run the entrywise raw-FD test used relative `1e-4` when `max(|analytic|,|FD|)>1e-5`, otherwise absolute `1e-8`. The mandatory fine-step full-map gate was always `1e-6`. After seeing that clipped raw fine-step maximum entry-relative error was about `1.439e-6`, the controller required an additional strict `1e-6` entrywise Richardson check. This added check uses only the same prescribed three finite-difference matrices; it does not replace any raw step/error or lower any production/test threshold. The original `1e-4` raw entry gate is not claimed as strict `1e-6` entry verification. Both extrapolations must also agree within `1e-6` entrywise (near-zero absolute `1e-8`). This timing and control refinement are deliberate and disclosed.

## Fixtures and numerical results

All are analytic native unit fixtures, not original-paper inputs. Historical sphere: coefficients proportional to `[.25,0,0,0,-1,-1,-1,0,0,0]`, subdivision 3, cube `[-1,1]`, outer cube `[-2,2]`, density 1. Generic: exactly `[-1.4,.2,-.4,.1,1.1,.8,1.3,.12,-.08,.17]`, normalized by the unchanged public constructor; both signs; fixed subdivision 1; design cube `[-2,2]` for containment and `[-1,1]` for clipping; outer cube extends one unit; wall setting `.5`; density `1.7`. The generic target frame has origin `(.31,-.24,.18)` and rotation `.43` radians about normalized `(1,2,-1)`; design/outer geometry is rigidly mapped to model coordinates. No hidden resolution increase or fixture replacement.

Full-map norm errors are `||J-FD||/||FD||`; entrywise relative maxima apply above `1e-5`, otherwise the table gives maximum absolute error. Values from retained native/final focused logs:

| Fixture | h | Full-map relative norm | Max entry-relative | Near-zero max absolute |
|---|---:|---:|---:|---:|
| Historical sphere | 3e-4 | 2.0414328021494352e-7 | 5.5125021166096028e-7 | 3.7006971716534333e-13 |
| Historical sphere | 1e-4 | 2.2685084364861124e-8 | 6.1249214344230276e-8 | 1.1102577868573119e-12 |
| Historical sphere | 3e-5 | 2.0419271028008771e-9 | 5.5130208748451277e-9 | 3.6998368055605564e-12 |
| Contained positive | 3e-4 | 9.5662311310986945e-7 | 2.0304317338751598e-6 | 0 |
| Contained positive | 1e-4 | 1.0628935379966245e-7 | 2.2560535067247654e-7 | 0 |
| Contained positive | 3e-5 | 9.5412543313942613e-9 | 2.040820629028855e-8 | 0 |
| Clipped positive | 3e-4 | 6.0850095062030915e-7 | 1.4389583703861323e-4 | 0 |
| Clipped positive | 1e-4 | 6.7611680544982892e-8 | 1.5991213943113645e-5 | 0 |
| Clipped positive | 3e-5 | 6.0879212151168996e-9 | 1.4392637411640616e-6 | 0 |
| Contained negative | 3e-4 | 9.5662428616136071e-7 | 2.0304386597657439e-6 | 0 |
| Contained negative | 1e-4 | 1.0629157998727681e-7 | 2.2560718345301841e-7 | 0 |
| Contained negative | 3e-5 | 9.566801468057228e-9 | 2.0313473964036129e-8 | 0 |
| Clipped negative | 3e-4 | 6.085013791302652e-7 | 1.4389582087789423e-4 | 0 |
| Clipped negative | 1e-4 | 6.7610700880704114e-8 | 1.5990634970634609e-5 | 0 |
| Clipped negative | 3e-5 | 6.0893569585839762e-9 | 1.4387008325733803e-6 | 0 |

Observed successive norm-error ratios are approximately `.11111` and `.0900`, consistent with second-order central-difference truncation over these three steps. This is an observed three-step trend, not a global/asymptotic-order guarantee. Sphere near-zero entries show roundoff growth as h decreases. No step was discarded.

Supplemental Richardson uses `R_coarse=(9 FD_1e-4-FD_3e-4)/8` and `R_fine=(100 FD_3e-5-9 FD_1e-4)/91`:

| Fixture | Fine extrapolation max entry-relative | Near-zero max absolute | Coarse/fine extrapolations max entry-relative difference |
|---|---:|---:|---:|
| Historical sphere | 6.2933926151204869e-11 | 4.010943257515957e-12 | 6.218008762983511e-11 |
| Contained positive | 1.4857376341911909e-9 | 0 | 1.6343993897530926e-9 |
| Clipped positive | 2.5970977415560908e-9 | 0 | 3.9081409239568119e-9 |
| Contained negative | 5.724118703217334e-11 | 0 | 6.9415572263143906e-11 |
| Clipped negative | 3.4204279656825174e-9 | 0 | 2.5376336442401809e-9 |

Both clipped signs construct 60 actual intersection vertices, including 8 on coplanar domain-face diagonals. Worst constraint residual is `1.2071229572633428e-15`; minimum singular ratio is `.040658789036594041`. Every baseline/perturbed output vertex has a unique tested source-feature correspondence; independent test geometry checks verify actual vertex, segment, or triangle membership. The source-feature connectivity hash matches across all 60 perturbations per fixture; no canonical coordinate-order assumption is used.

Tests also check exact repeated geometry/Jacobian/scientific metadata, unchanged input hash and opt-in vs legacy geometry/triangulation hashes, radial null `<2e-10`, independent signed-tetrahedron ten-moment integration, density scaling, complement moments and complement differential (two negations imply `J(-q)=J(q)`), unsupported plane without fallback, selected native solver rejection with reason, and actual radius-one sphere/cube input-vertex tangency without geometry discard. The historical radius-.5 sphere's analytic-volume triangulation relative error remains `.0086061573702495009`, distinct from its `2.04e-9` fine-step differentiation error.

## Commands and final verification

All commands ran in the strict worktree. Each red/green/build/check has its uniquely named `task-6d-*.log` here; no log was deleted. `/usr/bin/time -v` wraps native/final verification commands. Build concurrency never exceeds `-j4`; primary and isolation builds run sequentially.

```sh
cmake --build build/strict-v2-spin --target quadric-moment-evaluator-regression -j4
/usr/bin/time -v build/strict-v2-spin/quadric-moment-evaluator-regression --strict-discrete-negative-control
cmake -S . -B build/strict-v2-spin
cmake --build build/strict-v2-spin --target quadric-discrete-moment-regression -j4
/usr/bin/time -v build/strict-v2-spin/quadric-discrete-moment-regression
cmake --build build/strict-v2-spin --target quadric-discrete-moment-regression quadric-triangulator-regression quadric-inverse-gauss-map-regression quadric-moment-evaluator-regression quadric-hollowing-solver-regression -j4
/usr/bin/time -v cmake --build build/strict-v2-spin -j4
/usr/bin/time -v ctest --test-dir build/strict-v2-spin -R 'spin-faster-(discrete-moments|quadric-triangulator|inverse-gauss-map|quadric-moments|newton)$' --output-on-failure -V
/usr/bin/time -v ctest --test-dir build/strict-v2-spin -L spin-it-faster --output-on-failure
/usr/bin/time -v ctest --test-dir build/strict-v2-spin --output-on-failure
```

The initial/final configure emits the existing libigl old-CMake compatibility deprecation; compilation emits the existing Boost global-bind-placeholder deprecation note. These are disclosed, not hidden or fixed by changing dependencies. Primary and isolation caches both keep PyBullet OFF; the latter keeps Faster ON / 2014 OFF.

The final feature-isolation commands were:

```sh
cmake -S . -B build/task-6b-faster-only -DMAKE_IT_STAND_BUILD_SPIN_IT_FASTER=ON -DMAKE_IT_STAND_BUILD_SPIN_IT_2014=OFF
/usr/bin/time -v cmake --build build/task-6b-faster-only --target quadric-discrete-moment-regression -j4
/usr/bin/time -v ctest --test-dir build/task-6b-faster-only -R '^spin-faster-discrete-moments$' --output-on-failure -V
git diff --check
```

Final source verification completed with no failures:

| Check | Result | Wall time | Maximum RSS (KiB) | Evidence |
|---|---|---:|---:|---|
| Primary all-target build | exit 0 | 74.70 s | 3,345,056 | `task-6d-final-all-build.log` |
| Focused new/triangulator/inverse-map/legacy-moment/Newton | 5/5 | 5.73 s | 18,776 | `task-6d-final-focused.log` |
| Faster suite | 15/15 | 25.97 s | 75,868 | `task-6d-final-faster.log` |
| Full primary CTest | 183/183 | 194.95 s | 4,230,392 | `task-6d-final-full.log` |
| Faster ON / 2014 OFF build | exit 0 | 58.66 s | 3,343,700 | `task-6d-isolation-build.log` |
| Faster-only new native regression | 1/1 | 5.16 s | 17,788 | `task-6d-isolation-native.log` |

Source/test/report whitespace checks are clean. After staging raw evidence logs, unrestricted `git diff --cached --check` also reports eight preexisting-style trailing-space output lines in the two CMake configure logs (`BOOST_ROOT`/Boost discovery); those raw evidence bytes are preserved. The final focused log separately preserves the old regression's ordinary PASS at `derivative_error=0.00658869`; this is not relabeled as discrete verification. The new path's historical sphere geometry hash remains `0f00208452856b310c3735782290374aeac11b185b640fe0a157d102ef387101` in both modes.

## Exact final source/build identity and ownership

`task-6d-primary-sha256.log` records SHA-256 for all changed native source/header/test files and CMake, plus primary executable/core/cache. `task-6d-isolation-sha256.log` records the isolation executable/core/cache and reused unchanged Task6b source/header. These were measured after the final-source builds/tests, without source changes afterward.

| Artifact | SHA-256 |
|---|---|
| New helper source | `5e727428510818574a71257580155db46d8f2a70aa5ed55f1cc4bdd771c97e30` |
| New native regression source | `b7ae612c055ca4b37b9053046d5a7518804d75e847289ec37bfed0184979ce07` |
| Primary new native executable | `13a06629bcc51ed5a4b4a96c1ff4c7739641a4e8cfbc9f7bdb1d821ad02c2878` |
| Primary core archive | `0778661d61a7c14bf575052812523375f348c5af2ab8fc0c890a145a4b139553` |
| Primary CMake cache | `2c3f5a72c82d0ba5d8b5d0cd8869d13d89ea8a21e0daa405e88242e1778e517c` |
| Faster-only new native executable | `e6c283eff8bc29864e3ef9274bafb8a33e5a8869ca41d6d00532abd18d5345f8` |
| Faster-only core archive | `9570a4277763b58b2a2796241f9dfebf9d2bac45133145bffee79f641ae0f1b5` |
| Faster-only CMake cache | `bf9ce45dc0aee97188c0bd023182f18bf0b1d7fc927ee08ad54d38be72718127` |

All five fixtures' actual geometry, source-feature topology, and scientific hashes are printed in full in the final focused and isolation native logs. The scientific metadata is deterministic across repeated evaluations, but is explicitly not a substitute for these independently measured build artifacts.

The scoped implementation/evidence commit uses configured `volcano <xuxinyao307@gmail.com>` identity; the exact containing commit hash is returned to the controller separately (avoids a self-referential report hash). Source/configure/build ownership is handed back after this final verification and scoped commit. Controller independent task review remains the next gate; no merge/tag/push or release action was performed.

## Self-review

Checked transport orientation/density against independent signed-tetrahedron integrals and sign-complement identities; checked normalization/factor-two are confined to reused Task6b; checked edge constraints and the actual 8 coplanar-diagonal intersections; checked numerical rank and residual gates precede applicability; checked explicit vertex-contact rejection independent of rank; checked derivative selection and failure reason reach the native solver; checked legacy hash format and geometry remain unchanged; checked helper/CMake require no 2014 implementation; checked all prescribed FD steps, provenance and supplemental accuracy results remain visible. Controller task review is still independent and required.
