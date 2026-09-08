# Task6e: opt-in dual-monomial Newton and raw-equation stopping

Workspace `/data/xvxinyao/make-it-balance-paper-sim`, BASE `f433a149ad71e0283f3643596d6eb6711eb1efff`, branch `spin-it-strict-computational-v2`. This bounded implementation does not complete M4, any original experiment, or the strict M0–M8 plan. The dirty main checkout was not accessed or modified. No subagents, merge, tag, push, or historical scientific-output edits.

## Contract and implementation

`SpinSolverSettings::equationProfile` defaults to `LegacyNativeNormalizedV1`; explicit `DualMonomialRawV1` selects the independent literal-monomial reconstruction. Python selects `SpinSolverConfig(equation_profile="DUAL_MONOMIAL_RAW_V1")`. Invalid enum/profile values reject; strict budgets above 60 evaluations per attempt reject. Existing settings and result aggregate members keep their original order, with all new fields appended. `SpinSolverSettings{7}` retains its original meaning, verified by a compiled runtime regression. Python positional settings similarly retain their original order.

Native `QuadricCoefficients10` and final geometry remain unchanged: p stores actual symmetric matrix entries Ayz, Axz, Axy. For independent raw moments r, integral(beta)=p^T D r with D=diag(1,1,1,1,1,1,1,2,2,2). The focused `DualMonomialCoordinates10` adapter represents s=Dp/||Dp||, returns native geometry p=T s/||T s|| for T=D^-1, and differentiates dp/ds=(I-pp^T)T/||T s||. A coordinate object is never implicitly supplied as a geometry quadric.

Each strict evaluation transforms the native map into J_r_s=J_r_p dp/ds, projects multipliers afresh against [-J_g^T,s], and evaluates g_f+mu*s-J_g^T lambda. The bordered system uses H_L J_r_s+mu I, J_g J_r_s, multiplier column s, tangent row 2s^T, and tangent projection/geodesics on s. Initial sensitivity is checked after the coordinate transform. Both native and transformed numerical scales must be finite; projection cannot conceal a native norm overflow.

Raw residual is sqrt(||optimality||²+||g||²), separately exposed from the unchanged engineering normalization by max(1,||g_f||,||J_g||). Strict stopping and line-search decrease use the raw norm; default tolerance is 1e-12, alongside existing physical tolerance and positive-mu acceptance. Legacy stopping/decrease remain normalized. Restarts, failed-evaluation accounting, alpha=.5 and tmax=.5 remain intact. A failed solve has not been relabelled infeasible.

`SpinNewtonResult`, request-backed `QuadricHollowingResult`, every native event and trace carry the equation profile. Strict trace rows include native p, so their actual s is recoverable with the bound immutable profile; events already contain p when a quadric exists. Absent/invalid evaluations retain nullable raw metrics. Nonconverged/failed returned results retain strict identity.

Strict hashes are explicitly versioned: request v4, accepted trace v2, event stream/hash v2, output v3. They bind `DUAL_MONOMIAL_RAW_V1` and `literal_monomial_unit_sphere_v1`, as well as the raw residual where applicable. Legacy hash preimages retain request v3, trace v1, events v1 and output v2. Strict native JSON adds profile, coordinate identity, raw metrics and trace p; legacy JSON bytes omit the new fields. Python bindings expose enum/settings and result properties; the immutable facade consumes strict fields and defaults missing historical fields to the labelled legacy profile, without rewriting historical JSON. No Python derivative-method expansion was made.

This is an independently consistent coordinate convention, not recovered author source, a claim that the publication is wrong, or a claim of identical legacy Newton trajectories.

## TDD and mathematical evidence

All raw logs below are retained beside this report, including failures. Initial negative-control source is retained in `task-6e-red-tests.patch`.

1. Before production edits, built the existing Newton test with `--dual-pairing`, using the existing yo-yo mode and the controller's rational raw state r=(1,0,0,0,2,3,4,0,0,.5). Expected raw gradients at Qxx/Qyy/Qzz/Qxy are (-111/125,-101/125,26/25,2/25). The actual solver rejected the correctly paired stationary native p with normalized residual 0.017883. Logs: `task-6e-red-build.log`, `task-6e-red-pairing.log`. After opt-in, the exact same behavioral expectation passes.
2. Before production edits, `--raw-stopping` used a differentiable affine raw map for the supported lower-center objective with m=Lz=1e-5, p proportional to (1,0,0,-1,1e-14,0,0,0,0,0). Actual RED: `raw > 1e-12 accepted because engineering normalized residual < 1e-12`. Logs: `task-6e-red-raw-build.log`, `task-6e-red-raw.log`. GREEN proves raw>1e-12, normalized<1e-12, strict nonconvergence and unchanged legacy convergence.
3. During self-review, an added overflowing radial-native-J regression failed because projection concealed the overflowing native norm. Actual RED in `task-6e-guard-red.log`: `coordinate projection must not conceal overflowing native Jacobian norm`. The repair preserves native validation before transformation and checks the transformed norm as well. `task-6e-green-focused-2.log` confirms the repaired regression passes.

The native dual regression independently checks all three cross factors through polynomial evaluation, bidirectional conversion, sign, norm, adapter differential FD at h=3e-4/1e-4/3e-5, radial null, tangent projection and geodesic formula. It compares the actual solver's first trial against an independently assembled reduced 4-coordinate yo-yo bordered system with analytic quotient Hessian, rather than using the production adapter to construct the expected result. The nonlinear affine-map solve takes 4 evaluations and reaches raw 3.33356e-16; per-evaluation multipliers and accepted raw decrease are checked. Repeated strict event identity is `f797f4544ba273fb1aed5b55bb6dd1ea88292e657dc44fe2308cf8e5b083f16b`. Algebraic known-solution calls cover all five supported modes and are only unit validation, not geometry or author cases.

The raw line-search witness uses the supported objective f=Lz/m, p0=normalize((-0.96209441,0,0,-0.22572537,-0.15304383,0,0,0,0,0)), r0 mass=.54149167, Lz=-1.13932953, other moments zero, and r(p)=r0+A(p-p0). A is diagonal .1 except A00=-2.59280884, A33=-.93656388, A44=-3.68855061. Thus its native derivative is analytically A. These constants came from an independent NumPy search over this affine-map family using RNG seed 731; draw 161 separates the two decrease policies. At the actual solver's first trial, raw increases 1.11485→1.35535 while normalized decreases .259136→.150504. Strict rejects the trial and the legacy profile accepts it. This fixture is an algebraic line-search policy test, not a realizable geometry claim. The guarded radial-only map r0(p)=1+.5(p0²-1), at p=e0, separately proves native nonzero sensitivity cannot substitute for transformed nonzero sensitivity.

Task6d's geometry, steps, native full-map fine/entrywise/Richardson gates, signs and tolerances are unchanged. Added full 10×10 FD chains in s on the same contained/clipped asymmetric ellipsoid (both signs). Fine-step relative errors are 8.30084e-9 / 5.97380e-9 / 8.30816e-9 / 5.97266e-9 for contained-positive / clipped-positive / contained-negative / clipped-negative. Three-step errors, original Richardson results and provenance remain visible in `task-6e-green-focused-2.log`. Actual request-backed ellipsoid smoke repeats strict request/event/output identity; unsupported-plane differential rejection retains strict identity and null raw metric. No convergence was manufactured for those geometry smoke cases.

Python native tests exercise valid and invalid profile configuration, enum roundtrip, strict finite and failed evaluations, null metrics, strict trace/event identity, immutable trace coefficients, native result properties, repeated strict diagnostics, and unchanged historical serialization.

## Compatibility evidence

Before rebuilding production, executed the actual BASE Newton binary and captured its accepted-trace hash `4325d8f9b111a6925b47df38cd87e0f9023360b96efec3164a1fe20329e81d18`; it is now a regression expectation, not a comparison between two post-change defaults. BASE native Newton executable SHA was `e4c99ab66b7e4a98c504422196ba914dc8dff31925d00734ac5a7413357de181`, primary cache SHA `2c3f5a72c82d0ba5d8b5d0cd8869d13d89ea8a21e0daa405e88242e1778e517c`.

`task-6e-baseline-legacy.json` was emitted by the actual pre-change native binding for the failed constant-quadric request used by the Python test. Its exact SHA-256 (including print's extra terminal newline) is `c744e1daf10d8349b42801e787a9b54193464facfb62725dcc14b098fcf04dd8`. `cmp` with `task-6e-current-legacy.json` returned 0 after the new serializer/binding build, and the Python regression independently checks the same captured literal fingerprint. It binds input/event/output hashes and JSON bytes, not only status. Historical `experiments/` has no diff against BASE (`task-6e-historical-output-check.log`, exit 0).

## Verification commands and results

All work and builds were owned by this implementer in the strict worktree; build concurrency never exceeded -j4, primary and isolation builds were sequential, and no controller native launch overlapped them. Focused tests were used during iteration. The full primary suite is run once after the final all-target build.

```sh
cmake -S . -B build/strict-v2-spin
cmake --build build/strict-v2-spin --target quadric-hollowing-solver-regression -j4
build/strict-v2-spin/quadric-hollowing-solver-regression --dual-pairing
build/strict-v2-spin/quadric-hollowing-solver-regression --raw-stopping
build/strict-v2-spin/quadric-hollowing-solver-regression --dual-coordinates
cmake --build build/strict-v2-spin --target quadric-hollowing-solver-regression quadric-discrete-moment-regression _make_it_stand -j4
ctest --test-dir build/strict-v2-spin -R '^spin-faster-(newton.*|discrete-moments|python-api)$' --output-on-failure -V
/usr/bin/time -v cmake --build build/strict-v2-spin -j4
/usr/bin/time -v ctest --test-dir build/strict-v2-spin --output-on-failure
cmake -S . -B build/task-6b-faster-only -DMAKE_IT_STAND_BUILD_SPIN_IT_FASTER=ON -DMAKE_IT_STAND_BUILD_SPIN_IT_2014=OFF
/usr/bin/time -v cmake --build build/task-6b-faster-only --target quadric-hollowing-solver-regression quadric-discrete-moment-regression -j4
/usr/bin/time -v ctest --test-dir build/task-6b-faster-only -R '^spin-faster-(newton.*|discrete-moments)$' --output-on-failure -V
git diff --check
```

| Check | Result | Wall time | Maximum RSS (KiB) | Raw evidence |
|---|---|---:|---:|---|
| Final focused native/Python | 6/6 | 8.87 s | not timed | `task-6e-green-focused-2.log` |
| Final primary all-target build | exit 0 | 13.22 s | 1,615,128 | `task-6e-final-all-build.log` |
| Full primary registered CTest, run once | 186/186, no skipped/Not Run | 192.88 s | 4,253,840 | `task-6e-final-full.log` |
| Faster ON / 2014 OFF isolation build | exit 0 | 61.84 s | 3,343,544 | `task-6e-isolation-build.log` |
| Isolation relevant native tests | 5/5 | 8.65 s | 17,908 | `task-6e-isolation-native.log` |

The suite grows from baseline 183 to 186 through three explicitly registered Newton contract invocations; existing test targets remain registered. `task-6e-primary-sha256.log` records all changed source/header/test SHA-256 values and final binaries/core/cache. `task-6e-primary-flags.log` records actual cache options: Faster ON, 2014 ON, exact geometry ON, Python ON, reachability ON, paper Python tests ON, viewer OFF, PyBullet tests OFF, compiler observation OFF, BUILD_TESTING ON.

| Final primary artifact | SHA-256 |
|---|---|
| Newton executable | `bd77b54c83736c6ce20ff602e87620e0f17e9409999caf00f674d732d52452bf` |
| Discrete-map executable | `43e368d278df21b752efb99c3c107a569eb3a872778d3e25fb796bc7b42a7731` |
| Native Python extension | `ed48560b3ffe9e7a11fb0d94b2fa3fc312893a5d7694f12895d5d1811cfb59b2` |
| Core archive | `f300541613592815ff47af249365499eca9790f8c176eb10d98048cfc3bb69a9` |
| CMake cache | `2c3f5a72c82d0ba5d8b5d0cd8869d13d89ea8a21e0daa405e88242e1778e517c` |

`task-6e-isolation-sha256.log` records Newton executable `1155da93833b6cd28eb73be0292bd0b5ecb428f4334fb3a8b28756f59d4002df`, discrete executable `7f8dc970f85908aa422f511096e2eefc7651af04bbcd88fd117992b96970cef6`, core archive `e10debb806d08c4fb691edce78c89c76ae71f2f4074770f9033ec7839d1ae506`, and cache `bf9ce45dc0aee97188c0bd023182f18bf0b1d7fc927ee08ad54d38be72718127`. Actual isolation flags are Faster ON, 2014 OFF, exact geometry checks OFF, Python OFF, reachability OFF, viewer OFF, paper Python tests OFF, PyBullet tests OFF, compiler observation OFF, BUILD_TESTING ON (`task-6e-isolation-flags.log`). Native coordinate/discrete tests therefore require no 2014 implementation or Python feature.

Existing warnings: configure emits libigl's old-CMake compatibility deprecation; rebuild of the existing Boolean evaluator emits Boost's global bind-placeholder deprecation note. They are retained in raw logs, not hidden or fixed by dependency changes. The final incremental all-target build itself has no warning lines. No missing/disabled optional capability is presented as a paper-method success; PyBullet stays OFF as in the approved baseline.

After final verification, `sha256sum -c` confirms every source/primary artifact and isolation artifact still matches the recorded identity (`task-6e-final-source-identity-check.log`). Source/test/report whitespace checks are clean. The unrestricted staged whitespace check reports 19 raw-output trailing-space lines (8 CMake configure output lines, 9 verbose CTest/Python output lines, 2 literal context lines in the retained RED patch), plus the intentional extra EOF newline in each of the two byte-comparison JSON captures. Those raw evidence bytes are preserved, not reformatted.

## Changed files and self-review

Production: `CMakeLists.txt`; new `include/spin/DualMonomialCoordinates.h`, `src/spin/DualMonomialCoordinates.cpp`; `include/spin/SpinDesignRequest.h`, `include/spin/QuadricHollowingSolver.h`; `src/spin/SpinDesignRequest.cpp`, `src/spin/QuadricHollowingSolver.cpp`; `src/api/BalanceModel.cpp`; `bindings/pybind_module.cpp`; `python/make_it_stand/spin.py`. Tests: new `tests/dual_monomial_newton_regression.cpp`; `tests/quadric_hollowing_solver_regression.cpp`; `tests/quadric_discrete_moment_regression.cpp`; `tests/python/test_spin_it_faster.py`. This report and the explicitly scoped raw Task6e evidence are force-added individually under ignored SDD.

Self-review read the complete diff and new files: checked native p versus s at all evaluator, bordered-system, projection, geodesic, trace and output boundaries; native and transformed numerical guard ordering; raw stopping and decrease through behavioral adversaries; nullable failed metrics; strict identity on early nonconvergence; legacy hash/JSON literal fingerprints; and append-only aggregate/Python source compatibility. The existing solver and BalanceModel files are large; edits remain scoped and the coordinate adapter is focused. Independent reviewer belongs to the controller after handoff.

Limitations: no full-reproduction or feasibility certification follows; geometry smoke need not converge; analytic maps and known-solution fixtures are not original experiments; discrete ellipsoid differentiation retains Task6d applicability/topology restrictions, and default legacy evaluator selection is unchanged. ABI binary compatibility is not asserted for C++ structs extended with public fields; native source/aggregate and Python compatibility are tested. External preexisting warnings and final exact build identity are recorded below.
