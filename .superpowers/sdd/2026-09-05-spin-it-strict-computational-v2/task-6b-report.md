# Task 6b report: analytic finite inverse-Gauss point differential

## Status and bounded scope

Implemented the fixed-normal analytic coefficient differential of the existing
finite inverse-Gauss point map on base
`328d4c174c444109eff7386d74afbb08186a7054`. The actual ellipsoid and
hyperboloid triangulator now evaluates finite inverse-Gauss points through the
same prepared helper. Point evaluation retains its positive-square-root branch,
coefficient sign behavior, denominator cutoff `1e-15`, `wSquared` cutoff
`1e-14`, and bit-compatible representative geometry.

This is only the finite, full-rank, fixed-normal point-map increment. It does
not differentiate ideal vertices, truncation/cap vertices, planar,
paraboloidal, or cylindrical maps, Boolean intersection vertices, changing
topology, surface moments, or the complete discrete material-moment map. It
therefore does not fix or supersede the existing whole-map diagnostic error
`0.00658869`, and it does not complete Task 6/M4 or any original-paper case.

## Source changes and public contract

- Added `QuadricInverseGaussMap`, prepared once per normalized quadric. Its
  ordinary `point()` path caches the inverse quadratic matrix and denominator;
  the triangulator does not solve derivative systems while generating points.
- Added independent point and differential status enums, diagnostics, and
  `std::optional` point/Jacobian payloads. A valid point remains available when
  a requested differential is ill-conditioned; undefined states never expose
  a zero matrix as a successful Jacobian.
- Added the opt-in analytic normalized-coefficient `3x10` Jacobian and a native
  Faster regression. The source and test are registered only in the Faster
  CMake feature block and have no Spin-It 2014 dependency.
- Removed the triangulator's private duplicate point formula. A prepared common
  map is reused for all icosphere vertices and truncation-edge bisection points.

## Independent derivation implemented

For the current finite point `x`, fixed normal `n`, quadratic matrix `A`, and
linear vector `b`, each normalized coefficient direction solves

```text
[ A            -n ] [dx] = [-dA*x - db/2]
[(2*A*x+b)^T    0 ] [dt]   [-(x^T*dA*x + db^T*x + dc)].
```

The implementation uses a unit-scaled copy of the fixed normal only to make
the bordered-system conditioning independent of positive normal rescaling;
the authoritative point path still uses the caller's normal and historical
formula exactly. Native coefficient columns are
`[c,bx,by,bz,Axx,Ayy,Azz,Ayz,Axz,Axy]`. Thus each cross column changes both
symmetric matrix entries and its bottom equation contains `2yz`, `2xz`, or
`2xy`. After the ten raw bordered solves, right multiplication by
`I-p*p^T` is applied exactly once for perturb-and-renormalize derivatives of
the normalized coefficient vector `p`.

The derivative reciprocal-condition cutoff `1e-12`, solve residual cutoff
`1e-10`, and invalid-normal norm cutoff `1e-15` are implementation-default
engineering safeguards. They are not paper-published tolerances or rigorous
forward-error certificates. The existing triangulator full-rank classifier
threshold and historical point-map cutoffs remain unchanged.

## TDD and focused evidence

RED command after adding only the native regression and its CMake target:

```text
cmake -S . -B build/strict-v2-spin
cmake --build build/strict-v2-spin \
  --target quadric-inverse-gauss-map-regression -j2
```

It failed at the expected missing production header
`spin/QuadricInverseGaussMap.h`. The complete failure is retained in
`task-6b-red.log`.

Initial focused GREEN built the new regression and the existing triangulator
and moment-evaluator regressions with `-j2`, then ran:

```text
ctest --test-dir build/strict-v2-spin \
  -R '^(spin-faster-inverse-gauss-map|spin-faster-quadric-triangulator|spin-faster-quadric-moments)$' \
  --output-on-failure
```

Observed `3/3` passed. Logs are `task-6b-green-initial-build.log` and
`task-6b-green-initial-ctest.log`.

The final direct native regression reported:

```text
ellipsoid h=3e-4: 4.31344e-07
ellipsoid h=1e-4: 4.79273e-08
ellipsoid h=3e-5: 4.31273e-09
hyperboloid h=3e-4: 1.69238e-07
hyperboloid h=1e-4: 1.88043e-08
hyperboloid h=3e-5: 1.69223e-09
explicit analytic-zero entries: 23
maximum analytic-zero absolute error: 0
normalized radial-null norm: 2.61846e-16
```

These are relative Frobenius errors using `max(1, ||FD||)` for the smooth
asymmetric translated/cross-term ellipsoid and admissible hyperboloid fixtures.
All three preset errors are below `1e-6` and show the expected central-
difference convergence trend. The separate sphere oracle checks all ten
columns analytically, including nonzero linear and native cross terms, without
using finite differences or the bordered solver to construct expectations.
Repeated sphere eigenvalues, point/radius response, positive normal rescaling,
normal sign and complemented-quadric behavior are covered.

Invalid zero/NaN/Infinity normals, tiny-normal historical cutoff behavior,
zero/nonfinite coefficient construction, singular quadratic matrices,
degenerate denominators, ideal/inadmissible hyperboloid normals, and a
near-ideal valid point with an undefined ill-conditioned differential are
covered. The successful finite-difference fixtures establish these test cases,
not a global conditioning guarantee.

## Geometry compatibility

Before extraction, the representative truncated-hyperboloid triangulator hash
was
`e71f44edb04749919a3b2a3d2201ea1e21fc39b57bbef768fa2f768ba4cce847`.
After extraction and the final build it is byte-identical. The contained-sphere
moment-evaluator geometry hash remains
`0f00208452856b310c3735782290374aeac11b185b640fe0a157d102ef387101`.
No historical result or scientific output was regenerated or rewritten.
The direct output and executable-mode check are retained in
`task-6b-final-direct.log`.

One controller read-only launch overlapped the owner's triangulator relink and
temporarily returned exit 126/permission denied while the executable was being
replaced. The completed build restored normal mode `775`; subsequent direct
execution passed with the same hash above. This is recorded as a concurrent
relink race, not as a geometry-test failure.

## Faster-only and final suite evidence

A separate exact-source build at `build/task-6b-faster-only` configured
`MAKE_IT_STAND_BUILD_SPIN_IT_FASTER=ON`,
`MAKE_IT_STAND_BUILD_SPIN_IT_2014=OFF`, and `BUILD_TESTING=ON`. The initial
fresh configure's own DownloadProject performed a CGAL checkout/setup; that
observed setup is retained without claiming the tree was previously absent.
Current CGAL HEAD was then verified as the required pin
`f7c3c8212b56c0d6dae63787efc99093f4383415`, and the final configure reused it
with `LIBIGL_SKIP_DOWNLOAD=ON`; the project's dependency audit remained
enabled. The map and triangulator tests passed `2/2`. See
`task-6b-faster-only.log` and `task-6b-faster-only-final.log`.

Final current-source verification used `build/strict-v2-spin` with Faster,
Spin-It 2014, Python, reachability, exact geometry, and testing enabled. Local
pinned dependency reuse was explicitly selected without changing pins or
features.

```text
cmake --build build/strict-v2-spin -j2
ctest --test-dir build/strict-v2-spin -L spin-it-faster --output-on-failure
ctest --test-dir build/strict-v2-spin --output-on-failure
```

The all-target build exited zero with `182` registered tests. The Faster suite
passed `14/14` in 20.64 seconds. The full registered suite passed `182/182` in
188.55 seconds with no failed or Not Run tests. Logs are
`task-6b-full-build.log`, `task-6b-faster-suite.log`, and
`task-6b-full-suite.log`.

## Self-review and remaining gaps

Reviewed coefficient order, symmetric cross-term factors, bordered-system
signs, positive-root point semantics, the single normalization tangent chain,
rank/conditioning/residual failure paths, optional payload contracts, fixed-
normal scaling, common-path consumption, cached ordinary-map cost, CMake
feature isolation, and representative geometry hashes. `git diff --check`
was clean for production, tests, CMake, and this report before commit; raw
command logs retain their emitted trailing spaces verbatim.

Ideal point derivatives described by the Faster paper remain unimplemented in
this finite-only task, as do caps, intersections, topology changes, full
triangulation derivatives, and the complete discrete moment Jacobian. The
separate legacy overflow false-convergence issue from Task 6a review is also
unchanged; it was not repaired as scope creep. No Eq. 14 constant, acceptance
rule, solver behavior, original case, release evidence, merge, push, or tag is
part of this change.
