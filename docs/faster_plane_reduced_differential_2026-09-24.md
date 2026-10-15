# Faster plane-stratum discrete differential — 2026-09-24

## Scope and status

This change adds a **reduced, fixed-plane-stratum** source-motion and
Boolean-to-moment differential. It is not a ten-dimensional coefficient
Jacobian, does not enter the Eq. 14 full-rank Newton solver, and is not an
author-case result.

The chart has three coordinates (two plane-normal directions and one offset).
The source mesh is the triangulator's four-corner plane patch plus its
extrusion. The returned `10 x 3` Jacobian differentiates the same fixed
topology Boolean material map used by the existing restricted cylinder and
paraboloid paths.

## Differential and branch handling

For normalized plane coefficients `(c, l)`, the source motion differentiates
the native plane point `-c l / ||l||^2`, the local orthogonal frame, all four
patch corners, and the translated extrusion vertices. The result includes a
reduced coefficient basis and per-vertex provenance. The Boolean differential
rebuilds and byte-compares the source triangulation before transporting motion
through the fixed-provenance Boolean boundary.

Plane patches now share an explicit two-chart tangent frame with their
derivative: chart 0 uses the x-axis reference for `|n_x| < 0.9`; chart 1 uses
the y-axis reference otherwise. The chart seam is tested and returns an
explicit undefined-derivative diagnostic. This avoids the old
`unitOrthogonal()` frame discontinuity at the z-axis while retaining a
well-conditioned frame near coordinate-axis normals. The chart seam is an
implementation-coordinate limitation, not a claim that the underlying plane
geometry is singular.

## Verification

The regression was written first and failed for the expected missing plane
stratum source motion, then again for the missing plane Boolean differential.
After implementation, source and Boolean-map Jacobians were checked against
independent central differences at three step sizes while requiring unchanged
mesh topology, radius branch, and Boolean provenance signature.

Observed relative errors:

| Check | Step sizes | Relative error |
|---|---|---:|
| Plane source vertices | `1e-4`, `3e-5`, `1e-5`, `3e-6` | `8.16e-9`, `7.34e-10`, `8.30e-11`, `2.74e-11` |
| Plane Boolean material moments | `1e-4`, `3e-5`, `1e-5` | `2.00e-8`, `1.81e-9`, `2.18e-10` |
| Axis-aligned plane source vertices | `1e-4`, `3e-5`, `1e-5` | `5.00e-9`, `4.49e-10`, `4.85e-11` |

`ctest -R '^spin-faster-(source-motion|plane-family|quadric-triangulator|discrete-moments)$' -j1`
passed **4/4** on build
`build/strict-algorithms-exact-current-20260922` (21.11 seconds). The plane
family regression additionally checks axis-aligned source moments at a
`1e-5` finite-difference step and verifies that the chart seam is rejected.

The broader `ctest -L spin-it-faster` attempt was interrupted after its
registered unknown-seed sweep was found to run all five seeds, contrary to the
current plan's explicit deferral of the five-seed batch. Its partial run is
not counted as passed; no process remains running. The bounded six-test result
above is the accepted verification for this change.

## Remaining gaps

- Plane reduced derivatives do not yet participate in a stratum-aware
  optimizer.
- Fixed-plane Boolean derivatives still require unchanged topology and a
  locally stable source-frame chart.
- Hyperbolic cylinders, cones and other degenerate surface families remain
  unsupported or explicitly rejected.
- Cherchi/CGAL entity discrepancies remain unresolved.
- Faster author inputs, critical parameters, original case matrix, CAD cases,
  and paper-scale statistics remain unverified or blocked.
