# Parabolic-cylinder finite geometry and reduced Boolean derivative — 2026-09-24

## Scope

Implemented the rank-one real parabolic-cylinder path from classification
through finite triangulation, fixed-stratum source motion, and the actual CGAL
Boolean-to-moment differential. The equation `x^2 + 2y = 0` is a synthetic
method fixture. Its normalized coefficient stratum has six tangent dimensions;
the implementation does not pad its Jacobian to ten columns or claim a
rank-deficient Equation 14 optimizer.

Rank-one quadrics without a transverse linear term remain explicitly rejected
as a separate parallel/repeated/empty-plane family. This avoids silently
labeling those different degeneracies as parabolic cylinders.

## Implementation

- Added `ParabolicCylinder` surface/subtype classification and the
  `parabolicCylinderGridResolution` control.
- Triangulated a graph patch in principal-curvature and cylinder-axis
  coordinates. It is closed by translation along the transverse linear axis;
  extruding along the cylinder axis caused degenerate side faces and failed the
  strict boundary audit.
- Added the six-dimensional tangent constraints for a normalized rank-one
  quadratic block.
- Derived source-vertex velocities for curvature-axis motion, transverse
  linear direction, graph height, and translated closure vertices.
- Enabled the reduced discrete Boolean differential and public evaluator
  result path while retaining the true `10 x 6` Jacobian.

## Verification

The closed finite boundary passes strict manifold, winding, and intersection
audits. Source-motion finite differences under fixed-stratum retraction have
relative errors:

| Step | Relative error |
|---:|---:|
| `1e-4` | `2.34218496e-8` |
| `3e-5` | `2.10795372e-9` |
| `1e-5` | `2.34441311e-10` |

The clipped, offset design-domain Boolean output contains 30 constructed
intersection vertices. Its reduced ten-moment Jacobian agrees with finite
differences:

| Step | Relative error |
|---:|---:|
| `1e-4` | `6.56565715e-8` |
| `3e-5` | `5.90798518e-9` |
| `1e-5` | `6.52202591e-10` |

The regression also verifies the public unified evaluator exposes the reduced
result separately and does not mark an ambient `10 x 10` derivative as
defined. Its requested grid resolution is now validated and included in the
evaluator request identity (together with the hyperbolic-cylinder grid
resolution); both controls are exposed in the native Python binding and the
immutable Python settings object. A regression first reproduced the stale
identity failure, then passed after the fix. Legacy-default requests retain
their historical base-request hash; non-default extended identities use the
versioned `QuadricEvaluatorRequest-v2` encoding.

## Remaining limits

- No rank-deficient stratum optimizer or full ambient Newton derivative.
- No transitions across surface type, Boolean topology, or principal-axis
  chart; unsupported events remain fail-closed.
- Parallel/repeated/empty planes, intersecting planes, cones, repeated
  eigenspaces, and other degenerate quadrics remain separate gaps.
- No original Faster input or paper experiment is represented by this synthetic
  method validation.
