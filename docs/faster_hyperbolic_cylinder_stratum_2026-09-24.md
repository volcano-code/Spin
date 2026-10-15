# Hyperbolic-cylinder finite geometry and reduced Boolean derivative — 2026-09-24

## Scope

Implemented the simple-spectrum hyperbolic-cylinder branch across native
triangulation, fixed-stratum source motion, and the actual CGAL Boolean-to-moment
differential. The input stratum is rank two, so its coefficient chart has
dimension 7; this work does not fabricate a 10-column Jacobian or claim a
general Equation 14 optimizer for rank-deficient quadrics.

The finite geometry is a regular grid over both real branches, capped outside
the chosen truncation radius. It is an auditable finite construction for the
project's Boolean pipeline, not an exact representation of an unbounded
hyperbolic cylinder or a paper asset.

## Implementation

- Added `hyperbolicCylinderGridResolution` and native routing for the classified
  `HyperbolicCylinder` surface type.
- Constructed both branches with a common material-side orientation and
  exterior cap fans. Reversing one branch caused signed-volume cancellation
  and discontinuous global winding under coefficient perturbations; the
  regression now guards stable face orientation.
- Derived source vertex velocities in the true 10-by-7 fixed-stratum chart,
  including eigenframe motion, centered level/scales, longitudinal branch
  motion, and cap-centroid motion.
- Enabled the reduced discrete Boolean moment differential for this type.
- Corrected Boolean source-feature handling: dimension 3 is the no-feature
  sentinel. A single-source output vertex now uses the barycentric velocity of
  its source feature (or zero when supported only by the fixed design), while
  genuine two-input edge-edge contacts remain fail-closed.

## Verification

The source mesh has 340 vertices, 672 faces, and two boundary loops before
closing; its closed result passes the strict boundary audit. Independent source
finite differences give relative errors:

| Step | Relative error |
|---:|---:|
| `1e-4` | `4.4194469e-8` |
| `3e-5` | `3.9775333e-9` |
| `1e-5` | `4.4196425e-10` |

For the clipped non-symmetric design domain, the actual CGAL Boolean output has
32 constructed intersection vertices. The ten-moment reduced Jacobian agrees
with central finite differences:

| Step | Relative error |
|---:|---:|
| `1e-4` | `2.6447713e-8` |
| `3e-5` | `2.3816459e-9` |
| `1e-5` | `2.7074272e-10` |

The registered CTest selection `spin-faster-(quadric-triangulator|source-motion|discrete-moments|stratum-tangent|quadric-moments)` passed **5/5** after the implementation. The `spin-faster-source-motion` executable also passed, including the hyperbolic-cylinder source and Boolean finite-difference checks.

The current linked-worktree build repeated that same CTest selection on
2026-09-24 and passed **5/5**. This exposed and corrected a stale triangulator
regression: the implementation already produced the capped hyperbolic-cylinder
boundary, while the old assertion still expected an exception. The regression
now checks the generated subtype, two branch loops, interface residual, and
strict closed-boundary audit. The geometry remains a finite capped
method-validation construction, not an exact unbounded surface.

The corresponding Faster Python API suite also passes **15/15** on the current
build. Its historical serialized-JSON digest was refreshed only after two
independent evaluations produced the same bytes; the output changed because
the unsupported constant/zero-quadric diagnostic was clarified.

## Remaining limits

- No rank-deficient stratum optimizer or full ten-dimensional Newton derivative.
- No certification across surface-type, topology, truncation, or eigenframe-chart
  transitions; those events still fail closed.
- Repeated eigenvalues, intersecting/parallel plane families, cones, and other
  degenerate quadrics remain separate gaps.
- No Faster Table 1/Figure result or author input is covered by this synthetic
  method test.
