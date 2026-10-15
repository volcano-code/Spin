# Reduced paraboloid Boolean-moment differential — 2026-09-24

## Result

Added a fixed-stratum differential for the native finite paraboloid source
mesh, including the motion of its translated closure, and connected it to the
existing fixed-topology Boolean-to-moment differential. It returns a
10-by-8 Jacobian: ten geometric-moment rows and eight coordinates in the
paraboloid's normalized coefficient tangent. It does not pad those eight
columns to ten and does not enter the full-rank Equation-14 solver.

The source-motion fixture uses a rank-two paraboloid with distinct nonzero
curvatures, a fixed radius, and the native 16-by-16 parameter grid. Its source
vertex velocity was compared with independently retracted, retriangulated
central differences:

| Step | Relative source-vertex Jacobian error |
|---:|---:|
| `1e-4` | `6.0024e-8` |
| `3e-5` | `5.4019e-9` |
| `1e-5` | `6.0108e-10` |
| `3e-6` | `5.9271e-11` |

The actual CGAL Boolean fixture uses the paraboloid
`x^2 + 2 y^2 - z - 0.2 = 0`, density 1, a frozen translated cube design
domain, and 46 constructed source/design intersection vertices. The reduced
10-by-8 mass-moment Jacobian is checked against the actual Boolean result under
independent fixed-stratum retractions; topology is compared by sorted
birth-feature provenance rather than output vertex ordering:

| Step | Relative Boolean-moment Jacobian error |
|---:|---:|
| `1e-4` | `1.4382e-8` |
| `3e-5` | `1.2953e-9` |
| `1e-5` | `1.4425e-10` |

The test rejects topology changes, unavailable differentials, wrong chart
dimensions, and non-deterministic source identities. Repeated-eigenvalue or
other ill-conditioned charts continue to fail closed.

The source-motion adapter also checks the deterministic eigenvector sign-pivot
chart. If the two largest absolute components differ by at most `1e-10` and
have opposite signs, the reduced derivative is rejected with an explicit
sign-pivot-boundary diagnostic. Same-sign ties remain admissible because they
do not change the canonical vector sign; a rotated paraboloid control checks
this case. Regressions cover the unstable pivot in both a paraboloid null axis
and an elliptic-cylinder axis. This is a local branch guard; it neither
handles repeated eigenspaces nor certifies a neighborhood around a
nondegenerate point.

## Verification

- `quadric-source-motion-regression`: passed.
- `spin-faster-source-motion`, `spin-faster-discrete-moments`, and
  `spin-faster-stratum-tangent`: **3/3 passed**.
- The expanded triangulation/source-motion/discrete-map/Newton integration
  CTest group passed **11/11**.
- Both paraboloid and elliptic-cylinder sign-pivot-boundary rejection checks
  passed in the rebuilt `quadric-source-motion-regression` executable.
- A same-sign equal-magnitude pivot control remains accepted.
- Strict evidence:
  `build/faster-paraboloid-stratum-evidence-20260924-c/`, status
  `SYNTHETIC_REDUCED_PARABOLOID_BOOLEAN_DIFFERENTIAL_VALIDATED`, freshness
  `CURRENT`.

## Limits

This validates local derivatives on one fixed, simple-spectrum paraboloid
stratum and one fixed Boolean topology. It does not provide a ten-dimensional
ambient derivative at the rank-deficient point, optimize Eq. 14 over the full
coefficient sphere, handle transitions into/out of the paraboloid stratum,
cover planes, CAD, or an author input. Therefore it is not a complete Faster
solver or a reproduced paper case. Full-rank solver, planar/other rank-
deficient solver integration, Cherchi/CGAL cross-backend comparison, and
author experiments remain outstanding.
