# Faster paraboloid CAD bridge diagnostic — 2026-09-27

## Scope

This is a method-path diagnostic, not a paper-case result. It exercises one
sealed, reconstructed paraboloid candidate through the analytic OCCT CAD
material bridge. The pre-registered CAD/material acceptance limit remains
`1e-6`; no tolerance was relaxed.

## Findings

The first candidate exposed a sign-handling defect in the analytic paraboloid
carrier. The bridge canonicalized a negative-definite quadratic block, but the
carrier did not consistently preserve both (a) the input polynomial's negative
halfspace and (b) the curvature orientation in its carrier frame. The resulting
CAD material had a normalized mass-moment error of about `0.295` and the wrong
sign for the first z moment. The carrier now keeps transverse eigenvalues in
the canonical positive-quadratic frame, derives its axial direction from the
actual input sign, and flips the meridian curvature when the polynomial is
globally negated.

After this fix, the retained side agrees qualitatively with the sealed candidate,
but finite-resolution triangulation is not equivalent to the exact analytic CAD
surface. With the same frozen geometry and coefficients, the maximum
normalized material moment differences were:

| Paraboloid grid resolution | Maximum normalized moment error | Result |
|---:|---:|---|
| 32 | `1.2755738841e-2` | rejected |
| 64 | `3.3774490626e-3` | rejected |
| 128 | `8.2780459796e-4` | rejected |

The approximately fourfold reduction per grid doubling is consistent with a
second-order surface-discretization error. Extrapolating this trend suggests
that reaching `1e-6` by uniform refinement alone would require a grid in the
thousands and is not yet an acceptable validated path.

OCCT's adaptive volume integration reported `1.7174260632e-10` when requested
at `1e-10`. Its requested integration tolerance is now `1e-12`; the independent
acceptance check remains `<=1e-10`. This did not change the CAD/material result
criterion.

## Evidence and status

- The focused regression passes by asserting that a mismatched coarse
  paraboloid CAD result is rejected and that no B-rep, STEP, or identity sidecar
  is published.
- The analytical paraboloid carrier now reaches the common multi-material
  moment gate, but this does **not** pass CAD equivalence at the locked `1e-6`
  limit.
- `tests/python/test_occt_step_roundtrip.py`: 20 passed, 1 failed. The failure
  occurs before CAD export while sealing a reduced-stratum elliptic-cylinder
  candidate: the exact partition checker reports a missing whole-domain region
  and a positive-volume overlap between retained and carved parts. This is
  outside the OCCT integration path changed here and remains a separate
  regression. A later same-day follow-up isolated and fixed that regression;
  see `faster_cylinder_partition_sign_canonicalization_2026-09-27.md`.

## Remaining work

Do not loosen the moment threshold or emit a trusted CAD sidecar for the
paraboloid candidate. To close the CAD route, the solver/geometry contract must
either provide a sealed boundary whose representation agrees with the exact
CAD quadric within the existing tolerance, or define and independently verify
separate analytic-CAD and discrete-solver moment records. Until then,
elliptic-paraboloid CAD material equivalence remains blocked.
