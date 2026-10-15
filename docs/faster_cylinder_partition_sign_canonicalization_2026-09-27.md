# Faster elliptic-cylinder exact partition follow-up — 2026-09-27

## Failure and root cause

The sealed reduced-stratum cylinder test evaluated the same design cut using
coefficients `p` and `-p`. The exact material partition audit reported both a
positive-volume overlap and missing domain content. The geometric point sets
were numerically close, but not the same represented binary64 polyhedron: the
symmetric eigensolver received quadratic blocks `Q` and `-Q`, whose signed
eigenvalue order differed. Consequently, the independently generated radial
bases differed at roundoff scale. This is material to the exact audit, which
interprets each binary64 coordinate as an exact represented value.

## Repair

For elliptic-cylinder triangulation only, canonicalize the coefficient sign
from the quadratic-block trace before computing the eigenbasis and source
boundary. The original polynomial remains authoritative for the independent
positive/negative Boolean-side decision. This keeps the retained and carved
boundaries byte-identical without weakening the exact partition checker or
changing the requested material side.

The C++ triangulator regression now requires `p` and `-p` to produce identical
vertex and face arrays, in addition to its prior rounded geometric-facet
comparison.

## Verification

- The previously failing
  `test_reduced_stratum_cylinder_candidate_seals_after_exact_partition` now
  passes; it exports the cylinder candidate to B-rep and STEP and verifies the
  material-moment gate.
- `quadric-triangulator-regression` passes, including the new binary64
  representation-equality assertion.
- `tests/python/test_occt_step_roundtrip.py`: **21 passed**.
- `tests/python/test_faster_material_bridge.py`: **27 passed**.

This closes the specific sign-invariant cylinder partition regression. It does
not close the separate elliptic-paraboloid CAD-versus-discrete-mesh moment
discrepancy documented in
`faster_paraboloid_cad_diagnostic_2026-09-27.md`, nor does it establish the
complete Faster algorithm/application suite or either paper's original
experimental reproduction.
