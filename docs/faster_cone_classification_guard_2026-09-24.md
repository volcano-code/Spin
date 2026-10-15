# Faster singular-quadric classification guards — 2026-09-24

## Change

The quadric classifier now returns an explicit `Cone` type for a full-rank
indefinite quadratic whose translated constant level is zero. Such an input
has an apex, so the regular hyperboloid inverse-Gauss truncation routine is
not applicable. Triangulation fails closed with an `invalid_argument`
diagnostic identifying the unsupported degenerate cone instead of eventually
reporting that two truncation loops could not be constructed. The fixed-stratum
tangent evaluator likewise rejects the cone type rather than returning a
normalization-only basis that would falsely imply a valid cone tangent.

The same classifier now distinguishes a rank-two indefinite quadric with no
linear component along its null axis as either a `HyperbolicCylinder` (nonzero
centered level) or `IntersectingPlanes` (zero centered level). Both types have
explicit triangulation diagnostics. The coefficient-only fixed-stratum tangent
now supports a simple-spectrum hyperbolic cylinder using the rank-two null-axis
and axis-orthogonal-linear constraints; cones and intersecting planes still
fail closed there.

Regressions use `x^2 + y^2 - z^2 = 0` for a cone,
`x^2 - y^2 - 1 = 0` for a hyperbolic cylinder, and `x^2 - y^2 = 0` for two
intersecting planes. They confirm typed classification, explicit triangulation
rejection, and refusal to fabricate unsupported tangent bases. The standard
one-sheet-hyperboloid truncation and existing plane/elliptic-cylinder/
paraboloid/ellipsoid tangent fixtures remain positive controls.

## Verification

- The initial cone regression was first run against the old routing and escaped as
  the unrelated-looking runtime error `hyperboloid truncation did not produce
  exactly two closed loops`.
- After adding typed classification and fail-closed consumers, the targeted
  CTest selection covering triangulation, source motion, discrete moments,
  stratum tangent, and unified moment evaluation passed **5/5**.
- `git diff --check` passed for the implementation and regression test.

## Follow-up: hyperbolic-cylinder coefficient tangent

The fixed-stratum tangent test now includes `x^2-y^2-1=0`. It verifies the
classified subtype, three independent constraint rows, seven-dimensional
tangent chart, finite orthonormal basis, and translation tangent control. The
focused CTest `spin-faster-stratum-tangent` passes **1/1**.

## Scope limit

This initially implemented classification, safe-failure semantics, and a
coefficient-only hyperbolic-cylinder tangent chart. A later same-day increment
adds finite hyperbolic-cylinder triangulation and reduced source/Boolean-moment
derivatives; see
[`faster_hyperbolic_cylinder_stratum_2026-09-24.md`](faster_hyperbolic_cylinder_stratum_2026-09-24.md).
It still does not add a rank-deficient optimizer or reproduce a Faster paper case.
A later increment also adds the simple parabolic-cylinder finite geometry and
reduced Boolean derivative; see
[`faster_parabolic_cylinder_stratum_2026-09-24.md`](faster_parabolic_cylinder_stratum_2026-09-24.md).
Parallel/double planes, other singular quadrics, repeated eigenspaces, and
transitions between strata remain separate algorithm gaps.
