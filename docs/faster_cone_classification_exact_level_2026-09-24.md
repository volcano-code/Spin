# Faster exact cone versus near-cone classification — 2026-09-24

## Scope

This change corrects the algebraic classification of full-rank quadrics near a
cone. It is algorithm-kernel validation on synthetic coefficients, not an
author input or a paper-case result.

## Failure and correction

The old classifier used `abs(shiftedConstant) <= 1e-12 * scale` to label a
full-rank indefinite quadric as a cone. That conflated exact degeneracy with a
regular quadric whose nonzero centered level is small. A regression with
`A = diag(1,1,-1)` and centered level `+1e-14` failed before the change; the
opposite `-1e-14` side was also added to verify subtype selection. The same
threshold also rejected a nonempty definite quadric with centered level
`-1e-14` (and its sign-inverted equivalent), although a mathematically real
ellipsoid exists.

The classifier now uses the floating-point calculation as a fast sign filter
away from the degeneracy tolerance, then an exact rational predicate on the
binary64 input coefficients near zero:

`sign(shifted level) = sign((4*c*det(A) - b^T adj(A) b) / det(A))`.

Thus the exact cone remains classified as `CONE` and the triangulator continues
to reject its singular apex, while the two nonzero sides are classified as
`TWO_SHEET_HYPERBOLOID` and `ONE_SHEET_HYPERBOLOID`. The same signed predicate
now distinguishes nonempty and empty positive-/negative-definite level sets
without the arbitrary level cutoff. The exact predicate interprets the stored
binary64 coefficients as the actual polynomial; it does not claim intent
recovery for rounded, externally transformed cone data.

## Verification

- The new near-cone assertion failed against the previous tolerance-only
  classifier. A follow-up test for a tiny nonempty ellipsoid also failed before
  the definite-quadric branch used the exact sign predicate; both passed after
  the fix. Empty ellipsoids remain rejected.
- Rebuilt the current project configuration (**131/131 build steps**).
- The full Faster algorithm CTest selection passed **44/44**; the only omitted
  test was the long five-seed unknown-answer sweep, replaced in this batch by
  two independent runs of the frozen seed-0 diagnostic.
- The focused geometry/source-motion subset also passed **6/6**:
  `spin-faster-quadric-coefficients`, `spin-faster-quadric-triangulator`,
  `spin-faster-source-motion`, `spin-faster-stratum-tangent`,
  `spin-faster-quadric-moments`, and
  `spin-faster-quadric-degenerate-boolean-witness`.
- The retained degenerate-Boolean witness still rejects its actual material
  boundary; it was not waived or smoothed.
- The fixed seed-0 solver diagnostic was rerun in two independent processes.
  Both remained `RESOURCE_LIMIT_REACHED` after 120 function evaluations and
  produced identical event hash
  `734ddf59e81d6ae440c54a698d35187875c48029ee81751ba559b3a2bf392044`, trace
  hash `aec4f6c5f96a6f0269373c631095a8bd3fe221a80e130e9de864755bb5541183`, and
  output hash
  `bfd46a26e0802ba98ef0b2240bd580c2c4921849dcecf928c8b1aed9be50427b`.
  No verified final geometry was produced. A later audit found that the
  root-isolation events were hex-encoded and had been missed by the first text
  filter; see the corrective follow-up in
  `docs/faster_root_isolation_pruning_2026-09-24.md`. The isolated-root branch
  sign has since received a separate exact Sturm–Tarski fix and focused
  regression, but this historical solver run predates that correction. Boolean
  degeneracy and the fixed solver budget remain unresolved.
- `git diff --check` was clean after implementation.

Project status remains `PARTIAL_BLOCKED`.
