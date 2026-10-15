# Faster truncation-root isolation pruning — 2026-09-24

## Scope

This is a strict-profile numerical-kernel repair and synthetic solver diagnostic.
It is not an author-input experiment or a completed Spin-It Faster reproduction.

## Confirmed defect and repair

The exact rational Sturm isolator applied its depth limit before pruning
subintervals that were already certified to contain no roots. A second ordering
issue applied the depth limit before the positive-branch interval tests could
exclude a one-root interval. Near-cone hyperboloid coefficients reached depth
65 after only 92–93 visited nodes; one interval had zero Sturm roots, and the
adjacent one-root interval was excluded by its exact interval bounds.

The isolator now:

1. rejects an invalid negative Sturm count;
2. prunes a zero-root interval;
3. applies exact positive-branch interval exclusions; and
4. enforces the unchanged depth limit only if further subdivision is still
   required.

The node cap (4096), depth cap (64), Sturm arithmetic, geometric tolerances, and
solver budgets are unchanged. The regression uses the captured near-cone
coefficients from the unknown-answer seed-0 trajectory and requires source
motion to be defined.

## Verification

- Before the implementation change, the new assertion failed at depth 65 with
  `sturm_roots=0`.
- After the first ordering correction, the same input reached a depth-65
  one-root interval; its positive-branch interval bounds proved it inadmissible.
- After moving the depth guard behind both certified pruning cases,
  `spin-faster-source-motion` and `spin-faster-stratum-tangent` passed (2/2).
- The frozen unknown-answer seed-0 solver diagnostic was run twice in separate
  processes with the same request, starts, two-restart limit, and 60-evaluation
  per-restart cap. Both runs produced the same event hash
  `abb2b55652f2ec6568028838b625737c3e1b5a83d1e07b4faee70b1dfd210ad9`, trace
  hash `aec4f6c5f96a6f0269373c631095a8bd3fe221a80e130e9de864755bb5541183`, and
  output hash
  `1d840e865922c3753a22ddc6ecd2b70aa14d502f9a047771f9ba0d2563bc98f1`.
- Neither run reported a root-isolation exhaustion. Both remained
  `RESOURCE_LIMIT_REACHED` after 120 total evaluations; no final material was
  verified. The trajectory still encounters invalid/degenerate Boolean
  geometry and reaches its fixed evaluation ceiling. These are unresolved
  solver/geometry issues, not infeasibility evidence.
- `git diff --check` was clean after the source and test changes.

The source-motion kernel regression passes, but this repair alone does not
establish full solver convergence or paper reproduction. Project status remains
`PARTIAL_BLOCKED`.

## Corrective follow-up — exact branch sign at an isolated root

The seed-0 diagnostic claim above that root-isolation exhaustion was absent was
premature: the event details were hex-encoded, while the first summary filter
searched only plain text. Decoding the native events showed a separate unresolved
case: a depth-65 interval with exactly one Sturm root whose `u*v` interval bound
still straddled zero. This was not one of the zero-root or interval-excludable
cases fixed above. The solver diagnostic therefore remains
`RESOURCE_LIMIT_REACHED`; the old event hashes describe that earlier build and
must not be used as evidence that root isolation was fully resolved.

The follow-up keeps the depth cap at 64 and adds an exact Sturm–Tarski query for
`sign((u*v)(root))` when the interval contains exactly one simple root and the
positive-branch scale is certified. The signed remainder chain determines the
algebraic sign without further bisection. The captured case now returns a defined
source-motion result (`defined=1`), and the regression requires that result rather
than merely checking that the depth-budget message disappeared. The geometric
strict-negative branch condition, root count, node budget, and all solver limits
are unchanged. This addresses the remaining isolated-root sign ambiguity only;
it does not establish that the full seed-0 solver converges or yields valid final
material geometry.

The follow-up regression was rebuilt in the exact-geometry configuration, then
the full configured build completed (115 relink/build steps). The Faster-labelled
CTest selection, excluding only the separately identified long unknown-seed
sweep, passed **44/44**. `git diff --check` passed after the final documentation
update. The five-seed campaign and the full solver diagnostic have not been
rerun for this narrow certificate change.

## Exact-query validation follow-up — 2026-09-24

Exposed the internal exact Sturm–Tarski sign-sum helper through a
`spin::detail` test header and added known-rational-root cases for positive,
negative, zero-at-a-root, and degree-greater-than-root-polynomial queries. An
initial expectation for `1-x^5` was incorrect because both roots lie in `(0,1)`;
the corrected negative case is `-x^5`. The production query implementation
passed the corrected cases without modification, so the hypothesis that this
helper requires a remainder reduction was rejected rather than shipped as an
unverified change.

The exact-geometry build rebuilt the source-motion regression, which passed,
and the Faster CTest label passed **44/44** with only the separately identified
long unknown-seed sweep excluded. This validates the current local algorithm
tests only; it does not change the seed-0 result, which remains
`RESOURCE_LIMIT_REACHED` with no verified final material. The tiny-face
Boolean audit failures and hyperboloid loop-construction failures remain the
next concrete solver/geometry blockers.
