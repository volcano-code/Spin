# Actual source cage material moments

`SourceCageMomentFunctional` is an immutable Task 5/M3 derivative layer. It
retains an actual `SourceConformingBoundary`, positive finite density, and an
explicit `SourceCageMomentFunctionalConfig`. It accepts no supplied field,
mesh, shell moments, or caller success flag. Both 2014 and exact geometry build
options are required; Faster is independent.

## State and result

`beta` has one entry per retained adaptive leaf, in domain order. Its entries
are hollowed fractions in `[0,1]`. Parameters have exactly `4H` entries in the
original cage handle order: `[scale, tx, ty, tz]`. Initial scales are one and
translations zero. Scales must be finite and strictly positive; translations
must be finite. This interface explicitly fixes all rotations to identity.
It does not claim implementation of the author's rotation variables.

For rest point `p` and actual fixed handle center `c_h`,

```
p' = sum_h w_h(p) [c_h + scale_h (p-c_h) + t_h]
dp'/dscale_h = w_h(p) (p-c_h)
dp'/dt_h = w_h(p) I
```

The weights interpolate the successful upstream **final** BBW field using
eight-node rest-space embeddings for every original/material point, including
face centers and hanging points. No moved point is re-embedded. No weights are
clamped, normalized, dropped, or substituted.

`evaluate(beta, parameters)` returns `SourceCageMomentEvaluation` with density
weighted `finalMoments`, `10 x B` beta and `10 x 4H` handle Jacobians, and moved
points in exactly the boundary's original-plus-material order. Raw moment row
order is `M,LX,LY,LZ,QXX,QYY,QZZ,QYZ,QXZ,QXY`.
`evaluateValues` runs the same values/geometry path with `hasJacobians=false`
and both Jacobians explicitly `0 x 0`. It requires no Jacobian output capacity.

## Oriented material calculation

The functional integrates the moved original outer faces and subtracts each
actual conforming primitive boundary weighted by its beta owner. A shared
triangle is stored once, in positive-axis orientation. Its negative-side
owner's outward orientation is positive; its positive-side owner's outward
orientation is negative. Therefore the coefficient **added to the outer
surface sum** is `beta[positiveOwner] - beta[negativeOwner]`, treating absent
owners as zero. Beta derivatives use the same positive/negative signs. Density
multiplies the accumulated result once. Equal beta on both sides cancels the
shared contribution exactly at the coefficient level.

The private `SurfaceTriangleMoments.h` kernel factors the existing supplement
Algorithms 1/2 polynomial values and local `10 x 9` coordinate derivatives.
The existing whole-surface public APIs still validate closed oriented chains;
the private per-triangle kernel makes no closure claim. Handle derivatives
contract each triangle's local derivative with all fixed-rest handle weights.
Production does not use finite differences or signed tetrahedral formulas.

## Checked resources and identities

Defaults are `maximumWeightValues=1000000` and
`maximumJacobianValues=10000000`. Budgets must be positive and representable in
both Eigen indices and bytes. The Jacobian cap bounds the **combined** number
of output values `10B + 40H`. Weight scratch is exactly one `H`-entry row;
the cap must permit every actual handle and also obey the upstream binding and
grid interpolation budgets. A budget cannot select fewer handles. Retained
geometry uses the upstream boundary point cap; cached embeddings are bounded
by actual point count, container capacity, and checked byte arithmetic.

Storage is `O(P)` embeddings/points, `O(H)` weight scratch, and, only for
derivative evaluation, `O(10B+40H)` output derivatives. There is no global
point-by-handle or point-by-parameter matrix. Invalid/overflowing budgets and
allocation failure throw `length_error` containing `RESOURCE_LIMIT` without
partial output. Malformed state/density and degenerate triangles throw
`invalid_argument`. Nonfinite computed arithmetic throws `overflow_error`.
Every triangle is checked even at zero coefficient. Nonpositive final mass is
rejected. Empty carveable material is valid and contributes no cavity terms.

Version-one SHA-256 identities encode IEEE-754 doubles and fixed-width integers
in deterministic order. Configuration identity includes version, method,
layout, identity-rotation convention, density, all budgets, counts and upstream
boundary configuration. Output identity also includes actual boundary and BBW
output identities, both input vectors, evaluation mode and every produced
scientific array (row-major digest traversal). Timing is excluded. Value-only
and derivative modes have different output hashes; their values and geometry
are identical. Distinct input states remain distinguished even if their
motions cancel. Returned results can be modified without mutating retained
input snapshots.

Positive mass and finite/nondegenerate triangles do **not** establish geometric
legality. Self-intersection, inversion relative to rest, material overlap,
minimum shell thickness and final exact partition remain separate precommit
gates. This layer provides neither the full joint solver nor the paper
regularizer, SLQP, physical certificate, original-paper case acceptance, or
paper completion. Global status remains `PARTIAL_BLOCKED`.
