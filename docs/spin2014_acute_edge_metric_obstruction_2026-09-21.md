# Acute edges: a material ball-cut obstruction, not a wall measurement

## Why this check changes the next decision

The G7 optimizer, original rounding and independent mass checks have passed for
two identical reconstructed runs. Complete same-component feature-thickness
admission remains open. Earlier vertex/nearest-chord diagnostics cannot supply
all-surface coverage. Before implementing a complete epsilon-map as a possible
extra checker, it is necessary to ask whether its criterion can admit this
unchanged input at any positive threshold.

The reference [Epsilon-shapes, Sections 3 and 5](https://arxiv.org/html/1704.08049)
describes arbitrarily small topology-changing modifications near acute
configurations. Its formal signed-map definitions and later material-side
terminology conflict. This implementation therefore names the operator
explicitly: **material minus an open ball**, never ball union. It does not
claim to implement the full reference algorithm or Spin-It's wall rule.

## Exact continuum certificate

`experiments/teapot_spin2014/acute_wedge_query.py` uses the original represented
binary64 coordinates as rationals. All edge predicates are exact, without a
floating-point candidate filter. The supplied mesh must separately pass an
embedded, outward-material audit; the local query's link checks alone cannot
establish that prerequisite.

For two incident faces at a convex acute edge, let `m` be its midpoint and `u`
the inward direction in the first face perpendicular to the edge. Let `w` be
the projection of `u` onto the second face. Define

```
x(t) = m + t u
q(t) = m + t w
a = |u-w|²
b = |u|²
r(t)² = t²(a+b)/2
```

The checker verifies `0 < a < (a+b)/2 < b`, strict interior containment of
`x(t_max)` and `q(t_max)` in their triangles, the negative material halfspace
conditions, and `4 t_max² b < d²`. Here `d²` is the exact minimum squared
distance from `m` to **every other triangle**, not a sample or bounding-box
estimate. It additionally ensures `r(t_max)²` is below the explicit diagnostic
threshold squared. An exhausted halving budget returns `RESOURCE_LIMIT`, not a
certificate.

These conditions hold throughout `0 < t <= t_max`. The whole ball lies closer
to `m` than `2t|u|`, so it cannot encounter another triangle. Its radius is
smaller than the distance to the shared edge but larger than the distance to
the opposite face. It therefore intersects the two faces in disjoint disks.
The convex material between them meets the sphere in an annulus; ball removal
changes the material Euler characteristic by -1. Such radii approach zero.

This is a continuum counterexample to a positive **uniform ball-cut-neutrality
radius**, not a claim that a physical wall has zero measured thickness. It
does not assign a thickness value at the edge, does not exempt sharp parts,
does not authorize smoothing, and does not change the original outer-wall
clearance contract. Release wall thickness remains `UNKNOWN`.

## Regression and review evidence

Six tests in `tests/python/test_acute_wedge_query.py` pass. They exercise an
acute prism, ordinary right-angle/coplanar edges, reversed cavity orientation,
budget exhaustion, invalid requests and exact plane-preserving subdivision.
At four rational scales, independent whole-boundary cap filtration and segment
parity confirm two disk caps and an interior material chord. Those samples
cross-check the certificate; the inequalities above provide its all-scales
coverage.

Independent read-only review accepted the geometric argument and requested
the explicit operator naming. The updated state/operator tests failed on the
old naming before the implementation was corrected and all six passed again.

## G7 source-bound experiment

The reviewed experiment package is
`build/spin2014-acute-wedge-audit-20260921-reviewed/`. It uses the unchanged
latest G7 seal (manifest SHA-256
`7aa3f73d13fd9157322ddb79a18aeefb7a5fb8ba91f17b6b6f5f304fad6d4bc7`).
It reauthenticates original, final and allowed design-region geometry, scans
all edges, chooses the earliest indexed original convex-acute edge, constructs
its certificate separately on original and final geometry, and repeats each
twice. Independent fixed-center material-cut checks supplement the continuum
proof. An exact family-ball bounding cube is compared with all design-domain
bounding boxes: separation proves the whole certified family is outside those
domains; lack of separation must not be interpreted as a proof of overlap or
immutability.

The earlier `build/spin2014-acute-wedge-audit-20260921/` attempt was deliberately
stopped after review found the ambiguous operator name. Its files and the
pre-review checker/test source copies are retained. The explicitly stopped
unit was confirmed `not-found/inactive`, with MainPID 0. The controller's
second cleanup reported `JOB_CLEANUP_FAILED` after the unit had already been
removed; there is no successful run receipt and no scientific acceptance for
that attempt. This is not an optimizer failure or a recovered complete run.

Final counts, resources and test results must come from the reviewed job's
receipts. Until then, this section describes the experiment, not its outcome.

## Completed receipts (addendum; earlier attempt retained)

Both the reviewed earliest-index run and the separately registered
`build/spin2014-acute-wedge-audit-20260921-strongest/` run have now exited with
code 0. The latter took 56.570210303 s, peaked at 353,701,888 bytes, used no
swap, and executed all 13 geometry-diagnostic CTest groups successfully.
Each checks original and final geometry twice; the repeat observations match.

The exact scan finds 18 convex acute edges in each mesh. Seventeen occur at
the shaft cap and differ from right angles by less than 6e-16 degrees in the
represented binary64 geometry. These are not evidence of a substantial
physical thin feature. The separate strongest-angle selection instead chooses
spout edge [3854, 3858], with an interior dihedral of approximately 84.98077
degrees. Its exact moving-ball family is disjoint from every allowed design
domain, in both the unchanged original and final material. The supplemental
selection is `MAXIMUM_EXACT_NORMAL_COSINE_SQUARED_THEN_EDGE_ID`; it does not
replace the preserved earliest-index result.

The strongest-run observation digests are
`86b695a8f86bd3ddf953a9829eb5af931e66c559f561db62b1f9f182a5197672`
(original) and
`081c616551fdfd7015036759d4ede08c1cfb114fbdb283c45b7a6258ab7aca9a`
(final). The authoritative report is `job/payload/report.json` within that
batch. It explicitly retains `complete_wall_thickness=UNKNOWN`,
`epsilon_metric_adopted_for_release=false`, and paper/physical admission
`NOT_GRANTED`.

Consequently repeated hollowing cannot remove this fixed exterior obstruction
to the supplementary uniform ball-cut criterion. This does not define or
measure the printing wall thickness. Do not change the original wall contract,
smooth the mesh, or count a further ball-cut diagnostic as complete thickness
validation. A physically appropriate full-feature thickness contract and its
checker remain unresolved; independent simulator verification on analytic
fixtures can proceed without admitting this G7 material to physical tests.
