# Candidate complete feature-thickness contract v1 — review draft

**State:** `DRAFT_NOT_APPROVED`; this document does not alter a seal, a
manufacturing admission, or dynamics admission.

## Purpose and scope

This is a proposed engineering contract for the requirement
`spin2014.hollowing.minimum_wall`.  It is deliberately stronger than the
Spin-It fixed-shell constraint: the latter limits newly removed material, but
does not establish the thickness of old thin features or of material left
between boundaries.

The candidate applies only to a bounded, regular-closed, orientable,
self-intersection-free polyhedral material set `M` with a certified boundary
representation and an explicit *boundary-feature map*.  It is not a generic
triangle-mesh metric.  A mesh which lacks that map must return `UNKNOWN`, not
infer a map from face adjacency, triangle order, or a chosen sampling rate.

## Required boundary-feature map

The certified input is a tuple `(M, P, E, U)`, where:

- `P` is a finite collection of source-boundary patches.  Every boundary germ
  belongs to exactly one patch or to an explicitly declared exceptional set.
- `E` is an equivalence relation on local boundary germs which declares the
  two sides of a single exterior feature (including sharp creases) to be the
  *same local sheet*.  It is source geometry, not a triangle-index rule.
- `U` gives the physical unit conversion and the manufacturing threshold
  `tau > 0`, including the material/process authority for that threshold.

The map must survive vertex renumbering, triangle subdivision, face merging
and any exact re-triangulation of the same source patch.  A CAD B-rep may
provide it directly.  A tessellated source may provide it through stable
source-face IDs plus a certified patch correspondence.  Bare OBJ/STL input,
including the reconstructed G7 tea-pot, currently does not provide it.

## What counts as a material wall

Let `g_p` and `g_q` be boundary germs at `p` and `q`, with positive distance.
They are an **opposing-sheet witness** only if all of the following hold:

1. `g_p` and `g_q` are not equivalent under `E`.
2. Their source patches and material-side orientations are certified.
3. The open chord `(p,q)` lies in the interior of `M`; it must not travel
   through exterior air or a void.
4. The endpoint germ normals satisfy the declared opposing-angle bound.  The
   bound and its tolerance are part of `U`; it cannot be silently chosen by a
   checker.
5. The witness is classified as exactly one of: outer-to-cavity wall,
   inter-cavity web, or original thin feature.  A sharp exterior crease is an
   exceptional feature, not an opposing wall merely because nearby triangles
   approach at an acute angle.

This separates material sheets from triangle facets.  In particular, the two
sides of an acute convex exterior ridge can be equivalent germs even though
their triangulation contains two non-coplanar faces.  Conversely, opposite
germs on a thin neck can be distinct even when they live on one connected
outer boundary component.

## Complete all-surface certificate

For each admissible ordered pair of feature patches, the checker must solve
the global constrained infimum of `|p-q|` over its compact, stratified patch
domains and the five witness conditions above.  It must also certify coverage:

```
all non-exceptional boundary germs = union of checked patch strata
all possible opposing-sheet witnesses = union of checked pair domains
```

The certificate records each pair domain, whether it is empty, its global
lower bound, witness location when below threshold, and its numerical proof
method.  A pass requires every applicable lower bound to be at least `tau`.
It must separately report outer/cavity, cavity/cavity and original-feature
classes.  No finite point sample, ray sample, BVH nearest-pair query, sphere
topology test, or unlabelled non-incident-triangle scan can prove this
coverage.

For rational polyhedra the preferred proof uses exact predicates and a finite
stratified decomposition.  Otherwise, each reported bound `b` must include a
verified error interval `delta`; only `b - delta >= tau` is a pass.  Ambiguous
tangencies, non-manifold input, self-intersections, unclassified patch germs,
or an exhausted resource budget produce `UNKNOWN` (or `FAIL` when a certified
counterexample exists), never a positive certificate.

## Required adversarial validation before use

An implementation may be considered for independent review only after it has
fixtures with analytic expected outcomes for all of:

1. parallel thin plate;
2. narrow neck on one connected exterior component;
3. nested shell;
4. two cavities separated by a web;
5. curved shell with an inner cavity;
6. acute convex exterior wedge;
7. tangency, non-manifold, self-intersecting, and missing-feature-map inputs.

Every fixture must test source-patch re-triangulation invariance and the
boundary values immediately below, at, and above `tau`.  The wedge must not
fail merely because its local geometric feature size tends to zero; the thin
plate and neck must not pass because an implementation only inspected cavity
boundaries.

## Current consequence

No current sealed artifact supplies `(P, E, U)` with the above provenance and
coverage proof.  The existing exact fixed-shell, outer/cavity, and
inter-cavity checkers remain useful *partial* evidence, but do not implement
this contract.  Until this draft receives an independent technical review and
an implementation passes the adversarial suite, the G7 tea-pot remains:

```
same_component_feature_thickness = UNKNOWN
manufacturing_admission          = BLOCKED_UNKNOWN_WALL_THICKNESS
qualified PyBullet/OpenGL/video  = BLOCKED
```

