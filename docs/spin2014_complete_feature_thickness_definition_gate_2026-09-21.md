# Complete feature-thickness definition gate

## Decision

The current reconstructed G7 tea-pot remains
`same_component_feature_thickness=UNKNOWN` and
`manufacturing_admission=BLOCKED_UNKNOWN_WALL_THICKNESS`.  No new checker, seal
promotion, simulation admission, geometry repair, smoothing, remeshing, or
optimization rerun is authorized by this document.

This is a definition gate for the stronger engineering requirement in the
strict reproduction plan.  It is deliberately separate from Spin-It (2014)'s
fixed outer-shell constraint, which only limits where hollowing may occur.

## Facts already established

The existing sealed-material chain independently checks all of the following:

1. Exact fixed-region/design-region partition and preservation of the original
   outer surface.
2. Clearance from the original outer boundary to every newly carved boundary.
3. Exact separation between distinct bounded voids when there is more than one.
4. Reintegrated mass, center of mass, and full inertia from final geometry.

These prove neither a minimum thickness throughout a single connected material
component nor printability.  A nearest-triangle query can cross exterior air;
a finite surface sample has no all-surface coverage; and a distance from a
cavity to the outer boundary does not inspect pre-existing thin features.

## Rejected candidate definition

Do not use the infimum radius at which subtracting an open ball changes the
material topology as a wall-thickness measurement.  The exact acute-wedge
certificate on the unchanged tea-pot has a zero material-ball-cut infimum at a
fixed sharp exterior feature while remaining disjoint from every allowed
carving domain.  It is a valid topology observation, but it would incorrectly
declare ordinary sharp exterior detail to have zero printable wall thickness.

Likewise, do not adopt minimum normal-ray exit distance without a separate
sharp-feature convention: near a convex acute crease its infimum also tends to
zero although the observation concerns local boundary geometry, not two
opposing material sheets.

## Required definition before implementation

A future positive checker must publish a contract containing all of these
items before it can write `VERIFIED` into a `SealedSpinArtifact`:

1. A set-valued definition of which boundary sheets count as opposing material
   sheets, including treatment of faces meeting at an edge or vertex.
2. A representation-invariance rule: triangle subdivision, face ordering and
   canonical equivalent triangulations must not change the measurement.
3. An all-surface coverage argument or an explicit conservative unknown state;
   point sampling alone is insufficient.
4. Separate semantics for outer-to-cavity walls, inter-cavity webs, original
   thin parts, sharp exterior creases, tangencies and non-manifold or
   self-intersecting input.
5. A physical threshold unit and a declared relationship to the material and
   printing process.  A model-space clearance is not automatically a PLA
   manufacturing threshold.
6. Exact or bounded numerical error evidence, plus adversarial tests for a
   thin plate, thin neck, nested shell, multiple cavities, curved shell and
   acute convex wedge.

Until a contract satisfying all six items is independently reviewed, the only
valid outcome is `UNKNOWN`.  The existing fail-closed dynamics admission is
therefore correct: it may run analytic fixtures, but it must not claim a
qualified tea-pot simulation or physical validation.

## Consequence for the execution plan

Work may continue in parallel on source-bound optimization, mass moments,
derivatives, cage deformation, dual material, Faster geometry and strict asset
provenance.  The tea-pot's qualified PyBullet/OpenGL/video chain begins only
after this definition produces a real checker and a fresh seal passes it.
