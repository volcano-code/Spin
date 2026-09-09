# Original-surface boundary subtrees (Task 5c foundation)

`spin2014::SourceSurfaceOctree` is an immutable native geometry oracle. It owns
copies of the original outer and explicit design-region vertices and indexed
triangles. It accepts no occupancy reconstruction, optimizer state, or supplied
validation flags. The API and its linked implementation are available only when
both `MAKE_IT_STAND_BUILD_SPIN_IT_2014` and
`MAKE_IT_STAND_BUILD_EXACT_GEOMETRY_CHECKS` are enabled. Faster is independent.

This implements the geometry foundation of the boundary-cell subtrees described
in Spin-It (2014), Section 5.4. It does not implement adaptive optimization,
split/merge optimization history, the complete cage, or level-9 paper experiments.
Overall strict reproduction status remains `PARTIAL_BLOCKED`.

## Geometry and arithmetic contract

Both meshes must pass the existing exact indexed boundary and signed material
region audits, including manifoldness, orientation, self intersection, referenced
vertices and alternating component nesting. Valid disconnected components and
nested cavities are preserved. The original design material must be an exact
subset of the outer material. These tests use the existing Nef region operations
only as validation; they do not supply solver moments or a second moment checker.
No repair, reorientation of input, or Boolean fallback occurs.

`integrateSurfaceRawMoments` implements the supplement's Algorithm 1 polynomial
surface formula. Its ten moments integrate the unapproximated original polyhedron
in binary64 arithmetic. “Original surface moments” does not mean exact rational
moment arithmetic. The value-only and differential APIs share the same polynomial
accumulation, while only the differential API allocates and evaluates a dense
10-by-3N Jacobian. Both retain signed closed components and positive finite density.
Geometric legality is separate from the integrator's closure/finite/degeneracy
validation. Coordinates keep caller units: a degree-d moment has units
`density * length^(3+d)`. The octree fixes density to one for every moment.

The configured root has a finite origin, positive finite side, depth 0 through 20,
and a default maximum of 5,000,000 stored nodes. Both original boundaries must fit
the root; there is no clipping. Dyadic plane coordinates use
`fma(ldexp(integer, -level), side, origin)` so the same plane has identical
represented coordinates at every level. Numerical collapse, nonfinite arithmetic,
and zero/overflowed cell volume fail construction. Extremely small coordinate
spacing is conservatively rejected when the requested step is below the largest
binary64 spacing over the root. Moment sums remain subject to binary64 rounding.

## Closed boundary exclusion

The finest-cell policy is explicitly `IMPLEMENTATION_DEFAULT`, not a recovered
author sampling convention or exact cut-cell integration:

1. Each cell intersects original design triangles using pinned CGAL revision
   `f7c3c8212b56c0d6dae63787efc99093f4383415`,
   `Simple_cartesian<Gmpq>`, and `CGAL::do_intersect(Triangle_3, Bbox_3)`.
   This revision's closed-box SAT converts its represented box endpoints to the
   exact kernel before plane and axis projections. Touching counts as intersection.
   Inclusive source-coordinate bounding boxes and inherited intersecting triangle
   candidates accelerate the predicates without losing triangles enclosed by cells.
2. Only after proving no boundary intersection, the exact CGAL AABB-backed
   `Side_of_triangle_mesh` parity query determines actual material at the exact
   rational midpoint. Alternating signed nesting has already been validated.
   A full or empty terminal then represents the entire boundary-free cell.
3. Intersecting cells subdivide into all eight children until maximum depth.
   Finest intersecting cells contribute zero carveable moments and their whole
   cuboid volume as unresolved uncertainty.
4. Each internal node sums its actual children's ten moments and unresolved
   volumes. There are no dense grids, uniform upsampling, or synthetic descendants.

The accepted cuboids form a conservative subset of the design material. In exact
arithmetic, excluded design volume is no greater than the total unresolved cell
volume; reported binary64 quantities have rounding error. This is a volume bound,
not a relative-error certificate for all ten moments or a continuous infeasibility
proof. Cross and first moments can change sign or approach zero.

For root `[0,1]^3` and design `[1/4,3/4]^3`, touching excludes a layer on each side
of every boundary plane. At level L the accepted cube is
`[1/4 + 2^-L, 3/4 - 2^-L]^3` for L >= 3. At level 9 its side is `254/512`.
The original surface volume is `1/8`; the approximation has volume
`2048383/16777216`, and unresolved volume is `49153/8388608`.
These are analytic engineering fixtures, not original-paper cases.

## Queries, ownership, limits and identity

`cells()` exposes const stored nodes. Internal children occupy eight contiguous
indices starting at `firstChild`, ordered by x as bit 0, y as bit 1, z as bit 2.
Every stored node has its index, parent, level and integer coordinates. Stats
separate internal, unresolved boundary, full and empty nodes and expose actual
counts per visited level. A dyadic query walks this hierarchy. Descendants of a
coarser full/empty terminal return analytic cuboid/zero moments without allocating;
their index is `noIndex` and parent identifies the covering stored terminal.
Invalid levels or coordinates throw.

The node budget is checked before adding the root or any eight-child allocation.
Storage growth is capped by the budget. Exceeding it throws a `RESOURCE_LIMIT`
diagnostic; no partial oracle is returned. Allocation failure also propagates
instead of yielding valid empty geometry. Callers receive no mutable source,
configuration, tree or moment references; copying the object shares const state.

SHA-256 encoding uses explicit fixed-width big-endian integers and IEEE754 bits,
with row/column iteration rather than matrix storage bytes or object padding:

- Source identity binds both copied original V/F matrices and their dimensions.
- Configuration identity binds root, maximum depth, **maximum node guard**,
  versioned boundary policy, and pinned exact backend. The guard is an input
  control even though it limits resources.
- Partition identity binds both identities, original surface moments, and every
  actual cell's index, parent, children, level, coordinates, kind, ten moments and
  unresolved volume.

Timings, memory measurements and timestamps are excluded. Changing a node budget
changes config and partition identity even when the resulting tree is identical.
Identical inputs and configuration produce identical scientific hashes. Hashes
identify the calculation and do not grant feasibility, sealing, or reproduction
status. Existing lattice `AdaptiveVoxelDomain::source()` and solver/cage contracts
are unchanged.
