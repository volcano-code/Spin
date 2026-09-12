# Conforming source material boundary

`SourceConformingBoundary` consumes an actual immutable `SourceCageBinding`.
It builds closed, consistently oriented indexed surfaces for every included
Full cuboid, preserving the binding's primitive order and adaptive beta owners.
It requires both 2014 and exact geometry, independently of Faster.

This is an **IMPLEMENTATION_DEFAULT** surface discretization. The sparse
dyadic overlay and shared center fan are not recovered author triangulation
metadata. The overall strict reproduction remains **PARTIAL_BLOCKED**.

## Coordinates and public topology

Let D be `binding().source().config().maximumDepth` (at most 20).
`latticeLevel()` is D+1. `points()` starts with the byte-preserved original
outer rows. Its remaining rows correspond to lexicographically sorted unique
`materialPointKeys()` at lattice D+1. Coordinates use
`fma(ldexp(key,-(D+1)),rootSide,rootOrigin[axis])`. Original surface vertices
are never deduplicated with material vertices. Original outer triangles are
available unchanged through `binding().source().outerFaces()`.

All internal binding corner IDs are reindexed: double a binding corner's D-key,
binary-search `materialPointKeys()`, then add `outerPointCount()`. Every binding
corner is present. Do not reuse `binding().primitives()[i].pointIndices` to
index this boundary's points.

Each `Patch` records `axis`, `plane`, `level`, tile `coordinates`,
`negativeOwner`, `positiveOwner`, `triangleBegin`, and `triangleCount`.
The plane is in the **source D lattice**, so its output point-key coordinate
is `2*plane`. Tile coordinates are at `level`, with cyclic axes
u=(axis+1)%3 and v=(axis+2)%3. The output-key side length is
`1 << (D+1-level)`. Patch order is ascending (axis, plane), followed by
quadtree child order `uBit + 2*vBit` at each level.

`triangles()` stores each material triangle once as an array of three point
indices, oriented toward +axis. `primitiveBoundaries()[i]` lists
`OrientedTriangle{triangle,orientation}` for binding primitive i. An orientation
of +1 uses the stored triangle; -1 swaps its last two vertices. Upper faces
belong to `negativeOwner` and use +1; lower faces belong to `positiveOwner`
and use -1. An absent owner is `SourceConformingBoundary::noOwner`.
Shared patches reference identical triangles and points from both owners.

## Sparse construction and validation

Each primitive contributes six aligned dyadic squares. Each occupied plane has
a sparse quadtree, rooted at the root square. Face insertion takes expected
O(D * input faces) work with hashed plane roots; only face paths allocate nodes.
Traversal inherits owners and tiles covered missing quadrants without expanding
empty space. Same-side overlap is rejected. Plane roots are sorted before
traversal for deterministic output.

All patch corners and centers enter shared integer axis-line dictionaries.
Sorted range queries retrieve every point on each patch edge, including points
introduced by perpendicular faces. One center fan triangulates the resulting
boundary polygon. This avoids both independent coarse/fine face triangulations
and perpendicular edge T-junctions, without a finest-root-grid allocation or
primitive-pair/point-edge scanning. Point and line dictionary sorting adds the
usual logarithmic indexing cost. Every primitive must have a nonempty boundary
whose undirected edges occur exactly twice in opposite directions.

The represented D+1 spacing must strictly exceed the maximum endpoint ULP of
the root. Unsupported frames throw `invalid_argument`, with no snapping,
clamping, coarsening, perturbation, or fallback repair. Each rest triangle must
have distinct finite vertices and positive canonical area. Nonfinite arithmetic
throws `overflow_error`. An empty carveable binding produces no material
points, patches, or triangles and preserves its outer arrays.

## Resource and identity contract

All four configuration budgets are positive: `maximumFaceTreeNodes` (default
10,000,000), `maximumPatches` (2,000,000), `maximumPoints` (4,000,000, including
outer points), and `maximumTriangles` (16,000,000). Insertion is checked before
growth. Validation bounds Eigen/int indices, bytes, and derived counts of
two references per triangle, six edge entries per triangle, and three line
entries per point. Upstream source and BBW budgets retain their meanings.
These counts are not an exact combined resident-memory promise.
Malformed or overflowing budgets and allocation failures throw `length_error`
containing `RESOURCE_LIMIT`; no partially accepted object is returned.

Versioned SHA-256 identities bind the exact upstream binding output, all four
budgets, method and order conventions, lattice level, original geometry,
point keys/coordinate bytes, patch ownership and tiling, shared triangle indices,
and per-primitive oriented references. Hash iteration uses ordered arrays;
timing and caller-supplied pass flags are absent.

## Scope of the evidence

The native regression uses the actual outer cube (.3,.7), design cube
(.34,.66), source D5, eight adaptive leaves, 496 Full primitives, and a real
26-handle BBW field (coarse level 2/fine level 3). It independently enumerates
included unit faces and integrates signed tetrahedral moments, including all
ten moments and nonuniform beta totals. It applies an injective nonaffine map
and nonuniform all-handle BBW motion with fixed shared topology. These are
**SYNTHETIC_REGRESSION** checks, not original-paper experiments.

Indexed primitive closure and exact shared-triangle cancellation do not prove
global manifoldness at edge-only contacts, prevent deformed self-intersections,
or enforce shell thickness. Full source beta/phi/all-handle SLQP, paper
regularizers, analytic handle derivatives, geometric precommit rejection,
final exact material partition/seal, and original-paper assets remain required.
