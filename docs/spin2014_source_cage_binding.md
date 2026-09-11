# Source material / full BBW rest binding

`SourceCageBinding` is the immutable rest-geometry connection between a
`SourceAdaptiveDomain` and an actual successful `SourceCageBbw`. It is available
only with both Spin-It 2014 and exact geometry enabled, independently of Faster.
Global strict reproduction status remains **PARTIAL_BLOCKED**.

## Inputs and immutable identity

The constructor takes a const adaptive domain, a real BBW run and an explicit
`SourceCageBindingConfig`. It has no supplied-weight, success-flag or supplied
moment constructor. BBW must report success and an available, passing final
primal assessment. Both BBW protocols are accepted; independent normalization
retains its raw and final assessments and is not promoted to coupled KKT.

The original outer vertex and triangle arrays must match the cage's input
arrays in dimensions and scalar bytes. Root origin and root side must also
match in scalar bytes, including signed zero. The comparison does not compare
source digests: the material oracle has a separate design surface, while the
cage's source identity uses outer-as-design. Both upstream identities remain
available through the immutable objects.

The binding copies the adaptive partition, including original leaf order,
moments, adjacency, Laplacian and partition identity. Later caller split/merge
operations cannot change the binding's domain. The geometry oracle and BBW
retain their immutable shared upstream data. Accessors expose const references;
copies of the binding share private const Data.

## Actual material and point layout

The implementation calls `sourceBinaryCuboids(snapshot, ones, maximumPrimitives)`.
It retains actual included Full cuboids, including descendants of boundary
adaptive leaves. An adaptive leaf's bounding cube is never substituted for
these source subcells.

Each `Primitive` has its original `sourceKey`, canonical `bounds`, one
`betaIndex` in the copied domain's leaf order, and eight `pointIndices` with x
bit first. Owner lookup probes the primitive's ancestors in a hash map, giving
expected O(depth) lookups rather than scanning all leaves for every primitive.
Every primitive must have exactly one owner and every retained leaf must own
material. Primitive order follows the upstream sorted source keys.

`points()` has three columns and this layout:

1. The original outer vertices, with their original order and bytes.
2. Unique material corners in lexicographic integer `(x,y,z)` order on the
   source `maximumDepth` lattice.

`materialCornerKeys()[i]` identifies point `outerPointCount()+i`. Integer keys
deduplicate only material corners. Even a coincident original outer vertex
retains a separate point index. Coordinates come from actual canonical source
cuboid endpoints, with byte-consistency checks for coincident keys. They are
never snapped to the coarser weight grid. Original triangles remain available
as `source().outerFaces()` and index the unchanged initial point rows.

## Bounded full-handle transfer

`embeddings()` contains one actual `grid.embed(point)` result per point: the
selected cell, its integer key, eight node indices and eight coefficients.
No global dense point-by-handle matrix or `3N x 4H` LBS matrix is allocated.

`pointWeights(begin,count)` returns `count x handleCount()` values by summing
the cached eight node contributions from the actual final BBW field. It keeps
all handles in upstream order, without clamping, normalization or subsampling.
`begin==points().rows(), count==0` is valid and returns zero rows with all handle
columns; `begin>points().rows()` is invalid even when count is zero.

The three budgets default to 1,000,000 primitives, 2,000,000 points and 1,000,000
returned weight values per chunk. All must be positive and representable by
their relevant containers, Eigen indices and byte allocations. These are count
bounds rather than a combined process RSS guarantee. Existing upstream objects
also retain their own storage and resource contracts. Corner insertion stops
at the point cap; it does not first assemble an unbounded corner list. Actual
allocation exhaustion is translated to `length_error` with `RESOURCE_LIMIT`.

Malformed/overflow budgets and exceeded primitive, point or output budgets
throw `length_error` with `RESOURCE_LIMIT`. Invalid or mismatched inputs throw
`invalid_argument`; invalid chunk ranges, including overflowing requested
counts, throw `out_of_range`. Nonfinite interpolation arithmetic throws
`overflow_error`. Construction either returns a complete object or throws:
there is no partial binding, silent coarsening or dense root-lattice expansion.

## Scientific hashes

`configSha256()` uses the versioned `source-cage-binding-config-v1` convention
and binds upstream source/BBW configurations and every binding budget.
`outputSha256()` uses `source-cage-binding-output-v1` and additionally binds:

- Source content and source/adaptive partition identities, plus cage source
  and actual BBW output identities.
- Original outer arrays, root bytes and the complete final nodal weight field.
- All retained leaves in order, their keys and ten moment entries.
- Original point count, point coordinates, material integer keys and handle count.
- Every primitive's source key, canonical bounds, beta owner and connectivity.
- Every cached embedding's cell identity, node indices and coefficients.

Integers use fixed-width big-endian encoding; binary64 values use their IEEE
bits, and strings/matrices have explicit lengths/dimensions. Timings and
diagnostic text are excluded. No unrelated Faster data enters these hashes.

## Regression scope and remaining gates

The native synthetic regression uses outer cube `[.3,.7]^3`, design cube
`[.34,.66]^3`, root `[0,1]^3`, source depth 5, optimization initial level 2 /
maximum level 5, cage level 2, fine grid level 3, actual FEM and all 26 handles.
The immutable material contains 8 adaptive leaves, 496 Full primitives and 997
points. Its volume is `(10/32)^3 = 0.030517578125`, compared with `0.125` from
the eight adaptive bounding cells. It validates independent cuboid ten-moment
integrals with nonuniform owner beta, canonical points, manual eight-node
interpolation, repeat identities, split/merge snapshots, both BBW protocols,
resource limits, byte mismatches including signed zero, and transformed input.
A 182-value chunk budget accepts seven complete 26-handle rows but rejects the
25,922-value global point field.

This binding evaluates no deformed material or objective. Under nonaffine LBS,
independently triangulated cuboids can have hanging/shared-face conformity
problems. A conforming deformed material boundary, geometry validation, full
joint SLQP, the paper regularizer and original-case acceptance remain mandatory
future work. The binding does not seal a material artifact, grant feasibility,
recover author settings or remove source approximation uncertainty. Synthetic
examples are `SYNTHETIC_REGRESSION`, never paper reproduction evidence.
