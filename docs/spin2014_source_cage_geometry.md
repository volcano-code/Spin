# Actual-source represented geometry admission

`SourceCageGeometryChecker` owns immutable snapshots of the actual
`SourceCagePaperFunctional` and explicit settings. `check(fullState)` calls that
functional's `evaluateValues`; callers cannot supply geometry, evaluations or
pass flags. The returned `SourceCageGeometryAudit` has a private constructor and
const-only facts. Copying it preserves evidence; copying its facts does not make
another audit. There is no solver hookup in this component.

`CHECKS_PASSED` means all declared represented-polyhedral checks completed.
`REJECTED_GEOMETRY`, `RESOURCE_LIMIT`, and `UNKNOWN` never admit geometry. These
outcomes do not establish a feasible witness, convergence, stable spin, a sealed
artifact, the continuous BBW-map Jacobian throughout any cell, or original-paper
reproduction. Synthetic source fixtures are engineering regressions.

The world-space deformed points retain every original outer row and every
conforming material point. Original outer faces are unchanged. Each primitive
uses every stored oriented triangle reference and a compact sorted local index
map; IDs and beta owners remain in binding order. The existing pinned exact
boundary and signed-material-region implementation checks closedness, edge and
vertex manifoldness, nondegeneracy, orientation, exact self-intersection, and
alternating signed nested shells. Outer holes are not filled or reoriented. Each
primitive must have exactly one positive-volume component. Internal checker-only
Nef component orientation copies are the existing backend's representation
convention, and never alter source arrays.

Every primitive must be an inclusive subset of the actual outer material set.
Every unordered primitive pair is visited in lexicographic ID order. Strictly
disjoint binary64 AABBs can prove an impossible intersection; all remaining pairs
use the existing budget-aware exact Nef relation. Face, edge, and point sharing
are allowed; positive-volume overlap is rejected.

Exact outer/domain boundary separation uses finite nonnegative
`minimumOuterDomainClearance` in model units, default zero. Contact always
rejects, including when the threshold is zero. Each completed applicable
clearance record exposes the exact rational squared distance and witness face
indices. This is conservative clearance of the entire optimization domain from
the outer material boundary, independent of beta. It is not final binary-cavity
web thickness or all-feature shell thickness. Fractional beta does not realize
a half-cell cavity and is not binarized here.

## Deterministic evidence

The scientific SHA-256 binds raw binary64 full state (including beta, phi, and
all handles), actual source/frame/config identities, all settings, actual source
output, ordered attempted stages and their outcomes, resource plan and Nef
observations. Complete exact boundary, region, relation and clearance facts are
hashed, including incomplete outcomes and diagnostic vectors. Bounded records
retain primitive identities, beta owners, evidence and full-reason digests,
short diagnostic text and measured clearance. Timestamps are excluded. Changed
states never reuse a cached geometric pass. Invalid API dimensions/settings may
throw before an audit. Source failures such as invalid scales, nonfinite
arithmetic or invalid paper moments can precede geometry and retain UNKNOWN.

## Resource contract

Defaults are 100,000 points; 400,000 total outer-plus-primitive triangle
references; 4,096 primitives; 1,000,000 unordered pairs; 1,100,000 records;
512 MiB explicit checker accounting; 2,000,000 live explicit Nef complexity;
64 connected shells per mesh; 8,192 bytes per retained exact rational string.
Caps are explicit implementation defaults and do not recover an author protocol.
No cap is silently reduced. Addition/multiplication and native signed-index
narrowing are checked before size-dependent mesh, index or record allocations.

Preflight sums actual primitive reference counts without allocating a dense
primitive-pair table. For a local mesh with P upper-bound points and F faces,
accounting includes both its local mesh and the region's retained mesh copy,
index/DSU scratch, worst-case exact boundary diagnostic vectors including every
possible intersecting face pair, component-volume strings, and a conservative
F-by-F byte allowance for the existing signed-shell hierarchy. Primitive local
point upper bound is 3F before deduplication. All simultaneously retained region
meshes and known diagnostics are included. Record capacity is exactly
`3 + 3*primitiveCount + pairCount`; retained rational distances are separately
bounded. The source's immutable BBW and topology are shared; source evaluation
point/state/weight, beta-Laplacian and optional shape residual scratch is included.
Known diagnostic payloads receive a factor-four allowance for vector capacity
growth and simultaneous audit-to-region facts copies. Component dyadic volume
strings receive at least 8,192 bytes each in the transient accounting even when
the configured retained-string cap is smaller; the smaller cap still rejects
before retaining that region. Exact Nef complexity counts live explicit
application objects according to the existing backend budget definition.

The component-volume string allowance follows from the represented input bound:
every finite binary64 coordinate can be written as an integer of magnitude less
than 2^2098 over 2^1074. The existing auditor sums F origin-based tetrahedral
determinants, each with six cubic products, then divides by six. Here F is less
than 2^30 by the native-index cap. An unreduced numerator therefore needs fewer
than 6,328 bits and its common denominator fewer than 3,226 bits. Sign, slash and
decimal numerator/denominator together use fewer than 2,900 characters. The
8,192-byte allowance is conservative even before the factor-four copy/capacity
margin. This bound is for the auditor's signed volume strings; it is not a bound
on arbitrary intermediate Nef rational arithmetic. Clearance strings are
returned by the opaque exact-distance operation and checked before retention;
its temporary arithmetic/string storage remains within the external process
guard, not the claim of bounded retained evidence.

These are explicit application accounting limits, not an exact native allocator
or RSS proof. CGAL/GMP predicate and Nef temporary allocations, BVH internals,
allocator bookkeeping and rational arithmetic can allocate beyond the counted
payloads. Nef complexity is observed around existing operations; opaque native
peak memory inside a Boolean cannot be predicted. Run heavy checks sequentially
with OMP/BLAS/MKL threads set to one, `ulimit -v 25165824` and `timeout 3600`
(24 GiB virtual-memory and 60-minute process limits). A process killed by that
external guard has no completed audit. Resource exceptions and recognized
backend resource errors latch RESOURCE_LIMIT and stop further geometry. An
unknown or incomplete exact operation stays UNKNOWN; absent facts never count
as true. Sustained allocator failure may prevent allocating a returned digest
and propagate an exception; it cannot manufacture admission.

## Regression boundary

The registered regression uses the real source octree, conforming primitive
boundaries and all 26 actual BBW handles. The principal fixed source is the cube
[5/16,11/16]^3, design cube [.34,.66]^3, source depth 3, adaptive levels 2..3,
cage level 2 and FEM grid level 3. Its represented rest domain is [3/8,5/8]^3.
The independent analytic rest clearance is 1/16. Actual BBW summation rounds the
deformed binary64 coordinates even at identity: an independent rational
triangle-distance oracle gives minimum squared clearance
`19807040628565802923409276929/5070602400912917605986812821504`.
Its exact square root is the binary64 `0x1.fffffffffffc0p-5`; that threshold
passes and `0x1.fffffffffffc1p-5` rejects. No point is snapped to make 1/16 pass.

Coverage includes repeated scientific identity, all five objectives, changed
beta/phi and frame/config identities, actual global positive scale/translation,
a valid nonaffine motion, actual outer and primitive self-intersection, primitive
negative orientation, protrusion, exact cap boundaries, and a depth-4 signed
inner-void source retaining all 16 original rows and all 26 handles. A tiny
positive-scale state produces actual degenerate triangles rejected by the source
integrator before geometry; its audit is UNKNOWN at SourceEvaluation.

Explicit test-only `--probe-motion`, `--probe-collapse`, and `--probe-overlap`
arguments run bounded diagnostics; the no-argument registered test always runs
its assertions. The probes did not attain an exact outer-contact rejection or
a positive-volume pair overlap after all preceding checks passed. Those two
negative branches therefore lack an actual-source end-to-end fixture in this
task; the implementation composes the existing exact predicates and does not
claim those omitted cases were exercised. This limitation must remain visible
when considering subsequent solver integration.
