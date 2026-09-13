# Actual final boundary-component separation

`validateSourceCageComponentClearance(regionsAudit, settings)` accepts only the
immutable validated Task5z region product. The new immutable audit exposes
`facts()` and optional `validatedClearance()`; the latter retains the entire
original `regionsAudit()` and `regions()`. That chain preserves the original
binary material arrays, states, functional and world/spin moment evidence after
all originating audit variables have been destroyed. Copied facts and mutable
meshes cannot construct an audit or product.

This gate establishes **BOUNDARY_COMPONENT_SEPARATION** only. Its finite
nonnegative `minimumComponentSeparation` is in model units and is converted by
the existing exact backend from binary64 to rational without a tolerance. It
measures gaps between every two distinct connected components of the actual
final material boundary, including disconnected outer-root gaps through air.
Contact always rejects, even at threshold zero. Threshold equality passes.
`fullManufacturingThickness` always remains `NOT_ESTABLISHED`.

## Geometry and evidence

The input boundary hash is recomputed from the actual final mesh using the
backend's `exact-boundary-input-v1` / `AnyNonzero` protocol and matched against
both Task5z exact material and boundary input digests. Exact upstream
`faceComponents` supplies canonical component IDs. Faces retain their original
global order; each component's vertices are compacted in ascending global ID.
Coverage is independently checked: every original face and referenced vertex
belongs to exactly one extracted component, no component is empty, and no
vertex is shared across different components.

Each component retains source vertex/face IDs, original sign and parent/depth,
and original/checker input digests. Negative inward components are reversed
only in a local checker mesh by swapping triangle indices 1 and 2. The original
material stays byte-for-byte unchanged. Each whole extracted component—not
one selected voxel primitive—is validated as a standalone bounded positive
exact region. Its complete region and boundary facts are retained.

Every unordered pair `i < j` is measured in lexicographic order. Each record
retains complete exact clearance facts, the rational squared distance and exact
threshold, local witness face IDs and their mapped original material face IDs.
A measured rejection does not stop subsequent pair measurements. A multi-
component completed result requires exactly `n*(n-1)/2` complete measurements;
unknown or partial evidence cannot publish. The global minimum is compared
using exact rationals, with the earliest component pair winning ties. The
backend independently chooses its canonical local-face witness.

A single boundary component has explicit `NOT_APPLICABLE`, zero pairs and
absent distance, witness, contact and threshold-result optionals. This is
distinct from `INCOMPLETE` caused by unavailable measurements or resource
failure. Multiple components use aggregate state `MEASURED`; individual native
clearance records keep the backend state `COMPLETE`.

## Budgets and identities

Settings bound components, source points, source triangles, unordered pairs,
known owned bytes, retained exact diagnostics and shared native Nef complexity.
Zero, oversized or overflowing caps reject before variable allocation. Count
caps fit signed native indices; byte caps fit `ptrdiff_t`. Pair counts use
checked multiplication. Preflight covers the input face-component array,
component maps, source coverage and local-index scratch, mesh/native-wrapper
copies, retained pair containers and conservative edge/triangle-pair diagnostic
storage. Actual variable diagnostic strings and arrays, including rational
distance strings, are cumulatively checked before retained evidence copies.

All standalone component regions share one local `ExactNefBudget`. Live exact
regions and the mutable budget never escape; the audit receives a copied final
snapshot after the local operands release their charges. CGAL/GMP allocator
internals and rational intermediate allocations remain subject to the caller's
external 24 GiB / 60-minute guard. These counters are not an RSS guarantee.

Scientific SHA-256 binds every setting, upstream audit identity, actual material
identity, maps, reversal, full native region/boundary and pair facts, aggregate
diagnostics, outcome/reason, publication presence and copied budget statistics.
The hash contains no timestamp. `NotAdmitted`, measured `Rejected`,
`ResourceLimit` and `Unknown` are distinct and never publish a validated product.

## Limits

No separation between features of the same connected boundary, medial-axis
thickness, full wall/web manufacturing thickness, continuous BBW Jacobian,
author continuous design, sealing, feasibility authority, physics or paper
convergence is established. No deformed `SealedSpinArtifact` or Python override
is provided. The native real-pipeline regression uses fixed synthetic
`KNOWN_SOLUTION_VALIDATION` fixtures with zero objective weights, full handles
and nonzero phi; it is not an original-paper experiment.
