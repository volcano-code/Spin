# Actual deformed binary material validation

`validateSourceCageBinaryMaterial` consumes only a privately constructed
`SourceCageAdaptiveSolverAudit`. Admission requires `RelaxedAdaptiveComplete`,
an owned `lastSolved()`, no pending checkpoint, and the immutable terminal
solve/source/state/partition identity chain. A pending, resource-limited,
nonconverged or otherwise incomplete adaptive run returns `NOT_ADMITTED`.
There is no copied-facts, arbitrary-state or checkpoint overload.

This is a bounded **pre-seal** object. It is not `SealedSpinArtifact`. It does
not establish the exact final region partition, a continuous BBW Jacobian,
full cavity-web thickness, a feasible witness, stable spin, physical behavior,
or paper completion. Exact region partition/thickness and Python sealing
remain subsequent distinct gates.

## Rounding and emitted surfaces

Only beta changes: beta >= 0.5 becomes hollowed=1; otherwise it becomes 0.
The exact-half decision is an explicit implementation convention matching
the existing source rounding, not a recovered author tie rule. Every angle
and handle parameter is copied bit-for-bit. There is no reoptimization,
repair or pose recentering. The relaxed and binary full states and their
delta are retained. Because adaptation may have merged leaves, these vectors
use the terminal endpoint's layout, not necessarily the original layout.

The binary functional's `evaluateValues(binary)` provides the actual deformed
conforming points. The outer mesh retains the original triangle connectivity
and these actual deformed outer points. For every selected primitive, its
oriented triangle references are added; opposite internal references cancel.
Surviving triangle order follows the immutable conforming triangle order.
Referenced vertices are compacted by ascending original point index. No rest
cuboid or undeformed outer geometry is used for emission.

`carvedBoundary()` has the positive cavity orientation. `materialBoundary()`
appends its reversed faces to the outer mesh. An empty cavity has zero rows
and three columns in its vertex and face matrices, and `emptyCavity=true`.
The nonempty outer material is still independently integrated and exactly
audited; no empty-mesh audit is fabricated. Edge/corner-contact unions can
produce a nonmanifold cavity boundary and fail. No mesh repair is attempted.

## Independent numerical and geometry gates

The final emitted material triangles are integrated by signed tetrahedra
`(0,a,b,c)` using barycentric first/second moments with long-double
accumulation. This implementation does not call the source functional's
surface integrator. All ten density-weighted world moments are compared
with the actual binary functional, never with the relaxed functional.

For a degree-d moment, the residual is divided by `density * rootSide^(3+d)`.
Density and root side are fixed from the immutable INPUT source. The gate is
`max(abs(normalizedMomentResiduals)) <= 1e-10`. It is never normalized with a
deformed bounding box, moving center or output-derived length. Arithmetic
must remain finite and the mass must be strictly positive.

The independently integrated moments are transformed by explicit first- and
second-moment identities into the endpoint's fixed paper frame. The raw
Givens constraints are recomputed using the retained phi. First-moment
constraints are divided by `density * rootSide^4`; second-moment constraints
are divided by `density * rootSide^5`. Both `rawConstraints` and
`normalizedConstraints` are retained. The normalized infinity norm must not
exceed `constraintTolerance`, a finite positive `IMPLEMENTATION_DEFAULT` of
1e-6. `centralInertia` and its eigenvalues are recomputed in these fixed paper
coordinates. The minimum eigenvalue divided by `density * rootSide^5` must
be at least `-128 * epsilon(binary64)`, the explicit implementation allowance
for numerical PSD error matching the earlier source binary moment check.
A constraint failure remains `BINARY_VALIDATION_FAILED` and retains raw
residuals, including cases where small density would mask the raw violation.

The actual `SourceCageGeometryChecker` runs again on the binary state with
the explicit saved geometry settings. The final emitted material mesh is
then passed to `ExactBoundaryAudit` with `AnyNonzero`, which checks indexed
closure, winding, edge and vertex manifoldness and intersections. Both
audits run even when rounding has already failed the constraint gate.
`AnyNonzero` alone does not establish nesting or a final exact region partition.

## Ownership, diagnostics and bounded resources

The audit and material constructors are private; accessors return const data
or shared pointers to const owned data. References survive moves while an
owner remains; moved-from access throws. `candidate()` exposes raw mesh and
state evidence after emission even on subsequent validation failure. It is
not a validation grant. Only `validatedMaterial()` publishes bounded validated
material, and it is null on every failed, unknown or resource-limited result.
Copied scalar facts cannot change either accessor's result. Candidate data
retains the immutable functional after the originating audit is destroyed.

Evidence includes adaptive-input, functional, partition, configuration,
relaxed/binary-state, binary-functional-output, emitted-mesh-output and full
scientific digests. The scientific digest includes every raw numerical gate,
geometry audit, final boundary audit, status, resource reservation and reason.
No elapsed time is included. Repeated evaluation is deterministic.

Mesh counts, Eigen/int index capacities and checked size products are bounded
before size-dependent allocation. `maximumOwnedBytes` covers known state,
evaluation, mesh and assembly scratch buffers. `maximumExactDiagnosticBytes`
conservatively reserves four diagnostic payloads, including all possible
face pairs and 8192 bytes per possible component volume string, before the
final exact audit. Geometry checker caps are separate explicit settings.
These bounds do not promise native CGAL allocator bytes or RSS. Run under
the external 24 GiB virtual-memory and 60 minute process guards, with
OMP/BLAS/MKL threads set to one.

The regression uses actual source geometry, full BBW and guarded adaptive
solves. Its stationary fixtures are labeled `KNOWN_SOLUTION_VALIDATION`:
zero objective weights and a fixed principal frame make a known stationary
endpoint. They demonstrate validation behavior, not an optimization or
paper achievement.
