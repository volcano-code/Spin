# R8 final manufacturing-geometry acceptance record

Specification audit date: 2026-08-12 (Asia/Shanghai)

Acceptance status: **accepted at the documented R8 unit-contract ceiling**

## Scope and claim ceiling

R8 implements the final computational geometry and export layer requested by
Phase 8. It adds three deliberately separate contracts:

1. a fail-closed pre-slicer geometry gate and deterministic OBJ/STL export;
2. an opt-in, explicitly pinned local Prusa-style slicer smoke-test adapter;
3. an audit for an externally authored print-part declaration.

Acceptance means that these interfaces and their constructed fixtures pass.
It does **not** mean that an official paper case was reproduced, a real slicer
was run, a model was printed, a structure carries load, or a physical object
stands. Those claims require R9 artifacts and process-specific evidence.

## Final coherent boundary

`ManufacturingGeometry` compacts synchronized `M_O` and `M_I` into one
boundary. A ready report requires:

- finite vertices, valid indices, no degenerate or duplicate faces;
- closed edge- and vertex-manifold triangle components with consistent edge
  directions;
- exactly one outward outer component and zero or more inward cavity
  components;
- deterministic BVH broad phase followed by triangle-level
  self-intersection rejection;
- every cavity inside the outer shell, with no nested cavity components;
- exact closest triangle-feature clearance at least the requested model-space
  thickness;
- positive sampled Jacobians and nondegenerate boundary triangles for every
  active deformed hexahedral cell;
- exactly one six-neighbor material component; and
- positive signed material volume and agreement between surface and independent
  deformed-cell `(m,Q,c)` calculations.

The 27 Jacobian samples per trilinear cell are a conservative rejection
heuristic, not a continuous injectivity proof between samples. Minimum
clearance is geometric and expressed in model units. STL carries no physical
unit metadata; R8 does not convert the result into millimetres or validate a
printer's minimum feature capability.

The overload accepting an `AlternatingBalanceOptimizationResult` additionally
checks terminal status, reached targets, target generation, support mode,
occupancy, geometry epoch and content fingerprint, parameter state, mass/CoM,
and terminal trace consistency. It rejects a positive support band that only
previews hypothetical flattening. Because the result structure is publicly
constructible, `sourceBound` denotes supplied-result internal consistency—not
authenticity or proof of optimizer origin. The private export seal detects
post-preparation mutation; it is not a cryptographic attestation.

R8 performs no arbitrary mesh repair. Unreferenced-vertex compaction and
deterministic serialization do not change topology. A failed geometry gate
cannot be bypassed by calling the guarded exporter.

## Deterministic export

Ready geometry can be atomically exported as OBJ, ASCII STL or binary STL.
The exporter uses the classic locale, deterministic face order and stable
numeric formatting. Binary STL coordinates/normals are revalidated when
converted to single precision. Round-trip tests compare triangle geometry,
orientation, byte-stable repeated output and failure transactions. No partial
destination replaces a pre-existing file after an export error.

## Pinned local slicer adapter

`validateWithPinnedSlicer` is opt-in and fail-closed. It requires absolute
executable/profile paths, exact SHA-256 pins, exact full version output and a
positive timeout. Commands are passed as argv to `execv`; PATH and a shell are
never used. The output is staged in a mode-0700 private directory. The direct
child and same-process-group descendants are signalled on timeout or leader
exit; diagnostic streams and input/output artifacts are hashed, and the output
is atomically committed only after its hash remains stable. Known auto-repair
diagnostics are rejected.

The G-code audit is intentionally a smoke test: it requires text, a recognized
layer marker, and at least one finite positive-extrusion G0/G1 command with X
or Y displacement while tracking absolute/relative coordinate and extrusion
modes plus G92 resets.
It is not a full G-code parser, collision simulation, printer compatibility
check or build-success proof.

All slicer tests use the repository's fake executable. No real PrusaSlicer
binary, exported profile or generated process artifact is pinned in this
repository. Consequently, the adapter is `VERIFIED_UNIT`; actual slicer
compatibility is `NOT_RUN`. The child inherits the host environment and
dynamic-library closure, so this is not a hermetic or hostile-input sandbox.
The caller must retain normal, exclusive child ownership during the call: it
must not set `SIGCHLD` to ignored/`SA_NOCLDWAIT` or concurrently reap/fork in a
way that races this adapter. The fork/exec call must also be made while the
process is single-threaded. A descendant that deliberately leaves the process
group is outside the lifecycle guarantee.

## External partition declaration

`ExternalPartitionManifest` records a source-boundary hash, a declarative
model-unit-to-millimetre scale, external generator identity, per-part mesh and
evidence-file hashes, and rigid assembly transforms. The audit requires an
absolute manifest root, lexical-relative references whose resolved regular
files remain below that root, unique part identifiers, matching hashes and
finite proper rigid transforms. A canonical-field-stream SHA-256 fingerprint
binds the aggregate audit to the supplied declaration.

This is a local evidence audit over a trusted, non-mutating filesystem
snapshot. The path containment check is not an adversarial filesystem sandbox:
files are reopened by path for hashing and a concurrent actor could replace a
path after resolution. The fingerprint is a content digest, not a signature,
origin attestation or guarantee that files remain unchanged after the audit.

This interface does not generate a partition or interpret the referenced
geometry/slicer reports. It does not prove source coverage, non-overlap,
boolean reconstruction, cut quality, connector design, support/drainage
design, assembly path, tolerance, structural strength or printability. Those
evidence files are treated as opaque, hash-bound external artifacts.

## Registered R8 test matrix

| Group | Count | Covered contracts |
|---|---:|---|
| Geometry/export | 26 | Boundary topology/orientation/self-intersection; outer/cavity coherence; containment; wall threshold; deformed cell inversion/collapse; one material component; result binding/rejection; deterministic transactional OBJ/STL export; seal and scale behavior. |
| Slicer adapter | 18 | Not-configured state; executable/profile pins; full-version mismatch; stable invocation fingerprint; success; repair rejection; invalid/prime-only toolpaths; nonzero/signal/timeout; SIGCHLD-ownership failure; closed-stdio launch classification; transactional capture-limit failure; shell-safe names; deterministic JSON; path collisions. |
| Partition interface | 1 | Valid multi-part declaration; unique IDs; source/part/report hashes; rigid transforms; traversal, missing artifact and symlink escape rejection; manifest fingerprint binding. |

## Final local verification gate

The accepted record requires all of the following on the final committed tree:

```text
RelWithDebInfo focused R8 suite                         45/45 passed
RelWithDebInfo complete headless suite                 100/100 passed
R8 suite repeated until-fail:20                        passed
Independent Debug complete headless suite              100/100 passed
Viewer configure/build                                 passed
Changed R8 translation units with -Wall -Wextra
  -Wpedantic -Werror                                   passed
Recursive clean-clone headless configure/build/CTest   100/100 passed
Recursive clean-clone viewer configure/build           passed
```

These commands establish deterministic local implementation contracts only.
R9 remains responsible for official models, paper metrics/figures, performance
alignment, real process profiles and physical outcomes.
