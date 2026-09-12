# Shared source-material boundary verification — 2026-09-12

`SourceConformingBoundary` now gives adjacent coarse/fine material cuboids the
same indexed interface. Hanging edge points are included across perpendicular
faces as well. Shared triangles are stored once and consumed with opposite
orientations, preventing inconsistent per-cell triangulation during deformation.
Original outer geometry, actual material primitives, fill owners and full BBW
handles are preserved. This discretization is `IMPLEMENTATION_DEFAULT`, not
recovered author triangulation metadata.

## Verified in this checkpoint

- The baseline retains 496 material primitives and all 26 handles, producing
  1,968 shared patches, 2,965 points and 7,944 triangles.
- Independent unit-face coverage and closed-boundary checks passed on three
  inputs, including a transformed nonsymmetric shape and small genuine depth-9
  source geometry. The latter is not a level-9 optimization experiment.
- Independent ten-moment integration agreed before and after full-handle BBW
  motion: maximum relative error `4.163e-12`, within `1e-10`. Independently
  reconstructed point motion differed by at most `4.441e-16`.
- Each input ran twice, with byte-identical scientific exports. All 18
  deliberately corrupted result records were rejected by the independent checker.
- The scoped 2014 suite passed **34/34**, with no skipped or unrun tests, in
  226.63 seconds. This is not an all-project test result. The actual exact-OFF
  dependency build also excludes the new source and target.

Implementation commit: `ca9f976`. Fresh independent review approved the task
with no blocking findings. Constructor readability and existing dependency
CMake notices are recorded as maintenance items. Source, binary and retained
log hashes were checked against the tested implementation.

## Remaining work

Shared indexed closure is not a proof of global manifoldness, absence of
self-intersections or adequate shell thickness. The next integration must
connect these real source surfaces to analytic mass-moment derivatives and
joint fill/angle/full-handle optimization, then enforce geometric commit gates.
Paper regularizers, final material sealing, original assets/parameters and
original-paper experiment matrices remain mandatory.

These results are **SYNTHETIC_REGRESSION**, not original-paper experiments or
physical validation. Overall strict reproduction remains **PARTIAL_BLOCKED**.

See [implementation contract](spin2014_source_conforming_boundary.md) for topology,
coordinate, resource and identity semantics.
