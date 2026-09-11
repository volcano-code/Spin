# Source material / BBW binding verification — 2026-09-11

The new native `SourceCageBinding` connects actual source-octree material to the
complete BBW field. It does not replace boundary leaves with full bounding cubes,
reduce handles, or allocate a global point-by-handle matrix. Original outer
arrays and root coordinates must match byte-for-byte, and split/merge operations
cannot modify an already-created binding.

## Actual results

- The fixed synthetic fixture retains 8 fill variables, 496 actual material
  subcells, 997 bound points and all 26 handles. Material volume is
  0.030517578125, not the incorrect bounding-cell volume 0.125.
- Independent checking covered both BBW protocols and a transformed nonsymmetric
  model. Each complete set ran twice; exported scientific results were identical.
- Independently recomputed ten-moment error was at most 6.939e-18. All chunked
  weights matched independently calculated eight-node interpolation. Twenty-four
  deliberately corrupted records were rejected by the independent checker.
- The final scoped 2014 regression run passed 33/33 tests in 223.40 s, including
  source-material Python, native source/cage modules, dependency compilation and
  retained-workspace behavior. This is not an all-project suite.

Fresh independent code review approved this scope, with no blocking findings.
The final header, source and linked library hashes still match the independently
tested build. Implementation: `dd94566`; implementation report: `e7ac737`.
Existing CMake dependency deprecation warnings remain disclosed.

## What this does not complete

This is rest-material/interpolation integration, not the complete joint
deformation solver. Deformed shared faces must still form a consistent material
boundary; joint beta/angle/handle optimization, paper regularizers and final
geometry/material sealing remain separate mandatory gates. Original-paper
assets, parameters and computational case matrices retain their existing blockers.

Overall status remains `PARTIAL_BLOCKED`; these synthetic tests do not count as
original-paper experiments, physical validation or complete reproduction.

Implementation contract: [source binding](spin2014_source_cage_binding.md).
