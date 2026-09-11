# Task5o — immutable actual source-material/full-BBW binding

This implements the next prerequisite of the approved strict Spin-It 2014 M3 plan. It is NOT the full joint SLQP, paper regularizer, deformed material, or original-case acceptance. Global status remains PARTIAL_BLOCKED. Execute only after Task5n review passes and controller assigns source/build ownership.

## Scope / public contract

Add `include/spin2014/SourceCageBinding.h`, `src/spin2014/SourceCageBinding.cpp`, one native regression file, scoped CMake registration, and `docs/spin2014_source_cage_binding.md`. Use an immutable shared-private-Data value object, matching SourceCageGrid/Bbw. Accept a const SourceAdaptiveDomain snapshot, a real SourceCageBbw object and explicit resource config; no constructor taking supplied weights, pass flags, or precomputed moment claims.

Binding must:

1. Require actual BBW success and final primal validity. Keep protocol/raw/final assessments through the immutable upstream object; independent-normalized is valid but must not be mislabeled coupled KKT.
2. Verify byte-identical original outer vertices/faces and root origin/side between source oracle and cage. Dimension + scalar-byte equality, including signed zero. Compare actual input arrays, NOT equality of the source digests: oracle includes its separate design domain, cage uses outer-as-design. Preserve both upstream digests.
3. Copy the adaptive partition so later caller splits/merges cannot mutate this binding. Enumerate *actual included Full source cuboids* with existing sourceBinaryCuboids(domain, ones, maximumPrimitives). Boundary leaf bounding boxes MUST NOT substitute for their included subcells.
4. Each primitive stores its original SourcePartitionKey, canonical lower/upper box, exactly one owning adaptive beta index, and eight bound point indices (x-bit-first corners). Find owner by O(depth) ancestor map lookup, not O(primitives*leaves). Every primitive must have exactly one ancestor owner; bind all leaves and retain original order. Snapshot exposes source/domain/BBW read-only.
5. Point layout: original outer vertices first, byte-identical original order, then unique material-corner points sorted by integer coordinates on source.maximumDepth lattice. Deduplicate material corners by exact canonical keys, not float tolerance. Do not deduplicate/reindex the original surface into material corners. All points use actual source canonical endpoints; no snapping to BBW grid.
6. Cache sparse SourceCageGrid::Embedding (8 weights/node indices per point), using actual grid.embed. No global dense point-by-handle or 3N-by4H allocation. Preserve ALL handles and expose count. A bounded `pointWeights(begin,count)` API may return count-by-all-handles using cached eight-node embeddings and actual finalWeights. Handle zero count, out-of-range and overflow explicitly; no clamping/renormalization/handle subsampling.
7. Checked config budgets for maximumPrimitives, maximumPoints, maximumWeightValues (default conservative values). Guard size/Eigen/index overflow before allocation; incremental point maps stop at budget. RESOURCE_LIMIT length_error for resource exhaustion, invalid_argument for invalid/mismatched source/BBW, out_of_range for chunk ranges. No partial binding object, silent coarsening or dense root-grid expansion.
8. Versioned config/output SHA binds source config/content/partition, BBW output/field identity, all budgets, original geometry identity, point arrays, primitive ownership/connectivity/keys/bounds and embeddings. No timestamps/durations. Unrelated external Faster work excluded.

## Tests / required actual behavior

Use TDD: record missing-interface RED before implementation; then actual focused GREEN and the assigned safe scoped2014 suite once. Simple bodies are SYNTHETIC_REGRESSION only.

- Main fixture: original outer cube(.3,.7), design cube(.34,.66), root [0,1]^3, sourceDepth5, initialOpt2/max5; actual source coarse cage2/fine grid3/FEM/full26handle coupled BBW. Controller preflight establishes 8 adaptive leaves, 496 included primitives, actual volume (10/32)^3=0.030517578125 versus bounding-cell volume0.125. Do not tune inputs to reduce primitive count or handles.
- Every primitive exactly one owner and all primitive moment sums per leaf match independently integrated cuboids (not solver source function called twice). Check all ten moments/nonzero cross terms and nonuniform per-owner beta. Outer arrays unchanged.
- Verify bound point coordinates/connectivity/canonical ordering and sparse embedding reconstruction. Chunked weights equal independent eight-node manual interpolation, all26handles retained; use multiple chunks and a memory-limited configuration that cannot allocate the full point-by-handle field but accepts bounded chunks. Empty chunk, invalid range and finite-result behavior.
- Actual repeated complete binding produces identical arrays and scientific digests. Caller mutates partition via split/merge afterward: old snapshot unchanged, new binding changed and beta ownership remapped correctly. Original source/BBW objects unchanged.
- Different design with same exact outer is accepted; different outer ordering/bytes, root mismatch, invalid/failed actual BBW rejected. No public arbitrary-success constructor.
- Resource caps for primitive count, point insertion, weight chunk output; malformed zero/overflow budgets fail before large allocations. No fake upsampling or level9 optimization claim.
- Include one transformed nonsymmetric coordinate fixture if tractable using existing validated FEM/grid range; do not relax numerical tolerances to make it pass.

## Architectural boundaries / explicit ruling

This binds rest geometry and interpolation only. Under nonaffine LBS, independent cuboid triangulations can have hanging/shared-face conformity issues; this object is NOT a final material boundary or a geometry-valid deformation. It does not evaluate joint objective, deform cells, call SLQP, update rotations/regularizers, seal objects or grant feasibility. Those remain next mandatory tasks rather than silently using Task5m's legacy dense adapter. Source approximation uncertainty and implementation-default BBW protocols stay explicit.

Resource/collision rules: work only `/data/xvxinyao/make-it-balance-paper-sim`; preserve main project and six unrelated Faster paths. Use only existing FasterOFF `build/task-5c-no-faster-20260909`, max-j3, address-space cap25165824KiB, no main build. No broad suite containing unrelated destructive dependency scripts. Exact-OFF compatibility must not compile the new source binding. Controller will own independent verification and fresh review. Never spawn agents. Apply_patch edits; scoped commit only owned files under volcano <xuxinyao307@gmail.com>; no tag/push/merge/delete.

Report `.superpowers/sdd/2026-09-05-spin-it-strict-computational-v2/task-5o-report.md` with actual RED/GREEN commands/results, source/binary/input hashes, remaining limits, scoped commit and returned build ownership. Keep noisy dependency warnings visible and never claim full reproduction.
