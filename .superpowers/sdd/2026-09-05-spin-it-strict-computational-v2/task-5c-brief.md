# Task 5c — immutable original-surface boundary subtree foundation

BASE: c07663bba4bb8b82b3f8cde9e31e4eac302dafe0. This is a bounded part of approved Task 5, not its completion.

## Binding constraints

- Work only in `/data/xvxinyao/make-it-balance-paper-sim`, branch `spin-it-strict-computational-v2`. Preserve dirty main checkout and all historical results. No merge/push/tag/release.
- Analytic meshes are unit/engineering tests only, never original paper experiments. Overall `PARTIAL_BLOCKED` remains.
- Test-first; preserve RED and GREEN logs, no weakened tolerances or retuned inputs. Compile at most 4 concurrent jobs (prefer 3); experiments sequential, 24 GiB / 60 minutes, free-space requirements from approved plan.
- You are sole production source/configure/build writer. Do not spawn agents. Ask controller about substantive ambiguity. Commit only scoped files using configured identity. No broad staging.

## Scope and architectural ruling

Create a native immutable geometry oracle directly from original outer and explicit design-region triangle meshes. It must not accept or reconstruct its source from occupancy. It is a prerequisite for later adaptive optimization, not that optimizer. Keep existing AdaptiveVoxelDomain source() contract, solver/cage paths and golden behavior unchanged. No null-source adapter, dense export, Python API expansion, mesh repair, Boolean fallback, or new cage implementation this increment.

Paper Section 5.4 explicitly describes boundary-cell geometry subtrees. It does not specify every finest-cell sampling convention. This increment uses an explicit independent `IMPLEMENTATION_DEFAULT`: any finest-level closed cell intersecting the design boundary is excluded from the carveable approximation, and its whole volume is reported as unresolved boundary uncertainty. This is conservative geometry discretization, not exact cut-cell integration, author's recovered rule, continuous infeasibility proof, or level-9 optimization completion.

## 1. Value-only published surface integration

Add a public value-only unit/density-aware integrator next to SurfaceRawMomentDifferential. Factor the existing Algorithm 1 polynomial accumulation into shared code so value-only calls allocate/compute no dense vertex Jacobian. Existing differentiateSurfaceRawMoments semantics/results must remain compatible. Keep validation and signed-component integration contract explicit. Geometry legality is separate.

Tests: independent analytic translated asymmetric tetrahedron, box and reversed inner cavity; cross moments; density/translation covariance; invalid/nonfinite input; agreement with old differential is compatibility only, not independent proof. Include an actual allocation-sensitive regression (resource-limited subprocess or measured allocator), not a hand-written claimed allocation count. Never let value-only merely call the full differentiator.

## 2. SourceSurfaceOctree geometry oracle

Use focused new `include/spin2014/SourceSurfaceOctree.h` and `src/spin2014/SourceSurfaceOctree.cpp` (a private helper file is allowed if required to isolate CGAL template compilation; explain). Public immutable PIMPL/value-object API owns copied outer/design V/F, root cube, geometry maximum depth and node budget. Expose const source surfaces, source/config/partition SHA-256, exact original outer/design surface moments at density 1, approximate carveable moments, unresolved boundary volume, level/cell counts and queryable hierarchical cells. No caller-provided pass fields.

Root cube: explicit finite origin and positive finite side; both original boundaries must fit, no silent clipping. Level coordinates use integers. Supported depth 0..20, default max-node budget 5,000,000; reject numerically collapsed dyadic coordinates and overflow/nonfinite moment/volume arithmetic. Node budget is checked before allocation; exceeding it fails construction explicitly (`RESOURCE_LIMIT` diagnostic/exception), never yields a valid partial/empty result. Depth 9 must be genuinely traversed in a test.

Validate actual meshes with existing exact boundary/material-region audits and exact design-subset-outer relationship. Reject malformed/open/self-intersecting/unreferenced/misnested/outside input; preserve valid nested cavities and disconnected components. No reorientation or repairs. Reuse exact material region as validation only, not an independent solver moment checker.

Classify cells against original design triangles with pinned CGAL exact triangle–closed-box intersection (including touching). Spatial acceleration must not miss triangles entirely enclosed by a cell. With no boundary intersection, determine actual material inside/outside including nested cavities; then store a full/empty terminal. Boundary cells subdivide to the requested depth. At maximum depth, store unresolved boundary terminal with zero carveable moments and whole-cell uncertainty. Parent moments/uncertainty sum actual children. Keep genuine hierarchy/lineage; children in fixed deterministic order. No dense uniform grid loops/512-cubed arrays, no upsampling, no center/corner-only boundary test.

Queries: expose stored nodes/children and lookup a dyadic cell at any level <= maximum depth. For a descendant of a coarser full/empty terminal, return analytic cuboid/zero contribution without dense subdivision. Invalid coordinates fail. Query moments must partition consistently. Stats distinguish internal/boundary/full/empty nodes and actual maximum visited level. All moments geometrical (rho=1), fixed before optimization. Hashes bind copied source bytes, configuration and actual tree scientific contents; timestamps/measured resources excluded. The maximumNodes input guard IS included in configuration hash, since it governs success/failure; actual RSS/time are not. Explicit docs on approximate moments and uncertainty (volume only, not a magical ten-moment relative-error certificate).

Build this API only when BOTH 2014 and exact-geometry checks are enabled; existing 2014 without exact stays supported. It must not secretly depend on Faster. Add actual native linked/runtime feature-isolation evidence with Faster OFF. Use existing vendored CGAL, no dependency upgrade/download.

## 3. Fixed geometry acceptance tests

- Root [0,1]^3, outer unit box, design box [1/4,3/4]^3. Closed boundary-touch exclusion gives approximate box [1/4+h,3/4-h]^3 at depth L, h=2^-L; at L=9 side=254/512 and volume=(254/512)^3. Verify all ten analytic moments, actual depth-9 nodes, sparse node count/memory/time, exact surface volume 1/8, and exact-minus-approx volume <= computed boundary uncertainty. Native analytic normalized error <=1e-10 for these nondegenerate dyadic fixtures.
- Add a reversed inner cavity [7/16,9/16]^3 to the design shell; cells containing its triangles cannot be classified full solely from a sample point. Test disconnected small islands and a non-axis-aligned/asymmetric fixture; independently check analytic surface moments/nonzero cross terms and conservative volume enclosure.
- Parent/children and descendant-of-terminal queries consistent. Finite root containment and resource/node failures fail closed. Source-copy immutability: mutate caller buffers after construction and verify saved source/hash/query unchanged; separately modified new mesh changes source hash. Two identical builds give identical scientific hashes.
- Bad geometry (open, flipped/misnested, self intersection, outside design, NaN), invalid depth/coordinate/root and numeric collapse reject. No missing dependency skips.

## 4. Verification and deliverables

Add focused native CTest(s), standalone value-only allocation test as appropriate, public algorithm/provenance documentation. Record TDD commands/output, actual level9 geometry stats and resource limits/observations, actual linked executable/library hashes, feature isolation and main all-target build/full CTest once after focused iteration. All registered tests must run; retain failures. Do not rerun full suite per edit. No new original-paper case claims or sealing status.

Write full report to task-5c-report.md beside this brief. Include implementation, RED/GREEN, exact commands/results, scoped files, self-review, resource/build identities and remaining limitations. Reply under 15 lines: status, commits, test summary, concerns, report path. Controller will provide fresh independent review; do not self-dispatch a reviewer.
