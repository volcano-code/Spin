# Task5u — actual-source represented geometry admission foundation

Continue approved M2/M3 Task5 precommit geometry work, BASE9341e5c. This task
implements the checker first; solver line-search integration is the following
dependent task, not falsely claimed here. No original-case/physical/release grant.

## Scope and ownership

Canonical worktree /data/xvxinyao/moment-guided-spin-v1/code/make-it-balance-paper-sim,
branch spin-it-strict-computational-v2. Protected main remains read-only.
Own include/spin2014/SourceCageGeometry.h, src/spin2014/SourceCageGeometry.cpp,
tests/spin_it_2014_source_cage_geometry_regression.cpp,
docs/spin2014_source_cage_geometry.md and minimal CMakeLists.txt gated registration.
Report task-5u-report.md beside this brief separately committed. No generic
solver, existing geometry library, Faster/Canonical or foreign untracked edits.
If upstream change needed, ask controller with concrete evidence first.

## Required behavior

Typed immutable checker accepts actual SourceCagePaperFunctional plus explicit
settings and full raw state, then calls its actual evaluateValues. Do not accept
caller-provided geometry/evaluation/pass booleans as evidence. All original
outer rows, all conforming primitive boundary references and all handles stay.
Keep the five objectives compatible, but use world material geometry for checks.
Hash complete raw state, source/frame/settings, actual output, every attempted
geometry stage/outcome and resource facts. No timestamps in scientific digest.
Public facts must not grant FEASIBLE_WITNESS, convergence, stable spin, sealed
artifact or original-paper success; distinguish CHECKS_PASSED for the declared
represented-polyhedral checks, REJECTED_GEOMETRY, RESOURCE_LIMIT and UNKNOWN.
No result object with mutable public constructor capable of forging admission.

Checks on represented binary64 geometry using existing pinned exact primitives:
1. Deformed original outer mesh with unchanged faces: closed, oriented signed
   nested material, no degenerate faces, edge/vertex manifold, no self intersection.
   Preserve any original holes/nested components; do not force every shell to
   positive orientation or silently repair/reorient input arrays.
2. Each conforming primitive boundary: construct compact indexed local mesh from
   ALL stored OrientedTriangle references and actual deformed points; correct
   orientation, positive nonzero volume, closed/manifold/nondegenerate and no
   self-intersection. Keep primitive IDs and beta owners. This is a represented
   polyhedral-cell validity check, not a proof of the continuous BBW-map Jacobian
   throughout a cell. Never use rest-normal dot products that reject a valid
   rigid rotation as a flip.
3. Every deformed primitive is a subset of the actual outer material region,
   not merely inside its AABB. No positive-volume primitive-pair overlap;
   face/edge/point sharing is allowed and cannot be misclassified as overlap.
   Pair identity/order and counts must be checkable. Exact-safe AABB broad phase
   may avoid impossible overlaps, but cannot replace the actual set predicate
   for candidates. Use existing budget-aware Nef functions, no alternate backend.
4. Exact boundary separation between outer and each primitive: explicit finite
   nonnegative minimumOuterDomainClearance in model units, default0; outer contact
   is rejected even when threshold0. This protects the entire optimization domain,
   independently of beta. Name this conservative domain/outer clearance, not
   final binary cavity-web or all-feature shell thickness. A fractional beta is
   not a geometrically realized half-cell; no binarization here.

All geometric checks must complete before CHECKS_PASSED. Failure/incomplete
Boolean/source/resource evidence stays non-admitted; absence is never true.
Invalid API dimensions/settings may throw before an attempted audit; actual
source/geometry failures retain stage and deterministic failure category where
possible. Resource exceptions cannot become ordinary infeasibility. Do not keep
running geometry after resource latch. Document native allocator/opaque CGAL
limits honestly; no exact RSS claim. Do not swallow evidence into success.

No vertex snapping, smoothing, mesh repair, topology edits, primitive dropping,
partial handles, changed objective/tolerances or release shortcuts. Source input
is immutable. Cache immutable rest topology only if bounded; deformed checks
cannot reuse stale pass facts for changed state. This is an engineering guard
layer, not a recovered author-protocol assertion.

## Resource design proposal required before production

Propose finite positive caps for point/triangle/primitive counts, checked pair
counts, retained evidence and explicit mesh scratch, exact Nef live complexity.
Compute all integer/byte products with overflow checks before size-dependent
allocation, including index narrowing into existing int APIs and worst-case
pair/evidence storage. Avoid dense N*N allocations; sequential primitive/pair
work is acceptable with explicit time/memory cost. Existing ExactBoundaryAudit
BVH/CGAL internals do not expose exact native allocator caps, so state that
limitation and require process24GiB/60min. Do not invent a resource-proof claim.
Submit short API/accounting proposal to controller; then TDD production.

## Verification (real independent-source fixtures only)

Use existing actual SourceSurfaceOctree -> conforming -> full26handleBBW pipeline
from Task5t small fixtures, not mocks or direct fabricated success structs.
Write behavioral tests before production, actual registered RED then GREEN:
- identity and common positive scale/translation pass with analytically predicted
  volume and domain separation; repeated scientific output identical.
- nonuniform beta and phi changes do not alter geometry verdict, but bind state.
- real nonaffine valid motion, actual outer or primitive degeneracy/inversion/
  self-intersection failures, outer contact/protrusion and primitive overlap when
  actual attainable with the full controls. Do not fake omitted failure coverage;
  report any case that cannot be realized and resolve with controller.
- preserve holes/void semantics where feasible via existing actual-source inputs.
- exact requested clearance vs one-above boundary, caps exact/one-below, zero,
  overflow, bad inputs; small cap failure cannot give CHECKS_PASSED.
- typed result construction/copy semantics and changed input/settings identities;
  independent direct integration and analytic clearance checks for actual output.
No known-solution test or synthetic model counts as original paper reproduction.
Controller will separately export/check geometry/facts and mutate actual outputs.

Existing build build/task-5c-no-faster-20260909 only, FasterOFF,2014ON,exactON,
LIBIGL_SKIP_DOWNLOAD=ON. At most3 compiler jobs to leave controller one; all
heavy runs sequential, OMP/BLAS/MKL1, ulimit-v25165824, timeout3600.
Baseline source-cage-slqp and source-conforming tests are run by controller.
Run focused new tests while iterating, final related2014 registered group once
after frozen source. Verify current source-Python evidence and actual exactOFF
dependency-child graph excludes new component. No foreign builds or installs.
Preserve logs/failures/source+binary hashes. Scoped commits only, no merge/tag/push.
Read implementer prompt/TDD instructions, no subagents or worker-owned reviewer.

## Approved API/resource supplement before production

Immutable SourceCageGeometryChecker(functional,settings).check(fullState), private
result construction with const getters. Defaults proposed/approved:100kpoints,
400ktriangle references,4096primitives,1M unorderedpairs,1.1Mstages,512MiB known
checkerbuffers,2M liveNef complexity,64shells. All are IMPLEMENTATION_DEFAULT.
Retain ordered fixed-size stage evidence with IDs/owners/predicate results and
complete fact digests. Account actual output mesh copies, largest local boundary
diagnostic worst face-pair list/edges/components, upstream shell n*n scratch,
known region mesh copies and all retained records. Opaque rational/BVH/CGAL
allocator internals remain explicitly outside exact RSS accounting. Independently
exportable exact clearance strings require finite per-string/total limits and
resource rejection, not truncation-as-success. All integer widths checked.

Reuse region.boundaryFacts() to avoid duplicate audits if its full properties
suffice; each primitive requires one positive component, outer retains signed
nested semantics. AABB skip only when strictly separated. No per-pair dense table.

Task5t's depth5 fixtures contain7944triangles/1968patches; new synthetic geometry
regressions may explicitly use source depth3 with the same full26handle pipeline
and8beta, fixed before results, to bound Nef validation cost. Dyadic outer
5/16..11/16 and full interior cells3/8..5/8 give exact identity separation1/16.
These are newly declared unit fixtures, not replaced prior experiments or
reduced-resolution paper acceptance. Moderate refinement check if affordable,
otherwise explicit resource rejection. Source failure before geometry (e.g.
negative mass) stays UNKNOWN; do not manufacture downstream orientation facts.
