# Actual-source adaptive rebind transaction

`rebindSourceCageAdaptive` connects the immutable actual-source full-BBW paper
functional to one adaptive partition change. It is a native, bounded transaction,
not an adaptive optimization loop, binary material validator, sealed artifact,
paper acceptance criterion, or extension of the legacy fill-only adaptive solver.

```cpp
#include "spin2014/SourceCageAdaptiveRebind.h"
#include "spin2014/SourceCageGuardedSlqp.h"

spin2014::SourceCageAdaptiveRebindSettings settings;
settings.fractionalEpsilon = .05; // IMPLEMENTATION_DEFAULT
// Set settings.geometry and settings.budget explicitly for the application.
const auto audit = spin2014::rebindSourceCageAdaptive(functional, fullState, settings);
if (const auto checkpoint = audit.checkpoint()) {
    spin2014::SourceCageGuardedSlqpSettings solveSettings;
    solveSettings.geometry = settings.geometry;
    const auto solve = spin2014::solveSourceCageGuardedSlqp(
        checkpoint->functional(), checkpoint->state(), solveSettings);
    // Inspect the new solve's independent numerical, geometry and material facts.
}
```

The input consists only of the actual `SourceCagePaperFunctional`, its full state,
and settings. No external partition, BBW field, passing flag, checkpoint builder,
or old geometry audit is accepted. A checkpoint is constructed privately only
after every stage succeeds. `audit.checkpoint()` returns a
`shared_ptr<const SourceCageAdaptiveCheckpoint>` and is null for all failures.
Both its `functional()` and `state()` accessors return const references.
Audit facts have shared const storage; copies retain ownership and moves transfer
ownership without moving facts out of references already obtained from that
storage. Such references still require a live owner. Moved-from audit accessors
throw `logic_error`. A copied, edited facts struct cannot construct a checkpoint.

The eight staged records distinguish preflight, actual input geometry, input full
material evaluation, split, merge, rebuild, actual output geometry, and output
full material evaluation. The split and merge operate on a local copy of the
binding domain and beta. The implementation reconstructs
`SourceCageBinding -> SourceConformingBoundary -> SourceCageMomentFunctional ->
SourceCagePaperFunctional` with every original config, density, frame, source,
and full BBW retained. Original BBW weights remain shared and unchanged. The
reconstruction recalculates point samples, eight-node embeddings, primitive
owners, conforming boundaries, moment and constraint derivatives, and the fill
Laplacian. Static states have no angle; non-static phi and every handle's four
parameters are copied exactly using each functional's actual offsets. Even an
unchanged partition is rebuilt and checked again.

`facts()` exposes the old/new ordered leaf identities and beta, both upstream
transition maps, both actual geometry audits, old/new state/functional/partition/
binding hashes, full material evaluation hashes, and before/after world and spin
ten-moments. `deformedMomentDifference` is after minus before in world coordinates
and is meaningful only when `deformedMomentDifferenceAvailable` is true.
The geometry audits bind the value-only paper evaluations; material hashes bind
separate full derivative evaluations. They intentionally need not be equal.

Transition `removedBefore`/`removedAfter` are **rest geometric** moments, without
material density. Their conservation does not establish conservation of deformed
material moments. Refinement samples the same non-affine BBW map at different
points and can change the represented polyhedron. A subsequent guarded solve
therefore begins with empty numerical history and independently checks actual
geometry again. Prior convergence, KKT, feasibility, or geometry conclusions are
never transferred. Solver timeout, iteration limits and stalls do not establish
convergence or feasibility.

The scientific digest binds settings, old/new identities and states, ordered
transition maps and rest moments, geometry/checker digests, material results and
staged dispositions. It contains no elapsed time. Upstream config hashes can
include partition/dimension-dependent inputs; unchanged config fields do not
imply identical upstream config hashes after refinement.

`ResourceLimit`, `RejectedGeometry`, and `Unknown` outcomes never carry a
checkpoint. Finite-input arithmetic failures are unknown, not geometry proofs.
Invalid state dimensions, nonfinite fields, beta outside `[0,1]`, nonpositive
handle scales, invalid epsilon, or invalid API settings may throw. Allocation
failures are resource dispositions when diagnostic storage remains available;
allocation failure while constructing or finalizing diagnostic evidence may
propagate. Original functional and state remain untouched in all paths.

## Resource scope

The additional `maximumOwnedBytes` default is 256 MiB. Checked arithmetic bounds
the transaction's leaf/beta/state snapshots, split/merge evidence, retained
checker tapes and strings, fixed audit metadata, and simultaneous full-evaluation
arrays. The preflight uses conservative bounds of eight descendants per old leaf
and eight output cuboids per old cuboid, plus source-depth-bounded merge chains;
it can reject a transaction whose eventual allocation would have been smaller.
Vector growth is covered conservatively. Newly known evaluation dimensions are
checked again before allocating the remapped state or evaluating the output.
No configured budget is silently raised. A small fixed failure-diagnostic object
is constructed even when a requested cap is too small to admit a transaction.

The copied adaptive domain keeps its original leaf cap. Its internal temporary
copies and immutable rebuild layers use the existing domain, binding, boundary,
moment, paper and shape caps. Geometry execution uses the unchanged checker
caps, separately from the retained transaction tape. These payload/accounting
caps are **not** a native RSS guarantee, do not account for all allocator metadata
or opaque CGAL arithmetic, and do not replace external process guards:

```sh
ulimit -v 25165824
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
timeout 3600 ctest --test-dir build/task-5c-no-faster-20260909 \
  -R '^spin-it-2014-source-cage-adaptive-rebind$' --output-on-failure
```

## Synthetic regression evidence and remaining negative gap

The regression uses original source cubes, an actual coupled full 26-handle BBW
solve, and a depth-four source with one full `[.375,.5]^3` cuboid. Refinement from
level three to four creates eight cuboids. Independent signed-tetra integration
checks before/after deformed moments; analytic boxes and an explicit global
affine transform check all ten conservation identities in the applicable cases.
Tests cover complete state remapping, recursive exact-common-beta merges,
unequal near-bound non-merges, unchanged paths, rebuilt fill adjacency, static
and angle layouts, repeated scientific hashes, and a new guarded solve.

Actual insufficient-clearance input rejection, owned-byte/leaf/rebuild limits,
output checker exhaustion, input/output derivative-cap exhaustion, and unknown
arithmetic failures cannot publish. An output-only **geometry rejection** was
not attained in the bounded 44-state synthetic cosine-translation search. Its
input rejections and successful input/output checks do not close that gap;
output resource failures are not geometry rejection evidence. The optional
`--probe-clearance` regression argument reruns that bounded search. No clearance
threshold was tuned to floating-point differences. Exact-contact and positive
primitive-overlap negative fixtures from prior work remain open. None of these
synthetic cases is a paper case or optimizer validation.
