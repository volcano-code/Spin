# Task5aa — actual boundary-component separation, a thickness prerequisite

Approved parent: Task4 shared final geometry/thickness, user strict specification.
BASE69717a182d32b85812b0d2363d762c3bd26293d3. Task5z complete; do not redo.
Work only canonical isolated repo /data/xvxinyao/moment-guided-spin-v1/code/make-it-balance-paper-sim.
Protected main sibling and all foreign Faster dirty/untracked paths stay untouched.
No dependencies/downloads/merges/pushes/tags/newworktrees/deletions or agents.
Use apply_patch,TDD, preserve raw failures, user identity volcano <xuxinyao307@gmail.com>.
Only one implementer; you own build/task-5c-no-faster-20260909 -j3 until release.
External24GiB/60min, single-thread runtime. Controller owns independent checking,
remaining-work inventory and final broader suite; do not modify its files.

## Deliverable

Add include/spin2014/SourceCageComponentClearance.h,
src/spin2014/SourceCageComponentClearance.cpp,
tests/spin_it_2014_source_cage_component_clearance_regression.cpp,
docs/spin2014_source_cage_component_clearance.md and scoped CMake registration.
No changes to Task5z/backend/Python seal semantics.

Entry validateSourceCageComponentClearance(const SourceCageMaterialRegionsAudit&,settings).
Immutable privately constructed audit and optional immutable validated component
product retaining the entire original regions audit (and transitively original
binary audit/material/moments). Only genuine validatedRegions() may admit. No
candidate-only/arbitrary mutable meshes or copied facts can grant publication.

From the ACTUAL final material mesh and Task5z exact boundary faceComponents,
extract every closed indexed boundary component. Canonical order = existing
component ID; component faces in original face order, vertices compacted in
ascending original global vertex ID. Retain source vertex IDs, source face IDs,
original sign and parent/depth facts; independently confirm all source faces and
vertices are covered exactly and hash actual parent mesh against upstream exact
inputDigest. Original actual material is never altered. For inward shell only,
reverse a checker-only component copy, explicitly recording reversal. Construct
one standalone bounded positive ExactMaterialRegion per component using the
existing budgeted backend. Preserve each exact input/region digest and boundary
facts. Do not assume a cavity is a selected voxel primitive; components may have
many triangles/selected primitives.

For every unordered component pair i<j, compute existing exact boundary
clearance on those standalone regions using a finite nonnegative
minimumComponentSeparation in model units. Retain complete exact rational
squared distance, threshold, witness component-local face IDs AND corresponding
parent material face IDs. Contact rejects even threshold0. Noncomplete or unknown
never passes. Record all pairs if resource permits even after measured rejection;
never emit a successful summary based on only a prefix. Pair counts checked with
overflow and budgets before allocation/work. A single boundary component (e.g.
empty cavity solid) has an explicit NOT_APPLICABLE result, zero pairs, no invented
global distance/witness. Multiple components require exactly n*(n-1)/2 measured
pairs; save global minimum and canonical first-pair witness for ties.

Use one local ExactNefBudget shared by standalone component regions; only copied
final stats escape. Bound components/points/faces/pairs/known-owned bytes/exact
diagnostic bytes; check cap validity/overflow before variable allocations. Account
input faceComponents, component maps and stored per-pair diagnostics. Resource
failures/invalid inputs/unknown results remain distinct and never publish. Native
CGAL/rational allocator internals remain under external limit, not an RSS promise.
Deterministic complete scientific hash binds settings/upstream/material/component
maps/reversal/region+boundary facts/pair and aggregate diagnostics/outcome/budget;
no timestamps. Retain original bytes and moments across originating-audit lifetime.

## Required scope labels

This is BOUNDARY_COMPONENT_SEPARATION only. It measures geometry gap across any
two distinct material boundary components (including root-root gaps through air),
NOT a universal material-wall thickness definition. It does not check different
features of the SAME connected component, medial axes, continuous BBW Jacobian,
author continuous design, full thickness, sealing, feasible authority, physics or
paper convergence. fullManufacturingThickness MUST remain NOT_ESTABLISHED even
when this bounded gate passes. No deformed SealedSpinArtifact grant or Python
hollow-only override. Later full feature-thickness and sealing stay open.

## TDD cases

Real immutable actual source/BBW/adaptive/binary/regions fixtures: empty cavity
(one component); one cavity (two); two separated actual cavities (three). All
full handles and at least one nonidentity deformation / nonzero phi. Native
Task5z fixture construction available in its regression: duplicate fixture setup
only if needed, do not duplicate production extraction/integration logic.
Verify all component provenance/maps and all unordered pairs; a thin inter-cavity
gap must be rejected when outer gap passes (choose predetermined actual source
fixture and engineering thresholds, not paper input). Test representable equality
and nextafter when supported by actual identity-handle fixture literal expectation.
Reject upstream region rejection, budgets0/1/SIZE_MAX, invalid threshold, missing
or copied authority; lifetime/deterministic repeat; distinguish one-component N/A
from missing measurement. Pure analytic fixtures unit-only, zero objective known
stationary fixtures explicitly KNOWN_SOLUTION_VALIDATION. Do not claim original
paper experiments. If exact whole general component contract cannot be supported,
escalate without silently weakening it.

Report task-5aa-report.md (same SDD dir) with raw RED/GREEN commands/log paths,
API, units/scope/budget details, selfreview and limitations. Scoped commit only the
five authorized files after focused tests; report finalSHA and release build.
Send finalized public header/API early to controller for independent raw export.
