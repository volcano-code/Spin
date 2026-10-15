# Distinct bounded-void boundary separation (2026-09-21)

## Scope

This adds `make_it_stand.verification.inter_cavity_clearance.audit_inter_cavity_clearance`
to the shared Python verification package. It supplies a missing geometric
observation for multi-cavity Spin-It objects. It does **not** finish the current
G7 same-component feature audit, change the optimizer, modify an input, or grant
paper/manufacturing/physics admission.

The threshold must be supplied explicitly. The frozen `minimum_outer_wall`
contract concerns the outer surface and newly carved region; it is not silently
redefined to cover other distances. No existing sealed bundle format or old
facts are changed. A future full wall checker must explicitly bind its chosen
criterion and threshold, then consume this observation in that stated scope.

## Geometric construction

1. Reauthenticate the complete represented material using the existing exact
   CGAL region checker: closed embedded boundaries, signed orientation, strict
   nesting, and Nef material/void component counts.
2. Every negative boundary shell identifies one bounded void. Its **full**
   boundary contains that shell plus its immediate positive child shells,
   which bound material islands inside the void. Deeper negative shells belong
   to other voids, not the original one.
3. Reverse winding and compact unused vertices only in temporary void copies.
   Reauthenticate each copy as a connected finite exact region. Original
   material arrays are never edited.
4. Enumerate every distinct bounded-void pair and use the existing exact
   triangle-distance engine. Compare squared distances and the represented
   threshold squared as rationals. Threshold equality passes.
5. Report the global minimum and map its witness triangles back to **original
   material face IDs**, including material-island faces when appropriate.

This checks all surface points on the separate void boundaries, not just
vertices or sampled rays. It does not cover pairs on the same void boundary,
same-component local thin features, outer-surface features, or CAD error.

All-pairs separation at least the threshold provides a sufficient Euclidean
separation bound for distinct voids. A smaller separation alone is **not** a
material-side thickness defect: the connecting segment can cross exterior air
between separate material bodies. Both outcomes retain
`complete_wall_thickness=UNKNOWN` and manufacturing/physics `NOT_GRANTED`.

## Failure and evidence rules

- Zero or one bounded void yields `NOT_APPLICABLE`, with no invented distance,
  witness or successful threshold comparison.
- `maximum_pairs` is mandatory; exceeding it returns `RESOURCE_LIMIT` before
  pair enumeration. It does not purport to bound native mesh authentication
  costs; actual experiments still use the 24 GiB/swap-0/60-minute job wrapper.
- Failed native operations do not publish a partial minimum or successful
  comparison. `pairs_checked` is diagnostic only.
- The report binds the source material identity, threshold, budget, native
  region digests, checker file observations and all returned facts. Disk
  observations are not a claim of loaded-code/build provenance.
- Old sealed facts remain immutable; this is a supplemental observation, not
  a shortcut around `load_spin_artifact` and complete artifact revalidation.

## Tests and review

The initial seven tests were first run against the missing checker and failed
at its absence. After implementation all seven passed. Two further real native
tests exercise reordered face/vertex indices with an independent rational
triangle-distance oracle, and an exterior air gap that cannot be upgraded to a
wall failure. All nine passed through the newly registered CTest test
`exact-inter-cavity-clearance-python` (8.24 seconds).

The nested fixture catches a specific dangerous omission: checking negative
shells alone measures squared distance `25/16`, whereas the material-island
boundary gives the true `1/16`. The sheared fixture gives the analytic squared
distance `1/20`. Three cavities require all three pairs; a cap of two cannot
return the partial minimum.

A separate read-only review found no blocking defects and independently tested
native failure after one completed pair, shuffled nested-shell face order and
an exterior-air connection. Those failure injections are review evidence, not
yet a retained automated regression test.

## Source-bound G7 run

The separate run package is
`build/spin2014-inter-cavity-audit-20260921/`. It binds the latest G7 seal
manifest SHA-256 `7aa3f73d13fd9157322ddb79a18aeefb7a5fb8ba91f17b6b6f5f304fad6d4bc7`,
checks its material NPY bytes and represented geometry, repeats the supplemental
observation twice, and runs the targeted exact-geometry/sealed-object CTest
groups. It does not recompute the complete seal's mass integration and does not
replace that seal.

The job completed with `PROCESS_EXITED`, exit code 0, in **193.391948160 s**.
The enforced memory maximum was 24 GiB; observed peak memory was **433,983,488
bytes**, with swap maximum and observed peak both 0. Both G7 observations are
byte-identical with audit digest
`95d049d2993217f168b0dd08986dc2444a0d04529f10afdf0b147c477548509a`.
The native recheck confirms one bounded void: `NOT_APPLICABLE`, no distance or
threshold comparison. This does not change the remaining same-component
feature-thickness blocker.

All **9 registered targeted CTest groups** actually ran and passed in 112.35 s:
three native exact-geometry groups, four Python exact-geometry groups (including
the new nine-case separation suite), the sealed artifact/bundle group, and the
sealed-view data group. This is not a full-project, graphical desktop, or
PyBullet test claim. Required groups were checked against the JUnit report;
missing, skipped, or failed groups are rejected by the runner.

Evidence:

- `build/spin2014-inter-cavity-audit-20260921/job/result.json`
- `build/spin2014-inter-cavity-audit-20260921/job/payload/report.json`
- `build/spin2014-inter-cavity-audit-20260921/job/payload/tests.xml`
- `build/spin2014-inter-cavity-audit-20260921/job/payload/ctest.log`
- `build/spin2014-inter-cavity-audit-20260921/job/payload/g7-repeat-{0,1}.json`

Checker SHA-256: `6921346a48dffd85e5f25ed8e0fd0917894d145444e2df6c6bd3d38b69607f60`.
Test SHA-256: `f0049cd732b3b533a68ef66d67a4c02c1be7211be4faf7c59c67ff3fac9bec8f`.
The job binds source/checker files, the native extension, test executables,
CTest registration and original seal payloads. The original seal manifest
hash was independently rechecked after the run and remains unchanged.

## Remaining gate

This closes the **distinct bounded-void separation implementation and analytic
verification** subtask, not the complete wall-thickness requirement. The
current G7 object has only one void and gains no new wall-thickness admission.
The next material gate remains continuous same-component thin-feature coverage
with an explicit defensible metric; vertex-only and one-nearest-point checks
are not sufficient. Original-asset, level-9, full simulation, paper case matrix
and final release requirements remain open.
