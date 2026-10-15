# Faster prescribed-waterline input adapter — 2026-09-25

This is a synthetic method-integration result. It is not an author floating
case, a free-waterline solution, a fluid simulation, or a reproduction-complete
claim.

## Added path

`with_measured_buoyancy_geometry` now builds a `SpinItFasterRequest` from two
closed meshes:

- a caller-designated displacement boundary, whose sealed/floodable semantics
  remain the caller's explicit modeling choice;
- the actual fixed-shell mesh, integrated at the request's material density
  to derive fixed-shell moments and volume.

The displacement mesh is transformed into the request target frame, clipped by
the measured prescribed waterline, and used to fill displaced volume, center
of buoyancy, waterline area, and central second moments. The adapter rejects
left-handed/non-orthonormal target frames and waterline normals not aligned
with target-frame positive vertical. It does not infer cavity flooding or
validate a full material partition.

`solve_spin_it_faster_buoyant_plane` calls the existing restricted plane-family
solver and returns its trace, plane coefficients, material-boundary mesh,
constraint results, and numerical Eq. 7 flags. Equilibrium-constraint status
and metacenter stability remain separate outputs. The binding uses the
existing `solveBuoyantStabilityByPlane`; it does not pad its four-coefficient
plane family into the general ten-coefficient Newton system.

## Verification

- The new adapter tests were first run against the missing API and failed at
  the expected availability assertion.
- The boolean-density validation test failed before its guard was added, then
  passed after the change.
- The nontrivial box fixture constructs waterline data from geometry, runs the
  plane solver from a non-solution initialization, and checks convergence,
  nonzero accepted steps, numerical equality flags, plane coefficients, and
  analytic total mass / horizontal first moments.
- Current build rebuilt `_make_it_stand` successfully.
- `tests/python/test_spin_it_faster.py`: **22 passed**.
- CTest `spin-faster-buoyant-plane-balance`: **1/1 passed**.
- `py_compile` and `git diff --check`: passed.

The fixed-shell adapter currently requires the fixed shell to use the request's
single material density, matching `SpinDesignRequest`'s buoyancy contract.
Multiple fixed-material densities and chained Floating Island A→B provenance
are not implemented by this adapter; A→B still requires the separate approved
interface design. The native result is not a `SealedSpinArtifact`, and the
test does not constitute independent full-geometry, manufacturing, or author
asset validation.

## Historical-status correction — 2026-09-27

The paragraph above records the adapter state on 2026-09-25 and is superseded
for the later synthetic A→B path. The current
`solve_floating_island_b_buoyant_from_sealed_a` chain composes A's sealed
retained geometry and original fixed regions into B, measures B's displacement
geometry from its closed hull, runs and seals B's prescribed-waterline plane
solve, and binds the source A scientific/evidence digests. The current material
bridge regression exercises A solve/seal → B solve/seal → recomputed tilted and
reverse-pose hydrostatic checks. This remains synthetic method validation; it
does not supply the authors' Floating Island geometry, pose, or protocol.
