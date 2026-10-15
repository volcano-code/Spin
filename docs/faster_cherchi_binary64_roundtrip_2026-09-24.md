# Faster Cherchi Boolean: binary64 round-trip audit (2026-09-24)

## Scope

This diagnostic investigates the geometric discrepancy between the pinned Cherchi Boolean reference and the independent CGAL comparison path. It does not change the vendor source, solver, scientific tolerances, Boolean backend selection, or formal experiment inputs. The probe's optional exact-point sidecar and the 17-significant-digit OBJ serialization are diagnostic/infrastructure changes only.

The earlier Spin-It joint-deformation and dual-material findings remain documented in [the strict geometry and material audit](../experiments/spin_it_reproduction/strict/GEOMETRY_FAILURE_AND_DUAL_MATERIAL_AUDIT_2026-09-24.md): the deformation fixture reaches a legal-domain boundary and stalls without proving infeasibility; reconstructed Sumone's empty high-density insert is caused by the binary material result, not by export omission.

## Serialization regression and correction

The adapter originally wrote coordinates with `%.16g`. A round-trip regression using arbitrary binary64 box endpoints found 11 coordinate components that did not reload to the original binary64 values. The adapter now uses `%.17g`, which round-trips finite binary64 coordinates exactly. Standard and exact-diagnostic probe modes produce byte-identical OBJ output; the sidecar is separate and does not alter the ordinary output.

The focused regression in `experiments/cherchi_reference_probe/test_probe.py` passes. The previously recorded `.16g` comparison is historical evidence and must not be substituted for the current run.

## Exact-point diagnostic

The pinned Cherchi implementation's final-output path calls `getApproxXYZCoordinates(...)` in `external/reference-sources-20260917/sources/InteractiveAndRobustMeshBooleans-7bd6c2697695f555f9a15d7e154e084d841316f0/code/booleans.cpp` (around lines 1340–1350). Its implicit-point type also exposes `getExactXYZCoordinates` in `code/implicit_point.hpp` (around lines 868–900). The probe records, for each retained output point, its exact rational coordinates and the corresponding in-memory binary64 coordinates. OBJ output vertices are sorted, so matching must use the recorded coordinate key rather than sidecar row order.

For the current synthetic fixtures, the exact rational coordinates and the in-memory point differ as follows:

| Surface / operation | Different vertices | Maximum coordinate difference |
|---|---:|---:|
| Ellipsoid | 0 / 170 | 0 |
| One-sheet hyperboloid | 109 / 232 | 8.22e-15 |
| Two-sheet hyperboloid | 65 / 120 | 7.77e-16 |
| Plane | 2 / 12 | 6.66e-16 |
| Paraboloid | 20 / 31 | 2.22e-16 |

This establishes that approximate implicit-point evaluation contributes to coordinate variation in these fixtures. It does not establish that it is the only source of the Boolean discrepancy.

## Current exact-solid comparison

Current evidence is in `build/faster-cherchi-binary64-roundtrip-20260924-c`. The runner completed with `SYNTHETIC_REFERENCE_BOOLEAN_GEOMETRY_DISCREPANCY`; the project remains `PARTIAL_BLOCKED`.

| Surface / operation | Cherchi mesh | Exact set relation vs CGAL | Symmetric-difference volume | Normalized difference |
|---|---:|---|---:|---:|
| Ellipsoid / subtraction | 170 V / 332 F | equal | 0 | 0 |
| One-sheet hyperboloid / subtraction | 232 V / 464 F | different | 5.360275232887973e-16 | 6.700344041109961e-17 |
| Two-sheet hyperboloid / intersection | 120 V / 232 F | different | 1.3883993807219497e-16 | 1.735499225902435e-17 |
| Plane / intersection | 12 V / 20 F | different | 5.181040781584064e-16 | 6.476300976980073e-17 |
| Paraboloid / subtraction | 31 V / 58 F | different | 1.1031594163489501e-16 | 1.547689807951172e-17 |

An additional diagnostic-only mesh was produced by converting each retained exact rational point to its correctly rounded binary64 value. This is not the ordinary reference output and is not accepted as a solver result. Against CGAL, ellipsoid and plane become exactly equal; one-sheet hyperboloid, two-sheet hyperboloid, and paraboloid remain unequal, with symmetric-difference volumes respectively `1.7385484433342e-17`, `2.852118249060038e-18`, and `6.346711659150132e-18`. Thus exact-point rounding removes a substantial part of the measured volume discrepancy in these fixtures but does not eliminate the set mismatch for all surface types.

### Residual triangulation diagnosis

The reusable diagnostic `experiments/cherchi_reference_probe/analyze_triangulation_discrepancy.py` then compared the rational-rounded Cherchi meshes directly with the CGAL outputs. For each remaining mismatch, the two meshes have exactly the same binary64 vertex-coordinate set, but different unoriented triangle connectivity:

| Surface | Same vertex set | Shared triangles | CGAL-only / Cherchi-only | Exact set relation after point restoration |
|---|---:|---:|---:|---|
| One-sheet hyperboloid | 232 / 232 | 365 / 464 | 99 / 99 | different; symmetric difference `1.7385484433342e-17` |
| Two-sheet hyperboloid | 120 / 120 | 129 / 232 | 103 / 103 | different; symmetric difference `2.852118249060038e-18` |
| Paraboloid | 31 / 31 | 34 / 58 | 24 / 24 | different; symmetric difference `6.346711659150132e-18` |
| Plane control | 12 / 12 | 14 / 20 | 6 / 6 | equal; symmetric difference 0 |

Every unmatched triangle in those matrices lies within `1e-14` (maximum measured distance `4.45e-16`) of at least one input triangle's plane. This is evidence consistent with different triangulations of nearly coplanar Boolean patches after binary64 point conversion. It is not an exact coplanarity proof, and it does not establish that triangulation is the sole cause. The plane control shows that different triangulations can still represent the same exact solid; exact relation remains the deciding check. The machine-readable output is `build/cherchi-exact-point-probe-20260924-m/triangulation-discrepancy.json` (SHA-256 `07c7c14dec6ff33b319a11969959646f80b678b8c552f53c4487ae5bdd353f42`), generated by diagnostic source SHA-256 `43a76a71badc649c99f09b28ee73c8eafb8d90786152f421d4a9d4d9d465c9b1`. A second independent diagnostic execution reproduced the same JSON hash exactly.

## Verification and evidence identity

- `experiments/cherchi_reference_probe/test_probe.py`: passed; ordinary and exact diagnostic OBJ outputs are byte-identical, sidecar point mapping is checked, and binary64 round-trip is exact.
- `experiments/cherchi_reference_probe/test_safety.py`: passed 3/3.
- Focused strict integration tests: passed 2/2 when launched with `LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6` to select the system C++ runtime required by this build. The first launch without the preload failed at runtime-library loading and is retained as an environment failure, not counted as a test pass.
- Cherchi probe SHA-256: `78edb8cdcf0eae5ccb7014c65215b25d1ed73031de562a716bfc5ed95843e8a0`.
- Probe source SHA-256: `d280b6fd0779fa512522db488b5c8dc41292469093ded35f8268bf3ef5dd6242`.
- Probe regression test SHA-256: `2da28be3d63b8618db86ce7b63d4e336252374b2ab1060265ee6e2e892bd691d`.
- Current runner evidence SHA-256: `bb22346ea35836e70962fd246a86d3b2a277db1a8855f17fb4ea05269acd6cff`.
- Current runner scientific SHA-256: `029c83fdba3df5d36ff50a462bf82558ef50147f623037ac13e1c1bf0de6dc75`.
- The triangulation-diagnostic JSON and utility hashes are recorded above; all four diagnostic cases completed with exact set-relation checks.

## Decision and next work

The 17-digit serialization issue is fixed and regression-tested. The exact-point and face-connectivity investigation localizes the remaining differences to meshes with identical vertex sets but different triangle connections; the measured differences are tiny and near input-face planes, but exact-solid comparison still reports inequality for the hyperboloid and paraboloid cases. Do not add tolerance snapping, silently change Boolean backends, or reuse a Jacobian across topology/closure changes. The nonlinear solver's strict full-rank path uses the CGAL value mesh and a CGAL-derived discrete differential; focused solver/derivative regressions pass 6/6. This establishes internal method-path consistency for that backend, not Cherchi/CGAL equivalence or original-paper reproduction. This synthetic comparison is method-validation evidence, not a reproduction of the paper's original application cases.
