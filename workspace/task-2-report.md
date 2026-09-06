# Task 2 implementation report: independent signed-tetrahedron mass moments

## Status and scope

Implemented the M2-subset foundation on `spin-it-strict-computational-v2`, based on
`62427a8c91c873b32a5fc3b5df73f0291a604386`.  This change does **not** grant a
strict-release pass and does not modify historical scientific JSON.

Owned changes are limited to:

- `python/make_it_stand/verification/__init__.py`
- `python/make_it_stand/verification/mass_moments.py`
- `tests/python/test_independent_mass_moments.py`
- the focused `independent-mass-moments-contract` registration in `CMakeLists.txt`

## Implementation

The new public `make_it_stand.verification` package provides immutable
`TriangleMesh`, `DensityRegion`, `GeometricMoments`, and `MassProperties` values.
It exposes density-free signed geometric integration separately from
density-weighted single-region and multi-region mass-property aggregation.

The integrator uses only the Python standard library (`dataclasses`, `math`, and
`typing`) and does not import or call `_make_it_stand`, `spin.integrate_raw_mass_moments`,
or any solver implementation.  Each oriented face contributes a signed tetrahedron.
The checker accumulates signed volume, first moments, and every entry of the symmetric
raw second-moment matrix.  It derives central inertia as
`trace(Q_c) I - Q_c`, retaining all off-diagonal entries.

Conditioning is separate from acceptance: vertices are shifted by a bounding-box
midpoint and divided by the maximum bounding-box span before integration.  Moments
are rescaled by the appropriate powers of that span.  There is no absolute geometric
volume cutoff that changes merely with units.  Central moments are formed in a local
anchor frame, avoiding catastrophic world-coordinate subtraction under large
translations.  `math.fsum` is used for face and region aggregation.

Validation rejects empty inputs, malformed/non-integer/out-of-range indices,
repeated face vertices, duplicate triangles, exactly degenerate triangles,
unreferenced vertices, nonfinite coordinates, boundary/nonmanifold edges,
inconsistent edge orientation, multiple face fans at a shared vertex, nonpositive
or unrepresentable total oriented volume, nonpositive/nonfinite density, and
nonpositive/nonfinite aggregate mass.  Entire consistently reversed connected
components remain signed, so inner cavity surfaces subtract from outer material.

The focused CTest is registered under
`MAKE_IT_STAND_ENABLE_PAPER_PYTHON_TESTS` but outside
`MAKE_IT_STAND_BUILD_PYTHON`; therefore its registration and execution do not
depend on the native extension being enabled.

## TDD evidence

### RED 1: checker absent before production code

Command:

```text
PYTHONPATH=python:tests/python python3 -m unittest \
  test_independent_mass_moments.IndependentMassMomentOracleTests.test_hand_derived_tetrahedron_moments_and_inertia -v
```

Observed exit 1 before adding the checker:

```text
ModuleNotFoundError: No module named 'make_it_stand.verification'
Ran 1 test in 0.000s
FAILED (errors=1)
```

The first test used the brief's literal tetrahedron volume, first moment, full raw
second moment, and full central-inertia matrix.  No solver output generated its
expectations.

### GREEN 1: hand-derived tetrahedron oracle

The same command after the minimal implementation produced:

```text
test_hand_derived_tetrahedron_moments_and_inertia ... ok
Ran 1 test in 0.000s
OK
```

### RED 2: translation-conditioning regression

Self-review identified that centralizing already-translated world moments could lose
shape-scale information.  A test translated the analytic tetrahedron by
`(1e8, -2e8, 3e8)` and retained literal expected inertia values.

Command:

```text
PYTHONPATH=python:tests/python python3 -m unittest \
  test_independent_mass_moments.IndependentMassMomentOracleTests.test_large_translation_does_not_cancel_central_inertia -v
```

Observed exit 1 against the first implementation:

```text
AssertionError: 0.0 != 0.4875 within 10 places (0.4875 difference)
Ran 1 test in 0.001s
FAILED (failures=1)
```

After centralizing in the mesh-local anchor frame, the same command produced:

```text
test_large_translation_does_not_cancel_central_inertia ... ok
Ran 1 test in 0.000s
OK
```

## Final focused verification

Direct full module suite:

```text
PYTHONPATH=python:tests/python python3 -m unittest test_independent_mass_moments -v
Ran 13 tests in 0.073s
OK
```

The 13 tests cover:

1. hand-derived tetrahedron geometric and mass moments;
2. rigid rotation/translation with non-diagonal inertia;
3. large-translation covariance without central-moment cancellation;
4. analytic box raw moments;
5. concentric hollow box with reversed inner component (`V=7`, `Qii=31/12`, `Iii=31/6`);
6. offset reversed cavity and its shifted COM/inertia;
7. two disjoint density regions with literal aggregate moments;
8. density and unit-scale covariance at scales `1e-6` and `1e6`, normalized to `1e-10`;
9. public import while native extension module names are explicitly unavailable;
10. malformed indices, repeated indices, duplicate and degenerate triangles;
11. NaN and positive/negative infinity vertices;
12. open boundaries and inconsistent winding;
13. pinched shared-vertex fans, empty geometry, reversed/zero signed volume, bad density, and empty aggregation.

Syntax and whitespace verification:

```text
python3 -m py_compile python/make_it_stand/verification/__init__.py \
  python/make_it_stand/verification/mass_moments.py \
  tests/python/test_independent_mass_moments.py
git diff --check
```

Both exited 0 with no diagnostics.

Narrow configured CTest verification (no native build):

```text
cmake -S . -B build/strict-v2-paper -DLIBIGL_SKIP_DOWNLOAD=ON
ctest --test-dir build/strict-v2-paper \
  -R '^independent-mass-moments-contract$' --output-on-failure
```

Result:

```text
1/1 Test #9: independent-mass-moments-contract ... Passed
100% tests passed, 0 tests failed out of 1
```

CMake emitted only the pre-existing libigl minimum-version deprecation warning.

## Explicit limitations

This checker certifies only finite oriented triangle 2-manifold topology and signed
mass-moment integration at the represented floating-point scale.  It does **not**
detect or certify self-intersections, containment/nesting of reversed components,
wall thickness, material overlap between density regions, or general solid validity.
Geometrically coincident but separately indexed components are outside its topology
check.  Consequently this foundation must remain outside automatic strict-release
success until complete sealed-geometry checks exist.
