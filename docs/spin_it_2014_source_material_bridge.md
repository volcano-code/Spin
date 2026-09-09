# Source Spin2014 material bridge

This is a `SYNTHETIC_REGRESSION` bridge from the actual source adaptive solver to
an independently checked `SealedSpinArtifact`. It requires both native Spin2014
and exact geometry support; it is independent of the Faster feature. The project
remains `PARTIAL_BLOCKED`. This is not an original paper case, CAD, deformation,
manufacturing, physical validation, or a target feasibility certificate.

```python
from make_it_stand.verification.source_spin_material import (
    SourceGeometrySettings, SourceSolverSettings, AuditedSlqpSettings,
    SourceBoundarySettings, solve_and_seal_source_spin2014,
)

result = solve_and_seal_source_spin2014(
    locked_problem,
    source_settings=SourceGeometrySettings(
        root_origin=(-1, -1, 0), root_side=2, maximum_depth=4),
    solver_settings=SourceSolverSettings(
        maximum_optimization_level=4, laplacian_weight=1,
        slqp=AuditedSlqpSettings(maximum_iterations=300,
            equality_tolerance=1e-10, projected_gradient_tolerance=1e-8)),
    boundary_settings=SourceBoundarySettings(),
)
```

The public settings are frozen, validate complete numeric types and ranges, and
reject unknown keyword settings. Their defaults are `IMPLEMENTATION_DEFAULT`.
The direct native entrypoint requires complete mappings, rejects the same lossy
scalar/mesh encodings, and returns an owned result with no public constructor.
The native result contains the actual parsed settings, source/partition hashes,
ordered per-pass SLQP records, split/merge transitions, final beta/phi, cuboids,
moment checks, exact geometry facts and resource observations. Unavailable
diagnostics use JSON null with enclosing availability/status fields.

Exactly one design region and uniform density across all regions are supported.
Density and the strictly positive outer-wall request come from the locked problem.
The original, fixed and design surfaces must form the complete authenticated
partition. Original model coordinates are preserved throughout the native solve
and output. The target must have axis `(0,0,1)` and origin `(0,0,0)`. Static mode
requires exactly that origin as its sole support point; rotational modes require
an empty support list. `STATIC_LOWERING` maps to `STATIC_BALANCE`, `YOYO` to `YO_YO`,
and the three `TOP_*` objectives to `SPINNING_TOP`. Unsupported targets are rejected.

The source octree is constructed from copied original outer/design surface arrays.
The private terminal wrapper executes once, without an optimization retry.
The exported original outer bytes, including signed zero, must equal the locked
input. This export requires the cavity boundary not to contact the original
design boundary. Retained design is the unchanged original design surface plus
the reversed signed cavity surface. Existing independent exact partition checks
verify both retained/carved design and fixed/retained final material. Python's
separate signed tetrahedron integration checks the actual native moment claim.

For model moments `(m,b,Q)`, normalization origin `o`, length `L` and reference
density `R`, the claim is transformed by
`m'=m/(R L^3)`, `b'=(b-mo)/(R L^4)`, and
`Q'=(Q-ob^T-bo^T+m oo^T)/(R L^5)`. Geometry is never normalized or roundtripped.
The existing Python moment tolerance is `1e-6`; native geometry checking stays
at `1e-10`.

Inspect `result.status`, `result.seal_error`, and `result.artifact` separately.
`BINARY_MOMENT_VALIDATED` may coexist with a geometry resource/wall failure and
no artifact. `BINARY_VALIDATION_FAILED` and `*_UNKNOWN` never grant artifacts.
Returned arrays and nested identity/evidence objects are copies. Scientific
identity binds the locked input, complete settings, actual native run and material
seal; observations of paths, files and timing are separate.

The reproducible runner is
`experiments/spin_it_reproduction/strict/run_source_spin2014_material.py`.
Run it with the built extension on `PYTHONPATH`, a new `--output-directory`, the
extension's actual `--build-directory`, and one or more
`--dependency-root NAME=/absolute/dependency/root` arguments. These are declared
local roots, not compiler-discovered dependency closure. The runner writes actual
NPY inputs and configuration before optimization, `begin.json`, complete
`native-run.json`, path-independent `scientific.json`, bridge file observations,
the independently persisted/reloaded `sealed/` bundle, `evidence.json` and summary.
Existing output directories are rejected; failures and partial files are retained.
`--maximum-intermediate-complexity 4` exercises the actual native resource failure
without changing the fixture or numerical settings.

Use `PaperRunEvidenceV2.from_json(...)` and `check_freshness(...)` to re-read the
observed files. `CURRENT` is a freshness observation, not solver authentication.
Two output directories have identical bridge, native and artifact scientific
hashes, but provenance digests may differ because they include actual paths.
There is no generic loader that trusts solver-run JSON. Saving/loading
`result.artifact` through `spin_bundle` persists and rechecks its geometry seal
only; it does **not** authenticate a solver execution.

Evidence is `OBSERVED_FILES_ON_DISK_ONLY`, with build/source attestation
`NOT_VERIFIED`, original evidence `BLOCKED_UNKNOWN`, physical validation
`NOT_IN_SCOPE`, and solver-target feasibility `UNKNOWN`. Native Nef complexity
facts cover the native source-material operation chain; they exclude the extra
independent Python seal audits, CGAL allocator internals, and triangle-distance
storage. They do not measure total process memory.

## Retained Task5f execution evidence

Implementation `e1db0ee` passed independent task-scoped review. Full configured
CTest suite:198/198,445.61s; independent controller bridge rerun:10/10.
The source bridge also executed10/10 tests with Faster disabled.

Actual post-review runs are retained under
`build/task-5f-demo-run-1` and `build/task-5f-demo-run-2` (not original paper assets).
Both were independently reloaded and reintegrated: unchanged original outer
bytes,1actual cavity from72cuboids,39ordered optimization records, exact squared
outer-wall clearance25/64. Maximum normalized moment claim error:
`1.3877787807814457e-17`. Inputs, complete native trace, material and cavity geometry
are identical across the two runs. Both evidence checks were CURRENT at acceptance.

Scientific digest: `780431ce4e1970da972488483a575b2b566a61fa072fdc6d09b65b5399c6e6a7`.
Reintegrated artifact digest: `778f8dfb1dd0dc601152e2ffb78d5756efe8c42b1e61c1b3c9c8c03af9f0a97c`.
Full-suite log: `build/task-5f-full-ctest-verbose-01.log`.
Controller independent checks: `build/task-5f-controller-demos-02.log`.
The earlier controller check log is retained: it rejected freshness before the
controller explicitly loaded the required native module; no solver change was needed.

These are sealed synthetic numerical/material results, not physical simulations
or newly implemented OpenGL views. Geometry success does not certify the requested
support/rotation target. Original elephant/heart experiments and full paper
reproduction remain unverified; the overall project remains `PARTIAL_BLOCKED`.
