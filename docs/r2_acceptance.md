# R2 mass-property and projected-energy acceptance record

Implementation audit date: 2026-08-12 (Asia/Shanghai)

Acceptance status: **accepted**

## Scope

R2 replaces the inherited center-of-mass update with a paired mass-property
state:

```text
m = signed uniform-density mass
Q = integral rho*x dV = first moment
c = Q/m
```

`MassProperties` stores only `m` and `Q`; center of mass is derived and is not
maintained as an independently mutable quantity. The implementation provides:

- signed-tetrahedron surface integration for a closed, consistently
  edge-oriented triangular boundary;
- positive outer-surface and negative inward-cavity contributions;
- an exact uniform axis-aligned-box contribution;
- paired subtraction of each carved voxel's mass and first moment;
- a full grid-state recomputation after every carving step;
- normalized-gravity projected distance and paper Eq. 2 energy;
- explicit failure for undefined CoM, invalid box extents, invalid density,
  zero gravity, and solver inputs that violate its mesh/grid preconditions.

The production carving solver uses density `rho=1`. Its reported mass is
therefore signed volume in model-coordinate units, not calibrated SI mass in
kilograms. `computeSurfaceMassProperties` and the box helper support a uniform
positive density multiplier, but R2 does not support heterogeneous density.

## Equation interpretation

For an oriented triangle `(vi,vj,vk)`, R2 uses the standard
signed-tetrahedron contribution (shown below with a zero reference point for
clarity):

```text
mf = rho * vi dot (vj cross vk) / 6,
Qf = mf * (vi + vj + vk) / 4
```

up to an arbitrary reference-point translation used to reduce cancellation.
Summing these contributions is equivalent to the volume-integral definition
and divergence-theorem construction in the paper Appendix. Reversing every
face negates `m` and `Q` while preserving a defined `c=Q/m`.

There is a printed-factor discrepancy in the paper's symmetric Eq. 3 form. If
the raw cross product `((vj-vi)x(vk-vi))` is used, the displayed coefficient
`rho/6` produces three times the physical signed volume because the dot product
with `(vi+vj+vk)` repeats the same scalar triple product three times. The
physically consistent coefficient for that symmetric expression is `rho/18`.
The paper's volume-integral definition, its Appendix, and the beta source's
one-coordinate divergence form agree with the standard signed-tetrahedron
result used in R2. This discrepancy is recorded rather than reproduced.

For nonzero gravity `g`, R2 evaluates

```text
ghat = g / ||g||
delta_perp = (I - ghat*ghat^T) * (c - c*)
d = ||delta_perp||
E_CoM = 0.5 * d^2
```

The gravity direction is normalized after scaling by its largest absolute
component. Thus every finite nonzero direction—including magnitudes near
`1e-300` and `DBL_MAX`—has the same projection without norm overflow or
underflow. The solver continues to compare distance `d`, rather than squared
energy, with its inherited voxel-length target tolerance so dimensional units
remain consistent. Target construction and safe-region geometry remain R3
work.

## Solver preconditions

`CarvingSolver` now requires:

- a finite, closed, consistently edge-oriented triangular outer boundary;
- outward/global orientation yielding positive signed mass;
- a built, initially filled axis-aligned voxel grid with positive extents;
- finite target data and nonzero gravity.

The topology check requires every undirected triangle edge to occur exactly
twice with opposite directions. R2 rejects invalid solver input; it does not
repair open, non-manifold, inconsistently oriented, or inward outer meshes.

The incremental solver state subtracts

```text
mv = rho * dx * dy * dz
Qv = mv * voxel_center
```

for each carved axis-aligned voxel. The per-step oracle independently scans
the complete voxel state and accumulates these contributions in `long double`,
starting from the separately unit-tested initial outer-surface properties. It
is independent of the production incremental deletion accumulator; it is not a
second independent implementation of the initial surface integral.

## Registered test matrix

R2 registers four analytic physics tests:

| CTest | Contract | Verification status |
|---|---|---|
| `mass-surface-analytic` | Cuboid, tetrahedron, density scaling, global winding reversal, translation, rotation, uniform and tiny scaling, inward non-concentric cavity, and undefined zero-volume CoM. | Passed |
| `mass-voxel-analytic` | Anisotropic axis-aligned box, density scaling, aggregate/remove/add behavior, deletion prefixes, zero remaining mass, and invalid extent rejection. | Passed |
| `mass-deletion-exhaustive` | All `8! = 40320` deletion orders of an eight-voxel fixture; every prefix compares incremental `m`, `Q`, and `c` with a fresh sum. | Passed |
| `balance-energy-analytic` | Axis-aligned, oblique, parallel, tiny and near-`DBL_MAX` gravity scaling; invalid-gravity rejection; `E=0.5*d^2`. | Passed |

It also registers four carving regressions:

| CTest | Contract | Verification status |
|---|---|---|
| `carving-initial-state` | Zero-step state, initial `(m,Q,c)`, gravity-scale invariance, header-only CSV and valid solver initialization. | Passed |
| `carving-one-step` | Ordinary positive-candidate step, state/count continuity, exact voxel subtraction, oracle agreement and CSV round trip. | Passed |
| `carving-full-mass-consistency` | Production `carveAll()` with every-step long-double `m/Q/c` oracle checks; no hard-coded R0/R1 divergent floating-point target. | Passed |
| `carving-no-candidates` | One empty search, no mutation, `no_candidates` terminal reason and idempotent repeated `carveAll()`. | Passed |

The active carving CSV has exactly 39 columns. They cover iteration; before,
after, and recomputed CoM; before, after, and recomputed first moment;
projected distance and Eq. 2 energy; CoM/first-moment discrepancies; initial,
before, after, and recomputed mass plus mass discrepancy; candidate/removal and
inside-voxel counts; carving ratios; and stop reason. Every serialized numeric
field is round-tripped to the corresponding in-memory diagnostic value.

## Numerical tolerances

The focused analytic tests use absolute and relative tolerances of `1e-12`.
The Bunny carving regressions use scale-aware absolute-plus-relative checks:

| Quantity | Absolute tolerance | Relative tolerance |
|---|---:|---:|
| Mass and first moment | `1e-12` | `1e-10` |
| Derived CoM | `1e-10` | `1e-10` |
| Projected distance and energy | `1e-12` | `1e-10` |

These tolerances test numerical agreement rather than cross-platform bitwise
identity.

## Final clean-verification gate

The final gate used the pinned R0 toolchain and completed on 2026-08-12:

```text
RelWithDebInfo headless build and CTest                  8/8 passed
RelWithDebInfo full suite repeated until-fail:10         10/10 complete suites passed
Independent Debug build with assertions enabled         8/8 passed
Viewer configure/build                                  passed
Committed-tree clean-clone headless build and CTest      8/8 passed
Committed-tree clean-clone viewer build                  passed
Main repository and pinned libigl clean-state checks     passed
```

The historical R0 defect trace remained byte-identical at SHA-256
`3f9108649fc329a136a5371597335a99275021b6ddd9a3984198ce2b8e8c1856`.

## Conformance status ceiling

- signed mass is at most `VERIFIED_UNIT` at R2;
- signed first moment is at most `VERIFIED_UNIT` at R2;
- normalized-gravity Eq. 2 is at most `VERIFIED_UNIT` at R2;
- derived CoM reaches `VERIFIED_REGRESSION` only because every solver step is
  checked against the separately accumulated long-double voxel-state oracle;
- none of these rows reaches `VERIFIED_PAPER` at R2.

The R0/R1 CSV remains unchanged under `tests/baselines/` as historical defect
evidence. R2 deliberately retires its divergent numbers from executable
assertions rather than relaxing a tolerance to preserve them.

## Remaining non-conformances

R2 does not change the fixed batch of at most 1000 removed voxels, the
whole-grid empty-ratio denominator, target-point construction, the best-prefix
policy, temporary-energy-increase handling, or the paper's too-few-removals
termination rule. It also does not implement support polygons/safe regions,
multiple targets, wall-thickness validation, connectivity guarantees, BBW/LBS
deformation, Laplacian preservation, analytic gradients, alternating
optimization, mesh repair, or fabrication validation.

Accordingly, R2 does not establish objective monotonicity, successful balance,
static stability, printability, paper-faithful carving, or a complete
reproduction of *Make It Stand*.
