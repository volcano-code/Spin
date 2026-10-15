# Independent torque-free simulation audit — 2026-09-21

## Scope and decision

The approved plan requires independent energy, angular momentum and time-step
checks before using a simulator as evidence for optimized objects. This batch
verifies the existing ideal rigid-body simulator on an analytic fixture. It
does **not** simulate the G7 material, whose complete wall-thickness gate remains
`UNKNOWN`. It is `ANALYTIC_REGRESSION`, not an original-paper experiment, printing
validation, contact experiment, or paper admission.

No objective, material, optimizer setting, inertia calculation, differential
equation or integration tolerance was changed. The sole simulator change in
`experiments/reconstructed_heart_top/heart_top.py` adds the already integrated
body angular velocities to the returned trace. Previously only attitude and
self-reported invariant errors were returned, so the independent checker could
not recompute invariants directly from the actual integrated state.

## Independent oracle

The fixture has mass 1, center at the pivot, zero gravity, principal inertia
`diag(0.2, 0.2, 0.3)`, and principal angular velocity `(10,14,40*pi)`. The basis
rotation vector `(0.41,-0.29,0.23)` produces a body-frame inertia with all three
off-diagonal terms nonzero. This is a physically admissible synthetic inertia,
not a measured or reconstructed paper model.

For the axisymmetric torque-free body, the exact orientation is the product
of rotation about conserved spatial angular momentum and an opposing rotation
about the symmetry axis in body coordinates. The checker uses that closed form
and the analytic transverse angular velocity, without calling the numerical
Euler/quaternion right-hand side. It independently recomputes
`E = w^T I w / 2` and `L_world = R I w` from the numerical trace, rather than
trusting the simulator's summary metrics. The oracle deliberately supports the
axisymmetric case only; this is not analytic validation of arbitrary triaxial
motion, gravity or contact.

## Registered tests and review correction

`tests/python/test_torque_free_analytic.py` contains five tests. They cover the
full-tensor trajectory, a spherical closed-form limit, a deliberately dropped
cross-inertia negative control, NaN/mismatched truncation, and consistently
truncated or irregular time grids. CTest registers the suite as
`spin-torque-free-analytic-python`, using the required-suite adapter that rejects
missing, skipped, expected-failure and empty suites. Missing SciPy is a test
failure, never a skip or a pass.

Two tests initially failed because the angular-velocity trace was absent. They
passed after the additive trace field. Independent review then demonstrated
that consistently shortening all returned columns could escape the checker.
Two new negative cases failed before the grid validation was added; all five
tests now pass. The grid is supplied independently by the caller, not inferred
from the reported trajectory. The audit uses `arange(1921)/240`, checks every
timestamp, and evaluates the reference on those requested times. Timestamp
comparison allows 32 binary64 epsilons times the fixed grid scale for equivalent
grid arithmetic; it does not permit missing samples or shorter duration.

The reviewer independently rechecked all 48 saved trajectories and reproduced
all maximum errors below. The original checker is preserved at
`build/spin-torque-free-analytic-20260921/checker-v1/`.

## Actual complete batch

`build/spin-torque-free-analytic-20260921/job/` ran eight azimuths separated by
45 degrees, each with a 2-degree initial body-frame tilt, at maximum integration
steps `1/240`, `1/480`, and `1/960` seconds. Each run lasts 8 seconds, with 1,921
output samples. All 24 combinations were run twice: **48 actual simulations**.
The 24 paired raw trajectory JSON files have identical SHA-256 values.

The solver is adaptive DOP853: these are **maximum-step caps**, not fixed-step
integration. We verify errors against an analytic solution and refinement
consistency, but do not infer a fixed-step convergence order. Below a declared
`1e-11` orientation noise floor, the report makes no order estimate.

All listed acceptance errors have an engineering bound of `1e-8` (radians for
orientation; dimensionless for relative errors and quaternion norm). This is
not claimed to be a paper-specified tolerance.

| Maximum step | Maximum attitude error (rad) | Relative energy drift | Relative spatial momentum drift |
|---|---:|---:|---:|
| 1/240 s | 1.680e-9 | 1.632e-12 | 3.702e-10 |
| 1/480 s | 1.176e-11 | 1.422e-14 | 3.456e-13 |
| 1/960 s | 1.719e-12 | 7.395e-15 | 1.561e-14 |

The complete source-bound batch exited successfully in 106.810145396 seconds,
peak memory 165,896,192 bytes, swap 0. The envelope retained the 24-GiB/60-minute
limits and single-thread numerical libraries. NumPy 2.1.3 and SciPy 1.15.3 were
used; Python/package/native files were fingerprinted, not only version strings.
The five existing heart simulator regression tests also executed successfully.

## Fresh checker re-audit, not a rewritten historical run

The original batch used checker v1. After review, the stronger current checker
revalidated all 48 stored trajectories in a separately source-bound job. Its
authoritative current report is:

`build/spin-torque-free-analytic-20260921/reaudit-v3/payload/report.json`

This job completed in 1.717356580 seconds, peak memory 139,718,656 bytes, swap 0.
It authenticates the original successful job receipt and original input files
(using the preserved v1 checker for that historical identity), binds the current
checker separately, rechecks the entire independent time grid, recomputes the
analytic/invariant errors, and verifies both repeated trajectories. It also
executes the registered five-test suite, including its short simulations. It
does not rerun the full-duration batch or any optimization.

Current checker SHA-256:
`7414c0e0b3600f9863becffa25a3829e6ef8ad5e560a2fc1f55e9b6493316676`.
Repeated re-audit scientific digest:
`056ec2b51dad2a5ba13bde0007dc3a8e7644f2d669ef62e133897d27229b05dd`.

The earlier `reaudit-v2` receipt remains intact. Its `new_optimization_or_simulation_runs=0`
field was overly broad because CTest runs short regression simulations. Version
3 narrows the field to `new_full_duration_batch_simulation_runs=0` and explicitly
records the regression simulations. The earlier script is preserved as
`recheck_v2_preserved.py`; no historical receipt was rewritten.

## Remaining gates

- Complete same-component wall/feature thickness on the G7 object remains open.
- This axisymmetric full-tensor test does not cover arbitrary triaxial motion.
- Fixed-pivot gravity/critical-speed validation needs independent parameter and
  invariant checks; the existing simple slow/fast test alone is insufficient.
- Contact, friction and collision proxies require their separate protocol.
- Qualified G7 PyBullet, same-seal visualization/video and original-paper case
  admission remain blocked. Existing engineering videos are not upgraded.

The next simulation work must distinguish free rotation about the center of
mass from ideal rotation about a fixed pivot; the corresponding inertia origins
must not be interchanged. No material admission follows from this kernel test.
