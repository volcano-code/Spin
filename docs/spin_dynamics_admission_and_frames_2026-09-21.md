# Material admission and rotational reference frames — 2026-09-21

## Decision and scope

The G7 material's full feature-thickness gate remains `UNKNOWN`. No G7 physics
job was launched in this batch. Work addresses two prerequisites for the
approved same-material simulation stage, without changing the input, optimizer,
rounding, tolerances, or historical experiment outputs.

The legacy `experiments/teapot_spin2014/reconstructed_g7_dynamics.py` only checked
`MESH_MATERIAL_VALIDATED`, which does not establish full wall/feature thickness.
Its torque-free report section used pivot inertia under zero gravity. That
historical protocol explicitly disclosed a fixed ideal joint and suppressed
translation; it must not be reinterpreted as free rotation about the mass
center. The new code separates the two problems.

## Full material gate

The worker first rejects incomplete overall wall status before expensive exact
integration. This metadata preflight can only reject. For a potentially admitted
bundle it then calls `load_qualified_dynamics_bundle`, which uses the existing
`load_spin_artifact` to reconstruct and independently revalidate the actual
geometry and moments, checks the expected scientific digest, and requires:

- validated material geometry;
- `same_component_feature=VERIFIED` and `overall=VERIFIED`;
- outer-to-newly-carved clearance verified, or not applicable for empty carving;
- inter-cavity web verified, or not applicable for fewer than two bounded voids.

The existing seal checker cannot yet produce that complete result. This change
therefore does **not** create a newly admitted object. In particular a sealed
analytic box with a verified outer/cavity gap remains blocked when its full
feature thickness is unknown. Promoting JSON fields and recomputing their hashes
does not work: the independent geometry reconstruction contradicts the forged
facts. This is a fail-closed admission test, not a positive end-to-end qualified
material test; the latter still requires a real complete thickness checker.

## Explicit rotational origins

`simulate_rotational_protocol` always takes body-frame *central* inertia `Ic`.

| Protocol | Inertia supplied to Euler's equations | Gravity lever | Gravity |
|---|---|---|---|
| `FREE_CENTER_OF_MASS` | `Ic` | zero | zero |
| `IDEAL_FIXED_PIVOT_GRAVITY` | `Ic + m*((c dot c)*E-c*c^T)` | `c` from pivot to mass center | 9.81 m/s² |

Both are rotational models; translational motion is explicitly not integrated.
The kernel does not grant material admission. It rejects invalid central inertia
before a parallel-axis addition could conceal it, and requires an explicit
supported protocol. The existing worker now calls the appropriate branch and
labels free-rotation results with `CENTER_OF_MASS`.

Tests compare the fixed-pivot branch with a hand-computed full tensor, including
the three nonzero cross terms. Free orientation is invariant under changing the
model-coordinate origin while retaining central inertia. A synthetic offset
body demonstrates that the old zero-gravity pivot motion and free center motion
differ by approximately 0.0542651 radians over 0.125 s. This is a regression
fixture, not a new measurement of the tea-pot experiment.

## Native runtime integration

Independent review identified a launch dependency: a clean worker environment
did not contain the native geometry extension's directory. The launcher now
selects the binary explicitly from the frozen run spec, checks its ABI filename,
size and SHA-256 against the sealed native observation, and supplies its directory
in `PYTHONPATH`. The observation's `path` value is a label (`native`), not a
filesystem location; it is never used as a path to load.

The selected extension and resolved Linux C++ runtime are added to source-job
identities; the command explicitly preloads that C++ runtime. Additional loader
source dependencies are included. A clean-environment child-process test imports
the actual binary, and a mismatched native digest is rejected before import.
This runtime preparation alone does not authenticate geometry or permit physics.

## Executed evidence

Tests were added before implementation. Three initial material-admission tests
failed for the absent gate, then passed. Four frame tests failed for the absent
protocol interface, then passed. Review-driven runtime tests initially failed,
then passed after the launcher was fixed.

Two source-bound batches are retained separately:

1. `build/spin-dynamics-admission-20260921/job/`: four registered CTest groups,
   including existing sealed artifact/bundle regressions, passed before the
   runtime fix. Elapsed 79.603719036 s; peak 104,722,432 bytes; swap 0.
2. `build/spin-dynamics-runtime-20260921/job/`: the current three dynamics groups
   passed, containing five analytic-oracle tests, four frame tests and five
   native admission/runtime tests (**14 individual tests**). Elapsed
   21.265540405 s; peak 135,757,824 bytes; swap 0. Both jobs retained the
   24-GiB/60-minute resource limits and single-thread numerical libraries.

The latest independent reviewer reproduced the five native admission/runtime
tests and found no remaining actionable issue in this bounded change. This does
not replace a full-project release audit.

Current file hashes:

- dynamics script: `b5afd97694ebdc061fa12cd527beeb45b7c017625c8153bbe0fbf4b8ea15bc40`
- admission tests: `9259c07ad090f5c459c0320522d163a6452ccdfa883cfec397b291f406efa82f`
- reference-frame tests: `79f8d47874074526408cf42fc804d3276bc433307de7b8a85a4c1493878e1808`

Pre-change and pre-runtime-fix source copies are retained under
`build/spin-dynamics-admission-20260921/before/` and `before-runtime-fix/`.
The first batch is historical evidence for that earlier source, not a claim that
the later source already passed the same full group set.

## Remaining work

The legacy worker still points to historical frozen G7 batch paths. It is not a
prepared launch package for the new G7 seal. Its original stale-source checks
were not disabled to make it run. Preparing the next real material experiment
must bind the current admitted seal and sources in a new, independent directory.

Complete thickness, qualified material dynamics, the prescribed full perturbation
matrix, applicable critical-speed verification, contact/collision validation and
same-object video are still required. No original-paper experiment, fabrication
validation, or complete reproduction is granted by this batch.
