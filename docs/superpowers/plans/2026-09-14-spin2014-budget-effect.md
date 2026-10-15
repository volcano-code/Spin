# Spin-It 2014 budget control and reconstructed-effect audit

Approved by user continuation: keep the frozen native metric-tangent binary,
increase only SLQP iteration capacity, and check whether the reconstructed
teapot achieves the effects discussed by Spin-It 2014.

## Global constraints

Use the existing isolated worktree. Preserve dirty work, failed records and
all prior evidence. No commits, merges, deletion, tags or author-original claims.
Do not change solver/source/geometry/thresholds during running experiments.
Same G7/O7 input, 65536 leaf capacity, objective TOP_INERTIA_ONLY, no warm start.
Only increase maximum_iterations from960 to1920. Two sequential fresh runs share
a3600-second case-batch cap;24GiB address-space cap, one thread, minimum100GiB free.
If time/resource cap prevents completion, retain the attempted run and mark it
RESOURCE_LIMIT. No completion or physical-effects claim from partial moments.

### Task 1: Budget control

Verify baseline hashes and focused tests. Run the unchanged stream_probe.py with
the separately named1920 configuration twice under the shared deadline.
Rehash input/source/binary and scientific outputs, compare actual archives.
Record all statuses, KKT, binary checks, geometry admission and resource limits.

### Task 2: Effect audit

Read official2014 paper tea-pot result and equations, not2017 overview.
Extract initial solid, initial feasible-domain mixture, final continuous and
binary moments where available. Distinguish these material states explicitly.
Compute COM, full inertia about pivot and COM, target-axis alignment, inertia
ratios and applicable theoretical critical speed with units/assumptions.
Do not interpret inertia-only results as lowering/full-objective experiments.
Any diagnostic checker must have analytic/nonfinite/state-admission regressions
and independent review. Missing sealed geometry blocks final physical claims.

### Task 3: Evidence report

Combine actual results and an independent scoped review. Separate qualitative
effect agreement, numerical constraints, valid final geometry, simulation,
original-paper numerical reproduction. Do not derive a progress percentage.
Update remaining-work report with the actual next blocker.
