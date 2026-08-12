# Characterization baselines

`bunny_r0_known_failure.csv` is the exact 21-step diagnostic trace captured
before R1 added candidate-exhaustion termination. It is intentionally retained
as historical evidence of the R0 defect:

- candidates first reach zero at iteration 17;
- the solver then records four no-op iterations;
- it stops at iteration 21 only because of the legacy maximum-iteration rule;
- geometric projected error grows from about `0.0289962` to `0.811397`;
- the incremental/voxel-proxy CoM discrepancy reaches about `0.800743`.

This file is not a target output and must not be interpreted as evidence of
correct balancing. R1's live full-run regression stops at iteration 17 with
`no_candidates` while preserving the same known-divergent numerical state for
R2 to replace with mass/first-moment correctness tests.

R2 has now retired every hard-coded value from the active regression suite.
The live full-run test compares mass, first moment and derived center of mass
with an independent full recompute after every step. This CSV remains unchanged
only as historical defect evidence and retains its original SHA-256.
