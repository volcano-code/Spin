# G7 v2 frozen-prefix diagnostic

This is a read-only diagnostic of `build/spin2014-v2-g7-20260919`. It did not run, resume, retry, or alter the nonlinear optimizer, checkpoint carrier, or source artifacts. The job ended at its outer wall limit; stage 20 is partial and excluded.

## Authenticated historical evidence

- Outer receipt: `RESOURCE_LIMIT_WALL`, return code `-1`, elapsed `3600.334457750` s.
- Original source identity: rehashed and verified against `spec.json` (5,242 files).
- Checkpoint carrier: all 19 published boundaries, parent links, commit receipts, and referenced artifacts verified with the v2 read-only reader.
- Complete traces: all 19 sealed SLQP files were streamed through the frozen native reader and reported `ALL_NATIVE_RECORD_CHECKSUMS_VERIFIED`; their record counts run from 10 (stage 1) to 381 (stage 19).
- Stage 20 has no completed-attempt receipt and is not a stage result.

The analyzer records the stage-1 input only as an authenticated source hash because no initial partition snapshot was persisted. For stages 2--19 it validates the saved input/output snapshots, inclusive `epsilon <= beta <= 1-epsilon` fractional counts, and ten-moment partition deltas. The largest absolute stage-19 moment delta is `6.437450399132683e-20`.

The carrier has no authenticated split/merge transaction log. Parent/child relationships derived from adjacent snapshots are therefore `INFERRED_NET_PARTITION_CHANGE`; every operation order and repeated refinement/coarsening history is `UNAVAILABLE`. Nothing in this report treats the inferred relation as an exact execution sequence.

## Historical timing classification

The mechanically selected slowest completed stages are 19 and 18, by recorded native-call seconds (tie policy: earlier stage):

| Stage | Native/stage-worker seconds | Whole-job timing | Historical inner subtimers |
|---:|---:|---|---|
| 19 | 538.347963435 | 3600.334457750 s applies only to the outer job | unavailable |
| 18 | 503.343840861 | same | unavailable |

Historical objective/derivative, QP, line-search, split/merge, restore, and trace-I/O times were not recorded and remain unavailable. Fixed-QP measurements below are replay-only timings; they are not full-stage measurements or speedups.

## Fixed-QP hotspot replay

The saved traces were streamed to select the largest-peak working-set QP among `attempted && checked` records (strictly greater peak retains the earliest tie). A temporary replayer was compiled outside this batch only after the frozen G7 hashes matched the live source files:

- `SlqpTangentQp.cpp`: `3f9412551e4c311d93fa86d8a16eb8b9f025f7cd3935885a6a133e03018c0d56`
- `SlqpIterationStore.cpp`: `31c3ad4fc8a67919ae2b1ce81285365ce8b12e7a9dc13f1a72252569b9b2a90d`
- `SlqpLbfgsSolver.cpp`: `74144a5851fb587e82287229cf08ff0c45db0108722d1e70f9495a9c45b445c5`
- `SlqpAuditJson.cpp`: `6da7b8904bd53bf4a191026973af1eeb4c69cbe0fab99681497c9e8ace1e3b07`

The baseline binary SHA-256 was `b20b29ce956ea86950921333ee4092c36354f5e924694e20888acd27d479e328`. The timing binary SHA-256 was `4777b149cb09059bd694b6471f0b9b03dfefd70bc3728d8920e686cd5f569aba`; it is the same frozen QP source plus only timing probes (instrumented copy SHA-256 `f67bcaf33395a71448975bfdfd49567fc498358c18a92151bba4f89851f02dcf`). Both used one numerical thread, 24 GiB address-space limit, and a 180-second per-call limit.

| Stage | Trace SHA-256 | QP iteration | Peak faces | Reads (s), baseline / timing | QP (s), baseline / timing |
|---:|---|---:|---:|---|---|
| 18 | `547e63bf…17a2061f` | 279 | 412 | 2.424544, 2.438525 / 2.448618, 2.462142 | 0.714621, 0.706872 / 0.702762, 0.705225 |
| 19 | `e47575bf…db5bbe65` | 351 | 414 | 2.364494, 2.374669 / 2.408914, 2.386920 | 0.797605, 0.796057 / 0.795720, 0.800556 |

Each baseline and timing mode ran twice per selected QP (eight fixed-QP calls total). Within each stage, baseline/timing and repeat outputs had the same science SHA-256: stage 18 `8dce31e9b0b640784352ad8856cdb23015adf35d64bed452b98572eb01a8a5f1`; stage 19 `6cd012c13272ad665856a3e3400bc1dcc320a82b4c993806c58f9cf246a8f49b`.

Independent fixed-QP KKT checks passed in every call. Stage 18: stationarity `9.75781955236954e-18`, equality `3.3589012114030475e-17`, and zero box/face/dual/complementarity errors. Stage 19: stationarity `1.9434323941802667e-17`, equality `1.070038358856759e-16`, and zero box/face/dual/complementarity errors. These do not establish nonlinear convergence, binary validity, geometry, physics, or paper reproduction.

Timing-build sub-timers (repeat 0 / repeat 1, seconds) were:

| Stage | RHS | Metric columns | Schur/factor | Candidate/certificate | Other |
|---:|---|---|---|---|---|
| 18 | .000211289 / .000207821 | .024693549 / .024651344 | .190292480 / .194062388 | .429328378 / .428294209 | .058187556 / .057961475 |
| 19 | .000228784 / .000215786 | .058580840 / .058749946 | .194925961 / .195517898 | .483256906 / .487500188 | .058674303 / .058523344 |

`setup_phase` was zero for both selected already-feasible QPs. The six timing categories sum to the timing-replay QP wall time up to clock-boundary overhead.

## Limits

This diagnoses a reconstructed level-7 teapot and an incomplete adaptive prefix. It grants no final artifact, binary geometry, simulation, author-case, or paper-reproduction acceptance.
