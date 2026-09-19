# G7 v2 frozen-prefix diagnostic

This is a read-only diagnostic of `build/spin2014-v2-g7-20260919`. It did not run, resume, retry, or alter the nonlinear optimizer, checkpoint carrier, or source artifacts. The job ended at its outer wall limit; stage 20 is partial and excluded.

## Authenticated historical evidence

- Outer receipt: `RESOURCE_LIMIT_WALL`, return code `-1`, elapsed `3600.334457750` s.
- Original source identity: rehashed and verified against `spec.json` (5,242 files).
- Checkpoint carrier: all 19 published boundaries, parent links, commit receipts, and referenced artifacts verified with the v2 read-only reader.
- Complete traces: all 19 sealed SLQP files were streamed through the frozen native reader and reported `ALL_NATIVE_RECORD_CHECKSUMS_VERIFIED`; their record counts run from 10 (stage 1) to 381 (stage 19).
- Stage 20 has no completed-attempt receipt and is not a stage result.

The analyzer records the stage-1 input only as an authenticated source hash because no initial partition snapshot was persisted. For stages 2--19 it validates the saved input/output snapshots, inclusive `epsilon <= beta <= 1-epsilon` fractional counts, and ten-moment **partition-cell geometry** deltas. The largest absolute stage-19 cell-geometry moment delta is `6.437450399132683e-20`. These are sums of the `moments` carried by the partition cells, not carve/material moments, relaxed-solution moments, or constraint residuals.

The carrier has no authenticated split/merge transaction log. Parent/child relationships derived from adjacent snapshots are therefore `INFERRED_NET_PARTITION_CHANGE`; every operation order and repeated refinement/coarsening history is `UNAVAILABLE`. Nothing in this report treats the inferred relation as an exact execution sequence.

## Historical timing classification

The mechanically selected slowest completed stages are 19 and 18, by recorded native-call seconds (tie policy: earlier stage):

| Stage | Native/stage-worker seconds | Whole-job timing | Historical inner subtimers |
|---:|---:|---|---|
| 19 | 538.347963435 | 3600.334457750 s applies only to the outer job | unavailable |
| 18 | 503.343840861 | same | unavailable |

Historical objective/derivative, QP, line-search, split/merge, restore, and trace-I/O times were not recorded and remain unavailable. Fixed-QP measurements below are replay-only timings; they are not full-stage measurements or speedups.

All 19 recorded native-call values are retained in the frozen `stage.json`
receipts and are reproduced here for inspection (seconds):

| Stage | Native seconds | Stage | Native seconds | Stage | Native seconds |
|---:|---:|---:|---:|---:|---:|
| 1 | 27.749557874 | 2 | 27.499523885 | 3 | 27.794620535 |
| 4 | 27.792089042 | 5 | 29.176579019 | 6 | 31.471876552 |
| 7 | 37.084074772 | 8 | 38.406928557 | 9 | 40.602390074 |
| 10 | 53.205433611 | 11 | 58.624551204 | 12 | 352.109190847 |
| 13 | 211.586756191 | 14 | 275.739037547 | 15 | 256.520036633 |
| 16 | 331.868528843 | 17 | 341.123600418 | 18 | 503.343840861 |
| 19 | 538.347963435 |  |  |  |  |

## Readable late-stage lineage summary

The table reports each stage's **saved final partition**, not a transaction
trace. `Net split` and `net merge` are only adjacent-snapshot set relations,
thus all are `INFERRED_NET_PARTITION_CHANGE`; their order, any intervening
operation, and repeated refine/coarsen order remain unavailable. `f@L` means
the inclusive fractional-cell count at level `L`; entries have the form
`L3/L4/L5/L6/L7`. `fractional below max` is the number of saved fractional
cells at levels below the maximum level 7 that could, by level alone, still be
split. It is not a claim that a particular split was attempted or admissible.
Those counts use each checkpoint's saved final-partition `relaxed_fill` vector.
They do not substitute the current trace's terminal, pre-transition optimizer
state, which was not authenticated here as a matching pre-split/merge snapshot.

| Final stage | cells by level (L3/L4/L5/L6/L7) | f@L (L3/L4/L5/L6/L7) | total fractional | fractional below max | net split / net merge (inferred) | max abs cell-geometry 10-moment delta | terminal constraint residual change |
|---:|---|---|---:|---:|---:|---:|---|
| 12 | 30/46/104/426/8876 | 0/0/80/312/7742 | 8134 | 392 | 54 / 26 | 8.209866841262306e-19 | unavailable (not persisted) |
| 13 | 30/44/56/630/11244 | 0/0/16/508/10078 | 10602 | 524 | 354 / 2 | 2.2578722000247443e-19 | unavailable (not persisted) |
| 14 | 28/40/102/288/15166 | 0/4/64/160/14057 | 14285 | 228 | 512 / 2 | 3.210784265685832e-19 | unavailable (not persisted) |
| 15 | 28/34/96/534/16226 | 0/0/40/382/14967 | 15389 | 422 | 182 / 0 | 2.168404344971009e-19 | unavailable (not persisted) |
| 16 | 28/32/68/402/18868 | 0/0/4/200/17482 | 17686 | 204 | 348 / 0 | 1.8973538018496328e-19 | unavailable (not persisted) |
| 17 | 28/30/74/270/20326 | 0/0/16/66/18785 | 18867 | 82 | 174 / 0 | 2.456395547037471e-19 | unavailable (not persisted) |
| 18 | 26/42/64/266/20828 | 0/12/0/66/19074 | 19152 | 78 | 66 / 0 | 3.2356658585114273e-19 | unavailable (not persisted) |
| 19 | 26/40/72/234/21286 | 0/0/10/32/19260 | 19302 | 42 | 60 / 0 | 6.437450399132683e-20 | unavailable (not persisted) |

Stage 1 has an authenticated input *source hash*, but its initial partition
snapshot was not saved. Its input fractional census, input/output delta, and
initial operation history are therefore unavailable. No table row supplies an
unauthenticated constraint-residual change: the saved terminal diagnostics
only state that execution paused after solve and split/merge, and do not
persist a comparable relaxed constraint residual. Material/carve ten-moment
changes and pre/post-transition constraint shifts are likewise unavailable;
none is inferred from partition-cell geometry.

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

The complete machine-readable receipt, including the full 64-character trace
hashes, all eight commands/results, independent KKT values, binary/build
hashes, and the exact CMake recipe, is
`fixed-qp-replay-receipt.json` in this diagnostic batch. Its SHA-256 is in
`SHA256SUMS`. It is a replay receipt only: neither it nor its stored command
invokes a nonlinear optimizer.

`lineage.json` is the corresponding machine-readable stage-1--19 audit. It
contains every stage-12--19 inferred parent-to-children and children-to-parent
coordinate relation, level/fractional counts, saved partition-cell geometry
moment deltas, and recorded native duration. Cross-transition coordinate-key
overlap is explicitly `INFERRED_NET_RELATION_ONLY` and is never interpreted as
an operation history; absent initial/predecessor evidence is `UNAVAILABLE`.

## Evidence-bounded conclusion

1. The frozen evidence verifies source/checkpoint/trace integrity, 19 complete
   boundaries, and small partition-cell-geometry deltas. It does **not** show
   an authenticated implementation defect in those records, nor does it
   provide material-moment, constraint-residual, or historical inner-timing
   evidence needed to assess one.
2. A split/merge-looking net relation or reappearance of a fractional cell is
   insufficient to diagnose an algorithm defect: without an authenticated
   transaction trace, it remains evidence-insufficient rather than an exact
   split/merge or repeat-refinement history.
3. In the fixed-QP replay only, candidate/certificate time is the largest
   measured timing category and Schur/factor time is next. This is not a
   whole-stage profile and cannot be extrapolated to the nonlinear solver.
4. Any stage-2 work should begin with TDD and a controlled profile or
   optimization of this evidence-supported fixed-QP hotspot only. It must not
   alter thresholds or infer missing transition history.

Independent fixed-QP KKT checks passed in every call. Stage 18: stationarity `9.75781955236954e-18`, equality `3.3589012114030475e-17`, and zero box/face/dual/complementarity errors. Stage 19: stationarity `1.9434323941802667e-17`, equality `1.070038358856759e-16`, and zero box/face/dual/complementarity errors. These do not establish nonlinear convergence, binary validity, geometry, physics, or paper reproduction.

Timing-build sub-timers (repeat 0 / repeat 1, seconds) were:

| Stage | RHS | Metric columns | Schur/factor | Candidate/certificate | Other |
|---:|---|---|---|---|---|
| 18 | .000211289 / .000207821 | .024693549 / .024651344 | .190292480 / .194062388 | .429328378 / .428294209 | .058187556 / .057961475 |
| 19 | .000228784 / .000215786 | .058580840 / .058749946 | .194925961 / .195517898 | .483256906 / .487500188 | .058674303 / .058523344 |

`setup_phase` was zero for both selected already-feasible QPs. The six timing categories sum to the timing-replay QP wall time up to clock-boundary overhead.

## Limits

This diagnoses a reconstructed level-7 teapot and an incomplete adaptive prefix. It grants no final artifact, binary geometry, simulation, author-case, or paper-reproduction acceptance.
