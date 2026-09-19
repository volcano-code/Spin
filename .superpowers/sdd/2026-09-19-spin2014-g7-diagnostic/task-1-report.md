# Task 1 TDD report

## RED

Created `experiments/teapot_spin2014/test_g7_transition_lineage.py` before the implementation. Its first focused run failed during collection with `ModuleNotFoundError: No module named 'g7_transition_lineage'`, which is the expected missing-feature failure.

The tests cover inclusive fractional classification at both source boundaries, exact eight-child split and equal-binary sibling merge, ten-moment conservation, repeated split/merge ancestry, malformed keys, duplicate cells, hash mismatch, non-finite moments, conservation failure, stage-20 partial evidence, and the required exact/inferred/unavailable distinction.

## GREEN

Implemented the minimal read-only analyzer in `g7_transition_lineage.py` and controller in `build/spin2014-g7-diagnostic-20260919/diagnose.py`. The analyzer rejects malformed or unauthenticated synthetic inputs and cannot label net snapshot differences as exact transition history. The controller uses no optimizer, resume, retry, checkpoint publication, or source-run write path.

Focused verification command:

```text
python -m pytest -q -p no:cacheprovider experiments/teapot_spin2014/test_g7_transition_lineage.py
10 passed in 0.02s
```

The controller also rehashed the frozen source identity, validated all 19 checkpoint boundaries through the read-only v2 path, and streamed all 19 sealed native traces. It excluded the stage-20 partial attempt. `python -m py_compile` succeeded for both new Python modules.

## Fixed-QP replay

The two slowest completed stages (19 then 18) were selected mechanically from recorded native time. The saved trace reader selected the largest peak checked QP, with earliest-iteration tie behavior: stage 18 iteration 279 / 412 faces and stage 19 iteration 351 / 414 faces.

After verifying the G7 frozen source hashes against the current compile sources, the original fixed-QP implementation and an instrumentation-only timing build each ran twice per selected stage: eight calls total, single-threaded, 24 GiB address-space cap, and 180-second per-call cap. Each output was science-equal within its stage across both modes and repeats; independent first-order QP KKT checks passed. The full evidence and sub-timers are in `build/spin2014-g7-diagnostic-20260919/REPORT.md`.

No fixed-QP timing was represented as whole-stage timing or a nonlinear speedup.

## Postflight evidence packaging

The diagnostic batch also contains `fixed-qp-replay-receipt.json`: the
machine-readable eight-call receipt has complete trace hashes, selections,
commands, measured times/subtimers, science hashes, independent KKT values,
and the temporary build recipe and hashes. The readable report now includes
all 19 authenticated native timing values and a stage 12--19 saved-final-
snapshot lineage table. The latter explicitly labels every split/merge count
as inferred net state only and leaves residual changes and transaction order
unavailable where the frozen evidence did not persist them. Its ten-moment
delta is explicitly partition-cell geometry only, not a material/carve moment
or constraint-conservation assertion.
