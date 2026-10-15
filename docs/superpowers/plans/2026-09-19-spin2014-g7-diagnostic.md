# Spin-It G7 Post-Run Diagnostic Implementation Plan

> **For agentic workers:** Follow the subagent-driven implementation workflow for this task. Use TDD for executable behavior. Keep all prior experiment evidence immutable.

**Goal:** Independently diagnose the verified 19-stage G7 prefix, separating measurable runtime costs from unrecoverable historical timings and auditing adaptive partition transitions without advancing the optimization.

**Architecture:** A new read-only diagnostic batch consumes the completed v2 checkpoint chain and sealed native traces. It emits transition lineage, stage-cost accounting, two preselected fixed-QP profiles, and a report; it never mutates the source run or grants scientific acceptance.

**Tech Stack:** C++17 native SLQP trace reader and fixed-QP replay; Python 3.13 audit/analysis and pytest; existing source checkpoint v2 verifier.

**Spec:** User-approved Spin-It-first plan, Stage 1 (2026-09-19).

## Global Constraints

- Root: `/data/xvxinyao/moment-guided-spin-v1/code/make-it-balance-paper-sim`; preserve the main project and all existing dirty work.
- Input batch `build/spin2014-v2-g7-20260919` is immutable; analyze only its 19 complete stages. Never treat stage 20's partial trace as a completed solve.
- No new nonlinear solve, checkpoint resume, parameter search, scientific tolerance change, or automatic retry in this task.
- All new evidence belongs in `build/spin2014-g7-diagnostic-20260919`; bind every output to its source files, native executable, replay executable, configuration, and checker hashes.
- Keep native/Python source-run provenance separate from diagnostic-tool provenance. Historical stages do not contain internal sub-timers; mark those historical subcosts unavailable and measure only explicitly replayed QPs.
- Preserve fixed QP gates: independent KKT `1e-9`; box/dual `1e-10`. Use per-replay wall limit 180 seconds, address-space cap 24 GiB, single numerical thread.
- Fractional means the source predicate `epsilon <= beta <= 1-epsilon`; `epsilon` is the frozen input configuration.
- This task diagnoses a reconstructed level-7 teapot only. It cannot grant binary geometry, simulation, author-case, or paper-reproduction acceptance.

---

### Task 1: Audited historical adaptive-transition lineage and stage costs

**Files:**
- Create: `experiments/teapot_spin2014/g7_transition_lineage.py`
- Create: `experiments/teapot_spin2014/test_g7_transition_lineage.py`
- Create: `build/spin2014-g7-diagnostic-20260919/diagnose.py`
- Create: `build/spin2014-g7-diagnostic-20260919/REPORT.md`
- Create: `build/spin2014-g7-diagnostic-20260919/SHA256SUMS`

**Interfaces:**
- The lineage analyzer accepts the frozen v2 run directory and emits one JSON-compatible result for stages 1–19. Each stage includes verified input/output partition identities, split and merge parent/child keys, fractional counts by level using the inclusive source predicate, exact recorded ten-moment conservation deltas, and any repeated cell refinement/coarsening lineage.
- Any transition detail not recoverable from authenticated checkpoint/trace evidence is reported as unavailable; inferred net partition differences must be labeled inferred and cannot be presented as exact operation order.
- The controller verifies the existing terminal outer job receipts, original source identity, all 19 checkpoint boundaries, and every complete SLQP trace before analysis. It does not modify these files.
- The report distinguishes whole-job, stage-worker, native-call, and replay-only QP timing. Historical objective/derivative, QP, line-search, split/merge, restore, and trace I/O sub-times are marked unavailable unless directly measured by a separate replay.

- [x] **Step 1: Add failing behavior tests**
  - Inclusive fractional-boundary classification at `epsilon` and `1-epsilon`; values just outside are binary.
  - Parent-to-eight-child split, equal-binary sibling merge, conservation of all ten moments, and repeated refinement/coarsening ancestry.
  - Reject malformed keys, duplicate cells, mismatched partition hashes, non-finite moments, incorrect moment conservation, and stage-20 partial evidence.
  - Test that outputs distinguish exact recorded transitions from inferred net partition changes and unavailable history.
- [x] **Step 2: Run tests and confirm expected failures**
- [x] **Step 3: Implement the minimal read-only lineage auditor**
- [x] **Step 4: Re-run focused tests; validate stages 1–19 and source/checkpoint bindings**
- [x] **Step 5: Profile the two slowest completed stages, selected mechanically by recorded native time (ties: earlier stage). In each selected trace, choose the checked QP with the largest peak working set (ties: earliest iteration). Run current replay and timing build twice each; capture stage identity, QP identity, read time, QP sub-timers, science-output equality, and independent KKT results. Do not call a replay result a full-stage speedup.** The current-source reruns and fine Schur split are bound in `build/spin2014-current-fixed-qp-20260919/current-fixed-qp-schur-fine-profile-v3.json`.
- [x] **Step 6: Write the report and SHA-256 inventory; run the independent postflight checks and `git diff --check`**. The fine-grained supplement is `docs/spin2014_g7_schur_fine_profile_2026-09-19.md`; final profile/source hashes are retained under `build/spin2014-current-fixed-qp-20260919/`.
- [ ] **Step 7: Commit only this task's new plan, analyzer, tests, and diagnostic batch files; do not stage or commit any pre-existing dirty path.**

## Execution addendum: current-source fixed-QP replay

The earlier `build/spin2014-g7-diagnostic-20260919/fixed-qp-replay-receipt.json` was later found to bind an older `SlqpTangentQp.cpp` and `SlqpIterationStore.cpp`. It remains historical evidence only. A supplemental diagnostic was therefore run in `build/spin2014-current-fixed-qp-20260919/`, without modifying the frozen G7 run:

- The current QP was bound to live source SHA `7e25529452956f3cc8773ec0449c312a2bedc76ddb009e60362a2057cd26de62`; the current reader SHA is `1f0e7fbb94f662e994d46437694b28f0e331df4218c487f5eadd1015b9ec9b8a`.
- Stages 18 and 19 were replayed twice each. Both return checked first-order QP results, independently pass the fixed-QP KKT audit, and are deterministic within the current source. Their direction vectors differ from the historical serialized QPs only at about `1e-16`–`1e-15` maximum absolute error; this is not bitwise equality and is not evidence of full nonlinear rerun success.
- Temporary source instrumentation measured the current fixed-QP categories. Instrumented and uninstrumented scientific JSON compared exactly. The temporary patch was removed; the live QP source hash was checked back to the value above and the uninstrumented replayer hash returned to its receipt value.
- The current largest measured group is Schur matrix assembly/cache solve (`0.163`–`0.171` s per selected QP); candidate/certificate work is next (`0.045` s for stage 18 and `0.102` s for stage 19). The Schur group still needs finer decomposition before selecting a safe optimization.
- Details and current-source receipts are in [the supplemental report](../../spin2014_g7_current_source_qp_diagnostic_2026-09-19.md). This only completes the current-source replay supplement to diagnostic Step 5; it does not rebind the historical G7 stages or authorize a new nonlinear run.

## Execution addendum: fine-grained current-source Schur profile

The follow-up split Schur assembly, cache validation (including the exact cached-matrix scan), cache preparation, incremental Cholesky update, stability guard, triangular solve/residual check, and COD fallback. The temporary instrumented sources and compiler/link inputs were snapshotted before restoration; their hashes and those of the profile executable and four runs are recorded in `build/spin2014-current-fixed-qp-20260919/PROFILE_SHA256SUMS` and `current-fixed-qp-schur-fine-profile-v3.json`.

Both fixed QPs matched the same-source uninstrumented scientific JSON exactly and passed independent KKT. In 413/415 cache calls, respectively, there was one cold start and then incremental factor growth; there were no cache rejections, releases, matrix rebuilds, or COD fallbacks. Schur assembly and cache checks are measurable fixed-QP costs, but the historical stage traces lack corresponding sub-timers. No safe algorithm change is justified from this evidence; do not weaken input, matrix, factor, residual, KKT, box, or dual checks. The rejected candidate-tile change remains rejected. Production QP/cache source hashes and the uninstrumented executable were verified restored byte-for-byte.

**Stage 2 ruling:** no new solver repair and no new nonlinear run. Preserve `RESOURCE_LIMIT_WALL` / `INCOMPLETE_JOB`; the fixed-QP profile does not authorize a complete G7 retry or material sealing. Further G7 work needs either stage-wide attribution/remaining-budget evidence or a concrete, reviewable semantics-preserving repair with regression coverage. The reconstructed teapot remains outside original-author case counts.
