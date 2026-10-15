# G7 plan execution — diagnostics and direct-audit pilot (2026-09-21)

## Decision

The first from-initial-state G7 attempt remains **`INCOMPLETE_JOB`**: its expanded JSONL export was truncated at the 3,600 s deadline, and it is preserved without retroactive completion. A second, separately fingerprinted G7 run completed its native terminal check and direct trace audit in **3,583.948 s**, within the same 3,600 s limit. This closes the evidence-pipeline pilot, but not the geometry seal, deterministic repeat, or paper-reproduction gates.

The evidence supports an evidence-pipeline bottleneck, not a demonstrated adaptive split/merge or QP correctness defect. No scientific parameters or partition semantics were changed. Direct streaming audit of native traces with independently verifiable per-trace receipts has now completed one full G7 pilot. No sealed artifact, repeatability claim, or paper-reproduction claim is granted.

## Run and integrity state

- Run directory: `build/spin2014-g7-instrumented-20260921/repeat-0`.
- Started from the configured initial state; no old checkpoint was imported.
- Limits: 24 GiB memory, swap disabled, numerical libraries single-threaded, 3,540 s inner and 3,600 s outer wall.
- At launch, all 5,112 frozen source/input/build identity entries matched.
- Attempts 0–19 have durable stage-boundary receipts. Attempt 20 has a terminal native header and complete pass-20 `.slqptrace`, but no `stage.json`, terminal result receipt, stage audit, or root `results.json`.
- The final header reports `BINARY_MOMENT_VALIDATED`, normalized binary-constraint residual `2.51736511197818e-7`, and normalized moment error `8.93e-15`. Its geometry state is `NOT_ATTEMPTED_UNKNOWN`; this is not final-geometry validation.
- The attempt-20 `native-run.jsonl` reached 5,692,028,642 bytes (about 5.30 GiB / 5.69 GB) and ends in a truncated record. It is incomplete evidence and is not read as an archive. The required archive readback and `stage_audit_v2` audit never ran.
- The 20 completed paused workers account for 2,706.813 s of native time. Attempt 20 records 481.436 s in the native material stage; native work therefore reached the terminal binary check at roughly 3,188 s, leaving about 412 s of the outer envelope for terminal evidence processing. This is inferred from worker stage timers, not a new runtime claim.

The launch-time identity is preserved in `build/spin2014-g7-instrumented-20260921/spec.json`. This report and the diagnostic tool edits below changed three files in that identity afterward (`native_qp_replay.cpp` and the stage auditor under two labels). Accordingly, the historical G7 batch now fails a fresh identity check and must be treated as `STALE_RESULT` for any current-source acceptance. Its immutable traces remain useful diagnostic inputs; the new QP replay has its own explicit source, binary, and trace hashes.

### Direct-audit candidate run

- Run directory: `build/spin2014-g7-direct-audit-candidate-20260921/run/repeat-0/`.
- Started from the configured initial state. It completed 21 worker attempts: 20 durable stage boundaries plus the terminal stage. All 21 trace files were independently rehashed; total native trace size was `6,994,511,530` bytes.
- Outer elapsed time was `3,583.948 s`; peak memory was `8,546,852,864` bytes; accounted disk use was about `7.38 GB`. The job receipt records 24 GiB memory maximum and swap maximum 0.
- Solver status: `BINARY_MOMENT_VALIDATED`; normalized binary constraint residual `2.5173651121e-7`; normalized moment error `8.93e-15`. Header geometry state: `MATERIAL_BOUNDARY_VALIDATED`, with native geometry present. This is **not** the independent `SealedSpinArtifact` geometry/topology/wall-thickness audit, which remains `NOT_GRANTED`.
- The 21 passes reached level 7; terminal partition has 22,056 cells, of which 21,696 are level 7. The configuration remains G7/O7, not an author-resolution level-9 result.
- Original result: `build/spin2014-g7-direct-audit-candidate-20260921/run/results.json`. Trace receipt SHA-256: `2d872452a0595a4b274e06900aba3398c3ec54194393b572824a142d68819ae9`; stage-audit SHA-256: `4beb574ef5c8856745e76922a3669ed47e509bb5d65cb4ffea920158aee4f33b`.
- Offline verification recomputed all raw trace hashes and accepted the receipt (`VERIFIED`). The original report verifier had a JSON key-type bug in two display-only level-7 summary fields; an unbound post-hoc verifier normalizes those keys and independently reports `level7_partition_entered=true` and `level7_trace_recorded=true`. Its SHA-256 is `14ee451657d50d6b5e02569f16cd571eec94e250d1cb1ea12597789211d9005f`. No bound run file was changed after completion.
- This is a single pilot: determinism is `NOT_VERIFIED`. The original attempt and this candidate remain separate; the failed historical archive was not modified.

## Timing diagnosis

The following phase totals overlap and must not be added as if they were disjoint:

| Instrumented phase | Total |
|---|---:|
| Material-binding stages | 2,018.299 s |
| Adaptive SLQP solve | 1,273.598 s |
| Checkpoint decode | 568.314 s |
| Source octree build | 563.464 s |
| Reopening/auditing prior iteration traces during checkpoint decode | 558.166 s |
| Endpoint checkpoint restore | 178.544 s |
| Endpoint checkpoint capture | 177.970 s |

Each of the 20 complete paused workers only publishes a checkpoint. In the historical run, the terminal attempt expanded the full 21-pass history into JSONL and missed the deadline. In the direct-audit candidate, the terminal attempt traversed the same native records and emitted only compact stage-audit payloads, while recording SHA-256 for every raw trace; it completed inside the unchanged deadline. This demonstrates closure of the evidence pipeline for one G7 pilot, not a measured speedup for the optimizer itself. Split/merge lineage and moment checks show no conservation defect: the maximum observed per-transition relative moment discrepancy is approximately `3.2e-16`. No split/merge semantic change is justified by this run.

## Fixed-QP hotspot replay

The two slowest complete native stages were selected: one-based stages 19 and 20 (trace passes 18 and 19). A freshly built scanner selected the earliest checked QP attaining each maximum working-set size. The selected QP was then replayed twice with the current replay source and twice with a separately instrumented copy of the same QP kernel. Configuration was taken from the run: L-BFGS memory 8, `maximumFaces=-1` (one possible face per bounded variable), and 1,028 QP iterations.

| Stage | Trace records | Selected iteration | Peak faces | QP iterations | KKT stationarity | Equality residual | Result |
|---|---:|---:|---:|---:|---:|---:|---|
| 19 | 452 | 98 | 417 | 418 | `6.51e-18` | `9.90e-16` | independently verified |
| 20 | 313 | 270 | 414 | 415 | `4.12e-18` | `6.58e-16` | independently verified |

For both selected records, the replayed working set exactly equals the recorded one; the independent box, face, dual, and complementarity residuals are zero. Repeated scientific result fields are identical. These are fixed local QPs only: nonlinear convergence, binary material, geometry, and physical behavior were not checked by this replay.

The timed-kernel repetitions locate the dominant work inside the QP as follows:

- Stage 19: Schur solve `0.1575–0.1581 s`; candidate certification `0.0459–0.0460 s`; active-set work `0.0604–0.0617 s`.
- Stage 20: Schur solve `0.1562–0.1567 s`; candidate certification `0.1043–0.1047 s`; active-set work `0.0618–0.0623 s`.

Thus the QP is a measurable late-stage cost, but each hotspot replay is below half a second and its KKT evidence is clean. The QP is not the immediate reason the final campaign omitted its result receipt; the terminal expanded-trace pipeline is.

## Reproducible QP evidence

Raw traces:

- pass 18 SHA-256: `a23f8002185c3377f0898772d37694a0c7b78faf6c8aa188e00dfcb8432099f3`
- pass 19 SHA-256: `e03e241ef9dda1b43ae1b107a700541eb1b0914d902c24abadb7792bf244bf0c`

Diagnostic identities:

- Current replay source `experiments/teapot_spin2014/native_qp_replay.cpp`: `22bb306026b1c5d9206d1d3414698ad5f088ace138e91c2f01cc5b7026eee0e7`
- Current core archive `build/spin2014-g7-instrumented-20260921/native/libmake_it_stand_core.a`: `6720a34d0e631e413ebf76cd556e18249c3409164738ac9fa0ac85c4ce87210f`
- Current baseline replay executable: `bca4d5783a523f0d72ed916288943d6a167f08dd21fcb73a29642eefaa554fdf`
- Timed replay executable: `da40556ce981cfc347f795e929463048ef54fb30fcd693a1054f233dfe73efe9`
- QP replay outputs and timing stderr: `build/spin2014-g7-instrumented-20260921/qp-replay-current/`
- Partition lineage audit: `build/spin2014-g7-instrumented-20260921/g7-lineage-pass12-19.json`

The timing variant differs from `src/spin2014/SlqpTangentQp.cpp` only by RAII phase timers and counters; both source files and resulting executable are retained. The standalone replay test helper is `experiments/teapot_spin2014/test_selected_native_qp_replay.py`.

## Phase-2 preparation and test results

Added `audit_native_diagnostics`, a no-expanded-JSONL audit entry point which strictly parses each native event line, checks duplicate keys and non-finite values, independently enforces pass/iteration ordering, exhausts the stream, and runs the existing stage audit over the full events. The pybind iterator emits each pass and only its terminal iteration payload; it still reads and checks every SLQP record. A durable receipt binds the header, input/output/trace hashes, compact audit, and each referenced native trace file's path, size, and SHA-256. The receipt clearly does not grant geometry, solver, or paper acceptance.

Validation: the current suite `tests/python/test_stage_audit_v2.py`, `test_native_qp_replay_audit.py`, `test_terminal_trace_receipt.py`, and the full `test_source_spin_material_bridge.py`, plus the candidate campaign worker test, passed **72 tests** under the freshly compiled candidate binding. A separate worker smoke test exercises a tiny synthetic TOP objective through the real native extension, writes no expanded JSONL, and independently rehashes/verifies the resulting native trace receipt. Candidate binding provenance covers the recompiled binding object and headers, the frozen core archive and prior build inputs; **4,990** build paths were verified. The G7 spec contains **5,123** source/input/build identity entries; those fingerprints matched at launch and in offline verification. The real G7 QP replays above also passed the independent KKT checker.

The candidate implementation and all run evidence are isolated at `build/spin2014-g7-direct-audit-candidate-20260921/`; the earlier prepared-only spec is preserved at `run-preflight-v1/`. The completed run retains the 24 GiB / 3,600 s limits and original G7 scientific configuration. The old 5.69 GB partial JSONL remains untouched. The direct auditor eliminated that expanded copy, but the candidate still reread and hashed `6.99 GB` of binary traces; only one run is available, so no deterministic repeatability or optimizer-speedup claim follows.

## Identical-current-code repeat addendum (2026-09-21)

The earlier `NOT_VERIFIED_SINGLE_RUN` statements above describe the state when that report was written. A separate from-initial-state repeat was subsequently launched under the same scientific configuration and the same 24 GiB / swap-0 / 3,600 s outer and 3,540 s inner wall limits. It is preserved at `build/spin2014-g7-direct-audit-candidate-20260921/determinism-repeat-current-20260921/`; it is not folded into or substituted for the original candidate run.

- The repeat completed with a verified process receipt in `3,581.556148107 s`, peak memory `8,481,501,184 bytes`, swap peak `0`, and accounted disk `7,383,752,704 bytes`. It finished inside the fixed wall limit, with about `18.4 s` remaining.
- Both runs returned `BINARY_MOMENT_VALIDATED` after 21 worker stages. Their final `output_hash` is `b73fa0200b36fcc15741de219797c5ccef013b3b90f07a8375227da78ee83b20`; `trace_hash` is `a4aea4438689eac293898926603a672d26aa8036221036f8d51d1084ea4c4af2`; and `stage-audit.json` SHA-256 is `4beb574ef5c8856745e76922a3669ed47e509bb5d65cb4ffea920158aee4f33b` in both runs.
- All 21 per-stage `.slqptrace` archives were independently SHA-256 compared and match byte-for-byte. The final `outer.npz`, `material.npz`, and `cavity.npz` also match byte-for-byte (SHA-256 respectively `5622b4c16d0e467081efee982639b9f232aa3012c5548395df86de6963dfd825`, `8771887bc7b424dce3525a2840660e19b5d3c5fb17843393c924e3c66962ad60`, and `d52992a2070923da4f3b19589f7a6b80c9e4511ba8d73750353eb701550eba9b`).
- Both runs report the same normalized binary-constraint residual `2.51736511197818e-7`, normalized independent-moment error `8.931127442540148e-15`, and terminal source-cell distribution. The repeat completed in `3,416.413508363 s` of accumulated worker time versus `3,421.710218681 s` for the candidate (difference `-5.296710318 s`); this small timing difference is not a performance claim.
- Conclusion: deterministic equality is verified for these two executions of this **reconstructed G7/O7 configuration in this environment**. It does not upgrade the input to an author-original asset, establish level-9 resolution, prove continuous-domain feasibility, or grant paper reproduction.

The exact current-code geometry postprocess remains a separate, incomplete gate at `build/spin2014-g7-direct-audit-candidate-20260921/geometry-seal-current-20260921/`. Its independent mass-moment comparison is within `1e-6` and exact outer-to-carved clearance completed, but `same_component_feature` and `inter_cavity_web` thickness are `UNKNOWN`; overall thickness and manufacturing admission therefore remain blocked. The new deterministic repeat confirms that the sealed candidate geometry is the exact repeat geometry, but does not resolve those geometry checks. Do not run or label the object as a qualified physical/PyBullet artifact until that gate is resolved.

## Current checker clarification (2026-09-21)

The paragraph above records the original seal report without rewriting it. A narrowly scoped checker correction now uses the audited material topology: with exactly one bounded void, an *inter-cavity* web does not exist and is reported `NOT_APPLICABLE`; with two or more bounded voids it remains `UNKNOWN` because there is no inter-cavity thickness checker. Regression coverage passes for no cavity, one cavity, and multiple cavities. Current source SHA-256 for `python/make_it_stand/verification/spin_artifact.py` is `fddd6a70e61963b23e888e64fecdf2182d90050b9f3713511ecf879ae1bcf048`.

This semantic correction does **not** verify the remaining `same_component_feature` thickness. Consequently `wall_thickness_audit.overall` remains `UNKNOWN` and manufacturing admission remains blocked; no PyBullet run or qualified physical artifact is authorized by this result. The old seal manifest/report remain immutable historical evidence and still show their original conservative `inter_cavity_web=UNKNOWN` field. No replacement sealed bundle was published: a fresh exact-rational seal would repeat the approximately 24.8-minute mass-integration pass while not changing the remaining acceptance outcome.

The new checker and regression tests were exercised in the targeted 50-test suite (`test_sealed_spin_artifact`, `test_spin_bundle`, and the `SealedSpinViewTests` subset): all passed. A source-bound postprocess attempt was retained but did not grant a new seal: one initial attempt stopped at a too-strict run-config metadata equality check; the corrected attempt reached artifact writing but exposed a missing output-directory setup; a third attempt was intentionally stopped at 462.6 s during exact-rational reintegration, because a complete prior pass took 1,478.9 s and this semantic-only correction cannot clear the remaining wall-thickness blocker. These are audit-runner outcomes, not optimizer failures, and none changes the original G7 or deterministic-repeat evidence.

## Superseding scope decision

The deterministic-repeat gate in the historical scope decision is now satisfied for two executions of reconstructed G7/O7 in this environment. The next blocker is the independent same-component feature-thickness audit; until it has a sound implemented checker and passes on the sealed candidate, keep manufacturing admission and physics `BLOCKED`. This reconstructed, level-7 case still cannot be counted as an original Spin-It case or as a level-9 paper experiment.

## Fresh source-bound geometry seal (2026-09-21)

A new postprocess job now supersedes the historical checker report above without touching either optimizer run or the old seal directory:
`build/spin2014-g7-direct-audit-candidate-20260921/geometry-seal-repeat-audit-20260921-retry3/`.

- The source-bound job exited successfully in `1,494.999586826 s` (24m 55.0s), peak memory `702,787,584` bytes, swap `0`; its input identity binds the historical G7 spec/config/result, both terminal outputs, 42 stage result receipts, the previous seal, reconstructed inputs, native binding, current checker modules, and this worker.
- The frozen candidate config matches all 5,123 identities from its spec. Both completed runs have 21 stages and identical output/trace/stage-audit hashes. The fresh job revalidated all 21 paired stage-receipt digests and trace `{sha256,size}` entries; bytewise hashes of all raw traces and final geometry were also compared in the preceding deterministic-repeat audit.
- The new bundle was written to `payload/sealed-artifact/`. `save_spin_artifact` revalidates before serialization and reloads/revalidates the persisted bundle before returning. Its scientific digest is `6cb3c587fb44d9aac1306e8f0110242c3eba9e408087c919d877922ef2265765`; manifest SHA-256 is `7aa3f73d13fd9157322ddb79a18aeefb7a5fb8ba91f17b6b6f5f304fad6d4bc7`; current `spin_artifact.py` checker SHA-256 is `fddd6a70e61963b23e888e64fecdf2182d90050b9f3713511ecf879ae1bcf048`.
- Recomputed normalized ten-moment comparison is `WITHIN_ENGINEERING_THRESHOLD` (tolerance `1e-6`; largest listed absolute component error is approximately `2.88e-15`). Topology reports one material component and one bounded void. Outer-to-newly-carved clearance is `VERIFIED`; inter-cavity web is `NOT_APPLICABLE`; `same_component_feature` and overall wall thickness remain `UNKNOWN`, so manufacturing admission is still `BLOCKED_UNKNOWN_WALL_THICKNESS`.
- This is a fresh `RECONSTRUCTED` G7/O7 geometry/material seal, not an author-original case, level-9 run, continuous-domain feasibility certificate, physical validation, or paper-reproduction completion. The source-job wrapper's generic `scientific_validation` remains `NOT_GRANTED`; the seal's scoped geometry and moment facts are the evidence reported above.

The authoritative files are [geometry-postprocess-report.json](</data/xvxinyao/moment-guided-spin-v1/code/make-it-balance-paper-sim/build/spin2014-g7-direct-audit-candidate-20260921/geometry-seal-repeat-audit-20260921-retry3/payload/geometry-postprocess-report.json>) and [sealed-artifact/manifest.json](</data/xvxinyao/moment-guided-spin-v1/code/make-it-balance-paper-sim/build/spin2014-g7-direct-audit-candidate-20260921/geometry-seal-repeat-audit-20260921-retry3/payload/sealed-artifact/manifest.json>). Do not start PyBullet or represent this as print-ready until same-component thickness has an implemented, tested, source-bound checker and the seal admits it.

## Scope decision

1. Do not modify objective, tolerance, line-search, adaptive split/merge, or binary-rounding semantics.
2. Keep this original run as an incomplete historical diagnostic; do not patch its missing receipts into completion.
3. Preserve this single-pilot receipt and run a second, independently prepared identical G7 experiment to establish deterministic scientific output; do not count repeated artifacts from the same run as a repeat.
4. Before further full-resolution work, integrate the independently audited sealed-material geometry checks and correct the candidate display-key normalization in a new verifier revision; the existing run's source identity must not be rewritten.
5. Only after the repeated G7 evidence and geometry gate pass should the work advance to higher resolution and official case matrices. This G7 pilot alone is not paper admission.
