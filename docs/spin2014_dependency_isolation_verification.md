# Dependency isolation verification — 2026-09-11

Implementation commit `596cde5` fixes the Spin-It 2014 dependency contract's
destructive source-tree scratch-directory reuse. Runs now use exclusive retained
directories below the invoking build tree; live logs and running-phase status
remain available after interruption. No solver equations or tolerances changed.

Verification, with distinct source snapshots preserved:

| Check | Actual outcome |
|---|---|
| Final real dependency configure/build/audited runtime | 1/1 PASS, 89.25 s |
| Final retention/concurrency/failure/live-log behavior | 1/1 PASS, 6.32 s |
| 2014 label suite before live-log-only correction | 29/29 PASS, 217.40 s |
| Final three unlabelled 2014 integration/audit tests | 3/3 PASS, 0.02 s |
| Fresh independent task review | Approved, no blocking findings |

These are not one final all-project suite. Existing unrelated contract scripts
and unbuilt targets were not run. Source and retained result hashes were checked;
the main project's tracked changes and the six unrelated Faster files remain
unchanged. Historical interrupted runs are preserved, not reclassified as passes.

Read-only source-binding preflight also confirmed why the next deformation adapter
must preserve actual material subcells: a synthetic 8-leaf partition contains 496
Full subcells, with volume 0.030517578125. Substituting leaf bounding boxes would
yield 0.125. Independent per-leaf ten-moment sums agree within 1.09e-18, and two
runs have identical output hashes. This is a synthetic correctness check, not an
original paper model, optimized geometry or physical validation.

Overall status remains `PARTIAL_BLOCKED`. Source material/full-BBW binding, the
complete joint optimization/regularizer/geometry chain, and original-paper
experiments still require their own verification.
