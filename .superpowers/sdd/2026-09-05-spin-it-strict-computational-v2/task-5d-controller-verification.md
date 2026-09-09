# Task5d controller verification

Implementation under review: `45e47eac6cbbee481c1cad162fd997a77ca56e9a`, base `777c8f77270b3793478c5222a0e9adbad319bdf4`. This is a bounded native source-adaptive increment, not full Task5/M3 or original-paper experiment completion. Overall status remains PARTIAL_BLOCKED.

## Evidence read and independently checked

- Read the entire implementer report. Independently counted196 actual PASS records in `build/task-5d-full-ctest-01.log`; its summary is196/196,0failures,359.77seconds, without NotRun/Skipped records. This was the implementer's one full run for this code, not a controller rerun of the whole suite.
- Actual `sha256sum -c build/task-5d-source-bin-cache-hashes-01.log`: all17 entries matched current files, including10source/build/test inputs, both new executables, both archives, both caches and the main Python extension.
- Actual `sha256sum -c build/task-5d-exact-off-hashes-01.log`: both linked exactOFF audited executable/cache entries matched. The full registered dependency-contract test is evidence of actual2014 ON/exact OFF/Faster OFF configure/build/run, not merely a configure flag assertion.
- Actual `sha256sum -c build/task-5d-evidence-hashes-01.log`: all7 preserved RED/build/GREEN/full-suite logs matched. Full-suite log SHA256 `0867273f0651fc639952066e3a4803dcbb0d414ba0a2544439d761fbab0e5c6b`.
- At handoff, `git status --short` was empty. The entire task diff has no changes under `experiments/`. Dirty main tracked-diff SHA256 still `93e32af47ab1957a1efcb6716437656be129053f4269d2169378fcd95b05a35e`. No merge/push/tag/release/main edits.

## Independent mathematical check

See `task-5d-controller-math.md`. A separately derived analytic8cuboid constrained problem and SciPy solve (no native solver/column reuse) predict fractional values close to.25, final objective about.24987788473640274, and binary first-moment residual exactly1/131072 under the fixed normalization. Native observed final objective.24987788473640282, own KKT9.2934395364696739e-10, normalized binary residual7.6293945312495663e-6 and ten-moment disagreement6.9388939039072284e-18 agree with that expectation. Trajectory identity across different solvers is not required or claimed. The native fixed start.5 is not a known solution.

This is an intentional example of continuous convergence followed by BINARY_VALIDATION_FAILED; it is neither an original paper case nor an infeasibility certificate or physical tipping claim.

## Postcommit verification and review

Controller fresh focused native rerun used24GiB per-process address-space guard, one test process at a time,3600second external timeout; complete output captured by CTest in `build/task-5d-controller-focused.log` and `/usr/bin/time` in `build/task-5d-controller-focused-resource.log`. Actual5/5passed in21.82seconds, exit0, peak352224KiB, no swaps. Tests: source-adaptive, source-adaptive-geometry9, source-surface-octree, fill-functional and audited.

Fresh main executable output has exactly the same fixed request/trace/output SHA values as the reported FasterOFF run. It independently re-observed the raw/scaled binary failure and four actual solves. Fresh geometry9 has2,097,993nodes,8coarse optimization leaves and48face keys. This does not imply full level9 optimization. The previous functional and audited/Givens regression cases also passed; this focused run is not a second entire196test suite.

Fresh independent read-only task reviewer `/root/review_spin2014_source_adaptive_solver` received complete79,675byte review package BASE..HEAD, binding brief and implementer report. Review1 returned Needs fixes for twoImportant findings (empty angular source status; noncanonical emitted bounds), noCritical. Controller confirmed both; full review and adjudication in `task-5d-independent-review.md`. Original implementer resumed fixround1/5. The passing runs above are preserved pre-fix evidence, not a claim those newly identified cases were covered or that the task passed review.

## Evidence freshness

Actual calls on both prior `build/task-5c-refresh-faster-section-6-6-run-{1,2}/evidence.json` return STALE_RESULT because CMake content binding changed. Their files are preserved unchanged. Two new Task5d runs will be generated after code/review freeze. Refreshing unchanged mathematical results is not newly implementing Faster or completing its original3D experiments.

## Fix round1 postcommit verification

Fix head `6102533779a8940c0c808def1aad2a76bff4b488`. Read the entire appended fix report. Independently counted196actualPASS in final `build/task-5d-fix1-full-ctest-01.log`,0failures/noNotRun/noSkipped,351.02seconds. Its SHA256 is68887c883904896c9108562a1df4294a6491b801d9eab3bdaeab017278949ab9. Actual verification of fix1 manifests matched all15source/binary/cache,2exactOFF and12evidence-file entries. Old manifests correctly describe pre-fix binaries; they were not overwritten or represented as current amended identities.

Controller freshly ran4focused native tests after fix commit: source-adaptive, source-surface-octree, fill-functional and audited. Actual4/4passed in2.87seconds,exit0; logs `build/task-5d-fix1-controller-focused.log` and `build/task-5d-fix1-controller-focused-resource.log`. Directly observed all5objectives on empty source UNKNOWN/binary_check=0; actual512cuboids (224stored/288virtual),0canonicalendpointmismatches; fixed four-pass counterexample still BINARY_VALIDATION_FAILED with identical request/trace/output SHA and raw/scaled residuals. No production changes followed this run.

The same reviewer received the fix-only45e47ea..6102533 package (31,376bytes), report appendix, brief and two findings. Scoped review approved both fixes, no newCritical/Important breakage. Fullverbatimverdict in task-5d-independent-review.md. This is not a whole-branch review or authorization to release.

## Final evidence refresh and bounded handoff

After code/review freeze, actually reran the unchanged FasterSection6.6 runner into new directories `build/task-5d-refresh-faster-section-6-6-run-1` and `-2`. Both independent freshness calls return CURRENT/no mismatches; oldTask5c runs remain STALE_RESULT and unchanged. Source specSHA053f896c2b10ece89b136bfba24cd480bf44cd0da99c0661e6e9bd34495c15fa and paperSHAe7fc6395362ef7e63a703851176afd81d8e953bae4e8080279d7e244b3c51fcc are unchanged.

Both fresh runs have scientific SHA4a8e15c902415d85e381c9290f56398cd4f8fbd89d0753360d708e380d3a4754, result-byte SHA4f6e763cb1b681babe4dfddfb2dcae1743a7493897a70234a6ba9e1ce2ab72fb, report-byte SHA4a1cd2efc0febb98974b5d70a2a6e9961f3dea2527f6587b77c9a503856f5d98. Their evidenceJSON digests differ normally because paths/runtime provenance differ: run1effd415a41e924f2291d9cb804e2543027fb43028f4cf69ff9879b8804218498; run28f9f592db93bc3c28bf83961ccb177f4025c550c1f34a7ba166be6991bcefe73. This is a provenance refresh, not new Faster capability/original3D case completion.

Final controller focused raw log SHA b07917a95b22d655af5e8fed2bbad22287eb2b3e777c15829a4a74351e177d4b; resource log SHA6b654795071e122420b8fbd57f16a62aee5f3e8eff84cc11f8bf0c5e788fcea3. Pre-fix controller log/resource remain d50810d3c13ee63cc064ec9f2e83a6ed20ce0d04e3cb4b183e91a557e9fc7604 / a2c53139350c97d6d56fc92a3d4c57b5daefff1a24eed45db7595195f1ddfd0e. All raw logs stay preserved in build/, not rewritten for presentation.

Task5d is complete within its brief; program statusPARTIAL_BLOCKED. Remaining work includes complete source-native final material geometry/sealing, fullcage/BBW/jointdeformation, dual-density products, genuine high-resolution optimization experiments, missing originalassets/criticalparameters, remainingFaster/2D/hydrostatic/CAD and full original-case/simulation matrix. No release/merge/tag/push/printing or originalmain writes. Final main tracked-diff SHA still93e32af47ab1957a1efcb6716437656be129053f4269d2169378fcd95b05a35e.
