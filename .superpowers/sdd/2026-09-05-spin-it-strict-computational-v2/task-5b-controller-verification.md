# Task5b controller verification observations

Implementation commit: `c624c845578ddb6c738caf1c9a4fca988ae0f843`.
This records actual controller checks, not a whole-paper completion grant.

## Observed final implementation checks

- Read the implementer's final report, distinct RED stages and final-2 logs. The first 190-pass full suite predates the final-snapshot metadata fix and is retained as prior-candidate evidence.
- Final-code full run records 190 actual `Passed` entries, zero failures, 311.74 seconds, peak RSS 4,253,072 KiB and exit zero. The controller counted entries and read completion/resource output; it did not unnecessarily repeat the whole suite.
- Ran `sha256sum -c task-5b-final-2-sha256.log` twice, including after the full suite: all 26 entries match. These bind changed source/header/test/doc files, build cache, core archive, three executables and actual Python extension; this scoped manifest is not a replacement for the final release's complete dependency attestation.
- Independently reran postcommit `ctest --test-dir build/strict-v2-spin -R '^spin-it-2014-(audited|slqp|givens|python)$' -V --output-log <task-5b-controller-focused-20260909.log>`: 4/4 pass in .14 seconds, including 9 Python unittest methods.
- Python tests recompute KKT components from actual JSON state, gradient, constraints, Jacobian and multiplier vectors with ordinary arithmetic, without invoking the native residual checker. The controller's fresh run exercises that independent recomputation and repeated diagnostic identity.
- Fixed restoration case reaches x=.75 with zero equality residual; inconsistent constraints retain L1=.5; the singular example remains RESTORATION_STALLED; 24 finite rejected trials and 24 invalid trials survive; a real capped witness LP returns UNKNOWN. Actual Givens test reports full KKT 6.3137399822221774e-11 and phi 1.5707963268031664.
- Read actual 2014 ON/reachability OFF linked executable evidence (1/1, 91.26 seconds). This checks real HiGHS linking rather than just creation of a static archive; the contract also runs within the final full suite.
- Legacy SLQP golden hash remains ba7371831d27b970aeab54b8c1c14f8b411023e16bd06c7fc0c8a9d27f9d9aae. BASE/current default-Givens solver comparison gives identical 13fcb7597f2bc6a0275921674dd5f109f63abeb0013487e5216bfba14aa99eb5 and exact beta/phi. This is the scoped solver comparison described in the report, not a full historical application rebuild or a universal byte-equivalence claim.
- Directly queried the loaded native module: `/data/xvxinyao/make-it-balance-paper-sim/build/strict-v2-spin/python/_make_it_stand.cpython-313-x86_64-linux-gnu.so`, SHA256 5998f82afe48a005ae5ad065f7a7744cc80a6433d3c66b432156bffb0aaf3c35.
- `git diff --check` succeeds and tracked `experiments/` changes from the task base are empty. Dirty main's tracked diff SHA is unchanged at 93e32af47ab1957a1efcb6716437656be129053f4269d2169378fcd95b05a35e. One attempted report read used the main directory and failed without changing files; it was repeated successfully in the isolated tree.

## Raw-log identities

- Controller focused log: 814da5cb15909fcdfc5fb40efdfada6f16e62edcbaae4c0a49fb5c7b67e4fa0d.
- Final-2 full raw log: 9f6eae400aa85da602777dfa82412b9723e23972bdbcc28ca702150002e368d6.
- Final-2 full timed console: d30aa0fe66f90ea2ec07892c32a2797fee9a66e23b559ddeb797f00788476e16.

## Evidence freshness and review

Both previous `build/task-5a-refresh-faster-section-6-6-run-{1,2}` bundles are now STALE_RESULT because their CMake content binding changed. They remain untouched. An unchanged-spec refresh will be recorded only after code review/freeze; no old bundle is silently relabelled. The earlier material candidate remains historical and already stale under prior native-module changes.

The initial review found a reciprocal-curvature overflow; the next section records its correction and final scoped approval. First-order KKT is not proof of a minimum, continuous-domain infeasibility, physical stability, full SLQP equivalence, or full 2014/Faster reproduction. Overall status remains PARTIAL_BLOCKED.

## Final fix1 verification at 166c526

- Read actual RED: finite curvature 9.9998886718268301e-321 was accepted, evicting valid history and yielding inf/NaN. The fix validates the reciprocal before any mutation. Fresh final controller execution reports rejection, retained reciprocal .125 and finite inverse action 1.
- Final fixed-code full run: 190/190 actual passes, 314.01 seconds, peak RSS 4,253,644 KiB, exit zero. All 26 entries in task-5b-fix1-sha256.log matched after that run. Earlier full logs and manifests remain historical and were not overwritten.
- Postcommit controller command uses the same four-test regex and new log `task-5b-fix1-controller-focused.log`: 4/4 pass in .13 seconds, including the nine Python methods and their independent KKT recomputation. Normal legacy/default Givens behavior remains covered.
- Independently imported the real extension again: same resolved build path, new SHA256 54d47b3f8ee834d4ceaac145d7491a692f07ae9038687a38103166a0d980519c. This matches the worker's loaded-module observation and scoped manifest.
- Controller focused log SHA256: 51eb6999a1a4c45414c5c1c08d665a0d0320e0666095189dcb47567e2ac72bce. Fixed full raw SHA256: d20d48d23328cd4cbc519883f89407618e4a2cda5173eb4b223c5e797dc65173. Fixed full console SHA256: 8b3265c1c073d0bd3011ebf4b3ea0a15d48f210d5a1af2713a639552153e9e03.
- Same fresh reviewer accepted the scoped fix with no new Critical/Important breakage. Initial CannotVerify source/runtime/main observations are covered above and by the fix checks. Minor invalid-witness active-count metadata and inherited CMake warning remain explicitly deferred in task-5b-independent-review.md; no whole-method perfection claim follows.
- Dirty main's tracked diff SHA remains 93e32af47ab1957a1efcb6716437656be129053f4269d2169378fcd95b05a35e. No merge, push, tag, release or original-asset mutation occurred.

## Refreshed unchanged Faster Section6.6 evidence

After fix review/code freeze, executed the existing runner sequentially twice, with the unchanged exact published specification and paper:

```
env PYTHONPATH=python /home/lin/anaconda3/bin/python3 experiments/spin_it_reproduction/strict/run_faster_section_6_6.py --spec experiments/spin_it_reproduction/sources/faster_section_6_6_spec.json --paper build/strict-paper-sources/sif-final.pdf --output-directory build/task-5b-refresh-faster-section-6-6-run-N
```

New directories (N=1,2) were checked absent before execution. Both independent freshness checks return CURRENT with no mismatches. Spec SHA256 remains 053f896c2b10ece89b136bfba24cd480bf44cd0da99c0661e6e9bd34495c15fa; paper SHA256 remains e7fc6395362ef7e63a703851176afd81d8e953bae4e8080279d7e244b3c51fcc.

Both scientific hashes equal 4a8e15c902415d85e381c9290f56398cd4f8fbd89d0753360d708e380d3a4754. `cmp` confirms result.json and report.md are separately byte-identical between runs. Result SHA256 is 4f6e763cb1b681babe4dfddfb2dcae1743a7493897a70234a6ba9e1ce2ab72fb; report SHA256 is 4a1cd2efc0febb98974b5d70a2a6e9961f3dea2527f6587b77c9a503856f5d98. All are unchanged scientific outputs, not newly completed 3D paper cases. Old Section6.6 runs and the stale historical material candidate were preserved.
