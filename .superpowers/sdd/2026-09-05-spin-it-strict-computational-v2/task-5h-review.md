# Task5h independent task review

Initial range59ed96c12e8f24582f14459056c05f74d9c15492..4c6f89b8ef5a72a47e4ca6ac872ecf22dd122990. Fresh read-only reviewer `/root/review_spin2014_source_cage_grid`, brief/report/review-package supplied. No helpers or duplicate suite. Review applies production commit, not later controller docs.

Initial verdict: **SpecFAIL, qualityNeedsFixes**. One Important, no Critical.

Important finding, `SourceCageGrid.cpp:63`: `step < spacing` accepts equality, which does not guarantee strictly increasing FMA planes. For x-origin `nextafter(2^48,0)`, rootSide1, coarse3/fine4, step=ULP=.0625; represented fine planes(2,3),(4,5),…,(14,15) coincide although coarse planes are distinct. Thus fine cells can have zero width and distinct keys identical positions, defeating node identity and embedding. Reviewer independently reproduced the FMA numbers through read-only libm, inspected upstream root-validation semantics, but did not run the full source constructor. Original implementer must first obtain actual constructor regression RED, then fix representability without modifying original geometry or upstream legacy contracts.

Reviewed strengths: selected-only expansion and pre-count (`:79`), canonical complete identities (`:129`), bounded sparse lookup (`:164`), finite/budgeted generic interpolation (`:204`), 2014/exact conjunction independent of Faster (`CMakeLists.txt:262,496`).

Minor: existing CMake dependency deprecations, disclosed in report. Deferred separately; no dependency refresh.

Controller checks outside diff: main tracked diffSHA unchanged, no main mutation/old-output deletion; actual RED original tee logs read and independently exported connectivity/field checks completed. Those initial checks and201/201 full suite are **pre-fix evidence**, not proof against the newly identified numerical counterexample. Full build provenance closure and original experiments remain separate unfulfilled requirements.

Fix round1 dispatched to original implementer. Current task acceptance blocked until genuine regression/correction and scoped re-review. No initial approval inferred from passing ordinary-coordinate tests.

## Fix1 re-review — addressed

Fix range4c6f89b8ef5a72a47e4ca6ac872ecf22dd122990..da1bd34741c776275f420d274d4479f1291de00c. Original reviewer resumed for scoped fix review, supplied fix package/report/sole finding. No duplicate whole-task review, helpers or test reruns.

**Equal-ULP FMA plane collapse — ADDRESSED. SpecPASS / qualityApproved.** No new Critical/Important breakage; no new out-of-scope observation.

Reviewer evidence: strict `step > spacing` condition (`SourceCageGrid.cpp:66`) exceeds the maximum rounding-preimage width for either sign, including binade crossings. Header`:25` explicitly documents conservative rejection of otherwise valid aligned equal-ULP roots. Tests`:179` construct valid upstream cages outside the rejection catch, independently check coarse-plane order and collapse counts, and reject numerically unsupported grid construction. Actual positive-binade fixture establishes seven collapsed intervals before the fix; other two fixtures document conservative admission policy, not two additional collapse bugs.

Controller read actual main/FasterOFF4/4 covering-test logs and the appended fix report. Separate amended-binary probe again produces two byte-identical ordinary-coordinate outputs, with independent connectivity/field checks passing. Initial201/201 is retained as pre-fix historical evidence only. Controller final all-target/current-source full-suite run recorded separately after completion; re-review approval itself is not a test execution or full reproduction grant.

Deferred outside this fix: the upstream source-octree's related spacing guard remains a separate numerical-admission audit item; no claim that every existing geometry path has been fixed. Existing dependency deprecation warnings remain deferred. Dirty main and original-output preservation were controller checks, not inferred from the diff.

Final controller verification after fix/re-review: current-source all-target build exit0, full CTest201/201 passed482.35s, no duplicate/missing IDs or NotRun/Skipped markers in the outer capture. Production source/archive/executable digests unchanged. Logs and SHA bindings are in task-5h-controller-verification.md. This is execution evidence separate from the review judgment.
