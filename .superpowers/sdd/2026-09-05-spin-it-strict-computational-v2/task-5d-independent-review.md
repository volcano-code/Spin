# Task5d independent review and disposition

Review1 at45e47eac6cbbee481c1cad162fd997a77ca56e9a, reviewer `/root/review_spin2014_source_adaptive_solver`, full BASE777c8f7..HEAD package. Verdict Needs fixes, two Important, no Critical. No whole-branch approval.

## Reviewer report (verbatim)

### Spec Compliance

❌ Issues found: empty-partition handling and exact source-cuboid bounds do not fully satisfy Task5d. The remaining reviewed architecture, adaptation, fixed regression, and limited binary-status requirements are implemented.

⚠️ Pre-existing oracle correctness, full audited-optimizer correctness, and historical-result preservation are outside this task diff. Existing validation supports them; this was not a whole-branch review.

### Strengths

- `SourceAdaptiveDomain.cpp:40`: coverage is checked structurally, including nonempty branches and overlapping ancestor/descendant keys.
- `SourceAdaptiveDomain.cpp:96`: recursive merging requires identical beta values and oracle-confirmed child coverage, preserving fractions without averaging or snapping.
- `AdaptiveFillFunctional.cpp:14`: the new constructor preserves original outer moments and uses actual source contributions with density.
- `SourceAdaptiveSpinSolver.cpp:78`: adaptation is transactional, and changed partitions require another audited solve before completion.
- `SourceBinaryMoments.cpp:35`: independent integration streams original triangles and outward cuboid triangles through the tetrahedron path.
- `tests/spin_it_2014_source_adaptive_regression.cpp:96`: the fixed counterexample executes the real solver and preserves its binary constraint failure.

### Issues

#### Critical

None found.

#### Important

1. **Empty partitions can receive binary validation for angular objectives.**
   `src/spin2014/SourceAdaptiveSpinSolver.cpp:58`, `:84`

   The runner relies on solver rejection to handle empty beta. For YoYo/top objectives, the shared solver adds phi, producing a nonempty one-variable state even when there are no optimization leaves. Neither the evaluator nor KKT code rejects zero bounded variables. A centered symmetric outer can therefore converge, satisfy the unchanged-partition gate, and reach `BINARY_MOMENT_VALIDATED`.

   This violates the explicit requirement that empty cases remain UNKNOWN. Add an objective-independent empty-partition gate and regressions covering angular objectives. The existing empty-source test at `tests/spin_it_2014_source_adaptive_regression.cpp:162` exercises only StaticLowering.

2. **Emitted cuboid bounds differ from the source geometry’s canonical bounds.**
   `src/spin2014/SourceBinaryMoments.cpp:17`

   The emitter uses separate multiplication and addition; the source oracle uses fused multiply-add in `SourceSurfaceOctree.cpp:53`. For `origin=-0.3`, `side=1.4`, level 4, plane 7, the source boundary is exactly `0.3125`, while emission produces the preceding representable double. A selected cuboid can consequently extend outside the accepted source region, even though the moment discrepancy passes the numerical tolerance.

   Construct emitted bounds using the oracle’s canonical coordinate calculation. Add a shifted-root regression comparing emitted endpoints against accepted source bounds; current shifted-root tests cover queries but not emission.

#### Minor

- `src/spin2014/SourceAdaptiveSpinSolver.cpp:39` and `src/spin2014/SourceAdaptiveDomain.cpp:102`: long compound statements compress hashing and topology mutations into difficult-to-audit lines. Expand these into conventional blocks.
- `build/task-5d-all-target-build-01.log:1`: an existing libigl CMake deprecation warning remains. This is a nonblocking inherited warning, not a new compiler issue.

### Verification

- Read the complete supplied diff once; no files, index, HEAD, builds, or native tests were changed or executed.
- Confirmed the recorded `196/196` result at `build/task-5d-full-ctest-01.log:396` and inspected the fixed-case trace in `build/task-5d-green-focused-02.log`.
- Focused outside-diff checks covered nonfinite constraint rejection, audited trace/state handling, zero-bounded-variable KKT support, and canonical source bounds. The controller supplied the focused unchanged-evaluator check confirming no empty-beta rejection.
- Exact-rational arithmetic confirmed the shifted-root endpoint discrepancy.

### Assessment

**Task quality: Needs fixes.**

The principal integration is sound and meaningfully tested, but the empty angular case can grant a prohibited status, and cuboid emission must preserve the source’s actual geometric bounds.

## Controller adjudication

Both Important findings accepted for fix round1/5. Focused unchanged-evaluator inspection shows empty beta satisfies dimension/finite/bounds guards and angular objectives produce a phi column; no hidden rejection protects the source runner. Empty-source UNKNOWN is an explicit brief requirement; returning a numerical status for a zero-design phi-only problem would violate this contract even if the unchanged original object's moments satisfy constraints.

Controller independently inspected canonical source bounds versus emitter and reproduced the exact endpoint mismatch with Python `math.fma` and `Fraction`: canonical/once-rounded exact=0x1.4000000000000p-2, emitted=0x1.3ffffffffffffp-2, delta=-5.551115123125783e-17. A moment tolerance cannot substitute for consistent boundary geometry. Fix must preserve exact canonical source endpoints and maintain independent integration; no tolerance/fixture retuning.

Original implementer resumes sole source/build ownership for failing all-objective empty-source and non-dyadic shifted-root emission regressions, then scoped fixes and covering verification. Main does not patch production. Minors retained for later readability/final whole-branch work; no dependency upgrade/wholesale formatting in this fix round. Prior196/196 result remains preserved pre-fix evidence, not proof these newly identified cases passed. All project gates remain PARTIAL_BLOCKED.

Controller read actual fixround1 RED logs before production-fix report: `build/task-5d-fix1-red-empty-01.log` shows StaticLowering UNKNOWN but all four angular objectives BINARY_MOMENT_VALIDATED with binary_check=1, then assertion failure. `build/task-5d-fix1-red-endpoints-01.log` reports512cuboids (224stored/288virtual) and1664endpoint/literal-plane mismatches, then assertion failure. These are real behavior reproductions, not inferred-only bugs. Fix/final covering evidence and scoped re-review remain pending.

## Scoped re-review at6102533 (verbatim)

- **Empty partitions can receive binary validation for angular objectives — ADDRESSED.** `src/spin2014/SourceAdaptiveSpinSolver.cpp:59` now rejects empty partitions before constructing or calling the solver, independent of objective or phi. Tests at `tests/spin_it_2014_source_adaptive_regression.cpp:166` cover all five objectives and require UNKNOWN, preserved state, unavailable KKT, and absent binary artifacts.

- **Emitted cuboid bounds differ from the source geometry’s canonical bounds — ADDRESSED.** `src/spin2014/SourceBinaryMoments.cpp:17` now obtains endpoints from `SourceSurfaceOctree::cellBounds`. The accessor at `src/spin2014/SourceSurfaceOctree.cpp:259` validates keys and reuses the original canonical FMA calculation. Tests at `tests/spin_it_2014_source_adaptive_regression.cpp:191` check exact endpoints across stored and virtual FULL cells, literal regression anchors, and independent tetrahedron moment agreement.

### New Breakage in the Fix Diff

None found. Extracted query-key validation preserves existing checks; the independent moment checker and shared solver are unchanged.

### Verification

- Recorded RED evidence reproduces both defects: four angular objectives granted binary validation, and 1,664 endpoint mismatches occurred.
- `build/task-5d-fix1-green-focused-01.log:19` confirms all five objectives return UNKNOWN; line 24 confirms zero endpoint mismatches across 512 cuboids.
- `build/task-5d-fix1-full-ctest-01.log:396` records 196/196 passing.
- Read the complete fix diff once. No edits, builds, or native test reruns.

### Out-of-Scope Observations

None newly identified. Previous nonblocking readability and inherited deprecation-warning observations remain deferred.

### Verdict

**Fix round: All findings addressed, no new Critical/Important breakage. Approved.**

## Final controller disposition

Task5d bounded increment complete at6102533 after fixround1/5; independent review approved, allImportant/Critical closed. Actual final full196/196 and fresh controller4/4, hash checks, failure provenance and limitations recorded in task-5d-controller-verification.md. Prior Task5c shifted-root query coverage minor is addressed by this increment's query checks; Task5d adds exact emitted endpoint checks. Open nonblocking readability/CMake/Boost observations and priorTask5b invalid-witness activeVariables trace minor remain visible for later cleanup/whole-branch review.

Cannot-verify-from-diff items are not elevated to proof: existing oracle/shared-optimizer tests and prior task reviews support the selected contracts, not universal numerical correctness. Historical experiment directory unchanged and original main diffSHA preserved by actual controller checks. No original experiment, fulllevel9optimizer, fullgeometryseal, physical validation, completeM3/M4 or publication claim is granted. OverallPARTIAL_BLOCKED.
