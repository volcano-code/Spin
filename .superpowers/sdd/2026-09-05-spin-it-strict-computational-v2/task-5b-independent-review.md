# Task5b independent review and disposition

Reviewed base080b1f1d1f275cc873b1ef9317af0e0153e85d50..c624c845578ddb6c738caf1c9a4fca988ae0f843 by fresh read-only `/root/review_spin2014_kkt_restoration`.

Initial verdict: specification issues; quality Needs fixes. Critical: none.

Final scoped disposition: fix round1 is accepted at `166c526958a4d07326a5d69b3aaf74b47f2289c5`. The Important finding is addressed, with no new Critical/Important breakage. Task5b is accepted within its bounded first-order scope; the two deferred minor findings below remain visible and full Task5/M3 is not complete.

## Important finding (verbatim)

`src/spin2014/SlqpHistory.cpp:6–10`: `update()` checks finite positive curvature but stores `1./curvature` without checking that reciprocal. For finite one-dimensional `step=gradientChange=1e-160`, curvature is approximately `1e-320`, the cutoff underflows to zero, and the accepted inverse curvature becomes infinity. Subsequent two-loop recursion can generate NaN/Inf from an accepted history entry. Compute and validate the reciprocal before evicting or appending history; add this boundary to the actual history test.

Controller disposition: confirmed from current source. Reciprocal validation before any mutation is required by the explicit finite-history contract; it is not a changed model, tolerance, or optimization objective. Fix round1 will resume the original implementer. Normal legacy golden identities must remain verified; the invalid reciprocal cannot be retained for compatibility.

## Deferred minor findings (verbatim)

- `src/spin2014/SlqpLbfgsSolver.cpp:390–397`: An invalid witness returns before counting `activeVariables`. Such a row can contain nonzero active membership but report zero active variables; `finish()` retains that same-state row. Populate the count before this early return.
- `.superpowers/sdd/2026-09-05-spin-it-strict-computational-v2/task-5b-build-core.log:1`: The retained build output contains the libigl CMake compatibility warning. This is inherited dependency noise, not a reason to expand Task5b into a dependency upgrade.

These remain visible for subsequent trace cleanup/final whole-branch review. Neither is an authorization to edit historical logs or upgrade dependencies in this fix.

## External-evidence disposition

The review cannot independently establish current source/binary/cache or dirty-main state. Controller observations in task-5b-controller-verification.md cover the current reviewed commit with 26 matching identities, actual loaded extension, fresh focused tests and unchanged main tracked diff. These checks will be refreshed for the fix.

No original-paper, high-resolution, physical-validation or full-method acceptance follows from this bounded review. That limitation is confirmed, not waived. Overall PARTIAL_BLOCKED remains.

Review method: complete diff/brief/report through EOF; focused helper/consumer checks for active selection and dual/cage compatibility; a read-only arithmetic probe confirmed the subnormal-curvature counterexample. No suite, native binary, build or file mutation by reviewer.

## Scoped fix-round review (verbatim verdicts)

**Reciprocal-curvature overflow — ADDRESSED.** `src/spin2014/SlqpHistory.cpp:9–12` checks the reciprocal before any eviction or insertion and stores the validated value.

**Regression — sufficient.** `tests/spin_it_2014_audited.cpp:75–94` exercises the fixed `s=y=1e-160` case, verifies its overflow precondition, and checks rejection preserves all history entries and the finite inverse action.

**Evidence checked:** `task-5b-fix1-red-console.log:20` records acceptance with `inf`/`NaN`; `task-5b-fix1-final-focused-console.log:20` records rejection with retained inverse `.125` and metric `1`. Read the complete fix package through EOF; no tests, builds or mutations performed.

**New breakage in fix diff:** None.

**Out-of-scope observations:** No new observations; previously deferred minor findings remain outside this round.

**Fix round:** All findings addressed, no new Critical/Important breakage.
