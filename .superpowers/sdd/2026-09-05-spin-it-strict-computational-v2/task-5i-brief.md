# Task5i: upstream source-octree numerical admission

This bounded correction continues approved Task5 / M3, not a new reproduction scope. Workdir `/data/xvxinyao/make-it-balance-paper-sim`, branch `spin-it-strict-computational-v2`, BASE `a22de6b69b0dd0d550ad57157fbd9bac834cac60`. Read this first; it is your requirements. No subagents.

## Binding constraints

- Preserve all historical scientific outputs; old synthetic runs do not prove paper reproduction. Do not modify the dirty main checkout.
- No failed solve is an infeasibility certificate. No reconstructed asset or small analytic test counts as an original-paper case.
- Use test-first development and fresh verification. Commit only explicitly scoped changes, using the already configured user identity.
- Compile with at most 4 jobs; heavy experiments sequentially, 24 GiB / 60 minutes per case batch. Never delete old files to make space.

## Required work

Investigate `SourceSurfaceOctree.cpp::validateConfig` spacing admission. Last task corrected the finer cage grid's equal-ULP guard, but upstream uses `step < spacing`; checking the first fine cell alone cannot prove every dyadic plane is distinct. Establish an actual native regression before changing production. Candidate root x-origin `nextafter(2^48,0)`, side1, maximumDepth4: FMA planes(2,3),(4,5),...,(14,15) coincide, while first cell is noncollapsed. Use genuinely closed valid source/design polyhedra; independently record plane order/count and distinguish constructor's late error from early invalid_config rejection. An empty/invalid mesh rejection or merely duplicated floating formula is not enough. Include sign/alignment boundary checks, accepted strictly wider spacing, and ordinary-output compatibility. Valid equal-ULP aligned inputs may be conservatively rejected, but document this numerical scope honestly; do not call them additional collapse bugs.

Once root cause is observed, adopt a sufficient fail-closed admission bound consistent with SourceCageGrid: finest step strictly greater than maximum root-endpoint ULP, finite and positive. Document conservative numerical rejection in the public header and diagnostic; no recentering, geometry edits, coarsening, exhaustive root scan, or tolerance change. Do not silently relabel numerical unsupported as physical infeasible. Keep accepted ordinary source/config/partition outputs unchanged unless a concrete necessary identity change is justified in report. Scope source/header/tests plus report; don't refactor unrelated code or alter legacy solver behavior. Check SourceSurfaceCage and SourceCageGrid integration using existing covering tests.

Keep all RED/GREEN log captures in a fresh unique build directory (tee outer output; nested CTest can replace default LastTest.log). You own production edits and project build/configure/test commands. Run focused tests during iteration; after final changes run relevant source octree/level9, cage/grid, adaptive/material integration coverage and FasterOFF coverage. Then one main all-target build and full CTest before commit, not after every change. Existing main build `build/strict-v2-spin`; discover actual FasterOFF directory in build and verify cache flags. Use -j3 by default, ulimit -v25165824 and timeout3600. Report exact test counts, commands, exit codes, log paths and SHA bindings. No downloads or dependency updates.

Write report `.superpowers/sdd/2026-09-05-spin-it-strict-computational-v2/task-5i-report.md`: root cause, actual RED and GREEN, files, self-review, hashes, scope limits and residual risks. Commit only source/header/test/report. Never commit cumulative progress or controller docs. Return only DONE/concerns, commits, test summary, report path. Independent review is controller-owned after report.
