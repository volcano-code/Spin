# Task5o controller independent verification

Baselinebf925db611195942a06d119862d61a36d383d49e. Sole implementer owns native binding/build; controller separately compiled an exporter only after worker confirmed archive ready and no writes. Controller files do not import project/native functions for independent checking. Exporter calls actual SourceSurfaceOctree/SourceAdaptiveDomain/full SourceCageBbw/SourceCageBinding APIs and emits actual arrays, geometry, stored moments and upstream nodal weights.

## Actual checks

Two synthetic geometries × both COUPLED_PARTITION and INDEPENDENT_NORMALIZED protocols; entire executable run twice. Fixture0 retains original cube(.3,.7), designcube(.34,.66), level5source, initial2/max5partition, coarse2/fine3cage with26handles. Fixture1 uses nonsymmetric outerbox(.31,.32,.3)–(.69,.65,.7), designbox(.35,.36,.38)–(.61,.6,.65), scale1.7 and translation(-.2,.15,.6). No model selection based on solver success; these are binding tests only.

The independent NumPy/stdlib checker verifies complete finest-cell coverage without overlaps, original outer bytes/order, one ancestor beta owner per primitive, analytic ten box moments per owner, nonuniform-beta moment sums, exact integer corner deduplication and canonical FMA positions. It separately searches all actual fixture grid cells for the closed lexicographic incident cell, derives eight trilinear coefficients, and multiplies upstream actual nodal fields by those coefficients to check every returned chunk/all26handles. It does not call project moment or interpolation implementations.

| Geometry | Leaves | Full primitives | Bound points | Handles | Maximum leaf moment error | Total ten-moment error |
|---|---:|---:|---:|---:|---:|---:|
| Fixed cube, each protocol | 8 | 496 | 997 | 26 | 1.084e-18 | 0 |
| Transformed nonsymmetric box, each protocol | 8 | 147 | 324 | 26 | 5.204e-18 | 6.939e-18 |

Max nonuniform-beta error6.939e-18; every weight discrepancy0 in binary64. Two complete exports byte-identical. Six altered record types (primitive owner, outer point, returned weight, leaf moment, embedding coefficient, integer corner key) rejected for each case,24 actual negative checks. These checks validate the checker, not release-gate forgery defenses.

## Retained bindings

Controller directory `build/task-5o-controller-jOBBE6`; compile.log is empty, compile exit0, each actual export exits0 with empty stderr; check-1.log exits0. Resource environment: ulimit-v25165824KiB for compile/native runs, single-thread OMP/BLAS/MKL; no project archive writes.

- Native header SHA993887d711560770368331d314156ceaf57fd6162d2161a865a22de143aad991.
- Native source SHA69f09a52a0536bccd7903be8a592aad921a5588c01e7edbcf82d32bc6f5793fd.
- Linked archive SHA4c7b75d1620190594c06e9bee10bb0abf37cf8e6f6edd91865b7f84142404d43.
- Export executable SHA8e36c2bb77e54d2b1ac865ef6639284750777a3fb8fb6295413f72764bd60177.
- Each export SHAa26451a8e3f32644e86bb114922837338fa4edc72bf8e4b451ca5837a9b9b00a.
- Check log SHA7a5d8a2bcfb6a5ad1e0be867db9887e9380817b875070785941f504dacf97a78.

Final source/header/archive hashes were rechecked after implementationdd94566929fc3beb8e5a843a8fa773f43f186ffe and reporte7ac737: all unchanged from the independently executed probe. Final regression source d92cf6ec45a16aeabbb76aa36c3faf34f81a65d3c631d1070b96d6d8f62293ae, executable128306b8bf6b182f51533750acad3411cdbcc4d66022b0edc38ffbf5f181b560 and documentation1a918967ed136cbf29cdcf38a72bcbaadbf4b75d4186ee7fab0923049e22c186 matched the worker report. Controller matched33uniquePassedrecords toactualCTestselection; logSHA4baef868d43fe8ea47b15b3d0d38447818f3f02be647df76547f286cd0cc944c. Actual exact-OFF retainedbuildgraph SHAcec2bb7e4ce179404d0badcdb87cc1a4e9246c0ece7d4979387db6a4e33b49b4 matchesreport; inspectedcache2014ON/exactOFF and absentSourceCageBinding confirmed. MaintrackedSHA93e32af47ab1957a1efcb6716437656be129053f4269d2169378fcd95b05a35e unchanged.

Fresh taskreview SpecPASS/Approved, noCritical/Important findings; existingdeprecationwarning deferred. See task-5o-review.md. Task5o accepted for rest binding only. Any later changedproductioncontent requires revalidation/binding; no stale scientific result is silently upgraded.

## Explicit next geometry gate

`task-5o-controller-conformity.py` is an exact Fraction counterexample, NOT a project solver experiment: injective map(x,y,z)->(x,y,z+xy/16), Jacobian determinant1, maps a coarse triangulated shared face's center toz=1/32, while a newly evaluated fine-face center maps toz=1/64. Thus rest-point binding alone cannot certify a deformed conforming partition. Future joint material requires shared-face topology and actual geometry checks, not a false assumption that positive local mapping guarantees matching independently triangulated faces. No current SourceCageBinding deformation bug is claimed because this adapter intentionally does not deform or seal materials.

All cases here remain SYNTHETIC_REGRESSION. No source-jointSLQP, complete paper regularizer, original model, actuallevel9optimization, CAD, physics, or VERIFIED_COMPUTATIONAL_REPRODUCTION claim. Whole planPARTIAL_BLOCKED.
