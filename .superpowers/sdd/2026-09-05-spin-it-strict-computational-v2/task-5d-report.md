# Task 5d implementer report

Base: `777c8f77270b3793478c5222a0e9adbad319bdf4`, branch `spin-it-strict-computational-v2`, checkout `/data/xvxinyao/make-it-balance-paper-sim`.
Overall scientific program status remains **PARTIAL_BLOCKED**. This increment adds original-source adaptive fill and a deliberately limited numerical binary-cuboid check; it does not complete Task 5/M3 or any original experiment. All new numerical regressions are `SYNTHETIC_REGRESSION`, with zero additions to original experiment counts.

Implementation outcome: Task 5d bounded native increment implemented and verified; main all-target build passed and full registered CTest passed **196/196**, zero failures/NotRun, in **359.77 seconds**. Awaiting controller postcommit checks/fresh review; no integration beyond the scoped commit is performed here.

## Implementation

- `SourceAdaptiveDomain` owns a shared const original-surface oracle and a separate private optimization partition. It starts by recursively querying/pruning the root, and refines nonempty cells to the requested initial optimization level independently of source geometry depth. Keys, moments, levels, source tree identity and configured leaf budget bind the partition SHA.
- Split inherits each exact fractional beta; recursive merge requires exact common beta below/above epsilon and actual nonempty direct-child coverage, consulting the oracle for every missing child. It preserves near-bound values rather than snapping/averaging them. Candidate state and beta commit only after successful validation, adjacency construction, moment calculation and hashing. Rootward merging below initial level and split/merge cycle conservation are tested.
- Sparse face matching stores six dyadic face keys per leaf, then looks up opposite faces and ancestors. Positive-area coarse/fine interfaces give one unique undirected edge; corners/edges alone do not. Sparse Laplacian entries are unit weights. This graph describes bounding-cell adjacency, not proved material connectivity.
- A conditional `AdaptiveFillFunctional` overload consumes `SourceAdaptiveDomain` directly. Its base is density times original unapproximated outer surface moments; columns are density times actual oracle intersection moments. The legacy constructor/evaluation code, `AdaptiveVoxelDomain::source()`, lattice export, cage and Python contracts are unchanged. There is no occupancy/null-source adapter.
- `SourceAdaptiveSpinSolver` validates and hashes every request/control, rejects legacy algorithms, uses the existing `solveSlqpLbfgsPaper` implementation, and stores its full result, trace and KKT witness for each pass. It records partition/beta/phi before and after each pass, full transitions, fractional counts and level distributions. Any changed partition requires another actual solve. Each call resets history and carries phi; neither a changed Laplacian nor a passed caller flag can reuse an old KKT as completion.
- The runner returns explicit UNKNOWN/incomplete states for pass, leaf, output and numerical failures. Empty-variable or invalid initial numerical evaluation cannot grant a KKT or binary status: the actual solver is attempted, its unavailable evaluation is recorded, and the last valid state is retained.
- `SourceBinaryMoments` emits only the selected intersections with source FULL cuboids. It never removes the whole bounding optimization box when source accepted geometry covers only part of it. Output caps apply before publishing the list. Round-to-one is exactly beta >= 0.5 with no projection or reoptimization; pre-round beta and signed per-variable deltas remain in the result.
- The independent checker streams outward twelve-triangle box boundaries through `computeSurfaceMassMoments`, separately reintegrates original outer triangles through that signed-tetrahedron path, and subtracts the results. Its expected geometric integration never calls Algorithm 1, source cell moments or functional columns. Ten-moment errors use density times root-side powers 3/4/5, constraints use fixed first/second-moment powers 4/5, and central inertia is checked for PSD within 128 machine epsilons of fixed degree-five scale. No output-fitted scaling is used.
- Versioned deterministic hashes bind request, partition transitions, complete existing audited trace hashes, objective/regularizer decomposition, numerical output and configured guards. Diagnostic nonfinite/unavailable values are tagged explicitly, and measured timing/RSS is excluded.

## Fixed native counterexample (unchanged brief values)

Root origin `(-1,-1,0)`, side 2, source geometry depth 5. Original outer unit box x-shift is exactly `1/8192`. Design bounds are `[1/8,3/8] x [-1/8,1/8] x [5/8,7/8]`. Conservative accepted bounds are `[3/16,5/16] x [-1/16,1/16] x [11/16,13/16]`, volume `1/512`, first x moment `1/2048`.

Static native z-axis/origin; initial/max optimization levels 2/5; initial fill 0.5; epsilon 0.05; Laplacian, center and density weights 1; max iterations 300; equality tolerance `1e-10`; projected-gradient tolerance `1e-8`; all other audited defaults unchanged.

| Pass | Optimization level solved | Variables | Actual solver result | Adaptation |
|---|---:|---:|---|---|
| 0 | 2 | 2 | accepted step from `(0.5,0.5)` to `(0.25,0.25)`; converged | 2 split parents to 8 level-3 leaves |
| 1 | 3 | 8 | two accepted steps; converged | 8 actual parents split to 8 level-4 nonempty leaves |
| 2 | 4 | 8 | new audited solve, converged at inherited state | 8 actual parents split to 8 level-5 nonempty leaves |
| 3 | 5 | 8 | new audited solve, unchanged partition, converged | no split/merge; all 8 fractional variables at maximum optimization level |

Final fractions range from `0.24999284403976074` to `0.25000811027521025`. Final KKT infinity norm is `9.2934395364696739e-10`. The final partition is genuinely level 5, despite equal leaf counts at intermediate levels. All eight values round to zero, so the actual removed cuboid list is empty; the original shifted outer remains.

Outcome: **BINARY_VALIDATION_FAILED**. Independent raw first-x constraint residual is `0.00012207031249999306`; normalized residual is `7.6293945312495663e-06` (analytically `1/131072`), exceeding the fixed `1e-6` gate. Ten-moment normalized max error is `6.9388939039072284e-18`. Thus independent moments agree while binary paper constraints fail. This observed outcome agrees with, but was not replaced by, the controller's independent analytic/SciPy derivation. No fixture was retuned.

Final FasterOFF scientific hashes:

- Request: `cd904916aa65fb83fec62313b5302323027367d29d99c2ed7ba60973d34e0d80`
- Adaptive trace: `688accd3a09fc6e7865657c25c0234ec1b38c56b427269057015da35d6e82e43`
- Final partition: `23e1cf1eb46227b3ea171e954507aae6ebdc56fb3302261f1833e89d3332d5d8`
- Output: `8da963bb5f27c4f8c059ef429cf4510ebb624ed44c7019316628f65352a640f0`

## Tests and TDD chronology

All paths below are relative to this checkout. All failure logs remain preserved; no old log was overwritten.

1. Before production edits, wrote the full source-adaptive behavior test using the intended native API and a missing-API availability branch. Command: `c++ -std=c++17 -Iinclude tests/spin_it_2014_source_adaptive_regression.cpp -o build/task-5d-red-availability`, then ran it. Exit 1 with `RED: SourceAdaptiveSpinSolver source-native adaptive integration is unavailable`. Evidence: `build/task-5d-red-availability-01.log`. This is an API availability RED; it is not represented as a pre-implementation behavioral failure of already available code.
2. Initial linked build, `cmake --build build/task-5c-no-faster-20260909 --target spin-it-2014-source-adaptive-regression -j3`, failed because bundled Eigen lacks range-based vector iteration. Evidence: `build/task-5d-build-focused-01.log`. Fixed only the new indexed serialization loops; no dependency upgrade. Successful build: `build/task-5d-build-focused-02.log`.
3. First native behavioral run passed, including the fixed fractional trajectory: `build/task-5d-behavior-01.log`. Expanded cavity, independent tetrahedron, known symmetric full-carve, nonzero-phi, output-cap and nonconvergence tests also passed: `build/task-5d-build-focused-03.log`, `build/task-5d-behavior-02.log`.
4. During implementer self-review, added a real numerical behavior regression before fixing production: an empty conservative source caused `solveSlqpLbfgsPaper` to throw `initial fill field must lie in [0,1]` instead of retaining the valid state with explicit UNKNOWN. It failed with exit 1: `build/task-5d-red-numerical-failure-01.log`, build log `build/task-5d-build-focused-04.log`.
5. Fixed that runner boundary and tested both empty-variable and nonfinite-objective cases; both now return explicit UNKNOWN and no binary artifact/KKT grant. GREEN: `build/task-5d-build-focused-05.log`, `build/task-5d-green-focused-01.log`.
6. Added exhaustive control-hash mutation checks and exact raw/scaled fixed residual checks. Final FasterOFF linked focused GREEN: `build/task-5d-build-focused-06.log`, `build/task-5d-green-focused-02.log`.
7. Genuine depth-9 source/coarse-partition run: `timeout 3600 build/task-5c-no-faster-20260909/spin-it-2014-source-adaptive-regression --geometry9`. GREEN: `build/task-5d-green-geometry9-01.log`. Actual geometry nodes 2,097,993; optimization leaves 8; face keys 48; peak process RSS 351,580 KiB. This is not a complete level-9 optimization.
8. Main all-target build: `timeout 3600 cmake --build build/strict-v2-spin -j3`. Passed, `build/task-5d-all-target-build-01.log`, including the native extension and all registered executables. Existing CMake deprecation warnings remain; there are no new compiler warnings/errors in this successful build.
9. Full registered suite: `timeout 3600 ctest --test-dir build/strict-v2-spin -j1 --output-on-failure`, captured in `build/task-5d-full-ctest-01.log`. Passed **196/196**, zero failures/NotRun, in **359.77 seconds** on its single full run. Its `spin-2014-dependency-contract` passed in 93.33 seconds: it freshly configures, builds and runs `spin-it-2014-audited` with 2014 ON / exact OFF / Faster OFF / reachability OFF. That fixture's executable/cache hashes are preserved separately in `build/task-5d-exact-off-hashes-01.log`. This full run also passed legacy lattice golden checks, existing dual-material/cage APIs and Python contracts. No suite failure was suppressed or retried.

Heavy runs are sequential; compilation uses at most 3 jobs. Build/test batches use `ulimit -v 25165824` (24 GiB per-process address-space guard) and 3600-second external batch timeouts, except short focused build commands which finish promptly. Starting and finishing free space was approximately 1.4 TiB, above the 100 GiB start threshold and 20 GiB stop threshold. No observed batch approached the 24 GiB guard; the explicit measured geometry-9 resource test was below 0.4 GiB RSS. The process guard is not misrepresented as a measured aggregate cgroup memory cap.

## Implementer self-review

The complete SHA-256 source/executable/library/extension/cache manifest is preserved in `build/task-5d-source-bin-cache-hashes-01.log`, generated with `sha256sum` over all ten changed native/build/test files, both linked new regression executables, both core archives, both CMake caches and the main Python extension. Principal binary identities are:

| Artifact | SHA-256 |
|---|---|
| FasterOFF new regression executable | `08fc64fe41f7fe98134243033f1645f6738fee1215d893f86228d003277ccf34` |
| FasterOFF CMake cache | `e57a7b85de312d5a27cb6979c88b67c1b2806d43dcb8f884a4f711c2650b01e2` |
| Main new regression executable | `ffb6717dcafb6738b73aa36ad63195ab38a91faa4f3244cf3ed2e07e890987d1` |
| Main Python extension | `a018f6cc152844f56a264a809ce3cffdbfd3aa9117478c28cfbb0bf0a36a4b3a` |
| Main CMake cache | `2c3f5a72c82d0ba5d8b5d0cd8869d13d89ea8a21e0daa405e88242e1778e517c` |

After the full suite, both source/binary manifests were rechecked with `sha256sum -c`; all entries matched. RED/failure/GREEN/build/full-suite evidence hashes are in `build/task-5d-evidence-hashes-01.log`. The full-suite log hash is `0867273f0651fc639952066e3a4803dcbb0d414ba0a2544439d761fbab0e5c6b`. `git diff --check` passed before the scoped commit. Build logs are retained in the checkout's ignored `build/` evidence area, and are not represented as source code.

Checked the requirement boundary against production code and the actual tests. Structural source coverage cannot be replaced by nearly-equal mass sums: coverage uses dyadic keys and recursively queries nonempty branches; merge inspects all direct-child keys. Scalar values are never averaged or snapped during merge. Sorted leaf order is applied together with beta remapping. The sparse graph has no dense pairwise search or finest-face expansion; exhaustive pairwise geometry exists only in small tests. All mutation methods construct candidates and finish throwing work before committing.

The source functional never substitutes design or accepted source volume for original outer base moments. Analytic box expectations independently check all ten moments, all five paper objectives and gradient/Jacobian finite differences at `1e-4`, `1e-5`, `1e-6`; signed cavity and density cases independently check original material retention. Shifted non-default root queries check stored/virtual parent lineage and direct-child moment conservation. The zero-removal symmetric solve has the analytic beta-zero outcome; the centered upper carve has the analytic full-removal optimum. Both execute the audited solver. Nonzero initial phi survives two actual solves with new histories.

The output-cap test retains the exact continuous partition and beta but exposes no partial cuboid list. Independent tetrahedron integration agrees with hand-derived solid/cavity cases and rejects a wrong expected moment vector. Incomplete and numerical-failure paths cannot grant binary success. Repeated real runs have identical request, partition, event and output identities. Every numerical SLQP control, request weight/tolerance and resource guard has explicit request-hash mutation coverage.

No subagent/reviewer was spawned. The controller independently derived expected mathematics without editing implementation. Controller brief/math/notes are not staged. Commit uses the already configured `volcano <xuxinyao307@gmail.com>` identity. No main-checkout edits, historical experiment edits, merge, push, tag, release, dependency upgrade, Python API addition or CAD/manifold output are included.

## Remaining limitations

`BINARY_MOMENT_VALIDATED` means only numerical agreement of an emitted disjoint cuboid-list subtraction and paper equality residuals, with finite positive mass and a roundoff PSD inertia check. It is **not** SEALED, FEASIBLE_WITNESS, full manifold/self-intersection/wall-thickness validation, stability validation, a `SealedSpinArtifact`, paper success, recovered author trajectory, or an original-case result.

The original source oracle still excludes unresolved closed boundary cells conservatively. Source geometry and optimization depths are different quantities. Uniform adjacency is over bounding cells. Exact common-value recursive merging and mandatory post-merge resolve are explicit implementation defaults because paper terminal details are incomplete. All numerical bounds/guards/default controls are implementation choices, not published author parameters. Only single density is integrated; dual material, deformation/cage integration, Python exposure, CAD/manifold generation and original-case completion remain pending.

Empty-variable or nonfinite initial source evaluations currently terminate as explicit numerical UNKNOWN, since the shared audited solver has no zero-variable solve contract; no alternative dummy optimizer is substituted. The report does not claim convergence for those cases.

## Fix round 1/5 — accepted independent review findings

Fix base: `45e47eac6cbbee481c1cad162fd997a77ca56e9a`. The preceding report/evidence remains the historical pre-review record. Its empty-source claim was incomplete: an angular objective creates an unbounded phi coordinate even with zero beta variables, and the shared general solver can legitimately solve that phi-only problem. That behavior must not grant a source-adaptive binary status. The independent review also identified emitter multiply/add rounding that differed from the source oracle's canonical FMA endpoints. Both Important findings were accepted by the controller; this round changes only these two issues and their regression coverage.

### Reproduction and fixes

1. Wrote all-five-objective empty-source tests before editing production, using a centered symmetric original outer, a tiny design with zero conservative accepted cells, and nonzero initial phi for angular objectives. `--empty-source` failed at the review base: StaticLowering returned UNKNOWN, but YoYo and all three top objectives returned `BINARY_MOMENT_VALIDATED` with binary checks. The test requires UNKNOWN, no binary moment/list/check/fill, unchanged empty partition, carried phi, an invalid/unavailable KKT and an empty SLQP trace for every objective.
2. Wrote the shifted non-dyadic-root emitter test before editing production. Root origin `(-.3,-.2,-.1)`, side 1.4, source depth 5, optimization level 4. It checks every lower/upper endpoint on all axes against independently evaluated once-rounded dyadic planes and also anchors both level-4 x-plane-7 endpoint roles to the exact literal `0.3125`. The actual selected FULL list contains 512 cuboids: 224 stored cells and 288 virtual descendants. The RED run found 1,664 mismatches (including the additional literal checks), even though approximate moment tolerance could have accepted the discrepancy.
3. Added an objective-independent empty-partition guard before constructing/calling the shared paper solver. The existing source numerical-UNKNOWN recording path preserves the valid partition, beta and phi and records unavailable KKT/trace explicitly. This supersedes the earlier implementation's reliance on an empty general-solver state being rejected. No shared general SLQP semantics, objective definitions or callbacks were changed.
4. Added `SourceSurfaceOctree::cellBounds`, with the same query-key validation and the same existing canonical `bounds` FMA calculation used for geometry classification and source cell moments. The cuboid emitter obtains lower/upper endpoints from this method. The source FMA arithmetic itself is unchanged; the emitter no longer separately multiplies and adds. Stored and virtual accepted cells now emit exactly the geometry that the source oracle represents.
5. Kept the independent signed-tetrahedron checker unchanged. The canonical-endpoint regression additionally reintegrates the emitted list through that checker and confirms its ten-moment agreement. Exact endpoint identity is a separate gate and was not replaced by, or weakened to, the moment tolerance.

### Fix-round evidence

All commands ran in `/data/xvxinyao/make-it-balance-paper-sim`, using the existing caches, at most `-j3`, sequential heavy batches, `ulimit -v 25165824` and `timeout 3600`. Initial free space was approximately 1.4 TiB. No old evidence path was overwritten.

- Test-only build before production fixes: `build/task-5d-fix1-build-red-01.log`.
- Actual empty-source RED, exit 1: `build/task-5d-fix1-red-empty-01.log`, from `build/task-5c-no-faster-20260909/spin-it-2014-source-adaptive-regression --empty-source`.
- Actual canonical-endpoint RED, exit 1: `build/task-5d-fix1-red-endpoints-01.log`, from the same executable with `--canonical-cuboids`.
- Empty guard implemented and separately verified GREEN before the endpoint change: `build/task-5d-fix1-build-empty-01.log`, `build/task-5d-fix1-green-empty-01.log`.
- Endpoint source-boundary integration build: `build/task-5d-fix1-build-endpoints-01.log`.
- Final focused FasterOFF build and verbose CTest: `build/task-5d-fix1-build-focused-01.log`, `build/task-5d-fix1-green-focused-01.log`. Source adaptive (including binary/functionals), source geometry, and legacy fill functional all passed **3/3 in 2.77 seconds**. All five empty objectives were UNKNOWN, all 512 emitted cuboids had **zero exact endpoint mismatches**, and independent tetrahedron moments agreed.
- The fixed quarter-fill regression remains `BINARY_VALIDATION_FAILED` with four actual solves, raw residual `0.00012207031249999306`, normalized residual `7.6293945312495663e-06` and moment error `6.9388939039072284e-18`. Request, adaptive-trace, partition and output scientific hashes are unchanged from the historical fixed-case values recorded above. No fixture/tolerance retuning occurred.
- Final main all-target build passed: `build/task-5d-fix1-all-target-build-01.log`. The full registered suite passed **196/196**, zero failures/NotRun/Skipped, in **351.02 seconds**, on its single fix-round full run: `build/task-5d-fix1-full-ctest-01.log`. Both depth-9 source tests, legacy lattice/cage/Python tests and the exact-disabled linked fixture passed. There were no failed full-suite cases to suppress or retry.

The amended source, executable, core-library, extension and cache SHA-256 manifest is `build/task-5d-fix1-source-bin-cache-hashes-01.log`. It includes the changed source-oracle header/implementation, source runner/emitter, regression tests, unchanged CMake content, both native source-adaptive/source-geometry executables, both core archives, both caches and the main Python extension. Current principal identities are:

| Fix-round artifact | SHA-256 |
|---|---|
| FasterOFF adaptive executable | `85f8b94f635b92304ea8f919e6c6b5ca046c11080507b7489bf0507a6fcf84aa` |
| FasterOFF source-geometry executable | `0d0a036f82f12e72fca8722d88d490018e1be030ce084d2afefbda565d5f5bdb` |
| FasterOFF core library | `ff51c3c93c4df4fb27ad3b46c7257352a630286bf077f470a3f6421f328adaef` |
| Main adaptive executable | `28c9be2977564bda931959987130df00cd3e3e5ad9f07c6601f4cc6bee31d2b4` |
| Main core library | `71911fd761334402723ac04b8314d444ad366085d75c47a68157abc171c54772` |
| Main Python extension | `55768270bc1f0282fb468a6c88f2fa7cb10217ed08126e7163ab0d33b6eb56b9` |

The two existing cache hashes and CMake content hash remain unchanged. The earlier report's source/binary manifests identify the pre-fix build, not the amended binaries.

The full-suite exact-disabled 2014 dependency fixture freshly compiled and ran successfully in 90.69 seconds. Its current executable/cache hashes are retained in `build/task-5d-fix1-exact-off-hashes-01.log`; this fix did not alter shared general SLQP or the exact-disabled API.

Post-suite `sha256sum -c` verified all amended source/binary/cache manifest entries and both exact-disabled entries. The immutable fix-round RED/build/GREEN/full-suite evidence hash list is `build/task-5d-fix1-evidence-hashes-01.log`. Full-suite log SHA-256: `68887c883904896c9108562a1df4294a6491b801d9eab3bdaeab017278949ab9`. Final free space remained approximately 1.4 TiB. `git diff --check` passed. No broader test rerun or source change followed this green full run.

Implementer self-review confirmed that the new guard is independent of objective/phi and never calls the shared solver for an empty source partition; its UNKNOWN pass retains explicit invalid KKT and absent binary artifacts. The source bounds accessor validates the same keys as `cell()` and forwards to the original FMA bounds function without changing source classification, geometry moments or source hashes. The emitter alone consumes those endpoints; independent tetrahedron integration remains unchanged. Both accepted findings are addressed; fresh controller verification/re-review remains separate from this implementer disposition.

### Review disposition and deferred ledger

| Finding | Disposition |
|---|---|
| Important: angular phi-only empty source grants binary status | Fixed by the source-only objective-independent guard; real RED/GREEN across all five objectives. |
| Important: emitter endpoints differ from source canonical planes | Fixed using a source-owned validated canonical bounds method; real RED/GREEN covers lower/upper endpoints, all axes and stored/virtual FULL cells. |
| Minor: compressed hashing/topology statements | Deferred as directed to later readability/whole-branch work; no broad reformat in this round. |
| Minor: inherited libigl CMake deprecation warning | Deferred/nonblocking; no dependency upgrade. A pre-existing Boost bind deprecation pragma is also visible when rebuilding the source-oracle translation unit. |

The independent checker still does not certify manifoldness, self-intersections, wall thickness, stability, a sealed artifact or original-case paper success. Overall status remains **PARTIAL_BLOCKED**. No agents, controller-note edits, Section 6.6 refresh, historical experiment mutation, merge, push, tag or release is part of this fix round.
