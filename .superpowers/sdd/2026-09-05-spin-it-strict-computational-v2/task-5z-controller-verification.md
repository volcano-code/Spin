# Task5z controller verification — actual finite material regions

Scope: the bounded `SourceCageMaterialRegions` pre-seal layer only.
Whole-project status remains **PARTIAL_BLOCKED**. No new original-paper case,
nontrivial optimizer convergence, sealed deformation, full manufacturing thickness,
continuous BBW Jacobian, feasible witness or physical simulation is claimed.

Production commit: `73f52e2150db9d16fe6ac9559b911ec18e403fe8`, based on
`f2dfc71ec762d263f103db1a08ca53d50f3e6aed`, user-authorized identity volcano
<xuxinyao307@gmail.com>. Five scoped production/docs/test/CMake paths only.

## Actual computation and independent checking

Protocol is fixed in `task-5z-controller-protocol.md`. Four actual Task5y binary
inputs and six native region calls, with all26BBW handles and actual nonidentity
affine output. These are **KNOWN_SOLUTION_VALIDATION**, not paper experiments.

- Empty/full/tie inputs pass the bounded region gate at minimumOuterWall=.05.
- Full input at .13 is rejected by the measured actual outer clearance.
- Full input at Nef complexity1 returns RESOURCE_LIMIT without a product.
- Actual mixed-beta rounding failure is NOT_ADMITTED; no mutable facts grant.
- Native separate tests additionally recover two real disjoint cavities, verify
  exactly representable distance equality and nextafter-above rejection, retain
  original moments after both originating audits die, and exercise generic exact
  backend overlap/missing/nesting negatives. Those are not renamed controller cases.

Independent Python imports no native evaluator/checker. It rebuilds BBW point
interpolation, selected oriented triangle chains and compact primitive mesh
identities; recomputes signed rational tetrahedral moments and fixed-frame
constraints/inertia; checks containment for this convex fixed fixture; and
recomputes exact triangle-feature minimum distance and the actual witness pair.
It does not pretend its fixed-fixture convex check certifies arbitrary meshes.

Maximum checked normalized moment error: `4.3190471830795106e-16`.
All26 actual-export corruption negatives rejected, including false partition,
wrong primitive owner/key, mismatched mesh/evaluation, fabricated empty distance,
wrong clearance witness and publication of failed/resource-limited results.

Two new processes under24GiB virtual-memory/60min limits produced byte-identical
scientific JSONL (no time/resource fields in that digest):

`b6c39c8b7ef36de51fc166f03616270c3ccc84ab5d856675b3809f828716180a`

Actual runner directory: `build/task-5z-controller-replay-oI78iQ`.
Runtime30.06/29.38s; maximum RSS30940/30976KiB; both exit0. Independent checks
and negatives passed on final outputs; the prior preview is retained separately.

## Source, binary and dependency binding

Actual exporter: `build/task-5z-controller-build-5HBGSQ/export`.

| Loaded product | SHA-256 |
|---|---|
| Core archive | `267a1b9923ea5ded080565d13f69f5135daa6eacd369f6307496e5c965ec6bd7` |
| Region regression | `afbd3ef2323b29cc21ffc3b277d8c830dfc15c929a5d2587842e8e8f16ad2800` |
| Controller exporter | `b549e9338c777402ebfe954fc68d30f1e0b73622531278fb6545a0e385202ab1` |

`task-5z-inputs.sha256` records3642 entries including3423 actual Ninja core
dependency paths, compile configuration, link inputs, actual compiler/runtime,
runner/checker/protocol,170 observed Python-module paths and ELF dependencies.
The manifest passes before/after repeat and after the source commit. The exporter
was linked while final production bytes were frozen but the scoped commit was
being finalized; its build-time parent commit/diff record is retained, not changed
afterwards. The effective byte manifest, not a guessed commit-only provenance,
binds the actual loaded binary to final production. No reports/results enter the
source-input hash dependency closure.

The old Task5y frozen manifest now fails specifically for CMakeLists.txt,
build.ninja and core archive: **STALE_RESULT** for current-source replay. Its
historical run remains preserved. An initial lookup used the wrong old manifest
filename; that failure log remains separate and is NOT staleness evidence.
Initial evidence copy used `cp -n` and emitted portability warnings; later copies
use `--update=none`. This archival warning is not a scientific/test failure.
The staged source/script whitespace check reported an extra EOF blank line in
the already-frozen controller exporter. That original compiled input is retained
byte-for-byte rather than normalized after its repeats; this check is not claimed
clean. Production code/documentation checks are separate. No runtime/build error
or tolerance change is associated with the preserved blank line.

## Final test and review gate

Final registered2014-related batch: **47/47 actually passed**,563.34s, no
NotRun/skips/failures. All-target build completed62link steps. The new registered
region regression passed80.91s within that batch. Exact command/output is in
`task-5z-controller-related-tests.log`; registration inventory is
`task-5z-registered-tests.json`, including all47 distinct tests. The40entry
test-program/Python-extension manifest is checked again after the batch.
This profile has Faster disabled; no full-project or foreign Faster pass is implied.

Fresh reviewer `review_spin2014_binary_exact_regions` approved the scoped
production diff with no Critical/Important/Minor findings. It separately checked
Full primitive ownership and value-only allocation contracts. Its raw verdict is
archived in `task-5z-review.md`; no broad branch/merge review is claimed.
Controller checked the unchanged Python gate: `_check` requires outer canonical
bytes equal original and reports `deformation_validation=UNSUPPORTED_NOT_VERIFIED`
(`python/make_it_stand/verification/spin_artifact.py:85`, `:143`). Its file SHA
`3aaa4605d7e58f6e9a09c09de0bf47a5847702036557ed140ca0a1bfdf55cb68` matches BASE.
No whole-project pass implied.

## Build Graph and preservation

`task-5z-evidence.sha256` indexes the scoped archived raw outputs, scripts,
reports, review package and updated gap document; it deliberately excludes itself
and the cumulative working ledger. All input/test-program manifests were rechecked
after the complete suite, and archived output copies compare byte-for-byte.

Build Graph incrementally refreshed after production commit:442files,5880nodes,
62808edges, bash/c/cpp/python, zero update errors. It aids navigation only.
The header importer query returned0 despite known direct includes: the C++ path
resolution is incomplete, so this is not absence or correctness evidence.

Main tracked diff SHA remains
`93e32af47ab1957a1efcb6716437656be129053f4269d2169378fcd95b05a35e`;
the four foreign Faster tracked files remain
`a2229f437e553cdd4c584dcd65959e5ac92eb2e574db0b0e07d296d35eeffd75`.
No main edits, merge, push, tag, dependency download, deletion or foreign staging.
The whole worktree is not clean because protected foreign work remains present.

## Next boundary

Outer-to-carved separation does not prove inter-cavity webs or same-component
thin features. Full feature-thickness policy/validation and deformation-aware
sealing remain separate work; the existing Python hollow-only seal was not
weakened. Original case inputs/poses/parameters, nontrivial multiresolution
convergence,2014 dual-material/cutting, Faster all-family geometry/derivatives,
2D/hydrostatics/CAD, complete Table/Figure matrices and same-artifact independent
physics remain open. No complete-computational-reproduction declaration.
