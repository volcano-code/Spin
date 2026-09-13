# Task5aa controller verification — boundary-component separation

Status: bounded Task5aa verified; no full-thickness, seal or paper grant.
BASE: `69717a182d32b85812b0d2363d762c3bd26293d3`.
Whole project remains **PARTIAL_BLOCKED**.

## Scope

Actual immutable binary material → Task5z exact material regions → every
distinct final boundary-component pair. Original component maps, orientations,
exact distance witnesses, the original material and integrated moments remain
bound together. This is `BOUNDARY_COMPONENT_SEPARATION` only. It does not establish
same-component feature thickness, full manufacturing wall/web thickness,
continuous BBW deformation validity, deformation sealing or physical stability.

Seven actual controller calls use the fixed protocol in
`task-5aa-controller-protocol.md`. These are `KNOWN_SOLUTION_VALIDATION`, not
original-paper cases or evidence of nontrivial optimizer convergence. The
independent Python checker imports no native solver/validator. It reconstructs
face adjacency and complete component maps, checker orientation and input hashes,
actual BBW positions, rational tetrahedral moments, all component-pair distances,
local and parent witness faces and the exact global minimum.

## Preserved development attempts

- The native synthetic wide-outer fixture stopped at upstream BBW before reaching
  the new gate. Its raw log and parameters are preserved by the implementer;
  a 300s BBW limit was initially a hypothesis, not an infeasibility certificate.
- A test-only diagnostic expression failed compilation (`require` expects
  `const char*`). An accidentally launched older regression executable was
  terminated and is not counted as validation. Final build4 successfully
  relinked the updated core and test before further verification.
- Private exporter `build/task-5aa-controller-build-nZX38X/export` linked an
  earlier archive before that relink, and its source later lost one EOF blank
  line. It was never run or frozen and is superseded by a fresh exporter build.
- Initial focused old baseline `spin-it-2014-source-conforming-boundary` passed
  1/1 in0.60s. Independent analytic oracle self-tests passed. Neither is the
  new gate's final acceptance.

## Remaining scope inventory

`docs/strict_remaining_inventory_2026-09-13.md` and the source-derived
`task-5aa-catalogue-inventory.json` distinguish 88 current catalogue rows
(86 computational,2 physical out of scope),40 overlapping original-experiment
rows and8 remaining engineering work packages. The current catalogue is not a
proof that all paper details have exactly88 units. The M0/M1 release aggregator
still lacks integrated positive requirement checkers; its zero verified count
must not be misrepresented as zero algorithm implementation.

## Actual repeated computation

Production commit: `98e48bda5a7c22d322165799c2c5d9f55161dc99`, only five
authorized production/test/docs/CMake paths. Final exporter:
`build/task-5aa-controller-build-g2Y7y9/export`.

Independent seven-call check passed, with maximum checked normalized moment
error `4.6725865421339353e-16`. All29 actual-output corruptions were rejected.
Two fresh processes in `build/task-5aa-controller-replay-282Ilc` produced the
same complete scientific JSONL SHA-256:

`365ae1ecec6f45b269d5b8225a1a52d817c70d3ff16cba1fa6f3681246cef68d`.

Actual elapsed times44.74/46.40s, maximum RSS32572/32916KiB, both exit0.
No timestamps or performance measurements enter that scientific digest.
The native regression additionally measured a finer actual two-cavity negative
control where both outer gaps pass but the cavity-cavity gap rejects; inclusive
exact threshold equality and one-ULP-above rejection also passed. That is a
separate synthetic engineering test, not an additional original-paper case.

## Provenance

`task-5aa-inputs.sha256` contains3646 frozen entries:3425 actual Ninja core
dependency paths plus configuration, linked products, tools, checker/protocol
and170 observed Python-module paths and their ELF dependencies. Input checks
passed before and after the final repeats. The exporter build's actual Git
context predates the scoped commit; its content manifest binds the unchanged
compiled production bytes, not a retroactively invented build commit. Reports
and results are not effective source inputs and do not create circular hashes.

The prior Task5z manifest now fails exactly CMakeLists.txt, build.ninja and the
core archive. It is `STALE_RESULT` for a current-source replay, not evidence of
historical geometry failure. The raw expected-mismatch log is preserved in
`task-5aa-prior-result-freshness.log`; old results were not regenerated or edited.

| Loaded product | SHA-256 |
|---|---|
| Core archive | `544c4485e1d0383e8d66ffc4f87252f30f43ae14db36e5f8b56e6c1fe33cade4` |
| New native regression | `08cb6563833ad03a50e5b0686c44875a0374e79abb2e7c79052e1ae0072021b5` |
| Final controller exporter | `0776f3d6c87da9199d1438726960951eda0c580a795b18c3ae0d9684f7abda8b` |

## Review, broader regression and preservation

Final all-target build completed63link steps. The registered related batch
actually passed **48/48**,645.20s, including the new native gate85.25s and both
actual dependency build/isolation tests. No NotRun, skips or failures occurred.
The raw command is `task-5aa-related-tests.sh`; full output is retained in
`task-5aa-controller-related-tests.log`. Registration has48 distinct IDs.
Its41-entry actual-program/Python-extension manifest passed after the suite.
This build profile disables Faster and cannot establish correctness of the
protected foreign Faster work or imply whole-project all-required-group PASS.

Fresh scoped reviewer `review_spin2014_component_clearance` approved the
implementation with zero Critical/Important findings. Its raw verdict is in
`task-5aa-review.md`. Two Minor items are parked: readability of compressed
scientific code and preexisting libigl deprecation noise in early configure
logs. No frozen source is reformatted after execution. The unchanged native
region backend catches some allocator exceptions into untyped incomplete facts;
these can surface as `Unknown`, not `ResourceLimit`. This remains fail-closed,
but universal allocator-failure classification is **not established** and needs
a separately scoped backend change. It is not paper or feasibility authority.

The controller also checked the unchanged Python seal gate: outer canonical
bytes must match the original (`spin_artifact.py:85`) and deformation remains
`UNSUPPORTED_NOT_VERIFIED` (`:143`). File SHA matches BASE:
`3aaa4605d7e58f6e9a09c09de0bf47a5847702036557ed140ca0a1bfdf55cb68`.

Build Graph updated at production commit:453files,5965nodes,63620edges,
bash/c/cpp/python, zero update errors. Navigation only, not mathematical proof.
Main tracked diff SHA remains
`93e32af47ab1957a1efcb6716437656be129053f4269d2169378fcd95b05a35e`;
the four protected foreign Faster tracked files remain
`a2229f437e553cdd4c584dcd65959e5ac92eb2e574db0b0e07d296d35eeffd75`.
No final release, merge, push, tag, deletion or original-asset upgrade.

## Finalization and next step

Final source-input and test-program manifests are rechecked after the entire
suite. Archived preview/repeat outputs compare byte-for-byte with their actual
run directories. `task-5aa-evidence.sha256` indexes the scoped raw outputs,
reports, scripts, review package and remaining-work inventory, excluding itself
and the cumulative historical ledger. The main worktree and foreign changes
remain preserved; the entire isolated tree is not clean because foreign work
is intentionally retained.

Next: same-component thin-feature/full wall-web policy and validation,
continuous deformation legality and deformation-aware sealing; then actual
nontrivial source optimization and same-artifact simulations. The8 remaining
work packages in the inventory stay open. No original-paper experiment count
or physical-validation status increased in this turn.
