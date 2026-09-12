# Task 5r implementer report

## Scope and implementation

Base `86234ff489451521f4fce4d14ac43677171a2276`, canonical branch
`spin-it-strict-computational-v2`. This implements the bounded actual-source
joint paper functional after Task 5q. Global `PARTIAL_BLOCKED` is unchanged.

New immutable `SourceCagePaperFunctional` owns the actual source moment snapshot,
an explicit fixed proper world-to-spin frame, all five existing Givens objective
choices, and independent inertia/center/fill weights. Its state is beta followed
by optional phi and all original `[scale,tx,ty,tz]` handles. Static has no phi
coordinate; its angle offset is `size_t::max()`. The actual fixtures have 8 beta,
26 handles, 112 static or 113 nonstatic coordinates. The initial builder uses
the source's actual identity-scale/zero-translation handle state.

A constant exact 10x10 affine raw-moment map transforms the actual world moments
and beta/handle Jacobians, including mixed tensor entries. No density is applied
again, and no moving center-of-mass recentering occurs. Existing eleven-coordinate
Givens mathematics supplies objectives/equalities and analytic partials; all
moment/phi partials are chained to the joint state. All world outputs and moved
points are retained. There is no joint Hessian claim.

Sparse actual fixed-rest partition Laplacian regularization contributes
`0.5*fillWeight*beta^T*L*beta` and `fillWeight*L*beta`. Its angle/handle derivatives
are zero. It is explicitly **IMPLEMENTATION_DEFAULT source unit adjacency**;
this is not the complete paper shape/displacement regularizer. Reconstruct the
functional after a partition or BBW change. No adaptive rebind, joint SLQP,
geometry-commit gate, inversion/thickness/intersection completion, or original
paper asset/result authority is claimed.

The default positive joint derivative budget is 10,000,000 doubles. Checked
combined allocation is `(11+K)*N+B` (retained raw/constraint Jacobians, objective
gradient, and fill scratch). Upstream budgets remain independent. No P×H,
P×4H or dense B×B global matrix is created. Value-only results contain no joint
or upstream global derivative arrays and work with insufficient derivative
caps. Existing Givens uses bounded eleven-coordinate local workspace.
Configuration/state checks reject malformed values, frames and dimensions;
arithmetic fails closed. Capacity/allocation errors carry `RESOURCE_LIMIT`.

SHA256 identities bind the immutable actual upstream snapshot/output, full state,
frame, objective, all weights, budgets, value/derivative mode and results. No
timestamp or pass flag enters identity. Production mathematics stayed unchanged
from the first successfully linked core through all independent numerical checks.

## TDD and numerical investigation

The complete missing-interface regression and four scoped test registration
lines were written before production. RED command (all commands from canonical
worktree, 24 GiB virtual-memory limit, maximum -j3, downloads already disabled):

```
ulimit -v 25165824
cmake --build build/task-5c-no-faster-20260909 \
  --target spin-it-2014-source-cage-paper-regression -j3
```

`task-5r-red.log` exited 1 with the expected missing public interface:
`fatal error: spin2014/SourceCagePaperFunctional.h: No such file or directory`.
The new implementation then compiled/linked successfully (`task-5r-green-build-1.log`).

The first ordinary double two-point FD run failed honestly, with baseline YoYo
column 94 objective derivative -0.000191306, h=3e-4 and relative error 1.3276e-5
(`task-5r-focused-1.log`). A diagnostic measured 1.82015e-6 at h=1e-4,
1.08118e-6 at h=3e-5, and 4.62155e-7 at an additional diagnostic h=1e-5.
That last probe was not substituted for the original steps. All original
three-step errors were then collected without early exit in
`task-5r-focused-all-fd.log`; the test still exited 1.

The **predeclared** rule was 1e-6 relative for |analytic derivative|>1e-7 and
2e-9 absolute otherwise, in fixed native fixture units. No threshold, fixture,
weight or original step was changed. Ordinary double FD had 41 finest-step
discrepancies across the cases. Largest original coarse relative discrepancy
was 1.80011e-4 (baseline TopInertiaOnly); largest finest relative discrepancy
was 1.80343e-5 (baseline StaticLowering). Every near-zero absolute check passed.
All original failed logs remain preserved; ordinary double FD is not called
passing at 1e-6.

Root-cause investigation used an independent complex-step point/tetra/inertia
oracle. It takes actual rest embeddings/final BBW fields, independently computes
point motion, transforms points, integrates signed tetrahedra, and rotates the
inertia tensor. It never consumes production source Jacobians, the production
frame map or production Givens. All entries met the unchanged threshold:
worst relative 6.02589e-11, worst absolute 1.28786e-14. The controller independently
corroborated all 10 cases and rejected 20 data corruptions (controller-owned
evidence remains disjoint). No production defect was evidenced.

Per the explicit numerical verification addendum, the independent reference
was then carried in long double throughout BBW interpolation, handle motion,
point transform, integration and objective arithmetic. At least 64 mantissa
bits are required and the final native run records the actual value. A five-point
centered stencil uses **all three original** h=3e-4, 1e-4, 3e-5. All entries
must still meet 1e-6 relative / the original 2e-9 near-zero rule; nothing is
selected from a successful subset. `task-5r-focused-longdouble.log` preserved
its exit 1 from the deliberately retained old double-FD final counter even
though all high-precision checks passed. Only after collecting this evidence
and the addendum ruling was that final counter made an explicit diagnostic.

All long-double five-point worst relative errors (columns include objective,
every equality and all ten spin raw derivatives):

| Fixture / objective | h=3e-4 | h=1e-4 | h=3e-5 |
| --- | ---: | ---: | ---: |
| Baseline YoYo | 3.59417e-11 | 7.55166e-11 | 2.31589e-10 |
| Baseline TopInertiaOnly | 2.42445e-9 | 4.74852e-9 | 5.01335e-9 |
| Baseline TopLoweringOnly | 7.44706e-11 | 2.9515e-10 | 3.55072e-10 |
| Baseline TopMassAndLowering | 1.82410e-10 | 7.44531e-10 | 6.85686e-9 |
| Baseline StaticLowering | 2.77424e-9 | 4.65326e-9 | 5.43886e-9 |
| Transformed YoYo | 3.04919e-11 | 2.14669e-10 | 8.88286e-10 |
| Transformed TopInertiaOnly | 1.09125e-9 | 7.47078e-10 | 5.00733e-9 |
| Transformed TopLoweringOnly | 3.04919e-11 | 2.14669e-10 | 8.88286e-10 |
| Transformed TopMassAndLowering | 3.04919e-11 | 2.14669e-10 | 8.88286e-10 |
| Transformed StaticLowering | 1.67985e-10 | 2.75862e-9 | 3.16705e-9 |

Worst absolute long-double FD error was 7.79168e-13. Per-step absolute maxima
and the 10 complex-step maxima are in the full raw logs. Across the 10 cases,
1,128 joint coordinates and 17,602 scalar derivative entries are checked by
complex step, and 52,806 entries by each three-step FD path.

GREEN command:

```
ulimit -v 25165824
cmake --build build/task-5c-no-faster-20260909 \
  --target spin-it-2014-source-cage-paper-regression -j3
timeout 3600 build/task-5c-no-faster-20260909/spin-it-2014-source-cage-paper-regression
```

`task-5r-final-focused-build.log` and `task-5r-final-focused.log` exited 0.
Final output explicitly says `ordinary double finest-step discrepancies=41
(retained diagnostic, not claimed passing FD)` followed by
`SYNTHETIC_REGRESSION PASS: independent LD5 and complex-step accuracy gates`.
The final suite enforces at least 64 long-double mantissa bits. The compiler
observation in `task-5r-exact-off.log` records `__LDBL_MANT_DIG__ 64`.

## Other regression coverage

Both actual Task5q source/BBW fixtures retain source depth 5, cage 2, grid 3,
nonuniform beta, distinct all-handle scales/translations, nonzero fixed origin,
proper rotation, and nonzero nonstatic phi. World and spin moments independently
match signed tetrahedra of actual moved/transformed points. Explicit inertia
tensor expectations test all five objectives and equality ordering. Actual edge
differences independently check fill energy/gradient, constant-beta nullspace,
and zero angle/handle fill derivatives. Density covariance distinguishes
invariant inertia/pure-height terms, density-squared mass-height, and unchanged
fill. Frame behavior is checked against independently transformed geometry,
including mixed entries; no moving recentering is used.

Value-only geometry/moments/objectives/constraints match the derivative path
exactly, with all global derivatives absent. Tests cover exact combined budget,
one-below rejection, an insufficient upstream derivative budget, invalid budget
tokens, wrong lengths, invalid beta/scale/phi/objective/weights, malformed and
overflowing frames, upstream and fill arithmetic overflow, deterministic output,
and state/frame/weight/budget/mode identity sensitivity.

## Full build, scoped suite and retained evidence

The full assigned build relinked 55 steps with `cmake --build
build/task-5c-no-faster-20260909 -j3` under the 24 GiB limit and exited 0
(`task-5r-full-build.log`). It remained Faster OFF, exact ON, 2014 ON,
`LIBIGL_SKIP_DOWNLOAD=ON`. No cache repair/deletion/download/pin change occurred.

Scoped suite command (serial, 60-minute limit):

```
ulimit -v 25165824
timeout 3600 ctest --test-dir build/task-5c-no-faster-20260909 \
  --output-on-failure -j1 \
  -R '^(spin-it-2014-|spin-2014-|source-spin-material-python$)'
```

Suite result: **36/36 passed, 0 failed, 0 skipped, 253.59 seconds**, exit 0.
This includes source-moment, Givens, sparse fill, legacy, actual source material,
Python, and both dependency tests. No suite retry was needed. Begin snapshots
are `task-5r-suite-begin.diff`, `.status`, `.sha256`; end `.diff` and `.status`
match their beginnings byte-for-byte, including active foreign tracked changes.

The actual exact-OFF child is
`build/task-5c-no-faster-20260909/spin-2014-dependency-contract-runs/run-bfd2f0c4c905aef8/enabled`.
Its cache reports exact OFF, 2014 ON and `LIBIGL_SKIP_DOWNLOAD=ON`; `build.ninja`
has no `SourceCagePaperFunctional` or `source-cage-paper` matches. The retained
contract status is `PASSED` (`task-5r-exact-off.log`). No global unrelated suite
was run, no cache was deleted, and no dependency was downloaded.

The successful Python runner retained
`build/task-5f-runner-test-c7_vvver/success/evidence.json`. Its recorded source
begin/end SHA is identically
`f27a2a8d7488dfde7ff34b71350961b03f8a9988e8f3de73007aac01189de5f6`;
full source observations compare equal, and begin/end integrity is identically
`8fbe8b23325932490b18e68cf84b4d981729834e88fc85e5d8c4520891ca6b83`.
`task-5r-python-evidence.json` preserves the compact extraction. The original
test subsequently mutates its own inputs intentionally to verify stale-result
detection; the preceding successful recorded observation is the claim here.

Stable linked identities shared with the controller before full-suite work:

```
header  4c8d0f87615a72c1635762be9aa08d91329504f192922c99f1ba2e246f97d3de
source  cb1cb31c56d0ddc6f1287f5173c12cbd9bb765c50809e2658d5d3cc25dd1c87f
core    80947d9a9ac16045850f7131b61f034c70fc6b150e86e94bf9dbf7779a33621b
```

The core SHA remained unchanged through final relink. Additional final artifact
and log hashes are recorded in `task-5r-final.sha256`.

Durable report copy of key raw-log SHA256 values (raw files remain in the task
evidence directory; only implementation and this report are committed):

```
task-5r-red.log                 4537e32152f6026a8699ff527432a1c0c7ed018f3ca79a38061f43c557b5ce74
task-5r-focused-1.log           53491ad9c6e7fb7a73a55a3e282cd06367c18741572bdea628d9e78b3aeb08ab
task-5r-focused-diagnostic.log  5dc3c9562bd8b1a1ed471f617b3dc48cfcde34e6824b0e57c7a8a83770606a06
task-5r-focused-all-fd.log      b96d7561a06aeb3e74a76a3ade7d76a7fa4066006f7eb12174ed7c9610af90aa
task-5r-focused-complex.log     8280cc26f71ab60bd5e974dc6f54dc3b1423cb31692a6f229e71baa1706d2d45
task-5r-focused-longdouble.log  f297d60267b4bba75dbda6ce93fee73587f71f77aecf8fe6f6a52b7b394d0902
task-5r-final-focused.log       df7e22c7a0c13b779ba9c090f9fec9e53b00b39a038f6e3780f5edc9f870a25b
task-5r-full-build.log          c507ef6e029189b89f60bc8fd9cb7eff3359717ca363baa2e46ef2f95da63caa
task-5r-suite.log               4aa76e7041c1c7f0edf04e312febc334d7a8376938dc44bd005bddd4c6b9da19
task-5r-exact-off.log           20c7340bd8272fb9454ef06cccc6c697bc75fa1de29d45e9a1b3c80dd16c2b6e
task-5r-python-evidence.json    2a462272c1868acfc1bd4363d136289df7c0763bb9daf12046ccdb2f2fa31303
```

## Self-review and ownership

Reviewed every production/header line, scoped CMake diff, regression branches,
resource arithmetic, native mixed-moment order, phi contraction and immutable
source access. No production changes were needed after the initial implementation.
The numerical issue was a checker precision limitation; its original failures
and approved independent accuracy resolution remain visible. `git diff --check`
is clean for owned tracked changes. Public API is unchanged from the early
controller-approved proposal.

Owned production/test/doc files are `include/spin2014/SourceCagePaperFunctional.h`,
`src/spin2014/SourceCagePaperFunctional.cpp`,
`tests/spin_it_2014_source_cage_paper_regression.cpp`,
`docs/spin2014_source_cage_paper_functional.md`, and exactly six new CMake lines.
This report and explicitly named task-5r raw logs are task evidence. No existing
source-moment, Givens, solver, legacy deformation, source partition/BBW or
validator was edited. All active foreign Faster/Plane/Polyhedral/FullRank/
Truncation changes and experiments remain unstaged. Main is read-only.

Status: **DONE** for bounded Task 5r. Implementation commit
`4b2ec54c981e59c5b37540206d0ae11ffc79f94b` (`Add actual-source joint Givens paper
functional`) contains exactly the five assigned implementation paths. This
report is committed separately using authorized volcano
`<xuxinyao307@gmail.com>`. Remaining scientific limitations are intentional
scope boundaries, not waived reproduction gates. The ordinary double FD
discrepancies remain explicit; validated derivative accuracy comes from the
approved independent high-precision and complex-step checks.
