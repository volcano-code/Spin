# Task5r independent verification

BASE `86234ff`, canonical strict worktree at
`/data/xvxinyao/moment-guided-spin-v1/code/make-it-balance-paper-sim`.
Branch `spin-it-strict-computational-v2`; linked Git directory differs from
common Git directory and no superproject, so existing isolation retained.
No AGENTS found in scoped source or parent paths. Main is read-only; its
tracked-diff SHA at start is
`93e32af47ab1957a1efcb6716437656be129053f4269d2169378fcd95b05a35e`.
Foreign Faster/Plane/Polyhedral/FullRank/Truncation code and experiments are
excluded from this task. Their concurrent changes are not overwritten.

## Baseline and scope

Actual baseline command in existing assigned build:
`ctest --test-dir build/task-5c-no-faster-20260909 -R '^(spin-it-2014-source-cage-moments|spin-it-2014-givens)$' --output-on-failure -j1`.
2/2 passed, 1.82s; starting archive SHA
`1bdd1c8eb43e474359d6cbf6732f0c5faf18768fb7d4939ccaf4eff3306e86f2`.
2014/exact ON, Faster OFF, LIBIGL_SKIP_DOWNLOAD ON. No setup downloads needed.
Disk available1.4TiB. Owned builds max3 jobs/24GiB address space, tests serial.

The approved next step is bounded integration of source moments, fixed frame,
five explicit Givens objectives/constraints and source sparse fill term.
It does not add the still-required paper shape regularizer, joint SLQP,
deformed geometry admission, original cases or release authority.

Original local paper extract `build/strict-paper-sources/spin2014.txt` SHA
`5cf5ab047189575df48b01796de73c72ceadc602a3380b74f33d2dda81459161`:
Section4 explains user spin/support point and axis frame, Eq5–8 Givens
constraints, and the constrained axial shift `sz^2/s1` for tops. Section6.1
requires separate rotation-aware shape and displacement regularizers. These
cannot be silently replaced by the legacy uniform shape graph.

## Independent checker design

Helpers are `task-5r-controller-{export.cpp,check.py,negative.py}` in this
directory. Export actual two sourceD5/cage2/grid3/fullBBW fixtures, five modes
per fixture, all26handles/8beta, mixed handle scales/translations and beta.
Nonzero fixed origin and nontrivial proper rotation are independent of CoM.

NumPy/stdlib checker reconstructs eight-node rest interpolation itself, moves
all points with all handles, applies the frame to POINTS, then integrates signed
tetrahedra. It constructs and rotates the3x3 inertia tensor to evaluate all five
objectives and constraints. It does not import production integral/Givens code,
nor multiply the production10x10 frame map to obtain expected values.
Fill penalty uses explicit actual edge differences. Complex-step checks all113
nonstatic/112static columns of raw moments, objective and equality constraints.
Absolute scales come from the input rest bounding diagonal/density/weights,
not final mass or measured error. Derivative gate combines absolute1e-10 of
fixed scales and relative1e-7; values target1e-10. Native three-step FD remains
separate. Twenty corruption checks target values, gradients, equality/raw
Jacobians, value-only observations and fixed frame.

## Actual independent execution

Controller compiled the standalone exporter against stable actual archive SHA
`80947d9a9ac16045850f7131b61f034c70fc6b150e86e94bf9dbf7779a33621b`, using
`/usr/bin/c++ -std=c++17 -O2 -g -DNDEBUG -frounding-math`,2014/exact ON and
the assignedcore/HiGHS/Boost/GMP/MPFR/zlib. Helpers use canonical paths.
Directory `build/task-5r-controller-Bar1Rp` was created exclusively with mktemp.

Two actual native executions produced byte-identical complete exports:

- Exporter SHA `479f01c438ffe8c1240a0bc5f4427d67dc188880a9035b4edd9ecad530b05a29`.
- run1.txt and run2.txt SHA
  `d7863ca8bc8ad07e4ffdbe1d21627076cef913649b2e5c8dd79c109456188568`.

Both independent checker executions passed all10 fixture/mode combinations.
Largest normalized value error1.723e-12, joint derivative error3.426e-12.
Additional per-entry diagnostics and stricter explicit checks found maximum
relative error1.201e-10 for |reference derivative|>1e-7 (required<=1e-6);
near-zero absolute error0 on these fixtures (required<=1e-11). These extra
checks strengthen rather than replace the fixed-scale mixed gate. All20
mutated objective/constraint/Jacobian/value-only/frame observations were rejected.

Retained logs: initial `check1.log`, `check2.log`, `negative.log`;
additional diagnostics `check-detail1.log`, `check-detail2.log`; explicit
stronger-gate reruns `check-strict1.log`, `check-strict2.log`,
`negative-strict.log`. Existing outputs were not overwritten. All completed
commands exited0 and cmp confirmed identical scientific export bytes.

## Finite-difference investigation (not hidden)

The worker's initial ordinary double-precision centered FD test found actual
threshold failures at the three fixed step sizes. Changing fixture, increasing
the1e-6 threshold or enlarging its near-zero category is not authorized.
Controller's independent complex-step evidence above shows the actual analytic
derivatives agree with a separately integrated mathematical function, but it
does not make the failed ordinary FD runs pass.

Initial ruling: improve the independent FD reference's arithmetic
to long double throughout and use the five-point centered stencil at the SAME
three steps, keeping derivative thresholds unchanged. This targets subtractive
roundoff and O(h²) truncation, not the production function or input. Preserve
all original FD failures and do not silently claim they are green. Acceptance
requires actual new evidence and scoped review.

Controller inspected all30 `LD5` summaries in retained
`task-5r-focused-longdouble.log`: every fixture/mode/step's independent
five-point derivative comparison meets the unchanged limits. Largest relative
error6.85686e-9, largest absolute error7.79168e-13. This particular diagnostic
run still exited1 because the old double-FD-failure counter remained active;
it is NOT recorded as a passing overall test. The final test subsequently
distinguished those old diagnostics from the mandatory high-precision/complex-step
accuracy gates. Native complex-step independently reached relative6.02589e-11.

Helper SHAs after explicit stronger checks:

- Exporter source `c59c90c53005f123d52615fee07c0f934b7b30857c5a2adcca7d2ced14577e37`.
- Checker `c5e566faf24555804e05b87771916204db4a1309a5266d64e99ce5d757da7386`.
- Negative checker `6af35430c95e1c27acb46419a3236dedd9c4b966b69f10dfffdf5e3a58b70213`.

Final focused native run now exits0 while explicitly reporting41 ordinary
double finest-step discrepancies as retained diagnostics; all30 long-double
five-point and10 complex-step comparisons remain accuracy assertions.
Controller read that output. Its SHA is
`df7e22c7a0c13b779ba9c090f9fec9e53b00b39a038f6e3780f5edc9f870a25b`.
Old all-FD failure SHA
`b96d7561a06aeb3e74a76a3ade7d76a7fa4066006f7eb12174ed7c9610af90aa`;
intermediate long-double run still failing the original counter SHA
`f297d60267b4bba75dbda6ce93fee73587f71f77aecf8fe6f6a52b7b394d0902`.
Read-only compiler macro check confirms `__LDBL_MANT_DIG__=64`,
`__DBL_MANT_DIG__=53`. No fabricated higher precision is claimed.

After full assigned relink, production header/source and core hashes still
match the independently tested implementation:

- Header `4c8d0f87615a72c1635762be9aa08d91329504f192922c99f1ba2e246f97d3de`.
- Source `cb1cb31c56d0ddc6f1287f5173c12cbd9bb765c50809e2658d5d3cc25dd1c87f`.
- Core `80947d9a9ac16045850f7131b61f034c70fc6b150e86e94bf9dbf7779a33621b`.

## Final scoped suite

Actual final command completed **36/36 PASS**,253.59s, no Not Run, skip or
freshness recovery needed. Controller read the complete suite result;
`task-5r-suite.log` SHA
`4aa76e7041c1c7f0edf04e312febc334d7a8376938dc44bd005bddd4c6b9da19`.
This is the assigned2014/source-related subset, not the full project or all
foreign concurrent code. Production archive hash remains the one independently
linked/tested above.

Controller separately parsed actual Python runner evidence
`build/task-5f-runner-test-c7_vvver/success/evidence.json`: equal complete source
observations, begin=end SHA
`f27a2a8d7488dfde7ff34b71350961b03f8a9988e8f3de73007aac01189de5f6`.
It inspected actual retained dependency child
`build/task-5c-no-faster-20260909/spin-2014-dependency-contract-runs/run-bfd2f0c4c905aef8/enabled/`:
2014ON/exactOFF/FasterOFF/downloads disabled; generated Ninja graph contains
neither SourceCagePaperFunctional nor source-cage-paper. No dependency fallback
or download/pin change was used.

## Independent review and repair tracking

Initial implementation `4b2ec54`, report `a603876`. Fresh read-only review
confirmed the numerical integration but found an unbudgeted Eigen GEMM packing
workspace in the two dynamic Jacobian transforms. Controller verified the
resource contract and assigned a constant-workspace repair to the original
worker. Review details are in `task-5r-review.md`; initial verdict NEEDS_FIXES.
All results above remain evidence of the pre-repair implementation, not an
automatic verification of the repair. A changed archive requires independent
exporter relinking and new retained runs, followed by scoped re-review.

Controller rechecked the complete `task-5r-final.sha256` manifest: all17 entries
matched. Main tracked-diff SHA still equals the start hash; this says nothing
about completeness of unrelated untracked-file preservation. An accidental
read-only evidence lookup in main returned missing-file errors; no write ran
there. The subsequent canonical strict-worktree checks found the expected files.

Independent strict logs are additionally preserved byte-for-byte as
`task-5r-controller-check-strict1.log`, `task-5r-controller-check-strict2.log`
and `task-5r-controller-negative-strict.log` in this directory. Raw native
exports remain at the hashed build paths, not staged as source.

## Repair1 independent rebind

The production-used transform now copies each column into a fixed ten-vector,
then applies the fixed10x10 map. Its allocation probe rejected the former GEMM
route for both10x8 and10x104; the repaired helper passes the same guard and
independent scalar reference. Macro settings are local to the probe target;
normal functional tests still link the actual core. Controller read RED/GREEN
logs and the helper implementation. No geometry/objective/tolerance change.

Controller relinked the unchanged exporter source against the new archive
`c5fbb3112335fd83eb7f356461047ca47c51eed4d46c1a33cb3b302ae6fb5b25`.
New exclusively created directory `build/task-5r-controller-fix1-EuMfa7`;
same compile command and fixed single-thread environment as above.

- Export binary `97aa8d2f4497763af5ae042774098330822fe395e47daeab8bc6702e357c2e60`.
- Both new runs `0ecedced12e6cabcae8a91cac5d8525d5d083d8d21c03fe26e701b5780ac7599`.
- Both independent10-case checks passed, normalized value maximum1.723e-12,
  normalized joint derivative maximum3.432e-12, nonzero relative1.196e-10.
- All20 corruption checks rejected. Commands exited0; cmp confirmed equality.

Preserved copies `task-5r-controller-fix1-{check1,check2,negative}.log` bind the
new checks; old outputs remain intact. Scientific hashes differ from the old
implementation because floating-point contraction order changed; only the two
same-version runs are required and observed to be byte-identical. The values
and derivatives separately meet the unchanged independent numerical gates.

Repair commit `c79b9ed`, report-only commit `77e4f38`. Controller read the fresh scoped suite: **37/37 PASS**,
249.30s, no skipped/Not Run/failed tests or freshness recovery. Final focused
numerical log confirms64 mantissa bits, all30 LD5 comparisons (maximum relative
6.85705e-9, absolute7.79168e-13),10 complex-step comparisons (relative6.02589e-11),
and retained41 ordinary double discrepancies. The extra allocation test is
included in the37 count, not merely a manually run helper.

Controller jq verification of actual
`build/task-5f-runner-test-8yxf_px7/success/evidence.json` confirms identical
complete source observations and contentSHA
`7630a0025eb24e3904e9399d71473f0ba5cd10ebe115fb13a59502b39a706228`.
Begin/end integrity both
`ba236d3a45f1d8b103d40bc3e917a765d8f22bac5661534a3788cf8ef4ab770e`.
Actual dependency child `run-2a93d4cd902916b0/enabled` cache is2014ON,
exactOFF,FasterOFF,downloadsOFF; Ninja graph has no functional or paper probes.
Foreign untracked `tests/canonical_solver_process_replay.cpp` appeared during
the suite, was not built or modified here, and is excluded from Task5r. The
runner's actual selected-source observation remained equal; this is not a
whole-worktree immutability claim.

Controller validated all19 entries of `task-5r-fix1-final.sha256`, including
the unchanged historical manifest itself. No source/core mismatch is hidden by
the separate pre-repair and post-repair manifests. Builds are complete and
ownership released; scoped re-review package is `review-a603876..77e4f38.diff`.

## Bounded acceptance

Original independent reviewer re-reviewed the repair: all findings addressed,
no new Critical/Important breakage and no out-of-scope observation. See
`task-5r-review.md`. Task5r is accepted as actual-source joint paper-function
integration with independently checked first derivatives and bounded transform
scratch. It is not full joint SLQP, paper shape regularization, geometric
admission, original paper experiments or a release certification.

At final controller check, protected main tracked-diff SHA still matches
`93e32af47ab1957a1efcb6716437656be129053f4269d2169378fcd95b05a35e`.
Additional foreign Canonical/FullRank/Polyhedral files and experiments appeared
externally after the run; they remain unmodified and excluded from the scoped
commit/review/test claim. No whole-worktree-clean statement is made. Overall
strict reproduction status remains **PARTIAL_BLOCKED**.
