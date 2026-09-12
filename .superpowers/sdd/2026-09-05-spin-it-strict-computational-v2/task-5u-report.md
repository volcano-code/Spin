# Task5u — actual-source represented geometry admission foundation

Status: DONE_WITH_CONCERNS. Standalone checker implemented and all40registered
related tests pass; two actual-source negative-fixture gaps remain explicit.

## Scope and implementation

Implemented the approved standalone immutable `SourceCageGeometryChecker` around
the actual `SourceCagePaperFunctional::evaluateValues(fullState)`. No externally
supplied evaluation, mesh or successful Boolean is accepted. The typed audit has
a private constructor and const facts; ordinary copying preserves its evidence.
There is no solver/precommit integration in this task and no original-case,
physical, sealed-artifact, feasible-witness or convergence grant.

The original outer rows/faces and every conforming oriented primitive reference
are used. Signed nested outer material (including original voids), nondegenerate
closed manifold boundaries, exact self-intersection, one positive primitive
component, actual-material subset, exact outer/domain separation and all ordered
unordered-pair checks are required before CHECKS_PASSED. Strictly disjoint AABBs
can prove an impossible overlap; face/edge/point-sharing candidates still use the
existing budget-aware exact Nef relation. Beta never drops or binarizes cells.

Scientific hashes bind complete raw state, actual immutable source/frame/config,
all settings, actual evaluation output, every attempted stage/outcome, complete
exact fact digests, IDs/beta owners, measured rational clearance and resource
observations. No wall-clock fields are included. Resource errors stop further
geometry; missing/incomplete Boolean evidence stays non-admitted. Actual source
failures can occur before geometric checks and retain UNKNOWN rather than
inventing a later failure certificate.

Owned implementation paths only:

- `include/spin2014/SourceCageGeometry.h`
- `src/spin2014/SourceCageGeometry.cpp`
- `tests/spin_it_2014_source_cage_geometry_regression.cpp`
- `docs/spin2014_source_cage_geometry.md`
- six added registration lines in `CMakeLists.txt`, gated by 2014 and exact ON

No generic solver, upstream geometry, Faster/Canonical/Quadric, protected main,
or foreign untracked implementation was changed. This report is committed
separately from those five paths. Probe logs and hashes are preserved locally.

## Approved resources and limitations

The controller approved the API/resource proposal before production. Finite
defaults: 100,000 points; 400,000 total triangle references; 4,096 primitives;
1,000,000 pairs; 1,100,000 records; 512 MiB explicit checker accounting;
2,000,000 live explicit Nef complexity; 64 connected shells; 8,192 retained
exact-string bytes per string. All integer products/additions and existing
signed-index narrowing are checked before size-dependent allocations.

Preflight counts actual full references and worst-case `3+3P+P(P-1)/2` records.
It accounts source output/split-state/BBW row/beta-Laplacian/optional shape
residual scratch, local and retained region meshes, compact-index/DSU scratch,
all potentially retained audit diagnostics including F(F-1)/2 intersecting face
pairs, hierarchy storage and clearance strings. Diagnostic payloads have a
factor-four allowance for vector capacity growth and simultaneous facts copies.
No dense primitive-pair matrix is allocated. All primitive regions are retained
only under the explicit mesh/accounting preflight and one live Nef budget.

The documented component-volume string bound is explicit: finite binary64
coordinates share denominator 2^1074 with integer numerator magnitude <2^2098.
With F<2^30 faces, six cubic determinant terms per face and division by six,
unreduced signed volume uses <6,328 numerator bits and <3,226 denominator bits;
decimal sign/slash plus both integers needs <2,900 characters. Accounting keeps
at least 8,192 bytes per component even under a smaller retained-string cap.

This is not an exact allocator/RSS proof. Native CGAL/GMP arithmetic, BVH and
allocator internals, opaque Boolean peaks and temporary clearance arithmetic/
string storage remain subject to the required external 24 GiB virtual-memory /
60-minute process guard. Persistent allocation failure may prevent producing
an audit/digest and propagate; it cannot produce admission. All heavy runs were
sequential, OMP/OPENBLAS/MKL one thread, ulimit-v 25165824, timeout 3600, at most
three compilation jobs, using the assigned existing build only.

## TDD evidence and focused verification

The implementer prompt, TDD skill and writing-good-tests reference were read.
The test and gated CTest registration preceded production. Initial registered
RED in `task-5u-red.log`:

```
ctest --test-dir build/task-5c-no-faster-20260909 \
  -R '^spin-it-2014-source-cage-geometry$' --output-on-failure
FAIL: actual-source represented geometry admission is not implemented
0% tests passed, 1 tests failed
```

The test used the existing __has_include absent-component behavioral failure
pattern; it compiled and ran rather than failing to compile. Actual-source
fixture tests were present before the production implementation was added.

An additional registered RED caught an existing-backend limit mismatch in the
new settings validation: `INT_MAX+1` Nef complexity was not rejected by the
checker constructor. `task-5u-nef-cap-red.log` reports
`FAIL: unrepresentable exact backend cap rejected early`. Production validation
was corrected from LONG_MAX to the actual backend's INT_MAX limit.

An intermediate analytic scale fixture failure was a test parameter error,
not a production defect: handle scale is local about each handle center.
True common world scale s requires per-handle translation `(s-1)*center+T`.
The corrected test uses s=2 and T=(.125,-.25,.375), with volume scaled by eight.

The fixed source cube is [5/16,11/16]^3, design [.34,.66]^3, source depth3,
adaptive levels2..3, cage2, FEM grid3; actual full BBW has26handles and8primitives.
Its rest-domain clearance is analytically1/16, but actual BBW summation rounds
identity output. The attempted exact .0625 threshold correctly rejected and is
preserved in `task-5u-clearance-red.log`; no geometry was snapped or repaired.
The controller independently recomputed the actual triangle distances using
Fraction arithmetic and supplied the exact minimum:

```
19807040628565802923409276929/5070602400912917605986812821504
sqrt = 0x1.fffffffffffc0p-5 (exactly representable)
next = 0x1.fffffffffffc1p-5
```

The fixed literal equality threshold passes and its adjacent successor rejects.
Focused final GREEN: `task-5u-green.log`, 1/1 pass, 18.12 seconds. Coverage:

- Actual identity/repeated scientific digest; independently analytic material
  volume and common positive scale/translation volume.
- Changed nonuniform beta/phi retaining the geometry verdict but changing input
  identity; frame and optional shape/displacement regularizers; all five goals.
- Real valid nonaffine motion; all complete pair IDs/counts (8 cells, 28 pairs),
  shared-boundary candidates; exact rational distances and witness indices.
- Actual outer self-intersection (handle2, tx+1), primitive self-intersection
  (handle3, tx+1), negative primitive orientation (handle4, tx+1), protrusion
  (handle1, tx-1), using complete actual controls and beta0.
- Exact and one-below positive resource caps, every zero cap, unrepresentable
  settings, malformed state dimensions, explicit source failure UNKNOWN, tiny
  exact-string cap, live Nef exhaustion and terminal resource stage.
- Actual depth4 source with original inner void [.48,.52]^3 and off-center
  design, retaining16original outer rows and26handles, passes signed nesting.
- At common scale2^-49 plus translation(.5,.5,.5), real source integration
  rejects degenerate triangles before geometry: UNKNOWN at SourceEvaluation.
- Typed audit cannot be constructed from public facts; copies preserve evidence.

## Bounded probes and unresolved fixture coverage

Explicit diagnostic CLI arguments cannot be triggered by inherited environment
variables; no-argument registered CTest always runs assertions. The initial
environment-triggered probe form was corrected after controller review.

Preserved probe logs:

- `task-5u-motion-probe.log`:78 actual full-control single-handle motions,
  handles0..25, x translations+.5,+1,-1. Finds the invalid outer/primitive and
  protrusion cases fixed in the registered regression.
- `task-5u-collapse-probe.log`:11 centered common-scale states2^-45..2^-55,
  fixed translation(.5,.5,.5), static objective.45,46pass;47,48subset-reject;
 49..55fail the actual source triangle integrator and remain UNKNOWN.
- `task-5u-overlap-probe.log`:32 actual disconnected two-cell depth4 source
  motions. Directions preserve outer rows numerically via
  `e_h-Wouter^+ Wouter e_h`; handles0,4,10,12,13,15,21,25 and amplitudes±1,±2.
  Results are pass, primitive invalid, or subset rejection.

None of these121bounded cases reached an exact outer-contact rejection or a
positive-volume primitive-pair-overlap rejection after all prior checks passed.
Those two negative branches are composed from existing exact APIs but have not
been validated end-to-end by an actual-source fixture here. The controller
explicitly approved stopping the bounded search while retaining this limitation
for review and later solver integration. No omitted failure case is described
as tested, and no synthetic fixture is an original-paper reproduction.

## Frozen final group, hashes and commits

Final related group: **40/40 PASS,273.39seconds**, with no skipped, NotRun,
disabled or failed case counted as a pass. The complete group ran once after
source freeze, against a scoped relink of all existing2014CXX executable targets
and `_make_it_stand` at-j3. No foreign target was built.

```
ctest --test-dir build/task-5c-no-faster-20260909 --output-on-failure -j1 \
  -R '^(spin-it-2014-|spin-2014-|source-spin-material-python$)'
```

Source manifests `task-5u-source-begin.sha256` and `task-5u-source-end.sha256`
match byte-for-byte. Their manifest SHA256 is
`63531efe6c53b7a603e1c328c9bbdc509bd32a9755efb01026856e35e06d4e7f`.
Relink log is `task-5u-final-relink.log`; test group log is
`task-5u-final-group.log`, all beside this report. Group-log SHA256 is
`acd9a1e95e161c7533aa8df58e2dc15eb0d679f5bee9963078e0128857339480`.
Preserved complete raw CTest output is `task-5u-final-raw.log`, SHA256
`16bf2cde59473f2aeb7fc2aec323e38916ef48c60182d74728eb4c4be0d1c1a4`.

Current-source Python evidence in that raw log: lines475,477 report10tests OK in
23.338s for source-spin-material-python; lines991,993 report9tests OK in0.034s
for spin-it-2014-python; line1014 reports the existing miniature test OK.
The Python module SHA256 is
`5b1092d3baf5dafca9cfe8bd27df1469e583f1b28bef750aaa6302b4c0b34453`.
This is fresh test evidence, not a prior task's log or a paper-success claim.

The actual fresh exact-OFF dependency child is:

```
build/task-5c-no-faster-20260909/spin-2014-dependency-contract-runs/run-2a58452101d1cb98/enabled
```

The final raw log line1030 identifies its retained parent. Its cache records
2014ON, exactOFF, FasterOFF and downloads skipped. `task-5u-exact-off-graph.log`
records those cache lines and an actual graph search confirming no
`SourceCageGeometry`, `source_cage_geometry`, or `source-cage-geometry` source,
header or target name. Graph SHA256 is
`a708b48cc2a82450d7f094f32d81e5044bccc4dae623463e2d3dfac26cc6e3da`.
The registered dependency-contract and dependency-isolation tests both passed
(89.66s and7.30s). The exact-OFF graph was inspected, not guessed from gating.

Full artifact manifest: `task-5u-final-artifacts.sha256`. Frozen core SHA256 is
`c9b4a09c9b128e558feb2c13896ca3e074ec951a64c4eab79cb52484759fe72e`;
new native binary SHA256 is
`cc9388aaf9780a180414f861a99710ab91f47b55c3763ebeeb730eac868fe14f`.

Implementation commit: `5b518bc feat(spin2014): audit actual-source represented geometry`.
This report is the following separate scoped report commit. No merge/tag/push.

## Self-review

Reviewed complete implementation, header, test and CMake diff. Fixed the Nef
signed-index limit, strengthened simultaneous diagnostic and optional source
scratch accounting, documented a finite dyadic volume-string bound, and changed
probe activation to explicit argv. The controller's independent export/oracle
work is separate and its final evidence is not claimed as worker-owned tests.
Remaining concerns are the explicit native allocator limitation and the two
unreached actual-source negative cases above. No solver hookup was attempted.
