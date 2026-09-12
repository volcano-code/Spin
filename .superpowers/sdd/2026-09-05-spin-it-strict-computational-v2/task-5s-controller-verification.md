# Task5s controller verification — scoped acceptance

BASE `ce3e159a5b736951fe83daf94766c9bbd2aadfe8`; canonical strict linked worktree
`/data/xvxinyao/moment-guided-spin-v1/code/make-it-balance-paper-sim`, branch
`spin-it-strict-computational-v2`. Git/common directories differ, no superproject.
No relevant AGENTS found. Main read-only tracked-diff SHA
`93e32af47ab1957a1efcb6716437656be129053f4269d2169378fcd95b05a35e`; foreign
Faster/Canonical/etc source and experiments excluded. /data has1.4TiB free.

## Source and design check

This continues the already approved M3 shape-regularizer step. It does not
repeat Task5r or assert full SLQP/geometry admission. The binding full spec is
the user's strict M0–M8 plan; existing local plan/ledger and task brief inspected.
Using existing-worktree, SDD, TDD and fresh verification workflow. Worker
`spin2014_source_cage_shape` is sole scoped source/build writer; controller
owns these verification helpers and independent evidence.

Baseline actual command:
`ctest --test-dir build/task-5c-no-faster-20260909 -R '^(spin-it-2014-source-cage-paper|spin-it-2014-source-cage-paper-allocation)$' --output-on-failure -j1`
passed2/2,22.77s; archivec5fbb3112335fd83eb7f356461047ca47c51eed4d46c1a33cb3b302ae6fb5b25.
No dependency download/setup/rebuild in baseline.

Controller visually inspected original2014PDF page6 and cited Sorkine2004
pages3–4. Source hashes/URL in `task-5s-brief.md`; downloaded official8MiB PDF
only into new build/task-5s-sources-uIE9MG, no original replacement or repository
dependency modification. One initial page-number filename lookup failed, then
actual zero-padded filenames were read; no source content changed.

The published small-rotation shape energy uses the fixed linearized similarity
fit, not ordinary graph smoothing or iterative polar rotations. It penalizes
actual original outer vertices only. Separate displacement penalty carries no
density and no1/2 factor. Uniform-ring choice follows the cited2004 formula;
Spin-It author discretization/weights are not thereby recovered.

## Approved integration and resource decisions

Two classes: immutable mathematical rest-mesh operator and actual-source/BBW
chain. Existing SourceCagePaperFunctional adds optional independent weights,
both0 leaves shape disabled. Geometry, moments and constraints not changed.
State keeps full original handle order and centered-scale interpretation.

Worker proposed custom3x3 block CSR with S=sum(1+degree) blocks,9S scalar values,
sorted unique directed-edge construction. Controller approved no Gram/triplets/
denseV² or VH. Construction cap conservatively includes retainedgeometry,
allindices, coefficients, initial6F edgepairs and thinSVD local scratch
8*3*maxRing*7+2048 doubles, using checked double-equivalent byte accounting.
Evaluation stores3V residual and optional3V vertex gradient; wrapper bounds
4H returnedgradient andH row separately. Upstream immutable geometry/BBW/moment
storage is separately bounded, not falsely included. Exact bound tests passed
in the final registered shape test (construction/sparse/valence/evaluation and
104-gradient/26-row boundaries, including one-below rejection).

## Independent verification path prepared

`task-5s-controller-export.cpp` extends the independently generated Task5r
actual-source twofixtures/all26handles/8beta/fiveobjectives. Uses nonzero shape
and displacement weights, fixed nontrivial frame, actual source geometry;
compares enabled/disabled moments and constraints and exports both gradients.

`task-5s-controller-check.py` independently rebuilds all eight-node interpolants,
frame and signed tetrahedron moments. New shape expectations use per-ring raw
world-coordinate QR least-squares fits, not production's normalizedSVD/LD.
It checks weighted terms, residuals, vertex/handle/fulljointgradients with
complex-step, infinitesimal nullspaces, value-only behavior and additive effect.
No native/project numerical integrator imports. Thresholds remain fixed before
actual execution; inputscale-only normalizations. Companion negativechecker
must reject36 actual-observation corruptions. These are prepared checks, not
passing evidence until compiled/linked/executed against stable implementation.

## Actual independent execution against first stable core

Worker obtained registered expectedRED (missing implementation)1/1FAIL0.03s
at build/task-5s-red.log; controller read it. Initial actual-core focused
GREEN1/1PASS0.99s at build/task-5s-green-01.log names all10 combinations and
three-step QR finite differences. Final broader scope still pending.

Controller compiled exporter with /usr/bin/c++ -std=c++17 -O2 -g -DNDEBUG
-frounding-math,2014/exact macros ON, linking the actual assigned archive,
HiGHS/Boost/GMP/MPFR/zlib. No CMake/shared build was run by controller.
Archive SHA6450f4e751ca0b3a61c5c4a7ab9c4f22cc9193715f2fd96a6cb320a101ded854.
New directory build/task-5s-controller-9h3NTi created with mktemp; all executions
single-thread OMP/OPENBLAS/MKL, memory ulimit-v25165824, existing outputs never
overwritten. The exporter binary SHA
8db6501d22e45e35a949e93f02ab5d87be41bcceb6a326cc89e6fdf66a25841b.

Two actual exports run1.txt/run2.txt are byte-identical, SHA
7d28edb3ce9125a008d11d81a112e4b6ec6b0e72bb207f5bc3b4e1a21fd46563.
Both checker runs and36-corruption negative run exited0. All113/112 columns
checked; maximum normalized value error1.194e-12, normalized joint error
2.378e-12, nonzero joint derivative relative1.314e-10. Native moments and
constraints (including Jacobians) enabled/disabled are identical in exporter.
Independently reconstructed shape residuals, vertex/handle gradients, additive
effect and infinitesimal nullspaces passed, with value-only gradients absent.

Additional explicit shape-only error/near-zero gate added before final source
freeze (strengthening, no threshold relaxation). Both detailed reruns passed:
maximum shape vertex/handle nonzero relative1.836e-11, near-zero absolute0 on
these fixtures. All36 negative cases again rejected. Retained logs check1.log,
check2.log,negative.log,check-detailed1.log,check-detailed2.log,
negative-detailed.log. Same native exports reused, not relabelled as extra
independent runs. These are synthetic integration results, not paper cases.

Final helper SHAs before source freeze:
- exporter6415bd186cc2dd654bec6ebc0b124c8fb88ef54408d470ef0c1c63bf9c29e42d
- checkerdf7d20c5c18321bebcb8d2dbca309e2c26c59aacac4ac75e11c1c0d81742b083
- negativeafc9dfeac0717ffa0bb92d391e3ba16445ae6534591fd25aa8d437fecfcd51e7

## Final scoped execution and preservation

Implementation f6e5212, implementation report d8c640a. Final scoped suite
38/38 passed, 257.11s, no skipped/NotRun. Controller inspected raw Python
success evidence: source material10 tests, 2014 bridge9 tests, miniature1 test.
Final core remains6450f4e751ca0b3a61c5c4a7ab9c4f22cc9193715f2fd96a6cb320a101ded854.
Owned nine-source manifests match byte-for-byte (47b46625d6c054dd2def23feeae30c002719565f2d5dcad7606477e6d381a1d9).
The actual run-dcda002eade83cb8/enabled child cache has2014ON/exactOFF/FasterOFF;
controller checked its Ninja graph excludes both new implementations and test.
Main tracked diff hash remains93e32af47ab1957a1efcb6716437656be129053f4269d2169378fcd95b05a35e.
Foreign concurrent additions remain outside task and graphs; no whole-tree
quiescence, clean-tree release, or full-project test claim.

Scoped RED/GREEN/build/summary/raw/child-inspection/source-manifest logs are
preserved beside this report and cmp-checked against build originals. The
first controller raw-log copy contained tool-output truncation, was detected
by cmp, then replaced using complete output and cmp-verified. Original log
never changed; final raw SHA31d99093f735325b811646419f00da5955744051fc2b0da2ddffec4f59565965.
Historical configure deprecation/unused-option warnings remain in original
build logs; the actual correct Faster cache option was OFF throughout.

## Independent review and test-only correction

Reviewer review_spin2014_source_cage_shape reviewed ce3e159..d8c640a and
approved production mathematics/integration/resources with no Critical or
Important findings. Full report retained in task-5s-review.md. Its Minor
isolated-vertex fixture finding was verified by controller at test143–146:
the previous collapsed edge caused zero-area rejection before isolation.
The other Minor finding is historical configure noise, retained rather than
erased; later runs used the existing correct cache with no configure/dependency
change. No dependency warning cleanup was authorized or needed.

Original worker strengthened the diagnostic assertion first: expectedRED1/1
failed0.01s, then reset the source fixture and obtained GREEN1/1passed1.09s.
Test-only commit06f74a5, report20970ab. Controller read both logs and fixdiff;
sha256sum-c verified all eight production/doc/CMake inputs pluscore unchanged,
and the new test source/binary match greenmanifest. Testbinary now
adbd83579130fbe0ba7eae15e481b7b023501e578454bc407fe84a6e31c57629.
Fixlogs and before/RED/GREEN manifests preserved with cmp against originals.
The prior38/38 batch did not prove isolated-vertex coverage; the repaired
registered test supplies that missing coverage. No post-fix38batch is claimed.
Independent scientific exporter/check evidence remains bound to unchangedcore.
Scoped re-review of d8c640a..20970ab accepted both findings as addressed,
no new Critical/Important breakage; verbatim task-5s-re-review.md retained.
Controller accepts this source-bound shape-functional layer only. No full
joint SLQP, deformation geometry admission, original-case reproduction or
global publication status is granted. Overall remains PARTIAL_BLOCKED.
Next scoped algorithm work is actual-material beta/phi/full-handle joint SLQP
using these complete first-order terms, followed by geometry commit gates.
