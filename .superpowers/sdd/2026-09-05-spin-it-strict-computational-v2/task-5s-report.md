# Task 5s implementation report

Base: `ce3e159a5b736951fe83daf94766c9bbd2aadfe8`, branch
`spin-it-strict-computational-v2`, worktree
`/data/xvxinyao/moment-guided-spin-v1/code/make-it-balance-paper-sim`.
Status: DONE; all 38 covering tests passed. Scoped commit information is below.

## Scope and formulas

Implemented the immutable mathematical `LinearizedRotationLaplacian`, immutable
actual-binding `SourceCageShapeFunctional`, and optional composition in the
existing `SourceCagePaperFunctional`. No legacy solver/source moments or foreign
Faster/Canonical/FullRank/Plane/Polyhedral/Truncation files were edited or staged.
Only component/test registrations were changed in the 2014 AND exact CMake gate.

The mathematical operator follows the brief's Sorkine Section 3 Eq.2,5,6,8–12
coefficient order `(s,hx,hy,hz,tx,ty,tz)`, local order `{i}+sorted_unique(Ni)`,
and fixed sign `Ri=Di Ai^+-Li`. Equal-weight least squares uses every local
vertex, uniform one-ring `Li`, rest-only local centering and isotropic max-radius
normalization, and full-rank-7 thin Jacobi SVD. This normalization preserves the
dimensionless original-unit residual action. Relative rank tolerance defaults
to `1e-12`, is configurable/hashed, and rejects outside `(0,1)`. Planar,
noncollinear rings pass. No normal-equation inverse, damping, pin, omitted ring,
or singular-direction truncation is used.

The weighted energies are `gammaD,L*||LD vPrime||²` and
`gammaD,v*||vPrime-vRest||²` with no half or density. The vertex gradient is
`2 gammaD,L LD^T(LD vPrime)+2 gammaD,v(vPrime-vRest)`. All actual original
outer vertices and all actual cage handles participate. The wrapper uses the
unchanged centered scale chain `wih*(viRest-centerh)` and translation `wih I`,
computed from the binding's immutable actual BBW field/embeddings without
normalization/clamping. The paper passes original outer rows from its actual
moment result without re-deformation. Beta/phi additions are exactly zero;
moments, constraints, and their Jacobians remain unchanged.

Both weights default to zero; this disables construction/evaluation of LD.
Displacement-only remains an honest enabled source-operator mode. Value-only
has identical energies and no derivative arrays or derivative budget need.
Source/BBW/frame/state/config/operator/raw-output identities bind this mode;
no acceptance flags or scientific completion grants are added.

## Source authority and exclusions

Read the controller's complete brief first and followed the provided
implementer prompt and TDD skill. The controller supplied the visually checked
source equations. All four source hashes were independently checked locally:

| Source | SHA256 |
| --- | --- |
| `/home/lin/下载/SpinItPaper.pdf` | `3add70de03df9b0e173028d55c775e8520fa89a6b21d9b269f926ccbe8a76e4e` |
| `build/strict-paper-sources/spin2014.txt` | `5cf5ab047189575df48b01796de73c72ceadc602a3380b74f33d2dda81459161` |
| `build/task-5s-sources-uIE9MG/sorkine2004.pdf` | `246498341881005d40a27309c5aca0e3d6efd2de5d29a85f95bbcdcc978dfab4` |
| `build/task-5s-sources-uIE9MG/sorkine2004.txt` | `703ddc9364c235d6eb440920c04b04cd28f079588d258514a6bc460eb34b7095` |

Uniform weights are the published Sorkine variant but their selection for
Spin-It is formula reconstruction, not recovered author experiment parameters.
The mathematical vertex API does not certify input deformations. Broader
manifoldness/self-intersection checks, geometry legality, full joint Hessian,
beta/phi-dependent shape operator, large-rotation iteration, solver, original
experiment reproduction, and paper completion remain excluded. Tests are
`SYNTHETIC_REGRESSION`.

## Approved resources

Before coding, sent the public API and accounting to the controller. The
controller explicitly accepted the custom block CSR approach and the following
accounting without later changes. Let `V,F,d,S=sum(1+degree)` be original
vertices, faces, max unique degree, and block count; `m=3*(d+1)`.

Defaults: vertices 1,000,000; one-ring neighbors 4096; stored scalar
contributions 30,000,000; construction/evaluation workspace each 10,000,000
doubles; gradient/weight-row capacity each 1,000,000 doubles. All must be
positive. No cap silently overrides another.

Block CSR stores `V+1` size_t offsets, `S` size_t columns and all `9S`
row-major coefficient doubles, including structural zeros. Initial directed
edge capacity is `6F` pairs even after unique erasure. Simultaneous constructor
accounting is:

```text
ceil((3V*sizeof(double) + 3F*sizeof(int)
      + (V+1+S)*sizeof(size_t) + 9S*sizeof(double)
      + 6F*sizeof(pair<size_t,size_t>))/sizeof(double))
+ 8*m*7 + 2048
```

The initial geometry/edge bound is checked before edges; the full bound before
copied geometry, CSR and SVD allocations. Eight `m×7` buffers conservatively
cover A, thin U, scaled work, QR and packing. The fixed allowance covers small
SVD matrices/metadata and bounded stack work. Inspected the installed Eigen
JacobiSVD allocator: `m>7` only allocates the tall QR preconditioner, U is `m×7`,
V is `7×7`, and thin-U Householder workspace has 7 entries. Local coefficient
contractions are explicit scalar loops. No triplets/sparse conversion, global
Gram or global GEMM is present.

Operator evaluation arrays require `3V` for values or `6V` for gradient plus
returned residual; the wrapper moves them and adds separately capped `4H`
output and `H` row. No duplicate V-sized arrays. Existing upstream
source/binding/BBW/conforming/moment memory is separately bounded and excluded.
Constant evaluation control/hash objects are outside the variable array budget.
Integer/index/byte/cap failures and allocation failures are `RESOURCE_LIMIT`
length_error; invalid input is invalid_argument; nonfinite arithmetic is
overflow_error. None becomes infeasibility or silent zero energy.

## TDD and verification evidence

RED was a real registered test run before production code existed, using a
missing-header branch whose executable fails with the explicit missing-feature
assertion. Full regression expectations were authored before implementation.

```text
cmake -S /data/xvxinyao/moment-guided-spin-v1/code/make-it-balance-paper-sim
      -B /data/xvxinyao/moment-guided-spin-v1/code/make-it-balance-paper-sim/build/task-5c-no-faster-20260909
      -DMAKE_IT_STAND_BUILD_SPIN_IT_2014=ON
      -DMAKE_IT_STAND_BUILD_EXACT_GEOMETRY_CHECKS=ON
      -DMAKE_IT_STAND_BUILD_FASTER_SPIN=OFF -DLIBIGL_SKIP_DOWNLOAD=ON
cmake --build <same absolute B> --target spin-it-2014-source-cage-shape-regression -j3
ctest --test-dir build/task-5c-no-faster-20260909 --output-on-failure -j1
      -R '^spin-it-2014-source-cage-shape$'
```

The first configure spelling `MAKE_IT_STAND_BUILD_FASTER_SPIN` is not the
project option; the existing actual `MAKE_IT_STAND_BUILD_SPIN_IT_FASTER:BOOL=OFF`
cache entry remained OFF and was inspected. Thus Faster stayed excluded.
No dependency downloads or pin changes occurred. Existing libigl emits its
pre-existing CMake minimum-version deprecation warning on configure.

`build/task-5s-red.log`: **1/1 failed**, 0.03s, exact expected output
`FAIL: source linearized-rotation shape operator is not implemented`.
SHA256 `20d4cb9d4f0e7148463ae789f674071d13d642da50707ceec75446725ae60d38`.
Related logs: `task-5s-configure-red.log`, `task-5s-build-red.log`.

First full production/regression build logs:
`build/task-5s-build-green-01.log`, `build/task-5s-build-green-02.log`.
The first build had cached the initial missing-header conditional; the
subsequent test cap additions triggered recompilation of the complete test.
No passing claim was made for the initial executable.

GREEN command:

```text
cmake --build <same absolute B> --target spin-it-2014-source-cage-shape-regression -j3
ctest --test-dir build/task-5c-no-faster-20260909 --output-on-failure -V -j1
      -R '^spin-it-2014-source-cage-shape$'
```

`build/task-5s-green-01.log`: **1/1 passed**, 0.99s; all mathematical checks
and both actual source/BBW fixtures in all five modes, 112 static/113 angle
columns, passed independent raw-coordinate long-double QR finite differences
at `3e-4,1e-4,3e-5`. SHA256
`16e88f4ba0f629a2066d4c168dfe4f7fb15181b8871abcbb4d9ee1530756a810`.

The regression additionally covers rest/translation/uniform-scale/infinitesimal
rotation nullspaces, positive shear/local/finite-rotation energy, simultaneous
rest/current covariance and scale-squared units, every vertex derivative at
three steps, planar validity, rank failures, NaN/Inf, malformed faces/shapes,
zero area/isolated vertices, overflow, all zero budgets, exact-boundary/one-below
construction/sparse/valence/evaluation/104-gradient/26-row caps, linear sparse
fan growth, exact repeated outputs/copies, and changed actual BBW/source hashes
under zero motion. Near-zero FD rule is fixed `1e-10*inputScale`, smooth-region
relative `1e-6`; inputScale is declared from rest input before observations.
Tests also verify density independence, additive energies, unchanged moments/
constraints/Jacobians, zero beta/phi shape derivatives, independent weights,
disabled no-LD behavior under insufficient vertex cap, and value-only operation
under derivative shortages. No threshold was relaxed and no failure erased.

Final native relink/build log: `build/task-5s-build-final.log`; all 30 distinct
native test targets selected by CTest's registered command paths, plus
`_make_it_stand`, were built with `-j3`. Shell resource ceiling for builds and
the suite was `ulimit -v 25165824`. Final suite command:

```text
timeout 3600 ctest --test-dir build/task-5c-no-faster-20260909
  --output-on-failure -j1
  -R '^(spin-it-2014-|spin-2014-|source-spin-material-python$)'
```

Final suite: **38/38 passed, 0 failed, 0 skipped/NotRun, 257.11s**. New shape
test 1.07s; existing source paper test 22.77s; dependency contract 99.31s and
dependency isolation 6.37s. `build/task-5s-covering-final.log` SHA256
`1afdc02731861357a146fe8bb50403b9db994442c002998b9034c0c4fb7aa1aa`;
final build log SHA256
`9fb0ddf3fde3d7c27cdfb686626466d1ff1c80eb303ea06e43cb44cfed67aa66`.

The full current CTest output was copied to
`build/task-5s-covering-final-raw.log`, SHA256
`31d99093f735325b811646419f00da5955744051fc2b0da2ddffec4f59565965`.
Current Python success evidence there: lines 418–445 show
`source-spin-material-python`, **10 tests OK**, 23.288s unittest/23.38s CTest;
lines 935–962 show **9 tests OK** for `spin-it-2014-python`; lines 964–983 show
the paper miniature **1 test OK**. These logs are from this final run, not
reused prior-task evidence. The copy operation issued only the shell's
`cp -n` portability warning; the destination was verified by its content/hash.

Scoped begin/end manifests are `build/task-5s-source-begin.sha256` and
`build/task-5s-source-end.sha256`; `cmp` passed exactly and both have content
SHA256 `47b46625d6c054dd2def23feeae30c002719565f2d5dcad7606477e6d381a1d9`.
The stable core remained
`6450f4e751ca0b3a61c5c4a7ab9c4f22cc9193715f2fd96a6cb320a101ded854`.
Thus every owned source/test/doc/CMake input was frozen across verification.

Broader `include/src/tests/cmake/bindings` snapshots were also retained:
`build/task-5s-all-source-begin.sha256` SHA256
`f52ad154701fbeaebe1bc0d59feb4109ad7069c6164e9493b4dd234c0ad72836` and
`build/task-5s-all-source-end.sha256` SHA256
`154447965c63b876d88807c3852bb403256e464e44d54736e7d32a66174bbce0`.
All existing files matched. Concurrent foreign work added only
`tests/canonical_campaign_probe.cpp` and `tests/quadric_stratum_matrix.cpp`
between those snapshots. Both were verified absent from the assigned build
graph and actual exact-OFF child graph. No foreign file was edited, staged,
or compiled by this task. Entire-worktree quiescence is not claimed.

Actual exact-OFF child graph inspected after the registered dependency run:
`build/task-5c-no-faster-20260909/spin-2014-dependency-contract-runs/run-dcda002eade83cb8/enabled/build.ninja`,
SHA256 `35d278ac03b3f666727e257103514d258bf3aa98a00fab6122156595ea4df2f3`.
Its generated cache says 2014 ON, exact OFF, Faster OFF, downloads skipped.
Both new implementation names and the shape test are absent from that graph.
`build/task-5s-exact-off-inspection.log` records the actual cache/graph evidence,
SHA256 `bb160d6e7b086ad5251533705a2161ccadb8f3105ff82a19d6d2c44e73666bad`.

## Independent controller evidence

Controller independently authored/linked its exporter and raw-coordinate
NumPy QR/complex-step checker without production-operator reuse. As reported
by controller, both actual fixture exports repeated byte-identically and all
10 mode cases passed all 112/113 columns, residual/vertex gradients, handle
addition and value-only checks. Full joint max relative error `1.314e-10`;
shape-only vertex/handle max relative `1.836e-11`, near-zero absolute zero;
36/36 corrupted exports were rejected. Raw artifacts are retained at
`build/task-5s-controller-9h3NTi`; export SHA256
`7d28edb3ce9125a008d11d81a112e4b6ec6b0e72bb207f5bc3b4e1a21fd46563`.
This is controller-provided independent evidence, not implementer review.
The implementer did not edit or run the controller helpers.

Stable core (used by controller and final build):
`6450f4e751ca0b3a61c5c4a7ab9c4f22cc9193715f2fd96a6cb320a101ded854`.
No production source changed after this hash was notified to the controller.

## Owned files and self-review

Six new component/test/doc paths plus existing paper header/source and CMake
are the nine implementation paths in the brief. The report is separate.
Self-review checked formulas/sign/layout, local thin-SVD memory, source BBW
chain, value-only behavior, disabled/operator identities, overflow/error types,
and the scoped diff. No solver integration or new dependencies were introduced.
Remaining limitations are the explicit source-authority/model/legality
exclusions above; no additional correctness concern was found in self-review.

## Exact final build invocation

```sh
bash -o pipefail -c 'ulimit -v 25165824; cmake --build /data/xvxinyao/moment-guided-spin-v1/code/make-it-balance-paper-sim/build/task-5c-no-faster-20260909 --target spin-it-2014-source-cage-shape-regression spin-it-2014-source-cage-paper-regression spin-it-2014-source-cage-paper-allocation-probe spin-it-2014-source-cage-moments-regression spin-it-2014-source-conforming-boundary-regression spin-it-2014-source-cage-binding-regression spin-it-2014-source-bbw-joint-regression spin-it-2014-source-cage-bbw-regression spin-it-2014-source-cage-fem-regression spin-it-2014-source-cage-grid-regression spin-it-2014-source-surface-cage-regression spin-it-2014-source-material-regression spin-it-2014-source-adaptive-regression spin-it-2014-source-surface-octree-regression spin-it-2014-surface-value-regression spin-it-2014-adaptive-domain-regression spin-it-2014-fill-functional-regression spin-it-2014-audited spin-it-2014-slqp-regression spin-it-2014-end-to-end spin-it-2014-modes-regression spin-it-2014-givens-regression spin-it-2014-objectives-regression spin-it-2014-resolution-ladder spin-it-2014-dual-material-regression spin-it-2014-api-regression spin-it-2014-surface-differential-regression spin-it-2014-cage-functional-regression spin-it-2014-cage-solver-regression spin-it-2014-automatic-cage-regression _make_it_stand -j3 2>&1 | tee build/task-5s-build-final.log'
```

## Commits and handoff

Implementation commit: `f6e5212` —
`feat(spin2014): add source-bound linearized rotation shape functional`.
Its index was inspected immediately before commit and contained exactly the
nine owned paths, with `git diff --cached --check` clean. Identity was the
authorized `volcano <xuxinyao307@gmail.com>` for author/committer.
This report is committed separately with subject
`docs(spin2014): record Task5s shape functional verification`; no controller
helper, public status/gap document, foreign file, or whole folder is staged.
The controller owns independent review and any later acceptance decision.
