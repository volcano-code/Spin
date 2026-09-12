# Task5s: source-bound linearized-rotation shape regularization

## Pre-implementation API/resource ruling (approved)

LinearizedRotationLaplacian exposes immutable restVertices/faces/config,
vertexCount/coefficientCount/sparseValueCount, block-CSR rowOffsets,
columnVertices and3x3 row-major coefficients, config/operator hashes;
evaluate/evaluateValues(Ref<const MatrixXd>,shapeWeight,displacementWeight)
return weighted energies, Vx3 residual and optional Vx3 vertexGradient.
SourceCageShapeFunctional accepts actualbinding plusconfig, exposes laplacian,
binding andhashes; samevertex API adds actual4H handleGradient. hasGradient
andoutputSha256 identify each result. Paper exposes shapeFunctional pointer
(null when bothweights0), weighted energies andshapeOutputSha256.

Use custom blockCSR, S=sum_i(1+degree_i),9S stored doubles regardless ofzeros.
Unique sorted directed edges start from6F pairs. Explicit construction cap
counts simultaneously retainedrest/faces/indexarrays/coefficients/initialedgepairs
plus conservative thinSVD scratch8*3*maximumLocalRing*7+2048 doubles. Convert
integer storage to double-equivalent bytes with checked ceiling; no hidden
triplets, fullU or packing exception. maximumConstructionWorkspaceDoubles=10m,
maximumEvaluationWorkspaceDoubles=10m, maximumVertices=1m,
maximumOneRingVertices=4096, maximumSparseValues=30m. More permissive sparsecap
never overrides tighter constructioncap. Evaluation3V residual+optional3V
gradient; wrapper uses moves and separately bounds4H gradient/H row (1m each).
Fixed small block temporaries and separately boundedupstream storage are
explicitly excluded, not disguised as zero. These are implementation defaults.
Expose/check actualrequired capacities for exactcap/onebelow behavioral tests.

Paper adds shapeWeight/displacementWeight (default0), shapeOperatorConfig,
maximumShapeHandleGradientValues andmaximumShapeWeightValues. Existing state,
momentframe andgeometry unchanged. This section resolves proposed names and
resource representation below; no scientificscope or threshold changes.

BASE ce3e159a5b736951fe83daf94766c9bbd2aadfe8. This is the next approved M3
increment, not new design authority or full joint SLQP. Read this brief first.
Work only /data/xvxinyao/moment-guided-spin-v1/code/make-it-balance-paper-sim,
branch spin-it-strict-computational-v2. Main is read-only. All external Faster,
Canonical, FullRank, Plane, Polyhedral, Truncation edits/experiments are foreign.
Do not edit/stage them. No new dependencies, downloads, fallback geometry,
cache cleanup, merge, push, tag, other task or subagent. Preserve all failures.

## Source authority and scope

User's approved strict full-reproduction spec remains binding. Spin-It2014
Section6.1 adds gammaD,L * vPrime^T LD^T LD vPrime +
gammaD,v * ||vPrime-vRest||², with analytical handle derivatives. The paper
explicitly expects small rotations and cites Sorkine2004's linearized method.
Controller visually inspected original2014PDF page6 and SorkinePDF pages3–4.

- Original /home/lin/下载/SpinItPaper.pdf SHA
  3add70de03df9b0e173028d55c775e8520fa89a6b21d9b269f926ccbe8a76e4e.
- Existing extract build/strict-paper-sources/spin2014.txt SHA
  5cf5ab047189575df48b01796de73c72ceadc602a3380b74f33d2dda81459161.
- Sorkine official project PDF URL:
  https://igl.ethz.ch/projects/Laplacian-mesh-processing/Laplacian-mesh-editing/laplacian-mesh-editing.pdf
- Newly acquired build/task-5s-sources-uIE9MG/sorkine2004.pdf SHA
  246498341881005d40a27309c5aca0e3d6efd2de5d29a85f95bbcdcc978dfab4;
  text sorkine2004.txt SHA
  703ddc9364c235d6eb440920c04b04cd28f079588d258514a6bc460eb34b7095.

Implement Sorkine Section3 Eq2,5,6,8–12. Do NOT substitute ordinary graph
Laplacian, frozen differential target, ARAP polar rotations, large-rotation
iteration, cotangent variant, cage-surface regularization or cell regularization.
Uniform one-ring weights are published by Sorkine; their use in Spin-It is a
documented formula-reconstruction choice, not recovered author2014 parameters.
No original experiment or paper-completion grant, no geometry legality grant.
No solver yet; no full joint Hessian; no beta/phi-dependent shape operator.

## Required mathematical contract

Use every original source-surface vertex/triangle and all actual cage handles.
For vertex i, deterministic unique sorted one-ring Ni and local set {i}+Ni.
L_i v = vi - mean(vj, j in Ni), delta_i=L_i vRest. Fit seven coefficients
(s,hx,hy,hz,tx,ty,tz) using all vertices in that local set, equal LS weights.
For a local rest point (x,y,z), the3x7 block is
 [x,0,z,-y,1,0,0; y,-z,0,x,0,1,0; z,y,-x,0,0,0,1].
The differential transform block is the same3x4 linear part evaluated at
delta_i, followed by zero translation columns. Let D_i denote that3x7 block.
Construct residual rows R_i = D_i A_i^+ - L_i (scattered into3V coordinates).
Stack LD from these3-row blocks, xyz interleaved per vertex. Its sign is fixed
as above; energy does not care but exposed residual/digest must be consistent.

Rest local centering and isotropic normalization BEFORE factorization are
allowed and recommended: qj=(vj-vi)/restLocalScale, with delta from q; use
normalized A,D so LD remains dimensionless and is applied to original-unit
deformed vertices. Document this algebraic preconditioner; it cannot depend on
current deformation or alter physical energy units. Use rank-revealing thin
SVD/QR for the7-column problem, not explicit normal-equation inversion, not
full m-by-m U. Numerical rank must be7; reject rank deficiency/ill-conditioning
explicitly without damping, pinning, omitted rings or silent rank truncation.
Relative rank tolerance1e-12 is IMPLEMENTATION_DEFAULT, hashed/configurable;
invalid/nonfinite/out-of-range tolerances fail closed. Planar non-collinear
rings are valid and must not be rejected simply for being planar.

Compute Eshape=gammaD,L ||LD vPrime||² and Edisplace=gammaD,v ||vPrime-vRest||²
WITHOUT a half factor or density. Vertex gradient=2gammaD,L LD^T(LD vPrime)
 +2gammaD,v(vPrime-vRest). Sparse residual actions only: never form dense3V²,
LD^TLD, global VxH, Vx4H, or full state Hessian. Genuine large finite rotations
need not have zero energy. Infinitesimal fields sI+[h]cross and translations
are null modes of shape only; displacement still penalizes movement.

## Implementation surfaces and interfaces

Owned new files:
- include/spin2014/LinearizedRotationLaplacian.h
- src/spin2014/LinearizedRotationLaplacian.cpp
- include/spin2014/SourceCageShapeFunctional.h
- src/spin2014/SourceCageShapeFunctional.cpp
- tests/spin_it_2014_source_cage_shape_regression.cpp
- docs/spin2014_source_cage_shape_regularization.md
Owned existing edits:
- include/spin2014/SourceCagePaperFunctional.h
- src/spin2014/SourceCagePaperFunctional.cpp
- CMakeLists.txt: only component/test registration in2014ANDexact gate.
Do not touch SourceCageMomentFunctional, old CageSpinDeformation/solvers,
existing tests, core evidence grant code or foreign CMake sections.

Keep responsibilities separate: immutable mathematical rest-mesh operator;
immutable actual SourceCageBinding shape/BBW chain; existing paper functional
composition. Provide concise proposed public API to controller before coding
so its independent exporter uses exactly the production interfaces. A pure
vertex-evaluation API is mathematical only, never a source-validity certificate.
The source wrapper must accept actual immutable binding (not arbitrary supplied
weights); evaluate the original outer rows from the actual moment evaluation.
Chain vertex gradients using the existing centered scale convention
dvi/duh=wih*(viRest-centerh), dvi/dth=wih I, with full actual interpolated
BBW rows. Do not switch to a different parameter meaning. No material corners
contribute to surface regularization. Do not normalize/clamp weights again.

Append independently configurable shapeWeight/displacementWeight, default0,
and explicit shape/operator resource/numerical settings to paper config.
Both zero retains disabled shape behavior without constructing/evaluating LD.
Nonzero enables real operator; setting just displacement must still retain
honest mode identity. Paper result reports each weighted energy separately,
total includes both, all handle derivatives include both, beta/phi additions
are exactlyzero. Constraints/moments and their Jacobians remain unchanged.
Value-only returns identical energies with no derivative arrays or derivative
budget requirement. Do not recompute or mutate geometry to add regularization.
Source/data/frame/state/settings/rawoutput hashes bind the new mode and operator;
no post-hoc booleans or fabricated scientific acceptance. Expose enough source
geometry/operator identity for independent verification. Immutable copies and
changed source/BBW identity must remain distinguishable even under zero motion.

## Resource and invalid input contracts

Use positive explicit caps for vertices, one-ring size, sparse coefficient
contributions/nonzeros, construction workspace, evaluation workspace and
all-handle gradient/weight row. Proposed engineering defaults:1million vertices,
4096 ring vertices,30million sparse values,10million workspace doubles,
1million weight doubles. You may tighten/widen exact internal accounting before
implementation with a documented controller ruling, not silently later.
Validate integer/index/byte products BEFORE allocation. Account for simultaneous
triplets+sparse storage/factorization work, thinSVD temporaries and matrix
packing. Document what separately bounded upstream memory is excluded. Use
fixed-size local contractions or explicitly bounded local scratch. The
Task5r GEMM packing bug is the relevant prior finding: hidden temporary buffers
are not exempt. No O(V²), O(H²), O(VH) retained storage or full sparse Gram.
EvaluateValues must work when only derivative capacity is insufficient.
All invalid geometry/config/state dimensions/NaN/Inf rejected; repeated face
indices/zero-area triangles/isolated vertices and invalid local rank rejected.
Any broader manifold/self-intersection check not implemented must remain an
explicit exclusion. Allocation/size/cap issues use length_error RESOURCE_LIMIT;
nonfinite arithmetic overflow_error; invalid inputs invalid_argument. Do not
turn any of them into infeasibility or silently return zero energy.

## Tests first and independent verification

New regression and registration first, actual expected RED, then implementation.
Focused iteration only; final once all scoped related tests serial. Required:
1. Actual nontrivial rest mesh: rest, globaltranslation, uniformscale and
   infinitesimalrotation shape nullspaces. Independently computed LS values,
   no expected values from production LD. Nonuniform shear/local deformation
   positive energy; large rotation explicitly nonzero on suitable asymmetry.
2. Rotate/translate/scale rest and deformed meshes consistently: same normalized
   residual/expected scale² energy and covariant gradients. No density effect.
3. Planar noncollinear local patch valid; rank-deficient ring invalid. NaN/Inf,
   bad faces, zeroarea, isolated vertex, malformed shapes and bad budgets.
4. Both actual source fixtures from Task5r (D5,cage2,grid3,8beta,26handles), all
   fivepaper modes, nonzero weights, all113/112 columns. Independently reconstruct
   local transforms using LS on deformed points then energies/derivatives. All
   four parameter types/all handles tested, beta/phi unchanged by shape term.
   Explicit additive energy and unchanged moment/constraint tests. Zero weights
   preserve values/derivatives; changing just one weight affects only its term.
5. Analytical handle and vertex gradients: at least3 fixed FD steps, expected
   smooth-region relative<=1e-6 and nearzero absolute fixed<=1e-10 of inputscale;
   reference scales fixed from input, never final errors. SVD is construction
   only, so final energy is quadratic; centered derivative FD has no truncation
   term beyond numerical roundoff. Failures retained and diagnosed, no relaxing.
6. Sparse nnz/local-ring growth/budget fail-beforeallocation, value-only no
   gradients, caps boundary and highvalence rejection. Actual behavior not grep.
7. Identity/hashes different when appropriate; same inputs repeat2x exactly;
   source wrapper cannot accept fakeBBW; no state/geometry mutation.

Controller independently builds a NumPy LS-per-ring/checker/exporter, verifies
both actual fixtures and all joint variables and rejection of tampered evidence.
Do not write the controller checker or perform the controller review yourself.
Notify stable API/core hash early. Small fixtures remain SYNTHETIC_REGRESSION.

## Build, ownership and report

Assigned build build/task-5c-no-faster-20260909, canonical -S/-B,2014ON/exactON/
FasterOFF/LIBIGL_SKIP_DOWNLOAD=ON. max-j3, ulimit-v25165824, no downloads/pins.
Baseline source paper+allocation2/2PASS22.77s, corec5fbb311... atBASE.
Final covering command:
timeout3600 ctest --test-dir build/task-5c-no-faster-20260909 --output-on-failure
 -j1 -R '^(spin-it-2014-|spin-2014-|source-spin-material-python$)'
Maintain source quiescence during suite; do not rebuild unrelated foreign code.
Inspect actual exactOFF generated child graph to ensure component excluded.
Report exact registered count/time, no skip/NotRun, actual source begin/end.
Compiler jobs/test processes are yours; controller won't run shared builds.
Controller only independently links against stable core when notified.

Use apply_patch; preserve failurelogs and inputs. Before commit inspect index,
stage ONLY the nine owned paths (or explicit approved subset) and commit with
already authorized volcano <xuxinyao307@gmail.com>. Never stage whole folders.
Report .superpowers/sdd/2026-09-05-spin-it-strict-computational-v2/task-5s-report.md
with formulas/source boundaries, TDD RED/GREEN rawlogs, all commands/results,
hashes/resources/limits and concerns; separate report-only commit. No manual
completion grant. Return short DONE/DONE_WITH_CONCERNS/BLOCKED/NEEDS_CONTEXT,
commits/testsummary/concerns/reportpath. No subagents. Independent task review
will be dispatched by controller; fix loops resume you, not a fresh worker.
