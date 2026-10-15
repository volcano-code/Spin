# Spin-It model acquisition and algorithm-first case audit — 2026-09-22

## Acquired public source files

Files are retained under ignored `build/strict-author-assets-20260922/` and
are not copied into tracked source or redistributed by this project.

| Source package | URL | SHA-256 | Contents and scope |
|---|---|---|---|
| Elephant poses | `https://people.csail.mit.edu/sumner/research/deftransfer/data/elephant-poses.tgz` | `91ed4fa5c6b5a8c41a66b12c7af9e459577ed3e0d75b5eb1c73f4a09a782cd3c` | Eleven OBJ files: ten target poses and one reference mesh |
| Sculpture Generator | `https://people.eecs.berkeley.edu/~sequin/GEN/Sculpture%20Generator/SculptGen.zip` | `13c184b71e9d587d85334710d6a5dbcd198bcd8b1ec8eae0a0ee0c9511a8cc7a` | Windows generator executable, Tcl source/scripts, and 20 parameter demo files; no OBJ/STL mesh was included |

The Elephant package README attributes the meshes to the De Espona library and
describes them as the Elephant poses from Sumner and Popović's 2004 deformation
transfer work. This is a strong same-source lead for the 2014 Spin-It
acknowledgment, which credits De Espona Elephant geometry and Alec Jacobson's
posed Elephant. It does not identify the exact final Spin-It pose, design
domain, contact axis, shell, or optimization weights.

All eleven Elephant OBJ files load without repair and each has 42,321 vertices,
84,638 faces, one connected component, a watertight boundary, and consistent
winding under the local mesh audit. Their individual hashes and geometric
descriptors are recorded here:

| OBJ | SHA-256 | Volume (source units³) |
|---|---|---:|
| `elephant-reference.obj` | `75f2dae53d72c85c08714a566fc027fa7dadc3b115aea5aafb48ab165102bd6e` | 0.11091942637551805 |
| `elephant-01.obj` | `4706710f5490136fd5265ce88f26405fcba71a7f5cf65d7fce0344bd6dca449c` | 0.09800008399058689 |
| `elephant-02.obj` | `8d1254c05229a17badda821a00486246d2539d447439b6db19c920eec8ebbb02` | 0.10157440007926251 |
| `elephant-03.obj` | `99ceb10ab8b970748d35bacebffaf336d5af1546557ca1dba9cbfa5d3039cccc` | 0.09771362809185184 |
| `elephant-04.obj` | `b3e47bdd000427659012ae3dafaef3d5ed6880462f664933bfdb3b3c2249c209` | 0.10354915483744374 |
| `elephant-05.obj` | `cd371d6e8d5ab12d3996fef31af4838c1d0953173996d33408f75f114611a8d0` | 0.10420596510502929 |
| `elephant-06.obj` | `57c161de6f0753ba6ee202ff27e6c2fe63fcf8303bd6d99e13a653efe4362fea` | 0.10395435304706492 |
| `elephant-07.obj` | `745b11e77a06c3e96e194659d1b8a4a46d338ab141d3faa2bc5cce24bf0fcf4a` | 0.09486254362747433 |
| `elephant-08.obj` | `dce71f8078215bade4c3e75c799a7ba26926c7949b7ab918101982cdcb28738` | 0.10423532676532632 |
| `elephant-09.obj` | `b1a08762cd568ead0fd27cb3c8f618ebcab0e522fdb97e49633cb76536f510c5` | 0.10092126681762321 |
| `elephant-10.obj` | `3d88ddb5f1c0fbd447611c6345d27b21990bf58435af0bba78018082f4ff1bf2` | 0.10725746148389802 |

No asset is upgraded beyond `SOURCE_MODEL_ONLY` for Spin-It.

The native exact represented-boundary audit was run on all eleven OBJ files
without modifying them. Exact self-intersection counts were:

| OBJ | Exact intersecting face pairs | Boundary admission |
|---|---:|---|
| `elephant-01.obj` | 794 | rejected |
| `elephant-02.obj` | 1,369 | rejected |
| `elephant-03.obj` | 1,213 | rejected |
| `elephant-04.obj` | 1,302 | rejected |
| `elephant-05.obj` | 1,453 | rejected |
| `elephant-06.obj` | 1,263 | rejected |
| `elephant-07.obj` | 1,838 | rejected |
| `elephant-08.obj` | 918 | rejected |
| `elephant-09.obj` | 2,157 | rejected |
| `elephant-10.obj` | 1,247 | rejected |
| `elephant-reference.obj` | 658 | rejected |

For `elephant-01.obj`, 674 of the 794 reported pairs share no vertex index;
the remaining 120 share one index. The checker operates on the exact rational
interpretation of the input binary64 coordinates, not a distance tolerance.
Uniformly scaling the mesh to propose a reconstructed inner design region
preserves its intersections, so this is not caused by the octree or that
design-region construction. Although all files pass ordinary watertight and
winding checks, all eleven are rejected as direct solid inputs by the exact
material pipeline. No automatic repair, welding, smoothing, or remeshing was
applied. A repaired derivative could only be a separate reconstruction
experiment with the original rejected input retained and geometric deviation
reported.

## Armadillo same-source surface audit

The locally retained FAST checkout is pinned to commit
`74852c04200be5ff2423f7381f8d1ae6718c185c`. Its
`examples/armadillo-disconnected-ears/mesh.mesh` is a tetrahedral volume file
with an explicit surface-triangle section. The volume-vertex table has 47,162
rows, but the surface references only 43,243. Compacting that surface index
set removed 3,919 unreferenced *interior volume nodes* without changing any
referenced coordinate or surface triangle. The compact surface has 86,482
faces and passes the exact represented-boundary gate: one positive-volume
component; closed, edge- and vertex-manifold, consistently outward-oriented;
zero degenerate faces; zero exact self-intersection pairs; all vertices
referenced.

The original volume file SHA-256 is
`abcc4667c3237d26299b9a2c7907459ee62d10444194b66faa00229144c20f1a`; the
exact-boundary digest after lossless index compaction is
`a7b3a10a7541aa4209cf0210db9f6862f8897ed7768418686919cae252f88336`.
This is a valid complex surface-method input, but remains `SOURCE_MODEL_ONLY`:
FAST's disconnected-ears demonstration is not evidence of the particular
Spin-It hand/shell pose or its BBW/cage parameters.

A first explicitly reconstructed interior-domain proposal (scale `0.70`
about the mesh volume centroid) was rejected by the exact volume-containment
check: it is not an inclusive subset of this non-star-shaped outer solid. The
native source-octree/optimization path therefore did not run. This is an
input-domain finding, not an optimizer failure. No scale was searched after
the rejection and no alternative domain was silently substituted.

The FAST volume mesh contains 142,685 positively oriented tetrahedra. The
combinatorial boundary of those tetrahedra exactly matches all 86,482 explicit
surface triangles (no missing or extra boundary face; maximum tetrahedral
face incidence is two). This is useful provenance/consistency evidence, but
it does not itself prove non-overlap of every tetrahedron. A separate
reconstructed design proposal selected the largest face-connected component
among tets whose centroids lie within `0.12 * bbox_diagonal` of the surface
volume centroid. That proposal also did not pass the complete boundary gate.
A face-dual BFS proposal with fixed depth 32 likewise did not pass it. Both
were discarded before optimization; they do not justify changing the geometry
checker or tuning domain parameters until one passes.

## Current algorithm-first evidence state

The registered CTest inventory contains 287 tests. On 2026-09-22, two complete
serial executions (`run-8` and `run-9`) each ran all 287 tests with zero
failures; both test-ID inventories matched. The execution logs and captured
evidence are under `build/registered-ctest-v2-20260922-run-{8,9}/`. This is a
software regression gate, not paper-case evidence.

| Area | Current evidence | Still not established |
|---|---|---|
| Spin-It mass properties / derivatives | Native implementation, supplementary-algorithm tests, independent moment and derivative test groups in the 287-test runs | Full author geometry and paper-scale values on the original cases |
| Spin-It objectives / constraints / SLQP | Separate objective paths, Givens/SLQP/KKT tests registered and executed | Full paper case traces and author parameter matching |
| Adaptive source-surface octree / binary material | Source-derived adaptive implementation and transition tests; reconstructed teapot G7 history remains preserved | Paper Level-9 full optimization; G7 is a reconstructed teapot and its complete thickness gate remains unknown |
| Cage / BBW deformation | Implemented paths and synthetic geometry/FEM tests | Original Armadillo/Dancing Couple pose, full paper cage handles and paper outputs |
| Dual density / digital cuts | Method code and material tests exist | Original Sumone materials and author cut locations/results |
| Spin-It Faster moments / quadric solver | Multiple theorem, Newton, moment, planar and application-contract regression groups execute | Full Boolean-to-discrete-moment derivative chain and all original CAD inputs/applications |
| Spin-It Faster CAD / buoyancy | OCCT round-trip and reconstructed hydrostatic/planar method tests exist | Matching original CAD, waterline/design configuration and paper experiment outputs |
| Multi-case outcome | Reconstructed Cuboid/Heart low-level runner repeats deterministically; the earlier self-built Heart source-surface case has sealed, nontrivial material removal and independent checks | These are reconstructions, not original paper cases; current Cuboid/Heart low-level batch did not show material redistribution |

The historical strict audit snapshot `build/strict-current-baseline-audit-20260922-v26.json`
reported 49 of 98 catalogue rows verified, 47 blocked/not implemented, and two
physical measurements out of scope. It is not promoted as the current result:
source/build inputs have since changed, so a fresh complete audit with
current checker evidence is still required. The current base audit remains
`PARTIAL_BLOCKED`.

The generator archive can export STL interactively and contains demo parameter
files with branch counts from 1 to 7. It does not identify a Woven Ring
parameter file, provide a direct mesh export, or establish which demo produced
the Spin-It Woven Ring. The source attribution is therefore confirmed, but the
paper geometry remains un-recovered. The Windows executable was not run.

## Algorithm-first regression batch

`build/algorithm-first-cases-20260922/cuboid-heart-rerun.json` records a
two-repeat, maximum-level-4 run of the existing Cuboid yo-yo and Heart spinning
top reconstructed fixtures. Both report `passed=true` in this reconstructed
regression runner. The report explicitly says `RECONSTRUCTED_ASSET`,
`paper_resolution_reproduced=false`, and `paper_numerics_reproduced=false`.
These results validate that the existing two example paths still execute; they
do not count as author-case reproduction or Level-9 algorithm validation.

Each case's two repeats have matching trace, output, and geometry hashes.
Both cases report zero binary-constraint residual and zero axis-alignment
error. The Cuboid's outer and material OBJ hashes are identical, so this run
does not demonstrate a nontrivial mass redistribution. The Heart's run is also
already converged at the configured low resolution. These are regression
smokes, not sufficient optimization challenges; the next model run must report
material-change magnitude and iteration count before it can count as a
nontrivial end-to-end case.

The first launch failed before entering the solver because Conda's C++ runtime
did not expose the extension's required `GLIBCXX_3.4.32`. That failed attempt
was preserved. The identical scientific configuration then ran with the
project's system `libstdc++` preload and completed. Do not combine the failed
launch with scientific results.

Both reconstructed cases remain low-resolution smoke runs. The Cuboid's outer
and material OBJ hashes are identical (no demonstrated redistribution), and
the Heart was already converged under the configured low-resolution terminal
constraints. Neither is yet a nontrivial end-to-end optimization challenge.

## Next admission steps

1. Keep this Elephant package blocked from the direct material solver unless a
   valid unmodified pose is found; all ten target meshes plus the reference
   fail exact self-intersection admission. Search for the paper's separately
   published mesh or establish a deviation-audited reconstruction route
   without replacing the original rejected assets.
2. Match an Elephant pose to the Spin-It figure/video, then record the remaining
   shell, design-domain, axis, unit, and objective-parameter gaps.
3. For Armadillo, source a certified contained design region from its existing
   tetrahedralization or separately register a new reconstructed domain before
   another exact solver run; do not tune a global scale until it happens to
   pass.
4. Recover Woven Ring generator parameters or an exported mesh from a source
   that can be tied to the paper's specific design.
5. Add an OBJ-backed source-octree and material-solve fixture only after its
   fixed-shell/design partition and input settings are explicitly declared.
6. Keep all acquired objects in `SOURCE_MODEL_ONLY` or `RECONSTRUCTED` until
   their exact role in the Spin-It experiment is established.

The complete feature-thickness requirement remains separate from algorithm
validation and remains `UNKNOWN` for the reconstructed G7 teapot. It does not
block algebraic, derivative, solver, or source-domain method tests; it does
block claiming a fully manufacturing-qualified final artifact.
