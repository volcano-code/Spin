# Strict asset discovery checkpoint — 2026-09-06

Status: **PARTIAL_BLOCKED**. This is a source-discovery record, not experiment
readiness or a complete search of every public/private archive. The old
`NOT_PUBLICLY_RELEASED` assertions in `manifest.json` are superseded by the
narrower finding **NOT_OBTAINED_AND_AUTHENTICATED**. Failure to find an input
does not prove that its authors never released it.

## Sources checked

- [Spin-It 2014 official project](https://igl.ethz.ch/projects/spin-it/):
  paper, supplement, video and bibliography; model and posing credits.
- [Faster official record](https://research-explorer.ista.ac.at/record/17203):
  published paper, supplement and video.
- [Faster project page](https://visualcomputing.ist.ac.at/publications/2024/SIF/),
  reached from [Hafner's publication list](https://chrishafner.github.io/):
  the listed resources link to those documents/video, not an input CAD bundle.
- Local targeted filename inventories in the downloads directory, external
  Make-It-Stand reference directory, and the two project checkouts found
  historical generated fixtures and Make-It-Stand beta meshes. No final
  Spin-It posed/edited input or Faster CAD bundle was authenticated.

This is not an exhaustive inspection of differently named archives. Additional
model-source/CAD repository investigation remains open. No email has been sent,
no model purchased, and no non-public account accessed.

Follow-up inspection of the Faster PDF's embedded URL annotations (`pdfinfo -url`)
found ORCID and bibliographic links, but no Onshape document or input-CAD download
link. Its text describes modeling in Onshape without identifying a document
version. This finding is limited to that authenticated PDF; it is not a claim
that no public CAD document exists elsewhere.

The 2014 acknowledgments identify useful leads, not ready experimental inputs:
Carlo H. Sequin's Sculpture Generator 1 for Woven Ring; De Espona / Robert Sumner
for Elephant, with a posed model credited to Alec Jacobson; TurboSquid for T-Rex;
Stanford for Armadillo, with poses from the cited skinning demo; and Martin
Newell's Utah Teapot. No generator settings, licensed purchase, or final posed
surface was recovered merely from these credits. They cannot upgrade the rows
below to `AUTHOR_ORIGINAL`.

## Per-case readiness gaps

`MISSING` means not obtained and authenticated for this experiment, not
necessarily absent everywhere. An analytic case may instead qualify through
`PUBLISHED_EXACT_SPECIFICATION`, once all defining data are extracted and checked.

| Original computation | Current strict input status | Required unresolved evidence |
|---|---|---|
| 2014 asymmetric Ellipsoid; Cuboid yo-yo | Exact published specification not yet established | Dimensions, transforms, axes, shell/design domain and critical parameters |
| 2014 Heart; Teapot (three strategies); Woven Ring; Dancing Couple; Sumone; Elephant | MISSING final inputs | Exact topology/pose, units, shell/domain, density, support/axis and per-case weights |
| 2014 Armadillo hand and shell | MISSING posed inputs | Final posed surfaces, contact point/axis; an unposed source model is insufficient |
| 2014 T-Rex comparison | MISSING matched specification | Identical input, pose, support, shell, density and finest scale for both solvers |
| 2014 Table 1 cavities/cuts | MISSING original digital cuts/reference geometry | Manual planes, final parts and cavities; counts must be independently computed |
| 2024 vase; tea infuser; Heart; Noodles; Jellyfish; Flower | MISSING original inputs | CAD/mesh, fixed shell/domain, units, pose/axis, densities and application parameters |
| 2024 floating island A/B | MISSING original assembly | Geometry, floodability/waterline, transforms, and A's final material as B's fixed input |
| 2024 Figures 11–12 | MISSING egg/Heart design domains | Original layers, coordinates and exact objective parameter definitions |
| 2024 Figure 13 | Exact published specification not yet established | Literal cube/domain geometry, constraints, transforms and reference definitions |
| 2024 section 6.6 | PUBLISHED_EXACT_SPECIFICATION extracted; not yet computationally verified | Two literal one-dimensional problems recorded in `faster_section_6_6_spec.json`; algorithm/evidence checks remain open |
| 2024 Figure 10 | No authenticated original-case trajectory | Original input and complete residual/evaluation logs above the numerical noise floor |

The external beta directory
`/home/lin/xvxinyao/reference-only/make-it-stand/beta-extracted/make-it-stand_src_beta/data/meshes/`
contains Armadillo and Dinosaur meshes, but their identity with the final Spin-It
posed inputs has not been established. They are source leads, not
`AUTHOR_ORIGINAL` Spin-It evidence, and remain outside this repository.

Old generated ellipsoids, low-resolution occupancy, known-solution fixtures and
same-species meshes cannot fill these gaps. See [the correction](../STRICT_CORRECTION.md).

## Storage and document retention

Both checkouts now reside under `/data/xvxinyao`; about 1.4 TiB was available
at resumption, so the earlier storage blocker is lifted. Old `/home` paths are
compatibility symlinks, not additional copies. No full-resolution run is
authorized merely by having storage: input and algorithm gates still apply.

The prior `/tmp` paper cache no longer existed on resumption. New document
copies go in ignored `build/strict-paper-sources`; official URLs and known
SHA-256 values remain in `manifest.json`. Local 2014 reference remains
`/home/lin/下载/SpinItPaper.pdf`. PDF availability is not model availability.

The renewed Faster paper and supplement hashes match the existing manifest.
The 2014 supplement from
`https://igl.ethz.ch/projects/spin-it/SpinItSupMat.pdf` has SHA-256
`f0acc779410bacac74ce9fbce305fe6dd9e1ab570885141644abac805a1396e5`.
The section 6.6 extraction was checked against rendered pages 78:12–78:13,
not inferred from an optimization failure. Mathematical uniqueness there is
understood up to equality almost everywhere, as appropriate for integral data.

The next external coordination step, if authorized, is the
[author-material request draft](author_material_request.md).

## Public acquisition follow-up (2026-09-06)

The user confirmed that no separate original assets were available and requested
public acquisition. This authorizes public downloads, not contacting authors,
buying models or accessing private accounts. No such external action was taken.

Rechecked the [2014 project page](https://igl.ethz.ch/projects/spin-it/),
[2024 project page](https://visualcomputing.ist.ac.at/publications/2024/SIF/),
[ISTA record](https://research-explorer.ista.ac.at/record/17203), public author
homepages/repository lists, and title/model/CAD searches. No complete original
experiment model/pose/domain/CAD package was obtained or authenticated. This
does not establish that none exists publicly or that the authors cannot share it.

The 2014 acknowledgments explicitly credit the skinning demo for Armadillo
poses. Following the [official FAST page](https://igl.ethz.ch/projects/fast/)
recovered the following source leads in ignored `build/strict-author-assets/`:

| Item | Download/source identity | Result |
|---|---|---|
| FAST repository | `https://github.com/alecjacobson/FAST.git`, commit `74852c04200be5ff2423f7381f8d1ae6718c185c` | Shallow clone without submodules; only README, LICENSE and Armadillo example files materialized. No author code executed. |
| Armadillo mesh | `FAST-74852c0/examples/armadillo-disconnected-ears/mesh.mesh`, SHA-256 `abcc4667c3237d26299b9a2c7907459ee62d10444194b66faa00229144c20f1a` | Source mesh contains 47,162 vertices; not an authenticated Spin-It final pose. |
| Armadillo skinning weights | Same directory, `weights.dmat`, SHA-256 `3ab0ce560d8ea247443ef75b733a0f2de81454e8e20c58a16d635b9141f1351f` | Weights plus multiple `.bf` pose/bone files and `.rbr` settings recovered. No mapping to paper hand/shell cases established. |
| Selected-file checksum list | `FAST-armadillo-SHA256SUMS`, SHA-256 `236bff467d5539cd2bb70990430b630fd775821264efa3a4f196fe547bddcb04` | Paths relative to `FAST-74852c0`; generated from materialized non-hidden example files. The pinned Git tree additionally identifies all original blobs, including `.DS_Store`. |
| Historical official demo archive | `https://igl.ethz.ch/projects/fast/skinning.zip`, local `fast-skinning-2013.zip`, 61,952,724 bytes, SHA-256 `a2dec891e7755a0496748db2939040420a24098b805ab600997082a9e5d35d71` | ZIP CRC test passed, 505 entries; 252 have parent-traversal prefixes. **Not extracted or executed.** Listing contains no Armadillo/Elephant-named assets. Local filename is a label, not verified publication-date evidence. |

These files are `SOURCE_MODEL_ONLY` leads relative to Spin-It; **no original
case is upgraded or admitted**. The repository's top-level LICENSE identifies
MPL 2.0 for covered software; this is not an independently established blanket
redistribution license for every third-party model. Downloads stay separate
from tracked implementation and are not re-licensed or published here.

Missing information still includes the exact final hand/shell surfaces or pose
mapping, target frame, design domain and experiment parameters. Faster's original
Onshape document/version and CAD/domain inputs remain missing. The existing
author-material request is still **NOT SENT**.

## Official-resource recheck (2026-09-22)

The official pages were rechecked after the current strict-method evidence
refresh.  The [Spin-It 2014 project page](https://igl.ethz.ch/projects/spin-it/)
still lists only paper, supplemental material, video, and BibTeX under its
downloads.  The [Spin-It Faster project page](https://visualcomputing.ist.ac.at/publications/2024/SIF/)
still lists only paper, supplementary document, and submission video.  The
maintained [ISTA Visual Computing code index](https://visualcomputing.ist.ac.at/code/)
lists nearby projects with explicit source-code links but contains no Spin-It
Faster entry.  No author code, original model, CAD document, pose, or complete
case parameter package was obtained from these public official sources.

This is a time-bounded discovery result, not a claim that no public copy exists
elsewhere.  It does not alter any `MISSING` row, admit a reconstructed object,
or authorize author contact, purchase, or private-account access.

## Source-model geometry audit (2026-09-22)

The ignored local acquisition cache `build/strict-author-assets-20260922/` was
checked directly after the inventory above. Its Sumner Elephant archive was
also independently fetched over HTTPS and copied, byte-for-byte, into the
project source tree at
`sources/external/sumner_elephant_poses/elephant-poses.tgz`. Archive SHA-256 is
`91ed4fa5c6b5a8c41a66b12c7af9e459577ed3e0d75b5eb1c73f4a09a782cd3c`; all 13
tar entries were listed and contain no absolute paths, parent traversal,
symlinks or hardlinks. Its README identifies these as Elephant poses shown at
SIGGRAPH 2004, says the reference Elephant originated from De Espona, and
describes the target poses as deformations of that reference. This establishes
the archive's source lineage, not a match to Spin-It's final Elephant pose.

All 11 extracted OBJ files were read-only audited without repair or remeshing.
Each contains 42,321 vertices and 84,638 triangles, one connected component,
watertight and consistent winding, and Euler characteristic 2. The 10 posed
meshes and the reference have distinct geometry; units and Spin-It pose mapping
remain unknown. A separate signed-tetrahedron integrator was compared against
Trimesh's mass-properties path for all 11 files. Maximum relative differences
were `9.92e-16` in volume, `5.70e-16` in center of mass and `3.24e-15` in central
inertia, including its off-diagonal terms. This is an independent geometry /
moment-method validation of a source model, not a Spin-It optimization or paper
case.

The local official Séquin generator archive
`build/strict-author-assets-20260922/woven-ring/SculptGen.zip` has SHA-256
`13c184b71e9d587d85334710d6a5dbcd198bcd8b1ec8eae0a0ee0c9511a8cc7a`; its 202
members have no unsafe paths or duplicate names. It contains the old Windows
generator executable, Tcl UI and demo files, but not the WovenScherk mesh; the
executable was not run. The official [WovenScherk specification](https://people.eecs.berkeley.edu/~sequin/SFF/spec.wovenscherk.html) publishes
generator values (branches 2, storeys 7, height 1.25, flange 0.80, thickness
0.06, rim bulge 0, warp 720, twist 630, azimuth -45, mesh tiles 5, texture
tiles 1, detail 12) and a generated bounding box. The official FTP download was
attempted once without personal credentials but the connection failed while
reading; no partial STL was retained. These values concern the published
WovenScherk source model and do not yet prove that its surface is the same input
as Spin-It's Woven Ring case.

Readiness remains unchanged: Elephant is `SOURCE_MODEL_ONLY`; WovenScherk is a
published source-model specification lead; Spin-It Elephant and Woven Ring
original-case rows remain `MISSING` until the exact pose/design region and
critical experiment parameters are authenticated. The Elephant archive and
extracted meshes are now retained under the project source tree with hashes and
attribution in `sources/external/sumner_elephant_poses/SOURCE.md`; this is a
source-model lead, not an original Spin-It input.

## Cross-model source-surface algorithm check (2026-09-22)

The recovered models were next exercised at the source-surface geometry entry
point, without selecting a paper pose, support, design region, wall thickness,
or objective. Thus these are model-ingestion/moment checks only, not optimizer
runs or original-case evidence.

For `elephant-01.obj` (SHA-256
`4706710f5490136fd5265ce88f26405fcba71a7f5cf65d7fce0344bd6dca449c`), the
ordinary topology checks still report a single watertight, consistently wound
component with Euler characteristic 2. However, the strict exact-boundary
checker found **794 intersecting triangle pairs**. The source-octree entry
therefore rejects this mesh. No repair or remeshing was attempted; Elephant is
not accepted as a solver input on this evidence. The earlier independent
mass-moment comparison does not establish intersection-free geometry.

The FAST Armadillo file
`build/strict-author-assets/FAST-74852c0/examples/armadillo-disconnected-ears/mesh.mesh`
(SHA-256 `abcc4667c3237d26299b9a2c7907459ee62d10444194b66faa00229144c20f1a`)
is a tetrahedral volume-mesh file, not a surface-only OBJ. Its `Triangles`
section has 86,482 faces referencing 43,243 of the 47,162 listed vertices. The
3,919 non-surface volume vertices were excluded and the triangle indices
remapped; no surface face or vertex coordinate was altered. The extracted
surface has one component, is watertight, consistently wound, Euler
characteristic 2, and passed the strict exact-boundary test with zero
intersecting pairs. The independent signed-tetrahedron integrator measured
volume `0.06861869382742591` and center
`(-0.0026482301, 0.0384252246, -0.0396363874)` in the source mesh's unspecified
units. A source-surface octree/domain inspection at maximum depth 3 was accepted
(377 source cells, 2 adaptive leaves, 1 edge); its ten surface moments agreed
with the independent integrator to maximum relative error
`1.2351e-13`. The inspection took about 11 minutes end-to-end; this measured
cost argues against spending the next batch on a higher-depth example run before
profiling the shared exact-boundary/source-surface path.

The Armadillo result is `SOURCE_MODEL_ONLY`: FAST describes an Armadillo
skinning example, but this mesh/pose is not authenticated as either Spin-It
hand/shell input. The paper's support frame, cage layer, design region, and
optimization parameters remain missing. No optimization was run, and the
original Armadillo case readiness rows remain `MISSING`.

### Supplement surface-moment identity check on FAST Armadillo (2026-09-22)

Using the unmodified surface extracted from the FAST Medit mesh above, the
native Spin-It supplement surface-moment routine was compared against the
independent signed-tetrahedron integrator. The maximum scaled discrepancy over
the ten raw moments was `3.47e-16`. Analytic translation, dilation, and rotation
covariance checks had scaled residuals `4.72e-16`, `1.83e-15`, and at most
`3.08e-16`, respectively. The native Jacobian shape was `10 x 129729`.

This is useful evidence that the surface moment values and these analytic
identities agree on a large, non-axis-aligned source mesh. A subsequent
finite-difference check covered every Jacobian column. This still is not a
not a cage/BBW or optimization result. The loaded native module's SHA-256 is
recorded in `fast_armadillo_surface_algorithm_identity_audit_2026-09-22.json`,
and the extension was freshly rebuilt from the current worktree (99 build
steps, followed by a no-op rebuild). The worktree remains dirty, so this is not
an immutable clean-release attestation; therefore the result remains exploratory
algorithm evidence, not strict release evidence.
The asset remains `SOURCE_MODEL_ONLY`, and the original Spin-It hand/shell
readiness rows remain `MISSING`.

The centered finite-difference check covers all 129,729 coordinate columns
(43,243 vertices × 3 axes), comparing the native analytic Jacobian to an
independent signed-tetrahedron moment integrator. It computes only the changed
incident-triangle contributions for each perturbation; no columns are omitted.
At step sizes `1e-5` and `3e-6` times the bounding-box diagonal, maximum scaled
discrepancies were `1.84e-13` and `5.90e-13`; the two finite-difference results
differed by at most `6.16e-13` under the same scaling. This is a full Jacobian
check for this one source mesh, not a proof across every topology or solver
state. Exact scaling, distributions, and hashes are recorded in the linked
JSON audit.

## Full Sumner Elephant pose-family exact-boundary check (2026-09-22)

The official archive was staged under
`sources/external/sumner_elephant_poses/` without altering its bytes; the
archive SHA-256 matches the prior acquisition cache. All 11 source surfaces
(the reference and poses 01–10) were independently passed through the native
`audit_exact_triangle_boundary(..., "AnyNonzero")` checker using the pinned
strict build's `_make_it_stand` module
(SHA-256 `cbc6de4ea6cc27b81e17ea1315f70d699015b1461f305e2b90a798865895ba09`).
All audits completed. Every mesh remains closed, edge/vertex-manifold, and
consistently oriented, but every mesh also has exact self-intersections:

| Source mesh | Intersecting triangle pairs |
|---|---:|
| `elephant-reference.obj` | 658 |
| `elephant-01.obj` | 794 |
| `elephant-02.obj` | 1,369 |
| `elephant-03.obj` | 1,213 |
| `elephant-04.obj` | 1,302 |
| `elephant-05.obj` | 1,453 |
| `elephant-06.obj` | 1,263 |
| `elephant-07.obj` | 1,838 |
| `elephant-08.obj` | 918 |
| `elephant-09.obj` | 2,157 |
| `elephant-10.obj` | 1,247 |

The check made 6,335,169 exact candidate-pair tests in 182.734 seconds across
all 11 meshes. Full per-mesh records, including input hashes and timings, are in
[`elephant-source-exact-boundary-audit-retry1.jsonl`](../../../build/strict-faster-methods-refresh-20260922/elephant-source-exact-boundary-audit-retry1.jsonl)
(SHA-256 `2e5bc425f2fccde931721ff2fb159dfb1c35835a5f2513ebc6708e696bafffc9`).
No mesh was repaired or used for optimization. This closes the source-family
audit only: it does not prove these meshes match Spin-It's Elephant input, nor
does the rejection establish that the paper's algorithm fails on the original
input. The original Elephant case and its readiness remain blocked.

## Public asset recheck — 2026-09-23

The official Spin-It page still advertises only the paper, supplemental
material, video, and BibTeX under its downloads; the Faster project page still
advertises the paper, supplemental document, and submission video. The public
model leads were checked again and recorded in
[`public_asset_recheck_2026-09-23.json`](public_asset_recheck_2026-09-23.json).

The Sumner/Popovic Elephant archive remains a retained `SOURCE_MODEL_ONLY`
asset. Its archive hash is unchanged and its source page identifies the
reference Elephant as originating from De Espona, but this does not identify a
Spin-It pose or design domain. The Séquin WovenScherk page exposes a published
source-model specification and generator parameters, but not the Spin-It
edited input; the cited FTP mesh was not obtained. Neither result upgrades a
strict original-case row. No author contact, purchase, or private-account
access was attempted.

## Faster original-case admission baseline — 2026-09-30

The exhaustive machine-readable readiness inventory
`case_asset_readiness.json` contains 40 original-case rows across both papers.
Its Faster subset has 19 rows: 17 are `MISSING`, while the only two
`PUBLISHED_EXACT_SPECIFICATION` rows are the explicitly stated one-dimensional
Section 6.6 examples. Those two exact mathematical examples do not supply
geometry or critical parameters for any Table 1 case or Figures 10–13.

The 17 blocked Faster rows separately identify their missing artifacts and
critical parameters. They cover the eight Table 1 applications, the
figure-specific inputs, the chained floating-island material, CAD round-trip
assets, and the flower-spinner 2D/extrusion inputs. No reconstructed geometry
is promoted to fill these rows. The inventory SHA-256 at this review is
`07f88d0499b27a9a99a5d749c5382c8ae441a5a54b712af5540b370a60850211`.

This is an asset-admission baseline only. It does not change implementation or
verification status and does not assert that public searching is exhausted;
any newly acquired, authenticated author asset must update its own row with
source, bytes, and parameter evidence before the corresponding paper experiment
can be admitted.
