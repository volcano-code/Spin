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
