# Spin-It / Spin-It Faster public-material reproduction report

Computational public-complete gate: **PASS**.

- Formula coverage: `true`
- Experiment coverage: `true`
- Asset coverage: `RECONSTRUCTED_COMPLETE`
- Author-exact: `false`
- Physical validation: `false`

## Verified results

- Spin-It Faster: all eight reconstructed Table 1 cases passed the 50-start,
  at-most-eight-restart protocol. Successful-start coverage was 49/50 for the
  tilted vase, 48/50 for the tea infuser, and 50/50 for the other six cases.
  The worst successful normalized equation residual was `9.94e-13`; the worst
  physical relative residual was `4.34e-13`.
- Spin-It Faster Figures 10--13: convergence-pair data, two 41-point finite
  nondominated Pareto fronts, nine converged critical-distance rows, and both
  Section 6.6 one-dimensional degeneracy examples were regenerated.
- Spin-It 2014: all fourteen reconstructed paper-role cases passed twice at
  maximum octree level 9 (512 cubed output grid). Output hashes matched between
  repeats; every final material/cavity boundary was sealed and the maximum
  binary equality residual was zero. The dancing-couple cage deformation was
  executed and changed the silhouette while keeping valid geometry.
- The complete registered suite passed 164/164 tests with no Not Run entries.
- Representative Spin-It Faster and Spin-It 2014 bundles were opened on the
  real X11 desktop (`DISPLAY=:1`), including buoyancy waterline, displaced
  volume, metacenter and stability-margin overlays.

## Reproducibility hashes

- Faster Table 1: `0344190d59ebc11fceec32f0649e3bace2e8dc1e82150bffe44ebbccc39eee30`
- Faster Figures 10--13: `60d9fd032a5c9f2401780390322f9e43e8d2cc8f531ff6180c71a8bbce86bd83`
- Spin-It 2014 level 9: `970c342a58213c3558e64501ad28b1280b35b9a44ef72d16a84995323f5e8cbd`
- Aggregate gate: `e890ee650eae5682d0e7221daa74d80cde1c4445292dfc6606b68a9b91b11644`

## Claim boundary

This is the strongest claim supported by the published formulas and the
deterministic reconstructed assets in this repository. It is not bitwise
equivalence with unreleased author code, CAD/voxel assets, user interaction
histories, or hidden parameter files. It also does not reproduce the papers'
fabricated-object measurements. Those two stronger claims remain explicitly
false in the machine-readable release gate.

## Correction — 2026-09-05

The public-complete PASS claim above is superseded and is **not valid evidence
of complete computational reproduction**. Simplified assets, preconstructed
solutions, and an insufficient aggregate gate were used. Historical results
are preserved, with their corrected scope documented in
[STRICT_CORRECTION.md](STRICT_CORRECTION.md). The strict project status is
`PARTIAL_BLOCKED`; no original-paper experiment is completed by relabelling
these synthetic regression results.
