# Author-material request checklist — NOT SENT

This is a draft only. Sending it, purchasing models, using non-public accounts
or redistributing received assets requires separate authorization.

We are independently reproducing the computational methods of Spin-It (2014)
and Spin-It Faster (2024) in C++/Python. We distinguish computational
reproduction from author-code bitwise equivalence and physical experiments.
Could you share the following materials, identifying permitted usage and
redistribution terms?

## For every original case

1. Final posed/edited outer mesh or CAD document and downloadable version,
   with units and coordinate frame, rather than only the source-model link.
2. Fixed material and optimizable domain as geometry, shell thickness and
   density per material, with exclusion/overlap conventions.
3. Target axis/origin, contact point or support polygon, gravity direction,
   and all manually applied transforms.
4. Objective/constraint weights and application parameters, distinguished
   from implementation-specific solver defaults.
5. Final digital geometry, numerical references and plot/table data, including
   unsuccessful and degenerate comparisons, where available.

## Spin-It 2014

- Ellipsoid, Heart, Armadillo hand/shell, Teapot, Cuboid, Woven Ring,
  Dancing Couple, Sumone, T-Rex and Elephant final input versions.
- All three Teapot objective/weight configurations.
- Octree geometry and optimization levels, boundary subtree and shell rules,
  refinement/merge/rounding settings and SLQP controls.
- Complete cage extraction level/handles, finer-grid BBW settings,
  deformation regularizers and Dancing Couple hollow-only/joint results.
- Sumone's low/high-density components, densities and assembly transforms.
- Matched T-Rex Make-It-Stand input and configuration.
- Table 1 final cavities and manual cutting planes/parts. We do not presume
  that the paper used an automatic cutting algorithm.

## Spin-It Faster

- Original CAD/mesh inputs and domains for vase, infuser, floating island A/B,
  Heart, Noodles, Jellyfish and Flower, including discretization settings.
- Floating structures' sealed/floodable parts, prescribed waterline/fluid
  density and A's optimized material used while designing B.
- Flower's planar boundary loops, holes, bearing ring and layer dimensions.
- Figure 11 original three-stacked-egg CAD/final surface, assembly transforms,
  fixed shell/design domain, units, support/target origin and shell thickness;
  exact sampled weights and plotted solutions if available.
  We already have the published objective `-kappa*v+(1-kappa)*c3`, lateral
  first-moment constraints, kappa in [0,1], and both pure-objective endpoints.
  We also have the Section 6.1 printed-PLA density of 1187 kg/m^3 and do not
  assume three material layers or unpublished layer densities. Please flag
  only any case-specific departure from those published material conventions.
- Figure 12 original Heart specification and exact plotted sampling/output;
  Figure 13 cube/domain geometry, critical-distance numerical reference and
  shared unperturbed solution. We already have its published nine distances
  `d_m=d_crit-10^(-m)` and the instruction to share the unperturbed solution.
- Figure 10 trajectories and Table 1 full random-run/restart/failure logs,
  with function-evaluation definitions and available plotting data.
- Quadric discretization details and, if shareable, the modifications to the
  public Cherchi Boolean implementation used for the work.
- In particular, exact rank-two/paraboloid and cross-type source construction
  and derivative conventions: the supplement's explicit finite/ideal formulas
  use A inverse, which is unavailable for an exactly rank-two quadratic part.
- Any additional conventions required to reproduce the section 6.6 examples.

If some data cannot be shared, please identify them so those original
experiments remain explicitly blocked rather than being replaced by synthetic
inputs. Physical samples and manual spin/water-tank measurements are outside
this computational-only request.

Published-fact reconciliation: see
[the 2026-10-05 source correction](faster_published_protocol_correction_20261005.md).
The pre-correction draft is retained in
`snapshots/ef60354bb4ae69ba4ad2727819bf0f54cfd2c3ba2bf35b54074e60cf428d980b.md`.
