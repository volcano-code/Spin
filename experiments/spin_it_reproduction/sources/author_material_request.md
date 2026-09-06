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
- Figure 11 egg layers and parameter range; Figure 12 Heart specification;
  Figure 13 cube, critical-distance reference and shared initialization.
- Figure 10 trajectories and Table 1 full random-run/restart/failure logs,
  with function-evaluation definitions and available plotting data.
- Quadric discretization details and, if shareable, the modifications to the
  public Cherchi Boolean implementation used for the work.
- Any additional conventions required to reproduce the section 6.6 examples.

If some data cannot be shared, please identify them so those original
experiments remain explicitly blocked rather than being replaced by synthetic
inputs. Physical samples and manual spin/water-tank measurements are outside
this computational-only request.
