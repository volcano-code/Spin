# Spin-It Faster Figures 10--13

This directory is generated from the public formula implementation. The
authors' Onshape CAD models, Matlab program, random seeds, and raw plotting
tables are not public, so these rows remain `RECONSTRUCTED_ASSET` and
`author_exact=false`.

- Figure 10 consumes successful attempts from `run_paper_table1.py` and plots
  the Eq. 14 residual at iteration `n` against iteration `n+1` by application.
- Figure 11 sweeps `kappa` in `min -kappa*v+(1-kappa)*c3`, including both
  endpoints, and warm-starts each point from the previous quadric.
- Figure 12 sweeps the published spinning-top center-height weight and records
  every Euclidean quadric-type transition.
- Figure 13 first brackets and binary-searches the reconstructed tilted-cube
  feasibility boundary to `1e-9`, then evaluates `d_crit-10^-m` for `m=1..9`.
- Section 6.6 is reproduced analytically for the two one-dimensional
  `p=0` examples and never mislabeled as a numerical mesh experiment.

The runner emits data, definitions, evidence class, and independent claim
flags. A figure is not a passed reproduction until its result section is
populated by converged runs and its audit fields pass.
