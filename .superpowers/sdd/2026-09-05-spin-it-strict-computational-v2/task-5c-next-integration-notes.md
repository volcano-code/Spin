# Follow-on integration obligations (not implemented by Task5c)

Task5c supplies immutable mesh-sourced geometry, not a replacement for all lattice APIs. The following remain explicit blockers to full Task5:

1. Adaptive optimization leaves must query the new oracle while preserving original exterior mass. For an approximate carveable subset D_h of actual design D, `s_final=s_outer-sum(beta_k*s_cell_k)` leaves unresolved boundary material present; do not subtract exact D wholesale and then refill only D_h.
2. Separate geometry sampling depth, optimization leaf depth and binary exported geometry. Actual geometry depth9 does not establish any optimization variable ever reached depth9.
3. Split and merge must preserve a disjoint cover and mass, remap beta, rebuild sparse face adjacency and reset/update L-BFGS state. Existing lattice merge is a one-pass sibling grouping constrained by initialLevel and near-machine bound snapping; it is not evidence of full paper merge semantics.
4. Existing binaryVoidMask consumes source_->cellCount/grid/spacing and cannot accept null source. Maximum-resolution dense export is not an appropriate bridge. Need explicit sparse material boundary output with conformity at mixed-level cells and independent sealing.
5. A binary carved-cell union can create disconnected fragments or nonmanifold contacts, even if every cell is inside valid design geometry. Geometry oracle validity alone never grants SealedSpinArtifact or FEASIBLE_WITNESS.
6. Conservative design restriction may lose feasible candidates. Solver failure/LP infeasibility on the restricted model cannot imply infeasibility on the original continuous design.
7. The whole unresolved-cell volume bounds missed volume only. Moment-error bounds require fixed coordinate extrema per monomial, and optimizing with approximate moments still needs final geometry reintegration.
8. Full cage from all intersecting cells, finer-grid BBW, geometric derivatives and guard rollback stay separate. The value-only integrator saves memory but does not implement differentiated deformed cut cells.

These are future acceptance obligations, not an authorization to expand the current sole-writer task or claim original experiment progress.
