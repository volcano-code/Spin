# G7 evidence freshness after planar Faster additions — 2026-09-21

The existing two-run reconstructed G7/O7 optimizer evidence and its sealed
geometry remain preserved for historical diagnosis. They are not current
whole-tree source-bound acceptance after the planar Faster CMake additions.

The candidate G7 spec at
`build/spin2014-g7-direct-audit-candidate-20260921/run/spec.json` binds
`CMakeLists.txt` SHA-256
`e77cbf50d201a5127853110ec7dacc5a4adfe4b561d0f3acc0e6ba521995d528`.
The current file is
`6b4b9c8e15cfc68576796dba052ffb32065dd3195f3a6b61f46b1bba4ff97274`.
The strict identity comparison therefore fails. No historical receipt or
artifact was edited, and the old results are not relabeled as current.

This does not establish that the G7 numerical output changed. The planar
sources and test added in this session are gated on
`MAKE_IT_STAND_BUILD_SPIN_IT_FASTER`, but this whole-file hash comparison does
not isolate every difference from the historical CMake file. It does establish
that strict repository source identity differs. Until the runner's
dependency-closure semantics are independently verified or G7 is rerun from
the current source under its original resource envelope, use those results
only as preserved historical diagnostics. The separate same-component
feature-thickness audit remains `UNKNOWN` regardless.
