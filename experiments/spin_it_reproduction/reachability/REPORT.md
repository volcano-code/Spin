# Alternative-pose carving reachability report

Two independent full-resolution runs produced byte-identical JSON and report
artifacts.  All three reconstructed orientations are infeasible even under the
continuous relaxation that independently retains any eligible interior cell
fraction in `[0,1]`.

| Case | Evidence | Prior optimizer | Variables | Certificate residual | Input SHA-256 | Result SHA-256 |
|---|---|---|---:|---:|---|---|
| Dinosaur head contact | `RECONSTRUCTED_PARAMETER` + official cache | success with deformation | 13,526 | `6.106e-10` | `2c21079d3459c124a5a525b5a33c5ff42100a4109099fc910780e005e51b8a73` | `e0ed758e2d80ca6ce2fbfd02870669bd16807df79305c96426f9995bd579abc6` |
| Spheres 45° | `RECONSTRUCTED_PARAMETER` + official cache | failed | 51,724 | `7.782e-10` | `5cf8cc3583640bc47a62d700b628fbc28d69c0a28aaeb67f176ebe9ef5feadf0` | `998d8d2dbf6ede330dcca36009282496b6cdef0a969c5e6db083af75fe3470c0` |
| Toy head down | `RECONSTRUCTED_ASSET` | failed | 4,179 | `0` | `e7d8708522b1f76cff76d54358885ffaba6991d1db8f40fab942b510c782b364` | `a2943c4c425a4a99d6a3cab29b8cdf036d3ebcee40dafaa74c7ff34901deaec8` |

Each native status is `PROVEN_UNREACHABLE`; every independently checked Farkas
certificate has nonnegative weights, normalization residual below `3e-14`,
stationarity infinity norm at most `7.782e-10`, and a strictly negative
contradiction.  The corresponding report classification is
`PHYSICALLY_UNREACHABLE_UNDER_CONTINUOUS_CARVING_RELAXATION`.

“Physically unreachable” is deliberately scoped to the fixed support region,
zero-angle target, frozen external shape, and continuous carving model in
`cases.json`.  It does not prohibit deformation or a changed contact geometry;
that is why the prior Dinosaur run can succeed while its carving-only
relaxation is infeasible.
