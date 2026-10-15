# Fixed-pivot critical-speed oracle — 2026-09-21

## Source interpretation

Spin-It (2014) Equation 1 cites Lewis et al., *The heavy top: a geometric
treatment*, Nonlinearity 5 (1992), 1–48. Lewis Section 2.3, Equations 2.33–2.34,
gives a sufficient upright-stability condition for maximal-axis rotation;
Section 4.1 separately analyzes linear stability. Its inertia is about the fixed
point. The inertia of the center of mass must first be translated to that point.
Failure of a sufficient condition alone is not an instability proof.

Author-hosted reference listing:
<https://www.cds.caltech.edu/~marsden/bib/1992.html>.
The 1992 PDF and transport/hash record are retained in
`build/spin-heavy-top-reference-20260921/`. HTTPS requests failed certificate
verification; the public HTTP copy was retrieved without disabling TLS checks.
Its transport is explicitly not cryptographically authenticated. SHA-256 is
`397cb7681c9c691adebf42f4820d3f5643fe370afc876488d8c17343473dfe61`.
Pages 4, 13 and 28 were visually checked, as were Spin-It pages 3–5.

This work does not change Spin-It's published central-inertia optimization
objectives, assert an error in that paper, or equate an ideal pivot with contact
on a tabletop. It separates the simulation criterion from the objective.

## Implementation and mathematical scope

`python/make_it_stand/verification/fixed_pivot_critical.py` supplies an independent
exact-rational oracle. Inputs declare an upright body with center `(0,0,h)` and
central principal moments `(Ca,Cb,Cc)`, spinning about axis 3. The analyzed
pivot moments are `A=Ca+m*h*h`, `B=Cb+m*h*h`, `C=Cc`.

Writing `T=m*g*h`, the transverse characteristic polynomial is
`lambda^4 + a*lambda^2 + b`, with

```
a = ((2*A*B-C*(A+B)+C*C)*omega^2-(A+B)*T)/(A*B)
b = ((C-A)*omega^2-T)*((C-B)*omega^2-T)/(A*B)
```

The signs of `a`, `b` and `a*a-4*b` distinguish exponential growth, distinct
imaginary roots, and marginal cases needing further analysis. Equality never
becomes a strict pass. The asymmetric sufficient bound and axisymmetric linear
bound are returned in separate fields with separate applicability conditions.

No mesh is automatically diagonalized or projected onto the axis. Geometry
alignment is explicitly unchecked. This oracle cannot admit material, certify
global nonlinear/contact stability from a spectrum, or grant paper reproduction.
All results remain conditional on the declared ideal model. A qualified final
mesh must independently establish its frame, moments and alignment before these
criteria can be connected to it.

## Executed counterexamples and tests

For the analytic fixture `m=1`, `h=.2`, `Ic=diag(.02,.03,.04)`, `omega=15`,
inserting central moments into a fixed-point bound incorrectly gives
`225 > 196.2`. The actual pivot moments are `diag(.06,.07,.04)` and the spin
axis is not maximal. The independent Euler linearization has growing modes
with real part approximately `2.96383622 /s`. This is an invented regression
fixture, not a measurement of the tea-pot or a paper case.

A second fixture, `m=h=g=1`, `Ic=diag(4,5,7)`, `omega=.7`, lies below its
asymmetric sufficient threshold but has a simple imaginary spectrum. It guards
against treating sufficient as necessary. Additional tests cover thresholds at
equality, separate symmetric behavior, independent characteristic-polynomial
agreement, spin reversal, density/size transformations and invalid inputs.

The initial seven tests failed for the missing oracle, then passed. Read-only
review reproduced them and independently matched a symbolic determinant; two
review-suggested cases were added, bringing the new suite to nine tests.

The source-bound batch `build/spin-heavy-top-reference-20260921/job/` completed
in **1.945429313 s**, peak **114,028,544 bytes**, swap **0**. It ran six analytic
configurations twice with identical results, independently compared each
polynomial against the existing Euler linearization, and executed three
registered CTest groups: nine critical-oracle, five torque-free and four
reference-frame tests (**18 tests**, no skipped groups). This contains existing
short unit-test trajectories, not new full-duration or G7 physics simulations.

Current checker SHA-256:
`9b5d827b3ab32e998ad802b2618aaa5bfb050c5804a74d9c004d67fa604fda6e`.
Test SHA-256:
`1736ac04654d74c01a5d481af3e215ffb34d41075ebbcb94c615dc75e71d2cdc`.

## Remaining gate

The G7 fixed-shell check has now completed twice, but complete same-component
feature thickness remains `UNKNOWN`. This batch does not remove that blocker.
Actual qualified-object gravity trajectories, parameter-specific critical-speed
comparisons, full perturbation/time-step matrix and validated contact proxies
remain outstanding. Neither this oracle nor the historical engineering video
grants original-case or full-paper completion.
