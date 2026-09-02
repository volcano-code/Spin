# Spin-It Faster reconstructed OBJ report

All inputs are deterministic, license-free reconstructed assets. This is a
`RECONSTRUCTED_ASSET` formula reproduction, not author-CAD or author-code equivalence.
Fixed-shell moments are integrated as `outer - design domain`; the visible
outer surface is unchanged and both material and cavity boundaries are audited.

| Case | Mode | Faces | Interface | Eq. 14 residual | Physical residual | Repeat | Result |
|---|---|---:|---|---:|---:|---|---|
| static-faceted-shell | `STATIC_BALANCE` | 768 | `PLANE` | 8.361e-15 | 6.948e-16 | identical | PASS |
| fidget-faceted-shell | `FIDGET_SPINNER` | 768 | `ELLIPTIC_CYLINDER` | 2.793e-14 | 2.389e-15 | identical | PASS |
| yoyo-faceted-shell | `YO_YO` | 768 | `HYPERBOLOID` | 1.478e-13 | 8.262e-15 | identical | PASS |
| top-faceted-shell | `SPINNING_TOP` | 768 | `HYPERBOLOID` | 4.204e-14 | 2.200e-15 | identical | PASS |

Release gate: **PASS**.
Wall-clock durations are recorded but excluded from deterministic hashes.
