# Carving reachability classification

This report classifies reconstructed orientations against a fixed 
safe-region continuous-carving relaxation. `PROVEN_UNREACHABLE` does 
not rule out deformation, another support region, discrete manufacturing 
choices, or changing the external pose.

| Case | Prior optimizer outcome | Reachability | Cells | Certificate residual | Classification |
|---|---|---|---:|---:|---|
| dinosaur-head-contact | SUCCESS | PROVEN_UNREACHABLE | 13526 | 6.11e-10 | PHYSICALLY_UNREACHABLE_UNDER_CONTINUOUS_CARVING_RELAXATION |
| spheres-45deg | FAILED | PROVEN_UNREACHABLE | 51724 | 7.78e-10 | PHYSICALLY_UNREACHABLE_UNDER_CONTINUOUS_CARVING_RELAXATION |
| toy-head-down | FAILED | PROVEN_UNREACHABLE | 4179 | 0 | PHYSICALLY_UNREACHABLE_UNDER_CONTINUOUS_CARVING_RELAXATION |

Deterministic classification complete: `true`.
Local reachability release gate: `true`.
