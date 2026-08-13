# R10 evidence summary

> Scope: Validated local evidence from the modern implementation; not an exact paper reproduction.
> Paper-reference values are reference-only and are never pooled with local measurements.

## Integrity

Validated the integrity of 9 evidence case(s), their canonical manifests, complete run-record sets, canonical aggregates, and SHA-256-bound artifacts; this check does not establish scientific correctness.

## Local measurements

All rows below are `local_measurement` values from repository-defined synthetic or local protocols.

| Case | Variant | Executed | Algorithm success | Targets reached | Manufacturing ready | Target energy | Final cavity fraction | Material mass change | RMS displacement | Total time (s) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| synthetic-carving-only-r16-h0 | carving-only | 5/5 | 1.000 | 1.000 | 1.000 | 1.4251882e-32 | 0.021657689 | 0.021657689 | 0 | 0.11094113 |
| synthetic-deformation-only-no-laplacian-r16-h0 | no-carving-no-laplacian | 5/5 | 0.000 | 0.000 | 0.000 | 1.0303597e-16 | 0 | 0.018228411 | 0.012804142 | 0.12597819 |
| synthetic-deformation-only-r16-h0 | deformation-only | 5/5 | 0.000 | 0.000 | 0.000 | 1.1988563e-07 | 0 | 0.031488782 | 0.013394774 | 0.1467202 |
| synthetic-full-r12-h0 | resolution-r12 | 5/5 | 1.000 | 1.000 | 0.000 | 5.7007526e-32 | n/a | 0.022612936 | 3.2049378e-17 | 0.058082111 |
| synthetic-full-r16-h0 | full | 5/5 | 1.000 | 1.000 | 1.000 | 1.4251882e-32 | 0.021657689 | 0.021657689 | 0 | 0.11194625 |
| synthetic-full-r16-h1 | handle-high | 5/5 | 1.000 | 1.000 | 1.000 | 7.7037198e-32 | 0.021657689 | 0.021657689 | 2.2662333e-17 | 0.18431546 |
| synthetic-full-r20-h0 | resolution-r20 | 5/5 | 0.000 | 0.000 | 0.000 | 2.8582311e-13 | 0.020686498 | 0.020731303 | 1.9427951e-05 | 0.39536266 |
| synthetic-no-laplacian-r16-h0 | no-laplacian | 5/5 | 1.000 | 1.000 | 1.000 | 1.4251882e-32 | 0.021657689 | 0.021657689 | 0 | 0.11195496 |
| synthetic-suspension-r16-h0 | suspension | 5/5 | 1.000 | 1.000 | 1.000 | 0.00065705495 | 0.23268194 | 0.05242604 | 0.21633039 | 0.24327566 |

## Paper-reference values

No numerical paper-reference metric is declared by these manifests.

## Interpretation boundary

This report summarizes the supplied local regression, ablation, and robustness evidence only. `manufacturing_ready` denotes a sealed software geometry audit, not printing success. The report does not establish exact reproduction, object-class generality, or physical-print validation.
