# Task5h preflight

Plan: docs/superpowers/plans/2026-09-05-spin-it-strict-computational-v2.md, approved user M0–M8 specification remains binding. BASE59ed96c12e8f24582f14459056c05f74d9c15492, existing linked worktree on strict branch, not submodule, clean status. Main tracked diffSHA93e32af47ab1957a1efcb6716437656be129053f4269d2169378fcd95b05a35e unchanged. /data1.4TiB free. Relevant baseline source-cage/source-octree2/2 passed3.35s. No original-case/run/release gates changed.

| Shared interface / task | Producer vs consumer / finding |
|---|---|
| Task5g cage ->Task5h fine domain | Coarse selected cells/all handles define domain; finer source-only classification would discard outside-source handles. Explicit refine-cage-volume rule below. |
| Task5h ->later BBW | Exposes fine cells/nodes/edges/handle map; no graph-BBW claim. FEM or other audited discretization can consume cells later. |
| Task5h ->adaptive source solver | Query existing/new cell corner locations through same immutable interpolator; no fake dense upsampling optimization. |
| Task5h ->legacy GridDeformationDomain | Keep legacy untouched; no full root arrays/nearest-node substitutions in new path. |
| Task5h ->final seal/geometry moments | Computational cage volume is not material; interpolation outputs don't grant convergence or seal. |
| Task5h own tests/implementation | Analytic multiaffine exactness/closed-point lookup/sparse counts and independent bounds; no solved weights/original cases asserted. |

Ruling: refine the already validated coarse cage volume uniformly to the fine level rather than reclassifying only source-intersecting fine cells — keeps every coarse handle and matches the Figure4 fine grid covering the cage; precise domain choice is IMPLEMENTATION_DEFAULT, not recovered author metadata — if author details differ, adapt the domain and rerun downstream weights, without changing original-case acceptance requirements.

Ruling: implement domain plus general trilinear field transfer before BBW minimization — separates actual grid/lookup defects from discretization/solver correctness and enables independent polynomial tests — costs an additional solver-integration step; does not reduce full BBW/joint requirements.

Controller has read original2014 PDF complete page6 (Section6/Figure4/6.1), via PDF skill and existing rendered page. User's approved plan already authorizes these prerequisites; this is continued execution, not a new speculative design request. Implementation delegated per SDD, independent reviewer after commit. Controller parallel useful work: derive independent exact counts/lookup examples and audit subsequent BBW discrete contract without modifying production.
