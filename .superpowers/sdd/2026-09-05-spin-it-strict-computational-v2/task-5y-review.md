# Task5y scoped spec and quality review

Reviewed base `8f653b7eec496d3465162f255dce26d83d31baed` to head
`5acbfa8c02f452c11a568f52c4e6de7575d09008` using the supplied complete diff
package, read once. Read the brief first, including the final normalized-unit
and PSD clarification, then the implementer report. No git commands, builds,
test executions, production edits, or foreign Faster/main operations were run.
Only this report was written.

## Verdict

**Spec: PASS for Task5y. Quality: APPROVE within the specified pre-seal scope.**
No actionable P0/P1/P2/P3 findings were identified in the supplied change.
This review does not grant full partition/thickness/continuous-map validation,
SealedSpinArtifact status, feasibility, physical behavior, spin, or paper
completion. Controller-owned broader tests and independent raw-export checks
remain separate evidence; their results are not claimed here.

## Substantiated checks

- Admission consumes the private immutable adaptive audit, requires completed
  outcome/last solved/no pending checkpoint, and checks terminal object identity
  and state/functional/partition hashes. The upstream adaptive solver's completed
  path already validates the unchanged checkpoint and converged guarded solve.
- Binary construction copies the actual terminal full state and changes only
  beta at the explicit >=0.5 convention. It uses actual binary evaluated points,
  retains outer connectivity, sums oriented primitive triangle references,
  cancels internal references, compacts indices deterministically, and reverses
  the carved boundary when appending it to the outer material.
- Independent long-double signed-tetrahedron integration uses the correct
  first- and second-moment identities and all ten moment components. Fixed-frame
  transformation, Givens constraint ordering, rho L^4/rho L^5 normalization,
  positive mass, and the fixed numerical PSD allowance match the brief.
- The actual geometry checker and final emitted-material ExactBoundaryAudit
  AnyNonzero both run after numerical evaluation, including ordinary constraint
  failure and explicit empty cavity. Validation requires every numerical and
  geometry gate. Exceptions/resource limits cannot publish validated material;
  an established numerical failure retains BINARY_VALIDATION_FAILED precedence.
- Candidate evidence and validated publication remain separate. Private
  constructors and const shared ownership prevent copied facts from authorizing
  material. Count/index/product planning precedes variable mesh/evaluation
  allocation; known buffers and exact diagnostics have explicit separate caps.
  Native exact arithmetic/RSS remains honestly delegated to the external guard.
- Configuration, endpoint, states, output meshes, numerical facts, and both
  geometry audits feed deterministic evidence hashes. Documentation accurately
  limits the result to actual-deformed numerical/geometric pre-seal validation.

## Focused checks outside the diff

Read only the relevant existing code to settle these interface risks:

- `SourceCageAdaptiveSolver.cpp` and its header: private endpoint ownership,
  terminal hash domains, converged final-check conditions, unchanged-rebind
  identity, and pending checkpoint lifecycle.
- `SourceCagePaperFunctional.cpp` and its header: immutable shared functional
  ownership, beta/optional-angle/all-handle layout, 2/5/6 constraint counts,
  fixed frame convention, value-only evaluation and shape scratch.
- `Spin2014GivensProblem.cpp`: static, top, and YoYo constraint ordering and
  retained-phi mixed-moment formula.
- `SourceCageMomentFunctional.cpp`: actual deformed-point indexing, binary
  functional moment order/density, outer connectivity, signed patch convention,
  and value-only BBW scratch (distinct surface integration implementation).
- `SourceCageGeometry.cpp` and its header: explicit settings coverage, actual
  state evaluation, primitive/outer relation checks and resource outcome rules.
- `ExactBoundaryAudit.cpp`: AnyNonzero permits either nonzero component sign
  while still requiring indexed closure, winding, edge/vertex manifoldness,
  no unused vertices/degenerate faces, and completed intersection-free scanning.

## Test evidence and limits

Inspected the changed regression and implementer TDD narrative; directly read
the normalization RED log and final CTest log. The latter records **1/1 PASS**,
85.48 seconds, for the focused binary-material regression. The normalization
RED demonstrates that a low-density raw residual would falsely pass until the
fixed-input normalized gate is applied. The regression covers stationary
zero/fractional/tie/one endpoints, mixed rounding failure, actual deformation
and fixed frame, cancellation, nonmanifold rejection, geometry failure,
resource/invalid settings, incomplete endpoint rejection, determinism, and
immutable lifetime behavior.

No additional test run was justified by an unresolved code risk. The new
regression's end-to-end fixtures use TopMassAndLowering; static and YoYo branches
were checked against their upstream formulas by inspection. This is a stated
coverage limit, not a discovered defect or a request to broaden this task.
