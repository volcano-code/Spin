# Task5l: shared-edge FEM boundary cancellation fix

Scope: restore the actual non-dyadic transformed medium source cage that failed
in Task5k. Base `97d9f29`. Work in the existing isolated
`/data/xvxinyao/make-it-balance-paper-sim`; unrelated Faster edits and dirty
main checkout are preserved. No public API, mesh, tolerance, optimization
objective, resolution or handle changes. Overall `PARTIAL_BLOCKED` remains.

## Root cause and implementation

Source cube [.3,.7]^3, root origin (-.3,-.2,-.1), root side1.7, cageLevel2,
fineLevel3, vertices multiplied by1.7 then translated. The valid domain has
64 cells,125 nodes,26 handles and98 exposed boundary nodes.

The former third barycentric coordinate `1-w1-w2` evaluated to
`-2.710505431213761e-20` or `-5.421010862427522e-20` on real coarse edges.
The strict negative test excluded some nodes from both incident faces.
Diagnostic node123, key(6,6,5), was falsely excluded by triangles42 and45.
The direct opposite-edge area is zero there. The separate read-only diagnostic
found12 false-negative triangle tests; its complete coordinates are retained.

All three barycentric weights now use their own oriented opposite-edge areas.
Shared edge reversal uses the same reversed operands, and actual axis-aligned
edge zeros no longer arise from subtracting two rounded divisions from1.
Finite checks, strict nonnegative admission, reconstruction and sum checks,
duplicate-row consistency, all exposed-node completeness and budgets remain.
There is no epsilon widening, clamp, normalization, node snapping or mesh edit.
The numerical boundary-method identity is newly included in config hashes.
This remains floating-point numerical assembly, not an exact-predicate
certificate for arbitrary coordinates. Existing conservative coordinate and
resource applicability restrictions continue to apply.

## Tests and evidence

Fresh directory `build/task-5l-fem-boundary-8FCQG9/`:

- `red.log`: updated actual FEM regression against the old frozen core;
  exactly two new groups fail with incomplete boundary constraints.
- `diagnostic.log`:12 explicit barycentric false negatives using actual grid
  and cage positions; no production instrumentation.
- `green.log`: all FEM groups pass with only the new FEM object overriding
  the frozen static library. New regression checks98 rows, exactly two
  physical endpoint handles for node123, actual endpoint interpolation,
  independent topology/constraints, immutable geometry and repeat hash.
- `bbw-red.log`: updated BBW covariance regression linked to old core;
  all four1.7-scale combinations fail, preserving the original issue.
- `faster-off-build.log`, `faster-off-ctest.log`: normal project rebuild of
  four relevant targets with Faster disabled;4/4 registered tests pass7.26s.
  Existing FEM topology/mass/energy, holes, islands, uneven intervals and
  failure guards remain covered; BBW includes all18 groups plus both2and1.7
  transformed runs, both protocols and mass rules, two actual solves each.

The four1.7-scale actual BBW solves have maximal weight covariance error
`6.661338147750939e-16`, original normalized KKT at most
`3.330669073875470e-16`, and deterministic repeated scientific output.
All26handles and702coupled variables remain. Raw energies obey inverse scale.
No original-paper model or material/physics claim is added.

Commands use `ulimit -v 25165824`, compile at most-j3, timeout300 for standalone
probes and3600 for project build/CTest. The main build is intentionally not
rebuilt over ongoing unrelated Faster work; prior203/203 is historical, not a
new full-current-tree result. This task claims focused correctness only.

## Independent check and review completed

Implementation commit `dd5428c`. `task-5l-export.cpp` links the freshly rebuilt
Faster-disabled core and exports actual positions, integer keys, all handles
and sparse boundary rows. Each mass rule was exported twice; JSON bytes match.
`task-5l-exact-boundary-check.py` imports no project code: Python Fraction
constructs exact rational ratios from represented binary64 coordinates and
uses independent square-diagonal max/min formulas. It checks all98rows,
26completehandles,2548weight entries per rule and24shared-edge nodes with
exactly two endpoint handles. Maximum numerical-vs-exact error for each rule:
`1.1254693832033723e-20`. See `exact-check-{voronoi,barycentric}.log` and
`boundary-*-{1,2}.json`. The reference is exact arithmetic; production remains
floating-point and is not advertised as exact predicates.

Fresh reviewer `/root/review_spin2014_fem_edge_fix` returned scoped approval,
no Critical/Important/Minor requiring changes. It checked determinant
orientation, preserved gates, downstream method-hash propagation, actual red
and green evidence, and original-transform BBW coverage. No builds, tests or
mutations were performed by the reviewer. Review covered
`97d9f29..dd5428c`; independent Fraction evidence was completed separately
by the controller afterward.

`final-bindings.sha256` binds production source, both native regressions,
actual Faster-OFF core/test executables/HiGHS library/cache, and independent
export/checker source. All10 hashes were rechecked successfully. The main
checkout tracked-diff SHA remains
`93e32af47ab1957a1efcb6716437656be129053f4269d2169378fcd95b05a35e`.
Unrelated QuadricMomentEvaluator/QuadricHollowingSolver/test changes remain
uncommitted and untouched. No whole-current-tree clean/tested assertion.

This closes the specific Task5k non-dyadic shared-edge blocker and restores
that original fixture, not an arbitrary-coordinate proof. No whole-branch
release, original-case, joint-shape, sealed-material or physics grant. Next:
complete source-based beta/phi/handle joint optimization and its derivative
and geometry verification.
