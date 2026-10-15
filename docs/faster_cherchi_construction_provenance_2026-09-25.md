# Cherchi construction provenance sidecar — 2026-09-25

## Scope

The pinned Cherchi reference probe now emits a diagnostic
`<output>.construction.tsv` alongside its exact-coordinate and face-label
sidecars when `--dump-exact-rational` is selected. The new sidecar preserves,
for each canonically ordered output vertex, its Cherchi point kind and the
explicit point operands used to construct it (`EXPLICIT3D`, `LPI`, `TPI`, or
`LNC`; LNC also records its interpolation parameter).

This is provenance evidence only. It is not a derivative, a Boolean-backend
integration into the Faster evaluator, or evidence that the Cherchi and CGAL
outputs represent the same solid.

## Validation

Test-first development added assertions to
`experiments/cherchi_reference_probe/test_probe.py`. Before the implementation,
the regression failed because the construction sidecar was absent. After the
implementation, the pinned Cherchi probe rebuilt successfully and the complete
synthetic box intersection/subtraction regression passed:

- exact sidecar: 16 output vertices; 14 LPI and 2 explicit vertices;
- LPI line/plane operand coordinates exactly match original input vertices;
- each LPI plane's three defining points exactly match an input facet;
- output face-source recovery: 28 records, no ambiguous source faces in this
  box fixture;
- intersection volume and moments match the analytic box overlap;
- subtraction volume and moments match the analytic difference;
- two runs per operation produce identical output hashes;
- two exact-diagnostic runs produce byte-identical OBJ, exact-coordinate,
  construction, and face-label sidecars;
- ordinary output OBJ bytes are unchanged by diagnostic sidecar generation.

Observed hashes:

- intersection: `8dbd38f84340440199158cec9d9535ce7474748e9a8440a17cf598512257a037`
- subtraction: `fb3819e2c502eada4f494729f03aafb59823964b5c7cde5c85d5cd5aaf94985e`

The isolated component test ran with the single-threaded probe at
`build/cherchi-reference-singlethread-20260917/cherchi_probe`; evidence is under
`build/cherchi-construction-sidecar-determinism-20260925`.

## Remaining blocker

The sidecar records geometric construction operands, not stable source-triangle
IDs or a complete motion Jacobian. In particular, the box test's unique LPI
plane source is a favorable diagnostic fixture. Coplanar duplicate source
facets, topology changes, and the actual triangulated quadric inputs still
require explicit uniqueness and derivative checks. Cherchi/CGAL solid
discrepancies on several quadric types remain unresolved. The Faster solver
must therefore continue to reject this diagnostic output as a strict
Newton-value/Jacobian pair until the actual material map, provenance, derivative,
and independent solid-equivalence gates are all satisfied.
