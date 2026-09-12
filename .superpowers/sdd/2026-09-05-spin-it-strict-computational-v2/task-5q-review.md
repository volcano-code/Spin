# Task5q independent task review

Reviewer: `/root/review_spin2014_source_cage_moments`, fresh read-only agent.
Range: `6f61946cbb2b7d3329321f6741c4667ad616042f` through
`e2499187b5d7470fbe4eb454c8d12d26a4b92419`.
Package: `review-6f61946..e249918.diff`. Scope is Task5q, not the entire branch.

## Initial verdict: Needs fixes

Critical: none. Important: one legacy API compatibility regression.

At `src/spin2014/SurfaceRawMomentDifferential.cpp:39`, the factored triangle
kernel throws `std::overflow_error` when finite coordinates produce a nonfinite
cross product. The previous whole-surface implementation threw
`std::invalid_argument` for that condition. Both public surface APIs propagate
the changed exception. Preserve their legacy category while retaining the new
source API's arithmetic-overflow contract; add category-specific tests for both
legacy entry points. Controller independently confirmed BASE lines63–65.

The reviewer found actual shared-boundary signs, density application, all-handle
rest embeddings, local chain-rule contraction, value-only operation and bounded
storage consistent with the task. Independent tetrahedral/complex-step evidence
was inspected, as were final scoped tests and exact-OFF dependency graphs.
Unknown old dependency preservation is explicitly unverified, not inferred from
the current pinned sources. No reviewer tests, mutations or subagents were used.

Minor observations, deferred without extending this fix:

- A canceled-motion but distinct-state hash regression would strengthen coverage;
  current implementation hashes both complete parameter vectors, so no actual
  identity defect was found (`tests/spin_it_2014_source_cage_moments_regression.cpp:151`).
- Historical intended RED output includes inherited CMake deprecation warnings;
  retain that output, do not rewrite it as pristine.

Original implementer resumed for one bounded TDD compatibility fix.

## Fix round 1: Approved

Same original reviewer, range `e249918` through `82d7ad1`, package
`review-e249918..82d7ad1.diff`. All findings addressed; no new Critical/Important
breakage. Spec compliance and task quality approved for this bounded layer.

- `SurfaceRawMomentDifferential.cpp:40` applies the legacy invalid-argument
  category only to nonfinite normals; wrapper `:112` supplies that policy.
  Source API keeps the private kernel's default arithmetic category. No broad
  exception catch or polynomial change.
- `tests/spin_it_2014_surface_value_regression.cpp:11` checks both legacy APIs
  and distinguishes normal from later polynomial overflow.
  `tests/spin_it_2014_source_cage_moments_regression.cpp:142` checks both source
  modes retain arithmetic overflow.
- Reviewer inspected actual RED/GREEN, focused4/4, repaired34/35 plus separate
  Python1/1 recovery and source-change audit. It verified current core hash,
  two new scientific exports, independent112-column checks and12 corruptions.
- No new out-of-scope observations, tests rerun, mutations or subagents.

Two prior Minor observations stay deferred. This is not a complete branch,
joint solver, geometry-seal or paper-reproduction review.
