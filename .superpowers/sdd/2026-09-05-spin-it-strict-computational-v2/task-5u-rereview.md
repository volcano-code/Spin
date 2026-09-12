### Finding Verdicts

- **Invalid-motion tests did not distinguish primitive self-intersection from negative orientation** — **ADDRESSED**. The cases now declare expected categories at `tests/spin_it_2014_source_cage_geometry_regression.cpp:180`. The helper reconstructs the same SHA-matched actual output, verifies self-intersection through completed exact scans with nonempty intersecting pairs, verifies inversion as an intersection-free single negative-volume component with the exact nesting reason, and verifies protrusion by its deterministic subset reason at `tests/spin_it_2014_source_cage_geometry_regression.cpp:39`.

### New Breakage in the Fix Diff

None.

### Out-of-Scope Observations

- The deferred production-formatting Minor remains open.
- Exact-contact and positive-volume primitive-pair-overlap end-to-end coverage gaps remain open; this assertion-only fix does not address them.
- The native allocator limitation remains unchanged.

### Checks

- Retained discriminating RED confirmed category substitution fails at `.superpowers/sdd/2026-09-05-spin-it-strict-computational-v2/task-5u-minor2-red.log:4`.
- Retained focused GREEN passed 1/1 in 17.91 seconds at `.superpowers/sdd/2026-09-05-spin-it-strict-computational-v2/task-5u-minor2-green.log:4`.
- Artifact hashes match the current test, production sources, core library, binary, and logs. No tests were rerun.

### Verdict

**Fix round: All findings addressed, no new Critical/Important breakage.**
