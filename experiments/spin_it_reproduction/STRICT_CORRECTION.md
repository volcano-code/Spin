# Strict reproduction correction — 2026-09-05

The historical `PUBLIC_COMPLETE_REPRODUCTION`/PASS claim at commit
`e1190569294a8853bc12a86c6fd9d31b2428de05` is withdrawn as evidence of full
paper reproduction. Those tests establish properties of simplified local
fixtures, not the authors' complete computational experiments.

## Historical evidence retained, not rewritten

| File | SHA-256 of preserved bytes | Correct scope |
|---|---|---|
| release-gate.json | 6f4f9b9c3a4a20202485fa78b2afbb58e9a42cf50815f15c11561d4f01df4c27 | SUPERSEDED_CLAIM |
| faster/table1-results.json | 83d8d70b27bdb41927933a290ae8b76de9ebfa71dc8935544f6baf0524af113b | SYNTHETIC_REGRESSION |
| faster/figures/results.json | 8e152080dbf9b40dbf55fab7507b71d679c35e7406fd100b40cef44dd9c71d06 | SYNTHETIC_REGRESSION |
| spin2014/results-level9.json | a92e6ec38ec6e7197973a6366d931ca87749de99f5fa1421764ef35d57b49f15 | SYNTHETIC_REGRESSION; trex-static is KNOWN_SOLUTION_VALIDATION |

The full earlier tracked state remains recoverable from the historical commit.
In particular, the 16-cubed generated source lattice is not original level-9
geometry. The woven ring and dancing couple are ellipsoids. The synthetic
T-Rex cavity was constructed to balance on complete removal and initialized
at that known feasible solution. These facts prohibit promotion of those
results to original-paper evidence, even when residuals and hashes pass.

The legacy JSON and historical logs remain unchanged for auditability. New
reports must link this correction and use the strict requirements catalogue.
Current whole-project status: **PARTIAL_BLOCKED**. Physical validation:
**NOT_IN_SCOPE**. Author bitwise equivalence is not established.

## Working-tree preservation

Development occurs only in the linked `make-it-balance-paper-sim` worktree.
The main checkout's initial `git diff --binary` digest was
`93e32af47ab1957a1efcb6716437656be129053f4269d2169378fcd95b05a35e`.
That digest covers the tracked diff, not untracked files; it is not a claim
that a complete backup of the main checkout was made. No main-checkout files
are authorized for modification by this implementation checkpoint.
