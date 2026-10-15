# Faster restart-event ledger closure (2026-09-24)

This is a runner/evidence correction, not a paper-case result. The reconstructed
Table 1 runner previously retained Newton iteration traces but omitted the
native event stream for each externally scheduled restart. That lost failed
evaluation reasons and prevented an independent check that the recorded
function-evaluation count matched the actual solver calls.

The runner now serializes the complete deterministic native event fields into
each attempt and validates the attempt ledger before summarizing a case:

- every expected top-level run has recorded attempts;
- restart ordinals are unique, contiguous, and stop at the first convergence;
- each one-initial-guess native call has one start and one terminal event;
- `evaluation` event count equals the recorded function-evaluation count;
- each attempt remains within the strict 60-evaluation budget;
- unsuccessful runs retain all registered restart attempts.

Event details participate in the scientific result hash. Timing fields remain
excluded. A unit test first failed because the validator did not exist; after
implementation the complete runner unit/native file passed 8/8. Current-source
CTest also passed `spin-faster-python-api`, `spin-faster-table1-miniature`, and
`spin-faster-paper-figures-miniature` (3/3, 29.56 seconds).

The first independent-process replay then found that aggregate
`median_residual_seconds` and `median_derivative_seconds` were still included
in the scientific hash. A regression test reproduced that failure before the
hash filter was corrected to exclude all `*_seconds` fields. After the fix,
two separate processes ran the frozen `static-tilted-vase` miniature with one
start and one restart each. Both returned `CONVERGED`; each event ledger held
44 events, including 42 evaluation events matching the recorded 42 function
evaluations. Their scientific hashes matched:
`2502022eb005299c2e28e3beb9e753da1427574e2cd893ebc58c6392540cf402`.
The JSON file SHA-256 values differ because timing values remain in the full
payload; they are not scientific identity. Files are retained in
`build/faster-restart-ledger-replay-20260924/`.

Scope remains limited: these tests exercise a reconstructed miniature runner
and synthetic figure fixtures. They do not constitute the 800 original-input
Table 1 runs, author-source equivalence, or a completed Faster reproduction.
