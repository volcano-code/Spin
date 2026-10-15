# Spin-It 2014 online trace implementation

Approved design: the user approved compatible, opt-in online trace storage for
the reconstructed teapot/source-hollowing path, followed by cage integration.
This plan implements and verifies the source path without changing optimization.

## Global Constraints

- Work only in the existing isolated paper-sim worktree; preserve prior dirty
  changes, failed runs and frozen binaries. No commit, tag, push, cleanup or merge.
- Defaults and legacy diagnostic APIs remain compatible. Scientific values,
  event ordering and native scientific hashes must be identical across storage
  modes. Storage paths and timing must not enter scientific hashes.
- No tolerance, objective, rounding, restart or acceptance changes. This is not
  a claim that the teapot binary validation failure has been solved.
- Disk mode is explicit. Native code cannot call Python while the GIL is released.
  Use exclusive creation; write/read errors must fail closed and preserve partial
  evidence. No hidden fallback to full in-memory history.
- Keep the previous build frozen. New build: build/spin2014-online-trace-20260914/native.
  Compile at most four jobs; experiments single-threaded and <=24 GiB.
- Complete iteration records go to disk during optimization. The mutable last
  record must remain pending until finalState is known. Preserve the exact existing
  scientific hash preimage, including final state and KKT data.
- Account honestly for all remaining memory, particularly adaptive pass partitions
  and result vectors: bounding iteration history is not bounding all solver memory.

### Task 1: Native disk-backed iteration history and source-hollowing integration

Read these files first: include/spin2014/SlqpLbfgsSolver.h,
src/spin2014/SlqpLbfgsSolver.cpp, src/spin2014/SlqpAuditTrace.cpp,
src/spin2014/SlqpAuditJson.cpp, include/spin2014/SourceAdaptiveSpinSolver.h,
src/spin2014/SourceAdaptiveSpinSolver.cpp, bindings/SourceSpinMaterialBindings.cpp,
tests/python/test_source_trace_stream.py.

Implement an opt-in disk-backed SLQP iteration store, integrating it into the
actual source material binding and its lazy diagnostics iterator. Preserve the
old vector/default APIs. Use a shared record traversal/hash serializer so old
and new modes hash precisely the same records. Hold at most a pending last
record plus reader/writer buffers for iteration history in disk mode. Durable
history must include every field, not just reduced scalar traces. Re-reading
and full legacy diagnostics may intentionally materialize data only if explicitly
requested by the caller. Report trace counts and storage counters accurately.

The source solver's existing per-pass partitions can remain resident in this
first increment, but document and measure them separately. Do not claim the
whole adaptive solve has constant memory; pass-metadata spooling is subsequent
work. Cage paths do post-solve phi mutation, so reject unsupported disk mode
there rather than silently corrupting records; old cage default stays unchanged.

Storage must be created exclusively, validate read lengths and format versions,
detect corruption/truncation and missing completion. Make write failures
observable and preserve partial data. Prefer a narrow native record-store module
with RAII, bounded record allocations and explicit seal/count metadata. Do not
create a generic persistence framework. Scientific trace hashing may use a
streaming readback (bounded memory) if incremental hashing requires risky broad
changes; record that I/O cost. Never silently upgrade partial records to complete.

TDD: first add failing real tests, save RED output, implement and save GREEN.
Cover disk/memory exact numeric and hash equality on deterministic native/source
fixtures, final record finalState, failure result, multiple adaptive passes,
iterator ownership/repeat reading, corruption/truncation, existing path collision,
write/resource failure and bounded retained iteration memory. Register tests in
CMake. Preserve the old golden slqp trace hash test and old JSONL archive reader.

Use a new build, same options/dependencies as the frozen native-tangent build.
Own native source/header/binding/CMake and focused tests. Parent owns independent
experiment harness, provenance and documentation. Do not edit parent files.
Do not dispatch agents. No git commits. Return a report with exact RED/GREEN
commands/output paths, files changed, interface, remaining limitations and review
concerns. Snapshot current touched files before modifications for scoped diff.

### Task 2: Independent comparison and resource evidence

Parent prepares a separate experiment harness that runs the same small reconstructed
teapot/source fixture using frozen legacy and new disk modes in subprocesses.
Save source/binary/input hashes, raw records, RSS/time and exact comparison. Keep
the fixture identified as reconstructed/synthetic and never paper-original.
Run old golden hashing/archive tests and new focused tests, plus relevant source
solver regressions. Capture partial/error states instead of calling them success.
Do not launch the 1920-iteration large batch until small validation and review pass.

### Task 3: Independent review and handoff

Create a scoped diff against the pre-task file snapshots (not HEAD's dirty delta).
Review Task 1 with a fresh reviewer for spec and code quality, fix findings through
the implementer, then review the complete increment. Preserve review artifacts.
Publish a concise report separating engineering completion from unresolved binary
rounding, full pass-memory, original asset and paper reproduction gates.
