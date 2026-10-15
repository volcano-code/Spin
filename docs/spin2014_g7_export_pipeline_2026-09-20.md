# G7 terminal trace-export diagnosis and repair

Status: exporter change implemented and focused tests passed; no new G7 solve
was launched in this diagnostic step.

## Finding

The 2026-09-19 G7O7 run remains immutable and incomplete. Its native solver
returned `BINARY_MOMENT_VALIDATED`, but the worker then serialized a multi-GB
JSONL trace, reparsed every record while writing, and attempted another complete
archive readback before saving its mesh arrays and stage receipt. The one-hour
job ended during that post-solve path. The JSONL has no completion footer, so it
is not a complete archive and cannot be upgraded after the fact.

Each solve already writes an authenticated `SLQPTRC1` history: every iteration
record is framed and checksummed, and a footer binds the record count and end
offset. The existing native JSON iterator consumes this reader and emits all
iterations in order. The safe optimization is therefore to retain the complete
JSONL output and its v1 footer, but avoid reparsing every native-generated,
multi-megabyte JSON payload during the write and avoid an immediate redundant
full JSON readback in the time-critical native mode. The ordinary caller-data
writer and complete `read_archive` reader remain strict and unchanged.

## Change and evidence boundary

`write_trusted_native_archive` accepts only the compiled
`NativeSourceSpinMaterialResult`; it checks the complete native event order,
pass/iteration indices, declared counts against native trace-storage facts,
record-size and newline constraints, and recomputes payload/archive SHA-256
while writing. The native iterator verifies each sealed binary trace record as
it is consumed. The resulting JSONL retains every prior event and the same
`source-trace-archive-v1` completion footer.

The `COMPLETE_NATIVE_STREAM_FOOTER_DIGEST` run status is narrower than a full
JSON readback: offline verification binds the complete output-file hash to the
writer digest and checks the header, event count, and terminal footer. Full
record parsing remains available and required for independent semantic audits;
the new status does not grant stage-level rounding, sealed geometry, simulation,
or reproduction acceptance. The standard dictionary writer still performs a
full archive readback in the run path.

The existing run-specific `campaign.py` and its evidence under
`build/spin2014-g7-current-source-20260919` were not modified. A subsequent
from-initial-state batch must use a newly source-bound runner that adopts the
trusted native writer, persists terminal geometry before optional post-run
audits, and does not perform a second full parser pass inside the one-hour
solver envelope. It must still preserve the complete trace and supervisor
artifact hashing. Any deferred semantic audit must be completed and recorded
before geometry is sealed or simulation is run.

## Verification

Focused tests passed:

- 48 trace-stream, archive, and native disk-reader tests.
- 26 online-comparison and run-manifest tests.
- `py_compile` and `git diff --check` on the four touched implementation/test
  files.

The fast-path tests compare every decoded event against the native dictionary
iterator, verify complete record counts and the footer/payload/file digests,
and require the input, trace, and output scientific hashes to remain unchanged.
The footer-tamper test fails closed without invoking the full archive parser.

This is an engineering I/O-path repair only. The reconstructed teapot remains
`RECONSTRUCTED`; the old incomplete G7 job remains `RESOURCE_LIMIT_WALL` /
`INCOMPLETE_EXTERNAL_WALL_TIMEOUT`, with no standalone mesh seal, independent
geometry seal, determinism repeat, simulation, or video.
