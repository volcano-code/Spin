# R10 experiment evidence contract

R10 is a paper-alignment and experiment phase, not an exact-reproduction
claim.  Every experiment is represented by three immutable document types:

1. `make-it-stand-r10-experiment-manifest-v1` fixes the claim level,
   protocol, metric definitions, known deviations, commit, and SHA-256-bound
   inputs before execution.
2. `make-it-stand-r10-run-record-v1` binds one repeat to the manifest and
   records the build, non-identifying hardware description, timing scope,
   local measurements, and SHA-256-bound outputs.
3. `make-it-stand-r10-aggregate-v1` records median, first/third quartiles,
   IQR, minimum, and maximum across the complete repeat set.

The JSON Schemas live in `schemas/`.  Runtime validation is implemented in
`python/make_it_stand/r10.py`; it is intentionally stricter in two areas that
are awkward to express portably in JSON Schema: canonical bytes and rejection
of host-specific paths.

## Claim boundary

Choose exactly one level:

- `paper_governed`: the local protocol follows numerical rules stated in the
  paper, with every implementation deviation still disclosed.
- `configuration_anchored`: local inputs or parameters are anchored by an
  independently hashed configuration, but missing paper details prevent a
  stronger comparison.
- `local_robustness`: a repository-defined sensitivity or robustness study.
- `smoke`: execution evidence only.

`exact_reproduction_claim` is structurally required to be `false`.  Values
transcribed from a publication use `source: paper_reference` and
`comparability: reference_only`.  Measurements generated here use
`source: local_measurement`.  Aggregation never pools the two sources.

## Canonical JSON and binding

R10 JSON uses `utf8-json-sorted-keys-compact-v1`: UTF-8, lexicographically
sorted object keys, compact separators, finite JSON numbers, and no trailing
newline.  `document_sha256()` hashes those exact canonical bytes, so whitespace
cannot alter a manifest binding.  `write_canonical_json()` is the supported
writer and `read_canonical_json()` fails closed on any non-canonical encoding.

Paths stored in a record are canonical POSIX paths relative to the
caller-selected evidence root.  Absolute paths, parent traversal, symlinks,
hostnames, and working directories are rejected.  The root is used only while
hashing and is never serialized.  Input/output artifact entries include a
lowercase SHA-256 and byte count; provenance resources also require SHA-256.
The protocol names exactly one `native_extension` provenance resource, with a
byte count, and every run record copies that binding exactly. The local driver
populates it from the ordinary, non-symlink `_make_it_stand` extension file it
passes to the model wrapper. This binds the executed binary's bytes but is not
a cryptographic source-to-binary attestation: the declared commit, compiler,
and flags still depend on the recorded build provenance.

## Recording a repeat

```python
from make_it_stand.r10 import RunRecorder, bind_artifact, write_canonical_json

recorder = RunRecorder(
    manifest,
    run_id="spheres-full-r00",
    repeat_index=0,
    seed=31000,
    build={
        "build_type": "Release",
        "compiler": "GNU",
        "compiler_version": "recorded-by-driver",
        "commit_sha": "<git-object-id>",
        "dirty_worktree": False,
        "flags": ["-O3"],
    },
    timing_scope="model_setup_start-through_run_record_finalization",
    warmup_runs=1,
)
with recorder.phase("optimization"):
    result = run_one_fixed_configuration()

recorder.add_local_metric(
    "projected_com_error",
    result.projected_com_error,
    semantics="maximum normalized signed support violation at final state",
)
recorder.add_artifact(
    bind_artifact(
        evidence_root,
        "runs/r00/diagnostics.json",
        artifact_id="diagnostics-r00",
        role="diagnostics",
        media_type="application/json",
    )
)
write_canonical_json(
    evidence_root / "records/r00.json",
    recorder.finish(status="succeeded"),
)
```

Run-record timing uses a monotonic performance counter and the fixed scope
token `model_setup_start-through_run_record_finalization`: timing starts just
before model setup and ends when the run-record fields are finalized. For
publishable timing, the experiment manifest should additionally fix Release
mode, thread count, CPU policy, warm-up count, and repeat count; the record
captures hardware and build values.
The local driver defaults to five repeats and one discarded warm-up per repeat,
fixes OMP/BLAS thread-count variables to one, and discloses that affinity and
CPU-governor state are uncontrolled. Native phase timing is a separately
SHA-256-bound, successful-return-only sidecar; it is not inserted into the
outer non-nested phase list.

When a run record is checked against its manifest, its build commit and dirty
state must exactly match manifest provenance, its warm-up count must match
`protocol.algorithm.warmup_runs_per_measured_repeat`, and every local metric's
unit and semantics must exactly match the declared metric definition.

Aggregate a complete repeat set with:

```text
python python/r10_aggregate.py \
  --manifest experiments/r10/<experiment>/manifest.json \
  --run-record experiments/r10/<experiment>/records/r00.json \
  --run-record experiments/r10/<experiment>/records/r01.json \
  --output experiments/r10/<experiment>/aggregate.json
```

The CLI requires canonical inputs, validates every manifest hash and metric
definition, requires repeat indices `0..repeats-1`, and writes canonical JSON.
The quartiles use linearly interpolated Hyndman--Fan type 7 quantiles.

No third-party beta asset, external absolute path, or unverified paper result
belongs in a committed manifest.  Such assets can be selected through a local,
untracked evidence root and represented only by logical relative name and
SHA-256 when their license permits the intended use.

## Report layer

`python/r10_report.py` strictly revalidates an evidence tree and writes
`summary.csv`, `summary.md`, and `summary.svg`. Paper-reference rows remain
`reference_only` and never enter local aggregates or plots. The Markdown table
separates protocol execution from `algorithm_success`, `targets_reached`, and
the sealed software `manufacturing_ready` indicator. The convergence panel is
the first completed repeat's committed-event ledger, not a repeated-run
aggregate and not the paper/beta iteration unit. Report files are derived
presentation artifacts; the manifests, records, aggregates, and bound outputs
remain the evidence of record.

Two material quantities must not be conflated: `final_cavity_volume_fraction`
is cavity volume divided by the final outer envelope and is emitted only for a
strictly audited boundary; `material_mass_change_fraction` compares initial
and final material mass and may be negative after deformation.
