# R10 experiment workspace

This directory is reserved for license-safe, versioned R10 manifests and
derived evidence.  The schemas and recording workflow are documented in
`docs/r10_manifest.md`.

Do not commit third-party beta models, absolute local paths, fabricated paper
parameters, or placeholder measurements.  A publishable result directory must
contain a canonical manifest, one canonical run record for every declared
repeat, SHA-256-bound derived artifacts that may legally be redistributed, and
the canonical aggregate produced by `python/r10_aggregate.py`.

## License-safe local matrix

`python/r10_local_experiment.py` runs a repository-authored CC0 overhang across
five fixed profiles: full, carving-only, deformation-only, no-Laplacian, and
no-carving/no-Laplacian. The nine-case matrix also varies voxel resolution and
one deterministic handle placement and includes a suspension case. Multiple
bases are covered by the API/invariant regressions but are not a quantitative
axis of this matrix; the handle case is not described as a noise distribution.
The driver writes each manifest before native model loading, records only
measurements actually returned by diagnostics, and does not read the
third-party beta archive.

Example (write evidence outside the source tree):

```text
PYTHONPATH=python:build/python/python python3.12 python/r10_local_experiment.py \
  --output-root /tmp/make-it-stand-r10 \
  --repeats 5 \
  --warmup-runs 1 \
  --build-type Release \
  --compiler GNU \
  --compiler-version 13.3.0 \
  --build-flag=-O3 \
  --build-flag=-DNDEBUG
```

The CLI defaults to five measured repeats with one discarded warm-up per
repeat. It fixes the OMP/OpenBLAS/MKL thread-count environment to one; CPU
affinity and governor remain uncontrolled and are disclosed in each manifest.
The output root is operational input and is never serialized. A failed native
run produces a `failed` record with no invented metrics, allowing the remaining
matrix to finish while preserving an auditable failure count.  The matrix is
local robustness evidence; it is not an exact paper-reproduction claim.
Run-record status describes whether the protocol executed and exported its
artifacts; `algorithm_success`, `targets_reached`, and `manufacturing_ready`
remain separate measured indicators, so an unstable native result cannot be
silently promoted to a successful optimization result.

The native `performance.json` sidecar is available only after a successful
native return and is SHA-256-bound to its run. The outer run timer begins at
model setup and ends at run-record field finalization, including SHA-bound
artifact registration. Neither timing unit is the paper's beta iteration
count.

Each case also copies the ordinary, non-symlink `_make_it_stand` extension file
actually selected by the driver into `provenance/native_extension.bin`. Its
SHA-256/byte-count binding is declared by the manifest and copied into every
run record, so report validation fails if that binary is missing or changed.
This identifies the executed binary bytes; it does not cryptographically prove
that those bytes were built from the repository commit or compiler flags named
by the caller.

Render a deterministic presentation report after all cases are complete:

```text
PYTHONPATH=python python3 python/r10_report.py \
  --evidence-root /tmp/make-it-stand-r10 \
  --output-dir /tmp/make-it-stand-r10-report
```

The command revalidates canonical manifests, records, aggregates, and bound
artifacts, then emits `summary.csv`, `summary.md`, and `summary.svg`. These are
derived presentation files, not replacements for the canonical evidence. The
event curve uses only the first completed repeat and its x-axis is the modern
committed-event ordinal. `manufacturing_ready` is a sealed software geometry
audit indicator, not evidence that an object was printed or physically stood.

The archived `local-synthetic-release-v1` evidence binds clean implementation
commit `4b4303679bdca9091ca8ddcedd18095d5aec2de6` and a Release native extension
built with `-O3 -DNDEBUG`. Its 45/45 measured executions completed, while the
independent scientific indicators report algorithm success for six of nine
cases and manufacturing readiness for five of nine. See the checked-in
`reports/local-synthetic-release-v1/summary.md`; execution completion must not
be substituted for algorithm or manufacturing success.
