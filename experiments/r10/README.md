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
the full, carving-only, deformation-only, and no-Laplacian variants.  The fixed
matrix also varies voxel resolution and handle placement and includes a
suspension case.  It writes each manifest before native model loading, records
only measurements actually returned by diagnostics, and does not read the
third-party beta archive.

Example (write evidence outside the source tree):

```text
PYTHONPATH=python:build/python/python python3.12 python/r10_local_experiment.py \
  --output-root /tmp/make-it-stand-r10 \
  --repeats 3 \
  --build-type Release \
  --compiler GNU \
  --compiler-version 13.3.0 \
  --build-flag=-O3 \
  --build-flag=-DNDEBUG
```

The output root is operational input and is never serialized.  A failed native
run produces a `failed` record with no invented metrics, allowing the remaining
matrix to finish while preserving an auditable failure count.  The matrix is
local robustness evidence; it is not an exact paper-reproduction claim.
Run-record status describes whether the protocol executed and exported its
artifacts; `algorithm_success`, `targets_reached`, and `manufacturing_ready`
remain separate measured indicators, so an unstable native result cannot be
silently promoted to a successful optimization result.
