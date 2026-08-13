# R10 experiment workspace

This directory is reserved for license-safe, versioned R10 manifests and
derived evidence.  The schemas and recording workflow are documented in
`docs/r10_manifest.md`.

Do not commit third-party beta models, absolute local paths, fabricated paper
parameters, or placeholder measurements.  A publishable result directory must
contain a canonical manifest, one canonical run record for every declared
repeat, SHA-256-bound derived artifacts that may legally be redistributed, and
the canonical aggregate produced by `python/r10_aggregate.py`.
