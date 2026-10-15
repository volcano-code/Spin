# Spin-It 2014 online iteration storage — engineering evidence

Status: VERIFIED_STORAGE_ENGINEERING_INCREMENT. No paper or physics completion is granted; the overall paper reproduction remains partial.

## Scope

The source-hollowing/teapot path gains explicit online iteration storage while
the default retains its previous memory-backed behavior. Complete iteration
records, scientific hashes, solver objectives, rounding and acceptance tolerances
must agree. The previous experiment build remains frozen.

This increment bounds retained **iteration history**, not the entire adaptive
solver: per-pass partitions, fill vectors, source geometry, current numerical
state and L-BFGS history still use memory. Full pass-metadata spooling and cage
integration are separate remaining tasks. Calling the legacy full diagnostics
property explicitly materializes diagnostics and is not a low-memory operation.

## Preservation and reproducibility

- Worktree: `/data/xvxinyao/moment-guided-spin-v1/code/make-it-balance-paper-sim`.
- Main worktree is untouched; no commits, tags, cleanup or merges performed.
- Frozen native module SHA-256:
  `b3d63e5d6838a9d85e3098d84d9d221e3142ece51b9a723f6f815eea34138615`.
- New build and results: `build/spin2014-online-trace-20260914`.
- Scoped snapshots and review records:
  `.superpowers/sdd/2026-09-14-spin2014-online-trace`.
- A first baseline attempt was marked `STALE_RESULT` after the controller changed
  the checker during execution. Its archive, provenance, traceback and failure
  record are preserved under `controller/frozen-g6`; it is not accepted evidence.

## Validation completed before native integration

- Four frozen baseline CTests passed, 6.63 seconds; `controller/frozen-baseline.log`.
- Six independent comparison-checker tests passed; `controller/compare-green-v2.log`.
- The checker compares all canonical records and exhausts archive integrity checks,
  rejects identical-file comparisons, stale files and missing completion, and does
  not equate matching numerical failures with paper success. Signed zero is retained.
- The G6 engineering configuration was checked against the previous G6 experiment:
  solver, geometry, boundary, original input and minimum wall all match, with the
  same explicit `PAPER_POST_TRANSITION_V1` used by the prior runner.

## Numerical boundaries remain unchanged

The previously completed G7/1920-iteration reconstructed teapot run failed binary
validation (normalized residual approximately `4.60424e-5`, required `1e-6`).
Reducing trace memory does not repair that result. The model is reconstructed,
not the original paper asset, and no storage test upgrades its evidence grade.

## Final code verification

- Final integrated code review: READY, no remaining findings. Two earlier evidence
  binding defects and the relative-path reader ownership defect were fixed and
  independently re-reviewed. Review records are retained in the scoped SDD directory.
- Final stable-source selected CTests: **19/19 passed**, 231.26 seconds,
  `controller/final-ctest.log`. This is the selected regression group, not all
  registered project tests or a full-paper release gate.
- Independent comparison checker: **23/23 passed**, 0.30 seconds,
  `controller/final-checker-tests.log`. The online Python test group includes
  14 cases within the CTest total; do not add these as separate CTest runs.
- `git diff --check` passed on changed tracked files.

Opt-in native API:

```python
result = native.solve_source_spin2014_material(
    # existing geometry/settings arguments unchanged
    ...,
    iteration_store_directory="/absolute/new/run/iteration-store",
    iteration_store_max_bytes=8 * 1024**3,
)
for event in result.iter_diagnostics():
    # bounded-by-record traversal, without materializing result.diagnostics
    ...
storage_counters = result.trace_storage
```

The directory must be new. The budget applies cumulatively across passes.
Native record size is limited to 64 MiB; decoded allocation is separately bounded.
Files remain after result destruction or exceptions; no implicit memory fallback.
Path anchoring preserves readers after a working-directory change. Cage disk mode
is explicitly unsupported in this increment because cage post-processing changes
phi fields; legacy memory-backed cage behavior remains tested.

Scientific hashing uses bounded sequential readback rather than a new incremental
hash algorithm. Final JSONL export is still separate work and is measured below;
this increment does not eliminate its cost.

## Completed reconstructed teapot comparison

All four invocations completed with independent UUIDs and complete archives.
They used the same G6/O6 source geometry, pose, objective, tolerances and 960-step
per-pass cap. This is the reconstructed teapot storage regression, **not an
author-original or paper level-9 experiment**.

| Invocation | Full-process peak MiB | Native solve s | JSONL export s | Run total s |
|---|---:|---:|---:|---:|
| Frozen binary, memory | 869.6 | 68.9 | 40.0 | 125.1 |
| New binary, memory | 870.8 | 68.6 | 40.5 | 125.2 |
| New binary, disk 1 | 348.6 | 76.0 | 42.6 | 135.1 |
| New binary, disk 2 | 349.8 | 76.5 | 42.7 | 135.6 |

Disk mode reduced observed full-process peak RSS by **59.97% and 59.83%** relative
to the new memory mode. These are actual process measurements, not estimates
from the number of stored records. Solve-only peak RSS was 847.0 MiB in new
memory mode and 325.2/326.5 MiB in disk mode.

These timings are observations, not a controlled performance benchmark: read-only
verification of previous runs overlapped disk solves, and unrelated host jobs
continued. The disk runs did not demonstrate a speedup. Native spooling and
bounded hash readback add work; the separate JSONL export still costs about 43 s.

All three independently checked comparison kinds passed:

- `legacy-binary`: frozen memory versus new memory.
- `cross-mode`: new memory versus new disk.
- `deterministic-repeat`: two independent new disk invocations.

Each full archive has **2,073 records**: one header, ten pass records and **2,062
iteration records**. Every scientific record is equal, including signed-zero
representation. The ten native binary trace files are also byte-identical between
disk repeats. Each disk run wrote 614,589,497 native trace bytes, retained at most
one pending history record during solving and zero after sealing. Pass partitions
and beta vectors still retain 6,406,400 payload bytes; this is explicitly not all
remaining solver memory.

Shared final hashes:

- JSONL archive: `1c187791ab63f0cb6ebd4f6df1f63dd76ffae16d6914e1f7f1c0aebc7c4bd341`.
- Native scientific trace: `b0431da2d53518c2406bd994ebb0b724c06807e647fb51b0ff76b7bf29fe1f1b`.
- Native scientific output: `b8c52b05f03501f45b0da09b8b6d7a1827e89db9440c99ccd17c8f0f848a6e48`.
- New module: `7f84620395ec7a6ce1c2cba348a37a58ad71c774e24a06e3d823c8c51efe2146`.

The numerical status remains **BINARY_VALIDATION_FAILED in all four runs**.
The G6 normalized binary residual is `5.962606538739975e-6`, above the unchanged
`1e-6` gate. Independent moment-integration disagreement is about `3.95e-16`.
Storage compatibility and repeatability therefore passed; physical stability,
geometry admission and paper reproduction did not gain a pass from these runs.

All source/build/input/output hashes were rechecked after the comparisons by
`controller/collect_final.py`. Its generated `controller/final-observations.json`
binds the three comparison reports, the four invocations, test names/IDs, resource
measurements, unchanged frozen binary and all ten repeated native file hashes.
Each native invocation was single-threaded, limited to 24 GiB address space and
300 seconds; native experiments ran sequentially. Available storage remained above
1.3 TB. No old experiment files were removed.

## Remaining work

1. Spool adaptive pass metadata and integrate cage trace handling without changing
   the existing post-solve phi semantics or scientific evidence.
2. Remove redundant final export/readback costs through a separately verified
   archive/checker path, before relying on this to finish a two-run G7 batch budget.
3. Address and independently validate binary rounding residuals. Do not relax the
   `1e-6` threshold or reinterpret this numerical failure as infeasibility.
4. Original paper assets/parameters, full level-9 cases, geometry/physical gates
   and full release verification remain separate prerequisites.
