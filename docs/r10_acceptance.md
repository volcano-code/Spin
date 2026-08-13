# R10 local experiment infrastructure acceptance

Last audited: 2026-08-13 (Asia/Shanghai)

## Accepted scope

R10 accepts the repository's **local experiment infrastructure** at
`VERIFIED_REGRESSION`, not a complete reproduction of the paper's experiments.
The accepted surface comprises:

- independently switchable carving, deformation, and Laplacian terms, with
  the five Figure-14-style local profiles covered by focused regression;
- ordered committed-state events, explicit terminal states, initial/final
  geometry and material metrics, and distinct execution, optimization,
  target-satisfaction, and software-manufacturing indicators;
- multi-standing API and JSON configuration support with per-pose gravity and
  all-of target scheduling;
- successful-return-only native phase timing plus an outer monotonic run timer;
- fail-closed canonical manifests, run records, aggregates, SHA-256 artifact
  binding (including the actually loaded native extension bytes), portable
  paths, paper/local metric separation, and strict reports;
- a nine-case, repository-authored CC0 local matrix covering algorithmic
  profiles, three resolutions, two deterministic handle placements, standing,
  and suspension.

The event ledger binds occupancy, geometry, parameters, target snapshots,
epochs and the validated ablation profile at accepted terminals. Final cavity
fraction is kept distinct from net material-mass change and is reported only
when the corresponding boundary audit is valid.

## Evidence gates

The final source gate is:

- 103/103 headless CTests pass;
- 106/106 Python-preset CTests pass;
- the R10 label contains two C++ tests plus one Python evidence CTest;
- the Python evidence CTest contains the canonical-record, local-driver, and
  report contract suites;
- `git diff --check` passes.

A durable clean Release matrix was generated after the implementation and this
acceptance boundary were committed. Every manifest binds clean implementation
commit `4b4303679bdca9091ca8ddcedd18095d5aec2de6`, Release flags
`-O3 -DNDEBUG`, and the exact loaded native-extension SHA-256. All 45 measured
protocol executions completed; six of nine case aggregates report algorithm
success, and five of nine report sealed software-manufacturing readiness. The
three algorithm-level failures (both deformation-only profiles and resolution
20) and the additional resolution-12 manufacturing-audit failure remain
visible in the report. These are local robustness outcomes, not paper parity.

## Paper-alignment ceiling

The following are deliberately **not accepted**:

- EXP-01/Table 1: no five-case official-asset run, comparable beta/paper
  iteration unit, or completed performance discrepancy report;
- EXP-02/Figure 14: the same switch structure is implemented, but the local
  CC0 fixture is not the paper figure's input and no paper-trend parity is
  claimed;
- EXP-03: no quantitative multi-base matrix, handle-noise distribution, broad
  official-object study, or physical experiment;
- fabrication: `manufacturing_ready` is a sealed pre-slicer software geometry
  audit, not proof of slicing, printing, load bearing, or physical standing.

Accordingly EXP-01, EXP-02, EXP-03, and physical OUTPUT-01 remain `PARTIAL`.
The schema structurally requires `exact_reproduction_claim=false`. Modern
committed-event ordinals are not beta-2013 or paper iteration counts.

## External-oracle and license boundary

The official beta archive remains a read-only external oracle. Its archive
SHA-256 is recorded in `docs/provenance.md`; its research-only code,
configurations, `.vox`/`.bbw` files, binaries, and third-party model assets are
not copied into this MPL-2.0 repository or its evidence tree. The inspected
official meshes also require sanitation before the strict modern source-mesh
facade can accept them. This is an explicit reproduction deviation, not a
reason to substitute the CC0 matrix for official-paper evidence.

## Timing interpretation

The native sidecar covers a normal return from `BalanceModel.optimize()` after
configuration validation through diagnostics serialization, using
`steady_clock`. It does not measure exceptions or timeouts. The outer timer
covers model setup through run-record field finalization, including SHA-bound
artifact registration. The local CLI defaults to one discarded warm-up and five
measured repeats per case and fixes OMP/BLAS thread-count variables to one;
affinity and CPU-governor state remain uncontrolled and are disclosed.

The checked-in presets are `RelWithDebInfo`. Publishable timing requires a
separate clean Release build with actual compiler flags recorded. A CLI
`--build-type` string is provenance supplied by the caller, not automatic
binary introspection. The copied native-extension SHA-256 identifies the
executed bytes independently but is not source-to-binary attestation.

## Verdict

**Accepted:** R10 local experiment infrastructure, clean Release matrix, and
truthfully reported local robustness workflow.

**Not accepted:** complete R10 paper experiment reproduction or physical
fabrication reproduction.
