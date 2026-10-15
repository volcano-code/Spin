# Spin-It 2014 report-claim audit (2026-09-20)

## Scope

This audit preserves historical outputs and corrects how future outputs from
`experiments/spin_it_reproduction/spin2014/run_paper_cases.py` describe their
evidence. It does not change or upgrade any historical experiment result.

The preserved `results-level9.json` has SHA-256
`a92e6ec38ec6e7197973a6366d931ca87749de99f5fa1421764ef35d57b49f15`; its
bytes were not edited. The earlier [strict correction](../experiments/spin_it_reproduction/STRICT_CORRECTION.md)
already classifies it as `SYNTHETIC_REGRESSION`, with the synthetic T-Rex case
further limited to `KNOWN_SOLUTION_VALIDATION`.

## Root cause

The old aggregate fields exceeded what their inputs demonstrated:

| Historical field | What the old runner actually checked | Correct interpretation |
|---|---|---|
| `paper_resolution_reproduced=true` | The requested maximum-level setting was at least the paper's level 9 | A configured cap is not proof of a solved leaf reaching level 9; the old result records no observed maximum solved level. The inputs are 16³ reconstructed structural assets, not author assets. |
| `paper_numerics_reproduced=true` | The selected binary-constraint tolerance was no looser than the paper tolerance | Selecting a tolerance is not a numerical comparison against the paper's result. |
| `cage_deformation_executed=true` on every case | The old conditional returned true when cage deformation was not requested | Non-requested work is not executed work; only requested runs with successful cage results count. |
| `all_public_computational_cases_passed=true` | All reconstructed fixture runs met this runner's local checks | This is a synthetic regression outcome, not completion of the papers' public experimental cases. |

The relevant generator creates deterministic `RECONSTRUCTED_ASSET` voxel
inputs with dimensions 16×16×16 and sets their paper target level separately.
The historical level-9 batch overrode the configured cap to 9, but its JSON
does not contain observed solved-leaf depth evidence. The old report therefore
cannot establish either original-asset reproduction or actual level-9
optimization.

## Corrected output semantics

The current runner emits schema `spin-it-2014-reconstructed-run-v2`, records
`scope=RECONSTRUCTED_SYNTHETIC_REGRESSION`, and separates:

- `configured_level_cap_at_least_paper` from `observed_maximum_solved_level`;
- actual `paper_resolution_reproduced`, which requires eligible asset
  provenance and observed solved depth;
- selection of the published binary tolerance from a real numerical
  comparison against published output;
- requested cage runs from successfully executed cage runs;
- `all_reconstructed_regression_cases_passed` from
  `all_public_computational_cases_passed` (the latter remains false for these
  reconstructed inputs).

The existing native result type currently does not expose its deepest solved
leaf level, so the runner reports that evidence as unavailable rather than
inferring it from the configured cap. This limits claims; it does not prevent
the local reconstructed regression from passing.

## Verification performed

- Red/green tests cover configured-cap versus observed-depth semantics,
  rejection of depth evidence above its configured cap, reconstructed-asset
  claim rejection, cage request/execution semantics, and requested-case
  aggregation.
- `tests/python/test_spin_it_2014_paper_runner.py`: 10 passed using the
  source-bound G7 extension and its required `libstdc++` preload.
- A fresh two-repeat level-4 ellipsoid regression emitted
  `scope=RECONSTRUCTED_SYNTHETIC_REGRESSION`, reported local reconstructed
  checks passed, kept paper resolution/numerics false, and marked unrequested
  cage deformation unexecuted. This smoke run is not an author-paper case.
- The immutable level-9 JSON and its recorded SHA-256 were rechecked after the
  change; the historical file remains unchanged.

## G7 repeat remains incomplete

The first source-bound reconstructed G7 run reached `BINARY_MOMENT_VALIDATED`
after 21 passes in 3,544.3 seconds, with its terminal binary-trace audit and
sealed geometry handled in separate later jobs. The from-initial-state
determinism repeat reached only 20 complete passes and hit the 3,600-second
outer wall limit. Both runs used the same native extension SHA-256
`c85dccd584f128e8432add59aa8674cada40d1e5bc92ed7f4385e99aa311deb0` and
configuration SHA-256
`90c5582a71995119bd06d572a52887e45cc9b066b143a316af317b4d577294ec`.
All 20 common completed stage headers have identical trace and output hashes.
The repeat's pass 21 was in flight at timeout and has no completed stage
receipt, so it cannot be counted as a second complete solve. The evidence
supports a deterministic 20-stage prefix, not deterministic terminal
convergence. This audit did not resume or retry that run. Thus the exact
two-run deterministic acceptance gate remains open, and this reconstructed
teapot does not count as an original Spin-It case.
