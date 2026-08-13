# R9 Python and pybind11 acceptance record

Last audited: 2026-08-13 (Asia/Shanghai)

## Scope

R9 exposes the accepted C++ geometry and optimization pipeline through one
owning `api::BalanceModel` facade, a pybind11 extension, a strict Python
configuration layer, and five generic experiment utilities. It does not run or
claim the official paper cases, ablations, parameter sweeps, performance
alignment, figures, tables, or physical fabrication evidence; those remain R10.

The public workflow is:

```python
model = BalanceModel("horse.obj")
model.set_gravity(...)
model.set_support_polygon(...)
model.add_handle(...)
result = model.optimize()
result.save_outer_mesh(...)
result.save_inner_mesh(...)
result.save_printable_mesh(...)
result.save_diagnostics(...)
```

`BalanceModel` copies and owns its source state. `BalanceResult` owns an
immutable result snapshot, returns copied arrays, and permits printable export
only through the R8 manufacturing-ready seal. JSON configurations use an exact
field schema and reject duplicate, non-finite, unknown, or mistyped input.
Result consumers require the versioned native core fields used for summaries
while permitting additional forward-compatible diagnostic fields.

## Delivered interface

- `include/api/BalanceModel.h` and `src/api/BalanceModel.cpp`: owning C++ facade;
- `bindings/pybind_module.cpp`: `_make_it_stand` extension;
- `python/make_it_stand/`: validated configuration, API wrapper, deterministic
  atomic I/O, diagnostics summaries, CSV/SVG reporting, and CLI entry points;
- `python/run_balance.py`, `batch_experiment.py`, `evaluate.py`, `visualize.py`,
  and `compare_ablation.py`: source-tree command wrappers;
- `pyproject.toml`: isolated wheel metadata and installed console commands.

The comparison tool summarizes supplied diagnostics groups. Its schema and
rendered output explicitly state that it does not execute or establish an R10
paper ablation.

Configured outputs reject parent traversal, collisions, and symbolic links
observed below the caller-selected output root. This fail-closed experiment
guard is not claimed as a security sandbox against concurrent hostile
filesystem mutation.

## Verification evidence

The completed local audit records:

```text
Fresh headless CTest suite                        101/101 passed
Fresh Python-preset CTest suite                  103/103 passed
R9 focused CTest suite                                3/3 passed
Isolated Debug C++ R9 regression                      1/1 passed
Fresh inherited viewer configure/build                 passed
Direct PEP 517 wheel build                             passed
Source distribution and wheel-from-sdist build         passed
Isolated install/import/five-CLI/native-run gate        passed
```

The focused tests exercise the C++ facade, real pybind11 ownership and input
rejection, deterministic native optimization/diagnostics, guarded exports,
strict configuration parsing, batch behavior, evaluation, visualization, and
generic comparison. The complete suite also preserves all 100 R0--R8 tests.

The Python unit discovery contains 52 tests: 10 execute the real extension and
42 cover strict configuration and experiment-tool behavior. The isolated
installation audit ran outside the source and build trees with `PYTHONPATH`
unset. It verified the package-local extension, `--help` for all five installed
commands, functional execution of all five commands, and a real stable-cube
optimization. Two native runs produced identical diagnostics and geometry
hashes, and outer, inner, manufacturing-sealed printable, and diagnostics
exports all succeeded.

The audited wheel contained only the five package modules, the native
extension, and standard distribution metadata. The audited source distribution
contained no downloaded `libigl/external` tree, cache directory, nested Git
metadata, or host-local absolute path; a wheel built from that source
distribution succeeded.

## Completed acceptance gate

The wheel was installed into an isolated environment outside the source and
build trees, with no source-tree `PYTHONPATH`, and passed:

1. `import make_it_stand` and import of its package-local native extension;
2. `--help` for all five installed console commands;
3. one real in-memory or file-backed `BalanceModel.optimize()` call;
4. deterministic diagnostics plus outer, inner, and sealed printable export;
5. a check that the wheel contains the package and extension but no unintended
   libigl headers, caches, nested Git metadata, or host-local paths.

Accordingly, `PYAPI-01` and R9 are `VERIFIED_REGRESSION`. This status applies
only to the Python/API workflow. It does not advance any R10 experiment row,
claim equality to an official paper trajectory or asset, or establish
`VERIFIED_PAPER`.
