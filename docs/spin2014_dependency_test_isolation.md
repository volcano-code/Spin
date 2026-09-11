# Spin-It 2014 dependency-test isolation

The `spin-2014-dependency-contract` CTest keeps every invocation in a new,
retained directory below the build tree that invoked it. It no longer removes
or reuses `build/spin-2014-disabled-contract` or
`build/spin-2014-enabled-contract` in the source checkout.

CTest passes its current binary directory explicitly:

```sh
ctest --test-dir build/task-5c-no-faster-20260909 \
  -R '^spin-2014-dependency-contract$' -V
```

For a direct script invocation, the caller must provide an existing absolute
binary root that is not the source root:

```sh
cmake \
  -DMAKE_IT_STAND_DEPENDENCY_TEST_BINARY_ROOT=/absolute/path/to/build \
  -P tests/cmake/assert_spin_it_2014_dependency.cmake
```

Before any configure or build starts, the script prints
`SPIN_IT_2014_DEPENDENCY_RETAINED_DIR=<path>`. A short allocation lock protects
the random, collision-rejecting directory creation and is released before the
expensive child work. A pre-existing symlink at the reserved run namespace is
rejected rather than followed.

Each retained run contains separate `disabled` and `enabled` build trees.
Every configure, build, and audited-runtime phase has a `.command`,
`.stdout.log`, `.stderr.log`, and `.status` file. Phase status is written as
`RUNNING` before launch and replaced by the actual exit result afterward;
stdout and stderr are written directly so an interrupted CTest still leaves
usable partial diagnostics. `contract.status` is `RUNNING`, `FAILED`, or
`PASSED` according to the last durable contract state.

The contract still checks that 2014-OFF excludes its implementation, and that
2014-ON with Faster-OFF and reachability-OFF contains the audited implementation,
has no hidden CGAL linkage, builds the real `spin-it-2014-audited` target, and
runs that executable. The separate `spin-2014-dependency-isolation` test uses
only generated tiny repositories to exercise retention, validation,
concurrency, interruption-visible logs, and failure metadata; fixture success
is not evidence for the real project dependency contract.

Run directories and fixture directories are intentionally retained. This test
does not provide deletion, pruning, reuse, or passed-result shortcuts; storage
management is outside the contract.
