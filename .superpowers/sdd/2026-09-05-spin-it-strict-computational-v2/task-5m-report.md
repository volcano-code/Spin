# Task5m report: explicit Givens cage functional

Bounded implementation complete; overall scientific status **PARTIAL_BLOCKED**.
Base: `f78eac7aed4a8d47af99f3412092dfda2e4f5e42`.
Worktree: `/data/xvxinyao/make-it-balance-paper-sim`, branch
`spin-it-strict-computational-v2`. Implementation commit:
`ce733ac915f322d939e887f04dea40fc5e4a1f58` (six focused owned files).

## Delivered contracts

- Separate `evaluatePaper(beta, phi, handleParameters)` and
  `evaluatePaperJoint(state)` use `Spin2014GivensProblem`, not legacy-objective
  evaluation. Layout is beta, optional phi, every existing free handle column.
- Small count/layout API; static no-angle/phi-zero rule; all five modes and
  independent original inertia/center weights; existing fixed rotations and
  locked-handle basis preserved. An all-locked request is also tested.
- One shared geometry/moment and regularizer evaluation for legacy and paper
  paths. Raw signed beta and full handle Jacobians exposed; wrapper supplies
  complete joint gradient/Jacobian and zero raw angle column.
- Decomposed energies, total derivatives, actual final moments and all
  deformed meshes returned. Invalid sizes, beta bounds, NaN/Inf, static
  nonzero phi, and nonfinite computed output explicitly rejected.
- Existing uniform outer Laplacian and supplied symmetric fill Laplacian
  preserved. Section6.1 factors retain `gamma_D,L ||L_D v'||^2`,
  `gamma_D,v ||v'-v||^2`, and `.5 gamma_L beta^T L beta`. Documentation marks
  uniform graph discretization IMPLEMENTATION_DEFAULT, not Sorkine's
  rotation-aware operator or recovered author discretization.
- Existing core cage test extended; separate exact-feature-gated source-BBW
  integration test, with actual coupled solve and trilinear original/cell
  weights for all 26 handles (104 free handle columns, 107 joint columns).

## TDD and test evidence

All command paths below are relative to the worktree. Fresh retained logs:
`build/task-5m-logs-9xUnJr`. All directly requested builds used
`build/task-5c-no-faster-20260909` (Faster OFF, exact geometry ON, Python ON),
`ulimit -v 25165824`, `timeout 3600`, and `-j3`. The dependency-contract test's
unexpected nested build exception is documented below.

1. RED: declarations and temporary throwing explicit implementation made the
   new behavioral coverage compile. Core build succeeded, then core test
   exited 1 with `explicit cage Givens functional unimplemented`.
   Logs `red-build.log`, `red-test.log`.
2. RED integration: successful actual coupled source BBW solve reached the
   same intentional unimplemented exception. Build succeeded, test exited 1.
   Logs `red-integration-build.log`, `red-integration-test.log`.
3. First GREEN attempt: source integration passed, core failed on the initial
   legacy comparison oracle. Preserved `green-focused.log` and diagnostic
   build. YoYo diagonal-phi paper value 15.6941 versus legacy 8.42843.
   Root cause: legacy `SpinMomentProblem::centralInertia` centers all first
   moments, while paper Givens uses the centroid equality constraints. The
   controller independently confirmed and approved correcting this oracle.
   Original off-manifold geometry/value/derivative tests were retained;
   a separate posed-geometry request translated onto mode centroid equalities
   now checks diagonal-phi legacy VALUE equality. No tolerance was relaxed.
   Ambient gradients of the different extensions are not equated.
4. Corrected focused GREEN: 2/2 pass. `green-corrected-build.log` and
   `green-corrected-focused.log` retain all convergence measurements.
   Core per-mode maximum scaled objective/constraint errors at step 1e-5:
   YoYo 2.524e-10, TopInertiaOnly 2.291e-10, TopLoweringOnly 2.129e-10,
   TopMassAndLowering 2.569e-10, StaticLowering 6.579e-11.
   Errors converge from 1e-3 to 1e-4 to 1e-5; every column and equality is
   checked at final relative tolerance 1e-6 with explicit near-zero absolute
   allowance (core 2e-8, integration 2e-9). Raw derivatives are also checked.
5. Added bounded edge checks (all locked, finite beta out of bounds, finite
   inputs overflowing computed energy) pass in `final-focused-test.log`.
6. Full 2014 registered scope: all 30 runtime tests passed; final
   `spin-2014-dependency-contract` interrupted after its implicit scratch
   build was discovered. Thus **30 PASS, 1 interrupted**, not 31/31 success.
   `scoped-build.log` records successful builds for all 24 native 2014 targets
   and the Python extension. `scoped-tests.log` records all 31 attempts.
   CTest command returned 143 after deliberate TERM of its validated process
   group (3680723); no build descendants remained afterwards. Under the
   controller's subsequent bounded exception, the existing exact-OFF
   dependency build and audited binary completed separately (both exit 0).
   This is additional dependency evidence, not a 31/31 CTest result.

### Unexpected dependency-contract side effect and exact verification gap

The final registered test invokes
`tests/cmake/assert_spin_it_2014_dependency.cmake`, which internally calls
`file(REMOVE_RECURSE)` on `build/spin-2014-disabled-contract` and
`build/spin-2014-enabled-contract`, reconfigures both, and compiles
`spin-it-2014-audited -j3` in the enabled directory. This was discovered during
inspection of its long-running final test, after it had already regenerated
scratch build directories and begun compilation. The controller was informed
immediately, and the test/build process group was stopped. This indirectly
breached the assigned-build-only restriction; it must not be represented as
an authorized fresh variant verification or hidden behind the 30 passing
runtime tests. The script's removal affected only those exact scratch build
directories; prior build artifacts are not restored (rebuildable from source).
No main `build/strict-v2-spin` or protected sibling tree was touched. No
additional cleanup or deletion was performed.

The controller then inspected both existing generated directories and
authorized only finishing the already-configured
`build/spin-2014-enabled-contract` target `spin-it-2014-audited` and executing
it, with the same memory/time/parallel limits. No destructive script or
configure was rerun. This separate build and execution both exited 0.
Configuration: exact geometry OFF, reachability OFF, Faster OFF, 2014 ON.
The audited binary reports actual KKT/restoration/analytic/Givens tests PASS,
Givens KKT residual `6.3137399822221774e-11`. Logs:
`dependency-existing-build.log`, `dependency-existing-test.log`.
The named CTest contract remains interrupted; isolated, non-destructive
dependency-test scratch management is a mandatory followup before broad
suite/release. The exact-OFF cage functional executable itself was not run.

```bash
ulimit -v 25165824
timeout 3600 cmake --build build/spin-2014-enabled-contract --target spin-it-2014-audited -j3
timeout 3600 build/spin-2014-enabled-contract/spin-it-2014-audited
```

Exact command templates used for RED/focused builds and tests:

```bash
ulimit -v 25165824
timeout 3600 cmake --build build/task-5c-no-faster-20260909 --target spin-it-2014-cage-functional-regression -j3
timeout 3600 build/task-5c-no-faster-20260909/spin-it-2014-cage-functional-regression
timeout 3600 cmake --build build/task-5c-no-faster-20260909 --target spin-it-2014-source-bbw-joint-regression -j3
timeout 3600 build/task-5c-no-faster-20260909/spin-it-2014-source-bbw-joint-regression
timeout 3600 cmake --build build/task-5c-no-faster-20260909 --target spin-it-2014-cage-functional-regression spin-it-2014-source-bbw-joint-regression -j3
timeout 3600 ctest --test-dir build/task-5c-no-faster-20260909 -R 'spin-it-2014-(cage-functional|source-bbw-joint)$' -V
timeout 3600 ctest --test-dir build/task-5c-no-faster-20260909 -R '^spin-it-2014-cage-functional$' -V
```

Full build/test (union of registered 2014 names and 2014 labels, including
the otherwise unnamed source material Python test):

```bash
ulimit -v 25165824
mapfile -t task5m_targets < <(rg '^build spin-it-2014-[^/ ]+: CXX_EXECUTABLE_LINKER' build/task-5c-no-faster-20260909/build.ninja | cut -d' ' -f2 | tr -d ':')
timeout 3600 cmake --build build/task-5c-no-faster-20260909 --target "${task5m_targets[@]}" _make_it_stand -j3
timeout 3600 ctest --test-dir build/task-5c-no-faster-20260909 -R '2014|source-spin-material-python' --output-on-failure
git diff --check
```

## Independent controller evidence

Controller's disjoint `task-5m-controller-*` exporter/checker files sit beside
this report, with logs in `build/task-5m-controller-pBCwYb`. They do not form
part of the implementation commit. The controller reported two actual BBW
exports byte-identical and checked all five modes, all 106/107 columns,
nonidentity fixed rotations, separate nonsymmetric box/two-cell geometry,
independent signed-tetra moments and matrix-Givens complex-step derivatives.
I read `independent-check.log`: max raw moment error 5.55e-17, raw Jacobian
2.26e-17, joint derivative 4.79e-16, values 7.22e-16; three-step finite
differences also converge. No production oracle imports.

## Limitations and preserved scope

This is a dense legacy mesh adapter with unit density; independent arbitrary
positive densities are not implemented. Discrete regularizers remain the
legacy compatibility choice. Still mandatory: native source-adaptive
partition ownership/binding, audited SLQP geometry rejection, final material
integration, geometry seal, physical validation and original author cases.
The mesh integration test is not an ownership or final geometry certificate.
No Faster/upstream edits, Python bridge, SLQP integration, source-adaptive
extraction, geometry seal, pass-flag shortcut or original-case claim.

Aside from the interrupted dependency-contract side effect and controller's
bounded completion exception above, only the assigned build was used.
Exact-OFF core compilation/audited runtime succeeded separately; no Faster-ON
variant or exact-OFF cage regression execution is claimed. Core test source
and registration are independent of exact geometry, and the source-BBW test
is explicitly gated. No full-current-tree or unrelated-test success claim.
Unrelated dirty Quadric files and `quadric_certified_empty_cases.h` were
neither edited nor staged; no main strict-v2 build, protected sibling tree,
merge/push/tag/release, or subagent dispatch.

## Hashes / final status

SHA-256 (paths relative to worktree; archive is the final assigned-build
archive after the header-comment and extra-test build):

```text
cdc732e4467cdf62bfa8d65c9cdc18f7a63153d1334e36a476a144700b224647  include/spin2014/CageSpinDeformation.h
7aa2477cc5301df15722f05adea8a74a74776410c3aff2b3816f7a049a748f4f  src/spin2014/CageSpinDeformation.cpp
e186a34ca587d3d44d11a59ece09aec85541b17ae1fec3e56f8369e1709c9e07  tests/spin_it_2014_cage_functional_regression.cpp
d23bd3ff3167ef4d10e359796f508029f46449fb61a9396f487b8665c3d53bee  tests/spin_it_2014_source_bbw_joint_regression.cpp
869d42f31a6868d8ffee2394cccbbd2fc497479a48fee97b34d92b1ecf76cd1b  CMakeLists.txt
892846c15893ba4c69880ea42c298ce952de956c50f0906f3d1a5678f6fc566f  docs/spin2014_cage_givens_functional.md
85a80297bbb7d82dd176df016f821757f11d3e6b1ff4012a304de9d35f9fdc24  build/task-5c-no-faster-20260909/libmake_it_stand_core.a
699e486c08c731d792b2c92e220fd75682181b1b419998b521cdf52af0bd65f9  build/task-5m-logs-9xUnJr/red-test.log
5d322e72f67a5b89a80b0c79114fbf1b4dd4313f8c5d5b411bfc872df76ed58f  build/task-5m-logs-9xUnJr/red-integration-test.log
cd033bdceef28260f05c7d21c596da8ae86582309c732a29ac025d646d09f746  build/task-5m-logs-9xUnJr/green-corrected-focused.log
54159140bc3358f51520fee4fed56626fba2ddfa997344966b10a8da7b2b118d  build/task-5m-logs-9xUnJr/final-focused-test.log
c6ada08e6d5c1ed552b8d1a90443fbb4747a7652a195734e367e397fcf7e152e  build/task-5m-logs-9xUnJr/scoped-build.log
3ca2b59b7c68035c813ba913328fd257f4cc2c132f7048d4bcea187283fb6a61  build/task-5m-logs-9xUnJr/scoped-tests.log
23523bee62976eacdb2eec8780616c7dad3061bb6525a71a8147831edca39da5  build/task-5m-controller-pBCwYb/independent-check.log
35245677a7333961faf5162553a4d416c2157edca5e714fa2bc73cef349a6fca  build/task-5m-logs-9xUnJr/dependency-existing-build.log
01d25b5c038ccda81ff3d24055fce47769ba5a480d4d6e5d98946c1f5499994b  build/task-5m-logs-9xUnJr/dependency-existing-test.log
```

Source-BBW integration solve digest:
`0a56f039e35afe1e236135e4705a103d7929bed11add5e0c4f189ddafe8c6d8f`.
Its scaled all-column central errors were 1.09979e-6, 1.09993e-8, 1.11381e-10.
`git diff --check` passed. Final status contained only the five pre-existing
unrelated dirty/untracked Quadric paths. Report and logs reside in ignored
task/build directories and are not part of the focused source commit.
