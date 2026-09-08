# Task4c3a — actual Faster to complete material bridge

Base: `c5953093da80b84541bcbd025f4decd061663112` on `spin-it-strict-computational-v2`.
Implementation owner: this subagent only; no subordinate agents, merges, pushes,
asset downloads, geometry repairs, or historical-output overwrites.

## Scope and public behavior

Implemented the bounded Task4/6 prerequisite, not full M2/M4 or publication
admission. `HOLLOW_ONLY_MESH_MATERIAL_V1` accepts exactly one actual design
region and any supported actual fixed regions with their own positive densities.
The locked problem independently checks the initial fixed/design partition.
The bridge has no fixed-moment, final-moment, geometry-passed, or hydrostatic input.

New APIs:

- `make_it_stand.SpinEvaluatorConfig` / `make_it_stand.spin.SpinEvaluatorConfig`:
  immutable explicit derivative and triangulation controls, matching native
  legacy defaults when omitted. Both Python and native validation reject unknown
  derivative methods and invalid controls; native solve validates before starting.
- `evaluate_spin_it_faster_quadric(..., evaluator_settings=...)`: opt-in actual
  immutable geometry, fixed-normal/source provenance, derivative applicability and
  scientific hashes. Omitted native settings retain the historical six-key dict;
  old Python/native positional calls and legacy JSON/hash behavior remain intact.
- `solve_spin_it_faster(..., evaluator_settings=...)`: real native solver settings
  forwarding. The existing result additionally exposes copy-on-read derivative
  metadata and immutable converged solution coefficients/final target moments.
  The latter are `None` on nonconverged runs.
- `verification.faster_material.evaluate_and_seal_faster_candidate(problem,
  coefficients, *, evaluator_settings, solver_settings)`: native evaluation of
  actual p and -p, one call each, then independent complete material sealing.
  Result is `NOT_OPTIMIZED_CANDIDATE`, target feasibility `UNKNOWN`, project
  `PARTIAL_BLOCKED`, physical validation `NOT_IN_SCOPE`.
- `verification.faster_material.solve_and_seal_faster(problem, *,
  initial_quadrics, evaluator_settings, solver_settings)`: actual native solve,
  explicit initial sequence, immutable event/restart/status identity. Nonconverged
  runs return no artifact. Converged exported material would pass the same seal;
  a seal failure is separately recorded in `seal_error` with no artifact and
  without replacing the native run status. No converged geometry was fabricated.

Strict bridge settings require explicit `DISCRETE_ELLIPSOID` plus
`DUAL_MONOMIAL_RAW_V1`. Rotational modes require no support points. Per controller
ruling, static balance accepts exactly one locked support point exactly equal to
the locked target origin; other static supports are unsupported. Buoyancy,
multiple design regions, deformation, CAD and hand-set hydrostatics are rejected.

## Frame, moments, density, and complete geometry

Let locked normalization be `(c,L,R)`. Native model vertices are
`x_n=(x_model-c)/L`, each material density is `rho_n=rho/R`, and target origin is
`t=(target.origin-c)/L`. Derive target z by normalizing locked `target.axis`;
project the least-aligned Cartesian unit vector onto z's perpendicular plane to
obtain x, then y=z cross x. `A=[x y z]` is right-handed; this is deterministic
coordinate construction, not a pose optimization. Native target coordinates are
`x_t=A^T(x_n-t)`. Coefficients use native matrix-entry p in this normalized target
frame; the strict equation profile supplies the existing internal dual mapping.

Every fixed region is integrated from its actual normalized signed boundary by
the native published surface integrator, with density entering exactly once.
Their summed moments are passed in native model coordinates; native evaluation
itself transforms this fixed contribution into the target frame. Native retained
design moments receive only the design's normalized density.

To convert native target moments back, use
`m_n=m_t`, `l_n=A*l_t+m_t*t`, and
`Q_n=A*Q_t*A^T+(A*l_t)*t^T+t*(A*l_t)^T+m_t*t*t^T`.
The full symmetric Q is reconstructed, including QYZ/QXZ/QXY. The independent
artifact signed-tetrahedron integrator separately verifies every density-labelled
retained/fixed region and all normalized components at absolute error <=1e-6.
No native moments are passed through the independent checker to generate its
expected value.

Native retained/carved coordinates are mapped back via `x_model=L*x_n+c`.
Native outer output bytes are checked against the normalized locked request;
the separately sealed raw outer is the exact original snapshot. For internally
carved geometry the complete material boundary is the locked original signed
boundary plus the reversed actual carved boundary. Existing exact checks verify
that this is precisely fixed material plus retained design, each design splits
exactly into retained/carved, regions do not overlap, and new carve clears the
original exterior by the locked minimum wall. No remeshing, snapping or tolerance
relaxation is performed.

The non-dyadic axis `(1,2,3)` case is explicitly rejected by the exact partition
check after floating-point frame transport. This is a supported-scope limitation,
not repaired geometry or a successful candidate. Exact representable frame
positives include cyclic target axes, nonzero target origin, scale 2, and world
offsets `(2^20,-2^19,2^18)`. Both native complementary evaluations must be
applicable; unsupported/empty/contact constructions can fail explicitly.

## Scientific identity and persistence

Native base request identity remains separately available. The effective native
evaluation/solve request hash includes every nondefault evaluator control;
literal legacy defaults return the historical request hash. The discrete
derivative scientific identity now incorporates this effective request hash.
The immutable bridge identity additionally binds the locked input digest,
normalized request digest, frame, all evaluator/solver settings, input p or
initial sequence, native evaluation identities or actual event/status/count
stream, and successful artifact scientific digest. It contains no clocks.
All returned geometry is tuple-immutable or copy-on-read; result identity is
strict JSON copied on read, and result fields are frozen.

Existing `save_spin_artifact` / `load_spin_artifact(...,
expected_problem_digest=...)` persist and recheck the independent material seal.
That bundle alone does not persist or attest the separate solver/candidate run
binding. The export helper writes separate run JSON for inspection without
introducing a new authenticated bundle format.

## TDD and retained verification attempts

Working directory for every command: `/data/xvxinyao/make-it-balance-paper-sim`.
Python environment used:
`PYTHONPATH=python:build/strict-v2-spin/python:tests/python`
and `LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6`, interpreter
`/home/lin/anaconda3/bin/python3`.

1. Initial behavioral RED before production: `python3 -m unittest
   test_faster_material_bridge -v`. Actual output: `Ran 5 tests in 0.001s`,
   `FAILED (failures=5)`: four `locked Faster material bridge missing` assertion
   failures and one `immutable evaluator settings missing` assertion failure.
   This caught the missing usable operation/settings, not a mocked solver.
2. `task-4c3a-red-native.log`: same command after facade work but before native
   rebuild/bridge implementation. Actual `FAILED (failures=4, errors=1)`; settings
   correctly could not reach the old DSO (`native Faster evaluator settings
   support is unavailable`) and the locked bridge was still absent.
3. Native rebuilds `cmake --build build/strict-v2-spin --target _make_it_stand
   -j 3`: retained `task-4c3a-build-01.log`, `-02.log`, `-03.log`, all exit 0.
4. First GREEN: same five tests, `task-4c3a-focused-01.log`, `Ran 5 tests in
   30.683s`, `OK`. Prescribed geometry/coefficient/tolerance fixture unchanged.
5. Compatibility RED: `python3 -m unittest
   test_faster_material_bridge.FasterSettingsTests.test_public_settings_and_legacy_native_dictionary_compatibility
   -v`, `task-4c3a-red-compat.log`, public config assertion failed before export.
6. Native invalid-control RED: `python3 -m unittest
   test_faster_material_bridge.FasterSettingsTests.test_native_invalid_evaluator_controls_reject_before_solver
   -v`, `task-4c3a-red-invalid-config.log`, `AssertionError: ValueError not raised`.
   Corrected by shared native validation before solver/hash/evaluation.
7. Expanded GREEN: `ulimit -v 25165824; timeout 3600 env ... python3 -m unittest
   test_faster_material_bridge -v`, `task-4c3a-focused-02.log`, `Ran 13 tests in
   52.796s`, `OK`. Includes native invalid enum rejection, tampering, target scope,
   offdiagonal Fraction anchor, native frame transforms and explicit rounding fail.
8. Solver exposure RED: `python3 -m unittest
   test_faster_material_bridge.FasterSettingsTests.test_solver_exposes_immutable_actual_export_and_derivative_metadata
   -v`, `task-4c3a-red-solve-exposure.log`, missing derivative metadata assertion.
9. Actual Faster-disabled contract: `PYTHONPATH=python:build/task4c2c-exact-only-7Hs7nz/python`
   with the same interpreter/preload imports the real exact-only DSO and core
   package, then invokes the bridge. `task-4c3a-no-faster.log` records exit 0 and
   `strict Faster material bridge requires native Faster evaluator support`.
   No native-feature-missing condition is counted as positive geometry evidence.
10. Final all-target build: `ulimit -v 25165824; timeout 3600 cmake --build
    build/strict-v2-spin -j 3`, `task-4c3a-build-final.log`, exit 0. Existing Boost
    global-placeholder deprecation notes are present; no build error.
11. Final focused registered gates: `ulimit -v 25165824; timeout 3600 ctest
    --test-dir build/strict-v2-spin -R
    '^(spin-faster-material-bridge-python|spin-faster-python-api|sealed-spin-artifact-python)$'
    --output-on-failure -j 1`, `task-4c3a-focused-final.log`, 3/3 passed,
    116.64s. Legacy JSON SHA256 assertion
    `c744e1daf10d8349b42801e787a9b54193464facfb62725dcc14b098fcf04dd8` passed.

12. Full registered suite, once after final production source and all-target
    build: `ulimit -v 25165824; timeout 3600 ctest --test-dir build/strict-v2-spin
    --output-on-failure -j 1`, `task-4c3a-ctest-full.log`, exit 0:
    `100% tests passed, 0 tests failed out of 187`, total 253.93s. New bridge
    target passed in 57.03s; actual native discrete differential target passed in
    8.73s. The new test module now has 15 tests, including six fixed slabs with
    mixed densities. No feature-missing skip is used in this module.
13. Post-final-source export/reload: `ulimit -v 25165824; timeout 3600 env ...
    python3 .superpowers/sdd/2026-09-05-spin-it-strict-computational-v2/task-4c3a-export.py
    /data/xvxinyao/make-it-balance-paper-sim/build/task-4c3a-candidate-ae1OAozK`,
    `task-4c3a-export.log`, exit 0. The root was exclusively created with
    `mktemp -d .../build/task-4c3a-candidate-XXXXXXXX`; no existing bundle was
    overwritten. Actual save/reload used the external expected problem digest.
14. `git diff --check` passed. `sha256sum` of all nine production/test/build
    source files, primary DSO, discrete/solver binaries and primary CMake cache
    is retained as `task-4c3a-final-sha256.log`.

## Persisted engineering evidence

Inspection root:
`/data/xvxinyao/make-it-balance-paper-sim/build/task-4c3a-candidate-ae1OAozK`.
Contains a README, `material/` native-format existing sealed bundle (all original,
fixed, design, retained, carved and complete material snapshots), separate
`candidate-run.json`, `unsuccessful-solve.json`, and `summary.json`.
Generated build artifacts remain preserved locally; the export helper, summary
log, test/build logs and hashes are the scoped committed evidence.

- Candidate: `NOT_OPTIMIZED_CANDIDATE`; physical mass 108.42024156935989.
- Actual fixed shell mass is 95.2, analytic Q diagonal 2108/15. Its native
  normalized moment vector is approximately `(7/8,0,0,0,31/384,31/384,31/384,0,0,0)`.
- Maximum independent normalized moment error: `4.440892098500626e-16`.
- Problem digest:
  `33af8c98e4fac8ffe75b9d809e2f35f25da2c9f3c0198e367467d7639f0c8879`.
- Candidate scientific digest:
  `bd8c673fa2504fbd481e5b0394a54ba55b0cd4b7556467c3cf23e1b6b60e2741`.
- Artifact scientific digest:
  `9e73c79d39aadee8e15961d223c42c7eec290f341cb1f6e0a6fd77845ca1dbea`.
- Artifact observed-file evidence digest:
  `48e34bf22fca1cc8c5a64cfd2b609bd68d20ea4c4fb8ddff3c19a3b3cbbaa227`.
- Actual native solve: `RESOURCE_LIMIT_REACHED`, two real evaluations in two
  restarts `[1,1]`, no sealed final artifact, no infeasibility inference.
- Unsuccessful solve scientific digest:
  `b4931473c84828fac9fa7d62f2518b9d8f3a6d103e507c9db11c296aba760a3e`.
- Actual native event hash:
  `fb1d94dc88ec7863700a241beddb14b74fa652eeca27dc96019f14806ba1065f`.
- Final primary DSO SHA256:
  `00c14cc5f79daf9012b2d2c1ae5bbc30dae3db16032d25dc71babd86e499a24e`.
- Primary CMake cache SHA256:
  `2c3f5a72c82d0ba5d8b5d0cd8869d13d89ea8a21e0daa405e88242e1778e517c`.

## Self-review and remaining limits

Reviewed complete source diff and tests after the final build/full suite. Checked
full symmetric Q and translation terms against the controller's independent
Fraction anchor; density appears once per region; fixed native moments are
supplied in model-normalized coordinates before native target transformation;
complete material includes fixed shell and passes exact partition and wall
checks; original raw outer stays byte-identical. Checked explicit settings
reach native solve/evaluation, all controls bind effective request identity,
default JSON fingerprint stays unchanged, and native invalid enum/configs do
not become normal solver failures or silently select legacy behavior.

Result and settings data are frozen, returned geometry uses immutable tuples or
copies, copied JSON mutation cannot modify the scientific digest, and candidate
construction cannot label itself converged/feasible/publication-ready. All
positive geometry evidence uses actual native evaluation plus independent
recomputed sealing. The native binding and existing spin facade are relatively
large existing files; changes are confined to the relevant settings/exposure
sections and no unrelated restructuring was performed.

Bounded task is implemented with documented concerns: there is no naturally
converged strict geometric solve fixture, so the converged export/seal-error
branch has not been demonstrated by a successful optimization. Non-dyadic
frame rounding can fail exact partitions; unsupported geometric derivatives,
contacting boundaries and empty complementary material can fail explicitly.
The existing original-source admission, publication status, physical stability,
manufacturability, CAD and hydrostatics remain outside this material result.
`MESH_MATERIAL_VALIDATED` and connectivity are not physical feasibility claims.
Independent controller review remains pending.

Implementation commit: `4425e20` — `feat(spin): bridge locked Faster runs to
complete material seals`. Postcommit `sha256sum -c task-4c3a-final-sha256.log`
verified all 13 recorded source/DSO/binary/cache hashes; worktree was clean.
Source-only `git diff HEAD^ HEAD --check -- <the nine source paths>` passed.
The all-file staged whitespace check reported four pre-existing CMake/Boost
trailing-space lines in the raw `task-4c3a-build-03.log`; those original log bytes
were deliberately retained as evidence rather than edited. No source whitespace
error or test failure was involved.

## Source files

- `python/make_it_stand/verification/faster_material.py`
- `python/make_it_stand/spin.py`, `python/make_it_stand/__init__.py`
- `bindings/pybind_module.cpp`
- `include/spin/QuadricMomentEvaluator.h`
- `src/spin/QuadricMomentEvaluator.cpp`, `src/spin/QuadricHollowingSolver.cpp`
- `tests/python/test_faster_material_bridge.py`, `CMakeLists.txt`
- This scoped report, retained logs, and `task-4c3a-export.py` evidence helper.
