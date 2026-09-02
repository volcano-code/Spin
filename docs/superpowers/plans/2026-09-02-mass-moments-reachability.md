# Mass Moments and Carving Reachability Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Preserve the current Make-It-Stand work, add a canonical ten-component raw-moment API, and classify carving-only static-balance requests with a verifiable continuous-relaxation feasibility result.

**Architecture:** `RawMassMoments10` is a small value type layered over the existing `MassMoments` integrators. `CarvingReachabilityAnalyzer` converts the safe support polygon and per-cell mass contributions into a linear feasibility system, solves it with pinned HiGHS, and validates either the primal witness or an explicitly solved Farkas certificate before assigning a physical status.

**Tech Stack:** C++17, Eigen, CMake 3.28.3, libigl, HiGHS v1.14.0 at commit `7df0786de3088c832297e5ed821db236d8fab281`, CTest, Python 3.10 experiment tooling.

**Spec:** `docs/superpowers/specs/2026-09-02-spin-it-reproduction-design.md`

## Global Constraints

- Keep implementation and results in the Make-It-Balance Git repository.
- Preserve all pre-existing tracked and untracked work before feature edits.
- Do not commit generated multi-megabyte OBJ results; retain them in place and bind them through SHA-256 manifests.
- Use local Git identity `volcano <xuxinyao307@gmail.com>` only.
- `PROVEN_UNREACHABLE` requires an independently checked Farkas certificate.
- A feasible continuous relaxation is `RELAXATION_REACHABLE`, never a claim that a discrete manufacturable cavity exists.
- Numerical failure is `INDETERMINATE`, never physical infeasibility.
- Use test-driven development and make a focused commit after every task.
- Do not begin Spin-It Faster until this plan's complete milestone gate passes.

---

## File Map

**Create**

- `scripts/create_pre_spin_checkpoint.py`: deterministic current-worktree manifest and tracked binary patch generator.
- `tests/python/test_pre_spin_checkpoint.py`: provenance-script regression tests.
- `include/physics/RawMassMoments10.h`: canonical ten-component moment value type.
- `src/physics/RawMassMoments10.cpp`: conversions, transformations, validation, and derived physical properties.
- `tests/raw_mass_moments_regression.cpp`: analytic and transformation tests.
- `cmake/HighsDependency.cmake`: exact HiGHS dependency resolution.
- `include/reachability/CarvingReachability.h`: reachability inputs, statuses, certificates, and analyzer interface.
- `src/reachability/CarvingReachability.cpp`: polygon conversion, LP construction, HiGHS solve, witness and certificate validation.
- `tests/carving_reachability_regression.cpp`: analytic, exhaustive, and certificate tests.
- `python/make_it_stand/reachability.py`: strict JSON representation and experiment-summary helpers.
- `tests/python/test_reachability.py`: serialization and deterministic-summary tests.
- `experiments/spin_it_reproduction/reachability/README.md`: protocol and claim limits.

**Modify**

- `.gitignore`: ignore heavy generated result payloads while retaining reports, manifests, and configurations.
- `CMakeLists.txt`: register moment and reachability sources, pinned HiGHS, and regression tests.
- `include/physics/MassProperties.h`: explicit conversion bridge declarations only.
- `src/physics/MassProperties.cpp`: conversion bridge definitions only.
- `include/api/BalanceModel.h`: public reachability request/result entry point.
- `src/api/BalanceModel.cpp`: build fixed/variable per-cell moment inputs and serialize reachability diagnostics.
- `bindings/pybind_module.cpp`: expose statuses, witness summaries, and `analyze_carving_reachability`.
- `python/make_it_stand/__init__.py`: export the Python reachability API.
- `experiments/paper_sim_v1/REPORT.md`: link the new classifications without altering prior experiment outcomes.

---

### Task 1: Preserve the Current Worktree and Create the Spin Branch

**Files:**
- Create: `scripts/create_pre_spin_checkpoint.py`
- Create: `tests/python/test_pre_spin_checkpoint.py`
- Modify: `.gitignore`
- Create: `experiments/spin_it_reproduction/provenance/pre-spin/manifest.json` by running the script
- Create: `experiments/spin_it_reproduction/provenance/pre-spin/tracked.patch` by running the script

**Interfaces:**
- Consumes: repository root and output directory paths.
- Produces: `create_checkpoint(repository: Path, output: Path) -> dict[str, object]` and canonical manifest schema `spin-it-pre-checkpoint-v1`.

- [ ] **Step 1: Write the failing provenance tests**

```python
def test_checkpoint_records_tracked_diff_and_untracked_hashes(self) -> None:
    manifest = checkpoint.create_checkpoint(self.repository, self.output)
    self.assertEqual(manifest["schema"], "spin-it-pre-checkpoint-v1")
    self.assertEqual(manifest["git"]["branch"], "main")
    self.assertRegex(manifest["tracked_patch_sha256"], r"^[0-9a-f]{64}$")
    self.assertEqual(
        manifest["untracked_files"],
        sorted(manifest["untracked_files"], key=lambda item: item["path"]),
    )

def test_checkpoint_rejects_output_outside_repository(self) -> None:
    with self.assertRaises(ValueError):
        checkpoint.create_checkpoint(self.repository, self.outside)
```

- [ ] **Step 2: Run the tests and verify the missing module failure**

Run:

```bash
PYTHONPATH=python python3 -m unittest discover \
  -s tests/python -p 'test_pre_spin_checkpoint.py' -v
```

Expected: import failure for `create_pre_spin_checkpoint`.

- [ ] **Step 3: Implement deterministic checkpoint generation**

The script must:

```python
def create_checkpoint(repository: Path, output: Path) -> dict[str, object]:
    root = repository.resolve(strict=True)
    destination = output.resolve()
    destination.relative_to(root)
    patch = subprocess.run(
        ["git", "diff", "--binary", "--no-ext-diff"],
        cwd=root,
        check=True,
        stdout=subprocess.PIPE,
    ).stdout
    # Exclude destination itself, sort bytewise relative paths, stream SHA-256,
    # reject symlinks that resolve outside root, and atomically write outputs.
```

The manifest records branch, `HEAD`, local Git identity, UTC timestamp, tracked patch hash, and `{path,size,sha256}` for every untracked regular file. `.gitignore` must ignore raw files below `experiments/**/results*/` while negating `REPORT.md`, JSON manifests, CSV summaries, and SVG/PNG previews.

- [ ] **Step 4: Verify the script tests pass**

Run:

```bash
PYTHONPATH=python python3 -m unittest discover \
  -s tests/python -p 'test_pre_spin_checkpoint.py' -v
```

Expected: all tests pass.

- [ ] **Step 5: Generate and verify the real checkpoint**

Run:

```bash
python3 scripts/create_pre_spin_checkpoint.py \
  --repository /home/lin/xvxinyao/make-it-balance-paper-sim \
  --output experiments/spin_it_reproduction/provenance/pre-spin
sha256sum -c experiments/spin_it_reproduction/provenance/pre-spin/SHA256SUMS
```

Expected: every recorded file reports `OK`; the 1.6 GiB paper-simulation result directory remains in place.

- [ ] **Step 6: Commit the existing source checkpoint without raw generated meshes**

Review `git status --short`, then stage tracked source changes, new source/tests/docs/configurations, lightweight reports, the manifest, and `tracked.patch`. Do not stage ignored raw results.

```bash
git commit -m "checkpoint: preserve paper simulation state before Spin-It"
git switch -c spin-it-reproduction
```

Expected: branch `spin-it-reproduction`; no non-ignored pre-existing changes remain.

---

### Task 2: Add Canonical RawMassMoments10 Storage and Conversions

**Files:**
- Create: `include/physics/RawMassMoments10.h`
- Create: `src/physics/RawMassMoments10.cpp`
- Test: `tests/raw_mass_moments_regression.cpp`
- Modify: `include/physics/MassProperties.h`
- Modify: `src/physics/MassProperties.cpp`
- Modify: `CMakeLists.txt`

**Interfaces:**
- Consumes: existing `MassMoments { mass, firstMoment, secondMoment }`.
- Produces: `physics::RawMassMoments10`, `fromMassMoments`, `toMassMoments`, `mass`, `centerOfMass`, `inertiaAboutOrigin`, and `inertiaAboutCenter`.

- [ ] **Step 1: Write failing ordering and conversion tests**

```cpp
const MassMoments source = makeKnownMoments();
const physics::RawMassMoments10 raw =
    physics::RawMassMoments10::fromMassMoments(source);
require(raw.values() == Eigen::VectorXd::Map(expected.data(), 10));
require(raw.toMassMoments().secondMoment.isApprox(source.secondMoment, 1e-14));
require(raw.centerOfMass().isApprox(source.firstMoment / source.mass, 1e-14));
```

- [ ] **Step 2: Register and run the focused test**

Run:

```bash
cmake --preset headless
cmake --build --preset headless --target raw-mass-moments-regression -j2
ctest --preset headless -R '^raw-moments-(ordering|conversion)$' --output-on-failure
```

Expected: configure or compile failure because the type does not exist.

- [ ] **Step 3: Implement the ten-component value type**

```cpp
namespace physics {
class RawMassMoments10 {
public:
    enum Index { M = 0, LX, LY, LZ, QXX, QYY, QZZ, QYZ, QXZ, QXY, Count };
    static RawMassMoments10 fromMassMoments(const MassMoments& moments);
    MassMoments toMassMoments() const;
    const Eigen::Matrix<double, 10, 1>& values() const noexcept;
    double mass() const noexcept;
    Eigen::Vector3d centerOfMass(double minimumPositiveMass = 0.0) const;
    Eigen::Matrix3d inertiaAboutOrigin() const;
    Eigen::Matrix3d inertiaAboutCenter(double minimumPositiveMass = 0.0) const;
    bool isFinite() const noexcept;
private:
    Eigen::Matrix<double, 10, 1> values_ =
        Eigen::Matrix<double, 10, 1>::Zero();
};
}
```

The stored values are density-weighted raw moments. Density `1` therefore yields geometric volume moments, while multi-material sums remain direct.

- [ ] **Step 4: Run focused and existing mass tests**

Run:

```bash
cmake --build --preset headless --target raw-mass-moments-regression mass-properties-regression -j2
ctest --preset headless -R '^(raw-moments|mass-)' --output-on-failure
```

Expected: all selected tests pass.

- [ ] **Step 5: Commit**

```bash
git add CMakeLists.txt include/physics src/physics tests/raw_mass_moments_regression.cpp
git commit -m "feat: add canonical ten-component mass moments"
```

---

### Task 3: Add Raw-Moment Algebra and Affine Transformations

**Files:**
- Modify: `include/physics/RawMassMoments10.h`
- Modify: `src/physics/RawMassMoments10.cpp`
- Modify: `tests/raw_mass_moments_regression.cpp`

**Interfaces:**
- Consumes: `physics::RawMassMoments10` from Task 2.
- Produces: `operator+`, `operator-`, `scaledDensity(double)`, and `transformed(const Eigen::Matrix3d&, const Eigen::Vector3d&)`.

- [ ] **Step 1: Write failing analytic transformation tests**

```cpp
const Eigen::Matrix3d J = 2.0 * rotation;
const Eigen::Vector3d t(3.0, -4.0, 5.0);
const auto transformed = raw.transformed(J, t);
require(approximatelyEqual(transformed.mass(), 8.0 * raw.mass()));
require(transformed.centerOfMass().isApprox(J * raw.centerOfMass() + t, 1e-12));
require((raw + raw.scaledDensity(2.0) - raw).values().isApprox(
    raw.scaledDensity(2.0).values(), 1e-14));
```

- [ ] **Step 2: Run and verify failure**

Run:

```bash
cmake --build --preset headless --target raw-mass-moments-regression -j2
ctest --preset headless -R '^raw-moments-transform$' --output-on-failure
```

Expected: compile failure for missing algebra methods.

- [ ] **Step 3: Implement exact affine formulas**

For linear transform `J`, translation `t`, and `d = abs(det(J))`:

```text
m' = d m
l' = d (J l + m t)
Q' = d (J Q J^T + (J l)t^T + t(J l)^T + m tt^T)
```

Reject singular, non-finite transforms and non-finite density factors.

- [ ] **Step 4: Run all raw-moment and mass-property tests**

Run:

```bash
cmake --build --preset headless --target raw-mass-moments-regression mass-properties-regression -j2
ctest --preset headless -R '^(raw-moments|mass-)' --output-on-failure
```

Expected: all selected tests pass.

- [ ] **Step 5: Commit**

```bash
git add include/physics/RawMassMoments10.h src/physics/RawMassMoments10.cpp tests/raw_mass_moments_regression.cpp
git commit -m "feat: add raw moment transformations"
```

---

### Task 4: Pin HiGHS Behind a Reachability Build Option

**Files:**
- Create: `cmake/HighsDependency.cmake`
- Modify: `CMakeLists.txt`
- Test: `tests/cmake/assert_reachability_dependency.cmake`

**Interfaces:**
- Consumes: CMake option `MAKE_IT_STAND_BUILD_REACHABILITY`.
- Produces: exact target `highs::highs` from HiGHS tag and commit v1.14.0 / `7df0786de3088c832297e5ed821db236d8fab281`.

- [ ] **Step 1: Write the failing CMake contract test**

```cmake
execute_process(
  COMMAND "${CMAKE_COMMAND}" -S "${SOURCE_DIR}" -B "${BINARY_DIR}"
          -DMAKE_IT_STAND_BUILD_REACHABILITY=ON
          -DMAKE_IT_STAND_HIGHS_GIT_TAG=wrong
  RESULT_VARIABLE status)
if(status EQUAL 0)
  message(FATAL_ERROR "invalid HiGHS revision was accepted")
endif()
```

- [ ] **Step 2: Run and verify failure**

Run:

```bash
cmake -P tests/cmake/assert_reachability_dependency.cmake
```

Expected: failure because the contract test inputs and dependency module are absent.

- [ ] **Step 3: Implement pinned dependency resolution**

`cmake/HighsDependency.cmake` first accepts an installed package only if its reported version is exactly `1.14.0`; otherwise it uses `FetchContent` with:

```cmake
FetchContent_Declare(highs
  GIT_REPOSITORY https://github.com/ERGO-Code/HiGHS.git
  GIT_TAG 7df0786de3088c832297e5ed821db236d8fab281
  GIT_SHALLOW FALSE)
set(FAST_BUILD ON CACHE BOOL "" FORCE)
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(highs)
```

Restore the parent `BUILD_TESTING` value after dependency creation. With `MAKE_IT_STAND_BUILD_REACHABILITY=OFF`, no HiGHS configure or link step occurs.

- [ ] **Step 4: Verify enabled and disabled configurations**

Run:

```bash
cmake -S . -B build/reachability-off -DMAKE_IT_STAND_BUILD_REACHABILITY=OFF -DMAKE_IT_STAND_BUILD_VIEWER=OFF
cmake -S . -B build/reachability-on -DMAKE_IT_STAND_BUILD_REACHABILITY=ON -DMAKE_IT_STAND_BUILD_VIEWER=OFF
cmake -P tests/cmake/assert_reachability_dependency.cmake
```

Expected: all commands succeed and the enabled cache records the exact revision.

- [ ] **Step 5: Commit**

```bash
git add CMakeLists.txt cmake/HighsDependency.cmake tests/cmake/assert_reachability_dependency.cmake
git commit -m "build: pin HiGHS for reachability analysis"
```

---

### Task 5: Define Reachability Inputs and Safe-Polygon Halfspaces

**Files:**
- Create: `include/reachability/CarvingReachability.h`
- Create: `src/reachability/CarvingReachability.cpp`
- Create: `tests/carving_reachability_regression.cpp`
- Modify: `CMakeLists.txt`

**Interfaces:**
- Consumes: fixed and variable `MassProperties`, a `SupportRegion`, and a CCW safe polygon.
- Produces: `ReachabilityCell`, `CarvingReachabilityInput`, `ReachabilityStatus`, `ReachabilityCertificate`, `CarvingReachabilityResult`, and `safePolygonHalfspaces`.

- [ ] **Step 1: Write failing halfspace and validation tests**

```cpp
const Eigen::MatrixXd square = makeSquare(-1.0, 1.0);
const auto halfspaces = reachability::safePolygonHalfspaces(square);
require((halfspaces.A * Eigen::Vector2d::Zero()).maxCoeff() <=
        halfspaces.b.maxCoeff());
require(contains(halfspaces, Eigen::Vector2d(1.0, 0.0), 1e-12));
require(!contains(halfspaces, Eigen::Vector2d(1.01, 0.0), 1e-12));
requireThrowsInvalidArgument([&] {
    reachability::safePolygonHalfspaces(makeClockwiseSquare());
});
```

- [ ] **Step 2: Build and verify failure**

Run:

```bash
cmake --build build/reachability-on --target carving-reachability-regression -j2
ctest --test-dir build/reachability-on -R '^reachability-input$' --output-on-failure
```

Expected: compile failure because the reachability API is absent.

- [ ] **Step 3: Implement focused public types**

```cpp
namespace reachability {
struct ReachabilityCell {
    int cellIndex = -1;
    double mass = 0.0;
    Eigen::Vector2d firstMoment2D = Eigen::Vector2d::Zero();
};
struct CarvingReachabilityInput {
    double fixedMass = 0.0;
    Eigen::Vector2d fixedFirstMoment2D = Eigen::Vector2d::Zero();
    std::vector<ReachabilityCell> variableCells;
    Eigen::MatrixXd safePolygon2D;
    double feasibilityTolerance = 1e-9;
};
enum class ReachabilityStatus {
    AlreadySafe,
    ProvenUnreachable,
    RelaxationReachable,
    Indeterminate
};
}
```

Convert every CCW edge `e = v_next - v` to normalized outward halfspace `A_j = (e_y, -e_x)` and `b_j = A_j dot v`. Reject fewer than three vertices, clockwise/zero area, duplicate edges, non-finite values, non-positive fixed mass, negative cell mass, duplicate cell identifiers, and empty safe regions.

- [ ] **Step 4: Run input tests**

Run:

```bash
cmake --build build/reachability-on --target carving-reachability-regression -j2
ctest --test-dir build/reachability-on -R '^reachability-input$' --output-on-failure
```

Expected: pass.

- [ ] **Step 5: Commit**

```bash
git add CMakeLists.txt include/reachability src/reachability tests/carving_reachability_regression.cpp
git commit -m "feat: define carving reachability model"
```

---

### Task 6: Solve and Validate Reachable Relaxations

**Files:**
- Modify: `include/reachability/CarvingReachability.h`
- Modify: `src/reachability/CarvingReachability.cpp`
- Modify: `tests/carving_reachability_regression.cpp`

**Interfaces:**
- Consumes: `CarvingReachabilityInput` from Task 5.
- Produces: `CarvingReachabilityAnalyzer::analyze(const CarvingReachabilityInput&) -> CarvingReachabilityResult` with retained fractions and independently recomputed residuals.

- [ ] **Step 1: Write failing already-safe and reachable tests**

```cpp
const auto safe = analyzer.analyze(alreadySafeFixture());
require(safe.status == ReachabilityStatus::AlreadySafe);
require(safe.retainedFractions.size() == safe.inputCellCount);

const auto reachable = analyzer.analyze(oneMovableCellFixture());
require(reachable.status == ReachabilityStatus::RelaxationReachable);
require(reachable.retainedFractions.minCoeff() >= -1e-10);
require(reachable.retainedFractions.maxCoeff() <= 1.0 + 1e-10);
require(reachable.maximumPrimalViolation <= 1e-9);
```

- [ ] **Step 2: Run and verify failure**

Run:

```bash
cmake --build build/reachability-on --target carving-reachability-regression -j2
ctest --test-dir build/reachability-on -R '^reachability-(already-safe|feasible)$' --output-on-failure
```

Expected: compile failure for the missing analyzer.

- [ ] **Step 3: Implement the primal LP**

For each polygon halfspace, construct:

```text
sum_i ((A_j dot q_i) - b_j m_i) alpha_i
    <= b_j m_fixed - A_j dot q_fixed.
```

Use column bounds `[0,1]`, row lower bounds `-infinity`, row upper bounds from the right-hand side, zero objective, presolve enabled, one solver thread, and deterministic parallel settings. Recompute the final mass, first moment, projected center, every row residual, and bound residual without using HiGHS-reported residual fields.

- [ ] **Step 4: Run feasible tests twice**

Run:

```bash
ctest --test-dir build/reachability-on -R '^reachability-(already-safe|feasible)$' --repeat until-fail:2 --output-on-failure
```

Expected: both repetitions pass with identical retained fractions.

- [ ] **Step 5: Commit**

```bash
git add include/reachability/CarvingReachability.h src/reachability/CarvingReachability.cpp tests/carving_reachability_regression.cpp
git commit -m "feat: solve continuous carving reachability"
```

---

### Task 7: Produce an Independently Verifiable Infeasibility Certificate

**Files:**
- Modify: `include/reachability/CarvingReachability.h`
- Modify: `src/reachability/CarvingReachability.cpp`
- Modify: `tests/carving_reachability_regression.cpp`

**Interfaces:**
- Consumes: the primal inequality system `C alpha <= d`, including polygon rows and explicit `alpha <= 1`, `-alpha <= 0` rows.
- Produces: normalized non-negative Farkas vector `y`, certificate residuals, and `ProvenUnreachable` only after independent validation.

- [ ] **Step 1: Write failing infeasible and ambiguous-status tests**

```cpp
const auto result = analyzer.analyze(unreachableFixedShellFixture());
require(result.status == ReachabilityStatus::ProvenUnreachable);
require(result.certificate.weights.minCoeff() >= -1e-10);
require(result.certificate.stationarityInfinityNorm <= 1e-9);
require(std::abs(result.certificate.weights.sum() - 1.0) <= 1e-9);
require(result.certificate.contradiction < -1e-9);

const auto invalid = reachability::validateFarkasCertificate(
    fixture.C, fixture.d, Eigen::VectorXd::Ones(fixture.C.rows()));
require(!invalid.valid);
```

- [ ] **Step 2: Run and verify failure**

Run:

```bash
cmake --build build/reachability-on --target carving-reachability-regression -j2
ctest --test-dir build/reachability-on -R '^reachability-(infeasible|certificate)$' --output-on-failure
```

Expected: failure because no certificate path exists.

- [ ] **Step 3: Solve the explicit normalized Farkas LP**

For `C alpha <= d`, solve:

```text
minimize d^T y
subject to C^T y = 0,
           sum(y) = 1,
           y >= 0.
```

Accept the certificate only when all independently recomputed conditions hold:

```text
min(y) >= -tol
abs(sum(y)-1) <= tol
norm_inf(C^T y) <= tol
d^T y < -tol
```

If HiGHS reports primal infeasibility but this second LP or validation fails, return `Indeterminate` and retain both solver statuses in diagnostics.

- [ ] **Step 4: Exhaustively compare small fixtures**

Enumerate all binary states for two- through eight-cell fixtures. Assert:

- any binary reachable state prevents `ProvenUnreachable`;
- `ProvenUnreachable` implies the sampled continuous grid contains no feasible point;
- continuous-feasible/binary-infeasible fixtures are labeled `RelaxationReachable`.

Run:

```bash
ctest --test-dir build/reachability-on -R '^reachability-(infeasible|certificate|exhaustive)$' --output-on-failure
```

Expected: all selected tests pass.

- [ ] **Step 5: Commit**

```bash
git add include/reachability/CarvingReachability.h src/reachability/CarvingReachability.cpp tests/carving_reachability_regression.cpp
git commit -m "feat: certify unreachable carving requests"
```

---

### Task 8: Integrate Reachability with BalanceModel and Python

**Files:**
- Modify: `include/api/BalanceModel.h`
- Modify: `src/api/BalanceModel.cpp`
- Modify: `bindings/pybind_module.cpp`
- Create: `python/make_it_stand/reachability.py`
- Modify: `python/make_it_stand/__init__.py`
- Create: `tests/python/test_reachability.py`
- Modify: `tests/balance_model_regression.cpp`

**Interfaces:**
- Consumes: the same snapshot, voxel grid, deformation domain, cell moments, and standing target used by `BalanceModel::optimize`.
- Produces: `BalanceModel::analyzeCarvingReachability(std::size_t standingIndex = 0) const`, Python `BalanceModel.analyze_carving_reachability`, and strict schema `make-it-stand-carving-reachability-v1`.

- [ ] **Step 1: Write failing C++ adapter test**

```cpp
api::BalanceModel model = makeReachabilityFixtureModel();
model.addStandingConstraintPolygon(squareSupport(), gravity, radians(5.0));
const auto result = model.analyzeCarvingReachability(0);
require(result.status == reachability::ReachabilityStatus::RelaxationReachable);
require(result.variableCellCount > 0U);
require(result.fixedMass > 0.0);
```

- [ ] **Step 2: Write failing Python schema test**

```python
payload = result.to_dict()
self.assertEqual(payload["schema"], "make-it-stand-carving-reachability-v1")
self.assertIn(payload["status"], {
    "ALREADY_SAFE", "PROVEN_UNREACHABLE",
    "RELAXATION_REACHABLE", "INDETERMINATE",
})
self.assertRegex(payload["input_hash"], r"^[0-9a-f]{64}$")
```

- [ ] **Step 3: Run and verify both API failures**

Run:

```bash
cmake --build build/reachability-on --target balance-model-regression make_it_stand_native -j2
ctest --test-dir build/reachability-on -R '^balance-model-reachability$' --output-on-failure
PYTHONPATH=python python3 -m unittest discover \
  -s tests/python -p 'test_reachability.py' -v
```

Expected: missing method/module failures.

- [ ] **Step 4: Implement the adapter without duplicating physics**

Build the regular/deformed cell field exactly as `optimize` does. For each initially filled active cell:

- add non-carveable cells (`Voxel::isInner == false`) to fixed moments;
- add carveable cells to `ReachabilityCell` in ascending linear-cell order;
- project a cell first moment by computing its exact center of mass and applying `SupportRegion::toPlaneCoordinates`, then multiplying by cell mass;
- use `StandingTarget::safeRegion2D()` from the selected constraint;
- serialize only finite JSON numbers and stable, lowercase status tokens.

- [ ] **Step 5: Run C++ and Python API tests**

Run:

```bash
cmake --build build/reachability-on --target balance-model-regression make_it_stand_native -j2
ctest --test-dir build/reachability-on -R '^(balance-model-reachability|reachability-)' --output-on-failure
PYTHONPATH=python:build/reachability-on/python python3 -m unittest discover \
  -s tests/python -p 'test_reachability.py' -v
```

Expected: all selected tests pass.

- [ ] **Step 6: Commit**

```bash
git add include/api/BalanceModel.h src/api/BalanceModel.cpp bindings/pybind_module.cpp python/make_it_stand tests/balance_model_regression.cpp tests/python/test_reachability.py
git commit -m "feat: expose carving reachability analysis"
```

---

### Task 9: Classify the Existing Failed Orientation Experiments

**Files:**
- Create: `experiments/spin_it_reproduction/reachability/README.md`
- Create: `experiments/spin_it_reproduction/reachability/run_cases.py`
- Create: `experiments/spin_it_reproduction/reachability/cases.json`
- Create: `experiments/spin_it_reproduction/reachability/REPORT.md`
- Modify: `experiments/paper_sim_v1/REPORT.md`
- Test: `tests/python/test_reachability.py`

**Interfaces:**
- Consumes: existing orientation experiment inputs/configurations and `BalanceModel.analyze_carving_reachability`.
- Produces: canonical per-case JSON, aggregate report, and classification for Spheres 45 degrees, Toy head-down, and the successful Dinosaur head-contact case.

- [ ] **Step 1: Write failing deterministic report test**

```python
first = reachability.run_manifest(manifest, output_a)
second = reachability.run_manifest(manifest, output_b)
self.assertEqual(first["case_order"], second["case_order"])
self.assertEqual(first["result_hashes"], second["result_hashes"])
self.assertNotIn("PERFECT_REPRODUCTION", json.dumps(first))
```

- [ ] **Step 2: Run and verify failure**

Run:

```bash
PYTHONPATH=python:build/reachability-on/python python3 -m unittest discover \
  -s tests/python -p 'test_reachability.py' -v
```

Expected: failure for missing manifest runner.

- [ ] **Step 3: Implement the fixed three-case manifest and report**

The report must include original carving outcome, reachability status, primal or certificate residuals, cell counts, source/config/input hashes, and one of:

- `PHYSICALLY_UNREACHABLE_UNDER_CONTINUOUS_CARVING_RELAXATION`;
- `CURRENT_PLANE_CARVER_SOLVER_GAP`;
- `ALREADY_SAFE_OR_EXISTING_SUCCESS`;
- `INDETERMINATE_REQUIRES_NUMERICAL_INVESTIGATION`.

Do not tune carving parameters in this task.

- [ ] **Step 4: Run the three cases twice and compare hashes**

Run:

```bash
PYTHONPATH=python:build/reachability-on/python python3 \
  experiments/spin_it_reproduction/reachability/run_cases.py \
  --manifest experiments/spin_it_reproduction/reachability/cases.json \
  --output experiments/spin_it_reproduction/reachability/results
```

Expected: all three cases have a valid non-empty classification and repeat hashes match.

- [ ] **Step 5: Commit lightweight evidence**

```bash
git add experiments/spin_it_reproduction/reachability tests/python/test_reachability.py experiments/paper_sim_v1/REPORT.md
git commit -m "exp: classify carving reachability cases"
```

Do not add generated heavy OBJ/voxel payloads.

---

### Task 10: Run the Milestone Gate

**Files:**
- Modify: `docs/reproduction_environment.md`
- Create: `experiments/spin_it_reproduction/reachability/release-gate.json`

**Interfaces:**
- Consumes: all deliverables from Tasks 1-9.
- Produces: a machine-readable gate deciding whether the Spin-It Faster plan may begin.

- [ ] **Step 1: Build both supported configurations**

Run:

```bash
cmake --preset headless
cmake --build --preset headless -j2
cmake -S . -B build/reachability-off -DMAKE_IT_STAND_BUILD_REACHABILITY=OFF -DMAKE_IT_STAND_BUILD_VIEWER=OFF
cmake --build build/reachability-off -j2
```

Expected: both builds succeed.

- [ ] **Step 2: Run all registered C++ tests**

Run:

```bash
ctest --preset headless --output-on-failure
ctest --test-dir build/reachability-off --output-on-failure
```

Expected: zero failed and zero `Not Run` tests in both configurations.

- [ ] **Step 3: Run Python tests with required dependencies**

Run:

```bash
PYTHONPATH=python:build/headless/python python3 -m unittest discover -s tests/python -v
```

Expected: all enabled core and reachability tests pass; optional dataset tests are reported under their separate preset rather than skipped unexpectedly.

- [ ] **Step 4: Verify deterministic evidence and repository state**

Run:

```bash
git diff --check
git status --short
sha256sum -c experiments/spin_it_reproduction/provenance/pre-spin/SHA256SUMS
```

Expected: no source diff errors, only documented ignored/generated paths, and all preserved files verify.

- [ ] **Step 5: Write and commit the release gate**

`release-gate.json` must set `spin_it_faster_unblocked` to true only when every build/test/evidence check above passes and no reachability result is `INDETERMINATE`.

```bash
git add docs/reproduction_environment.md experiments/spin_it_reproduction/reachability/release-gate.json
git commit -m "docs: close reachability milestone"
```

Expected: the milestone ends with a clean tracked worktree and an explicit decision to begin or block the Spin-It Faster plan.
