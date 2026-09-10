# Task7v1 — independently sealed spin geometry in the native viewer

Scope: bounded first viewer slice from the approved Task7/M6 plan. Base `5b60db63f71dcd014faaf3397198dd4b76d924a2`, branch `spin-it-strict-computational-v2`, worktree `/data/xvxinyao/make-it-balance-paper-sim`. No helper/reviewer agents, GUI launch, host-setting edits, solver changes, physics/CAD/run-authentication claims, push, tag, or merge. This is not completion of Task7, M6, or the paper reproduction.

The controller corrected this increment's identifier from Task7a to Task7v1 because `task-7a-report.md` already described historical Faster work. That report was never edited; its observed SHA256 remains `49795bf4e416a43cd984826fe18bb9c9272306a00344d43009c3dc912f6fdc91`. Early new RED logs retain their original `task-7a-red-*` names. Historical bundles, failed attempts, and older viewer builds were not rewritten.

## Implementation

- `python/make_it_stand/verification/spin_view.py`: public `load_spin_scene(bundle_directory, *, expected_problem_digest)` calls the existing `load_spin_artifact` for independent reintegration, exact partition/clearance validation, original byte identity, expected-input binding, and implementation-staleness checks. It returns owned model-frame arrays, per-region densities, and facts; no solver rerun or JSON/OBJ transfer stage.
- Separate original and final outer boundaries, actual fixed-shell material regions, final whole material, retained design regions, design-domain regions, and newly carved cavity. Explicitly empty cavity remains two `(0,3)` arrays; no dummy geometry. Different densities enter the existing checker exactly once, and complete world mass/center/inertia are used.
- Full symmetric world central inertia goes through `numpy.linalg.eigh`. Eigenvalues are grouped at relative `1e-10` of their actual inertia scale; no unit-sized absolute floor. Only singleton groups produce axis lines. Multiplicities, the display tolerance, and conventional axis signs are exposed. Isotropic geometry produces no unique axes; a repeated xy pair preserves only the unique z line.
- `src/viewer/SealedSpinView.cpp` plus its focused header: `pybind11::embed` with the same configured Python interpreter and build-local actual extension. Isolated CPython ignores user `PYTHONPATH`/`PYTHONHOME`; source/build paths are explicitly prepended. Canonical actual interpreter, package, adapter, extension, and imported project-module paths are validated. All geometry is converted into owning C++ row-major matrices, and hashes of those actual C++ bytes must equal the admitted mesh identities. This is local file/path observation, not execution/build authentication.
- The same admission and C++ conversion are used by graphical and headless entrypoints. Python objects and interpreter are destroyed before display; scene values own their data. Headless output is computed JSON summary only, with counts/hashes, densities, statuses, world physical facts, target, eigenspaces, and actual loaded module paths. Invalid input exits nonzero before window initialization.
- Existing `make-it-stand` executable gains `--sealed-spin-bundle BUNDLE EXPECTED_PROBLEM_SHA256` and `--check-sealed-spin-bundle BUNDLE EXPECTED_PROBLEM_SHA256`. Optional graphical `--frames N` accepts integers `1..10000` for bounded real-display observation. It redraws without automatic rotation. Missing/unknown strict arguments fail explicitly. If Python/exact support is absent, the strict flags produce an explicit unavailable error.
- Focused native libigl viewports show original/result side by side. Shared oblique camera is applied after libigl's automatic initialization framing and synchronized across both viewports; camera transforms never alter scene coordinates. The controls/status area reserves screen space separately from geometry, including framebuffer scaling. Established gray-blue outer/fixed, gold material/retained, blue design, and red cavity colors remain.
- Independent layer toggles; default final whole material plus cavity with triangle-display cutaway enabled. Cutaway omits triangles crossing the chosen x plane and has **no section caps**. Cavity isolation and optional display-copy translation are labelled; there is no export operation for these display representations. Fixed shell is labelled separately from outer surface.
- Complete material center/unique axes and actual target axis/support points only. No bounding square masquerading as support, certified support hull, fabricated waterline, or solver-convergence display. Iteration replay remains deferred pending run admission.
- Persistent labels derive geometry status and artifact digest from the recomputed seal. Labels retain `PARTIAL_BLOCKED`, `original_evidence_admission: BLOCKED_UNKNOWN`, `physical_validation: NOT_IN_SCOPE`, `build_source_attestation: NOT_VERIFIED`, `solver_target_feasibility: UNKNOWN`, `numerical_solver_trace_convergence: NOT_AUTHENTICATED_BY_GEOMETRY_BUNDLE`, `HOLLOW_ONLY`, original posture, and display-not-simulation. Scene/session status changes no release gate. Legacy `--spin-result` remains available and is explicitly labelled legacy/unverified.
- CMake changes are local to Python embedding, viewer compilation/runtime path, and new exact+Python tests. No solver or dependency pins changed. The headless main build registers Python admission tests and requires no GUI; viewer builds additionally register native same-entrypoint admission tests. The older result-view regression remains.

## TDD and build evidence

Commands below ran from the worktree. Log paths are relative to it. Python commands use `/home/lin/anaconda3/bin/python3` with `PYTHONPATH=python:BUILD/python:tests/python` and the selected GCC runtime preload for the Python test driver. Native subprocess tests deliberately remove `LD_PRELOAD` and set `DISPLAY=''`.

1. Before production implementation, `python3 -m unittest test_sealed_spin_view.SealedSpinViewTests -v`: six expected missing-adapter assertion failures, `build/task-7a-red-python-01.log`. Additional native-availability and two-dimensional-eigenspace RED: two expected failures, `build/task-7a-red-native-01.log`.
2. Adapter GREEN: seven tests, 32.635 s, `build/task-7v1-green-python-01.log`. A later small-scale test initially failed because the reused *test-only* L-prism claim helper hardcodes the unit L; replaced that test claim with independent literal scale-invariant normalized moments. Both `task-7v1-small-inertia-01.log` (fixture failure) and `-02.log` (pass) are preserved. No production workaround was added.
3. Fresh configure command:

   ```sh
   cmake -S . -B build/task-7v1-viewer -DCMAKE_BUILD_TYPE=RelWithDebInfo \
     -DMAKE_IT_STAND_BUILD_VIEWER=ON -DMAKE_IT_STAND_BUILD_PYTHON=ON \
     -DMAKE_IT_STAND_BUILD_EXACT_GEOMETRY_CHECKS=ON \
     -DMAKE_IT_STAND_BUILD_SPIN_IT_2014=ON -DMAKE_IT_STAND_BUILD_SPIN_IT_FASTER=ON \
     -DFETCHCONTENT_SOURCE_DIR_HIGHS=/data/xvxinyao/make-it-balance-paper-sim/build/strict-v2-spin/_deps/highs-src \
     -DPython_EXECUTABLE=/home/lin/anaconda3/bin/python3 \
     -DPython3_EXECUTABLE=/home/lin/anaconda3/bin/python3
   cmake -S . -B build/task-7v1-viewer -DLIBIGL_SKIP_DOWNLOAD=ON
   ```

   Both configure logs are preserved as `build/task-7v1-configure-01.log` and `-02.log`. The first default libigl configure unexpectedly refreshed/cloned the GUI external source directories rather than only using the old cache stamps. The exact declared revisions are recorded below; no binary shim or numerical pin alteration was used. Subsequent configuration explicitly disables libigl downloads. The initial automatic dependency refresh is a setup deviation, not evidence of preservation of unobserved prior external working-tree contents.
4. Native build commands uniformly use `ulimit -v 25165824; timeout 3600 cmake --build BUILD --target make_it_stand_viewer -j3`. Build01 failed because pinned libigl `MeshGL.h` assumes the caller included `<cstdint>`; added that include in the new TU, matching existing main.cpp, without editing libigl. Build02 linked successfully. Logs preserve all failures.
5. Real native runtime RED: initial successful admission under the test driver's preload printed libigl constructor help before the summary, so JSON parsing failed (`task-7v1-native-red-runtime-01.log`). Direct unpreloaded startup also failed with `GLIBCXX_3.4.32` absent from Conda's older runtime (`task-7v1-viewer-ldd-red-01.log`, `task-7v1-native-red-rpath-01.log`). The viewer now sets libigl's existing quiet-constructor definition and puts the configured compiler runtime directory before the Python directory in BUILD_RPATH. Build03 passed. This is a narrow viewer runtime correction, not a global loader/host change.
6. Viewer GREEN `ctest --test-dir build/task-7v1-viewer -R '^sealed-spin-view-' --output-on-failure -V`: 2/2 CTest registrations passed, comprising 9 Python and 2 native tests, 54.85 s (`build/task-7v1-view-green-01.log`). These cover real signed cavities/density partitions, byte and copy ownership, empty cavity, full non-diagonal/small-scale inertia, repeated eigenspaces, actual target data, wrong input binding, unsupported scopes, stale/payload/fake-passed rejection, and native malformed CLI. Additional final native case exercises the unique z axis and actual target-point C++ conversion.

7. Self-review found libigl automatically reframes visible layers *after* `callback_init`. Final shared original framing is therefore applied after `launch_init`; projection also synchronizes between viewports. The control reservation accounts for framebuffer scaling, and the window has a readable minimum size. Clarified the old legacy footer. `build/task-7v1-viewer-build-final-01.log` records the actual final compile/link, exit 0.
8. Final main all-target command `ulimit -v 25165824; timeout 3600 cmake --build build/strict-v2-spin -j3` passed (`build/task-7v1-main-all-target-build-01.log`). Final full registered `timeout 3600 ctest --test-dir build/strict-v2-spin --output-on-failure -V` passed **199/199**, zero failures/NotRun/unexpected skips, **518.63 s** (`build/task-7v1-main-full-ctest-01.log`). This includes all nine new Python scene tests, actual source-material/geometry9 checks, and the older `result-view-preserves-sealed-posture` regression. The preexisting dependency-contract tests regenerated their own disposable contract-build directories as part of that suite; historical artifact bundles and prior viewer builds were not rewritten.
9. Final viewer `ctest --test-dir build/task-7v1-viewer -R '^sealed-spin-view-native$' --output-on-failure -V` passed **1/1** registration, **3/3** native cases, **14.74 s** (`build/task-7v1-native-final-01.log`). The new case converts the unique z axis and actual non-origin support point/nonunit target axis. Tests also verify the bundle bytes are unchanged after successful native admission. Native subprocesses have `DISPLAY=''` and remove `LD_PRELOAD`.
10. Re-ran the actual process-maps observer against the final binary in a new directory using `python3 build/task-7v1-observe-runtime.py /data/xvxinyao/make-it-balance-paper-sim/build/task-7v1-runtime-final-01`. It exited 0, with both loader override variables removed. `build/task-7v1-runtime-final-01.log` and that new directory contain the final viewer hash, extension hash, mapped runtime paths, computed JSON, and empty stderr. Nothing was overwritten in the first observation directory.
11. Unsupported configuration was compiled from the actual final source in a second NEW directory:

    ```sh
    cmake -S . -B build/task-7v1-legacy-viewer -DCMAKE_BUILD_TYPE=RelWithDebInfo \
      -DMAKE_IT_STAND_BUILD_VIEWER=ON -DMAKE_IT_STAND_BUILD_PYTHON=OFF \
      -DMAKE_IT_STAND_BUILD_EXACT_GEOMETRY_CHECKS=OFF \
      -DMAKE_IT_STAND_BUILD_REACHABILITY=OFF -DMAKE_IT_STAND_BUILD_SPIN_IT_2014=OFF \
      -DMAKE_IT_STAND_BUILD_SPIN_IT_FASTER=OFF -DBUILD_TESTING=OFF -DLIBIGL_SKIP_DOWNLOAD=ON
    timeout 3600 cmake --build build/task-7v1-legacy-viewer --target make_it_stand_viewer -j3
    ```

    Configure/build exited 0 under the same address-space bound. Both strict flags, invoked on that binary with `DISPLAY=''` and no preload/library-path overrides, exited **1** with `sealed-spin unavailable: viewer requires Python embedding and exact geometry support`. Logs: `task-7v1-legacy-configure-01.log`, `task-7v1-legacy-build-01.log`, `task-7v1-unavailable-check-01.log`, `task-7v1-unavailable-view-01.log`. The legacy constructor retains its old usage preamble; it does not create a window.
12. A direct final native check against preserved `build/task-5f-demo-run-1/sealed` with its recorded problem digest `2f45a44a9c74fd3d02bc68a7c36e11b8d7fdee0a0fe22a9857236338d6ed3929` exited **1**, specifically `stale implementation observation`, before geometry/window work (`build/task-7v1-historical-stale-01.log`). Its recorded native hash is the unchanged main build's `9308e1ed...`, whereas the separate viewer extension is `6fffa706...`. This negative staleness probe is not positive original-problem admission from a self-read digest. The controller must use an independently held binding and fresh new output for positive acceptance.

## Final hashes and self-review

`build/task-7v1-source-build-hashes-01.log` records the full source/build SHA256 list; `sha256sum -c` passed for every listed entry after the main suite and final viewer build. These are file observations, not build provenance authentication.

| File | SHA256 |
|---|---|
| `CMakeLists.txt` | `b65c27599f2aeffab13dd1ac0be620ca3894b454b5a4762757193fd86cabda52` |
| `main.cpp` | `6231ff6f6b5feba9bc8e1ad6a2a0d09b4b1e16fad33d99513dcb84a6aa1026e8` |
| `include/viewer/SealedSpinView.h` | `d8cc4dd8968f28b121ff1bb7c74612037a69d8424d9d40821ba3cb9e49d5d43f` |
| `src/viewer/SealedSpinView.cpp` | `485327b5e704614763e4e9ea77cbca7db0d3feef7fb9ef300944c87c20cb6c1e` |
| `python/make_it_stand/verification/spin_view.py` | `ba71340284e449a41159f443b969d374cf81d50a48226a67545e60eb02e574d7` |
| `tests/python/test_sealed_spin_view.py` | `1e088aa56d730d4bdc231fc5abe2b44d6cad908edea0ff1d4e952b9d1bfb0193` |
| `build/task-7v1-viewer/CMakeCache.txt` | `7e855a6251525c0bd85a68f924b0028613c4492367c056ac5d80ace1624f1256` |
| `build/task-7v1-viewer/make-it-stand` | `2bb5879543c94ea4c415b3ce74068f0f59a419f2e0eee32e259ea8270bebf218` |
| `build/task-7v1-viewer/python/_make_it_stand.cpython-313-x86_64-linux-gnu.so` | `6fffa706fcf8f409dc5f7a5ba1fc5bef5905512a77acb58208dfcc1194f4f7a7` |
| `build/strict-v2-spin/CMakeCache.txt` | `2c3f5a72c82d0ba5d8b5d0cd8869d13d89ea8a21e0daa405e88242e1778e517c` |
| `build/strict-v2-spin/python/_make_it_stand.cpython-313-x86_64-linux-gnu.so` | `9308e1edcc781814ba699ca1912bcb64f06e74d0951ce70011810104e5fa403f` |
| `build/task-7v1-legacy-viewer/make-it-stand` | `5d93867eba41439d2ebb5784e41fdf8a9e61b4709547f8367cd34a4108baa9a7` |
| `build/task-7v1-legacy-viewer/CMakeCache.txt` | `12c2a60d8adf72a5dca9a457a7dca3dcf4ed970aa3d3ca9fc5dc98aa6bb93319` |

Self-review checked the brief against the final diff, public mandatory digest binding, exact loader reuse, owning matrices and byte hashes, true shell/outer separation, complete density-weighted world inertia, degenerate-axis handling, actual support data, immutable display copies, no solver/physics promotion, all strict CLI branches, and conditional CMake registration. `git diff --check` passed. No helper reviewer was used. Source changes are limited to six implementation/test/build files plus this report; controller brief, acceptance notes, cumulative progress, and historical Task7a report are excluded from staging.

Worker implementation/headless verification is complete. Source/configure/build ownership is released to the controller with the scoped commit. The separate real-display/controller acceptance remains outstanding; no worker rendering-success claim is made.

## Runtime/dependency observations

`build/task-7v1-observe-runtime.py` is a test-only evidence script stored under ignored build output. It generated a fresh analytic sealed fixture in the new `build/task-7v1-runtime-observation` directory, launched only `--check-sealed-spin-bundle` with `DISPLAY=''`, removed both `LD_PRELOAD` and `LD_LIBRARY_PATH`, sampled `/proc/PID/maps`, and retained summary/stderr/mapped paths/hashes. It exited 0 and observed the actual build-local extension during admission. The generated artifact is not an original paper asset or solver-run seal.

Observed artifact `9abe6cc0316e130494affb4e4322dd05643f66d18d28f9fa1d8619aa4ba7bdfb`; expected problem `3b9abff17804998eebc1f0b725f6d9ea9310d39a9811baf8d4365b6a0de8c70a`. Native-extension SHA256 `6fffa706fcf8f409dc5f7a5ba1fc5bef5905512a77acb58208dfcc1194f4f7a7`. Runtime log: `build/task-7v1-runtime-observation-01.log`.

Actual mapped paths included system `/usr/lib/x86_64-linux-gnu/libstdc++.so.6.0.33`, `libgmp.so.10.5.0`, `libmpfr.so.6.2.1`, selected `/home/lin/anaconda3/lib/libpython3.13.so.1.0`, and Conda `/home/lin/anaconda3/lib/libz.so.1.2.13`. The C++ runtime startup failure is resolved for the observed default environment; GMP/MPFR used the selected system implementations. Python also loaded its Conda zlib. Standalone extension `ldd` resolves system zlib, so the live-process observation is more specific than standalone ldd. This mixed zlib observation and the still-emitted CMake safe-RPATH warning remain disclosed; no complete dependency closure/build authentication is asserted. A caller overriding loader environment variables can change runtime resolution.

GUI source observations under `libigl/external`:

| Directory | Resolved commit | Declaration |
|---|---|---|
| glad | `09b4969c56779f7ddf8e6176ec1873184aec890f` | exact libigl pin |
| glfw | `3327050ca66ad34426a82c217c2d60ced61526b7` | exact libigl pin |
| imgui | `5503c0a12e0c929e84b3f61b2cb4bb9177ea3da1` | libigl tag `v1.76` |
| libigl-imgui | `99f0643089b19f6daf5b3efd9544a65c9a851966` | exact libigl pin |
| imguizmo | `a23567269f6617342bcc112394bdad937b54b2d7` | exact libigl pin |
| cgal | `f7c3c8212b56c0d6dae63787efc99093f4383415` | existing exact numerical pin |

HiGHS uses existing `build/strict-v2-spin/_deps/highs-src`, commit `7df0786de3088c832297e5ed821db236d8fab281`. Python 3.13.5 and pybind11 3.0.1 were selected by CMake with Interpreter, Development.Module, and Development.Embed. Existing CMake deprecation and Boost.Bind placeholder warnings remain; the new embedding RPATH warning is separately documented above.

## Run interface and remaining boundary

```sh
build/task-7v1-viewer/make-it-stand --check-sealed-spin-bundle BUNDLE EXPECTED_PROBLEM_SHA256
build/task-7v1-viewer/make-it-stand --sealed-spin-bundle BUNDLE EXPECTED_PROBLEM_SHA256
# Only on an already-authorized actual display, for bounded observation:
build/task-7v1-viewer/make-it-stand --sealed-spin-bundle BUNDLE EXPECTED_PROBLEM_SHA256 --frames 120
```

These are build-tree-bound entrypoints; the configured source/interpreter/extension paths must exist. The external expected problem digest must come from the caller's locked original problem binding. A bundle sealed with another extension/checker observation can be stale and is rejected; historical evidence is not rewritten to make it load. Generate new evidence in a new directory with the viewer build's actual extension when required. Do not use `native-run.json` or `scientific.json` as convergence admission.

The worker did not open a graphical window. Real-display rendering, screenshot inspection, and interaction acceptance belong to the controller after source/build ownership handoff. Automated results establish headless admission and conversion, not rendering success. Cutaway is intentionally triangle filtering without caps; display transforms do not represent new material. Support remains target points only; there is no certified hull construction or new physics. Solver execution trace/convergence, physical stability, paper/original-evidence admission, CAD, hydrostatics, deformation, and full Task7/M6 completion remain outside this increment.
