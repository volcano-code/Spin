# Reproduction environment lock

Last audited: 2026-09-01 (Asia/Shanghai)

## Reference host toolchain

| Component | Locked/verified value |
|---|---|
| OS | Ubuntu 24.04.3 LTS, x86-64 |
| Kernel during R0 audit | Linux 6.17.0-40-generic |
| CMake | `/usr/bin/cmake` 3.28.3 |
| C/C++ compiler | GCC/G++ 13.3.0 |
| Generator | Ninja 1.13.0 |
| C++ language mode | C++17, extensions disabled |
| Build type | `RelWithDebInfo` |
| Paper/Python CPython | `/home/lin/anaconda3/bin/python3` 3.13.5 |
| pybind11 headers | 3.0.1 |
| NumPy runtime | 2.1.3 |
| Optional simulation | PyBullet build 2026-03-16; trimesh 4.9.0 |
| Reachability LP | HiGHS 1.14.0, commit `7df0786de3088c832297e5ed821db236d8fab281` |
| Spin-It Faster Boolean | CGAL commit `f7c3c8212b56c0d6dae63787efc99093f4383415` |
| Exact arithmetic | GMP 6.3.0; MPFR 4.2.1 |

The project rejects CMake 4.x. The pinned libigl revision invokes legacy
download subprojects whose old compatibility mode was removed in CMake 4.
Using the system CMake above avoids modifying the libigl submodule.

## Dependency lock

The parent repository pins libigl through its gitlink:

```text
libigl f06fb13014aadff7eb1d7423ac018620441a3657
```

On the audited host, that libigl revision resolves the following cached
dependency source commits:

| Dependency | Commit/tag |
|---|---|
| Eigen | `21ae2afd4edaa1b69782c67a54182d34efe43f9c` (`3.3.7`) |
| glad | `09b4969c56779f7ddf8e6176ec1873184aec890f` |
| GLFW | `3327050ca66ad34426a82c217c2d60ced61526b7` |
| ImGui | `5503c0a12e0c929e84b3f61b2cb4bb9177ea3da1` (`v1.76`) |
| ImGuizmo | `a23567269f6617342bcc112394bdad937b54b2d7` |
| libigl-imgui | `99f0643089b19f6daf5b3efd9544a65c9a851966` |
| CGAL (Faster only) | `f7c3c8212b56c0d6dae63787efc99093f4383415` |

These external trees are selected by the pinned libigl CMake scripts and are
not committed to the parent repository. The first clean configure can require
network access; hashes above provide an audit check after resolution.

## Reproduction commands

```bash
git submodule update --init --recursive
/usr/bin/cmake --preset headless
/usr/bin/cmake --build --preset headless --parallel
/usr/bin/ctest --preset headless
/usr/bin/cmake --preset viewer
/usr/bin/cmake --build --preset viewer --parallel
/usr/bin/cmake --preset python
/usr/bin/cmake --build --preset python --parallel
/usr/bin/ctest --preset python
/usr/bin/cmake --preset paper-python
/usr/bin/cmake --build --preset paper-python --parallel
/usr/bin/ctest --preset paper-python
/usr/bin/cmake --preset physics-python
/usr/bin/cmake --build --preset physics-python --parallel
/usr/bin/ctest --preset physics-python
```

Run the four CTest presets sequentially.  Some native regression executables
write atomic-diagnostics fixtures at fixed paths in the system temporary
directory.  Concurrent preset invocations use separate build trees but can
still race on those fixtures and produce a false corruption failure.  The
release gate deliberately invokes the presets one at a time.

Configure separate build trees sequentially as well.  The pinned libigl CMake
scripts use a shared `libigl/external` source/download directory; concurrent
first-time configurations can race while materializing Eigen even when their
binary directories differ.  This is an upstream dependency-layout constraint,
not a compiler or reachability failure.

## Carving reachability configuration

Reachability is controlled by `MAKE_IT_STAND_BUILD_REACHABILITY` and defaults
to `ON`.  The enabled configuration pins and statically links HiGHS; it sets
the deterministic solver controls to one thread, parallel mode off, random
seed zero, and primal/dual feasibility tolerances `1e-9`.  The disabled
configuration neither builds nor links HiGHS, and omits the reachability API
and its eight native tests while preserving all original targets.

The 2026-09-02 milestone gate built and tested both configurations:

| Configuration | Build | CTest |
|---|---|---:|
| `build/headless`, reachability enabled | passed | 129/129 |
| `build/reachability-off` | passed | 121/121 |

The extension-backed Python discovery run passed 140 tests.  Two additional
official N=80 optimizer smoke tests were intentionally skipped because they
require explicit `MAKE_IT_STAND_RUN_OFFICIAL_BETA_*` opt-in and are not part of
the enabled reachability gate.  The full three-case reachability experiment
was instead run twice at its declared resolutions; its two summary files are
byte-identical.

## Spin-It Faster geometry configuration

`MAKE_IT_STAND_BUILD_SPIN_IT_FASTER` defaults to `OFF`.  In that state CMake
forces `LIBIGL_WITH_CGAL=OFF`; the core has no CGAL/GMP/MPFR link dependency
and an invalid Faster-only CGAL override is deliberately ignored.  When the
option is enabled, configuration checks the CGAL revision written in pinned
libigl's download declaration, verifies the resolved checkout's Git HEAD, and
requires the audited system GMP/MPFR versions above.  The core then links
`igl::cgal` for exact mesh Booleans.  No alternate Boolean or voxel fallback
is accepted by the Faster reproduction path.

The presets intentionally use separate build trees. Headless testing does not
enable OpenGL, GLFW or ImGui; the viewer build does not register tests. The
third tree is also headless and enables only the optional CPython extension
and R9/R10 tests. The checked-in Python preset contains no host-specific path and
normally resolves Python and pybind11 from the active environment. On the
audited host, the preset selected Anaconda CPython 3.13 and pybind11 3.0.1.
Because that
Anaconda runtime carries a `libstdc++` older than the GCC 13 extension
requires, CMake asks the selected C++ compiler for its matching
`libstdc++.so.6` and applies it only to extension-backed CTest processes.
The pybind11 package is header-only and compiles against the selected CPython
development module. The extension is not
stable-ABI (`abi3`) and must be rebuilt for each CPython minor version and
platform. Python support is off in the other presets.

PEP 517 wheel builds obtain scikit-build-core and pybind11 from the exact
`pyproject.toml` build requirements, while NumPy is declared as a runtime
dependency. A source distribution includes the tracked pinned libigl source
tree needed by CMake but excludes all locally downloaded `libigl/external`
trees and nested Git metadata; its first wheel build therefore resolves Eigen
through the pinned libigl download machinery. Build source distributions only
from an initialized recursive checkout and inspect their member list and size
before release.

The R9 acceptance audit built both a direct wheel and a source distribution,
built a second wheel from that source distribution, and installed the wheel in
an isolated environment outside the checkout with `PYTHONPATH` unset. Package
and native-extension import, all five installed console commands, a real native
optimization, deterministic diagnostics, and guarded mesh exports passed. The
source archive audit found no downloaded external tree, cache, nested Git
metadata, or host-local absolute path. See `docs/r9_acceptance.md` for the
scope boundary: these checks verify the experiment interface, not R10 paper
experiments.

The active local R10 driver defaults to one discarded warm-up and five
measured repeats per case, fixes OMP/BLAS thread-count variables to one, and
records uncontrolled CPU affinity/governor policies in each manifest. The
driver also binds the exact loaded native-extension bytes into every case;
declared commit/compiler/flags remain caller-supplied build provenance rather
than cryptographic source-to-binary attestation. The
checked-in presets remain `RelWithDebInfo`; a publication timing study must use
an explicitly configured clean Release tree and report its actual flags rather
than merely passing a `Release` label to the driver.

## Clean-state checks

Run after building:

```bash
git status --porcelain
git -C libigl status --porcelain
git submodule status --recursive
git diff --check
```

A reproducible committed R0 checkout must leave both repositories clean. The
submodule status line must begin with a space, not `+`, `-`, or `U`.
