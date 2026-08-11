# Reproduction environment lock

Last audited: 2026-08-11 (Asia/Shanghai)

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
```

The two presets intentionally use separate build trees. Headless testing does
not enable OpenGL, GLFW or ImGui; the viewer build does not register tests.

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
