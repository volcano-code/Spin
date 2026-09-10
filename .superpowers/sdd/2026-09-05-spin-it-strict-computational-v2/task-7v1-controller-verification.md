# Task7v1 controller verification

Implementation under test: `e780f7acbcbe0cd0a2d1c928c72433c37764ce91`; worktree `/data/xvxinyao/make-it-balance-paper-sim`. This closes only the bounded verified-geometry native display increment. Full reproduction remains `PARTIAL_BLOCKED`.

## Tests and builds

Main all-target build and final viewer link completed. Controller inspected actual full CTest log: 199/199 passed, 518.63 s; no Failed, Not Run or unexpected skips. `build/task-7v1-main-full-ctest-01.log` SHA256 `5b179bab7852771da053b33f69f713ff4b25084701b07827f4d8c72144bc7897`. Post-commit controller native admission: three actual cases passed, 14.665 s (one CTest registration, 14.79 s), `build/task-7v1-controller-native-01.log`.

Final viewer SHA256 `2bb5879543c94ea4c415b3ce74068f0f59a419f2e0eee32e259ea8270bebf218`; actual viewer extension `6fffa706fcf8f409dc5f7a5ba1fc5bef5905512a77acb58208dfcc1194f4f7a7`. Eleven recorded source/build hashes were independently rechecked. Unsupported feature configuration rejects strict entrypoints. Old Task5f bundle is rejected by the different viewer extension; its files were not modified.

## Two actual source runs

After source review approval, ran the existing real source runner twice, sequentially, into exclusive new directories `build/task-7v1-source-run-1` and `-2`. No synthetic known-answer constructor was used in place of optimization; the input itself is deliberately classified `SYNTHETIC_REGRESSION`, not an original paper case.

Command for each run (replace N with 1 or 2):

```bash
ulimit -v 25165824
set -o pipefail
timeout 3600 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  PYTHONPATH=python:build/task-7v1-viewer/python \
  LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6 \
  /home/lin/anaconda3/bin/python3 \
  experiments/spin_it_reproduction/strict/run_source_spin2014_material.py \
  --output-directory build/task-7v1-source-run-N \
  --build-directory build/task-7v1-viewer \
  --dependency-root libigl=/data/xvxinyao/make-it-balance-paper-sim/libigl \
  --dependency-root highs=/data/xvxinyao/make-it-balance-paper-sim/build/strict-v2-spin/_deps/highs-src \
  --dependency-root pybind11=/home/lin/anaconda3/lib/python3.13/site-packages/pybind11/include
```

Both exited zero, `BINARY_MOMENT_VALIDATED`, sealed, freshness CURRENT. Runner elapsed observations were 24.58175419602776 s and 10.33681309298845 s; these are not hardware-comparable paper timings. The 24 GiB bound is process address space, not a cgroup/aggregate RSS guarantee. This small regression did not add peak-RSS measurement. Dependency manifests observe declared roots/files and loaded native modules, not authenticated build provenance.

Controller `task-7v1-verify-demos.py` independently reloads both seals using literal original digest `2f45a44a9c74fd3d02bc68a7c36e11b8d7fdee0a0fe22a9857236338d6ed3929`, checks analytical world/normalized ten moments, complete inertia, outer dtype/shape/bytes, exact clearance squared `25/64`, actual cavity count 1, all seven scene layers and unique z eigenspace. It observes the actual loaded extension before checking run freshness. Both runs have 39 ordered solve records, 72 input cuboids and matching locked native input/trace/output/mesh hashes.

`build/task-7v1-controller-pair-01.log` SHA256 `69f3f85e9fc7f47ee004f2c6594f3af91d31576a73235ca058dbc31a7d8c368d` records PASS. Maximum normalized mass-moment discrepancy `1.3877787807814457e-17`. Scientific material digest `778f8dfb1dd0dc601152e2ffb78d5756efe8c42b1e61c1b3c9c8c03af9f0a97c`; bridge scientific digest `780431ce4e1970da972488483a575b2b566a61fa072fdc6d09b65b5399c6e6a7`.

The following files are pairwise byte identical (paths/elapsed observations are separately retained):

| File | SHA256 |
|---|---|
| inputs.json | c415ea8f44faa23ee2c2e7e77ca3ed3ffbf8c8a77c56662893948097b1ba3978 |
| native-run.json | fa943e91287b47254d38a8a2f0f0fb4aedfe990f681ab4d21cc30c22cdb50dd9 |
| scientific.json | 80cc6eb362b1f46cbcfc7e3094f69860cdb5f79dbf4b8dfdc05d6fa3e3d10227 |

Direct native `--check-sealed-spin-bundle` on fresh run1 exited zero without DISPLAY, LD_PRELOAD or LD_LIBRARY_PATH; computed JSON agrees with the independently checked scene, full inertia, seven layers and scope labels.

## Real desktop observation

Actual X11 `DISPLAY=:1`; physical DP-0 1920x1080 at 144 Hz; direct NVIDIA RTX 5090 D OpenGL 4.6. No Xvfb, virtual desktop, system package install or host display setting change. Native process observed PID 2219798, X11 client `0x560000b`, exact title `Sealed spin | HOLLOW_ONLY | PARTIAL_BLOCKED`. The window manager maximized the requested window: client1854x1011 at(66,69). Shell session96220 owns the interactive viewer; no optimization continues in it.

System Pillow/XCB captured the named task-window rectangle. The initial `default-cutaway.png` was occluded by Codex and is NOT rendering evidence. An initial guarded retry refused because another window was active; no image was produced. The helper was tightened to verify active-window identity, then activation and capture were performed together. These are controller capture corrections, not production viewer changes.

Controller visually inspected:

- `activated-cutaway.png`: original left, exposed red cavity and retained gold face right; readable labels and degeneracy disclosure.
- `full-material.png`: actual checkbox disabled cutaway; both complete outer silhouettes agree under shared camera, cavity correctly hidden inside opaque material.
- `isolated-cavity.png`: actual checkbox hides retained material, red cavity alone remains with the CoM/axis/point overlay.
- `exploded-cavity.png`: isolated mode disabled, cutaway enabled, explode enabled; only the red display copy moves; source pose/CoM stay unchanged.
- `restored-cutaway.png`: explode disabled again. Host notification overlays partly occlude this capture; it is not used to prove every disclosure is visible. Earlier unobscured panel captures provide that evidence.

Clicks were scoped to the named active task window using XTest, at visually observed widget coordinates. No unrelated window controls were operated. Host notifications/desktop companion/IME may appear over screen captures and were not removed or silently retouched. A separate `task-window-native.xwd`/PNG attempt produced a partially black/occluded result and is explicitly NOT acceptance evidence. Those failed observations are retained, not presented as clean screenshots.

This is display triangle filtering with NO section caps, not a production Boolean section or simulation. Only a unique z axis is shown because transverse eigenvalues repeat. Input is one support point, explicitly not a certified support hull. Display does not grant solver/target/physics/original-model success.

Before/after-interaction snapshots cover all18 run1 seal files. `cmp` exited zero; both SHA256 maps equal `890a621a363f5214cb79a43d38692fc9dc968b48057ecaa35fa2b08076e03055`.

| Actual screen capture | SHA256 |
|---|---|
| full-material.png | 2f408a53d5e5dc8eaabfba1b2663c01bc89b5c211bfc1600c51aac9ff41c41d9 |
| isolated-cavity.png | bd9cfa49cd533a3230850532d921865dd09bd552c94e5a5053635098fbd8dd11 |
| exploded-cavity.png | bdc63bad8def5baab4ef96cf39a471108f2cf126f28fc67c61d557e33c75629d |

## Scope and preservation

Independent code review approved the bounded implementation, with no actionable findings. Controller real-desktop acceptance adds observation, not solver execution authentication. Runtime remains explicitly mixed system/Conda and only observed, not a complete attested dependency closure. GUI-dependency contents before the initial configure refresh were not recorded, so their historical preservation cannot be proven retrospectively.

Main tracked diff SHA256 remains `93e32af47ab1957a1efcb6716437656be129053f4269d2169378fcd95b05a35e`. No writes/merge into dirty main, no old experiment deletion, no push, release tag or complete-reproduction status. Historical Task7a report remains unchanged. Task7v1 controller docs/helpers are separately scoped; cumulative progress is not staged.

Remaining full cage/BBW deformation, dual-material, Faster2D/hydrostatic/CAD and remaining quadric cases, original assets/critical parameters, complete original experiment matrices, authenticated trace replay and independent simulation are not closed by this increment.
