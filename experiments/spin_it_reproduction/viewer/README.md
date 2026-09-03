# Spin-It reproduction viewer

Build the viewer with both solvers enabled, then open a five-file sealed
bundle with `--spin-result OUTER DESIGN MATERIAL CAVITY DIAGNOSTICS TITLE`.

The checked-in examples are:

```bash
./build/spin-faster-viewer/make-it-stand --spin-result \
  experiments/spin_it_reproduction/viewer/faster-buoyant/floating-island-cat/outer.obj \
  experiments/spin_it_reproduction/viewer/faster-buoyant/floating-island-cat/design.obj \
  experiments/spin_it_reproduction/viewer/faster-buoyant/floating-island-cat/material.obj \
  experiments/spin_it_reproduction/viewer/faster-buoyant/floating-island-cat/cavity.obj \
  experiments/spin_it_reproduction/viewer/faster-buoyant/floating-island-cat/diagnostics.json \
  "Spin-It Faster — Floating Island Cat"

./build/spin-faster-viewer/make-it-stand --spin-result \
  experiments/spin_it_reproduction/viewer/spin2014/ellipsoid-top/repeat-1/outer.obj \
  experiments/spin_it_reproduction/viewer/spin2014/ellipsoid-top/repeat-1/design.obj \
  experiments/spin_it_reproduction/viewer/spin2014/ellipsoid-top/repeat-1/material.obj \
  experiments/spin_it_reproduction/viewer/spin2014/ellipsoid-top/repeat-1/cavity.obj \
  experiments/spin_it_reproduction/viewer/spin2014/ellipsoid-top/repeat-1/diagnostics.json \
  "Spin-It 2014 — Ellipsoid Top (proxy)"
```

Gold is retained material, red is the removed cavity, transparent blue is the
fixed shell/design domain, and RGB lines are principal axes. Buoyant cases use
the target-frame plane as the waterline and show displaced volume, metacenter
height and stability margin. The 2014 viewer bundle is intentionally the small
level-4 proxy for interactive inspection; the release gate uses the separate
level-9, 512-cubed result matrix.
