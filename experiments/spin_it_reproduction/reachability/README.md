# Alternative-pose carving reachability

This directory classifies the three alternative-pose cases previously shown
through the OpenGL workflow.  The orientations are reconstructed parameters,
not recovered author interaction histories:

- `spheres-45deg`: official Spheres cache, reconstructed 45-degree gravity;
- `dinosaur-head-contact`: official Dinosaur cache, reconstructed headward gravity;
- `toy-head-down`: reconstructed watertight Toy surface and reconstructed inverted gravity.

These OpenGL comparisons used the prior zero-angle feasibility target; the
separate paper-simulation gate retains its 5-degree safe-region protocol.

The analyzer is a fixed-safe-region continuous relaxation.  An infeasibility
certificate proves only that no independent fractional removal of the eligible
interior cells can place the projected center of mass in that fixed region.
It does not rule out deformation, a changed contact region, a new external
pose, or a different manufacturing model.

Run from the repository root after building the Python preset with
reachability enabled:

```bash
LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6 \
PYTHONPATH=build/python/python:python python3 \
  experiments/spin_it_reproduction/reachability/run_cases.py
```

`cases.json` is the immutable case definition.  `results/summary.json`, the
per-case JSON files, and `results/REPORT.md` are generated deterministically.
`results-repeat/` is a separately computed byte-identical repeat.
