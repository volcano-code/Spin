# Spin-It 2014 formula-reproduction report

Status: **formula core complete; reconstructed-asset milestone PASS**.

This is not an author-code or author-asset claim.  The paper's original model
caches, interactive cage choices and exact per-model parameters are not public,
so the strongest honest evidence class for the experiments below is
`RECONSTRUCTED_ASSET`.

## Implemented paper components

- Adaptive octree fill field with fractional-cell splitting and binary-sibling
  merging.
- Exact ten raw mass moments, analytic objective/constraint derivatives, the
  paper's static, spinning-top and yo-yo objectives, and active-set
  SLQP/L-BFGS optimization.
- Deterministic binary rounding, sealed material/cavity extraction, independent
  geometry-to-moment validation and immutable Python/C++ artifacts.
- Two-density boundary/interior material optimization, including an exact
  equal-density reduction to the one-material formulation.
- Joint cage deformation with scale/translation handles, analytic surface
  moment differentials, fill Laplacian, mesh bi-Laplacian and displacement
  regularizers.
- Deterministic sparse cage sampling from the voxel boundary, bounded
  biharmonic weights and interpolation to the outer surface and void cells.

## Verified results

The 4^3 deterministic cross-mode oracle was executed twice per mode from source
commit `f95463e6397e885bcadc7302020a026fc7b80540`:

| Mode | Status | Binary equality residual | Repeat output hash |
|---|---|---:|---|
| Static balance | `CONVERGED` | 0 | identical |
| Spinning top | `CONVERGED` | 0 | identical |
| Yo-yo | `CONVERGED` | 0 | identical |

The complete machine-readable record is in `modes/summary.json`.  The combined
configuration (Make-It-Stand + reachability + Spin-It Faster + Spin-It 2014 +
Python) passed **159/159** registered tests.  The 2014-only label passed
**12/12**, and the combined label passed **13/13** including Python.

## Real-desktop OpenGL evidence

The viewer was executed on the host GNOME Xorg desktop (`DISPLAY=:1`), not an
Xvfb/virtual display.  The 16^3 static fixture is in `viewer-smoke/`; top and
yo-yo captures are stored beside their mode artifacts.

The new `Explode cavity for inspection` control translates only the displayed
cavity.  It does not alter the saved geometry.  Red is removed material, blue
is the design domain, and gold is retained material.  `outer.obj` and
`material.obj` have the same outer boundary in the symmetric fixture because
the cavity is completely internal; `cavity.obj` is the actual carving result.

Run the static result with:

```bash
DISPLAY=:1 build/spin-faster-viewer/make-it-stand --spin-result \
  experiments/spin_it_reproduction/spin2014/viewer-smoke/outer.obj \
  experiments/spin_it_reproduction/spin2014/viewer-smoke/design.obj \
  experiments/spin_it_reproduction/spin2014/viewer-smoke/material.obj \
  experiments/spin_it_reproduction/spin2014/viewer-smoke/cavity.obj \
  experiments/spin_it_reproduction/spin2014/viewer-smoke/diagnostics.json \
  "Spin-It 2014 Adaptive Hollowing"
```

## Remaining limits

- This milestone does not reproduce the paper's original examples or exact
  interactive edit histories; those remain blocked on unpublished assets.
- A 16^3 yo-yo stress fixture with adaptive refinement from level 2 reached the
  relaxation limit under default solver settings.  Failed valid requests now
  return a structured unsealed artifact instead of throwing during diagnostic
  serialization.  Larger/highly refined convergence therefore remains an
  engineering target, not a claimed result.
- No fabrication or physical validation is claimed.

See `release-gate.json` for hashes and exact gate outcomes.
