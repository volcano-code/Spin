# Make-It-Stand paper-simulation closeout

This directory is the fail-closed experiment layer for the isolated
`make-it-stand-paper-sim-v1` branch.  It never treats an old result as current
unless its source tree, executable, configuration, dependencies, and input
hashes still match.

## Evidence classes

- `EXACT_CACHE`: untouched official OBJ/MIS/VOX/BBW input and no reconstructed
  parameter.
- `RECONSTRUCTED_PARAMETER`: complete cache, but an unarchived parameter such
  as lambda was selected by the documented bounded search.
- `RECONSTRUCTED_ASSET`: an input artifact had to be regenerated.  Horse uses
  an occupancy-derived outer surface for simulation.  The rocking-toy OBJ is
  an open 29-component surface after indexed cleanup, so its runner first makes
  a deterministic voxel-union/marching-cubes shell, then rebuilds voxels, BBW,
  pose, support threshold, and handles.
- `BLOCKED_MISSING_ASSET`: the requested experiment cannot be run honestly.

The public archive remains external and read-only.  Its files and claim
ceilings are listed in `asset_ledger.json`.

## Execution order

Use the Python extension produced by the `python` preset:

```bash
export PYTHONPATH="$PWD/python:$PWD/build/python/python"
export LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6
export MAKE_IT_STAND_BETA_ROOT=/home/lin/xvxinyao/reference-only/make-it-stand/beta-extracted/make-it-stand_src_beta
```

1. Search lambda independently for spheres, dinosaur, armadillo1, and horse
   with `run_lambda_search.py`.  The fixed coarse set is
   `{5,10,20,40,80}`; a successful neighbor is refined geometrically.
2. Put the selected five values in a JSON object keyed by case, then run
   `run_table1_matrix.py`.  The enforced order is spheres, dinosaur,
   armadillo1, horse, toyrocking; each case runs twice at full resolution with
   24 GiB and 60 minute limits.
3. Run `run_figure14.py` on official N=80 spheres.
4. Run `run_static_oracle.py` on the five current sealed outputs.
5. Build one multi-case simulation manifest with
   `prepare_pybullet_manifest.py`, then run `run_pybullet_matrix.py`.  This is
   18 unique scenarios per case and two repeats: 90 unique case-scenes and 180
   simulations in total.
6. Invoke `release_gate.py` with the five case roots and the three reports.

`PAPER_PROFILE` and `RELEASED_BETA_PROFILE` are never combined.  Every
optimizer run emits `iteration-trace-v2.json`, which separately names paper
outer iterations, beta optimize calls, carving calls, deformation trials,
target generations, objective contexts, rollback decisions, and state hashes.

## Geometry and simulation contract

Only `SealedBalanceArtifact` may enter PyBullet.  It applies the archived
support projection, re-audits the exact flattened outer and cavity boundaries,
repairs only explicitly audited structured-cavity and collapsed planar-support
defects, measures outer/cavity clearance and containment, then reintegrates
mass, center of mass, and inertia from that final material geometry.  A raw
optimizer mesh, incomplete material-shell audit, or pre-flatten moment payload
is rejected.

`quasi_static_table_tilt` keeps the designated support polygon fixed and tilts
gravity.  `reseeded_initial_pose_drop` rotates and lifts the object before a
dynamic release.  They answer different questions and are reported separately.

The release gate additionally requires a clean tree and configured Git name
and email.  Until the user supplies that identity, no final commit or tag may
be created.  Even after every computational gate passes, the strongest claim
is `VERIFIED_PAPER_SIM`, never author-identical or bitwise reproduction of
unpublished assets and interaction history.
