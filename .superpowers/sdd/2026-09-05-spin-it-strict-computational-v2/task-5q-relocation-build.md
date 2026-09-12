# Observed external relocation / dependency regeneration

Both project aliases became symlinks at 2026-09-12 11:44 to the corresponding
directories under `/data/xvxinyao/moment-guided-spin-v1/code/`. This task did not
perform the move. HEAD remained `6f61946`; protected main tracked-diff SHA and
six external Faster-file hashes remained unchanged.

The worker's already-requested assigned target build regenerated CMake. Its
cache had `LIBIGL_SKIP_DOWNLOAD=OFF`. The libigl download helper therefore
invoked its generated CGAL ExternalProject clone step. Observed child chain:
149210 cmake gitclone ->149213 git clone. The regenerated gitinfo command and
stamp paths now use the canonical relocated directory. Its generated comment
states that changed download information reruns the download. This is
consistent with relocation invalidating the old download metadata; the old
gitinfo bytes were not separately captured, so no exact before/after metadata
diff is claimed.

The controller requested termination of the owned build. Worker terminated
149065/149064/149041 (session92877). By that point cloning/regeneration had
finished and 16/73 compilation actions had run. Subsequent process inventory
showed no related build, compiler or clone jobs. Raw `task-5q-build.log` remains
in this SDD directory; it is an interrupted attempt, never a passing build.

Important side effect: generated
`libigl/external/.cache/cgal/cgal-download-prefix/tmp/cgal-download-gitclone.cmake`
removes the old `libigl/external/cgal` directory before clone. This action had
already completed. Therefore **no deletion-free claim is made**. This task did
not create a backup of the old dependency working directory. Current CGAL is
clean at the same previously pinned commit
`f7c3c8212b56c0d6dae63787efc99093f4383415`; Eigen is clean at
`21ae2afd4edaa1b69782c67a54182d34efe43f9c`. This verifies present source identity,
not preservation of any unknown old untracked dependency data.

Follow-up: Eigen download stamp files were also dated 11:45. Re-creation of
Eigen during the same regeneration cannot be excluded; there was no separate
live clone observation for it. Current Eigen identity is verified as above.
Do not infer that only CGAL's directory was affected.

Ruling: reconfigure only the existing assigned build with
`-DLIBIGL_SKIP_DOWNLOAD=ON`, retaining all scientific/test feature settings.
The helper guards only downloading; it does not disable CGAL or checks. Keep
all logs/caches/partial objects; do not repair Git, delete caches, change pins,
or touch external project edits. Then resume the bounded build and tests.
Cost if wrong: diagnose remaining path/cache mismatch; no relaxed scientific
gate, no unverified dependency fallback. User was informed of actual paths,
dependency replacement and backup limitation.

Further bounded prevention: fresh registered dependency-contract child builds
must also pass SKIP_DOWNLOAD=ON. A real generated fixture rejects OFF first;
then two explicit arguments correct both child invocations. All scientific
feature/build/runtime and retained-directory assertions stay intact. This
addendum avoids relying on potentially invalidated stamps to prevent downloads.

## Alias rejection in evidence tests

First post-feature scoped suite passed the new native moment test, then
`source-spin-material-python` failed: generated CTest PYTHONPATH still used the
alias, so imported Python `__file__` passed to `implementation_observation()`
contained a symlink. `run_evidence._path` correctly rejected it. This is not a
reason to weaken provenance checks. Partial suite is preserved as
`task-5q-safe-scoped-suite.log`, never promoted to all-green.

Ruling: stop that owned suite and configure the existing build using canonical
absolute source and binary paths with download suppression still ON. Inspect
generated paths, rebuild affected outputs and rerun focused Python plus a new
scoped suite. If CMake refuses the source-path change, diagnose without editing
or deleting the existing cache. No scientific validation source is changed.
Cost if wrong: another recorded path/build attempt, not relaxed acceptance.
