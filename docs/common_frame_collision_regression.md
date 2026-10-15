# Common-frame collision construction regression

The versioned implementation is `experiments/heart_stem_contact/bullet_contact_common_frame.py`. The historical constructor is intentionally retained so old experimental source identities remain meaningful. This change is not automatically selected by existing runners.

Portable test command, from the repository root:

```sh
python tests/test_bullet_common_frame.py
```

Requires numpy, scipy, trimesh and pybullet. Missing dependencies cause an error; no conditional success/skip is provided. The test is registered as `pybullet-common-frame-contracts` when `MAKE_IT_STAND_ENABLE_PYBULLET_TESTS=ON`, including the existing `physics-python` preset. Configuration checks all four imports. Enabling that group and running it alone cannot satisfy the full release test gate.

CTest integration was verified in a separate build directory:

```sh
cmake --preset physics-python -B build/common-frame-ctest -DPython3_EXECUTABLE=/home/lin/anaconda3/bin/python
ctest --test-dir build/common-frame-ctest -L pybullet --output-on-failure
```

Both registered PyBullet CTest entries passed, actually executing rather than skipping. This run did not compile or execute the unrelated C++ tests registered in the same directory; it is not a full-project green declaration. Existing CMake dependency deprecation warnings were not changed.

The synthetic asymmetric convex fixture is generated into a temporary directory. It checks that three COM/principal-frame descriptions of a common model produce the same engine collision vertex set in world coordinates. It also checks full world inertia, mass, COM, angular/linear initial velocity and no external fixed constraint. A separate negative test rejects an empty collision-part list.

Regression witness:

```sh
CONSTRUCTOR_MODULE=bullet_contact python tests/test_bullet_common_frame.py
```

The witness failed with a 7.72595 micrometer vertex-set mismatch. The new constructor passed both tests. Synthetic fixture results are engineering verification only, never original-paper case evidence or a surface-fidelity certificate.

Actual reconstructed-teapot checks are preserved separately under `build/pla-glass-20260920`: 24 support/material-frame subcases, 13-part frame consistency, free-rotation invariance, drop initialization and short-time integration diagnostics. Those records do not establish full spin-down convergence or complete reproduction. Required next work includes the full release audit, review of collision surface fidelity and bounded full-duration convergence with explicit source selection.
