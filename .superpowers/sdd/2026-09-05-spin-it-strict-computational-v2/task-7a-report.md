# Task7a: exact Faster Section 6.6 computational checks

Workspace `/data/xvxinyao/make-it-balance-paper-sim`, BASE
`14e39ca16b1d4cb28d99adb9952d792f71553746`, branch
`spin-it-strict-computational-v2`. No subagents, main-checkout edits, merge, tag, push,
author contact, historical result edits, or historical specification status edits.

## Outcome and scope

Implemented a bounded standard-library `Fraction` checker for exactly the two published
one-dimensional Section 6.6 definitions and a thin exclusive-output runner. This is not
a 3D proxy, ordinary Newton run, CGAL/CAD/native-extension computation, physical
simulation, whole reproduction, M4 completion, release gate, or proof that the two cases
exhaust all problems with p=0. The historical specification remains byte-identical with
its `NOT_VERIFIED` / `NOT_RUN` fields.

Each scoped case reports `VERIFIED_ANALYTIC_PROPERTIES`; the generated data separately
reports global `PARTIAL_BLOCKED`, author bitwise exact `NOT_VERIFIED`, physical validation
`NOT_IN_SCOPE`, no aggregate status change, and no manual release gate.

## Exact computations and independent conditions

The checker parses only exact integer/rational inputs, rejects floats as defining
coefficients, and performs polynomial multiplication/integration with `Fraction`; it uses
no samples, floats, `eval`, native modules, solver statuses, expected-outcome labels, or
caller-provided pass booleans.

For `faster.section6_6.nonunique` -> `faster.section6_6.example.1`, it evaluates
`chi_tau(x)=1/2+tau(2x-1)` for real `tau in [-1/2,1/2]`:

- exact volume polynomial `(1/2,0)`;
- bilinear rectangle corners `(1,0,0,1)`, with the explicit multi-affine
  corner-extremum argument proving `0<=chi<=1` throughout the rectangle;
- `integral_0^1 (2x-1)^2 dx=1/3>0`, hence distinct real parameters have positive
  squared L2 distance `(t-s)^2/3` and give an uncountable family modulo a.e.;
- independent integrals `1/2` for indicators of `(0,1/2)` and `(1/2,1)`, plus the
  nonbinary constant witness `chi=1/2`;
- indicators of `(t,t+1/2)`, real `t in [0,1/2]`, stay in the domain and have integral
  `1/2`. For ordered `0<=s<t<=1/2`, overlap is `1/2-(t-s)` and symmetric-difference
  measure is `2(t-s)>0`, whose exact affine polynomial is `(0,2)`. Thus this is also an
  uncountable binary family modulo a.e., not a finite witness count.

For `faster.section6_6.full_boundary` -> `faster.section6_6.example.2`, `|Omega|=1/2`
and `chi=1` integrates to `1/2`. At any optimum the box bound makes `1-chi`
nonnegative, while its integral is `|Omega|-1/2=0`; the
nonnegative-integrable-zero-integral lemma gives `chi=1` almost everywhere. The result
therefore records `UNIQUE_ALMOST_EVERYWHERE`, never pointwise uniqueness.

For both cases the square objective has global lower bound zero, each checked witness
attains zero, `2(v-1/2)=0`, and `p=-2(v-1/2)=0`. No normalized nonzero quadric is
demanded at p=0.

## Validation and fail-closed behavior

The checker accepts only the fixed schema, objective, no-additional-constraint contract,
case IDs, domains, and exactly two unique cases. It rejects missing/duplicate cases,
altered target/objective/domain/constraints/IDs, corrupted affine coefficients, truncated
parameter intervals, zero-width binary parameter family, changed full-domain witness,
duplicate JSON keys, and nonfinite JSON. A fake top-level `passed=true` plus a mutated
definition is rejected. Expected outcome/value labels are present for source-schema
identity but are not read when deriving the scientific outcomes.

The runner verifies the immutable specification and published PDF hashes before checking,
then creates a new output directory with `exist_ok=False`. A pre-existing output is
rejected. Evidence binds actual runner, checker, CMake configuration, specification, PDF,
result, report, Git HEAD and three binary diff contexts. The scientific digest excludes
clock/timing and output paths. Freshness re-hashes checker/spec/PDF/output bytes and
recomputes the scientific digest; tests show mismatched checker/spec/PDF records cannot
remain `CURRENT`. Full build attestation is explicitly
`NOT_APPLICABLE_PURE_PYTHON_FORMULA_PROOF`, not generalized from this formula proof.

## TDD evidence

RED was run before either production file existed:

```sh
PYTHONPATH="$PWD/python:$PWD/tests/python" python3 -m unittest test_faster_section_6_6_exact -v
```

It failed 6/6: the checker tests/freshness test reported
`exact Faster Section 6.6 checker is missing`, and the real CLI test reported the runner
file was absent. This was the expected feature-missing failure. The initial invocation
from the wrong cwd produced five `ModuleNotFoundError` errors and was corrected before
accepting RED; it is retained too. Raw relevant output: `task-7a-red.log`.
The original tool stream was not redirected to a file, so the persisted RED log is an
exact transcription of its test lines, representative traceback, summary, and the
intermediate command error rather than a byte-for-byte capture of repeated tracebacks.

GREEN used the identical repository-root command after implementation: 6/6 passed in
0.564 s. Focused registered CTest also passed 1/1 in 0.58 s. Raw output:
`task-7a-green.log`.

Behavior coverage uses independently supplied literal anchors. It mutates actual input
objects and runs the real CLI in subprocesses; there are no mocks or source-text tests.

## Preserved formula-evidence runs and source identity

After final source commit `266d8352e7aab8580e30c6dab5a7f2e0c51d3681`, two real
exclusive runs were retained:

- `build/task-7a-faster-section-6-6-run-1`
- `build/task-7a-faster-section-6-6-run-2`

Both freshness assessments are `CURRENT`, Git HEAD is `266d835...`, and recorded changed
paths are empty. Both scientific-data hashes are
`0ef670f30e2c6b08a00be6dae4abb4f590a389ba9b2eb671a109977b91f908e3`.
Their deterministic `result.json` files are byte-identical at SHA256
`a0bd1420d2aa471b557ea957c54feda9f9191d641ae5456c744f407aeee57782`;
their deterministic `report.md` files are byte-identical at SHA256
`ba33d1749e2d3ea1af9b951a1c258901c06c02d3ae74b189ef360ec0ee5163df`.
Evidence JSON differs as intended because runtime measurements are outside scientific data.

Actual input bindings in both runs:

| Input | SHA256 |
|---|---|
| checker | `d3377443944d68469bc4f74312cbb7b14dd8fe25afc9f3cc644d61d6df5aafff` |
| runner | `30c73a10b10d39be4aaeb268ad1581fd82e89eac4c1de47094758c23b450bdfd` |
| CMake configuration | `8b280acf93d0b5dbe015be98a535e78f71ed7755905f46d0bfa8c04573d50d67` |
| immutable specification | `053f896c2b10ece89b136bfba24cd480bf44cd0da99c0661e6e9bd34495c15fa` |
| published PDF | `e7fc6395362ef7e63a703851176afd81d8e953bae4e8080279d7e244b3c51fcc` |

The generated reports make the source/local comparison explicit and equal. Original case
counts outside these two are neither computed nor claimed.

## Final verification

The existing primary `build/strict-v2-spin` was reconfigured and grew from the dispatched
187 registered tests to 188 through one small Python CTest. Final source was committed
before the preserved runs and final verification.

```sh
cmake -S . -B build/strict-v2-spin
cmake --build build/strict-v2-spin -j4
ctest --test-dir build/strict-v2-spin -R '^faster-section-6-6-exact-contract$' --output-on-failure
cmake --build build/strict-v2-spin -j4
ctest --test-dir build/strict-v2-spin --output-on-failure -j4
```

- configure: exit 0; existing libigl old-CMake compatibility deprecation retained;
- focused all-target build: exit 0 with concurrency limited to four;
- exact registered CTest: 1/1 passed;
- final all-target build: exit 0 with concurrency limited to four;
- full primary suite, run once after final source: 188/188 passed, 0 failed, 65.90 s.

Final verification summary: `task-7a-final-verification.log`. The all-target build console
was not redirected; its exact live tool output remains in the task transcript, while the
persisted summary records its observed terminal lines and exit status.

### Full-suite raw-capture correction and recovery

The earlier statement that
`build/strict-v2-spin/Testing/Temporary/LastTest.log` at SHA256
`b9a4599ad670c92fa8aba343c4e9202e2ed9ec244fb688d6367390315b1e3231`
was a complete per-test log is **superseded and incorrect**. Inspection showed that file
contained only three lines: its start timestamp, a separator, and its end timestamp. Those
exact observed bytes are retained without relabelling at
`.superpowers/sdd/2026-09-05-spin-it-strict-computational-v2/task-7a-prior-lasttest-three-line.log`
with the same SHA256. There is no `task-7a-final-full.log` claim or file.

No authentic full capture of the original parallel run was recoverable. The precise reason
CTest's mutable `LastTest.log` contained no per-test lines is unknown. A controller-authorized
recovery ran one fresh suite without a rebuild, enumeration, or source change, with explicit
capture, serial execution, 24 GiB virtual-memory limit, and 60-minute timeout:

```sh
timeout 60m bash -c 'ulimit -v 25165824; exec ctest --test-dir build/strict-v2-spin \
  --output-on-failure --output-log \
  .superpowers/sdd/2026-09-05-spin-it-strict-computational-v2/task-7a-recovery-full.log -j1'
```

The recovery passed 188/188 with zero failures in 253.50 s. Its explicit output log has
448 lines, contains 188 `Passed` result lines, and is committed byte-for-byte as
`.superpowers/sdd/2026-09-05-spin-it-strict-computational-v2/task-7a-recovery-full.log`,
SHA256 `299a086e7c7f2bc7275db26a55f1f92e9eb9787640df0c4acf845a1f018c8a78`.
Immediately after recovery, CTest's mutable `LastTest.log` was again only three lines (now
with the recovery timestamp), reinforcing that it is not the full capture; its omission
mechanism remains unknown. The explicit `--output-log` is the authoritative retained raw
full-suite evidence. This recovery supplements rather than hides the original successful
parallel run, whose non-redirected observation and capture gap remain disclosed.

## Changed files and self-review

Source commit `266d835` changes only:

- `CMakeLists.txt`
- `python/make_it_stand/verification/faster_section_6_6.py`
- `experiments/spin_it_reproduction/strict/run_faster_section_6_6.py`
- `tests/python/test_faster_section_6_6_exact.py`

This report and the scoped logs are separate documentation/evidence commits. Build
outputs remain under ignored `build/` and are preserved for controller inspection.

Self-review read every new file and the CMake diff; rechecked the exact formulas, a.e.
language, independent label handling, two-case boundary, source bindings, output
exclusivity, deterministic/non-deterministic separation, no-native scope, freshness
recomputation, and restricted statuses. `git diff --check` and `py_compile` passed.

Concern: the pre-existing configure step emits libigl's CMake deprecation warning. It is
unrelated to this pure-Python formula checker and was not suppressed or broadened into
dependency work. No correctness concern remains within the bounded Task7a scope;
independent review belongs to the controller.
