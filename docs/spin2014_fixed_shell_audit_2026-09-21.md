# Fixed boundary shell versus complete feature thickness — 2026-09-21

## Source and scope

Spin-It (2014), Section 5.1 and Figure 3, partitions the original solid into
fixed boundary shell and interior, and restricts hollowing to the interior.
Section 5.3 initialization only counts the interior contribution of cells that
overlap the shell. These are input-domain and material-preservation requirements,
not just a check that a particular final cavity happens to be distant from the
outer surface.

The relevant complete PDF pages 4–5 were extracted and visually inspected from
`/home/lin/下载/SpinItPaper.pdf`, SHA-256
`3add70de03df9b0e173028d55c775e8520fa89a6b21d9b269f926ccbe8a76e4e`.
The project separately requires a complete thin-feature/thickness engineering
audit. That stronger gate is retained. Neither fixed-shell preservation nor
outer/cavity clearance replaces it. The previously studied material-ball-cut
metric is not adopted as a release metric by this change.

## Implemented independent geometry sidecar

`python/make_it_stand/verification/fixed_boundary_shell.py` checks real exact
polyhedral regions, independently of the optimizer and its occupancy variables:

1. The locked original domain is partitioned by fixed and design regions.
2. Every entire design region has at least the locked minimum clearance from
   the original boundary, including when nothing was removed in the result.
3. Every design region partitions into its actual retained and removed regions.
4. Final material is exactly the union of fixed and retained regions; reported
   removal is exactly the union of removed regions.
5. The reported outer vertices and faces remain canonically byte-identical.

Inclusive exact material partitions establish containment before the distance
check; distances alone would not distinguish an interior region from an
exterior one. Every returned result includes input/output identities and native
partition/clearance facts. No geometric repair or parameter change is performed.

The narrow success state is `FIXED_BOUNDARY_SHELL_VALIDATED`. Full feature
thickness stays `UNKNOWN`; material admission, original-paper reproduction and
physics are not granted. This sidecar performs no new mass integration and does
not alter the frozen bundle format or upgrade historical sealed facts.

## Requirement coverage correction

The strict requirement catalogue previously omitted this distinct Section 5.1
constraint. Revision 2 adds `spin2014.hollowing.fixed_shell`, bringing the list
from 88 to 89. All original 88 rows are unchanged. The row-layout schema remains
v1, while its content revision is 2. Old evidence omitting the new row fails
validation; it is not auto-migrated to a pass. A caller-written `VERIFIED` also
fails. The new checker is not yet a registered strict original-paper requirement
checker: reconstructed G7 geometry alone cannot supply the required parameter
and original-case provenance.

The preceding catalogue, tests and CMake file are preserved under
`build/spin-fixed-shell-audit-20260921/before/`. No historical experimental log,
original input, optimizer, objective, threshold, or rounding rule was changed.

## Regression evidence

Two new strict-evidence tests initially failed because the requirement was
absent; the six initial geometric tests failed because the checker was absent.
After implementation, 22 tests passed (six native geometry plus 16 strict audit).
Read-only review found no actionable correctness issue and reproduced these
tests. Reviewer probes for a nontrivial retained/removed split, incorrect
aggregate removal and overlapping fill were then added to the native suite;
all eight native tests passed on the seal-bound extension.

Checks include equality at the exact clearance threshold, rejection at the next
binary64 value above it, no-removal input with an insufficient shell, missing
fixed material, changed outer geometry, missing/duplicate/unknown domain IDs,
out-of-domain removal, and overlapping or incorrectly aggregated results.

## Actual G7 sidecar execution

The first source-bound job stopped before geometry acceptance because its
launcher selected the current development native module, whose hash differs
from the module recorded in the historical seal. Its failed receipt and script
are retained in `job/` and `before/run-v1.py`. No hash check was disabled.

The corrected launcher selects the original seal's native module explicitly:
`build/spin2014-g7-direct-audit-candidate-20260921/native-audit/python/`, SHA-256
`bdc914b7eef0b4ac578dc116918bf4c37b7ce6c8227de5b4176fb9c5dfb2c01f`.
The development native module used by CTest is independently included in the
job identity as well. The new job is `job-native-pinned/`. At this report's
checkpoint it is still running: the observed service is
`spin2014-job-29148c3eb7094d1fb225e83d23adc3f7.service`, with geometry worker
PID 2287960. The first `fixed-shell-0.json` has now returned
`FIXED_BOUNDARY_SHELL_VALIDATED`, including the whole input design domain's
clearance at the unchanged 0.015 model-unit threshold. The second check and
required registered tests are still pending. There is no completed whole-job
receipt or two-repeat acceptance yet; a later completion addendum must inspect
the actual terminal receipt, not infer success here.

Both jobs retain 24 GiB memory, swap 0, a 60-minute complete-job deadline,
single-thread numerical libraries, and the existing disk-space safety limits.
No optimization or physical simulation is launched here.

## Final local regression checkpoint

After adding the review probes, the current eight native geometry tests and
16 strict-audit tests ran together: **24 passed in 3.986 s**, no skips. This is
a local regression execution, not the still-pending source-bound G7 job.
The native geometry suite also passed separately against the original seal's
native extension. Read-only review of the new batch runner found no actionable
provenance or correctness issue.

Current SHA-256 observations:

- fixed-shell checker: `ceaed1c4ef56471e6bc8f568ec86b10b7aa3032aabdf89851c7b2a21140ecede`
- fixed-shell tests: `6ae7049c82adc9a2bba784792da91f09809aae02bb90bbeea8fe96eaab76bc9f`
- requirement catalogue: `45e1dffdbfe97c04095e1d619b7349ada0c570468cbdb2301e5feadcb7db329e`

The G7 check uses the seal-native extension; registered CTest uses the separately
pinned development-build extension. File hashes are observed disk provenance,
not loaded-memory or complete build-source attestation.

## Completed source-bound G7 check

The previously live `job-native-pinned/` exited successfully in
752.820779065 s, peak memory 637,464,576 bytes, swap 0. Both fixed-shell outputs
are byte-identical. Both required CTest groups ran and passed: eight fixed-shell
tests and 16 strict-audit tests, with no skipped groups.

The actual complete input design region passes the locked 0.015 model-unit
clearance. Exact initial, per-design, final-material and removal partitions all
pass. This establishes the fixed-shell geometry fact for these reconstructed G7
inputs and outputs, not provenance of an author-original shell parameter.

The complete whole-job receipt, two outputs, checker digest and test logs are
under `build/spin-fixed-shell-audit-20260921/job-native-pinned/`. A later CMake
registration change for the separate critical-speed oracle does not retroactively
update this historical job's frozen CMake identity. Full feature thickness remains
`UNKNOWN`; physics/material and strict-original-case admission remain ungranted.
