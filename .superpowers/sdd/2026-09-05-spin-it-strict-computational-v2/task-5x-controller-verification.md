# Task5x controller verification — actual guarded adaptive loop

Status: VERIFIED_BOUNDED_ADAPTIVE_ORCHESTRATION; overall PARTIAL_BLOCKED.
This document is evidence navigation, never a pass flag.
Base: d81bc808806fbff24fd40fad4c6d5e16e5dc3dc0. Worktree:
/data/xvxinyao/moment-guided-spin-v1/code/make-it-balance-paper-sim.

## Scope and preservation

Only Task5x actual-source adaptive orchestration and its synthetic verification.
No original-paper case, binary material, sealed geometry, feasible witness,
physics, continuous-map legality, or complete Spin-It reproduction is granted.
Faster is disabled in the assigned build and foreign pending work is excluded.
Protected main tracked diff before:93e32af47ab1957a1efcb6716437656be129053f4269d2169378fcd95b05a35e.
Foreign four-file tracked diff before:a2229f437e553cdd4c584dcd65959e5ac92eb2e574db0b0e07d296d35eeffd75.
Authorized identity volcano <xuxinyao307@gmail.com>; no new worktree/dependency,
merge, push, tag, deletion of historical records, or changes to main checkout.

## Controller protocol (fixed before native loop export)

Outer cube [.3125,.6875]^3, centered design [.425,.575]^3, source depth4,
optimization levels2..4, cage2/fine3, all26 actual coupled BBW handles,
density1.17, beta=.31. Static frame origin(.5,.5,.3), identity rotation.
Known-stationary cases have zero center/inertia/regularizer weights and are
KNOWN_SOLUTION_VALIDATION only. The nontrivial control uses center weight1
and origin(.49,.52,.3), maximum one SLQP iteration/two line-search steps,
and is SYNTHETIC_REGRESSION. No tolerance is changed to make it converge.

Six fixed cases: stationary_multi (four-pass allowance), stationary_pass_limit
(one pass), nontrivial_stop, initial_geometry_reject (clearance.2),
rebind_geometry_reject (rebind-only clearance.2), aggregate_limit (320MiB
known-buffer cap). Normal clearances/tolerances use explicit implementation
defaults. Same beta dimension across levels deliberately tests partition identity
rather than dimension-only change detection. This differs from the implementer's
rho1.3 rotated stationary fixture; do not combine their scientific outputs.

Independent checks use raw callbacks/guard records, actual objective-gradient
re-evaluation in the native exporter, separate Python KKT arithmetic, independent
eight-node BBW interpolation, integer ancestry remapping and rational signed
tetrahedral integration. Native loop success booleans are not used as proof.

## Preflight evidence

Baseline CTest exact regex: ^spin-it-2014-source-(cage-adaptive-rebind|cage-guarded-slqp-move|adaptive)$.
3/3PASS23.25s, executed by controller (tool session21802). Assigned existing
Faster-OFF build only. /data available1.4TiB; external24GiB/3600s, threads1;
compilation at most4 total jobs and heavy native runs sequential.

Build Graph MCP available this turn: preflight incremental6parsedfiles,
17nodes/350edges updated; total419files5641nodes60344edges, languages
bash/c/cpp/python, noerrors. Last-build metadata advanced to baseline d81bc80.
Unqualified rebind query was ambiguous across declarations/definition; qualified
definition resolves54callees. This graph aids navigation only; it includes
foreign tracked source read-only and has incomplete C++ member/import resolution.

## Current verified implementation and suite

Production:97f7e12c4a8152779f5ad0691325e261e63b769e. Fresh scoped review
Approved (task-5x-review.md), no Critical/Important/Minor findings. Review is
limited to this adapter and its upstream contracts, not the whole paper.
Implementer focused1/1PASS35.49s. Controller all-target Faster-OFF build finished
60links, exit0; then the exact related regex
`^(spin-it-2014-|spin-2014-|source-spin-material-python$)` passed45/45 in400.50s.
All45 actually ran; zero Not Run/skip/failure. This is NOT the entire project,
Faster, or original-paper experiment matrix. Native and Python source checks and
fresh exact-OFF/Faster-OFF dependency builds are included in this related group.

Core archive SHA256:0aaac93b1d0be39473fed75cf2a8adb2b3c116fb062a7724bebb67a1efb9fa41.
Native regression SHA256:7b537275b077cf8f8c47e0e5f127bc735eccc87fa548b3de09974b334bb8186f.
GCC13.3.0, CMake3.28.3, Python3.13.5, NumPy2.1.3; Release,
2014/exact ON, Faster OFF, LIBIGL_SKIP_DOWNLOAD ON. libigl commit
f06fb13014aadff7eb1d7423ac018620441a3657 plus actual compiler dependency hashes.

Controller preview6cases independently passed; maximum normalized moment error
4.6904363622493164e-16. Actual3pass known-stationary chain preserves8beta values
while changing partition identity twice. Nontrivial control MAXIMUM_ITERATIONS;
initial geometry rejection has no numerical result; failed rebind retains solved
endpoint but no candidate checkpoint. Fixed320MiB case stops BEFORE first guarded
solve reservation, with no solved or pending endpoint. No input was changed.

## Preserved controller failures and corrections

- Preview checker initially assumed320MiB would allow one pass. Actual raw
  records showed pre-solve reservation failure; corrected only the expectation,
  retained fixed budget and original failing check log.
- A real-export negative test deleted the failed rebind record; the first
  checker mistakenly accepted it. Added terminal-specific required evidence
  checks; all16 actual-export corruptions now rejected. The first failed run
  remains task-5x-controller-preview-negative.log.
- Runtime dependency inventory first parsed absolute `ldd` file headers ending
  in colon as library paths. Raw ldd output identified the cause; exclude those
  headers. The incomplete manifest is retained as
  task-5x-controller-inputs-attempt1-failed.sha256; the failed run log remains.
  Corrected manifest has3633entries and sha256sum --check passes.
- Old Task5w replay input check fails on CMakeLists/core archive exactly as
  expected: STALE_RESULT, historical bytes preserved. It is not new evidence
  for the current core.

The current input manifest binds actual Ninja-discovered core compiler inputs,
build/cache options, linked archives/native binary/compiler/interpreter, controller
helpers,170 code modules actually imported during independent verification,
and ELF dependencies resolved from the trusted binary/interpreters/imported native
modules. Builtin Python modules are bound through the interpreter. This is an
observed build/checker code+runtime closure, not a claim to archive the entire OS.
Reports/results are excluded from scientific-input hashes to avoid cycles.

Post-production graph update:3newfiles parsed,53nodes/567edges updated,
total422files5693nodes60832edges, noerrors, build metadata97f7e12. See
task-5x-graph.json. Subsequent evidence-only commits need not alter source graph
identities; the recorded metadata names the exact production commit.

Main/foreign tracked diff hashes were rechecked after production and match the
preflight hashes above. Task-owned production paths are clean; the whole worktree
remains intentionally dirty with excluded foreign Faster work.

## Frozen replay and handoff

Final replay exit0, build/task-5x-controller-replay-fXtWIE. Six fixed cases each
ran twice in separate processes; complete JSONL output byte-identical:
64c3a3fcba0967088b6cf5a420f40c0fc1f44ff998e933370b7a914fc3ba284b.
Exporter executable SHA:d7f2cf3d6685000b0699bb8154cb0ca4a49e89427c1699e6f2ef68c91d2f6032.
Both final independent checks and all16real-export corruption tests passed.
Elapsed25.86/26.11s; maximumRSS28896/29056KiB, exits0. Resource text is separate
from scientific output. Raw copies preserved as task-5x-controller-run1/2.jsonl
and task-5x-controller-resource1/2.log, with final replay log retained.
Frozen3633-entry input check passed again after replay. No source/runtime change
was hidden by rewriting the output hashes. Build/test/review conclusions all
refer to production97f7e12 and core0aaac93b..., not excluded dirty Faster work.

This closes Task5x only. Next: binary material revalidation from the actual
deformed source, recomputed moments/constraints and sealed-object integration.
Nontrivial optimized convergence, continuous BBW internal legality, previous
actual geometry-negative gaps, original assets/parameters, original tables and
figures, Faster remaining algorithms/applications/CAD and physical-model
simulations remain open. No release tag or original-paper status is granted.
