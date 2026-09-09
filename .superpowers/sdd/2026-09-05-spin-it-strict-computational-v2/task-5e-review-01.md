# Task5e independent review — 20aa427

Reviewer `/root/review_spin2014_source_material_boundary`, read-only full BASE dba17bd..20aa427 diff review. Verdict: spec issues found; quality Needs fixes. No Critical findings.

Important: `include/verification/ExactCuboidUnion.h:20` excludes reused auditor explicit application-level Nef objects/Boolean results, not merely allocator internals. Component construction/nesting/signed assembly/relations/partitions in `src/verification/ExactMaterialRegion.cpp:139,145,176,230,277` do not enforce configured intermediate complexity. Zero carving with maximumIntermediateComplexity=4 passes empty-union accounting but subsequent ordinary cube region audits exceed the cap and can publish success. Propagate scoped budget through reused operations, including live operands/results; add real regression for this path.

Minor/report correction: report says full passing-test stdout retained in `build/task-5e-full-ctest-detailed-01.log`, but controller and reviewer read exactly three start/end lines. Full summary197/197 and focused verbose evidence remain valid. Preserve the incomplete file and correct claim; do not regenerate history or rerun a suite just to replace this artifact.

Reviewed strengths: immutable owned computed geometry, private actual solver terminal handoff, signed nested-shell/retained-island handling, exact inclusive partition with separate interior overlap, failure publication ordering, cap and clearance equality tests, actual success/failure solver regressions. Historical preservation and runtime identity require controller evidence. No review builds/reruns/writes.

Controller ruling: accepted the concrete resource gap; earlier instruction to avoid speculative broad auditor refactoring did not authorize cap bypass. A narrowly scoped opt-in budget propagated through existing exact auditor operations is within the approved task and preserves legacy behavior when absent. No tolerance/output-cap relaxation or hard-allocator claim. Fix round1 returns to original implementer. Report-only correction also required, despite Minor severity, because evidence claims must be accurate.

Controller postcommit checks at20aa427: main focused4/4 in58.10s (`build/task-5e-controller-focused-01.log`), all source/bin/cache/evidence manifests match. Actual Python-loaded extension path is `build/strict-v2-spin/python/_make_it_stand.cpython-313-x86_64-linux-gnu.so`, SHA0f58767d04085af22fdbab2e7e6bde03eb9dd045d364b0886a9c48482948ed9f. These precede the resource fix and must not be relabeled final amended-code verification.
