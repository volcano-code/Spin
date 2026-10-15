# Faster two-restart evidence parser and focused probe — 2026-09-24

## Scope

This is a development-only method diagnostic for one frozen reconstructed
box-yo-yo request, seed 2. It is not the five-seed repeated batch and is not an
author-paper case.

## Changes and checks

- The native diagnostic emits both deterministic start coefficient vectors,
  the per-restart evaluation counts, and every attempt-start, attempt-terminal,
  failed, or rejected solver event.
- The Python evidence parser enforces the frozen maximum of two starts, at most
  60 evaluations per attempted start, exact agreement between per-start and
  total evaluation counts, and exact cross-process equality of start and event
  records for the full-sweep runner.
- Focused evidence-parser tests passed **7/7**; the native target
  `spin-faster-unknown-seed-sweep` built successfully.
- A single-seed diagnostic completed in a 24 GiB memory scope with swap capped
  at zero and numerical libraries limited to one thread. The real transcript
  was consumed by the one-seed parser: two starts, 46 retained event records,
  two attempt starts, and two attempt terminals.

## Outcome

Seed 2's first start failed its initial Boolean boundary audit after one
evaluation. Its second start reached the 60-evaluation cap. The summary reports
`restart_evaluations=[1,60]`, 61 total evaluations, normalized residual
`9.1993847767608391e-06`, and no verified final geometry. The result is a
reproducible failure diagnostic, not convergence or an infeasibility proof.

Transcript: `build/faster-two-restart-seed2-20260924-a/seed2.stdout`

Transcript SHA-256:
`b66ee9af12c1fc9b9dc13e389da47b3b6da57e3f5dee6deaefda01f6adbae932`

Project reproduction status remains `PARTIAL_BLOCKED`.
