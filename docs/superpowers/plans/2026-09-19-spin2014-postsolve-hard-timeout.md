# Spin-It Post-Solve Worker Hard-Timeout Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Ensure the approved G7 diagnostic controller returns within its 180-second per-call wall budget even if its worker ignores graceful termination or sends a result but does not exit.

**Architecture:** Keep the historical diagnostic runner immutable. Harden only `experiments/teapot_spin2014/postsolve_replay.py`: reserve a bounded cleanup window inside the call deadline, use only timed joins, escalate from terminate to kill, and never report a successful diagnostic unless its worker has exited cleanly. Test timeout/cleanup behavior with deterministic fake process objects; do not invoke the native solver or G7 replay.

**Tech Stack:** Python 3.13, `multiprocessing`, pytest.

**Spec:** `docs/superpowers/plans/2026-09-19-spin2014-g7-postsolve-replay.md` and its global resource constraints; the 180-second diagnostic cap must include child cleanup.

## Global Constraints

- Work only in the existing isolated worktree `/data/xvxinyao/moment-guided-spin-v1/code/make-it-balance-paper-sim`.
- Preserve all existing dirty paths and historical experiment directories.
- Do not modify or execute the archived out-of-plan `build/spin2014-g7-poststage-timing-20260919/run_diagnostic.py`.
- Do not invoke the optimizer, resume a checkpoint, import the native solver for a real G7 replay, change scientific settings, or run a large build.
- Keep the total `_bounded_native_call` wall time at or below its supplied budget; address-space remains 24 GiB and numerical-library threads remain one.
- Tests must verify timeout cleanup without launching native code or touching frozen G7 artifacts.
- Do not commit mixed changes from an already dirty file; document why this task remains uncommitted if earlier hunks in the same files cannot be safely separated.

---

### Task 1: Bounded process cleanup for the supplemental diagnostic

**Files:**
- Modify: `experiments/teapot_spin2014/postsolve_replay.py`
- Modify: `experiments/teapot_spin2014/test_postsolve_replay.py`
- Create: `.superpowers/sdd/2026-09-19-spin2014-postsolve-hard-timeout/task-1-report.md`

**Interfaces:**
- Keep `_bounded_native_call(callback)` backward-compatible; an optional private-test wall-budget argument may be added, with production default `WALL_SECONDS`.
- A private worker-cleanup helper must return whether the process exited, must never call `join()` without a finite timeout, and must escalate `terminate()` to `kill()` before the total call deadline.
- Return `UNAVAILABLE_RESOURCE_LIMIT_WALL` when the worker is reaped after the call budget; return `UNAVAILABLE_NATIVE_CLEANUP_FAILED` with the child PID when it is still alive after SIGKILL and the bounded reap deadline. Never accept a payload from a worker that failed to exit cleanly.

- [x] **Step 1: Write failing tests**
  - Add a fake process that remains alive after `terminate()`, records every `join(timeout)` and exits only after `kill()`; assert cleanup escalates and all joins have finite timeouts.
  - Add a controller-path test with a fake pipe that produces no payload and a short injected wall budget; assert timeout status, escalation, endpoint closure, and no unbounded wait.
  - Add a success-path test in which a fake worker provides a valid compact result but does not exit until joined; assert that result acceptance occurs only after bounded worker exit.
- [x] **Step 2: Run the new focused tests and confirm they fail because the existing controller performs an unbounded `join()` after `terminate()` and after receiving a payload.**
- [x] **Step 3: Implement deadline accounting and bounded cleanup.** Reserve enough of the supplied wall budget for terminate/reap and kill/reap; replace every unbounded `join()` with deadline-bounded waits; close pipe endpoints on every branch; report cleanup failure instead of success if the worker remains alive.
- [x] **Step 4: Run `PYTHONDONTWRITEBYTECODE=1 python3 -m pytest -q -p no:cacheprovider experiments/teapot_spin2014/test_postsolve_replay.py`; confirm all tests pass and no native or optimizer invocation occurred.**
- [x] **Step 5: Run `git diff --check`; append test evidence and the reason no commit was created to the task report. Do not stage shared dirty paths.**

### Bound of this controller

The controller enforces deadline-derived waits and cleanup once the synchronous `os.fork()` call returns. No in-process implementation can preempt a stalled fork syscall or OS scheduling delay; an independent outer watchdog is still required for an absolute whole-process wall cap. This remaining boundary is recorded in `.superpowers/sdd/2026-09-19-spin2014-postsolve-hard-timeout/task-1-report.md` and must not be represented as a proven hard real-time guarantee.
