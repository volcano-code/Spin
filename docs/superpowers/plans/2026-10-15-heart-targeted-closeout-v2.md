# Heart Targeted Closeout Implementation Plan

> **For agentic workers:** Use superpowers:executing-plans to implement this plan task-by-task. Do not restore old queues or spawn agents automatically.

**Goal:** Diagnose the frozen Heart failure and deliver a validated, repeatable optimized material within the registered speed gates.

**Architecture:** Preserve the formal Cherchi producer and native Newton. Use saved same-call states for diagnosis; introduce only evidence-backed changes in independent candidates.

**Tech Stack:** Python, C++, Cherchi, GMP, MPFR, MuJoCo, Blender.

**Spec:** docs/heart-specialist-closeout-progress-20261015.md and the user-approved Heart closeout requirements.

## Global Constraints

- Frozen geometry, axis, target, density, units and 0.375 mm clearance remain unchanged.
- Raw residual <=1e-12; analytic moments <=1e-10; applicable original three-scale derivatives <=1e-6; independent terminal moments/constraints <=1e-6.
- Successful complete core run and repeat <=60 s each; cold preparation through seal/reload <=600 s, separately reported.
- Seed 0/1, eight attempts, sixty evaluations per attempt, alpha=tmax=0.5; no historical answer initialization.
- Large jobs serial, 24 GiB, swap 0, 60 minutes; compile <=4; numerical libraries single-threaded.
- Preserve all failures and uncommitted changes; no commits, pushes or deletions.

## Review Focus

- Attempt-final residual must not replace the best accepted residual in diagnosis.
- Real construction changes must not be relabeled as smooth retriangulation.
- A small-step decrease does not establish eventual convergence.
- Timing must include failed trials, backtracking and restarts.
- Prepared-state caches must not retain all trial materials indefinitely or reuse changed inputs.

## Task 1: Close the current run ledger before changing dependencies

- [ ] Read seed 1 result after its actual worker exits; preserve running bindings until then.
- [ ] Parse seed 0/1 events by attempt, evaluation, accepted residual, multiplier, step, restart and rejection stage.
- [ ] Compare the best accepted state, last accepted state and final returned state explicitly. Seed 0 currently reports OPTIMIZATION_DID_NOT_CONVERGE, normal 529.395436001 s, preparation 234.89493218 s; this is not a convergence or speed pass.
- [ ] Extract the first stagnating attempt and its first scientifically meaningful failure, with exact coefficients, same-call material and Jacobian identities.

## Task 2: Diagnose the actual stagnation, not model appearance

- [ ] Recompute Eq.14 residual, multipliers, tangent system singular values and linear-system backward error at the saved state.
- [ ] Compare predicted and actual residual changes along the actual Newton direction and registered backtracking steps. Reuse existing measurements rather than repeating the 81-point batch.
- [ ] Attribute changes to source construction, Boolean incidence, derivative conditioning or linear-system assembly using the same-call provenance.
- [ ] If a formula or coordinate defect is found, first add its saved-state failing regression, then make the smallest correction and rerun that regression plus unaffected ellipsoid/hyperboloid cases.
- [ ] If formulas are correct but the frozen initialization encounters a difficult basin, record that limitation. Do not silently add trust regions, change weights, shorten the Newton direction or relax convergence; an algorithm extension requires a separate decision.

## Task 3: Remove measured repeated cost

- [ ] Profile the saved stagnating state once cold and three times warm: source, Boolean, boundary, integration, derivative, Newton and persistence, with exclusive times.
- [ ] Inspect ownership of evaluate-cache and producer material snapshots; test that evicting rejected transient states preserves the accepted terminal attachment and scientific events.
- [ ] Keep the certified interval integral with exact fallback only where bitwise equivalence is established. The existing four-state whole-evaluation gain is about 12–18%, not enough for the final speed gate.
- [ ] Optimize the largest measured remaining phase in an independent build; verify material, source, moments, Jacobian and acceptance decisions before adoption.
- [ ] Run one complete seed only after a measurable candidate benefit or a demonstrated convergence correction; no unchanged hourly reruns.

## Task 4: Freeze and validate the successful path

- [ ] Execute the original seed 0/1 protocol and retain every failed attempt. Require real accepted updates, changed material, positive multipliers and raw residual <=1e-12.
- [ ] Seal the actual accepted same-call material, independently reload and recompute moments, objective, constraints and partition.
- [ ] Repeat the smallest successful seed in a new process from its original random initialization; compare scientific events, material and summaries, excluding clocks.
- [ ] Report paper-corresponding computation, complete normal solve, preparation and terminal processing separately. Speed failure remains explicit even if correctness passes.

## Task 5: Deliver only verified artifacts

- [ ] Export fixed, retained, cavity, shaft, components and whole assembly OBJ from the sealed material; independently reload and apply the registered geometric bounds and manifold/partition checks.
- [ ] Bind SI mass, center of mass, full inertia, roles and hashes in JSON.
- [ ] Only after collision-input admission run the approved MuJoCo 64-scenario comparison; Blender reads the actual trajectories.
- [ ] Publish separate correctness, repeatability, speed, OBJ and dynamics conclusions. Rendering cannot substitute for any gate.

## Budget and stopping rules

Run-ledger and failure attribution: 0.5–1 working day. Evidence-backed correction and same-state performance work: 1–3 days. Complete run/repeat and material delivery: 1–2 days after admission. Dynamics/rendering: 1–2 additional days after a valid material exists. These are conditional effort estimates, not a guarantee that the frozen input converges. If no implementation defect is found, deliver the precise limitation and request a decision before extending the paper algorithm.
