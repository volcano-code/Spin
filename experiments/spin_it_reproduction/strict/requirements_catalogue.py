"""Versioned, fixed requirement catalogue for the strict Spin-It audit.

This module contains catalogue data only.  Validation and evidence evaluation
belong in :mod:`audit`; callers cannot replace this list with selected rows.
"""

from __future__ import annotations

from typing import Any


CATALOGUE_SCHEMA = "spin-it-strict-requirements-v1"
CATALOGUE_VERSION = 1


def _row(
    requirement_id: str,
    paper_locator: str,
    requirement_kind: str,
    expected_outcome: str,
    required_evidence: str,
    dependencies: tuple[str, ...] = (),
) -> dict[str, Any]:
    return {
        "requirement_id": requirement_id,
        "paper_locator": paper_locator,
        "requirement_kind": requirement_kind,
        "expected_outcome": expected_outcome,
        "dependencies": list(dependencies),
        "required_evidence": required_evidence,
    }


_ROWS = [
    *[
        _row(
            f"spin2014.eq.{number}",
            f"Spin-It (2014), Equation {number}",
            "FORMULA",
            "MATCH_PUBLISHED_EQUATION",
            f"Independent symbolic and numerical checks of Equation {number} with published conventions and tolerances.",
        )
        for number in range(1, 9)
    ],
    _row("spin2014.objective.static", "Spin-It (2014), static balancing objective", "FORMULA", "MATCH_PUBLISHED_OBJECTIVE", "Independent value, gradient, and Hessian checks for the static objective."),
    _row("spin2014.fill.continuous", "Spin-It (2014), Sections 5.2-5.3", "METHOD", "MATCH_CONTINUOUS_FILL_SEMANTICS", "Trace demonstrating continuous fill variables and published bounds."),
    _row("spin2014.fill.binary", "Spin-It (2014), Section 5.3", "METHOD", "MATCH_BINARY_SEMANTICS", "Independent binary-tolerance and rounding audit at the published threshold.", ("spin2014.fill.continuous",)),
    _row("spin2014.octree.boundary_subtrees", "Spin-It (2014), Section 5.4", "METHOD", "MATCH_BOUNDARY_SUBTREE_CONSTRUCTION", "Source-to-tree audit proving boundary subtrees are constructed from the original surface."),
    _row("spin2014.octree.split_merge", "Spin-It (2014), Section 5.3", "METHOD", "MATCH_SPLIT_MERGE_REFINEMENT", "Deterministic refinement trace covering both split and merge operations.", ("spin2014.octree.boundary_subtrees",)),
    _row("spin2014.regularizer.laplacian", "Spin-It (2014), Section 5.3", "METHOD", "MATCH_PUBLISHED_LAPLACIAN", "Independent matrix/value audit using the published weight and neighborhood semantics."),
    _row("spin2014.solver.slqp", "Spin-It (2014), Section 5.4", "METHOD", "MATCH_SLQP", "Iteration trace and independent KKT audit for the published SLQP scheme."),
    _row("spin2014.solver.active_set", "Spin-It (2014), Section 5.4", "METHOD", "MATCH_ACTIVE_SET", "Trace demonstrating published active-set transitions and constraint handling.", ("spin2014.solver.slqp",)),
    _row("spin2014.solver.lbfgs", "Spin-It (2014), Section 5.4", "METHOD", "MATCH_LBFGS", "Trace and independent step audit of the published L-BFGS subproblem behavior.", ("spin2014.solver.slqp",)),
    _row("spin2014.cage.complete_octree", "Spin-It (2014), Section 6 and Figure 4", "METHOD", "MATCH_COMPLETE_OCTREE_CAGE", "Geometry audit proving every required intersecting cell contributes to the cage."),
    _row("spin2014.cage.bbw_interpolation", "Spin-It (2014), Section 6 and Figure 4", "METHOD", "MATCH_BBW_INTERPOLATION", "Independent BBW solve and interpolation comparison on the complete cage.", ("spin2014.cage.complete_octree",)),
    _row("spin2014.deformation.joint", "Spin-It (2014), Section 6.1", "METHOD", "MATCH_JOINT_DEFORMATION", "Joint beta/shape optimization trace with independent geometry checks.", ("spin2014.cage.bbw_interpolation",)),
    _row("spin2014.material.dual_density", "Spin-It (2014), Section 7", "METHOD", "MATCH_DUAL_DENSITY", "Independent material-label and density-weighted moment audit."),
    _row("spin2014.supplement.algorithm.1", "Spin-It (2014), Supplement Algorithm 1", "ALGORITHM", "MATCH_PUBLISHED_ALGORITHM", "Stepwise trace checked independently against Supplement Algorithm 1."),
    _row("spin2014.supplement.algorithm.2", "Spin-It (2014), Supplement Algorithm 2", "ALGORITHM", "MATCH_PUBLISHED_ALGORITHM", "Stepwise trace checked independently against Supplement Algorithm 2."),
]


def _case(
    requirement_id: str,
    locator: str,
    outcome: str,
    description: str,
    dependencies: tuple[str, ...] = (),
) -> dict[str, Any]:
    return _row(
        requirement_id,
        locator,
        "ORIGINAL_EXPERIMENT",
        outcome,
        "Checked AUTHOR_ORIGINAL or PUBLISHED_EXACT_SPECIFICATION asset plus explicit critical-parameter evidence; "
        + description,
        dependencies,
    )


_ROWS += [
    _case("spin2014.case.ellipsoid", "Spin-It (2014), Figure 6", "MATCH_ELLIPSOID_TOP", "independent sealed-geometry, moment, solve, and outcome checks."),
    _case("spin2014.case.heart", "Spin-It (2014), Figure 7", "MATCH_HEART_TOP", "independent sealed-geometry, moment, solve, and outcome checks."),
    _case("spin2014.case.armadillo.hand", "Spin-It (2014), Figure 8 hand", "MATCH_ARMADILLO_HAND", "hand pose, cavity, cuts, and outcome checked separately."),
    _case("spin2014.case.armadillo.shell", "Spin-It (2014), Figure 8 shell", "MATCH_ARMADILLO_SHELL", "shell pose, cavity, cuts, and outcome checked separately."),
    _case("spin2014.case.teapot.inertia_only", "Spin-It (2014), Section 8 rotational-stability comparison: inertia-only strategy", "MATCH_TEAPOT_INERTIA_ONLY", "the inertia-only objective and resulting geometry are independently checked."),
    _case("spin2014.case.teapot.lowering_only", "Spin-It (2014), Section 8 rotational-stability comparison: lowering-only strategy", "MATCH_TEAPOT_LOWERING_ONLY", "the lowering-only objective and resulting geometry are independently checked."),
    _case("spin2014.case.teapot.mass_and_lowering", "Spin-It (2014), Section 8 rotational-stability comparison: mass-and-lowering strategy", "MATCH_TEAPOT_MASS_AND_LOWERING", "the mass-and-lowering objective and resulting geometry are independently checked."),
    _case("spin2014.case.cuboid", "Spin-It (2014), Figure 10 top", "MATCH_CUBOID_YOYO", "yo-yo geometry, cavity/cuts, and outcome are independently checked."),
    _case("spin2014.case.woven_ring", "Spin-It (2014), Figure 10 bottom", "MATCH_WOVEN_RING_YOYO", "woven-ring geometry, cavity/cuts, and outcome are independently checked."),
    _case("spin2014.case.dancing_couple.hollowing", "Spin-It (2014), Figure 11 hollowing", "MATCH_HOLLOWING_INSUFFICIENCY", "hollowing is evaluated without deformation and must faithfully retain the observed insufficient dominant-axis alignment; ordinary converged/stable success is not the acceptance criterion."),
    _case("spin2014.case.dancing_couple.joint_deformation", "Spin-It (2014), Figure 11 joint deformation", "MATCH_JOINT_DEFORMATION_IMPROVEMENT", "joint deformation is evaluated separately and must demonstrate the published alignment improvement relative to the insufficient hollow-only result.", ("spin2014.deformation.joint",)),
    _case("spin2014.case.sumone.dual_density", "Spin-It (2014), Figure 5", "MATCH_SUMONE_DUAL_DENSITY", "both material regions and density-weighted moments are independently checked.", ("spin2014.material.dual_density",)),
    _case("spin2014.case.trex.matched_make_it_stand", "Spin-It (2014), Figure 12", "MATCH_TREX_COMPARISON", "Spin-It and Make-It-Stand use identical source inputs and discretization."),
    _case("spin2014.case.elephant", "Spin-It (2014), Table 1 elephant", "MATCH_ELEPHANT", "all reported cavities, cuts, and accessibility outcomes are checked."),
    _case("spin2014.table1.family.c1_cut1", "Spin-It (2014), Table 1: Cuboid, Ellipsoid, Heart, and T-Rex (1 void, 1 cut)", "MATCH_CAVITY_CUT_FAMILY", "Cuboid, Ellipsoid, Heart, and T-Rex each have one void and one cut; each named model is independently checked."),
    _case("spin2014.table1.family.c1_cut3", "Spin-It (2014), Table 1: Armadillo (hand) (1 void, 3 cuts)", "MATCH_CAVITY_CUT_FAMILY", "Armadillo (hand) has one void and three cuts and is independently checked."),
    _case("spin2014.table1.family.c5_cut2", "Spin-It (2014), Table 1: Armadillo (shell) (5 voids, 2 cuts)", "MATCH_CAVITY_CUT_FAMILY", "Armadillo (shell) has five voids and two cuts and is independently checked."),
    _case("spin2014.table1.family.c2_cut1", "Spin-It (2014), Table 1: Teapot (2 voids, 1 cut)", "MATCH_CAVITY_CUT_FAMILY", "Teapot has two voids and one cut and is independently checked."),
    _case("spin2014.table1.family.c1_cut4", "Spin-It (2014), Table 1: Woven Ring (1 void, 4 cuts)", "MATCH_CAVITY_CUT_FAMILY", "Woven Ring has one void and four cuts and is independently checked."),
    _case("spin2014.table1.family.c1_cut5", "Spin-It (2014), Table 1: Dancing Couple (1 void, 5 cuts)", "MATCH_CAVITY_CUT_FAMILY", "Dancing Couple has one void and five cuts and is independently checked."),
    _case("spin2014.table1.family.c3_cut2", "Spin-It (2014), Table 1: Elephant (3 voids, 2 cuts)", "MATCH_CAVITY_CUT_FAMILY", "Elephant has three voids and two cuts and is independently checked."),
]


_ROWS += [
    _row("faster.moments.raw", "Spin-It Faster, Sections 3-4", "FORMULA", "MATCH_RAW_MOMENTS", "Independent raw-moment values and convention audit."),
    _row("faster.moments.coefficient_duality", "Spin-It Faster, Sections 3-4", "THEOREM", "MATCH_COEFFICIENT_MOMENT_DUALITY", "Independent symbolic and numerical duality check."),
    _row("faster.theorem.quadric_reduction", "Spin-It Faster, main reduction theorem", "THEOREM", "MATCH_THEOREM_CONDITIONS", "Independent evaluation of hypotheses, conclusion, and numerical instances."),
    _row("faster.theorem.optimality", "Spin-It Faster, optimality theorem", "THEOREM", "MATCH_THEOREM_CONDITIONS", "Independent optimality and boundary-condition checks."),
    _row("faster.theorem.degenerate_conditions", "Spin-It Faster, degenerate conditions", "THEOREM", "MATCH_DEGENERATE_CONDITIONS", "Explicit audit of every stated degenerate precondition and outcome."),
    _row("faster.eq14.newton", "Spin-It Faster, Equation 14", "ALGORITHM", "MATCH_NEWTON_UPDATE", "Independent residual and Newton-step trace."),
    _row("faster.eq14.geodesic", "Spin-It Faster, Equation 14", "ALGORITHM", "MATCH_GEODESIC_UPDATE", "Independent sphere-geodesic and step-cap audit.", ("faster.eq14.newton",)),
    _row("faster.eq14.multiplier", "Spin-It Faster, Equation 14", "ALGORITHM", "MATCH_MULTIPLIER_UPDATE", "Independent per-step multiplier computation audit.", ("faster.eq14.newton",)),
    _row("faster.eq14.restart", "Spin-It Faster, Equation 14 and implementation notes", "ALGORITHM", "MATCH_RESTART_PROTOCOL", "All initializations, failed attempts, evaluation limits, and restarts retained.", ("faster.eq14.newton",)),
    _row("faster.discrete_moment_map.differentiable", "Spin-It Faster, Section 5", "METHOD", "MATCH_DIFFERENTIABLE_DISCRETE_MAP", "Independent finite-difference audit including intersection-vertex motion."),
    _row("faster.geometry.inverse_gauss_map", "Spin-It Faster, Section 5", "METHOD", "MATCH_INVERSE_GAUSS_MAP", "Independent construction and topology audit from fixed sphere normals."),
    _row("faster.geometry.ideal_points", "Spin-It Faster, Section 5", "METHOD", "MATCH_IDEAL_POINTS", "Independent audit of ideal-point construction and domain handling."),
    _row("faster.geometry.caps", "Spin-It Faster, Section 5", "METHOD", "MATCH_CAPS", "Independent cap construction, closure, and orientation audit."),
    _row("faster.geometry.degenerate_types", "Spin-It Faster, Section 5", "METHOD", "MATCH_DEGENERATE_TYPES", "Each published degenerate geometry type classified and checked independently."),
    _row("faster.geometry.boolean_independent", "Spin-It Faster, Section 5", "INDEPENDENT_CHECK", "MATCH_BOOLEAN_RESULT", "Boolean result checked by an implementation independent of the primary construction."),
    _row("faster.2d.dogleg", "Spin-It Faster, 2D dogleg example", "ORIGINAL_EXPERIMENT", "MATCH_2D_DOGLEG", "Actual 2D moments, topology, and result checked independently."),
    _row("faster.2d.polygon_circle", "Spin-It Faster, polygon/circle example", "ORIGINAL_EXPERIMENT", "MATCH_POLYGON_CIRCLE", "Polygon/disk Boolean and moment results checked independently."),
    _row("faster.2d.extrusion", "Spin-It Faster, 2D-to-3D extrusion comparison", "ORIGINAL_EXPERIMENT", "MATCH_EXTRUSION_EQUIVALENCE", "Independent 2D/3D moment and objective equivalence audit."),
    _row("faster.hydrostatics.actual", "Spin-It Faster, buoyancy examples", "METHOD", "MATCH_ACTUAL_HYDROSTATICS", "Geometry-derived displaced volume, buoyancy center, and stability checked independently."),
    _row("faster.floating_island.chained_assets", "Spin-It Faster, floating-island examples", "ORIGINAL_EXPERIMENT", "MATCH_CHAINED_FLOATING_ISLAND", "Author-original chained A-then-B assets and both optimization stages checked.", ("faster.hydrostatics.actual",)),
    _row("faster.cad.roundtrip", "Spin-It Faster, CAD example", "ORIGINAL_EXPERIMENT", "MATCH_CAD_ROUNDTRIP", "Original CAD input, pinned kernel, B-rep/STEP roundtrip, and independent mass audit."),
]


for case_id, label in (
    ("static_tilted_vase", "static tilted vase"),
    ("floating_tea_infuser", "floating tea infuser"),
    ("floating_island_cat", "floating island cat"),
    ("floating_island_dog_human", "floating island dog/human"),
    ("heart_spinning_top", "heart spinning top"),
    ("noodles_spinning_top", "noodles spinning top"),
    ("jellyfish_yoyo", "jellyfish yo-yo"),
    ("flower_fidget_spinner", "flower fidget spinner"),
):
    _ROWS.append(
        _case(
            f"faster.table1.{case_id}",
            f"Spin-It Faster, Table 1: {label}",
            "MATCH_TABLE1_CASE",
            "50 starts, failures, first-success definition, timings, geometry, residuals, and two full repeats are checked.",
        )
    )


_ROWS += [
    _case("faster.figure.10", "Spin-It Faster, Figure 10", "MATCH_FIGURE_10", "quadratic convergence is estimated independently above numerical noise."),
    _case("faster.figure.11", "Spin-It Faster, Figure 11", "MATCH_FIGURE_11", "all initial weights, refinements, endpoints, and dominated rows are retained."),
    _case("faster.figure.12", "Spin-It Faster, Figure 12", "MATCH_FIGURE_12", "all initial weights, refinements, endpoints, and dominated rows are retained."),
    _case("faster.figure.13", "Spin-It Faster, Figure 13", "MATCH_FIGURE_13", "near-feasibility bracket and nine common-initialization offsets are independently checked."),
    _case("faster.section6_6.example.1", "Spin-It Faster, Section 6.6 first example", "NONUNIQUE_ZERO_SENSITIVITY", "on Omega=(0,1), independently verify objective nonnegativity, volume 1/2 witnesses with objective 0 and p=0, and uncountably many binary and nonbinary optima; ordinary nondegenerate Newton convergence is not the acceptance criterion."),
    _case("faster.section6_6.example.2", "Spin-It Faster, Section 6.6 second example", "FULL_DOMAIN_ZERO_SENSITIVITY", "on Omega=(0,1/2), independently verify the unique optimum chi=1 almost everywhere, volume 1/2, objective 0 and p=0; this is a full-domain boundary solution, not infeasibility or an ordinary nondegenerate quadric interface."),
    _row("shared.sealed_geometry", "Both papers, all experiments", "EVIDENCE_CHAIN", "SEALED_GEOMETRY_MATCHES", "Independent topology, containment, wall-thickness, and final geometry audits."),
    _row("shared.independent_moments", "Both papers, all experiments", "EVIDENCE_CHAIN", "INDEPENDENT_MOMENTS_MATCH", "Independent mass, first moment, second moment, inertia, and density audit."),
    _row("shared.provenance", "Both papers, all experiments", "EVIDENCE_CHAIN", "PROVENANCE_COMPLETE", "Frozen source, build, runner, checker, asset, and parameter provenance."),
    _row("shared.determinism", "Both papers, all experiments", "EVIDENCE_CHAIN", "TWO_RUN_HASH_MATCH", "Two same-environment scientific runs with matching hashes and complete retained failures."),
    _row("shared.registered_tests", "Both papers, verification suite", "EVIDENCE_CHAIN", "REGISTERED_TESTS_RUN", "CTest inventory and execution record with no unexpected skip or Not Run."),
    _row("shared.independent_simulations", "Both papers, validation", "EVIDENCE_CHAIN", "INDEPENDENT_SIMULATIONS_MATCH", "Independent dynamics, fixed-pivot, tilt, and hydrostatic perturbation results."),
    _row("shared.real_desktop_visualization", "Both papers, visualization", "EVIDENCE_CHAIN", "REAL_DESKTOP_VIEWER_AUDITED", "Real desktop viewer trace bound to the same sealed artifact as computation."),
    _row("physical.2014.measurements", "Spin-It (2014), fabricated-object measurements", "PHYSICAL_MEASUREMENT", "NOT_IN_SCOPE", "Explicit declaration that fabricated-object measurements are not part of computational reproduction."),
    _row("physical.faster.measurements", "Spin-It Faster, fabricated-object measurements", "PHYSICAL_MEASUREMENT", "NOT_IN_SCOPE", "Explicit declaration that fabricated-object measurements are not part of computational reproduction."),
]


REQUIREMENTS = tuple(_ROWS)
REQUIRED_REQUIREMENT_IDS = tuple(row["requirement_id"] for row in REQUIREMENTS)


def catalogue() -> dict[str, Any]:
    """Return the immutable catalogue content as JSON-compatible values."""

    return {
        "schema": CATALOGUE_SCHEMA,
        "schema_version": CATALOGUE_VERSION,
        "requirements": [
            {
                **row,
                "dependencies": list(row["dependencies"]),
            }
            for row in REQUIREMENTS
        ],
    }
