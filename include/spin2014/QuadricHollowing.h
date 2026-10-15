#pragma once
#include "spin2014/QuadricMomentProblem.h"
#include "spin/CanonicalFullRankHollowing.h"

namespace spin2014
{
/** Copies geometry and computes the fixed shell from outer minus design in
 * MODEL coordinates at its own density. Exact inclusive containment is a
 * prerequisite. A nonzero supplied fixed moment vector must agree; never use
 * it as authoritative input. No material manufacture/seal authority here.
 */
spin::SpinDesignRequest prepareSpin2014QuadricRequest(
    const spin::SpinDesignRequest&,const QuadricObjectiveSettings&,double fixedMaterialDensity);

/** Versioned algorithm attachment around existing request/result types.
 * FullRankHollowingResult is the ordinary numerical/material audit; this
 * adapter does not create a second completion/acceptance system. */
struct QuadricMaterialRun
{
    std::string objectiveProfile, objectiveIdentity;
    QuadricObjectiveSettings objectiveSettings;
    spin::SpinDesignRequest measuredRequest;
    double fixedMaterialDensity = 0;
    spin::FullRankHollowingResult solve;
    std::optional<spin::FullRankMomentEvaluation> cavityEvaluation;
    /** Availability of both exports, not a partition/manufacturing pass. */
    bool materialExportComplete = false;
    std::string materialExportDiagnostic;
};

/** Explicit rational geometry export: DESIGN minus the SAME represented
 * retained material. Independently checks moments, containment and partition;
 * does not optimize, certify KKT, qualify floating display or grant a seal.
 * Returned cavity has no quadric Jacobian: its source is the retained tape. */
spin::FullRankMomentEvaluation exportSpin2014RationalCavity(
    const spin::SpinDesignRequest&,const spin::FullRankMomentEvaluation& retained);

QuadricMaterialRun solveSpin2014QuadricHollowing(
    const spin::SpinDesignRequest&,const QuadricObjectiveSettings&,double fixedMaterialDensity,
    const spin::QuadricMomentEvaluatorSettings&,
    const std::vector<spin::QuadricCoefficients10>& initialGuesses = {},
    const spin::SpinNewtonEventSink& eventSink = {});
} // namespace spin2014
