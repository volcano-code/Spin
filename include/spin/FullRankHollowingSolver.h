#pragma once
#include "spin/FasterTheoremConditions.h"
#include "spin/FullRankMomentEvaluator.h"
#include "spin/QuadricHollowingSolver.h"

namespace spin {
/** Immutable request-backed map; Jacobian columns are native p, not dual s. */
SpinMomentEvaluator makeFullRankMomentEvaluator(
    const SpinDesignRequest&, const QuadricMomentEvaluatorSettings& = {});
/** Strict positive native boundary audit, then independent signed-tetrahedron
 * integration of degree 0/1/2 polynomials. Does not use surface-moment code. */
physics::RawMassMoments10 integrateFullRankBoundaryMoments(
    const TriangleMesh&, double density);
physics::RawMassMoments10 integrateFullRankBoundaryMoments(
    const TriangleMesh&, double density, bool exactBinary64SelfIntersection);
physics::RawMassMoments10 integrateFullRankBoundaryMoments(
    const TriangleMesh&, double density, bool exactBinary64SelfIntersection,
    bool exactBinary64MaterialBoundary);
/** Independently audit the represented exact signed material and integrate
 * its density-weighted degree-0/1/2 tetrahedron moments without rounding.
 * Rational approximations are NOT exact analytic quadric roots. A failed
 * display boundary is neither consumed nor promoted by this entry point. */
std::vector<std::string> integrateFullRankBoundaryMoments(
    const ExactBooleanVertexCoordinates&, const Eigen::MatrixXi&, double density);
struct FullRankHollowingResult {
    SpinNewtonResult numerical;
    bool finalGeometryVerified = false;
    std::string verificationDiagnostic;
    std::optional<FullRankMomentEvaluation> finalEvaluation;
    std::optional<physics::RawMassMoments10> independentTotalMoments;
    std::optional<double> independentlyRecomputedRawResidual;
    /** Narrow Eq. (7) audit; never a theorem-completion declaration. */
    std::optional<FasterTheoremOneNumericalAudit> theoremOneNumericalAudit;
    double totalWallSeconds = 0, finalVerificationSeconds = 0;
    std::string inputHash, outputHash;
};
FullRankHollowingResult solveFullRankQuadricHollowing(
    const SpinDesignRequest&, const QuadricMomentEvaluatorSettings& = {},
    const std::vector<QuadricCoefficients10>& initialGuesses = {},
    const SpinNewtonEventSink& eventSink = {});
/** Same full-rank material path and independent verification, for an explicit
 * moment objective. Nonempty problemIdentity binds the caller's formula and
 * units to the scientific input/output hashes. */
FullRankHollowingResult solveFullRankQuadricHollowing(
    const SpinDesignRequest&, const MassMomentProblem&, const std::string& problemIdentity,
    const QuadricMomentEvaluatorSettings& = {},
    const std::vector<QuadricCoefficients10>& initialGuesses = {},
    const SpinNewtonEventSink& eventSink = {});
}
