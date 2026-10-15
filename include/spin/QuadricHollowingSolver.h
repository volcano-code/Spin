#pragma once

#include "spin/QuadricCoefficients.h"
#include "spin/QuadricMomentEvaluator.h"
#include "spin/SpinDesignRequest.h"
#include "spin/SpinMomentProblem.h"
#include "spin/MassMomentProblem.h"
#include "spin/LegacyPredictionTransaction.h"

#include <Eigen/Core>

#include <functional>
#include <limits>
#include <optional>
#include <string>
#include <vector>


namespace spin
{


struct SpinMomentMapEvaluation
{
    Eigen::Matrix<double, 10, 1> rawMoments =
        Eigen::Matrix<double, 10, 1>::Zero();
    /** Variable-material response before adding fixed-shell moments. */
    std::optional<Eigen::Matrix<double, 10, 1>> variableRawMoments;
    Eigen::Matrix<double, 10, 10> momentJacobian =
        Eigen::Matrix<double, 10, 10>::Zero();
    bool derivativeDefined = false;
    std::string derivativeDiagnostic;
    double residualEvaluationSeconds = 0.0;
    double derivativeEvaluationSeconds = 0.0;
    std::string geometryHash;
    std::string topologyHash;
    /** Per-call identity of the exact material/moment/Jacobian map. Empty for
     * synthetic or external evaluators that do not provide a bound identity. */
    std::string evaluationMapIdentity;
    std::string booleanBackendIdentity;
    std::string derivativeIdentity;
};


using SpinMomentEvaluator = std::function<SpinMomentMapEvaluation(
    const QuadricCoefficients10&
)>;

/** Optional synchronous observer for finalized solver events. It is intended
 * for durable diagnostics; an empty sink preserves the legacy in-memory path. */
using SpinNewtonEventSink = std::function<void(const SpinSolverEvent&)>;


/** Result of a reversible diagonal equilibration solve for one Eq. 14 linear
 * system. The original matrix is never regularized or otherwise modified. */
struct SpinNewtonLinearSolveResult
{
    bool solved = false;
    bool usedEquilibration = false;
    Eigen::Index originalRank = 0;
    Eigen::Index equilibratedRank = -1;
    Eigen::VectorXd solution;
    Eigen::VectorXd rowScales;
    Eigen::VectorXd columnScales;
    double residualInfinityNorm = std::numeric_limits<double>::infinity();
    double backwardError = std::numeric_limits<double>::infinity();
    std::string diagnostic;
};


/** Solve A x = b, retrying a rank-deficient double-precision factorization
 * only after reversible row/column max-norm scaling. The candidate must pass
 * a backward-error check against the original unscaled system. */
SpinNewtonLinearSolveResult solveSpinNewtonLinearSystem(
    const Eigen::MatrixXd& matrix,
    const Eigen::VectorXd& rightHandSide
);

/** Eq.14 linearization on an explicitly supported dual-coefficient subspace.
 * This is an assembly result, not an optimization/convergence certificate.
 * Unknown order: tangent direction, constraint multipliers, multiplier mu.
 */
struct SpinMomentSubspaceLinearization
{
    Eigen::MatrixXd solverTangent;
    Eigen::MatrixXd momentTangentJacobian;
    Eigen::VectorXd multipliers;
    double multiplierMu = 0.0;
    Eigen::VectorXd fullResidual;
    double rawResidual = std::numeric_limits<double>::infinity();
    Eigen::MatrixXd matrix;
    Eigen::VectorXd rightHandSide;
};

/** Consume a genuine directional material derivative without ambient padding.
 * supportBasis: 10 x k orthonormal basis in dual-monomial coordinates.
 * nativeTangentBasis: 10 x (k-1) orthonormal native geometry directions.
 * nativeDirectionalJacobian: 10 x (k-1), of rawMoments at this same quadric.
 * All omitted optimality equations and their local derivatives must vanish;
 * otherwise throws invalid_argument. Only a local supported-space contract:
 * callers must revalidate at every trial and retain the full raw residual.
 */
SpinMomentSubspaceLinearization linearizeSpinMomentSubspace(
    const SpinMomentProblem& problem,
    const QuadricCoefficients10& quadric,
    const Eigen::Matrix<double, 10, 1>& rawMoments,
    const Eigen::MatrixXd& supportBasis,
    const Eigen::MatrixXd& nativeTangentBasis,
    const Eigen::MatrixXd& nativeDirectionalJacobian
);


struct SpinNewtonResult
{
    SpinSolverStatus status = SpinSolverStatus::Indeterminate;
    std::string diagnostic;
    std::optional<QuadricCoefficients10> quadric;
    double normalizedEquationResidual =
        std::numeric_limits<double>::infinity();
    double physicalRelativeResidual =
        std::numeric_limits<double>::infinity();
    double multiplierMu = 0.0;
    Eigen::VectorXd multipliers;
    int functionEvaluations = 0;
    int restartCount = 0;
    double residualEvaluationSeconds = 0.0;
    double derivativeEvaluationSeconds = 0.0;
    double evaluatorWallSeconds = 0.0;
    std::vector<int> restartFunctionEvaluations;
    int eventStreamVersion = 1;
    std::vector<SpinSolverEvent> events;
    std::vector<SpinEvaluationTiming> evaluationTimings;
    std::string eventHash;
    std::vector<SpinIterationTrace> trace;
    /** Accepted states: legacy v1, or profile-bound dual-coordinate v2. */
    std::string traceHash;
    SpinEquationProfile equationProfile = SpinEquationProfile::LegacyNativeNormalizedV1;
    std::optional<double> rawEquationResidual;
};


/** Eq. 14 spherical Newton with an explicit equation/coordinate profile.
 * The evaluator always consumes native geometry coefficients p and returns
 * independent raw moments and their native-p Jacobian, in either profile.
 */
SpinNewtonResult solveSpinMomentSystem(
    const SpinMomentProblem& problem,
    const SpinSolverSettings& settings,
    const SpinMomentEvaluator& evaluator,
    const std::vector<QuadricCoefficients10>& initialGuesses,
    const SpinNewtonEventSink& eventSink = {}
);

/** Same spherical Newton kernel, with an owned general moment objective.
 * The caller still supplies a genuine geometry/moment derivative map; this
 * low-level interface does not admit or seal a material object on its own.
 */
SpinNewtonResult solveMassMomentSystem(
    const MassMomentProblem& problem,
    const SpinSolverSettings& settings,
    const SpinMomentEvaluator& evaluator,
    const std::vector<QuadricCoefficients10>& initialGuesses,
    const SpinNewtonEventSink& eventSink = {}
);

/** Timing-free scientific comparison. The historical capture-mode eventHash
 * includes payload timing and remains unchanged for backwards compatibility. */
std::string spinSolverScientificEventHash(const SpinNewtonResult& result);

struct LegacyTransactionResult
{
    SpinNewtonResult solver;
    LegacyTransactionAttachment attachment;
};

/** Native/testing adapter only. Public request-backed bindings must derive
 * fixedMomentsTarget from owned real material geometry, never caller numbers.
 * An empty predictor shares the original kernel without altering its events. */
LegacyTransactionResult solveSpinMomentSystemLegacyTransactional(
    const SpinMomentProblem& problem,
    const SpinSolverSettings& settings,
    const SpinMomentEvaluator& trueEvaluator,
    const std::vector<QuadricCoefficients10>& initialGuesses,
    const std::string& contextIdentity,
    const Vector10& fixedMomentsTarget,
    const LegacyPredictor& predictor,
    const LegacyTransactionSettings& transactionSettings
);

/** Directional material map. The legacy ambient momentJacobian must remain
 * zero; only the genuine native basis/Jacobian below is consumed. */
struct SpinMomentSubspaceMapEvaluation
{
    SpinMomentMapEvaluation material;
    Eigen::MatrixXd nativeTangentBasis;
    Eigen::MatrixXd nativeDirectionalJacobian;
};
using SpinMomentSubspaceEvaluator = std::function<SpinMomentSubspaceMapEvaluation(
    const QuadricCoefficients10&)>;

struct SpinMomentSubspaceEvaluationRecord
{
    int evaluationIndex = 0;
    std::array<double, 10> quadricCoefficients{};
    std::optional<SpinMomentSubspaceMapEvaluation> evaluation;
    std::optional<SpinMomentSubspaceLinearization> linearization;
    std::string error;
    std::string errorType;
};

/** Versioned attachment to the ordinary Newton events, not a new acceptance
 * system. Scientific hash binds support, all calls and decisions, excluding
 * timing. A solver status alone does not validate final material geometry. */
struct SpinNewtonSubspaceResult
{
    int attachmentVersion = 1;
    Eigen::MatrixXd supportBasis;
    SpinNewtonResult solver;
    std::vector<SpinMomentSubspaceEvaluationRecord> evaluations;
    std::string scientificHash;
};

/** Explicit supported-space entry into the SAME Newton loop. Strict dual/raw
 * profile only. Each trial is freshly evaluated and relinearized; unsupported
 * states fail through the normal backtracking/restart ledger. */
SpinNewtonSubspaceResult solveSpinMomentSubspaceSystem(
    const SpinMomentProblem& problem,
    const SpinSolverSettings& settings,
    const Eigen::MatrixXd& supportBasis,
    const SpinMomentSubspaceEvaluator& evaluator,
    const std::vector<QuadricCoefficients10>& initialGuesses,
    const SpinNewtonEventSink& eventSink = {}
);


/** Platform-stable, mode-aware restart sequence derived from the fixed seed. */
std::vector<QuadricCoefficients10> deterministicInitialQuadrics(
    SpinMode mode,
    const SpinSolverSettings& settings
);


/** Complete request-backed solver using the exact Boolean moment evaluator. */
QuadricHollowingResult solveQuadricHollowing(
    const SpinDesignRequest& request,
    const QuadricMomentEvaluatorSettings& evaluatorSettings =
        QuadricMomentEvaluatorSettings(),
    const std::vector<QuadricCoefficients10>& initialGuesses = {}
);


} // namespace spin
