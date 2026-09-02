#pragma once

#include "spin/QuadricCoefficients.h"
#include "spin/QuadricMomentEvaluator.h"
#include "spin/SpinDesignRequest.h"
#include "spin/SpinMomentProblem.h"

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
    Eigen::Matrix<double, 10, 10> momentJacobian =
        Eigen::Matrix<double, 10, 10>::Zero();
    bool derivativeDefined = false;
    std::string geometryHash;
    std::string topologyHash;
};


using SpinMomentEvaluator = std::function<SpinMomentMapEvaluation(
    const QuadricCoefficients10&
)>;


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
    std::vector<SpinIterationTrace> trace;
    std::string traceHash;
};


/** Paper Eq. 14 spherical Newton solve against an arbitrary moment map. */
SpinNewtonResult solveSpinMomentSystem(
    const SpinMomentProblem& problem,
    const SpinSolverSettings& settings,
    const SpinMomentEvaluator& evaluator,
    const std::vector<QuadricCoefficients10>& initialGuesses
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
