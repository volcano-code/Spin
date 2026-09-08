#pragma once

#include "spin2014/AdaptiveFillFunctional.h"

#include <Eigen/Core>

#include <string>
#include <vector>


namespace spin2014
{


enum class SlqpStatus
{
    Converged,
    MaximumIterations,
    LineSearchExhausted,
    InvalidProblem
};


struct SlqpLbfgsSettings
{
    int maximumIterations = 120;
    int maximumLineSearchSteps = 24;
    int lbfgsMemory = 8;
    double equalityTolerance = 1e-8;
    double projectedGradientTolerance = 1e-7;
    double activeSetTolerance = 1e-10;
    double meritPenalty = 100.0;
    double lineSearchReduction = 0.5;
    double armijo = 1e-4;
};


struct SlqpIteration
{
    int iteration = 0;
    double objective = 0.0;
    double equalityInfinityNorm = 0.0;
    double projectedGradientInfinityNorm = 0.0;
    double kktInfinityNorm = 0.0;
    double meritValue = 0.0;
    double acceptedStep = 0.0;
    int lineSearchEvaluations = 0;
    int activeVariables = 0;
    int lbfgsHistorySize = 0;
    bool usedFirstOrderFallback = false;
    double phi = 0.0;
};


struct SlqpLbfgsResult
{
    SlqpStatus status = SlqpStatus::InvalidProblem;
    std::string diagnostic;
    Eigen::VectorXd beta;
    AdaptiveFillEvaluation evaluation;
    std::vector<SlqpIteration> trace;
    std::string traceHash;
    Spin2014Formulation formulation = Spin2014Formulation::LegacyAngleEliminated;
    double initialPhi = 0.0;
    double finalPhi = 0.0;
};


const char* slqpStatusToken(SlqpStatus status) noexcept;

SlqpLbfgsResult solveSlqpLbfgs(
    const AdaptiveFillFunctional& functional,
    const Eigen::VectorXd& initialBeta,
    const SlqpLbfgsSettings& settings = SlqpLbfgsSettings()
);

// Joint [beta,phi], box bounds apply only to beta. Static requires initialPhi=0.
SlqpLbfgsResult solveSlqpLbfgsPaper(
    const AdaptiveFillFunctional& functional,
    const Eigen::VectorXd& initialBeta,
    double initialPhi,
    const SlqpLbfgsSettings& settings = SlqpLbfgsSettings()
);


} // namespace spin2014
