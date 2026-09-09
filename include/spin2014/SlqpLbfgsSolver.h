#pragma once

#include "spin2014/AdaptiveFillFunctional.h"
#include "spin2014/SlqpKkt.h"

#include <Eigen/Core>

#include <string>
#include <functional>
#include <vector>


namespace spin2014
{

// Differentiable equalities with a [0,1] bounded prefix; remaining variables
// are unbounded. Callback identity is not a code or input attestation.
struct EqualityBoxProblem
{
    Eigen::Index variableCount = 0;
    Eigen::Index betaCount = 0;
    std::function<AdaptiveFillEvaluation(const Eigen::VectorXd&)> evaluate;
};


enum class SlqpStatus
{
    Converged,
    MaximumIterations,
    LineSearchExhausted,
    InvalidProblem,
    RestorationStalled,
    NumericalUnknown
};

enum class SlqpAlgorithm { LegacyV3, AuditedKktRestorationV1 };
const char* slqpAlgorithmToken(SlqpAlgorithm) noexcept;


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
    SlqpAlgorithm algorithm = SlqpAlgorithm::LegacyV3;
    SlqpAuditSettings audit;
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
    KktResidual kkt;
    Eigen::VectorXd state, gradient, constraints;
    Eigen::MatrixXd constraintJacobian;
    // 0 free, 1 lower, 2 upper; unbounded coordinates always 0.
    std::vector<unsigned char> activeMembership;
    std::string membershipHash;
    std::string historyResetReason;
    int historyBefore = 0;
    int historyAfter = 0;
    bool curvaturePairAccepted = false;
    bool curvaturePairRejected = false;
    double linearizedConstraintResidual = 0.;
    RestorationEvidence restoration;
    std::vector<RestorationTrial> objectiveTrials;
    bool finalState = false;
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
    SlqpAlgorithm algorithm = SlqpAlgorithm::LegacyV3;
    KktResidual finalKkt;
    Eigen::VectorXd finalState;
    Eigen::Index boundedVariableCount = 0;
};


const char* slqpStatusToken(SlqpStatus status) noexcept;

SlqpLbfgsResult solveSlqpLbfgs(const EqualityBoxProblem& problem,
    const Eigen::VectorXd& initial, const SlqpLbfgsSettings& settings = {});

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
