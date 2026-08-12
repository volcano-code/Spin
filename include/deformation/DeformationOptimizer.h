#pragma once

#include "deformation/DeformationObjective.h"

#include <Eigen/Core>

#include <functional>
#include <vector>


namespace deformation
{


struct PaperLineSearchConfig
{
    double initialStep = 1.0;
    double reductionRatio = 0.8;
    double minimumStep = 0.4;
    double requiredRelativeDecrease = 0.01;
    double muDecrement = 0.05;
};


struct PaperLineSearchState
{
    double step = 1.0;
    double mu = 0.75;
};


enum class PaperTrialDecision
{
    Accepted,
    StepReduced,
    MuReduced,
    MuExhausted
};


struct PaperLineSearchResult
{
    PaperTrialDecision decision = PaperTrialDecision::StepReduced;
    bool accepted = false;
    bool rolledBack = true;
    Eigen::VectorXd parameters;
    double energyBefore = 0.0;
    double trialEnergy = 0.0;
    double energyAfter = 0.0;
    double attemptedStep = 0.0;
    double muBefore = 0.0;
    /** True when energyAfter uses nextState.mu rather than muBefore. */
    bool objectiveReweighted = false;
    PaperLineSearchState nextState;
};


using PaperEnergyCallback = std::function<double(
    const Eigen::VectorXd&,
    double
)>;


/** One transactional paper §5.3 gradient trial. */
PaperLineSearchResult paperGradientDescentTrial(
    const Eigen::VectorXd& parameters,
    const Eigen::VectorXd& gradient,
    const PaperLineSearchState& state,
    const PaperLineSearchConfig& config,
    const PaperEnergyCallback& evaluateEnergy
);


struct DeformationOptimizerConfig
{
    PaperLineSearchConfig lineSearch;
    int maximumIterations = 200;
    double gradientTolerance = 1.0e-8;
};


enum class DeformationOptimizationStatus
{
    GradientConverged,
    MuExhausted,
    MaximumIterations
};


struct DeformationIterationRecord
{
    int iteration = 0;
    PaperTrialDecision decision = PaperTrialDecision::StepReduced;
    bool accepted = false;
    double energyBefore = 0.0;
    double trialEnergy = 0.0;
    double energyAfter = 0.0;
    double gradientNorm = 0.0;
    double attemptedStep = 0.0;
    double muBefore = 0.0;
    bool objectiveReweighted = false;
    PaperLineSearchState stateAfter;
};


struct DeformationOptimizationResult
{
    DeformationOptimizationStatus status =
        DeformationOptimizationStatus::MaximumIterations;
    Eigen::VectorXd parameters;
    double mu = 0.0;
    DeformationEnergyEvaluation evaluation;
    std::vector<DeformationIterationRecord> iterations;
};


class DeformationOptimizer
{
public:
    DeformationOptimizer(
        const DeformationObjective& objective,
        const DeformationOptimizerConfig& config =
            DeformationOptimizerConfig()
    );

    DeformationOptimizationResult optimize(
        const Eigen::VectorXd& initialParameters,
        double initialMu
    ) const;

private:
    const DeformationObjective& objective_;
    DeformationOptimizerConfig config_;
};


} // namespace deformation
