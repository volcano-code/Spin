#include "deformation/DeformationOptimizer.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>


namespace deformation
{
namespace
{


void validateLineSearchConfig(
    const PaperLineSearchConfig& config
)
{
    if (
        !std::isfinite(config.initialStep) ||
        !std::isfinite(config.reductionRatio) ||
        !std::isfinite(config.minimumStep) ||
        !std::isfinite(config.requiredRelativeDecrease) ||
        !std::isfinite(config.muDecrement) ||
        config.initialStep <= 0.0 ||
        config.minimumStep <= 0.0 ||
        config.minimumStep > config.initialStep ||
        config.reductionRatio <= 0.0 ||
        config.reductionRatio >= 1.0 ||
        config.requiredRelativeDecrease < 0.0 ||
        config.requiredRelativeDecrease >= 1.0 ||
        config.muDecrement <= 0.0 ||
        config.muDecrement > 1.0
    )
    {
        throw std::invalid_argument(
            "paper line-search configuration is invalid"
        );
    }
}


void validateLineSearchState(
    const PaperLineSearchState& state,
    const PaperLineSearchConfig& config
)
{
    if (
        !std::isfinite(state.step) ||
        !std::isfinite(state.mu) ||
        state.step < config.minimumStep ||
        state.step > config.initialStep ||
        state.mu < 0.0 ||
        state.mu > 1.0
    )
    {
        throw std::invalid_argument(
            "paper line-search state is invalid"
        );
    }
}


double checkedEnergy(
    const PaperEnergyCallback& evaluateEnergy,
    const Eigen::VectorXd& parameters,
    double mu
)
{
    const double value = evaluateEnergy(parameters, mu);
    if (!std::isfinite(value) || value < 0.0)
    {
        throw std::domain_error(
            "paper line-search energy must be finite and non-negative"
        );
    }
    return value;
}


} // namespace


PaperLineSearchResult paperGradientDescentTrial(
    const Eigen::VectorXd& parameters,
    const Eigen::VectorXd& gradient,
    const PaperLineSearchState& state,
    const PaperLineSearchConfig& config,
    const PaperEnergyCallback& evaluateEnergy
)
{
    validateLineSearchConfig(config);
    validateLineSearchState(state, config);

    if (
        parameters.size() == 0 ||
        parameters.size() != gradient.size() ||
        !parameters.allFinite() ||
        !gradient.allFinite() ||
        !evaluateEnergy
    )
    {
        throw std::invalid_argument(
            "paper line-search vectors/callback are invalid"
        );
    }

    PaperLineSearchResult result;
    result.parameters = parameters;
    result.attemptedStep = state.step;
    result.muBefore = state.mu;
    result.nextState = state;
    result.energyBefore = checkedEnergy(
        evaluateEnergy,
        parameters,
        state.mu
    );
    result.trialEnergy = std::numeric_limits<double>::infinity();

    const Eigen::VectorXd trialParameters =
        parameters - state.step * gradient;
    bool validTrial = trialParameters.allFinite();

    if (validTrial)
    {
        try
        {
            result.trialEnergy = checkedEnergy(
                evaluateEnergy,
                trialParameters,
                state.mu
            );
        }
        catch (const InfeasibleDeformationState&)
        {
            validTrial = false;
        }
    }

    const double requiredMaximum =
        (1.0 - config.requiredRelativeDecrease) *
        result.energyBefore;

    if (validTrial && result.trialEnergy <= requiredMaximum)
    {
        result.decision = PaperTrialDecision::Accepted;
        result.accepted = true;
        result.rolledBack = false;
        result.parameters = trialParameters;
        result.energyAfter = result.trialEnergy;
        return result;
    }

    /* Every insufficient trial, including a small decrease, is rolled back. */
    result.accepted = false;
    result.rolledBack = true;
    result.energyAfter = result.energyBefore;

    if (state.step > config.minimumStep)
    {
        result.decision = PaperTrialDecision::StepReduced;
        result.nextState.step = std::max(
            config.minimumStep,
            state.step * config.reductionRatio
        );
        return result;
    }

    result.nextState.step = config.initialStep;

    /*
     * Repeated decimal decrements (for example 0.20 - 3*0.05) need not land
     * exactly on 0.05 in binary floating point. Treat tolerance-scale residue
     * as the paper's exact endpoint so termination is deterministic.
     */
    const double muTolerance =
        64.0 * std::numeric_limits<double>::epsilon() *
        std::max({1.0, std::abs(state.mu), config.muDecrement});

    if (state.mu > config.muDecrement + muTolerance)
    {
        result.decision = PaperTrialDecision::MuReduced;
        result.objectiveReweighted = true;
        result.nextState.mu = state.mu - config.muDecrement;
        if (
            std::abs(result.nextState.mu - config.muDecrement) <=
                muTolerance
        )
        {
            result.nextState.mu = config.muDecrement;
        }
        result.energyAfter = checkedEnergy(
            evaluateEnergy,
            parameters,
            result.nextState.mu
        );
        return result;
    }

    result.decision = PaperTrialDecision::MuExhausted;
    result.objectiveReweighted = true;
    result.nextState.mu = 0.0;
    result.energyAfter = checkedEnergy(
        evaluateEnergy,
        parameters,
        0.0
    );
    return result;
}


DeformationOptimizer::DeformationOptimizer(
    const DeformationObjective& objective,
    const DeformationOptimizerConfig& config
)
    :
    objective_(objective),
    config_(config)
{
    validateLineSearchConfig(config_.lineSearch);
    if (
        config_.maximumIterations <= 0 ||
        !std::isfinite(config_.gradientTolerance) ||
        config_.gradientTolerance < 0.0
    )
    {
        throw std::invalid_argument(
            "deformation optimizer configuration is invalid"
        );
    }
}


DeformationOptimizationResult DeformationOptimizer::optimize(
    const Eigen::VectorXd& initialParameters,
    double initialMu
) const
{
    if (!std::isfinite(initialMu) || initialMu < 0.0 || initialMu > 1.0)
    {
        throw std::invalid_argument(
            "initial deformation mu must lie in [0,1]"
        );
    }

    DeformationOptimizationResult result;
    result.parameters = initialParameters;
    result.mu = initialMu;
    result.evaluation = objective_.evaluate(
        result.parameters,
        result.mu
    );

    PaperLineSearchState state;
    state.step = config_.lineSearch.initialStep;
    state.mu = initialMu;

    for (int iteration = 0;
         iteration < config_.maximumIterations;
         ++iteration)
    {
        const double gradientNorm = result.evaluation.gradient.norm();
        if (!std::isfinite(gradientNorm))
        {
            throw std::overflow_error(
                "deformation optimizer gradient norm overflowed"
            );
        }

        if (gradientNorm <= config_.gradientTolerance)
        {
            result.status =
                DeformationOptimizationStatus::GradientConverged;
            return result;
        }

        const PaperLineSearchResult trial = paperGradientDescentTrial(
            result.parameters,
            result.evaluation.gradient,
            state,
            config_.lineSearch,
            [this](const Eigen::VectorXd& parameters, double mu)
            {
                return objective_.evaluate(parameters, mu).totalEnergy;
            }
        );

        DeformationIterationRecord record;
        record.iteration = iteration;
        record.decision = trial.decision;
        record.accepted = trial.accepted;
        record.energyBefore = trial.energyBefore;
        record.trialEnergy = trial.trialEnergy;
        record.energyAfter = trial.energyAfter;
        record.gradientNorm = gradientNorm;
        record.attemptedStep = trial.attemptedStep;
        record.muBefore = trial.muBefore;
        record.objectiveReweighted = trial.objectiveReweighted;
        record.stateAfter = trial.nextState;
        result.iterations.push_back(record);

        result.parameters = trial.parameters;
        state = trial.nextState;
        result.mu = state.mu;
        result.evaluation = objective_.evaluate(
            result.parameters,
            result.mu
        );

        if (trial.decision == PaperTrialDecision::MuExhausted)
        {
            result.status =
                DeformationOptimizationStatus::MuExhausted;
            return result;
        }
    }

    result.status = DeformationOptimizationStatus::MaximumIterations;
    return result;
}


} // namespace deformation
