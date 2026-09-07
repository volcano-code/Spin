#include "spin/QuadricHollowingSolver.h"

#include "manufacturing/SlicerValidation.h"

#include <Eigen/LU>
#include <Eigen/QR>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>


namespace
{


constexpr int momentDimension = 10;


struct NewtonState
{
    spin::SpinMomentMapEvaluation map;
    spin::SpinMomentEvaluation problem;
    Eigen::VectorXd lambda;
    double mu = 0.0;
    Eigen::Matrix<double, momentDimension, 1> optimalityResidual =
        Eigen::Matrix<double, momentDimension, 1>::Zero();
    double normalizedResidual = std::numeric_limits<double>::infinity();
    double physicalResidual = std::numeric_limits<double>::infinity();
};


class InvalidNumericalState : public std::runtime_error
{
public:
    using std::runtime_error::runtime_error;
};


void validateSettings(const spin::SpinSolverSettings& settings)
{
    if (
        settings.maximumIterations <= 0 || settings.maximumRestarts <= 0 ||
        settings.maximumFunctionEvaluationsPerRestart <= 0 ||
        !std::isfinite(settings.equationTolerance) ||
        settings.equationTolerance <= 0.0 ||
        !std::isfinite(settings.physicalRelativeTolerance) ||
        settings.physicalRelativeTolerance <= 0.0 ||
        settings.lineSearchReduction != 0.5 ||
        settings.maximumGeodesicStep != 0.5
    )
    {
        throw std::invalid_argument("invalid spherical Newton settings");
    }
}


NewtonState evaluateState(
    const spin::SpinMomentProblem& problem,
    const spin::QuadricCoefficients10& quadric,
    spin::SpinMomentMapEvaluation map
)
{
    NewtonState state;
    state.map = std::move(map);
    if (
        !state.map.rawMoments.allFinite() ||
        !state.map.momentJacobian.allFinite()
    )
    {
        throw InvalidNumericalState(
            "moment evaluator returned non-finite data"
        );
    }
    const double rawMomentNorm = state.map.rawMoments.norm();
    const double momentJacobianNorm = state.map.momentJacobian.norm();
    if (
        !std::isfinite(rawMomentNorm) ||
        !std::isfinite(momentJacobianNorm)
    )
    {
        throw InvalidNumericalState(
            "moment evaluator data overflowed its numerical scale"
        );
    }
    if (!state.map.derivativeDefined)
    {
        throw std::runtime_error("moment evaluator derivative is undefined" +
            (state.map.derivativeDiagnostic.empty() ? std::string() :
                ": " + state.map.derivativeDiagnostic));
    }
    state.problem = problem.evaluate(state.map.rawMoments);
    const Eigen::Index constraintCount = state.problem.constraints.size();
    if (
        state.problem.constraintJacobian.rows() != constraintCount ||
        state.problem.constraintJacobian.cols() != momentDimension ||
        state.problem.constraintHessians.size() !=
            static_cast<std::size_t>(constraintCount)
    )
    {
        throw std::runtime_error("moment problem derivative shape is invalid");
    }
    bool finiteConstraintHessians = true;
    for (const auto& hessian : state.problem.constraintHessians)
    {
        finiteConstraintHessians =
            finiteConstraintHessians && hessian.allFinite();
    }
    if (
        !std::isfinite(state.problem.objective) ||
        !state.problem.objectiveGradient.allFinite() ||
        !state.problem.objectiveHessian.allFinite() ||
        !state.problem.constraints.allFinite() ||
        !state.problem.constraintJacobian.allFinite() ||
        !finiteConstraintHessians
    )
    {
        throw InvalidNumericalState(
            "moment problem produced a non-finite objective or derivative"
        );
    }

    Eigen::MatrixXd projection(momentDimension, constraintCount + 1);
    if (constraintCount > 0)
    {
        projection.leftCols(constraintCount) =
            -state.problem.constraintJacobian.transpose();
    }
    projection.col(constraintCount) = quadric.values();
    const Eigen::VectorXd coefficients =
        projection.completeOrthogonalDecomposition().solve(
            -state.problem.objectiveGradient
        );
    if (!coefficients.allFinite())
    {
        throw InvalidNumericalState(
            "moment multiplier solve produced non-finite coefficients"
        );
    }
    state.lambda = coefficients.head(constraintCount);
    state.mu = coefficients(constraintCount);
    state.optimalityResidual =
        state.problem.objectiveGradient + state.mu * quadric.values();
    if (constraintCount > 0)
    {
        state.optimalityResidual -=
            state.problem.constraintJacobian.transpose() * state.lambda;
    }

    if (!state.optimalityResidual.allFinite())
    {
        throw InvalidNumericalState(
            "optimality residual produced non-finite coefficients"
        );
    }
    const double optimalitySquared =
        state.optimalityResidual.squaredNorm();
    const double constraintSquared = state.problem.constraints.squaredNorm();
    const double combinedSquared = optimalitySquared + constraintSquared;
    const double objectiveGradientNorm =
        state.problem.objectiveGradient.norm();
    const double constraintJacobianNorm =
        state.problem.constraintJacobian.norm();
    const double equationScale = std::max({
        1.0,
        objectiveGradientNorm,
        constraintJacobianNorm
    });
    const double physicalScale = std::max(1.0, rawMomentNorm);
    if (
        !std::isfinite(optimalitySquared) ||
        !std::isfinite(constraintSquared) ||
        !std::isfinite(combinedSquared) ||
        !std::isfinite(objectiveGradientNorm) ||
        !std::isfinite(constraintJacobianNorm) ||
        !std::isfinite(equationScale) ||
        !std::isfinite(physicalScale)
    )
    {
        throw InvalidNumericalState(
            "Eq. 14 residual normalization overflowed"
        );
    }
    state.normalizedResidual = std::sqrt(combinedSquared) / equationScale;
    state.physicalResidual = std::sqrt(constraintSquared) / physicalScale;
    if (
        !std::isfinite(state.normalizedResidual) ||
        !std::isfinite(state.physicalResidual) || !std::isfinite(state.mu)
    )
    {
        throw InvalidNumericalState("Eq. 14 evaluation overflowed");
    }
    return state;
}


spin::SpinIterationTrace traceEntry(
    int restart,
    int iteration,
    int evaluations,
    double step,
    const NewtonState& state
)
{
    spin::SpinIterationTrace entry;
    entry.restart = restart;
    entry.iteration = iteration;
    entry.functionEvaluations = evaluations;
    entry.normalizedEquationResidual = state.normalizedResidual;
    entry.physicalRelativeResidual = state.physicalResidual;
    entry.objective = state.problem.objective;
    entry.multiplierMu = state.mu;
    entry.acceptedGeodesicStep = step;
    entry.quadricHash.clear();
    entry.materialGeometryHash = state.map.geometryHash;
    entry.materialTopologyHash = state.map.topologyHash;
    return entry;
}


std::string traceHash(const std::vector<spin::SpinIterationTrace>& trace)
{
    std::ostringstream stream;
    stream.imbue(std::locale::classic());
    stream << "SpinNewtonTrace-v1;" << std::hexfloat;
    for (const spin::SpinIterationTrace& entry : trace)
    {
        stream << entry.restart << ';' << entry.iteration << ';'
               << entry.functionEvaluations << ';'
               << entry.normalizedEquationResidual << ';'
               << entry.physicalRelativeResidual << ';'
               << entry.objective << ';' << entry.multiplierMu << ';'
               << entry.acceptedGeodesicStep << ';'
               << entry.quadricHash << ';'
               << entry.materialGeometryHash << ';'
               << entry.materialTopologyHash << ';';
    }
    return manufacturing::computeContentSha256(stream.str());
}


std::array<double, momentDimension> coefficientArray(
    const spin::QuadricCoefficients10& quadric
)
{
    std::array<double, momentDimension> result{};
    for (int index = 0; index < momentDimension; ++index)
    {
        result[static_cast<std::size_t>(index)] = quadric.values()(index);
    }
    return result;
}


std::optional<double> finiteMetric(double value)
{
    return std::isfinite(value) ? std::optional<double>(value) : std::nullopt;
}


spin::SpinSolverEvent restartEvent(
    const char* type,
    int restart,
    int iteration,
    int evaluationIndex,
    int restartEvaluationIndex,
    const char* outcome,
    const char* reason
)
{
    spin::SpinSolverEvent event;
    event.type = type;
    event.restart = restart;
    event.iteration = iteration;
    event.evaluationIndex = evaluationIndex;
    event.restartEvaluationIndex = restartEvaluationIndex;
    event.outcome = outcome;
    event.reason = reason;
    return event;
}


std::string solverEventHash(const std::vector<spin::SpinSolverEvent>& events)
{
    std::ostringstream stream;
    stream.imbue(std::locale::classic());
    stream << "SpinSolverEvents-v1;" << std::hexfloat;
    const auto optionalDouble = [&stream](const std::optional<double>& value) {
        if (value) stream << "some:" << *value;
        else stream << "none";
        stream << ';';
    };
    const auto optionalString = [&stream](
        const std::optional<std::string>& value
    ) {
        if (value) stream << "some:" << value->size() << ':' << *value;
        else stream << "none";
        stream << ';';
    };
    for (const spin::SpinSolverEvent& event : events)
    {
        stream << event.type.size() << ':' << event.type << ';'
               << event.restart << ';' << event.iteration << ';'
               << event.evaluationIndex << ';'
               << event.restartEvaluationIndex << ';';
        optionalDouble(event.step);
        stream << event.outcome.size() << ':' << event.outcome << ';'
               << event.reason.size() << ':' << event.reason << ';'
               << event.detail.size() << ':' << event.detail << ';';
        if (event.quadricCoefficients)
        {
            stream << "some;";
            for (double coefficient : *event.quadricCoefficients)
            {
                stream << coefficient << ';';
            }
        }
        else
        {
            stream << "none;";
        }
        stream << event.quadricHash.size() << ':' << event.quadricHash << ';';
        optionalDouble(event.normalizedEquationResidual);
        optionalDouble(event.physicalRelativeResidual);
        optionalDouble(event.objective);
        optionalDouble(event.multiplierMu);
        optionalString(event.materialGeometryHash);
        optionalString(event.materialTopologyHash);
    }
    return manufacturing::computeContentSha256(stream.str());
}


bool converged(
    const NewtonState& state,
    const spin::SpinSolverSettings& settings
)
{
    return state.normalizedResidual <= settings.equationTolerance &&
        state.physicalResidual <= settings.physicalRelativeTolerance &&
        state.mu > 0.0;
}


std::uint64_t splitMix64(std::uint64_t* state)
{
    *state += 0x9e3779b97f4a7c15ULL;
    std::uint64_t value = *state;
    value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31U);
}


double deterministicSignedUnit(std::uint64_t* state)
{
    const std::uint64_t bits = splitMix64(state) >> 11U;
    const double unit = static_cast<double>(bits) / 9007199254740992.0;
    return 2.0 * unit - 1.0;
}


spin::SpinMomentProblem problemForRequest(
    const spin::SpinDesignRequest& request
)
{
    switch (request.mode)
    {
    case spin::SpinMode::StaticBalance:
        return request.staticObjective == spin::StaticObjective::ParetoTradeoff
            ? spin::SpinMomentProblem::staticBalancePareto(
                request.staticParetoWeight
            )
            : spin::SpinMomentProblem::staticBalance(
                request.staticObjective
            );
    case spin::SpinMode::BuoyantStability:
        return spin::SpinMomentProblem::buoyantStability(
            request.buoyancy.liquidDensity *
                request.buoyancy.displacedLiquidVolume,
            request.buoyancy.centerOfBuoyancyTarget.head<2>()
        );
    case spin::SpinMode::FidgetSpinner:
        return spin::SpinMomentProblem::fidgetSpinner();
    case spin::SpinMode::YoYo:
        return spin::SpinMomentProblem::yoYo();
    case spin::SpinMode::SpinningTop:
        return spin::SpinMomentProblem::spinningTop(
            request.centerOfMassWeight
        );
    }
    throw std::invalid_argument("unsupported Spin-It Faster mode");
}


std::string hollowingOutputHash(
    const spin::QuadricHollowingResult& result
)
{
    std::ostringstream stream;
    stream.imbue(std::locale::classic());
    stream << "QuadricHollowingResult-v2;" << result.inputHash << ';'
           << result.geometryHash << ';' << result.traceHash << ';'
           << result.eventHash << ';'
           << spin::spinSolverStatusToken(result.status) << ';'
           << std::hexfloat << result.normalizedEquationResidual << ';'
           << result.physicalRelativeResidual << ';';
    if (std::isfinite(result.buoyantMetacenterHeight))
    {
        stream << result.buoyantMetacenterHeight << ';'
               << result.buoyantStabilityMargin << ';';
    }
    else
    {
        stream << "not-buoyant;";
    }
    for (int index = 0; index < result.quadricCoefficients.size(); ++index)
    {
        stream << result.quadricCoefficients(index) << ';';
    }
    return manufacturing::computeContentSha256(stream.str());
}


} // namespace


namespace spin
{


SpinNewtonResult solveSpinMomentSystem(
    const SpinMomentProblem& problem,
    const SpinSolverSettings& settings,
    const SpinMomentEvaluator& evaluator,
    const std::vector<QuadricCoefficients10>& initialGuesses
)
{
    validateSettings(settings);
    if (!evaluator)
    {
        throw std::invalid_argument("spherical Newton evaluator is empty");
    }
    if (initialGuesses.empty())
    {
        throw std::invalid_argument(
            "spherical Newton requires deterministic initial guesses"
        );
    }

    SpinNewtonResult result;
    bool sawGeometryFailure = false;
    const auto finish = [&result]() {
        result.traceHash = ::traceHash(result.trace);
        result.eventHash = ::solverEventHash(result.events);
    };
    const int restartLimit = std::min(
        settings.maximumRestarts,
        static_cast<int>(initialGuesses.size())
    );
    for (int restart = 0; restart < restartLimit; ++restart)
    {
        result.restartCount = restart + 1;
        result.diagnostic.clear();
        QuadricCoefficients10 current = initialGuesses[restart];
        int evaluationsThisRestart = 0;
        result.restartFunctionEvaluations.push_back(0);
        SpinSolverEvent start = restartEvent(
            "attempt_start", restart, 0, result.functionEvaluations, 0,
            "started", "restart_started"
        );
        start.quadricCoefficients = coefficientArray(current);
        start.quadricHash = current.hash();
        result.events.push_back(std::move(start));

        struct EvaluationResult
        {
            std::optional<NewtonState> state;
            std::string error;
            bool invalidNumericalState = false;
        };
        const auto evaluateAttempt = [&](const QuadricCoefficients10& quadric,
                                         int iteration,
                                         std::optional<double> step) {
            ++evaluationsThisRestart;
            ++result.functionEvaluations;
            result.restartFunctionEvaluations.back() =
                evaluationsThisRestart;
            SpinSolverEvent event = restartEvent(
                "evaluation", restart, iteration,
                result.functionEvaluations, evaluationsThisRestart,
                "failed", "evaluation_exception"
            );
            event.step = step;
            event.quadricCoefficients = coefficientArray(quadric);
            event.quadricHash = quadric.hash();
            SpinEvaluationTiming timing;
            timing.evaluationIndex = result.functionEvaluations;
            timing.restart = restart;
            timing.restartEvaluationIndex = evaluationsThisRestart;
            const auto started = std::chrono::steady_clock::now();
            const auto recordFailure = [&](const std::exception& error,
                                           const char* reason,
                                           bool invalidNumericalState) {
                const auto failed = std::chrono::steady_clock::now();
                if (!timing.evaluatorReturned)
                {
                    timing.wallSeconds = std::chrono::duration<double>(
                        failed - started
                    ).count();
                    result.evaluatorWallSeconds += timing.wallSeconds;
                }
                event.reason = reason;
                event.detail = error.what();
                result.evaluationTimings.push_back(std::move(timing));
                result.events.push_back(std::move(event));
                return EvaluationResult{
                    std::nullopt,
                    error.what(),
                    invalidNumericalState
                };
            };
            try
            {
                SpinMomentMapEvaluation map = evaluator(quadric);
                const auto returned = std::chrono::steady_clock::now();
                timing.evaluatorReturned = true;
                timing.wallSeconds =
                    std::chrono::duration<double>(returned - started).count();
                timing.reportedResidualSeconds =
                    map.residualEvaluationSeconds;
                timing.reportedDerivativeSeconds =
                    map.derivativeEvaluationSeconds;
                result.evaluatorWallSeconds += timing.wallSeconds;
                result.residualEvaluationSeconds +=
                    map.residualEvaluationSeconds;
                result.derivativeEvaluationSeconds +=
                    map.derivativeEvaluationSeconds;
                if (!map.geometryHash.empty())
                {
                    event.materialGeometryHash = map.geometryHash;
                }
                if (!map.topologyHash.empty())
                {
                    event.materialTopologyHash = map.topologyHash;
                }
                NewtonState state = evaluateState(
                    problem, quadric, std::move(map)
                );
                event.outcome = "evaluated";
                event.reason = "evaluation_succeeded";
                event.normalizedEquationResidual = finiteMetric(
                    state.normalizedResidual
                );
                event.physicalRelativeResidual = finiteMetric(
                    state.physicalResidual
                );
                event.objective = finiteMetric(state.problem.objective);
                event.multiplierMu = finiteMetric(state.mu);
                result.evaluationTimings.push_back(std::move(timing));
                result.events.push_back(std::move(event));
                return EvaluationResult{std::move(state), "", false};
            }
            catch (const InvalidNumericalState& error)
            {
                return recordFailure(
                    error, "invalid_numerical_state", true
                );
            }
            catch (const std::exception& error)
            {
                return recordFailure(error, "evaluation_exception", false);
            }
            catch (...)
            {
                const auto failed = std::chrono::steady_clock::now();
                timing.wallSeconds = std::chrono::duration<double>(
                    failed - started
                ).count();
                result.evaluatorWallSeconds += timing.wallSeconds;
                result.evaluationTimings.push_back(std::move(timing));
                throw;
            }
        };
        const auto terminate = [&](int iteration, const char* reason) {
            result.events.push_back(restartEvent(
                "attempt_terminal", restart, iteration,
                result.functionEvaluations, evaluationsThisRestart,
                "terminated", reason
            ));
        };

        EvaluationResult initial = evaluateAttempt(current, 0, 0.0);
        if (!initial.state)
        {
            sawGeometryFailure = true;
            result.diagnostic = initial.error;
            terminate(
                0,
                initial.invalidNumericalState
                    ? "initial_numerical_state_invalid"
                    : "initial_evaluation_failed"
            );
            continue;
        }
        NewtonState state = std::move(*initial.state);
        result.events.back().outcome = "accepted";
        result.events.back().reason = "initial_evaluation_succeeded";
        result.trace.push_back(
            traceEntry(restart, 0, evaluationsThisRestart, 0.0, state)
        );
        result.trace.back().quadricHash = current.hash();
        if (state.map.momentJacobian.norm() <= 1e-15)
        {
            result.status = SpinSolverStatus::OptimizationDidNotConverge;
            result.diagnostic = "initial moment Jacobian vanishes";
            terminate(0, "initial_sensitivity_vanished");
            continue;
        }
        if (converged(state, settings))
        {
            result.status = SpinSolverStatus::Converged;
            result.quadric = current;
            result.normalizedEquationResidual = state.normalizedResidual;
            result.physicalRelativeResidual = state.physicalResidual;
            result.multiplierMu = state.mu;
            result.multipliers = state.lambda;
            result.diagnostic.clear();
            terminate(0, "converged");
            finish();
            return result;
        }

        for (int iteration = 1; iteration <= settings.maximumIterations;
             ++iteration)
        {
            const Eigen::Index constraintCount =
                state.problem.constraints.size();
            Eigen::Matrix<double, momentDimension, momentDimension> hessianL =
                state.problem.objectiveHessian;
            for (Eigen::Index constraint = 0;
                 constraint < constraintCount;
                 ++constraint)
            {
                hessianL -= state.lambda(constraint) *
                    state.problem.constraintHessians[
                        static_cast<std::size_t>(constraint)
                    ];
            }
            if (!hessianL.allFinite())
            {
                result.status = SpinSolverStatus::OptimizationDidNotConverge;
                result.diagnostic =
                    "Eq. 14 Lagrangian Hessian overflowed";
                terminate(iteration, "newton_system_invalid");
                break;
            }

            const Eigen::Index systemSize =
                momentDimension + constraintCount + 1;
            Eigen::MatrixXd system =
                Eigen::MatrixXd::Zero(systemSize, systemSize);
            Eigen::VectorXd rightHandSide = Eigen::VectorXd::Zero(systemSize);
            system.topLeftCorner(momentDimension, momentDimension) =
                hessianL * state.map.momentJacobian +
                state.mu *
                    Eigen::Matrix<double, momentDimension, momentDimension>::Identity();
            if (constraintCount > 0)
            {
                system.block(
                    0,
                    momentDimension,
                    momentDimension,
                    constraintCount
                ) = -state.problem.constraintJacobian.transpose();
                system.block(
                    momentDimension,
                    0,
                    constraintCount,
                    momentDimension
                ) = state.problem.constraintJacobian *
                    state.map.momentJacobian;
            }
            system.block(
                0,
                momentDimension + constraintCount,
                momentDimension,
                1
            ) = current.values();
            system.block(
                momentDimension + constraintCount,
                0,
                1,
                momentDimension
            ) = 2.0 * current.values().transpose();
            rightHandSide.head(momentDimension) =
                -state.optimalityResidual;
            if (constraintCount > 0)
            {
                rightHandSide.segment(momentDimension, constraintCount) =
                    -state.problem.constraints;
            }
            if (!system.allFinite() || !rightHandSide.allFinite())
            {
                result.status = SpinSolverStatus::OptimizationDidNotConverge;
                result.diagnostic = "Eq. 14 bordered Newton system overflowed";
                terminate(iteration, "newton_system_invalid");
                break;
            }

            Eigen::FullPivLU<Eigen::MatrixXd> factorization(system);
            if (factorization.rank() < systemSize)
            {
                result.status = SpinSolverStatus::OptimizationDidNotConverge;
                result.diagnostic = "bordered Newton system is rank deficient";
                terminate(iteration, "bordered_system_rank_deficient");
                break;
            }
            Eigen::Matrix<double, momentDimension, 1> direction =
                factorization.solve(rightHandSide).head(momentDimension);
            direction = current.projectToTangent(direction);
            const double directionNorm = direction.norm();
            if (!direction.allFinite() || !std::isfinite(directionNorm))
            {
                result.status = SpinSolverStatus::OptimizationDidNotConverge;
                result.diagnostic = "Newton direction is non-finite";
                terminate(iteration, "newton_direction_invalid");
                break;
            }
            if (directionNorm <= 1e-15)
            {
                result.status = SpinSolverStatus::OptimizationDidNotConverge;
                result.diagnostic = "Newton direction vanished before convergence";
                terminate(iteration, "newton_direction_vanished");
                break;
            }

            double step = std::min(
                directionNorm,
                settings.maximumGeodesicStep
            );
            bool accepted = false;
            NewtonState candidateState;
            QuadricCoefficients10 candidate = current;
            while (
                evaluationsThisRestart <
                    settings.maximumFunctionEvaluationsPerRestart &&
                step > std::numeric_limits<double>::epsilon()
            )
            {
                candidate = current.geodesicUpdate(direction, step);
                EvaluationResult trial = evaluateAttempt(
                    candidate, iteration, step
                );
                if (trial.state)
                {
                    candidateState = std::move(*trial.state);
                    if (candidateState.normalizedResidual <
                        state.normalizedResidual)
                    {
                        result.events.back().outcome = "accepted";
                        result.events.back().reason = "residual_decreased";
                        accepted = true;
                        break;
                    }
                    result.events.back().outcome = "rejected";
                    result.events.back().reason = "no_residual_decrease";
                }
                else
                {
                    sawGeometryFailure = true;
                    result.diagnostic = trial.error;
                }
                step *= settings.lineSearchReduction;
            }
            if (!accepted)
            {
                result.status =
                    evaluationsThisRestart >=
                        settings.maximumFunctionEvaluationsPerRestart
                    ? SpinSolverStatus::ResourceLimitReached
                    : SpinSolverStatus::LineSearchExhausted;
                if (result.diagnostic.empty())
                {
                    result.diagnostic = "Eq. 14 line search found no decrease";
                }
                terminate(
                    iteration,
                    evaluationsThisRestart >=
                            settings.maximumFunctionEvaluationsPerRestart
                        ? "maximum_evaluations_reached"
                        : "line_search_exhausted"
                );
                break;
            }

            current = candidate;
            state = candidateState;
            result.trace.push_back(traceEntry(
                restart,
                iteration,
                evaluationsThisRestart,
                step,
                state
            ));
            result.trace.back().quadricHash = current.hash();
            if (converged(state, settings))
            {
                result.status = SpinSolverStatus::Converged;
                result.quadric = current;
                result.normalizedEquationResidual = state.normalizedResidual;
                result.physicalRelativeResidual = state.physicalResidual;
                result.multiplierMu = state.mu;
                result.multipliers = state.lambda;
                result.diagnostic.clear();
                terminate(iteration, "converged");
                finish();
                return result;
            }
            if (iteration == settings.maximumIterations)
            {
                result.status = SpinSolverStatus::ResourceLimitReached;
                result.diagnostic = "maximum Newton iterations reached";
                terminate(iteration, "maximum_iterations_reached");
            }
        }
    }

    if (result.status == SpinSolverStatus::Indeterminate)
    {
        result.status = sawGeometryFailure
            ? SpinSolverStatus::GeometryEvaluationFailed
            : SpinSolverStatus::OptimizationDidNotConverge;
    }
    if (!result.trace.empty())
    {
        const SpinIterationTrace& final = result.trace.back();
        result.normalizedEquationResidual =
            final.normalizedEquationResidual;
        result.physicalRelativeResidual = final.physicalRelativeResidual;
        result.multiplierMu = final.multiplierMu;
    }
    finish();
    return result;
}


std::vector<QuadricCoefficients10> deterministicInitialQuadrics(
    SpinMode mode,
    const SpinSolverSettings& settings
)
{
    validateSettings(settings);
    std::vector<QuadricCoefficients10> guesses;
    guesses.reserve(static_cast<std::size_t>(settings.maximumRestarts));
    std::uint64_t state =
        static_cast<std::uint64_t>(settings.deterministicSeed) ^
        0x5350494e49544641ULL;
    for (int restart = 0; restart < settings.maximumRestarts; ++restart)
    {
        QuadricCoefficients10::Vector values;
        for (int coordinate = 0; coordinate < values.size(); ++coordinate)
        {
            values(coordinate) = deterministicSignedUnit(&state);
        }
        switch (mode)
        {
        case SpinMode::StaticBalance:
            values(0) = std::abs(values(0)) + 0.25;
            break;
        case SpinMode::BuoyantStability:
            values.tail<6>().setZero();
            values(0) = std::abs(values(0)) + 0.25;
            break;
        case SpinMode::FidgetSpinner:
            values(3) = 0.0;
            values(6) = 0.0;
            values(7) = 0.0;
            values(8) = 0.0;
            break;
        case SpinMode::YoYo:
            values(0) = 0.0;
            values(6) = -std::abs(values(6)) - 0.25;
            break;
        case SpinMode::SpinningTop:
            values(6) = -std::abs(values(6)) - 0.1;
            break;
        }
        guesses.push_back(
            QuadricCoefficients10::fromUnnormalized(values)
        );
    }
    return guesses;
}


QuadricHollowingResult solveQuadricHollowing(
    const SpinDesignRequest& request,
    const QuadricMomentEvaluatorSettings& evaluatorSettings,
    const std::vector<QuadricCoefficients10>& initialGuesses
)
{
    validateSpinDesignRequest(request);
    const std::vector<QuadricCoefficients10> guesses =
        initialGuesses.empty()
        ? deterministicInitialQuadrics(request.mode, request.settings)
        : initialGuesses;
    std::string lastDerivativeDiagnostic;
    std::string lastDerivativeScientificHash;
    const SpinMomentEvaluator evaluator =
        [&](const QuadricCoefficients10& quadric) {
            const QuadricMomentEvaluation evaluation =
                evaluateQuadricMaterialMoments(
                    request,
                    quadric,
                    evaluatorSettings
                );
            SpinMomentMapEvaluation map;
            lastDerivativeDiagnostic = evaluation.derivativeDiagnostic;
            lastDerivativeScientificHash = evaluation.derivativeScientificHash;
            map.rawMoments = evaluation.totalMomentsTarget.values();
            map.momentJacobian = evaluation.momentJacobian;
            map.derivativeDefined = evaluation.derivativeDefined;
            map.derivativeDiagnostic = evaluation.derivativeDiagnostic;
            map.residualEvaluationSeconds =
                evaluation.residualEvaluationSeconds;
            map.derivativeEvaluationSeconds =
                evaluation.derivativeEvaluationSeconds;
            map.geometryHash = evaluation.geometryHash;
            map.topologyHash = evaluation.topologyHash;
            return map;
        };
    const SpinNewtonResult solved = solveSpinMomentSystem(
        problemForRequest(request),
        request.settings,
        evaluator,
        guesses
    );

    QuadricHollowingResult result;
    result.derivativeMethod = quadricMomentDerivativeMethodToken(evaluatorSettings.derivativeMethod);
    result.derivativeDiagnostic = lastDerivativeDiagnostic;
    result.derivativeScientificHash = lastDerivativeScientificHash;
    result.status = solved.status;
    result.diagnostic = solved.diagnostic;
    result.inputHash = computeSpinDesignRequestHash(request);
    result.outerSurface = request.outerSurface;
    result.normalizedEquationResidual = solved.normalizedEquationResidual;
    result.physicalRelativeResidual = solved.physicalRelativeResidual;
    result.trace = solved.trace;
    result.traceHash = solved.traceHash;
    result.solveFunctionEvaluations = solved.functionEvaluations;
    result.solveRestartFunctionEvaluations =
        solved.restartFunctionEvaluations;
    result.solveEvaluatorWallSeconds = solved.evaluatorWallSeconds;
    result.eventStreamVersion = solved.eventStreamVersion;
    result.events = solved.events;
    result.evaluationTimings = solved.evaluationTimings;
    result.eventHash = solved.eventHash;
    result.residualEvaluationSeconds = solved.residualEvaluationSeconds;
    result.derivativeEvaluationSeconds = solved.derivativeEvaluationSeconds;
    if (solved.quadric.has_value())
    {
        const auto evaluateForMaterialExport =
            [&](const QuadricCoefficients10& quadric) {
                ++result.postSolveMaterialEvaluationCount;
                const auto started = std::chrono::steady_clock::now();
                try
                {
                    QuadricMomentEvaluation evaluation =
                        evaluateQuadricMaterialMoments(
                            request, quadric, evaluatorSettings
                        );
                    result.postSolveMaterialEvaluationWallSeconds +=
                        std::chrono::duration<double>(
                            std::chrono::steady_clock::now() - started
                        ).count();
                    result.postSolveMaterialResidualEvaluationSeconds +=
                        evaluation.residualEvaluationSeconds;
                    result.postSolveMaterialDerivativeEvaluationSeconds +=
                        evaluation.derivativeEvaluationSeconds;
                    return evaluation;
                }
                catch (...)
                {
                    result.postSolveMaterialEvaluationWallSeconds +=
                        std::chrono::duration<double>(
                            std::chrono::steady_clock::now() - started
                        ).count();
                    throw;
                }
            };
        result.quadricCoefficients = solved.quadric->values();
        try
        {
            const QuadricMomentEvaluation finalEvaluation =
                evaluateForMaterialExport(*solved.quadric);
            result.finalMoments = finalEvaluation.totalMomentsTarget;
            result.derivativeDiagnostic = finalEvaluation.derivativeDiagnostic;
            result.derivativeScientificHash = finalEvaluation.derivativeScientificHash;
            if (request.mode == SpinMode::BuoyantStability)
            {
                result.buoyantMetacenterHeight =
                    buoyantMetacenterHeight(request.buoyancy);
                result.buoyantStabilityMargin = buoyantStabilityMargin(
                    result.finalMoments.values(),
                    request.buoyancy
                );
                if (result.buoyantStabilityMargin <= 0.0)
                {
                    result.status = SpinSolverStatus::Infeasible;
                    result.diagnostic =
                        "equilibrium center of mass is not below the "
                        "prescribed-waterline metacenter";
                }
            }
            result.materialBoundary = finalEvaluation.materialBoundary;
            const QuadricCoefficients10 complement =
                QuadricCoefficients10::fromUnnormalized(
                    -solved.quadric->values()
                );
            const QuadricMomentEvaluation cavityEvaluation =
                evaluateForMaterialExport(complement);
            result.cavityBoundary = cavityEvaluation.materialBoundary;
            result.geometryHash = manufacturing::computeContentSha256(
                std::string("QuadricHollowingGeometry-v2;") +
                result.inputHash + ';' + finalEvaluation.geometryHash + ';' +
                cavityEvaluation.geometryHash
            );
        }
        catch (const std::exception& error)
        {
            result.status = SpinSolverStatus::InvalidMaterialArtifact;
            result.diagnostic = error.what();
        }
    }
    result.outputHash = hollowingOutputHash(result);
    return result;
}


} // namespace spin
