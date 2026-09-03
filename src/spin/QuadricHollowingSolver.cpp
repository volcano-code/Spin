#include "spin/QuadricHollowingSolver.h"

#include "manufacturing/SlicerValidation.h"

#include <Eigen/LU>
#include <Eigen/QR>

#include <algorithm>
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
    const spin::SpinMomentEvaluator& evaluator,
    const spin::QuadricCoefficients10& quadric
)
{
    NewtonState state;
    state.map = evaluator(quadric);
    if (
        !state.map.rawMoments.allFinite() ||
        !state.map.momentJacobian.allFinite()
    )
    {
        throw std::runtime_error("moment evaluator returned non-finite data");
    }
    if (!state.map.derivativeDefined)
    {
        throw std::runtime_error("moment evaluator derivative is undefined");
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
    state.lambda = coefficients.head(constraintCount);
    state.mu = coefficients(constraintCount);
    state.optimalityResidual =
        state.problem.objectiveGradient + state.mu * quadric.values();
    if (constraintCount > 0)
    {
        state.optimalityResidual -=
            state.problem.constraintJacobian.transpose() * state.lambda;
    }

    const double combinedSquared =
        state.optimalityResidual.squaredNorm() +
        state.problem.constraints.squaredNorm();
    const double equationScale = std::max({
        1.0,
        state.problem.objectiveGradient.norm(),
        state.problem.constraintJacobian.norm()
    });
    const double physicalScale = std::max(1.0, state.map.rawMoments.norm());
    state.normalizedResidual = std::sqrt(combinedSquared) / equationScale;
    state.physicalResidual =
        state.problem.constraints.norm() / physicalScale;
    if (
        !std::isfinite(state.normalizedResidual) ||
        !std::isfinite(state.physicalResidual) || !std::isfinite(state.mu)
    )
    {
        throw std::runtime_error("Eq. 14 evaluation overflowed");
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
        return spin::SpinMomentProblem::staticBalance(
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
    stream << "QuadricHollowingResult-v1;" << result.inputHash << ';'
           << result.geometryHash << ';' << result.traceHash << ';'
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
    const int restartLimit = std::min(
        settings.maximumRestarts,
        static_cast<int>(initialGuesses.size())
    );
    for (int restart = 0; restart < restartLimit; ++restart)
    {
        result.restartCount = restart + 1;
        QuadricCoefficients10 current = initialGuesses[restart];
        int evaluationsThisRestart = 0;
        NewtonState state;
        try
        {
            state = evaluateState(problem, evaluator, current);
            ++evaluationsThisRestart;
            ++result.functionEvaluations;
        }
        catch (const std::exception& error)
        {
            sawGeometryFailure = true;
            result.diagnostic = error.what();
            continue;
        }
        result.trace.push_back(
            traceEntry(restart, 0, evaluationsThisRestart, 0.0, state)
        );
        result.trace.back().quadricHash = current.hash();
        if (state.map.momentJacobian.norm() <= 1e-15)
        {
            result.status = SpinSolverStatus::OptimizationDidNotConverge;
            result.diagnostic = "initial moment Jacobian vanishes";
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
            result.traceHash = ::traceHash(result.trace);
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

            Eigen::FullPivLU<Eigen::MatrixXd> factorization(system);
            if (factorization.rank() < systemSize)
            {
                result.status = SpinSolverStatus::OptimizationDidNotConverge;
                result.diagnostic = "bordered Newton system is rank deficient";
                break;
            }
            Eigen::Matrix<double, momentDimension, 1> direction =
                factorization.solve(rightHandSide).head(momentDimension);
            direction = current.projectToTangent(direction);
            const double directionNorm = direction.norm();
            if (!direction.allFinite() || directionNorm <= 1e-15)
            {
                result.status = SpinSolverStatus::OptimizationDidNotConverge;
                result.diagnostic = "Newton direction vanished before convergence";
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
                try
                {
                    candidateState = evaluateState(
                        problem,
                        evaluator,
                        candidate
                    );
                    ++evaluationsThisRestart;
                    ++result.functionEvaluations;
                    if (candidateState.normalizedResidual <
                        state.normalizedResidual)
                    {
                        accepted = true;
                        break;
                    }
                }
                catch (const std::exception& error)
                {
                    sawGeometryFailure = true;
                    result.diagnostic = error.what();
                    ++evaluationsThisRestart;
                    ++result.functionEvaluations;
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
                result.traceHash = ::traceHash(result.trace);
                return result;
            }
            if (iteration == settings.maximumIterations)
            {
                result.status = SpinSolverStatus::ResourceLimitReached;
                result.diagnostic = "maximum Newton iterations reached";
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
    result.traceHash = ::traceHash(result.trace);
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
    const SpinMomentEvaluator evaluator =
        [&](const QuadricCoefficients10& quadric) {
            const QuadricMomentEvaluation evaluation =
                evaluateQuadricMaterialMoments(
                    request,
                    quadric,
                    evaluatorSettings
                );
            SpinMomentMapEvaluation map;
            map.rawMoments = evaluation.totalMomentsTarget.values();
            map.momentJacobian = evaluation.momentJacobian;
            map.derivativeDefined = evaluation.derivativeDefined;
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
    result.status = solved.status;
    result.diagnostic = solved.diagnostic;
    result.inputHash = computeSpinDesignRequestHash(request);
    result.outerSurface = request.outerSurface;
    result.normalizedEquationResidual = solved.normalizedEquationResidual;
    result.physicalRelativeResidual = solved.physicalRelativeResidual;
    result.trace = solved.trace;
    result.traceHash = solved.traceHash;
    if (solved.quadric.has_value())
    {
        result.quadricCoefficients = solved.quadric->values();
        try
        {
            const QuadricMomentEvaluation finalEvaluation =
                evaluateQuadricMaterialMoments(
                    request,
                    *solved.quadric,
                    evaluatorSettings
                );
            result.finalMoments = finalEvaluation.totalMomentsTarget;
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
                evaluateQuadricMaterialMoments(
                    request,
                    complement,
                    evaluatorSettings
                );
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
