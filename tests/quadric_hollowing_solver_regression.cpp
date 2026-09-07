#include "spin/QuadricHollowingSolver.h"

#include <Eigen/Core>

#include <cmath>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>


namespace
{


using Vector10 = spin::QuadricCoefficients10::Vector;


int fail(const std::string& message)
{
    std::cerr << "[FAIL] scenario=quadric-hollowing-solver "
              << message << '\n';
    return EXIT_FAILURE;
}


spin::SpinMomentMapEvaluation identityMoments(
    const spin::QuadricCoefficients10& quadric
)
{
    spin::SpinMomentMapEvaluation result;
    result.rawMoments = quadric.values();
    result.momentJacobian = Eigen::Matrix<double, 10, 10>::Identity();
    result.derivativeDefined = true;
    result.geometryHash = quadric.hash();
    result.topologyHash = "synthetic-identity";
    return result;
}


spin::SpinMomentMapEvaluation overflowedFidgetObjective()
{
    spin::SpinMomentMapEvaluation value;
    value.rawMoments(0) = 1.0;
    value.rawMoments(4) = 1e308;
    value.rawMoments(5) = 1e308;
    value.momentJacobian.setIdentity();
    value.derivativeDefined = true;
    value.geometryHash = "overflowed-fidget-geometry";
    value.topologyHash = "overflowed-fidget-topology";
    return value;
}


spin::SpinMomentMapEvaluation finiteFidgetSolutionMoments()
{
    spin::SpinMomentMapEvaluation value;
    value.rawMoments(0) = 1.0;
    value.rawMoments(4) = 1.0;
    value.rawMoments(5) = 1.0;
    value.momentJacobian.setIdentity();
    value.derivativeDefined = true;
    return value;
}


spin::SpinMomentMapEvaluation overflowedResidualNormalization()
{
    spin::SpinMomentMapEvaluation value;
    value.rawMoments(0) = 1.0;
    value.rawMoments(1) = 1.0;
    value.rawMoments(4) = 1e308;
    value.rawMoments(5) = 1e308;
    value.momentJacobian.setIdentity();
    value.derivativeDefined = true;
    value.geometryHash = "overflowed-normalization-geometry";
    value.topologyHash = "overflowed-normalization-topology";
    return value;
}


spin::QuadricCoefficients10 initialPoint()
{
    Vector10 values;
    values << 1.0, 0.12, -0.08, 0.09, -0.07,
        0.05, -0.04, 0.03, -0.02, 0.01;
    return spin::QuadricCoefficients10::fromUnnormalized(values);
}


spin::QuadricCoefficients10 analyticSolution()
{
    Vector10 values = Vector10::Zero();
    values(0) = 1.0;
    return spin::QuadricCoefficients10::fromUnnormalized(values);
}


bool sameScientificEvent(
    const spin::SpinSolverEvent& left,
    const spin::SpinSolverEvent& right
)
{
    return left.type == right.type && left.restart == right.restart &&
        left.iteration == right.iteration &&
        left.evaluationIndex == right.evaluationIndex &&
        left.restartEvaluationIndex == right.restartEvaluationIndex &&
        left.step == right.step && left.outcome == right.outcome &&
        left.reason == right.reason && left.detail == right.detail &&
        left.quadricCoefficients == right.quadricCoefficients &&
        left.quadricHash == right.quadricHash &&
        left.normalizedEquationResidual == right.normalizedEquationResidual &&
        left.physicalRelativeResidual == right.physicalRelativeResidual &&
        left.objective == right.objective &&
        left.multiplierMu == right.multiplierMu &&
        left.materialGeometryHash == right.materialGeometryHash &&
        left.materialTopologyHash == right.materialTopologyHash;
}


bool allPresentMetricsAreFinite(const spin::SpinSolverEvent& event)
{
    return (!event.step || std::isfinite(*event.step)) &&
        (!event.normalizedEquationResidual ||
         std::isfinite(*event.normalizedEquationResidual)) &&
        (!event.physicalRelativeResidual ||
         std::isfinite(*event.physicalRelativeResidual)) &&
        (!event.objective || std::isfinite(*event.objective)) &&
        (!event.multiplierMu || std::isfinite(*event.multiplierMu));
}


} // namespace


int main()
{
    spin::SpinSolverSettings settings;
    settings.maximumIterations = 60;
    settings.maximumRestarts = 1;
    settings.maximumFunctionEvaluationsPerRestart = 60;
    const spin::SpinMomentProblem problem =
        spin::SpinMomentProblem::staticBalance(
            spin::StaticObjective::MaximumMass
        );
    const std::vector<spin::QuadricCoefficients10> deterministic =
        spin::deterministicInitialQuadrics(
            spin::SpinMode::YoYo,
            settings
        );
    const std::vector<spin::QuadricCoefficients10> deterministicRepeat =
        spin::deterministicInitialQuadrics(
            spin::SpinMode::YoYo,
            settings
        );
    if (
        deterministic.size() != 1 ||
        deterministic.front().hash() != deterministicRepeat.front().hash() ||
        deterministic.front().values()(0) != 0.0 ||
        deterministic.front().values()(6) >= 0.0
    )
    {
        return fail("mode-aware deterministic restart sequence is invalid");
    }
    const std::vector<spin::QuadricCoefficients10> buoyantInitial =
        spin::deterministicInitialQuadrics(
            spin::SpinMode::BuoyantStability,
            settings
        );
    if (
        buoyantInitial.size() != 1 ||
        buoyantInitial.front().values().tail<6>().norm() != 0.0
    )
    {
        return fail("buoyant restart does not use the published planar family");
    }
    const spin::SpinNewtonResult solved = spin::solveSpinMomentSystem(
        problem,
        settings,
        identityMoments,
        {initialPoint()}
    );
    if (
        solved.status != spin::SpinSolverStatus::Converged ||
        !solved.quadric.has_value() ||
        solved.normalizedEquationResidual > 1e-12 ||
        solved.physicalRelativeResidual > 1e-12 ||
        solved.multiplierMu <= 0.0 || solved.trace.empty() ||
        solved.traceHash.size() != 64
    )
    {
        return fail("identity moment system did not converge");
    }
    Vector10 expected = Vector10::Zero();
    expected(0) = 1.0;
    if (!solved.quadric->values().isApprox(expected, 1e-11))
    {
        return fail("Newton solution is not the analytic coefficient vector");
    }
    const spin::SpinNewtonResult repeated = spin::solveSpinMomentSystem(
        problem,
        settings,
        identityMoments,
        {initialPoint()}
    );
    if (
        repeated.traceHash != solved.traceHash ||
        repeated.quadric->hash() != solved.quadric->hash()
    )
    {
        return fail("Newton solve is not deterministic");
    }

    if (
        solved.eventStreamVersion != 1 || solved.events.empty() ||
        solved.eventHash.size() != 64 ||
        solved.functionEvaluations != static_cast<int>(solved.evaluationTimings.size()) ||
        solved.restartFunctionEvaluations != std::vector<int>{solved.functionEvaluations}
    )
    {
        return fail("complete event stream/count contract is missing");
    }
    if (
        repeated.events.size() != solved.events.size() ||
        repeated.eventHash != solved.eventHash
    )
    {
        return fail("scientific event stream is not deterministic");
    }
    for (std::size_t index = 0; index < solved.events.size(); ++index)
    {
        if (!sameScientificEvent(solved.events[index], repeated.events[index]))
        {
            return fail("scientific event content changed with timing");
        }
    }
    const spin::SpinNewtonResult differentlyTimed =
        spin::solveSpinMomentSystem(
            problem,
            settings,
            [](const spin::QuadricCoefficients10& quadric) {
                spin::SpinMomentMapEvaluation value = identityMoments(quadric);
                value.residualEvaluationSeconds = 123.0;
                value.derivativeEvaluationSeconds = 456.0;
                return value;
            },
            {initialPoint()}
        );
    if (
        differentlyTimed.eventHash != solved.eventHash ||
        differentlyTimed.events.size() != solved.events.size() ||
        differentlyTimed.residualEvaluationSeconds ==
            solved.residualEvaluationSeconds
    )
    {
        return fail("timings leaked into deterministic scientific events");
    }
    for (std::size_t index = 1; index < solved.trace.size(); ++index)
    {
        if (solved.trace[index].normalizedEquationResidual >=
            solved.trace[index - 1].normalizedEquationResidual)
        {
            return fail("accepted residual trace is not strictly decreasing");
        }
    }

    spin::SpinMomentEvaluator undefined =
        [](const spin::QuadricCoefficients10& quadric) {
            spin::SpinMomentMapEvaluation result = identityMoments(quadric);
            result.derivativeDefined = false;
            return result;
        };
    const spin::SpinNewtonResult geometryFailure =
        spin::solveSpinMomentSystem(
            problem,
            settings,
            undefined,
            {initialPoint()}
        );
    if (geometryFailure.status !=
        spin::SpinSolverStatus::GeometryEvaluationFailed)
    {
        return fail("undefined geometry derivative was not fail-closed");
    }


    int restartCalls = 0;
    spin::SpinMomentEvaluator failThenSucceed =
        [&restartCalls](const spin::QuadricCoefficients10& quadric) {
            ++restartCalls;
            if (restartCalls == 1)
            {
                throw std::runtime_error("initial \\\"bad\\\"\\npath");
            }
            spin::SpinMomentMapEvaluation value = identityMoments(quadric);
            value.residualEvaluationSeconds = 0.02 * restartCalls;
            value.derivativeEvaluationSeconds = 0.01 * restartCalls;
            return value;
        };
    spin::SpinSolverSettings restartSettings = settings;
    restartSettings.maximumRestarts = 2;
    const spin::SpinNewtonResult recovered = spin::solveSpinMomentSystem(
        problem,
        restartSettings,
        failThenSucceed,
        {initialPoint(), analyticSolution()}
    );
    if (
        recovered.status != spin::SpinSolverStatus::Converged ||
        recovered.diagnostic != "" || restartCalls != 2 ||
        recovered.functionEvaluations != 2 ||
        recovered.restartFunctionEvaluations != std::vector<int>({1, 1}) ||
        recovered.events.size() != 6 ||
        recovered.events[1].type != "evaluation" ||
        recovered.events[1].outcome != "failed" ||
        recovered.events[1].reason != "evaluation_exception" ||
        recovered.events[1].quadricHash != initialPoint().hash() ||
        !recovered.events[1].quadricCoefficients.has_value() ||
        recovered.events[1].normalizedEquationResidual.has_value() ||
        recovered.events[1].materialGeometryHash.has_value() ||
        recovered.events[2].reason != "initial_evaluation_failed" ||
        recovered.events.back().reason != "converged" ||
        recovered.evaluationTimings.size() != 2 ||
        recovered.evaluationTimings[0].evaluatorReturned ||
        recovered.evaluationTimings[0].wallSeconds < 0.0 ||
        recovered.evaluationTimings[0].reportedResidualSeconds.has_value() ||
        !recovered.evaluationTimings[1].reportedResidualSeconds.has_value()
    )
    {
        return fail("initial failure/later restart evidence is incomplete");
    }

    int allFailureCalls = 0;
    spin::SpinMomentEvaluator alwaysFails =
        [&allFailureCalls](const spin::QuadricCoefficients10&) -> spin::SpinMomentMapEvaluation {
            ++allFailureCalls;
            throw std::runtime_error("always fails");
        };
    const spin::SpinNewtonResult allFailed = spin::solveSpinMomentSystem(
        problem,
        restartSettings,
        alwaysFails,
        {initialPoint(), analyticSolution()}
    );
    if (
        allFailed.status != spin::SpinSolverStatus::GeometryEvaluationFailed ||
        allFailureCalls != 2 || allFailed.functionEvaluations != 2 ||
        allFailed.restartFunctionEvaluations != std::vector<int>({1, 1}) ||
        allFailed.events.size() != 6 ||
        allFailed.events[2].reason != "initial_evaluation_failed" ||
        allFailed.events[5].reason != "initial_evaluation_failed"
    )
    {
        return fail("all-initial-failure accounting is incomplete");
    }

    spin::SpinMomentEvaluator zeroJacobian =
        [](const spin::QuadricCoefficients10& quadric) {
            spin::SpinMomentMapEvaluation result = identityMoments(quadric);
            result.momentJacobian.setZero();
            return result;
        };
    const spin::SpinNewtonResult singular = spin::solveSpinMomentSystem(
        problem,
        settings,
        zeroJacobian,
        {initialPoint()}
    );
    if (singular.status != spin::SpinSolverStatus::OptimizationDidNotConverge)
    {
        return fail("zero moment Jacobian was not rejected");
    }
    if (
        singular.functionEvaluations != 1 || singular.events.size() != 3 ||
        singular.events.back().reason != "initial_sensitivity_vanished"
    )
    {
        return fail("vanished sensitivity terminal evidence is incomplete");
    }

    Vector10 rankDeficientValues = Vector10::Zero();
    rankDeficientValues(1) = 1.0;
    const spin::SpinNewtonResult rankDeficient =
        spin::solveSpinMomentSystem(
            problem,
            settings,
            [](const spin::QuadricCoefficients10& quadric) {
                spin::SpinMomentMapEvaluation value = identityMoments(quadric);
                value.rawMoments = analyticSolution().values();
                return value;
            },
            {spin::QuadricCoefficients10::fromUnnormalized(
                rankDeficientValues
            )}
        );
    if (
        rankDeficient.status !=
            spin::SpinSolverStatus::OptimizationDidNotConverge ||
        rankDeficient.events.back().reason !=
            "bordered_system_rank_deficient"
    )
    {
        return fail("rank-deficient terminal evidence is incomplete");
    }

    int exceptionTrialCalls = 0;
    spin::SpinMomentEvaluator exceptionTrial =
        [&exceptionTrialCalls](const spin::QuadricCoefficients10& quadric) {
            ++exceptionTrialCalls;
            if (exceptionTrialCalls == 2)
            {
                throw std::runtime_error("trial failed");
            }
            return identityMoments(quadric);
        };
    const spin::SpinNewtonResult exceptionRecovered =
        spin::solveSpinMomentSystem(
            problem,
            settings,
            exceptionTrial,
            {initialPoint()}
        );
    bool foundFailedTrial = false;
    for (const spin::SpinSolverEvent& event : exceptionRecovered.events)
    {
        foundFailedTrial = foundFailedTrial ||
            (event.type == "evaluation" && event.iteration == 1 &&
             event.outcome == "failed" && event.step.has_value() &&
             event.quadricCoefficients.has_value() &&
             !event.normalizedEquationResidual.has_value());
    }
    if (
        exceptionRecovered.status != spin::SpinSolverStatus::Converged ||
        !exceptionRecovered.diagnostic.empty() || !foundFailedTrial ||
        exceptionRecovered.functionEvaluations != exceptionTrialCalls
    )
    {
        return fail("exception line-search trial was not retained");
    }

    int rejectedTrialCalls = 0;
    spin::SpinMomentEvaluator rejectedFiniteTrial =
        [&rejectedTrialCalls](const spin::QuadricCoefficients10& quadric) {
            ++rejectedTrialCalls;
            spin::SpinMomentMapEvaluation value = identityMoments(quadric);
            if (rejectedTrialCalls == 2)
            {
                value.rawMoments(1) = 100.0;
                value.rawMoments(2) = 100.0;
                value.geometryHash = "finite-rejected-geometry";
                value.topologyHash = "finite-rejected-topology";
            }
            return value;
        };
    const spin::SpinNewtonResult rejectedRecovered =
        spin::solveSpinMomentSystem(
            problem,
            settings,
            rejectedFiniteTrial,
            {initialPoint()}
        );
    bool foundRejectedFiniteTrial = false;
    for (const spin::SpinSolverEvent& event : rejectedRecovered.events)
    {
        foundRejectedFiniteTrial = foundRejectedFiniteTrial ||
            (event.type == "evaluation" && event.outcome == "rejected" &&
             event.reason == "no_residual_decrease" &&
             event.normalizedEquationResidual.has_value() &&
             event.materialGeometryHash == "finite-rejected-geometry" &&
             event.materialTopologyHash == "finite-rejected-topology");
    }
    if (
        rejectedRecovered.status != spin::SpinSolverStatus::Converged ||
        !foundRejectedFiniteTrial ||
        rejectedRecovered.functionEvaluations != rejectedTrialCalls
    )
    {
        return fail("finite rejected trial metrics/hashes were not retained");
    }

    bool nonStandardPropagated = false;
    try
    {
        spin::solveSpinMomentSystem(
            problem,
            settings,
            [](const spin::QuadricCoefficients10&) -> spin::SpinMomentMapEvaluation {
                throw 17;
            },
            {initialPoint()}
        );
    }
    catch (int value)
    {
        nonStandardPropagated = value == 17;
    }
    if (!nonStandardPropagated)
    {
        return fail("non-standard evaluator exception was silently converted");
    }

    Vector10 overflowQuadricValues = Vector10::Zero();
    overflowQuadricValues(4) = 1.0;
    overflowQuadricValues(5) = 1.0;
    const spin::QuadricCoefficients10 fidgetSolution =
        spin::QuadricCoefficients10::fromUnnormalized(
            overflowQuadricValues
        );
    const spin::SpinNewtonResult overflowObjective =
        spin::solveSpinMomentSystem(
            spin::SpinMomentProblem::fidgetSpinner(),
            settings,
            [](const spin::QuadricCoefficients10&) {
                return overflowedFidgetObjective();
            },
            {fidgetSolution}
        );
    bool sawRejectedOverflowEvaluation = false;
    bool numericalStateContractFailed = false;
    for (const spin::SpinSolverEvent& event : overflowObjective.events)
    {
        if (!allPresentMetricsAreFinite(event))
        {
            return fail("scientific event contains a non-finite metric");
        }
        if (event.type == "evaluation")
        {
            sawRejectedOverflowEvaluation =
                event.outcome == "failed" &&
                event.reason == "invalid_numerical_state" &&
                !event.normalizedEquationResidual.has_value() &&
                !event.physicalRelativeResidual.has_value() &&
                !event.objective.has_value() &&
                !event.multiplierMu.has_value();
        }
    }
    if (
        overflowObjective.status !=
            spin::SpinSolverStatus::GeometryEvaluationFailed ||
        overflowObjective.quadric.has_value() ||
        overflowObjective.functionEvaluations != 1 ||
        overflowObjective.trace.size() != 0 ||
        !sawRejectedOverflowEvaluation ||
        overflowObjective.events.back().reason !=
            "initial_numerical_state_invalid"
    )
    {
        std::cerr << "[DETAIL] overflowed initial objective was accepted\n";
        numericalStateContractFailed = true;
    }
    const spin::SpinNewtonResult overflowObjectiveRepeated =
        spin::solveSpinMomentSystem(
            spin::SpinMomentProblem::fidgetSpinner(),
            settings,
            [](const spin::QuadricCoefficients10&) {
                return overflowedFidgetObjective();
            },
            {fidgetSolution}
        );
    if (
        overflowObjectiveRepeated.eventHash != overflowObjective.eventHash ||
        overflowObjectiveRepeated.events.size() !=
            overflowObjective.events.size()
    )
    {
        std::cerr << "[DETAIL] invalid-state event hash is not deterministic\n";
        numericalStateContractFailed = true;
    }
    else
    {
        for (std::size_t index = 0; index < overflowObjective.events.size();
             ++index)
        {
            if (!sameScientificEvent(
                    overflowObjective.events[index],
                    overflowObjectiveRepeated.events[index]))
            {
                std::cerr << "[DETAIL] invalid-state event content changed\n";
                numericalStateContractFailed = true;
                break;
            }
        }
    }

    int overflowRestartCalls = 0;
    const spin::SpinNewtonResult overflowRestartRecovered =
        spin::solveSpinMomentSystem(
            spin::SpinMomentProblem::fidgetSpinner(),
            restartSettings,
            [&overflowRestartCalls](const spin::QuadricCoefficients10&) {
                ++overflowRestartCalls;
                return overflowRestartCalls == 1
                    ? overflowedFidgetObjective()
                    : finiteFidgetSolutionMoments();
            },
            {fidgetSolution, fidgetSolution}
        );
    if (
        overflowRestartRecovered.status != spin::SpinSolverStatus::Converged ||
        overflowRestartCalls != 2 ||
        overflowRestartRecovered.functionEvaluations != 2 ||
        overflowRestartRecovered.restartFunctionEvaluations !=
            std::vector<int>({1, 1}) ||
        overflowRestartRecovered.events.size() != 6 ||
        overflowRestartRecovered.events[1].outcome != "failed" ||
        overflowRestartRecovered.events[1].reason !=
            "invalid_numerical_state" ||
        overflowRestartRecovered.events[2].reason !=
            "initial_numerical_state_invalid" ||
        overflowRestartRecovered.events.back().reason != "converged"
    )
    {
        std::cerr << "[DETAIL] invalid initial state bypassed restart recovery"
                  << " status=" << spin::spinSolverStatusToken(
                         overflowRestartRecovered.status)
                  << " calls=" << overflowRestartCalls
                  << " events=" << overflowRestartRecovered.events.size()
                  << " last=" << overflowRestartRecovered.events.back().reason
                  << '\n';
        numericalStateContractFailed = true;
    }

    int overflowTrialCalls = 0;
    const spin::SpinNewtonResult overflowTrialRecovered =
        spin::solveSpinMomentSystem(
            problem,
            settings,
            [&overflowTrialCalls](const spin::QuadricCoefficients10& quadric) {
                ++overflowTrialCalls;
                return overflowTrialCalls == 2
                    ? overflowedResidualNormalization()
                    : identityMoments(quadric);
            },
            {initialPoint()}
        );
    bool sawRejectedOverflowTrial = false;
    for (const spin::SpinSolverEvent& event : overflowTrialRecovered.events)
    {
        if (!allPresentMetricsAreFinite(event))
        {
            std::cerr << "[DETAIL] recovery event contains non-finite metric\n";
            numericalStateContractFailed = true;
        }
        sawRejectedOverflowTrial = sawRejectedOverflowTrial ||
            (event.type == "evaluation" && event.iteration == 1 &&
             event.outcome == "failed" &&
             event.reason == "invalid_numerical_state" &&
             event.step.has_value() &&
             !event.normalizedEquationResidual.has_value() &&
             !event.physicalRelativeResidual.has_value() &&
             !event.objective.has_value() &&
             !event.multiplierMu.has_value() &&
             event.materialGeometryHash ==
                 "overflowed-normalization-geometry" &&
             event.materialTopologyHash ==
                 "overflowed-normalization-topology");
    }
    if (
        overflowTrialRecovered.status != spin::SpinSolverStatus::Converged ||
        overflowTrialCalls <= 2 || !sawRejectedOverflowTrial ||
        overflowTrialRecovered.functionEvaluations != overflowTrialCalls ||
        overflowTrialRecovered.restartFunctionEvaluations !=
            std::vector<int>({overflowTrialCalls}) ||
        overflowTrialRecovered.events.back().reason != "converged"
    )
    {
        std::cerr << "[DETAIL] line-search recovery status="
                  << spin::spinSolverStatusToken(overflowTrialRecovered.status)
                  << " calls=" << overflowTrialCalls
                  << " events=" << overflowTrialRecovered.events.size()
                  << " last=" << overflowTrialRecovered.events.back().reason
                  << '\n';
        for (const spin::SpinSolverEvent& event : overflowTrialRecovered.events)
        {
            if (event.type == "evaluation" && event.iteration == 1)
            {
                std::cerr << "[DETAIL] line-search numerical state outcome="
                          << event.outcome << " reason=" << event.reason
                          << " objective_present="
                          << event.objective.has_value() << '\n';
            }
        }
        numericalStateContractFailed = true;
    }

    const auto rejectedInitialNumericalState = [&](
        const spin::SpinMomentProblem& numericalProblem,
        const spin::SpinMomentEvaluator& numericalEvaluator,
        const char* scenario
    ) {
        const spin::SpinNewtonResult numericalResult =
            spin::solveSpinMomentSystem(
                numericalProblem,
                settings,
                numericalEvaluator,
                {initialPoint()}
            );
        const bool rejected =
            numericalResult.status ==
                spin::SpinSolverStatus::GeometryEvaluationFailed &&
            numericalResult.functionEvaluations == 1 &&
            numericalResult.trace.empty() &&
            numericalResult.events.size() == 3 &&
            numericalResult.events[1].outcome == "failed" &&
            numericalResult.events[1].reason == "invalid_numerical_state" &&
            numericalResult.events[2].reason ==
                "initial_numerical_state_invalid";
        if (!rejected)
        {
            std::cerr << "[DETAIL] invalid " << scenario
                      << " state was not rejected at evaluation\n";
        }
        return rejected;
    };

    numericalStateContractFailed =
        !rejectedInitialNumericalState(
            problem,
            [](const spin::QuadricCoefficients10&) {
                return overflowedResidualNormalization();
            },
            "residual normalization"
        ) || numericalStateContractFailed;
    numericalStateContractFailed =
        !rejectedInitialNumericalState(
            problem,
            [](const spin::QuadricCoefficients10& quadric) {
                spin::SpinMomentMapEvaluation value = identityMoments(quadric);
                value.momentJacobian(0, 0) = 1e308;
                value.momentJacobian(0, 1) = 1e308;
                return value;
            },
            "moment-Jacobian norm"
        ) || numericalStateContractFailed;
    numericalStateContractFailed =
        !rejectedInitialNumericalState(
            spin::SpinMomentProblem::staticBalance(
                spin::StaticObjective::LowerCenterOfMass
            ),
            [](const spin::QuadricCoefficients10&) {
                spin::SpinMomentMapEvaluation value;
                value.rawMoments(0) = 1e-309;
                value.rawMoments(3) = 1.0;
                value.momentJacobian.setIdentity();
                value.derivativeDefined = true;
                return value;
            },
            "derived objective/derivative"
        ) || numericalStateContractFailed;

    const spin::SpinNewtonResult overflowedNewtonSystem =
        spin::solveSpinMomentSystem(
            spin::SpinMomentProblem::staticBalance(
                spin::StaticObjective::LowerCenterOfMass
            ),
            settings,
            [](const spin::QuadricCoefficients10&) {
                spin::SpinMomentMapEvaluation value;
                value.rawMoments(0) = 3.162277660168379e-78;
                value.momentJacobian.setIdentity();
                value.momentJacobian(0, 0) = 1e154;
                value.derivativeDefined = true;
                return value;
            },
            {initialPoint()}
        );
    if (
        overflowedNewtonSystem.status !=
            spin::SpinSolverStatus::OptimizationDidNotConverge ||
        overflowedNewtonSystem.functionEvaluations != 1 ||
        overflowedNewtonSystem.trace.size() != 1 ||
        overflowedNewtonSystem.events.back().reason !=
            "newton_system_invalid"
    )
    {
        std::cerr << "[DETAIL] overflowed Newton system was not rejected before factorization"
                  << " status=" << spin::spinSolverStatusToken(
                         overflowedNewtonSystem.status)
                  << " evals=" << overflowedNewtonSystem.functionEvaluations
                  << " trace=" << overflowedNewtonSystem.trace.size()
                  << " last=" << overflowedNewtonSystem.events.back().reason
                  << '\n';
        numericalStateContractFailed = true;
    }
    if (numericalStateContractFailed)
    {
        return fail("non-finite numerical states crossed acceptance boundaries");
    }

    spin::SpinSolverSettings limited = settings;
    limited.maximumIterations = 1;
    limited.maximumFunctionEvaluationsPerRestart = 1;
    const spin::SpinNewtonResult exhausted = spin::solveSpinMomentSystem(
        problem,
        limited,
        identityMoments,
        {initialPoint()}
    );
    if (exhausted.status != spin::SpinSolverStatus::ResourceLimitReached)
    {
        return fail("iteration/evaluation limit did not return resource status");
    }
    if (
        exhausted.functionEvaluations != 1 || exhausted.events.back().reason !=
            "maximum_evaluations_reached"
    )
    {
        return fail("evaluation budget terminal evidence is incomplete");
    }

    spin::SpinSolverSettings iterationOnly = settings;
    iterationOnly.maximumIterations = 1;
    iterationOnly.maximumFunctionEvaluationsPerRestart = 60;
    const spin::SpinNewtonResult iterationExhausted =
        spin::solveSpinMomentSystem(
            problem,
            iterationOnly,
            identityMoments,
            {initialPoint()}
        );
    if (
        iterationExhausted.status != spin::SpinSolverStatus::ResourceLimitReached ||
        iterationExhausted.events.back().reason != "maximum_iterations_reached"
    )
    {
        return fail("iteration budget terminal evidence is incomplete");
    }

    spin::SpinSolverSettings searchExhaustion = settings;
    searchExhaustion.maximumIterations = 1;
    searchExhaustion.maximumFunctionEvaluationsPerRestart = 100;
    const spin::SpinNewtonResult searchExhausted =
        spin::solveSpinMomentSystem(
            problem,
            searchExhaustion,
            [](const spin::QuadricCoefficients10& quadric) {
                spin::SpinMomentMapEvaluation value = identityMoments(quadric);
                value.rawMoments(0) = 1.0;
                value.rawMoments(1) = 1e50;
                value.rawMoments(2) = 1e50;
                return value;
            },
            {initialPoint()}
        );
    if (
        searchExhausted.status != spin::SpinSolverStatus::LineSearchExhausted ||
        searchExhausted.events.back().reason != "line_search_exhausted" ||
        searchExhausted.functionEvaluations >=
            searchExhaustion.maximumFunctionEvaluationsPerRestart
    )
    {
        return fail("line-search exhaustion terminal evidence is incomplete");
    }

    std::cout << "[PASS] scenario=quadric-hollowing-solver iterations="
              << solved.trace.size() << " residual="
              << solved.normalizedEquationResidual << " trace="
              << solved.traceHash << '\n';
    return EXIT_SUCCESS;
}
