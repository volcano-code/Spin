#include "spin/QuadricHollowingSolver.h"

#include <Eigen/Core>

#include <cmath>
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


spin::QuadricCoefficients10 initialPoint()
{
    Vector10 values;
    values << 1.0, 0.12, -0.08, 0.09, -0.07,
        0.05, -0.04, 0.03, -0.02, 0.01;
    return spin::QuadricCoefficients10::fromUnnormalized(values);
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

    spin::SpinSolverSettings limited = settings;
    limited.maximumIterations = 1;
    limited.maximumFunctionEvaluationsPerRestart = 2;
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

    std::cout << "[PASS] scenario=quadric-hollowing-solver iterations="
              << solved.trace.size() << " residual="
              << solved.normalizedEquationResidual << " trace="
              << solved.traceHash << '\n';
    return EXIT_SUCCESS;
}
