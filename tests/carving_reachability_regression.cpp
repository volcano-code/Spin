#include "reachability/CarvingReachability.h"

#include <Eigen/Core>

#include <cstdlib>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>

namespace
{

int fail(const std::string& scenario, const std::string& message)
{
    std::cerr << "[FAIL] scenario=" << scenario << ' ' << message << '\n';
    return EXIT_FAILURE;
}

Eigen::MatrixXd square(bool clockwise = false)
{
    Eigen::MatrixXd polygon(4, 2);
    polygon << -1.0, -1.0,
                1.0, -1.0,
                1.0,  1.0,
               -1.0,  1.0;
    if (clockwise)
    {
        polygon.colwise().reverseInPlace();
    }
    return polygon;
}

Eigen::MatrixXd scaledSquare(double halfExtent)
{
    return halfExtent * square();
}

reachability::CarvingReachabilityInput alreadySafeFixture()
{
    reachability::CarvingReachabilityInput input;
    input.fixedMass = 1.0;
    input.fixedFirstMoment2D = Eigen::Vector2d::Zero();
    input.safePolygon2D = scaledSquare(0.25);
    input.variableCells = {
        reachability::ReachabilityCell{1, 1.0, Eigen::Vector2d::Zero()}
    };
    return input;
}

reachability::CarvingReachabilityInput oneMovableCellFixture()
{
    reachability::CarvingReachabilityInput input;
    input.fixedMass = 1.0;
    input.fixedFirstMoment2D = Eigen::Vector2d::Zero();
    input.safePolygon2D = scaledSquare(0.25);
    input.variableCells = {
        reachability::ReachabilityCell{2, 1.0, Eigen::Vector2d(2.0, 0.0)}
    };
    return input;
}

bool throwsInvalidArgument(const std::function<void()>& action)
{
    try
    {
        action();
    }
    catch (const std::invalid_argument&)
    {
        return true;
    }
    return false;
}

int runInput()
{
    const reachability::Halfspaces2D halfspaces =
        reachability::safePolygonHalfspaces(square());
    if (!reachability::contains(halfspaces, Eigen::Vector2d::Zero(), 1e-12))
    {
        return fail("input", "square does not contain its center");
    }
    if (!reachability::contains(
            halfspaces,
            Eigen::Vector2d(1.0, 0.0),
            1e-12
        ))
    {
        return fail("input", "square excludes its boundary");
    }
    if (reachability::contains(
            halfspaces,
            Eigen::Vector2d(1.01, 0.0),
            1e-12
        ))
    {
        return fail("input", "square contains an exterior point");
    }
    if (!throwsInvalidArgument([&] {
            static_cast<void>(reachability::safePolygonHalfspaces(
                square(true)
            ));
        }))
    {
        return fail("input", "clockwise polygon was accepted");
    }

    reachability::CarvingReachabilityInput input;
    input.fixedMass = 1.0;
    input.safePolygon2D = square();
    input.variableCells = {
        reachability::ReachabilityCell{7, 1.0, Eigen::Vector2d::Zero()}
    };
    reachability::validateInput(input);

    input.variableCells.push_back(input.variableCells.front());
    if (!throwsInvalidArgument([&] {
            reachability::validateInput(input);
        }))
    {
        return fail("input", "duplicate cell identifier was accepted");
    }
    return EXIT_SUCCESS;
}

int runAlreadySafe()
{
    const reachability::CarvingReachabilityResult result =
        reachability::CarvingReachabilityAnalyzer().analyze(
            alreadySafeFixture()
        );
    if (result.status != reachability::ReachabilityStatus::AlreadySafe)
    {
        return fail("already-safe", "safe full material was not recognized");
    }
    if (
        result.retainedFractions.size() != result.inputCellCount ||
        !result.retainedFractions.isApprox(
            Eigen::VectorXd::Ones(result.inputCellCount),
            1e-14
        )
    )
    {
        return fail("already-safe", "full-material witness is not all ones");
    }
    return EXIT_SUCCESS;
}

int runFeasible()
{
    const reachability::CarvingReachabilityAnalyzer analyzer;
    const reachability::CarvingReachabilityResult result =
        analyzer.analyze(oneMovableCellFixture());
    const reachability::CarvingReachabilityResult repeated =
        analyzer.analyze(oneMovableCellFixture());
    if (
        result.status !=
        reachability::ReachabilityStatus::RelaxationReachable
    )
    {
        return fail("feasible", "reachable relaxation was not classified");
    }
    if (
        result.retainedFractions.size() != 1 ||
        result.retainedFractions.minCoeff() < -1e-10 ||
        result.retainedFractions.maxCoeff() > 1.0 + 1e-10 ||
        result.maximumPrimalViolation > 1e-9
    )
    {
        return fail("feasible", "returned witness violates the primal LP");
    }
    if (!result.retainedFractions.isApprox(
            repeated.retainedFractions,
            1e-14
        ))
    {
        return fail("feasible", "LP witness changed between identical runs");
    }
    const reachability::Halfspaces2D halfspaces =
        reachability::safePolygonHalfspaces(
            oneMovableCellFixture().safePolygon2D
        );
    if (!reachability::contains(halfspaces, result.finalCenter2D, 1e-9))
    {
        return fail("feasible", "independently recomputed center is unsafe");
    }
    return EXIT_SUCCESS;
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 2)
    {
        std::cerr << "usage: carving-reachability-regression SCENARIO\n";
        return EXIT_FAILURE;
    }
    const std::string scenario(argv[1]);
    try
    {
        if (scenario == "input") return runInput();
        if (scenario == "already-safe") return runAlreadySafe();
        if (scenario == "feasible") return runFeasible();
    }
    catch (const std::exception& error)
    {
        return fail(scenario, error.what());
    }
    return fail(scenario, "unknown scenario");
}
