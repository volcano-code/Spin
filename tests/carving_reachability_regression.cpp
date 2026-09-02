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
    }
    catch (const std::exception& error)
    {
        return fail(scenario, error.what());
    }
    return fail(scenario, "unknown scenario");
}
