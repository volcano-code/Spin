#include "reachability/CarvingReachability.h"

#include <Eigen/Core>

#include <cstdlib>
#include <cmath>
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

reachability::CarvingReachabilityInput unreachableFixedShellFixture(
    int cellCount = 1
)
{
    reachability::CarvingReachabilityInput input;
    input.fixedMass = 1.0;
    input.fixedFirstMoment2D = Eigen::Vector2d(2.0, 0.0);
    input.safePolygon2D = scaledSquare(0.25);
    for (int index = 0; index < cellCount; ++index)
    {
        input.variableCells.push_back(reachability::ReachabilityCell{
            100 + index,
            1.0,
            Eigen::Vector2d(3.0 + index, 0.0)
        });
    }
    return input;
}

void fullInequalitySystem(
    const reachability::CarvingReachabilityInput& input,
    Eigen::MatrixXd& C,
    Eigen::VectorXd& d
)
{
    const reachability::Halfspaces2D halfspaces =
        reachability::safePolygonHalfspaces(input.safePolygon2D);
    const Eigen::Index n = static_cast<Eigen::Index>(
        input.variableCells.size()
    );
    C = Eigen::MatrixXd::Zero(halfspaces.A.rows() + 2 * n, n);
    d = Eigen::VectorXd::Zero(C.rows());
    for (Eigen::Index row = 0; row < halfspaces.A.rows(); ++row)
    {
        const Eigen::Vector2d normal = halfspaces.A.row(row).transpose();
        const double offset = halfspaces.b(row);
        d(row) = offset * input.fixedMass -
            normal.dot(input.fixedFirstMoment2D);
        for (Eigen::Index column = 0; column < n; ++column)
        {
            const reachability::ReachabilityCell& cell =
                input.variableCells[static_cast<std::size_t>(column)];
            C(row, column) = normal.dot(cell.firstMoment2D) -
                offset * cell.mass;
        }
    }
    for (Eigen::Index column = 0; column < n; ++column)
    {
        C(halfspaces.A.rows() + column, column) = 1.0;
        d(halfspaces.A.rows() + column) = 1.0;
        C(halfspaces.A.rows() + n + column, column) = -1.0;
    }
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

int runInfeasible()
{
    const reachability::CarvingReachabilityResult result =
        reachability::CarvingReachabilityAnalyzer().analyze(
            unreachableFixedShellFixture()
        );
    if (result.status != reachability::ReachabilityStatus::ProvenUnreachable)
    {
        return fail(
            "infeasible",
            std::string("unreachable shell lacks a proof; primal=") +
                result.primalSolverStatus + "; certificate=" +
                result.certificateSolverStatus + "; contradiction=" +
                std::to_string(result.certificate.contradiction) +
                "; stationarity=" +
                std::to_string(
                    result.certificate.stationarityInfinityNorm
                ) + "; weights=" +
                ([&] {
                    std::string text;
                    for (Eigen::Index index = 0;
                         index < result.certificate.weights.size();
                         ++index)
                    {
                        if (!text.empty()) text += ',';
                        text += std::to_string(
                            result.certificate.weights(index)
                        );
                    }
                    return text;
                })()
        );
    }
    if (
        !result.certificate.valid ||
        result.certificate.weights.size() == 0 ||
        result.certificate.minimumWeight < -1e-10 ||
        result.certificate.stationarityInfinityNorm > 1e-9 ||
        result.certificate.normalizationResidual > 1e-9 ||
        result.certificate.contradiction >= -1e-9
    )
    {
        return fail("infeasible", "reported Farkas certificate is invalid");
    }
    return EXIT_SUCCESS;
}

int runCertificate()
{
    Eigen::MatrixXd C;
    Eigen::VectorXd d;
    fullInequalitySystem(unreachableFixedShellFixture(), C, d);
    const reachability::ReachabilityCertificate invalid =
        reachability::validateFarkasCertificate(
            C,
            d,
            Eigen::VectorXd::Ones(C.rows()),
            1e-9
        );
    if (invalid.valid)
    {
        return fail("certificate", "invalid all-ones proof was accepted");
    }
    return EXIT_SUCCESS;
}

bool binaryStateIsFeasible(
    const Eigen::MatrixXd& C,
    const Eigen::VectorXd& d,
    unsigned int state
)
{
    Eigen::VectorXd retained(C.cols());
    for (Eigen::Index column = 0; column < C.cols(); ++column)
    {
        retained(column) = (state >> column) & 1U;
    }
    return (C * retained - d).maxCoeff() <= 1e-10;
}

bool ternaryGridHasFeasiblePoint(
    const Eigen::MatrixXd& C,
    const Eigen::VectorXd& d
)
{
    std::size_t stateCount = 1;
    for (Eigen::Index column = 0; column < C.cols(); ++column)
    {
        stateCount *= 3;
    }
    for (std::size_t state = 0; state < stateCount; ++state)
    {
        std::size_t encoded = state;
        Eigen::VectorXd retained(C.cols());
        for (Eigen::Index column = 0; column < C.cols(); ++column)
        {
            retained(column) = 0.5 * static_cast<double>(encoded % 3);
            encoded /= 3;
        }
        if ((C * retained - d).maxCoeff() <= 1e-10) return true;
    }
    return false;
}

int runExhaustive()
{
    const reachability::CarvingReachabilityAnalyzer analyzer;
    for (int cellCount = 2; cellCount <= 8; ++cellCount)
    {
        const reachability::CarvingReachabilityInput input =
            unreachableFixedShellFixture(cellCount);
        Eigen::MatrixXd C;
        Eigen::VectorXd d;
        fullInequalitySystem(input, C, d);
        bool binaryFeasible = false;
        for (unsigned int state = 0; state < (1U << cellCount); ++state)
        {
            binaryFeasible = binaryFeasible ||
                binaryStateIsFeasible(C, d, state);
        }
        const reachability::CarvingReachabilityResult result =
            analyzer.analyze(input);
        if (
            binaryFeasible &&
            result.status == reachability::ReachabilityStatus::ProvenUnreachable
        )
        {
            return fail("exhaustive", "binary witness contradicted proof");
        }
        if (
            result.status == reachability::ReachabilityStatus::ProvenUnreachable &&
            ternaryGridHasFeasiblePoint(C, d)
        )
        {
            return fail("exhaustive", "sampled continuous witness contradicted proof");
        }
    }

    reachability::CarvingReachabilityInput fractional;
    fractional.fixedMass = 1.0;
    fractional.safePolygon2D = scaledSquare(0.1);
    fractional.safePolygon2D.col(0).array() += 2.0 / 3.0;
    fractional.variableCells = {
        reachability::ReachabilityCell{1, 1.0, Eigen::Vector2d(2.0, 0.0)},
        reachability::ReachabilityCell{2, 0.0, Eigen::Vector2d::Zero()}
    };
    Eigen::MatrixXd C;
    Eigen::VectorXd d;
    fullInequalitySystem(fractional, C, d);
    bool binaryFeasible = false;
    for (unsigned int state = 0; state < 4; ++state)
    {
        binaryFeasible = binaryFeasible || binaryStateIsFeasible(C, d, state);
    }
    const reachability::CarvingReachabilityResult fractionalResult =
        analyzer.analyze(fractional);
    if (
        binaryFeasible ||
        fractionalResult.status !=
            reachability::ReachabilityStatus::RelaxationReachable
    )
    {
        return fail(
            "exhaustive",
            std::string(
                "fractional-only solution was not kept distinct from binary "
                "reachability; binary="
            ) + (binaryFeasible ? "true" : "false") +
                "; status=" + reachability::toString(fractionalResult.status) +
                "; primal=" + fractionalResult.primalSolverStatus +
                "; violation=" +
                std::to_string(fractionalResult.maximumPrimalViolation) +
                "; center=" + std::to_string(fractionalResult.finalCenter2D.x())
        );
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
        if (scenario == "infeasible") return runInfeasible();
        if (scenario == "certificate") return runCertificate();
        if (scenario == "exhaustive") return runExhaustive();
    }
    catch (const std::exception& error)
    {
        return fail(scenario, error.what());
    }
    return fail(scenario, "unknown scenario");
}
