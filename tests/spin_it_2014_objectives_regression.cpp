#include "spin2014/AdaptiveVoxelSpinSolver.h"
#include "spin2014/SpinIt2014MomentProblem.h"

#include "voxel/VoxelDomain.h"

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>


namespace
{


using Objective = spin2014::Spin2014Objective;
using Problem = spin2014::SpinIt2014MomentProblem;
using Vector10 = Problem::Vector10;


int fail(const std::string& message)
{
    std::cerr << "[FAIL] scenario=spin-it-2014-objectives " << message << '\n';
    return EXIT_FAILURE;
}


Vector10 moments()
{
    Vector10 value;
    // Native raw order: [M,Lx,Ly,Lz,Qxx,Qyy,Qzz,Qyz,Qxz,Qxy].
    value << 4.0, 0.4, -0.8, 2.0, 10.0, 14.0, 6.0, 0.3, -0.4, 0.2;
    return value;
}


std::shared_ptr<const voxel::VoxelDomain> fixture()
{
    voxel::GridSpec grid;
    grid.dimensions = Eigen::Vector3i(1, 1, 1);
    grid.minimum = Eigen::Vector3d(-1.0, -1.0, 0.0);
    grid.maximum = Eigen::Vector3d(1.0, 1.0, 2.0);
    return voxel::VoxelDomain::create(grid, {1U}, {}, {1U});
}


bool throwsInvalid(const spin2014::AdaptiveVoxelSpinRequest& request)
{
    try
    {
        spin2014::validateAdaptiveVoxelSpinRequest(request);
    }
    catch (const std::invalid_argument&)
    {
        return true;
    }
    return false;
}


double maximumGradientError(
    const Problem& problem,
    const Vector10& point,
    double step
)
{
    const auto exact = problem.evaluate(point);
    double error = 0.0;
    for (int coordinate = 0; coordinate < 10; ++coordinate)
    {
        Vector10 plus = point;
        Vector10 minus = point;
        plus(coordinate) += step;
        minus(coordinate) -= step;
        const double finiteDifference =
            (problem.evaluate(plus).objective - problem.evaluate(minus).objective) /
            (2.0 * step);
        error = std::max(
            error,
            std::abs(finiteDifference - exact.objectiveGradient(coordinate))
        );
        if (
            (coordinate == 7 || coordinate == 8) &&
            (std::abs(finiteDifference) > 1e-10 ||
             std::abs(exact.objectiveGradient(coordinate)) > 1e-14)
        )
        {
            return std::numeric_limits<double>::infinity();
        }
    }
    return error;
}


double maximumHessianError(
    const Problem& problem,
    const Vector10& point,
    double step
)
{
    const auto exact = problem.evaluate(point);
    double error = 0.0;
    for (int coordinate = 0; coordinate < 10; ++coordinate)
    {
        Vector10 plus = point;
        Vector10 minus = point;
        plus(coordinate) += step;
        minus(coordinate) -= step;
        const Vector10 finiteDifference =
            (problem.evaluate(plus).objectiveGradient -
             problem.evaluate(minus).objectiveGradient) /
            (2.0 * step);
        error = std::max(
            error,
            (finiteDifference - exact.objectiveHessian.col(coordinate))
                .cwiseAbs().maxCoeff()
        );
    }
    return error;
}


} // namespace


int main()
{
    const Vector10 raw = moments();
    constexpr double gammaI = 2.5;
    constexpr double gammaC = 0.7;
    struct Case
    {
        Objective objective;
        double expected;
        spin::SpinMode mode;
        int constraintCount;
        const char* token;
    };
    // Hand-derived from Ia=18.84, Ib=14.96, Ic=23.8, Iab=-0.28,
    // f_yo/gammaI=1.0220040957559493, ell=0.5 and ell*M=2.
    const Case cases[] = {
        {Objective::YoYo, 2.5550102393898735, spin::SpinMode::YoYo, 5, "YOYO"},
        {Objective::TopInertiaOnly, 2.5550102393898735,
         spin::SpinMode::SpinningTop, 4, "TOP_INERTIA_ONLY"},
        {Objective::TopLoweringOnly, 2.7300102393898733,
         spin::SpinMode::SpinningTop, 4, "TOP_LOWERING_ONLY"},
        {Objective::TopMassAndLowering, 5.355010239389873,
         spin::SpinMode::SpinningTop, 4, "TOP_MASS_AND_LOWERING"},
        {Objective::StaticLowering, 0.175,
         spin::SpinMode::StaticBalance, 2, "STATIC_LOWERING"},
    };

    for (const Case& item : cases)
    {
        const Problem problem = Problem::fromObjective(
            item.objective, gammaI, gammaC
        );
        const auto value = problem.evaluate(raw);
        if (
            std::abs(value.objective - item.expected) > 1e-12 ||
            problem.objective() != item.objective ||
            problem.mode() != item.mode ||
            problem.inertiaWeight() != gammaI ||
            problem.centerWeight() != gammaC ||
            value.constraints.size() != item.constraintCount ||
            std::string(spin2014::spin2014ObjectiveToken(item.objective)) !=
                item.token
        )
        {
            return fail(std::string("formula/API mismatch for ") + item.token);
        }
        const double steps[] = {1e-3, 3e-4, 1e-4};
        double bestGradientError = std::numeric_limits<double>::infinity();
        double bestHessianError = std::numeric_limits<double>::infinity();
        for (double step : steps)
        {
            bestGradientError = std::min(
                bestGradientError,
                maximumGradientError(problem, raw, step)
            );
            bestHessianError = std::min(
                bestHessianError,
                maximumHessianError(problem, raw, step)
            );
        }
        if (bestGradientError > 2e-7 || bestHessianError > 2e-6)
        {
            return fail(std::string("finite differences disagree for ") + item.token);
        }
    }

    const double density = 3.0;
    const Vector10 scaled = density * raw;
    const Problem inertiaOnly = Problem::fromObjective(
        Objective::TopInertiaOnly, gammaI, gammaC
    );
    const Problem loweringOnly = Problem::fromObjective(
        Objective::TopLoweringOnly, gammaI, gammaC
    );
    const Problem full = Problem::fromObjective(
        Objective::TopMassAndLowering, gammaI, gammaC
    );
    const Problem staticLowering = Problem::fromObjective(
        Objective::StaticLowering, gammaI, gammaC
    );
    if (
        std::abs(loweringOnly.evaluate(scaled).objective -
                 loweringOnly.evaluate(raw).objective) > 1e-12 ||
        std::abs(staticLowering.evaluate(scaled).objective -
                 staticLowering.evaluate(raw).objective) > 1e-12
    )
    {
        return fail("height-only penalty changed under uniform density scaling");
    }
    const double rawPenalty = full.evaluate(raw).objective -
        inertiaOnly.evaluate(raw).objective;
    const double scaledPenalty = full.evaluate(scaled).objective -
        inertiaOnly.evaluate(scaled).objective;
    if (std::abs(scaledPenalty / rawPenalty - density * density) > 1e-12)
    {
        return fail("mass-and-lowering penalty did not scale with density squared");
    }

    if (
        Problem::yoYo().objective() != Objective::YoYo ||
        Problem::yoYo().inertiaWeight() != 1.0 ||
        Problem::spinningTop(0.7).objective() != Objective::TopMassAndLowering ||
        Problem::staticBalance(0.7).objective() != Objective::StaticLowering ||
        std::abs(Problem::spinningTop(0.7).evaluate(raw).objective -
                 3.822004095755949) > 1e-12
    )
    {
        return fail("legacy factory defaults changed or weights are not independent");
    }

    spin2014::AdaptiveVoxelSpinRequest request;
    request.domain = fixture();
    request.initialLevel = 0;
    request.maximumLevel = 0;
    request.maximumAdaptivePasses = 1;
    const std::string legacyHash =
        spin2014::computeAdaptiveVoxelSpinRequestHash(request);
    std::set<std::string> selectorHashes;
    for (Objective objective : {
             Objective::TopInertiaOnly,
             Objective::TopLoweringOnly,
             Objective::TopMassAndLowering})
    {
        request.objective = objective;
        selectorHashes.insert(
            spin2014::computeAdaptiveVoxelSpinRequestHash(request)
        );
    }
    if (selectorHashes.size() != 3U || selectorHashes.count(legacyHash) != 0U)
    {
        return fail("selector or legacy resolution source is absent from request hash");
    }
    const std::string beforeWeight =
        spin2014::computeAdaptiveVoxelSpinRequestHash(request);
    request.inertiaWeight = 1.25;
    const std::string afterInertia =
        spin2014::computeAdaptiveVoxelSpinRequestHash(request);
    request.centerOfMassWeight = 1.25;
    const std::string afterCenter =
        spin2014::computeAdaptiveVoxelSpinRequestHash(request);
    if (beforeWeight == afterInertia || afterInertia == afterCenter)
    {
        return fail("independent weights are absent from request hash");
    }

    request.mode = spin::SpinMode::StaticBalance;
    request.objective = Objective::TopInertiaOnly;
    if (!throwsInvalid(request))
    {
        return fail("selector/mode mismatch was accepted");
    }
    request.objective.reset();
    request.inertiaWeight = -1.0;
    if (!throwsInvalid(request))
    {
        return fail("negative inertia weight was accepted");
    }
    request.inertiaWeight = std::numeric_limits<double>::quiet_NaN();
    if (!throwsInvalid(request))
    {
        return fail("nonfinite inertia weight was accepted");
    }

    std::cout << "[PASS] scenario=spin-it-2014-objectives five-objective API"
              << " derivatives/density/hash\n";
    return EXIT_SUCCESS;
}
