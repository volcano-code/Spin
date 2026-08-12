#include "physics/BalanceObjective.h"
#include "stability/BalanceTargetSchedule.h"
#include "stability/SupportRegion.h"

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>


namespace
{


constexpr double TOLERANCE = 1.0e-12;


int fail(const std::string& scenario, const std::string& reason)
{
    std::cerr
        << "[FAIL] scenario=" << scenario
        << " reason=" << reason << '\n';
    return EXIT_FAILURE;
}


bool close(double first, double second)
{
    return
        std::isfinite(first) &&
        std::isfinite(second) &&
        std::abs(first - second) <=
            TOLERANCE *
            std::max({1.0, std::abs(first), std::abs(second)});
}


Eigen::MatrixXd squareSupportVertices()
{
    Eigen::MatrixXd vertices(5, 3);
    vertices <<
        -2.0, 0.0, -2.0,
         2.0, 0.0, -2.0,
         2.0, 0.0,  2.0,
        -2.0, 0.0,  2.0,
         0.0, 0.2,  0.0;
    return vertices;
}


stability::SupportRegion makeSquareSupport()
{
    return stability::SupportRegion::fromMesh(
        squareSupportVertices(),
        Eigen::Vector3d(0.0, -1.0, 0.0),
        0.0
    );
}


Eigen::Vector3d pointAbove(
    const stability::SupportRegion& support,
    const Eigen::Vector2d& projection,
    double height
)
{
    return
        support.fromPlaneCoordinates(projection) -
        height * support.gravity();
}


stability::BalanceTargetSchedule oneStandingSchedule(
    const stability::SupportRegion& support,
    double gamma
)
{
    stability::BalanceTargetSpecifications specifications;
    specifications.emplace_back(
        stability::StandingTargetSpecification(support, gamma)
    );
    return stability::BalanceTargetSchedule(specifications);
}


int runRefresh(const std::string& scenario)
{
    const stability::SupportRegion support = makeSquareSupport();
    const Eigen::MatrixXd fixedPolygon = support.supportPolygon2D();
    stability::BalanceTargetSchedule schedule = oneStandingSchedule(
        support,
        std::atan(0.5)
    );

    if (schedule.hasSnapshot() || schedule.generation() != 0)
    {
        return fail(scenario, "a new schedule already has a snapshot");
    }

    const Eigen::Vector3d unsafeCenter = pointAbove(
        support,
        Eigen::Vector2d(1.5, 0.0),
        2.0
    );
    const stability::BalanceTargetSnapshot first =
        schedule.refresh(unsafeCenter);

    if (
        first.generation != 1 ||
        !first.contextChanged ||
        first.allReached ||
        first.objectives.size() != 1 ||
        first.constraints.size() != 1 ||
        first.constraints[0].type !=
            stability::BalanceConstraintType::Standing ||
        first.constraints[0].standingState !=
            stability::StandingState::FeasibleButUnsafe ||
        first.constraints[0].requestedMarginAttainable != true ||
        first.constraints[0].reached ||
        !support.toPlaneCoordinates(
            first.objectives[0].target()
        ).isApprox(Eigen::Vector2d(1.0, 0.0), TOLERANCE)
    )
    {
        return fail(scenario, "initial standing snapshot is incorrect");
    }

    const stability::BalanceTargetSnapshot repeated =
        schedule.refresh(unsafeCenter);
    if (
        repeated.generation != 1 ||
        repeated.contextChanged ||
        repeated.constraints[0].standingState !=
            stability::StandingState::FeasibleButUnsafe
    )
    {
        return fail(scenario, "identical refresh changed objective generation");
    }

    const Eigen::Vector3d stableCenter = pointAbove(
        support,
        Eigen::Vector2d(0.25, 0.0),
        2.0
    );
    const stability::BalanceTargetSnapshot stable =
        schedule.refresh(stableCenter);
    if (
        stable.generation != 2 ||
        !stable.contextChanged ||
        !stable.allReached ||
        stable.constraints[0].standingState !=
            stability::StandingState::Stable ||
        !stable.constraints[0].reached ||
        !support.toPlaneCoordinates(
            stable.objectives[0].target()
        ).isApprox(Eigen::Vector2d(0.25, 0.0), TOLERANCE) ||
        !close(
            aggregateProjectedCenterOfMassEnergy(
                stableCenter,
                stable.objectives
            ),
            0.0
        )
    )
    {
        return fail(scenario, "fresh stable target was not committed");
    }

    /* State/height may change without changing the objective target. */
    const stability::BalanceTargetSnapshot sameTarget = schedule.refresh(
        pointAbove(support, Eigen::Vector2d(0.25, 0.0), 1.5)
    );
    if (
        sameTarget.generation != 2 ||
        sameTarget.contextChanged ||
        !sameTarget.allReached ||
        !schedule.currentSnapshot().objectives[0].target().isApprox(
            sameTarget.objectives[0].target(),
            0.0
        ) ||
        !std::get<stability::StandingTargetSpecification>(
            schedule.specifications()[0]
        ).support.supportPolygon2D().isApprox(fixedPolygon, 0.0)
    )
    {
        return fail(scenario, "state-only refresh changed target/support context");
    }

    std::cout
        << "[PASS] scenario=" << scenario
        << " generations=" << schedule.generation() << '\n';
    return EXIT_SUCCESS;
}


int runStability(const std::string& scenario)
{
    const stability::SupportRegion support = makeSquareSupport();
    stability::BalanceTargetSchedule standing = oneStandingSchedule(
        support,
        std::atan(0.5)
    );

    struct StandingFixture
    {
        double projection;
        stability::StandingState state;
        bool reached;
    };
    const StandingFixture standingFixtures[] =
    {
        {2.5, stability::StandingState::OutsideSupport, false},
        {1.5, stability::StandingState::FeasibleButUnsafe, false},
        {1.0, stability::StandingState::Stable, true},
        {0.0, stability::StandingState::Stable, true}
    };

    for (const StandingFixture& fixture : standingFixtures)
    {
        const stability::BalanceTargetSnapshot snapshot = standing.refresh(
            pointAbove(
                support,
                Eigen::Vector2d(fixture.projection, 0.0),
                2.0
            )
        );
        if (
            snapshot.constraints[0].standingState != fixture.state ||
            snapshot.constraints[0].reached != fixture.reached ||
            snapshot.allReached != fixture.reached
        )
        {
            return fail(scenario, "standing terminal predicate is incorrect");
        }
    }

    stability::BalanceTargetSchedule unattainable = oneStandingSchedule(
        support,
        std::atan(1.5)
    );
    const stability::BalanceTargetSnapshot impossible = unattainable.refresh(
        pointAbove(support, Eigen::Vector2d::Zero(), 2.0)
    );
    if (
        impossible.constraints[0].requestedMarginAttainable != false ||
        impossible.constraints[0].standingState !=
            stability::StandingState::FeasibleButUnsafe ||
        impossible.constraints[0].reached ||
        impossible.allReached
    )
    {
        return fail(scenario, "unattainable standing margin was reported stable");
    }

    stability::BalanceTargetSpecifications suspensionSpecifications;
    suspensionSpecifications.emplace_back(
        stability::SuspensionTargetSpecification(
            Eigen::Vector3d::Zero(),
            Eigen::Vector3d(0.0, -1.0, 0.0),
            std::atan(0.25)
        )
    );
    stability::BalanceTargetSchedule suspension(
        suspensionSpecifications
    );

    const Eigen::Vector3d suspensionQueries[] =
    {
        Eigen::Vector3d(0.0, -2.0, 0.0),
        Eigen::Vector3d(0.5, -2.0, 0.0),
        Eigen::Vector3d(0.5001, -2.0, 0.0),
        Eigen::Vector3d(0.0, 2.0, 0.0),
        Eigen::Vector3d::Zero()
    };
    const bool suspensionExpected[] =
    {
        true, true, false, false, false
    };

    for (std::size_t index = 0; index < 5; ++index)
    {
        const stability::BalanceTargetSnapshot snapshot =
            suspension.refresh(suspensionQueries[index]);
        if (
            snapshot.constraints[0].reached !=
                suspensionExpected[index] ||
            snapshot.allReached != suspensionExpected[index] ||
            snapshot.constraints[0].type !=
                stability::BalanceConstraintType::Suspension ||
            snapshot.constraints[0].standingState.has_value() ||
            snapshot.constraints[0].requestedMarginAttainable.has_value()
        )
        {
            return fail(scenario, "suspension angular predicate is incorrect");
        }

        if (
            index == 3 &&
            snapshot.constraints[0].suspensionDownstream != false
        )
        {
            return fail(scenario, "upstream suspension state was accepted");
        }
        if (
            index == 4 &&
            (
                snapshot.constraints[0].suspensionDownstream.has_value() ||
                snapshot.constraints[0].suspensionDeviationRadians.has_value()
            )
        )
        {
            return fail(scenario, "zero suspension offset has a direction");
        }
    }

    /* Suspension state changes do not change its fixed objective context. */
    if (
        suspension.generation() != 1 ||
        suspension.currentSnapshot().contextChanged
    )
    {
        return fail(scenario, "fixed suspension target changed generation");
    }

    std::cout
        << "[PASS] scenario=" << scenario
        << " safe_boundary=closed suspension_boundary=closed\n";
    return EXIT_SUCCESS;
}


int runMulti(const std::string& scenario)
{
    const stability::SupportRegion support = makeSquareSupport();
    stability::BalanceTargetSpecifications specifications;
    specifications.emplace_back(
        stability::StandingTargetSpecification(
            support,
            std::atan(0.5)
        )
    );
    specifications.emplace_back(
        stability::StandingTargetSpecification(
            support,
            std::atan(0.75)
        )
    );
    stability::BalanceTargetSchedule schedule(specifications);

    const Eigen::Vector3d partiallyStableCenter = pointAbove(
        support,
        Eigen::Vector2d(0.75, 0.0),
        2.0
    );
    const stability::BalanceTargetSnapshot partial = schedule.refresh(
        partiallyStableCenter
    );
    if (
        partial.constraints.size() != 2 ||
        partial.objectives.size() != 2 ||
        !partial.constraints[0].reached ||
        partial.constraints[1].reached ||
        partial.allReached ||
        !close(
            aggregateProjectedCenterOfMassEnergy(
                partiallyStableCenter,
                partial.objectives
            ),
            0.03125
        )
    )
    {
        return fail(scenario, "multi-standing all-of predicate is incorrect");
    }

    const Eigen::Vector3d fullyStableCenter = pointAbove(
        support,
        Eigen::Vector2d(0.5, 0.0),
        2.0
    );
    const stability::BalanceTargetSnapshot complete = schedule.refresh(
        fullyStableCenter
    );
    if (
        !complete.constraints[0].reached ||
        !complete.constraints[1].reached ||
        !complete.allReached ||
        !close(
            aggregateProjectedCenterOfMassEnergy(
                fullyStableCenter,
                complete.objectives
            ),
            0.0
        )
    )
    {
        return fail(scenario, "multi-standing completion is incorrect");
    }

    /* A mixed schedule preserves specification order exactly. */
    const Eigen::Vector3d attachment = pointAbove(
        support,
        Eigen::Vector2d(0.5, 0.0),
        4.0
    );
    stability::BalanceTargetSpecifications mixedSpecifications;
    mixedSpecifications.emplace_back(
        stability::SuspensionTargetSpecification(
            attachment,
            support.gravity(),
            0.0
        )
    );
    mixedSpecifications.emplace_back(
        stability::StandingTargetSpecification(
            support,
            std::atan(0.75)
        )
    );
    mixedSpecifications.emplace_back(
        stability::StandingTargetSpecification(
            support,
            0.0
        )
    );
    stability::BalanceTargetSchedule mixed(mixedSpecifications);
    const stability::BalanceTargetSnapshot mixedSnapshot =
        mixed.refresh(fullyStableCenter);

    if (
        mixedSnapshot.constraints.size() != 3 ||
        mixedSnapshot.objectives.size() != 3 ||
        mixedSnapshot.constraints[0].type !=
            stability::BalanceConstraintType::Suspension ||
        mixedSnapshot.constraints[1].type !=
            stability::BalanceConstraintType::Standing ||
        mixedSnapshot.constraints[2].type !=
            stability::BalanceConstraintType::Standing ||
        !mixedSnapshot.objectives[0].target().isApprox(attachment, 0.0) ||
        !support.toPlaneCoordinates(
            mixedSnapshot.objectives[1].target()
        ).isApprox(Eigen::Vector2d(0.5, 0.0), TOLERANCE) ||
        !mixedSnapshot.allReached
    )
    {
        return fail(scenario, "mixed objective order/all-of state changed");
    }

    std::cout
        << "[PASS] scenario=" << scenario
        << " partial_energy=0.03125 constraints=3\n";
    return EXIT_SUCCESS;
}


} // namespace


int main(int argc, char** argv)
{
    if (argc != 2)
    {
        return fail(
            "unknown",
            "usage: balance-target-schedule-regression <refresh|stability|multi>"
        );
    }

    const std::string scenario = argv[1];
    if (scenario == "refresh")
    {
        return runRefresh(scenario);
    }
    if (scenario == "stability")
    {
        return runStability(scenario);
    }
    if (scenario == "multi")
    {
        return runMulti(scenario);
    }
    return fail(scenario, "unknown scenario");
}
