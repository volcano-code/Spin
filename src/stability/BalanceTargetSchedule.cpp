#include "stability/BalanceTargetSchedule.h"

#include "physics/BalanceEnergy.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>


namespace stability
{
namespace
{


constexpr double HALF_PI =
    1.57079632679489661923132169163975144;


void mixSnapshotWord(std::uint64_t& hash, std::uint64_t word) noexcept
{
    for (int byte = 0; byte < 8; ++byte)
    {
        hash ^= (word >> (byte * 8)) & UINT64_C(0xff);
        hash *= UINT64_C(1099511628211);
    }
}


void mixSnapshotDouble(std::uint64_t& hash, double value) noexcept
{
    std::uint64_t bits = 0;
    static_assert(sizeof(bits) == sizeof(value), "unexpected double size");
    std::memcpy(&bits, &value, sizeof(bits));
    mixSnapshotWord(hash, bits);
}


void validateAngle(
    double angle,
    const char* message
)
{
    if (
        !std::isfinite(angle) ||
        angle < 0.0 ||
        angle >= HALF_PI
    )
    {
        throw std::invalid_argument(message);
    }
}


bool targetsEquivalent(
    const BalanceObjective& first,
    const BalanceObjective& second
)
{
    const double scale = std::max({
        1.0,
        first.target().cwiseAbs().maxCoeff(),
        second.target().cwiseAbs().maxCoeff()
    });
    const double tolerance =
        64.0 * std::numeric_limits<double>::epsilon() * scale;

    return
        (first.target() - second.target()).norm() <= tolerance &&
        first.gravity().isApprox(second.gravity(), 0.0);
}


bool objectiveContextsEquivalent(
    const BalanceObjectives& first,
    const BalanceObjectives& second
)
{
    if (first.size() != second.size())
    {
        return false;
    }

    for (std::size_t index = 0; index < first.size(); ++index)
    {
        if (!targetsEquivalent(first[index], second[index]))
        {
            return false;
        }
    }

    return true;
}


struct SuspensionState
{
    bool hasDirection;
    bool downstream;
    double deviationRadians;
    bool reached;
};


SuspensionState evaluateSuspension(
    const Eigen::Vector3d& centerOfMass,
    const SuspensionTargetSpecification& specification
)
{
    const Eigen::Vector3d displacement =
        centerOfMass - specification.objective.target();

    if (!displacement.allFinite())
    {
        throw std::overflow_error(
            "suspension attachment displacement overflowed"
        );
    }

    const double scale = displacement.cwiseAbs().maxCoeff();
    if (!(scale > 0.0))
    {
        return SuspensionState{false, false, 0.0, false};
    }

    const Eigen::Vector3d scaled = displacement / scale;
    const double scaledNorm = scaled.norm();
    if (!std::isfinite(scaledNorm) || !(scaledNorm > 0.0))
    {
        throw std::overflow_error(
            "suspension attachment direction overflowed"
        );
    }

    const Eigen::Vector3d direction = scaled / scaledNorm;
    const double parallel = std::max(
        -1.0,
        std::min(
            1.0,
            direction.dot(specification.objective.gravity())
        )
    );
    const Eigen::Vector3d perpendicular =
        direction - parallel * specification.objective.gravity();
    const double deviation = std::atan2(
        perpendicular.norm(),
        parallel
    );

    if (!std::isfinite(deviation))
    {
        throw std::overflow_error(
            "suspension angular deviation overflowed"
        );
    }

    const bool downstream = parallel > 0.0;
    const double angleTolerance =
        64.0 * std::numeric_limits<double>::epsilon() *
        std::max({
            1.0,
            deviation,
            specification.maximumDeviationRadians
        });

    return SuspensionState{
        true,
        downstream,
        deviation,
        downstream &&
            deviation <=
                specification.maximumDeviationRadians + angleTolerance
    };
}


} // namespace


std::uint64_t balanceTargetSnapshotContentHash(
    const BalanceTargetSnapshot& snapshot
) noexcept
{
    std::uint64_t hash = UINT64_C(1469598103934665603);
    mixSnapshotWord(hash, snapshot.generation);
    mixSnapshotWord(hash, snapshot.contextChanged ? 1U : 0U);
    mixSnapshotWord(hash, snapshot.allReached ? 1U : 0U);
    mixSnapshotWord(hash, snapshot.objectives.size());
    for (const BalanceObjective& objective : snapshot.objectives)
    {
        for (int axis = 0; axis < 3; ++axis)
        {
            mixSnapshotDouble(hash, objective.target()(axis));
            mixSnapshotDouble(hash, objective.gravity()(axis));
        }
    }
    mixSnapshotWord(hash, snapshot.constraints.size());
    for (const BalanceConstraintSnapshot& constraint : snapshot.constraints)
    {
        mixSnapshotWord(hash, static_cast<std::uint64_t>(constraint.type));
        mixSnapshotWord(hash, constraint.reached ? 1U : 0U);
        mixSnapshotWord(
            hash,
            constraint.standingState.has_value()
                ? 1U + static_cast<std::uint64_t>(*constraint.standingState)
                : 0U
        );
        mixSnapshotWord(
            hash,
            constraint.requestedMarginAttainable.has_value()
                ? (*constraint.requestedMarginAttainable ? 2U : 1U)
                : 0U
        );
        if (constraint.suspensionDeviationRadians.has_value())
        {
            mixSnapshotWord(hash, 1U);
            mixSnapshotDouble(hash, *constraint.suspensionDeviationRadians);
        }
        else
        {
            mixSnapshotWord(hash, 0U);
        }
        mixSnapshotWord(
            hash,
            constraint.suspensionDownstream.has_value()
                ? (*constraint.suspensionDownstream ? 2U : 1U)
                : 0U
        );
    }
    return hash;
}


StandingTargetSpecification::StandingTargetSpecification(
    const SupportRegion& fixedSupport,
    double requestedTopplingAngleRadians
)
    :
    support(fixedSupport),
    gammaRadians(requestedTopplingAngleRadians)
{
    validateAngle(
        gammaRadians,
        "standing toppling angle must be finite and in [0, pi/2)"
    );
}


SuspensionTargetSpecification::SuspensionTargetSpecification(
    const Eigen::Vector3d& attachmentPoint,
    const Eigen::Vector3d& gravityDirection,
    double allowedDeviationRadians
)
    :
    objective(
        makeSuspensionObjective(
            attachmentPoint,
            gravityDirection
        )
    ),
    maximumDeviationRadians(allowedDeviationRadians)
{
    validateAngle(
        maximumDeviationRadians,
        "suspension deviation must be finite and in [0, pi/2)"
    );
}


BalanceTargetSchedule::BalanceTargetSchedule(
    const BalanceTargetSpecifications& specifications
)
    :
    specifications_(specifications),
    snapshot_{0, false, false, BalanceObjectives(), {}},
    hasSnapshot_(false),
    generation_(0)
{
    if (specifications_.empty())
    {
        throw std::invalid_argument(
            "balance target schedule requires at least one specification"
        );
    }
}


BalanceTargetSnapshot BalanceTargetSchedule::refresh(
    const Eigen::Vector3d& centerOfMass
)
{
    if (!centerOfMass.allFinite())
    {
        throw std::invalid_argument(
            "balance target refresh requires a finite center of mass"
        );
    }

    BalanceTargetSnapshot next{
        generation_,
        false,
        true,
        BalanceObjectives(),
        {}
    };
    next.objectives.reserve(specifications_.size());
    next.constraints.reserve(specifications_.size());

    for (const BalanceTargetSpecification& specification : specifications_)
    {
        if (
            const StandingTargetSpecification* standing =
                std::get_if<StandingTargetSpecification>(&specification)
        )
        {
            const StandingTarget target = makeStandingTarget(
                standing->support,
                centerOfMass,
                standing->gammaRadians
            );
            const bool reached =
                target.requestedMarginAttainable() &&
                target.state() == StandingState::Stable;

            next.objectives.push_back(target.objective());
            next.constraints.push_back(BalanceConstraintSnapshot{
                BalanceConstraintType::Standing,
                target.state(),
                target.requestedMarginAttainable(),
                std::nullopt,
                std::nullopt,
                reached
            });
            next.allReached = next.allReached && reached;
            continue;
        }

        const SuspensionTargetSpecification& suspension =
            std::get<SuspensionTargetSpecification>(specification);
        const SuspensionState state = evaluateSuspension(
            centerOfMass,
            suspension
        );

        next.objectives.push_back(suspension.objective);
        next.constraints.push_back(BalanceConstraintSnapshot{
            BalanceConstraintType::Suspension,
            std::nullopt,
            std::nullopt,
            state.hasDirection
                ? std::optional<double>(state.deviationRadians)
                : std::nullopt,
            state.hasDirection
                ? std::optional<bool>(state.downstream)
                : std::nullopt,
            state.reached
        });
        next.allReached = next.allReached && state.reached;
    }

    next.contextChanged =
        !hasSnapshot_ ||
        !objectiveContextsEquivalent(
            next.objectives,
            snapshot_.objectives
        );

    if (next.contextChanged)
    {
        if (generation_ == std::numeric_limits<std::size_t>::max())
        {
            throw std::overflow_error(
                "balance target generation overflowed"
            );
        }
        ++generation_;
    }

    next.generation = generation_;
    snapshot_ = next;
    hasSnapshot_ = true;
    return snapshot_;
}


bool BalanceTargetSchedule::hasSnapshot() const
{
    return hasSnapshot_;
}


const BalanceTargetSnapshot&
BalanceTargetSchedule::currentSnapshot() const
{
    if (!hasSnapshot_)
    {
        throw std::logic_error(
            "balance target schedule has not been refreshed"
        );
    }
    return snapshot_;
}


std::size_t BalanceTargetSchedule::generation() const
{
    return generation_;
}


const BalanceTargetSpecifications&
BalanceTargetSchedule::specifications() const
{
    return specifications_;
}


} // namespace stability
