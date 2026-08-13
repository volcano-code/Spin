#pragma once

#include <Eigen/Core>
#include <Eigen/StdVector>

#include "physics/BalanceObjective.h"
#include "stability/SupportRegion.h"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <variant>
#include <vector>


namespace stability
{


enum class BalanceConstraintType
{
    Standing,
    Suspension
};


/** One standing constraint whose support remains fixed during optimization. */
struct StandingTargetSpecification
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    SupportRegion support;
    double gammaRadians;

    StandingTargetSpecification(
        const SupportRegion& fixedSupport,
        double requestedTopplingAngleRadians
    );
};


/** One suspension constraint with a fixed attachment and angular tolerance. */
struct SuspensionTargetSpecification
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    BalanceObjective objective;
    double maximumDeviationRadians;

    SuspensionTargetSpecification(
        const Eigen::Vector3d& attachmentPoint,
        const Eigen::Vector3d& gravityDirection,
        double allowedDeviationRadians
    );
};


using BalanceTargetSpecification = std::variant<
    StandingTargetSpecification,
    SuspensionTargetSpecification
>;

using BalanceTargetSpecifications = std::vector<
    BalanceTargetSpecification,
    Eigen::aligned_allocator<BalanceTargetSpecification>
>;


/** Per-constraint state accompanying one ordered objective snapshot. */
struct BalanceConstraintSnapshot
{
    BalanceConstraintType type;

    /** Set only for a standing constraint. */
    std::optional<StandingState> standingState;
    std::optional<bool> requestedMarginAttainable;

    /** Set only for a suspension constraint with non-zero attachment offset. */
    std::optional<double> suspensionDeviationRadians;
    std::optional<bool> suspensionDownstream;

    bool reached;
};


/**
 * Ordered targets and satisfaction states evaluated at one center of mass.
 *
 * contextChanged means that this snapshot's objective targets differ from the
 * previous committed snapshot.  Consequently, energies across that boundary
 * are not a same-objective descent pair.  generation increases only at such
 * an objective-context boundary; state-only changes leave it unchanged.
 */
struct BalanceTargetSnapshot
{
    std::size_t generation;
    bool contextChanged;
    bool allReached;
    BalanceObjectives objectives;
    std::vector<BalanceConstraintSnapshot> constraints;
};


/** Deterministic content hash for a complete ordered target snapshot. */
std::uint64_t balanceTargetSnapshotContentHash(
    const BalanceTargetSnapshot& snapshot
) noexcept;


/**
 * Fixed-support target schedule for the R7 alternating optimizer.
 *
 * Standing supports and suspension attachments never move.  A refresh only
 * projects the current center of mass onto each standing safe region and
 * evaluates every constraint.  Objective order is identical to specification
 * order, and terminal success uses all-of semantics.
 */
class BalanceTargetSchedule
{
public:
    explicit BalanceTargetSchedule(
        const BalanceTargetSpecifications& specifications
    );

    BalanceTargetSnapshot refresh(
        const Eigen::Vector3d& centerOfMass
    );

    bool hasSnapshot() const;
    const BalanceTargetSnapshot& currentSnapshot() const;
    std::size_t generation() const;

    const BalanceTargetSpecifications& specifications() const;

private:
    BalanceTargetSpecifications specifications_;
    BalanceTargetSnapshot snapshot_;
    bool hasSnapshot_;
    std::size_t generation_;
};


} // namespace stability
