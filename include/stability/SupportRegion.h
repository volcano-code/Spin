#pragma once

#include <Eigen/Core>
#include <Eigen/StdVector>

#include "physics/BalanceObjective.h"

#include <vector>


namespace stability
{


enum class StandingState
{
    OutsideSupport,
    FeasibleButUnsafe,
    Stable
};


class StandingTarget;


/** A fixed convex support polygon expressed in its gravity-orthogonal plane. */
class SupportRegion
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    static SupportRegion fromMesh(
        const Eigen::MatrixXd& vertices,
        const Eigen::Vector3d& gravityDirection,
        double contactThreshold
    );

    const Eigen::Vector3d& gravity() const;
    const Eigen::Vector3d& basisU() const;
    const Eigen::Vector3d& basisV() const;
    double supportHeight() const;
    /** User-requested hypothetical flattening band; zero uses true contacts. */
    double contactThreshold() const;
    const std::vector<int>& contactVertexIndices() const;
    const Eigen::MatrixXd& supportPolygon2D() const;

    Eigen::Vector2d toPlaneCoordinates(
        const Eigen::Vector3d& point
    ) const;

    Eigen::Vector3d fromPlaneCoordinates(
        const Eigen::Vector2d& point
    ) const;

    bool contains2D(
        const Eigen::Vector2d& point
    ) const;

    double maximumInset() const;

    double topplingAngle(
        const Eigen::Vector3d& point
    ) const;

private:
    SupportRegion();

    Eigen::Vector3d gravity_;
    Eigen::Vector3d basisU_;
    Eigen::Vector3d basisV_;
    Eigen::Vector3d planeOrigin_;
    double supportHeight_;
    double contactThreshold_;
    std::vector<int> contactVertexIndices_;
    Eigen::MatrixXd supportPolygon2D_;
    double maximumInset_;
    double geometryTolerance_;

    friend class StandingTarget;
    friend StandingTarget makeStandingTarget(
        const SupportRegion& support,
        const Eigen::Vector3d& centerOfMass,
        double gammaRadians
    );
};


/** The safe-region projection and balance objective for one standing pose. */
class StandingTarget
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    StandingState state() const;
    bool requestedMarginAttainable() const;
    double height() const;
    double requestedInset() const;
    double maximumInset() const;
    const Eigen::MatrixXd& safeRegion2D() const;
    const Eigen::Vector2d& projectedCenterOfMass2D() const;
    const Eigen::Vector2d& target2D() const;
    const BalanceObjective& objective() const;

private:
    StandingTarget(
        StandingState state,
        bool requestedMarginAttainable,
        double height,
        double requestedInset,
        double maximumInset,
        const Eigen::MatrixXd& safeRegion2D,
        const Eigen::Vector2d& projectedCenterOfMass2D,
        const Eigen::Vector2d& target2D,
        const BalanceObjective& objective
    );

    StandingState state_;
    bool requestedMarginAttainable_;
    double height_;
    double requestedInset_;
    double maximumInset_;
    Eigen::MatrixXd safeRegion2D_;
    Eigen::Vector2d projectedCenterOfMass2D_;
    Eigen::Vector2d target2D_;
    BalanceObjective objective_;

    friend StandingTarget makeStandingTarget(
        const SupportRegion& support,
        const Eigen::Vector3d& centerOfMass,
        double gammaRadians
    );
};


StandingTarget makeStandingTarget(
    const SupportRegion& support,
    const Eigen::Vector3d& centerOfMass,
    double gammaRadians
);


} // namespace stability
