#include "stability/SupportRegion.h"

#include "physics/BalanceEnergy.h"

#include <Eigen/Geometry>
#include <Eigen/LU>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>


namespace stability
{
namespace
{


struct HalfPlane
{
    Eigen::Vector2d inwardNormal;
    double offset;
};


double cross2D(
    const Eigen::Vector2d& first,
    const Eigen::Vector2d& second
)
{
    return first.x() * second.y() -
        first.y() * second.x();
}


double geometryScale(const Eigen::MatrixXd& vertices)
{
    const Eigen::RowVector3d minimum =
        vertices.colwise().minCoeff();

    const Eigen::RowVector3d maximum =
        vertices.colwise().maxCoeff();

    return (maximum - minimum).norm();
}


double coordinateTolerance(
    double scale
)
{
    const double epsilon =
        std::numeric_limits<double>::epsilon();

    return 256.0 * epsilon * scale;
}


bool lexicographicLess(
    const Eigen::Vector2d& first,
    const Eigen::Vector2d& second
)
{
    if (first.x() < second.x())
    {
        return true;
    }

    if (first.x() > second.x())
    {
        return false;
    }

    return first.y() < second.y();
}


std::vector<Eigen::Vector2d> uniqueSortedPoints(
    std::vector<Eigen::Vector2d> points,
    double tolerance
)
{
    std::sort(
        points.begin(),
        points.end(),
        lexicographicLess
    );

    std::vector<Eigen::Vector2d> uniquePoints;

    for (const Eigen::Vector2d& point : points)
    {
        if (
            uniquePoints.empty() ||
            (point - uniquePoints.back()).norm() > tolerance
        )
        {
            uniquePoints.push_back(point);
        }
    }

    return uniquePoints;
}


std::vector<Eigen::Vector2d> convexHull(
    const std::vector<Eigen::Vector2d>& inputPoints,
    double coordinateToleranceValue
)
{
    std::vector<Eigen::Vector2d> points =
        uniqueSortedPoints(
            inputPoints,
            coordinateToleranceValue
        );

    if (points.size() <= 2)
    {
        return points;
    }

    Eigen::Vector2d minimum = points.front();
    Eigen::Vector2d maximum = points.front();

    for (const Eigen::Vector2d& point : points)
    {
        minimum = minimum.cwiseMin(point);
        maximum = maximum.cwiseMax(point);
    }

    const double crossTolerance =
        coordinateToleranceValue *
        std::max(
            (maximum - minimum).norm(),
            coordinateToleranceValue
        );

    std::vector<Eigen::Vector2d> lower;

    for (const Eigen::Vector2d& point : points)
    {
        while (
            lower.size() >= 2 &&
            cross2D(
                lower.back() - lower[lower.size() - 2],
                point - lower.back()
            ) <= crossTolerance
        )
        {
            lower.pop_back();
        }

        lower.push_back(point);
    }

    std::vector<Eigen::Vector2d> upper;

    for (auto iterator = points.rbegin();
         iterator != points.rend();
         ++iterator)
    {
        const Eigen::Vector2d& point = *iterator;

        while (
            upper.size() >= 2 &&
            cross2D(
                upper.back() - upper[upper.size() - 2],
                point - upper.back()
            ) <= crossTolerance
        )
        {
            upper.pop_back();
        }

        upper.push_back(point);
    }

    lower.pop_back();
    upper.pop_back();
    lower.insert(lower.end(), upper.begin(), upper.end());

    return lower;
}


Eigen::MatrixXd pointsToMatrix(
    const std::vector<Eigen::Vector2d>& points
)
{
    Eigen::MatrixXd matrix(
        static_cast<Eigen::Index>(points.size()),
        2
    );

    for (std::size_t index = 0;
         index < points.size();
         ++index)
    {
        matrix.row(static_cast<Eigen::Index>(index)) =
            points[index].transpose();
    }

    return matrix;
}


std::vector<HalfPlane> polygonHalfPlanes(
    const Eigen::MatrixXd& polygon,
    double tolerance
)
{
    std::vector<HalfPlane> halfPlanes;
    halfPlanes.reserve(
        static_cast<std::size_t>(polygon.rows())
    );

    for (Eigen::Index index = 0;
         index < polygon.rows();
         ++index)
    {
        const Eigen::Vector2d first =
            polygon.row(index).transpose();

        const Eigen::Vector2d second =
            polygon.row(
                (index + 1) % polygon.rows()
            ).transpose();

        const Eigen::Vector2d edge = second - first;
        const double length = edge.norm();

        if (!(length > tolerance))
        {
            throw std::domain_error(
                "support polygon contains a degenerate edge"
            );
        }

        HalfPlane halfPlane;
        halfPlane.inwardNormal =
            Eigen::Vector2d(-edge.y(), edge.x()) /
            length;

        halfPlane.offset =
            halfPlane.inwardNormal.dot(first);

        halfPlanes.push_back(halfPlane);
    }

    return halfPlanes;
}


bool satisfiesInsetConstraints(
    const Eigen::Vector2d& point,
    double inset,
    const std::vector<HalfPlane>& halfPlanes,
    double tolerance
)
{
    for (const HalfPlane& halfPlane : halfPlanes)
    {
        if (
            halfPlane.inwardNormal.dot(point) <
            halfPlane.offset + inset - tolerance
        )
        {
            return false;
        }
    }

    return true;
}


std::pair<double, Eigen::Vector2d> computeMaximumInset(
    const std::vector<HalfPlane>& halfPlanes,
    double tolerance
)
{
    double bestInset = -std::numeric_limits<double>::infinity();
    Eigen::Vector2d bestCenter = Eigen::Vector2d::Zero();

    for (std::size_t first = 0;
         first < halfPlanes.size();
         ++first)
    {
        for (std::size_t second = first + 1;
             second < halfPlanes.size();
             ++second)
        {
            for (std::size_t third = second + 1;
                 third < halfPlanes.size();
                 ++third)
            {
                Eigen::Matrix3d system;
                Eigen::Vector3d rightHandSide;

                const std::size_t indices[3] = {
                    first,
                    second,
                    third
                };

                for (int row = 0; row < 3; ++row)
                {
                    const HalfPlane& halfPlane =
                        halfPlanes[indices[row]];

                    system(row, 0) =
                        halfPlane.inwardNormal.x();

                    system(row, 1) =
                        halfPlane.inwardNormal.y();

                    system(row, 2) = -1.0;
                    rightHandSide(row) = halfPlane.offset;
                }

                Eigen::FullPivLU<Eigen::Matrix3d> decomposition(
                    system
                );

                if (decomposition.rank() < 3)
                {
                    continue;
                }

                const Eigen::Vector3d solution =
                    decomposition.solve(rightHandSide);

                const Eigen::Vector2d center = solution.head<2>();
                const double inset = solution.z();

                if (
                    !solution.allFinite() ||
                    inset < -tolerance ||
                    !satisfiesInsetConstraints(
                        center,
                        inset,
                        halfPlanes,
                        tolerance
                    )
                )
                {
                    continue;
                }

                if (
                    inset > bestInset + tolerance ||
                    (
                        std::abs(inset - bestInset) <= tolerance &&
                        lexicographicLess(center, bestCenter)
                    )
                )
                {
                    bestInset = inset;
                    bestCenter = center;
                }
            }
        }
    }

    if (!std::isfinite(bestInset) || bestInset <= tolerance)
    {
        throw std::domain_error(
            "support polygon has no positive-area inset"
        );
    }

    return std::make_pair(
        std::max(0.0, bestInset),
        bestCenter
    );
}


Eigen::MatrixXd insetRegion(
    const std::vector<HalfPlane>& halfPlanes,
    double inset,
    double tolerance,
    const Eigen::Vector2d& fallbackPoint
)
{
    std::vector<Eigen::Vector2d> intersections;

    for (std::size_t first = 0;
         first < halfPlanes.size();
         ++first)
    {
        for (std::size_t second = first + 1;
             second < halfPlanes.size();
             ++second)
        {
            const Eigen::Vector2d& firstNormal =
                halfPlanes[first].inwardNormal;

            const Eigen::Vector2d& secondNormal =
                halfPlanes[second].inwardNormal;

            const double determinant =
                cross2D(firstNormal, secondNormal);

            if (
                std::abs(determinant) <=
                128.0 * std::numeric_limits<double>::epsilon()
            )
            {
                continue;
            }

            const double firstOffset =
                halfPlanes[first].offset + inset;

            const double secondOffset =
                halfPlanes[second].offset + inset;

            Eigen::Vector2d point;
            point.x() =
                (
                    firstOffset * secondNormal.y() -
                    firstNormal.y() * secondOffset
                ) /
                determinant;

            point.y() =
                (
                    firstNormal.x() * secondOffset -
                    firstOffset * secondNormal.x()
                ) /
                determinant;

            if (
                point.allFinite() &&
                satisfiesInsetConstraints(
                    point,
                    inset,
                    halfPlanes,
                    tolerance
                )
            )
            {
                intersections.push_back(point);
            }
        }
    }

    if (
        intersections.empty() &&
        satisfiesInsetConstraints(
            fallbackPoint,
            inset,
            halfPlanes,
            tolerance
        )
    )
    {
        intersections.push_back(fallbackPoint);
    }

    const std::vector<Eigen::Vector2d> hull =
        convexHull(intersections, tolerance);

    if (hull.empty())
    {
        throw std::domain_error(
            "failed to construct the requested support inset"
        );
    }

    return pointsToMatrix(hull);
}


double pointSegmentDistance(
    const Eigen::Vector2d& point,
    const Eigen::Vector2d& first,
    const Eigen::Vector2d& second,
    Eigen::Vector2d* closestPoint
)
{
    const Eigen::Vector2d edge = second - first;
    const double squaredLength = edge.squaredNorm();

    double parameter = 0.0;

    if (squaredLength > 0.0)
    {
        parameter = std::clamp(
            (point - first).dot(edge) / squaredLength,
            0.0,
            1.0
        );
    }

    const Eigen::Vector2d closest =
        first + parameter * edge;

    if (closestPoint != nullptr)
    {
        *closestPoint = closest;
    }

    return (point - closest).norm();
}


bool containsClosedConvexSet(
    const Eigen::MatrixXd& region,
    const Eigen::Vector2d& point,
    double tolerance
)
{
    if (region.rows() == 1)
    {
        return (
            point - region.row(0).transpose()
        ).norm() <= tolerance;
    }

    if (region.rows() == 2)
    {
        return pointSegmentDistance(
            point,
            region.row(0).transpose(),
            region.row(1).transpose(),
            nullptr
        ) <= tolerance;
    }

    for (Eigen::Index index = 0;
         index < region.rows();
         ++index)
    {
        const Eigen::Vector2d first =
            region.row(index).transpose();

        const Eigen::Vector2d second =
            region.row(
                (index + 1) % region.rows()
            ).transpose();

        const Eigen::Vector2d edge = second - first;
        const double scaledTolerance =
            tolerance * std::max(edge.norm(), tolerance);

        if (
            cross2D(edge, point - first) <
            -scaledTolerance
        )
        {
            return false;
        }
    }

    return true;
}


Eigen::Vector2d closestPointInClosedConvexSet(
    const Eigen::MatrixXd& region,
    const Eigen::Vector2d& point,
    double tolerance
)
{
    if (containsClosedConvexSet(region, point, tolerance))
    {
        return point;
    }

    if (region.rows() == 1)
    {
        return region.row(0).transpose();
    }

    double minimumDistance =
        std::numeric_limits<double>::infinity();

    Eigen::Vector2d result = Eigen::Vector2d::Zero();
    const Eigen::Index edgeCount =
        region.rows() == 2 ? 1 : region.rows();

    for (Eigen::Index index = 0;
         index < edgeCount;
         ++index)
    {
        const Eigen::Vector2d first =
            region.row(index).transpose();

        const Eigen::Vector2d second =
            region.row(
                (index + 1) % region.rows()
            ).transpose();

        Eigen::Vector2d candidate;
        const double distance = pointSegmentDistance(
            point,
            first,
            second,
            &candidate
        );

        if (
            distance < minimumDistance - tolerance ||
            (
                std::abs(distance - minimumDistance) <= tolerance &&
                lexicographicLess(candidate, result)
            )
        )
        {
            minimumDistance = distance;
            result = candidate;
        }
    }

    return result;
}


double distanceToPolygonBoundary(
    const Eigen::MatrixXd& polygon,
    const Eigen::Vector2d& point
)
{
    double minimumDistance =
        std::numeric_limits<double>::infinity();

    for (Eigen::Index index = 0;
         index < polygon.rows();
         ++index)
    {
        minimumDistance = std::min(
            minimumDistance,
            pointSegmentDistance(
                point,
                polygon.row(index).transpose(),
                polygon.row(
                    (index + 1) % polygon.rows()
                ).transpose(),
                nullptr
            )
        );
    }

    return minimumDistance;
}


} // namespace


SupportRegion::SupportRegion()
    :
    gravity_(Eigen::Vector3d::Zero()),
    basisU_(Eigen::Vector3d::Zero()),
    basisV_(Eigen::Vector3d::Zero()),
    planeOrigin_(Eigen::Vector3d::Zero()),
    supportHeight_(0.0),
    contactThreshold_(0.0),
    contactVertexIndices_(),
    supportPolygon2D_(0, 2),
    maximumInset_(0.0),
    geometryTolerance_(0.0)
{
}


SupportRegion SupportRegion::fromMesh(
    const Eigen::MatrixXd& vertices,
    const Eigen::Vector3d& gravityDirection,
    double contactThreshold
)
{
    if (
        vertices.cols() != 3 ||
        vertices.rows() == 0 ||
        !vertices.allFinite()
    )
    {
        throw std::invalid_argument(
            "support vertices must be a non-empty finite N-by-3 matrix"
        );
    }

    if (
        !std::isfinite(contactThreshold) ||
        contactThreshold < 0.0
    )
    {
        throw std::invalid_argument(
            "contact threshold must be finite and non-negative"
        );
    }

    SupportRegion result;
    result.gravity_ = normalizedGravityDirection(
        gravityDirection
    );

    const double scale = geometryScale(vertices);

    if (!std::isfinite(scale) || !(scale > 0.0))
    {
        throw std::domain_error(
            "support vertices have zero geometric extent"
        );
    }

    result.geometryTolerance_ =
        coordinateTolerance(scale);

    Eigen::Index leastAlignedAxis = 0;
    result.gravity_.cwiseAbs().minCoeff(
        &leastAlignedAxis
    );

    Eigen::Vector3d referenceAxis =
        Eigen::Vector3d::Zero();

    referenceAxis(leastAlignedAxis) = 1.0;

    result.basisU_ =
        referenceAxis -
        referenceAxis.dot(result.gravity_) *
            result.gravity_;

    const double basisUNorm = result.basisU_.norm();

    if (!std::isfinite(basisUNorm) || !(basisUNorm > 0.0))
    {
        throw std::domain_error(
            "failed to construct the support-plane basis"
        );
    }

    result.basisU_ /= basisUNorm;
    result.basisV_ =
        result.gravity_.cross(result.basisU_).normalized();

    Eigen::VectorXd heights(vertices.rows());

    for (Eigen::Index index = 0;
         index < vertices.rows();
         ++index)
    {
        heights(index) =
            -result.gravity_.dot(
                vertices.row(index).transpose()
            );
    }

    result.supportHeight_ = heights.minCoeff();
    result.contactThreshold_ = contactThreshold;
    const Eigen::RowVector3d minimumVertex =
        vertices.colwise().minCoeff();
    const Eigen::RowVector3d maximumVertex =
        vertices.colwise().maxCoeff();

    result.planeOrigin_ =
        minimumVertex.transpose() +
        0.5 * (maximumVertex - minimumVertex).transpose();
    result.planeOrigin_ +=
        (
            -result.supportHeight_ -
            result.gravity_.dot(result.planeOrigin_)
        ) * result.gravity_;

    std::vector<Eigen::Vector2d> contactPoints;

    for (Eigen::Index index = 0;
         index < vertices.rows();
         ++index)
    {
        if (
            heights(index) <=
            result.supportHeight_ +
                contactThreshold +
                result.geometryTolerance_
        )
        {
            result.contactVertexIndices_.push_back(
                static_cast<int>(index)
            );

            const Eigen::Vector3d point =
                vertices.row(index).transpose() -
                result.planeOrigin_;

            contactPoints.emplace_back(
                point.dot(result.basisU_),
                point.dot(result.basisV_)
            );
        }
    }

    const std::vector<Eigen::Vector2d> hull =
        convexHull(
            contactPoints,
            result.geometryTolerance_
        );

    if (hull.size() < 3)
    {
        throw std::domain_error(
            "contact vertices do not define a support polygon"
        );
    }

    result.supportPolygon2D_ = pointsToMatrix(hull);

    const std::vector<HalfPlane> halfPlanes =
        polygonHalfPlanes(
            result.supportPolygon2D_,
            result.geometryTolerance_
        );

    result.maximumInset_ = computeMaximumInset(
        halfPlanes,
        result.geometryTolerance_
    ).first;

    return result;
}


const Eigen::Vector3d& SupportRegion::gravity() const
{
    return gravity_;
}


const Eigen::Vector3d& SupportRegion::basisU() const
{
    return basisU_;
}


const Eigen::Vector3d& SupportRegion::basisV() const
{
    return basisV_;
}


double SupportRegion::supportHeight() const
{
    return supportHeight_;
}


double SupportRegion::contactThreshold() const
{
    return contactThreshold_;
}


const std::vector<int>&
SupportRegion::contactVertexIndices() const
{
    return contactVertexIndices_;
}


const Eigen::MatrixXd& SupportRegion::supportPolygon2D() const
{
    return supportPolygon2D_;
}


Eigen::Vector2d SupportRegion::toPlaneCoordinates(
    const Eigen::Vector3d& point
) const
{
    if (!point.allFinite())
    {
        throw std::invalid_argument(
            "support-plane point must be finite"
        );
    }

    return Eigen::Vector2d(
        (point - planeOrigin_).dot(basisU_),
        (point - planeOrigin_).dot(basisV_)
    );
}


Eigen::Vector3d SupportRegion::fromPlaneCoordinates(
    const Eigen::Vector2d& point
) const
{
    if (!point.allFinite())
    {
        throw std::invalid_argument(
            "support-plane coordinates must be finite"
        );
    }

    return planeOrigin_ +
        point.x() * basisU_ +
        point.y() * basisV_;
}


bool SupportRegion::contains2D(
    const Eigen::Vector2d& point
) const
{
    if (!point.allFinite())
    {
        throw std::invalid_argument(
            "support-plane coordinates must be finite"
        );
    }

    return containsClosedConvexSet(
        supportPolygon2D_,
        point,
        geometryTolerance_
    );
}


double SupportRegion::maximumInset() const
{
    return maximumInset_;
}


double SupportRegion::topplingAngle(
    const Eigen::Vector3d& point
) const
{
    if (!point.allFinite())
    {
        throw std::invalid_argument(
            "toppling query point must be finite"
        );
    }

    const double height =
        -gravity_.dot(point) - supportHeight_;

    if (!std::isfinite(height) || !(height > 0.0))
    {
        throw std::invalid_argument(
            "toppling query must lie above the support plane"
        );
    }

    const Eigen::Vector2d projected =
        toPlaneCoordinates(point);

    if (!contains2D(projected))
    {
        return 0.0;
    }

    return std::atan2(
        distanceToPolygonBoundary(
            supportPolygon2D_,
            projected
        ),
        height
    );
}


StandingTarget::StandingTarget(
    StandingState state,
    bool requestedMarginAttainable,
    double height,
    double requestedInset,
    double maximumInset,
    const Eigen::MatrixXd& safeRegion2D,
    const Eigen::Vector2d& projectedCenterOfMass2D,
    const Eigen::Vector2d& target2D,
    const BalanceObjective& objective
)
    :
    state_(state),
    requestedMarginAttainable_(requestedMarginAttainable),
    height_(height),
    requestedInset_(requestedInset),
    maximumInset_(maximumInset),
    safeRegion2D_(safeRegion2D),
    projectedCenterOfMass2D_(projectedCenterOfMass2D),
    target2D_(target2D),
    objective_(objective)
{
}


StandingState StandingTarget::state() const
{
    return state_;
}


bool StandingTarget::requestedMarginAttainable() const
{
    return requestedMarginAttainable_;
}


double StandingTarget::height() const
{
    return height_;
}


double StandingTarget::requestedInset() const
{
    return requestedInset_;
}


double StandingTarget::maximumInset() const
{
    return maximumInset_;
}


const Eigen::MatrixXd& StandingTarget::safeRegion2D() const
{
    return safeRegion2D_;
}


const Eigen::Vector2d&
StandingTarget::projectedCenterOfMass2D() const
{
    return projectedCenterOfMass2D_;
}


const Eigen::Vector2d& StandingTarget::target2D() const
{
    return target2D_;
}


const BalanceObjective& StandingTarget::objective() const
{
    return objective_;
}


StandingTarget makeStandingTarget(
    const SupportRegion& support,
    const Eigen::Vector3d& centerOfMass,
    double gammaRadians
)
{
    constexpr double HALF_PI =
        1.57079632679489661923132169163975144;

    if (!centerOfMass.allFinite())
    {
        throw std::invalid_argument(
            "center of mass must be finite"
        );
    }

    if (
        !std::isfinite(gammaRadians) ||
        gammaRadians < 0.0 ||
        gammaRadians >= HALF_PI
    )
    {
        throw std::invalid_argument(
            "toppling angle must be finite and in [0, pi/2)"
        );
    }

    const double height =
        -support.gravity_.dot(centerOfMass) -
        support.supportHeight_;

    if (!std::isfinite(height) || !(height > 0.0))
    {
        throw std::invalid_argument(
            "center of mass must lie above the support plane"
        );
    }

    const double requestedInset =
        height * std::tan(gammaRadians);

    if (!std::isfinite(requestedInset))
    {
        throw std::invalid_argument(
            "requested stability inset is not finite"
        );
    }

    const std::vector<HalfPlane> halfPlanes =
        polygonHalfPlanes(
            support.supportPolygon2D_,
            support.geometryTolerance_
        );

    const std::pair<double, Eigen::Vector2d> maximum =
        computeMaximumInset(
            halfPlanes,
            support.geometryTolerance_
        );

    const bool attainable =
        requestedInset <=
        maximum.first + support.geometryTolerance_;

    const double effectiveInset = attainable ?
        std::min(requestedInset, maximum.first) :
        maximum.first;

    const Eigen::MatrixXd safeRegion = insetRegion(
        halfPlanes,
        effectiveInset,
        support.geometryTolerance_,
        maximum.second
    );

    const Eigen::Vector2d projectedCenterOfMass =
        support.toPlaneCoordinates(centerOfMass);

    const bool insideSupport =
        support.contains2D(projectedCenterOfMass);

    const bool insideSafeRegion =
        containsClosedConvexSet(
            safeRegion,
            projectedCenterOfMass,
            support.geometryTolerance_
        );

    StandingState state = StandingState::OutsideSupport;

    if (insideSupport)
    {
        state = attainable && insideSafeRegion ?
            StandingState::Stable :
            StandingState::FeasibleButUnsafe;
    }

    const Eigen::Vector2d target =
        closestPointInClosedConvexSet(
            safeRegion,
            projectedCenterOfMass,
            support.geometryTolerance_
        );

    const BalanceObjective objective(
        support.fromPlaneCoordinates(target),
        support.gravity_
    );

    return StandingTarget(
        state,
        attainable,
        height,
        requestedInset,
        maximum.first,
        safeRegion,
        projectedCenterOfMass,
        target,
        objective
    );
}


} // namespace stability
