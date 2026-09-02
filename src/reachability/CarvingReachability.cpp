#include "reachability/CarvingReachability.h"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <unordered_set>

namespace reachability
{
namespace
{

bool finite(const Eigen::MatrixXd& value)
{
    return value.array().isFinite().all();
}

double signedDoubleArea(const Eigen::MatrixXd& polygon)
{
    double area = 0.0;
    for (Eigen::Index index = 0; index < polygon.rows(); ++index)
    {
        const Eigen::Index next = (index + 1) % polygon.rows();
        area += polygon(index, 0) * polygon(next, 1) -
            polygon(next, 0) * polygon(index, 1);
    }
    return area;
}

} // namespace

Halfspaces2D safePolygonHalfspaces(const Eigen::MatrixXd& polygon)
{
    if (polygon.cols() != 2 || polygon.rows() < 3)
    {
        throw std::invalid_argument(
            "safe polygon must contain at least three 2D vertices"
        );
    }
    if (!finite(polygon))
    {
        throw std::invalid_argument("safe polygon contains non-finite values");
    }
    if (!(signedDoubleArea(polygon) > 0.0))
    {
        throw std::invalid_argument(
            "safe polygon vertices must be strictly counter-clockwise"
        );
    }

    Halfspaces2D result;
    result.A.resize(polygon.rows(), 2);
    result.b.resize(polygon.rows());
    for (Eigen::Index index = 0; index < polygon.rows(); ++index)
    {
        const Eigen::Index next = (index + 1) % polygon.rows();
        const Eigen::Vector2d edge =
            polygon.row(next).transpose() - polygon.row(index).transpose();
        const double length = edge.norm();
        if (!(length > 0.0) || !std::isfinite(length))
        {
            throw std::invalid_argument(
                "safe polygon contains a duplicate or invalid edge"
            );
        }
        const Eigen::Vector2d outward(edge.y(), -edge.x());
        result.A.row(index) = (outward / length).transpose();
        result.b(index) = result.A.row(index).dot(polygon.row(index));
    }

    for (Eigen::Index vertex = 0; vertex < polygon.rows(); ++vertex)
    {
        const Eigen::VectorXd violations =
            result.A * polygon.row(vertex).transpose() - result.b;
        if (violations.maxCoeff() > 1e-12)
        {
            throw std::invalid_argument(
                "safe polygon must be convex and counter-clockwise"
            );
        }
    }
    return result;
}

bool contains(
    const Halfspaces2D& halfspaces,
    const Eigen::Vector2d& point,
    double tolerance
)
{
    if (
        halfspaces.A.cols() != 2 ||
        halfspaces.A.rows() != halfspaces.b.size() ||
        !finite(halfspaces.A) ||
        !halfspaces.b.array().isFinite().all() ||
        !point.array().isFinite().all() ||
        !std::isfinite(tolerance) ||
        tolerance < 0.0
    )
    {
        throw std::invalid_argument("invalid halfspace containment query");
    }
    if (halfspaces.A.rows() == 0)
    {
        return false;
    }
    return (
        halfspaces.A * point - halfspaces.b
    ).maxCoeff() <= tolerance;
}

void validateInput(const CarvingReachabilityInput& input)
{
    if (!(input.fixedMass > 0.0) || !std::isfinite(input.fixedMass))
    {
        throw std::invalid_argument("fixed shell mass must be positive");
    }
    if (!input.fixedFirstMoment2D.array().isFinite().all())
    {
        throw std::invalid_argument(
            "fixed shell first moment contains non-finite values"
        );
    }
    if (
        !(input.feasibilityTolerance > 0.0) ||
        !std::isfinite(input.feasibilityTolerance)
    )
    {
        throw std::invalid_argument("feasibility tolerance must be positive");
    }
    static_cast<void>(safePolygonHalfspaces(input.safePolygon2D));

    std::unordered_set<int> identifiers;
    for (const ReachabilityCell& cell : input.variableCells)
    {
        if (
            cell.cellIndex < 0 ||
            cell.mass < 0.0 ||
            !std::isfinite(cell.mass) ||
            !cell.firstMoment2D.array().isFinite().all()
        )
        {
            throw std::invalid_argument("invalid variable cell moments");
        }
        if (!identifiers.insert(cell.cellIndex).second)
        {
            throw std::invalid_argument("duplicate variable cell identifier");
        }
    }
}

const char* toString(ReachabilityStatus status) noexcept
{
    switch (status)
    {
    case ReachabilityStatus::AlreadySafe: return "ALREADY_SAFE";
    case ReachabilityStatus::ProvenUnreachable: return "PROVEN_UNREACHABLE";
    case ReachabilityStatus::RelaxationReachable:
        return "RELAXATION_REACHABLE";
    case ReachabilityStatus::Indeterminate: return "INDETERMINATE";
    }
    return "INDETERMINATE";
}

} // namespace reachability
