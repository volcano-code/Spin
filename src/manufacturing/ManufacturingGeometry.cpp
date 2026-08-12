#include "manufacturing/ManufacturingGeometry.h"

#include "deformation/DeformedVoxelMassField.h"
#include "deformation/GridDeformationDomain.h"
#include "extraction/VoidMeshExtractor.h"
#include "physics/MassProperties.h"
#include "voxel/VoxelGrid.h"

#include <Eigen/Geometry>

#include <igl/facet_components.h>
#include <igl/is_vertex_manifold.h>
#include <igl/winding_number.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <cerrno>
#include <fcntl.h>
#include <unistd.h>


namespace manufacturing
{
namespace
{


struct EdgeUse
{
    int forward = 0;
    int reverse = 0;
};


struct Bounds
{
    Eigen::Vector3d minimum = Eigen::Vector3d::Zero();
    Eigen::Vector3d maximum = Eigen::Vector3d::Zero();
};


struct Triangle
{
    std::array<Eigen::Vector3d, 3> points;
    std::array<int, 3> vertices;
    Bounds bounds;
    Eigen::Vector3d centroid = Eigen::Vector3d::Zero();
    int face = -1;
};


struct BvhNode
{
    Bounds bounds;
    int begin = 0;
    int end = 0;
    int left = -1;
    int right = -1;

    bool leaf() const noexcept
    {
        return left < 0;
    }
};


constexpr int OUTWARD_CELL_TRIANGLES[12][3] =
{
    {2, 1, 0}, {3, 1, 2},
    {4, 5, 6}, {6, 5, 7},
    {0, 1, 4}, {4, 1, 5},
    {6, 3, 2}, {7, 3, 6},
    {4, 2, 0}, {6, 2, 4},
    {1, 3, 5}, {5, 3, 7}
};


double coordinateScale(const Eigen::MatrixXd& vertices)
{
    if (vertices.rows() == 0 || vertices.cols() != 3 || !vertices.allFinite())
    {
        return 1.0;
    }

    const Eigen::Vector3d minimum =
        vertices.colwise().minCoeff().transpose();
    const Eigen::Vector3d maximum =
        vertices.colwise().maxCoeff().transpose();

    return std::max(
        {
            1.0,
            vertices.cwiseAbs().maxCoeff(),
            (maximum - minimum).norm()
        }
    );
}


double resolvedTolerance(
    const Eigen::MatrixXd& vertices,
    double requested
)
{
    if (!std::isfinite(requested) || requested < 0.0)
    {
        throw std::invalid_argument(
            "manufacturing geometric tolerance must be finite and nonnegative"
        );
    }

    if (requested > 0.0)
    {
        return requested;
    }

    return 2048.0 * std::numeric_limits<double>::epsilon() *
        coordinateScale(vertices);
}


Bounds triangleBounds(const std::array<Eigen::Vector3d, 3>& points)
{
    Bounds result;
    result.minimum = points[0];
    result.maximum = points[0];

    for (int index = 1; index < 3; ++index)
    {
        result.minimum = result.minimum.cwiseMin(points[index]);
        result.maximum = result.maximum.cwiseMax(points[index]);
    }

    return result;
}


Bounds mergeBounds(const Bounds& first, const Bounds& second)
{
    Bounds result;
    result.minimum = first.minimum.cwiseMin(second.minimum);
    result.maximum = first.maximum.cwiseMax(second.maximum);
    return result;
}


bool boundsOverlap(
    const Bounds& first,
    const Bounds& second,
    double tolerance
)
{
    return
        (first.minimum.array() <= second.maximum.array() + tolerance).all() &&
        (second.minimum.array() <= first.maximum.array() + tolerance).all();
}


double squaredBoundsDistance(const Bounds& first, const Bounds& second)
{
    double result = 0.0;
    for (int axis = 0; axis < 3; ++axis)
    {
        const double gap = std::max(
            {
                0.0,
                first.minimum(axis) - second.maximum(axis),
                second.minimum(axis) - first.maximum(axis)
            }
        );
        result += gap * gap;
    }
    return result;
}


std::vector<Triangle> makeTriangles(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces
)
{
    std::vector<Triangle> result;
    result.reserve(static_cast<std::size_t>(faces.rows()));

    for (Eigen::Index row = 0; row < faces.rows(); ++row)
    {
        Triangle triangle;
        triangle.face = static_cast<int>(row);

        for (int corner = 0; corner < 3; ++corner)
        {
            triangle.vertices[static_cast<std::size_t>(corner)] =
                faces(row, corner);
            triangle.points[static_cast<std::size_t>(corner)] =
                vertices.row(faces(row, corner)).transpose();
            triangle.centroid +=
                triangle.points[static_cast<std::size_t>(corner)] / 3.0;
        }

        triangle.bounds = triangleBounds(triangle.points);
        result.push_back(std::move(triangle));
    }

    return result;
}


class TriangleBvh
{
public:
    explicit TriangleBvh(const std::vector<Triangle>& triangles)
        : triangles_(triangles), order_(triangles.size()), nodes_()
    {
        std::iota(order_.begin(), order_.end(), 0);
        if (!triangles_.empty())
        {
            build(0, static_cast<int>(order_.size()));
        }
    }

    const std::vector<Triangle>& triangles() const noexcept
    {
        return triangles_;
    }

    const std::vector<int>& order() const noexcept
    {
        return order_;
    }

    const std::vector<BvhNode>& nodes() const noexcept
    {
        return nodes_;
    }

    bool empty() const noexcept
    {
        return nodes_.empty();
    }

private:
    int build(int begin, int end)
    {
        const int nodeIndex = static_cast<int>(nodes_.size());
        nodes_.push_back(BvhNode());

        Bounds bounds = triangles_[static_cast<std::size_t>(order_[begin])].bounds;
        Bounds centroidBounds;
        centroidBounds.minimum =
            triangles_[static_cast<std::size_t>(order_[begin])].centroid;
        centroidBounds.maximum = centroidBounds.minimum;

        for (int index = begin + 1; index < end; ++index)
        {
            const Triangle& triangle =
                triangles_[static_cast<std::size_t>(order_[index])];
            bounds = mergeBounds(bounds, triangle.bounds);
            centroidBounds.minimum =
                centroidBounds.minimum.cwiseMin(triangle.centroid);
            centroidBounds.maximum =
                centroidBounds.maximum.cwiseMax(triangle.centroid);
        }

        nodes_[static_cast<std::size_t>(nodeIndex)].bounds = bounds;
        nodes_[static_cast<std::size_t>(nodeIndex)].begin = begin;
        nodes_[static_cast<std::size_t>(nodeIndex)].end = end;

        if (end - begin <= 4)
        {
            return nodeIndex;
        }

        Eigen::Index axis = 0;
        (centroidBounds.maximum - centroidBounds.minimum).maxCoeff(&axis);
        const int middle = begin + (end - begin) / 2;

        std::nth_element(
            order_.begin() + begin,
            order_.begin() + middle,
            order_.begin() + end,
            [this, axis](int first, int second)
            {
                const double a = triangles_[static_cast<std::size_t>(first)]
                    .centroid(axis);
                const double b = triangles_[static_cast<std::size_t>(second)]
                    .centroid(axis);
                return a == b ? first < second : a < b;
            }
        );

        const int left = build(begin, middle);
        const int right = build(middle, end);
        nodes_[static_cast<std::size_t>(nodeIndex)].left = left;
        nodes_[static_cast<std::size_t>(nodeIndex)].right = right;
        return nodeIndex;
    }

    std::vector<Triangle> triangles_;
    std::vector<int> order_;
    std::vector<BvhNode> nodes_;
};


double orient2d(
    const Eigen::Vector2d& a,
    const Eigen::Vector2d& b,
    const Eigen::Vector2d& c
)
{
    return
        (b(0) - a(0)) * (c(1) - a(1)) -
        (b(1) - a(1)) * (c(0) - a(0));
}


bool pointOnSegment2d(
    const Eigen::Vector2d& point,
    const Eigen::Vector2d& first,
    const Eigen::Vector2d& second,
    double tolerance
)
{
    if (std::abs(orient2d(first, second, point)) > tolerance)
    {
        return false;
    }

    return
        point(0) >= std::min(first(0), second(0)) - tolerance &&
        point(0) <= std::max(first(0), second(0)) + tolerance &&
        point(1) >= std::min(first(1), second(1)) - tolerance &&
        point(1) <= std::max(first(1), second(1)) + tolerance;
}


bool segmentsIntersect2d(
    const Eigen::Vector2d& a,
    const Eigen::Vector2d& b,
    const Eigen::Vector2d& c,
    const Eigen::Vector2d& d,
    double tolerance
)
{
    const double abc = orient2d(a, b, c);
    const double abd = orient2d(a, b, d);
    const double cda = orient2d(c, d, a);
    const double cdb = orient2d(c, d, b);

    if (
        ((abc > tolerance && abd < -tolerance) ||
         (abc < -tolerance && abd > tolerance)) &&
        ((cda > tolerance && cdb < -tolerance) ||
         (cda < -tolerance && cdb > tolerance))
    )
    {
        return true;
    }

    return
        pointOnSegment2d(c, a, b, tolerance) ||
        pointOnSegment2d(d, a, b, tolerance) ||
        pointOnSegment2d(a, c, d, tolerance) ||
        pointOnSegment2d(b, c, d, tolerance);
}


bool pointInTriangle2d(
    const Eigen::Vector2d& point,
    const std::array<Eigen::Vector2d, 3>& triangle,
    double tolerance
)
{
    const double first = orient2d(triangle[0], triangle[1], point);
    const double second = orient2d(triangle[1], triangle[2], point);
    const double third = orient2d(triangle[2], triangle[0], point);
    const bool negative =
        first < -tolerance || second < -tolerance || third < -tolerance;
    const bool positive =
        first > tolerance || second > tolerance || third > tolerance;
    return !(negative && positive);
}


Eigen::Vector2d projectedPoint(const Eigen::Vector3d& point, int droppedAxis)
{
    if (droppedAxis == 0)
    {
        return Eigen::Vector2d(point(1), point(2));
    }
    if (droppedAxis == 1)
    {
        return Eigen::Vector2d(point(0), point(2));
    }
    return Eigen::Vector2d(point(0), point(1));
}


bool coplanarTrianglesIntersect(
    const Triangle& first,
    const Triangle& second,
    const Eigen::Vector3d& normal,
    double tolerance
)
{
    Eigen::Index droppedAxis = 0;
    normal.cwiseAbs().maxCoeff(&droppedAxis);

    std::array<Eigen::Vector2d, 3> a;
    std::array<Eigen::Vector2d, 3> b;
    for (int index = 0; index < 3; ++index)
    {
        a[static_cast<std::size_t>(index)] = projectedPoint(
            first.points[static_cast<std::size_t>(index)],
            static_cast<int>(droppedAxis)
        );
        b[static_cast<std::size_t>(index)] = projectedPoint(
            second.points[static_cast<std::size_t>(index)],
            static_cast<int>(droppedAxis)
        );
    }

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            if (segmentsIntersect2d(
                    a[static_cast<std::size_t>(i)],
                    a[static_cast<std::size_t>((i + 1) % 3)],
                    b[static_cast<std::size_t>(j)],
                    b[static_cast<std::size_t>((j + 1) % 3)],
                    tolerance
                ))
            {
                return true;
            }
        }
    }

    return
        pointInTriangle2d(a[0], b, tolerance) ||
        pointInTriangle2d(b[0], a, tolerance);
}


bool segmentIntersectsTriangle(
    const Eigen::Vector3d& start,
    const Eigen::Vector3d& end,
    const Triangle& triangle,
    double tolerance
)
{
    const Eigen::Vector3d direction = end - start;
    const Eigen::Vector3d edge1 = triangle.points[1] - triangle.points[0];
    const Eigen::Vector3d edge2 = triangle.points[2] - triangle.points[0];
    const Eigen::Vector3d p = direction.cross(edge2);
    const double determinant = edge1.dot(p);
    const double scale = std::max(
        1.0,
        direction.norm() * edge1.norm() * edge2.norm()
    );

    if (std::abs(determinant) <= tolerance * scale)
    {
        return false;
    }

    const double inverse = 1.0 / determinant;
    const Eigen::Vector3d offset = start - triangle.points[0];
    const double u = offset.dot(p) * inverse;
    const Eigen::Vector3d q = offset.cross(edge1);
    const double v = direction.dot(q) * inverse;
    const double parameter = edge2.dot(q) * inverse;

    return
        u >= -tolerance &&
        v >= -tolerance &&
        u + v <= 1.0 + tolerance &&
        parameter >= -tolerance &&
        parameter <= 1.0 + tolerance;
}


bool trianglesIntersect(
    const Triangle& first,
    const Triangle& second,
    double tolerance
)
{
    if (!boundsOverlap(first.bounds, second.bounds, tolerance))
    {
        return false;
    }

    const Eigen::Vector3d firstNormal =
        (first.points[1] - first.points[0]).cross(
            first.points[2] - first.points[0]
        );
    const Eigen::Vector3d secondNormal =
        (second.points[1] - second.points[0]).cross(
            second.points[2] - second.points[0]
        );
    const double firstNorm = firstNormal.norm();
    const double secondNorm = secondNormal.norm();

    std::array<double, 3> secondToFirst{};
    std::array<double, 3> firstToSecond{};
    for (int index = 0; index < 3; ++index)
    {
        secondToFirst[static_cast<std::size_t>(index)] =
            firstNormal.dot(
                second.points[static_cast<std::size_t>(index)] -
                first.points[0]
            ) / firstNorm;
        firstToSecond[static_cast<std::size_t>(index)] =
            secondNormal.dot(
                first.points[static_cast<std::size_t>(index)] -
                second.points[0]
            ) / secondNorm;
    }

    const auto strictlyOneSide = [tolerance](const std::array<double, 3>& d)
    {
        return
            (d[0] > tolerance && d[1] > tolerance && d[2] > tolerance) ||
            (d[0] < -tolerance && d[1] < -tolerance && d[2] < -tolerance);
    };

    if (strictlyOneSide(secondToFirst) || strictlyOneSide(firstToSecond))
    {
        return false;
    }

    const bool coplanar =
        std::all_of(
            secondToFirst.begin(),
            secondToFirst.end(),
            [tolerance](double value)
            {
                return std::abs(value) <= tolerance;
            }
        ) &&
        std::all_of(
            firstToSecond.begin(),
            firstToSecond.end(),
            [tolerance](double value)
            {
                return std::abs(value) <= tolerance;
            }
        );

    if (coplanar)
    {
        return coplanarTrianglesIntersect(
            first,
            second,
            firstNormal,
            tolerance
        );
    }

    for (int edge = 0; edge < 3; ++edge)
    {
        if (
            segmentIntersectsTriangle(
                first.points[static_cast<std::size_t>(edge)],
                first.points[static_cast<std::size_t>((edge + 1) % 3)],
                second,
                tolerance
            ) ||
            segmentIntersectsTriangle(
                second.points[static_cast<std::size_t>(edge)],
                second.points[static_cast<std::size_t>((edge + 1) % 3)],
                first,
                tolerance
            )
        )
        {
            return true;
        }
    }

    return false;
}


double squaredPointTriangleDistance(
    const Eigen::Vector3d& point,
    const Triangle& triangle
)
{
    const Eigen::Vector3d a = triangle.points[0];
    const Eigen::Vector3d ab = triangle.points[1] - a;
    const Eigen::Vector3d ac = triangle.points[2] - a;
    const Eigen::Vector3d ap = point - a;
    const double d1 = ab.dot(ap);
    const double d2 = ac.dot(ap);
    if (d1 <= 0.0 && d2 <= 0.0)
    {
        return ap.squaredNorm();
    }

    const Eigen::Vector3d b = triangle.points[1];
    const Eigen::Vector3d bp = point - b;
    const double d3 = ab.dot(bp);
    const double d4 = ac.dot(bp);
    if (d3 >= 0.0 && d4 <= d3)
    {
        return bp.squaredNorm();
    }

    const double vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0)
    {
        const double v = d1 / (d1 - d3);
        return (point - (a + v * ab)).squaredNorm();
    }

    const Eigen::Vector3d c = triangle.points[2];
    const Eigen::Vector3d cp = point - c;
    const double d5 = ab.dot(cp);
    const double d6 = ac.dot(cp);
    if (d6 >= 0.0 && d5 <= d6)
    {
        return cp.squaredNorm();
    }

    const double vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0)
    {
        const double w = d2 / (d2 - d6);
        return (point - (a + w * ac)).squaredNorm();
    }

    const double va = d3 * d6 - d5 * d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0)
    {
        const Eigen::Vector3d bc = c - b;
        const double w = (d4 - d3) /
            ((d4 - d3) + (d5 - d6));
        return (point - (b + w * bc)).squaredNorm();
    }

    const double denominator = 1.0 / (va + vb + vc);
    const double v = vb * denominator;
    const double w = vc * denominator;
    return (point - (a + ab * v + ac * w)).squaredNorm();
}


double squaredSegmentDistance(
    const Eigen::Vector3d& firstStart,
    const Eigen::Vector3d& firstEnd,
    const Eigen::Vector3d& secondStart,
    const Eigen::Vector3d& secondEnd
)
{
    const Eigen::Vector3d first = firstEnd - firstStart;
    const Eigen::Vector3d second = secondEnd - secondStart;
    const Eigen::Vector3d offset = firstStart - secondStart;
    const double a = first.dot(first);
    const double e = second.dot(second);
    const double f = second.dot(offset);
    double s = 0.0;
    double t = 0.0;

    if (a <= std::numeric_limits<double>::epsilon())
    {
        t = std::clamp(f / e, 0.0, 1.0);
    }
    else
    {
        const double c = first.dot(offset);
        if (e <= std::numeric_limits<double>::epsilon())
        {
            s = std::clamp(-c / a, 0.0, 1.0);
        }
        else
        {
            const double b = first.dot(second);
            const double denominator = a * e - b * b;
            if (denominator != 0.0)
            {
                s = std::clamp((b * f - c * e) / denominator, 0.0, 1.0);
            }
            t = (b * s + f) / e;

            if (t < 0.0)
            {
                t = 0.0;
                s = std::clamp(-c / a, 0.0, 1.0);
            }
            else if (t > 1.0)
            {
                t = 1.0;
                s = std::clamp((b - c) / a, 0.0, 1.0);
            }
        }
    }

    return (
        (firstStart + s * first) -
        (secondStart + t * second)
    ).squaredNorm();
}


double squaredTriangleDistance(
    const Triangle& first,
    const Triangle& second,
    double tolerance
)
{
    if (trianglesIntersect(first, second, tolerance))
    {
        return 0.0;
    }

    double result = std::numeric_limits<double>::infinity();
    for (int index = 0; index < 3; ++index)
    {
        result = std::min(
            result,
            squaredPointTriangleDistance(
                first.points[static_cast<std::size_t>(index)],
                second
            )
        );
        result = std::min(
            result,
            squaredPointTriangleDistance(
                second.points[static_cast<std::size_t>(index)],
                first
            )
        );
    }

    for (int firstEdge = 0; firstEdge < 3; ++firstEdge)
    {
        for (int secondEdge = 0; secondEdge < 3; ++secondEdge)
        {
            result = std::min(
                result,
                squaredSegmentDistance(
                    first.points[static_cast<std::size_t>(firstEdge)],
                    first.points[static_cast<std::size_t>((firstEdge + 1) % 3)],
                    second.points[static_cast<std::size_t>(secondEdge)],
                    second.points[static_cast<std::size_t>((secondEdge + 1) % 3)]
                )
            );
        }
    }

    return result;
}


bool shareIndexedVertex(const Triangle& first, const Triangle& second)
{
    for (const int firstVertex : first.vertices)
    {
        for (const int secondVertex : second.vertices)
        {
            if (firstVertex == secondVertex)
            {
                return true;
            }
        }
    }
    return false;
}


bool trianglesHaveDisallowedSelfIntersection(
    const Triangle& first,
    const Triangle& second,
    double tolerance
)
{
    if (!shareIndexedVertex(first, second))
    {
        return trianglesIntersect(first, second, tolerance);
    }

    const double toleranceSquared = tolerance * tolerance;
    for (int firstVertex = 0; firstVertex < 3; ++firstVertex)
    {
        const int index = first.vertices[static_cast<std::size_t>(firstVertex)];
        if (std::find(
                second.vertices.begin(),
                second.vertices.end(),
                index
            ) == second.vertices.end() &&
            squaredPointTriangleDistance(
                first.points[static_cast<std::size_t>(firstVertex)],
                second
            ) <= toleranceSquared)
        {
            return true;
        }
    }
    for (int secondVertex = 0; secondVertex < 3; ++secondVertex)
    {
        const int index = second.vertices[static_cast<std::size_t>(secondVertex)];
        if (std::find(
                first.vertices.begin(),
                first.vertices.end(),
                index
            ) == first.vertices.end() &&
            squaredPointTriangleDistance(
                second.points[static_cast<std::size_t>(secondVertex)],
                first
            ) <= toleranceSquared)
        {
            return true;
        }
    }

    for (int firstEdge = 0; firstEdge < 3; ++firstEdge)
    {
        const int firstA =
            first.vertices[static_cast<std::size_t>(firstEdge)];
        const int firstB =
            first.vertices[static_cast<std::size_t>((firstEdge + 1) % 3)];
        for (int secondEdge = 0; secondEdge < 3; ++secondEdge)
        {
            const int secondA =
                second.vertices[static_cast<std::size_t>(secondEdge)];
            const int secondB =
                second.vertices[static_cast<std::size_t>((secondEdge + 1) % 3)];
            if (
                firstA == secondA || firstA == secondB ||
                firstB == secondA || firstB == secondB
            )
            {
                continue;
            }

            if (squaredSegmentDistance(
                    first.points[static_cast<std::size_t>(firstEdge)],
                    first.points[static_cast<std::size_t>((firstEdge + 1) % 3)],
                    second.points[static_cast<std::size_t>(secondEdge)],
                    second.points[static_cast<std::size_t>((secondEdge + 1) % 3)]
                ) <= toleranceSquared)
            {
                return true;
            }
        }
    }

    return false;
}


void testSelfLeafPair(
    const Triangle& first,
    const Triangle& second,
    double tolerance,
    IntersectionAudit& audit
)
{
    if (
        !boundsOverlap(first.bounds, second.bounds, tolerance)
    )
    {
        return;
    }

    ++audit.broadPhasePairCount;
    ++audit.narrowPhasePairCount;
    if (trianglesHaveDisallowedSelfIntersection(
            first,
            second,
            tolerance
        ))
    {
        if (audit.intersectingPairCount == 0)
        {
            audit.firstFace = first.face;
            audit.secondFace = second.face;
        }
        ++audit.intersectingPairCount;
    }
}


void auditCrossNodes(
    const TriangleBvh& firstBvh,
    int firstNodeIndex,
    const TriangleBvh& secondBvh,
    int secondNodeIndex,
    double tolerance,
    IntersectionAudit& audit
)
{
    const BvhNode& firstNode =
        firstBvh.nodes()[static_cast<std::size_t>(firstNodeIndex)];
    const BvhNode& secondNode =
        secondBvh.nodes()[static_cast<std::size_t>(secondNodeIndex)];

    if (!boundsOverlap(firstNode.bounds, secondNode.bounds, tolerance))
    {
        return;
    }

    if (firstNode.leaf() && secondNode.leaf())
    {
        for (int first = firstNode.begin; first < firstNode.end; ++first)
        {
            const Triangle& firstTriangle = firstBvh.triangles()[
                static_cast<std::size_t>(firstBvh.order()[first])
            ];
            for (int second = secondNode.begin; second < secondNode.end; ++second)
            {
                const Triangle& secondTriangle = secondBvh.triangles()[
                    static_cast<std::size_t>(secondBvh.order()[second])
                ];
                if (!boundsOverlap(
                        firstTriangle.bounds,
                        secondTriangle.bounds,
                        tolerance
                    ))
                {
                    continue;
                }

                ++audit.broadPhasePairCount;
                ++audit.narrowPhasePairCount;
                if (trianglesIntersect(
                        firstTriangle,
                        secondTriangle,
                        tolerance
                    ))
                {
                    if (audit.intersectingPairCount == 0)
                    {
                        audit.firstFace = firstTriangle.face;
                        audit.secondFace = secondTriangle.face;
                    }
                    ++audit.intersectingPairCount;
                }
            }
        }
        return;
    }

    if (
        secondNode.leaf() ||
        (!firstNode.leaf() &&
         firstNode.end - firstNode.begin >= secondNode.end - secondNode.begin)
    )
    {
        auditCrossNodes(
            firstBvh,
            firstNode.left,
            secondBvh,
            secondNodeIndex,
            tolerance,
            audit
        );
        auditCrossNodes(
            firstBvh,
            firstNode.right,
            secondBvh,
            secondNodeIndex,
            tolerance,
            audit
        );
    }
    else
    {
        auditCrossNodes(
            firstBvh,
            firstNodeIndex,
            secondBvh,
            secondNode.left,
            tolerance,
            audit
        );
        auditCrossNodes(
            firstBvh,
            firstNodeIndex,
            secondBvh,
            secondNode.right,
            tolerance,
            audit
        );
    }
}


void auditSelfNode(
    const TriangleBvh& bvh,
    int nodeIndex,
    double tolerance,
    IntersectionAudit& audit
);


void auditSelfCrossNodes(
    const TriangleBvh& bvh,
    int firstNodeIndex,
    int secondNodeIndex,
    double tolerance,
    IntersectionAudit& audit
)
{
    const BvhNode& firstNode =
        bvh.nodes()[static_cast<std::size_t>(firstNodeIndex)];
    const BvhNode& secondNode =
        bvh.nodes()[static_cast<std::size_t>(secondNodeIndex)];

    if (!boundsOverlap(firstNode.bounds, secondNode.bounds, tolerance))
    {
        return;
    }

    if (firstNode.leaf() && secondNode.leaf())
    {
        for (int first = firstNode.begin; first < firstNode.end; ++first)
        {
            for (int second = secondNode.begin; second < secondNode.end; ++second)
            {
                testSelfLeafPair(
                    bvh.triangles()[static_cast<std::size_t>(bvh.order()[first])],
                    bvh.triangles()[static_cast<std::size_t>(bvh.order()[second])],
                    tolerance,
                    audit
                );
            }
        }
        return;
    }

    if (
        secondNode.leaf() ||
        (!firstNode.leaf() &&
         firstNode.end - firstNode.begin >= secondNode.end - secondNode.begin)
    )
    {
        auditSelfCrossNodes(
            bvh,
            firstNode.left,
            secondNodeIndex,
            tolerance,
            audit
        );
        auditSelfCrossNodes(
            bvh,
            firstNode.right,
            secondNodeIndex,
            tolerance,
            audit
        );
    }
    else
    {
        auditSelfCrossNodes(
            bvh,
            firstNodeIndex,
            secondNode.left,
            tolerance,
            audit
        );
        auditSelfCrossNodes(
            bvh,
            firstNodeIndex,
            secondNode.right,
            tolerance,
            audit
        );
    }
}


void auditSelfNode(
    const TriangleBvh& bvh,
    int nodeIndex,
    double tolerance,
    IntersectionAudit& audit
)
{
    const BvhNode& node = bvh.nodes()[static_cast<std::size_t>(nodeIndex)];
    if (node.leaf())
    {
        for (int first = node.begin; first < node.end; ++first)
        {
            for (int second = first + 1; second < node.end; ++second)
            {
                testSelfLeafPair(
                    bvh.triangles()[static_cast<std::size_t>(bvh.order()[first])],
                    bvh.triangles()[static_cast<std::size_t>(bvh.order()[second])],
                    tolerance,
                    audit
                );
            }
        }
        return;
    }

    auditSelfNode(bvh, node.left, tolerance, audit);
    auditSelfCrossNodes(bvh, node.left, node.right, tolerance, audit);
    auditSelfNode(bvh, node.right, tolerance, audit);
}


IntersectionAudit auditSelfIntersections(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces,
    double tolerance
)
{
    IntersectionAudit result;
    const TriangleBvh bvh(makeTriangles(vertices, faces));
    if (!bvh.empty())
    {
        auditSelfNode(bvh, 0, tolerance, result);
    }
    result.completed = true;
    result.intersectionFree = result.intersectingPairCount == 0;
    return result;
}


void closestCrossNodes(
    const TriangleBvh& outer,
    int outerNodeIndex,
    const TriangleBvh& inner,
    int innerNodeIndex,
    double tolerance,
    double& bestSquared,
    WallClearanceAudit& audit
)
{
    const BvhNode& outerNode =
        outer.nodes()[static_cast<std::size_t>(outerNodeIndex)];
    const BvhNode& innerNode =
        inner.nodes()[static_cast<std::size_t>(innerNodeIndex)];

    if (squaredBoundsDistance(outerNode.bounds, innerNode.bounds) > bestSquared)
    {
        return;
    }

    if (outerNode.leaf() && innerNode.leaf())
    {
        for (int outerIndex = outerNode.begin;
             outerIndex < outerNode.end;
             ++outerIndex)
        {
            const Triangle& outerTriangle = outer.triangles()[
                static_cast<std::size_t>(outer.order()[outerIndex])
            ];
            for (int innerIndex = innerNode.begin;
                 innerIndex < innerNode.end;
                 ++innerIndex)
            {
                const Triangle& innerTriangle = inner.triangles()[
                    static_cast<std::size_t>(inner.order()[innerIndex])
                ];
                if (squaredBoundsDistance(
                        outerTriangle.bounds,
                        innerTriangle.bounds
                    ) > bestSquared)
                {
                    continue;
                }

                ++audit.broadPhasePairCount;
                ++audit.narrowPhasePairCount;
                const double distance = squaredTriangleDistance(
                    outerTriangle,
                    innerTriangle,
                    tolerance
                );
                if (distance < bestSquared)
                {
                    bestSquared = distance;
                    audit.outerFace = outerTriangle.face;
                    audit.innerFace = innerTriangle.face;
                }
            }
        }
        return;
    }

    if (
        innerNode.leaf() ||
        (!outerNode.leaf() &&
         outerNode.end - outerNode.begin >= innerNode.end - innerNode.begin)
    )
    {
        closestCrossNodes(
            outer,
            outerNode.left,
            inner,
            innerNodeIndex,
            tolerance,
            bestSquared,
            audit
        );
        closestCrossNodes(
            outer,
            outerNode.right,
            inner,
            innerNodeIndex,
            tolerance,
            bestSquared,
            audit
        );
    }
    else
    {
        closestCrossNodes(
            outer,
            outerNodeIndex,
            inner,
            innerNode.left,
            tolerance,
            bestSquared,
            audit
        );
        closestCrossNodes(
            outer,
            outerNodeIndex,
            inner,
            innerNode.right,
            tolerance,
            bestSquared,
            audit
        );
    }
}


WallClearanceAudit auditWallClearance(
    const Eigen::MatrixXd& outerVertices,
    const Eigen::MatrixXi& outerFaces,
    const Eigen::MatrixXd& innerVertices,
    const Eigen::MatrixXi& innerFaces,
    double requestedThickness,
    double tolerance
)
{
    WallClearanceAudit result;
    result.requestedThickness = requestedThickness;
    result.completed = true;

    if (innerFaces.rows() == 0)
    {
        result.intersectionFree = true;
        result.meetsRequestedThickness = true;
        return result;
    }

    result.applicable = true;
    const TriangleBvh outer(makeTriangles(outerVertices, outerFaces));
    const TriangleBvh inner(makeTriangles(innerVertices, innerFaces));
    double bestSquared = std::numeric_limits<double>::infinity();
    closestCrossNodes(
        outer,
        0,
        inner,
        0,
        tolerance,
        bestSquared,
        result
    );

    result.measuredMinimumClearance = std::sqrt(std::max(0.0, bestSquared));
    result.intersectionFree = result.measuredMinimumClearance > tolerance;
    result.meetsRequestedThickness =
        result.intersectionFree &&
        result.measuredMinimumClearance + tolerance >= requestedThickness;
    return result;
}


ContainmentAudit auditContainment(
    const Eigen::MatrixXd& outerVertices,
    const Eigen::MatrixXi& outerFaces,
    const Eigen::MatrixXd& innerVertices,
    const Eigen::MatrixXi& innerFaces,
    const BoundaryAudit& innerAudit,
    bool shellsIntersectionFree
)
{
    ContainmentAudit result;
    result.completed = true;

    if (innerFaces.rows() == 0)
    {
        result.allCavitiesInsideOuter = true;
        return result;
    }

    result.applicable = true;
    if (!shellsIntersectionFree || innerAudit.faceComponentCount <= 0)
    {
        return result;
    }

    Eigen::VectorXi components;
    const int componentCount = igl::facet_components(innerFaces, components);
    result.testedComponentCount = componentCount;
    result.allCavitiesInsideOuter = true;

    std::vector<Eigen::MatrixXi> componentFaces(
        static_cast<std::size_t>(componentCount)
    );
    std::vector<Eigen::RowVector3d> representatives(
        static_cast<std::size_t>(componentCount),
        Eigen::RowVector3d::Zero()
    );

    for (int component = 0; component < componentCount; ++component)
    {
        const Eigen::Index count =
            (components.array() == component).count();
        componentFaces[static_cast<std::size_t>(component)].resize(count, 3);
        Eigen::Index destination = 0;
        for (Eigen::Index row = 0; row < innerFaces.rows(); ++row)
        {
            if (components(row) == component)
            {
                componentFaces[static_cast<std::size_t>(component)].row(
                    destination++
                ) = innerFaces.row(row);
            }
        }

        if (count <= 0)
        {
            result.allCavitiesInsideOuter = false;
            result.failedComponent = component;
            break;
        }

        for (int corner = 0; corner < 3; ++corner)
        {
            representatives[static_cast<std::size_t>(component)] +=
                innerVertices.row(
                    componentFaces[static_cast<std::size_t>(component)](
                        0,
                        corner
                    )
            ) / 3.0;
        }

        const double winding = igl::winding_number(
            outerVertices,
            outerFaces,
            representatives[static_cast<std::size_t>(component)]
        );

        if (!std::isfinite(winding) || winding < 0.5)
        {
            result.allCavitiesInsideOuter = false;
            result.failedComponent = component;
            result.failedWindingNumber = winding;
            break;
        }
    }

    if (result.allCavitiesInsideOuter)
    {
        for (int containing = 0; containing < componentCount; ++containing)
        {
            for (int nested = 0; nested < componentCount; ++nested)
            {
                if (containing == nested)
                {
                    continue;
                }

                const double winding = igl::winding_number(
                    innerVertices,
                    componentFaces[static_cast<std::size_t>(containing)],
                    representatives[static_cast<std::size_t>(nested)]
                );
                if (std::isfinite(winding) && std::abs(winding) >= 0.5)
                {
                    if (result.nestedCavityPairCount == 0)
                    {
                        result.firstContainingComponent = containing;
                        result.firstNestedComponent = nested;
                    }
                    ++result.nestedCavityPairCount;
                }
            }
        }
    }

    if (result.nestedCavityPairCount > 0)
    {
        result.allCavitiesInsideOuter = false;
    }

    return result;
}


Eigen::Vector3d trilinearDerivative(
    const std::array<Eigen::Vector3d, 8>& corners,
    int derivativeAxis,
    double x,
    double y,
    double z
)
{
    const double coordinates[3] = {x, y, z};
    Eigen::Vector3d result = Eigen::Vector3d::Zero();

    for (int bit = 0; bit < 8; ++bit)
    {
        double coefficient = 1.0;
        for (int axis = 0; axis < 3; ++axis)
        {
            const bool high = (bit & (1 << axis)) != 0;
            if (axis == derivativeAxis)
            {
                coefficient *= high ? 1.0 : -1.0;
            }
            else
            {
                coefficient *= high
                    ? coordinates[axis]
                    : 1.0 - coordinates[axis];
            }
        }
        result += coefficient * corners[static_cast<std::size_t>(bit)];
    }
    return result;
}


DeformedCellAudit auditDeformedCells(
    const deformation::GridDeformationDomain& domain,
    const Eigen::MatrixXd& deformedCorners,
    double tolerance
)
{
    DeformedCellAudit result;
    result.completed = true;
    result.activeCellCount = static_cast<int>(domain.activeCellCount());
    if (
        deformedCorners.rows() != domain.fullRestCorners().rows() ||
        deformedCorners.cols() != 3 ||
        !deformedCorners.allFinite()
    )
    {
        result.invalidCellCount = result.activeCellCount;
        result.firstInvalidCell = result.activeCellCount > 0 ? 0 : -1;
        return result;
    }
    const double restJacobian = std::abs(domain.voxelSize().prod());
    const double jacobianFloor = std::max(
        tolerance * tolerance * tolerance,
        4096.0 * std::numeric_limits<double>::epsilon() * restJacobian
    );
    const double areaFloorSquared =
        tolerance * tolerance *
        std::max(1.0, domain.voxelSize().squaredNorm());
    const double samples[3] = {0.0, 0.5, 1.0};

    for (int cell = 0; cell < domain.cellCount(); ++cell)
    {
        if (!domain.isActiveCell(cell))
        {
            continue;
        }

        const std::array<int, 8> global = domain.cellGlobalCorners(cell);
        std::array<Eigen::Vector3d, 8> corners{};
        for (int index = 0; index < 8; ++index)
        {
            corners[static_cast<std::size_t>(index)] =
                deformedCorners.row(
                    global[static_cast<std::size_t>(index)]
                ).transpose();
        }

        bool valid = true;
        for (const auto& face : OUTWARD_CELL_TRIANGLES)
        {
            const Eigen::Vector3d normal =
                (corners[static_cast<std::size_t>(face[1])] -
                 corners[static_cast<std::size_t>(face[0])]).cross(
                    corners[static_cast<std::size_t>(face[2])] -
                    corners[static_cast<std::size_t>(face[0])]
                );
            if (!normal.allFinite() || normal.squaredNorm() <= areaFloorSquared)
            {
                ++result.degenerateBoundaryTriangleCount;
                valid = false;
            }
        }

        for (const double x : samples)
        {
            for (const double y : samples)
            {
                for (const double z : samples)
                {
                    Eigen::Matrix3d jacobian;
                    jacobian.col(0) = trilinearDerivative(
                        corners, 0, x, y, z
                    );
                    jacobian.col(1) = trilinearDerivative(
                        corners, 1, x, y, z
                    );
                    jacobian.col(2) = trilinearDerivative(
                        corners, 2, x, y, z
                    );
                    const double determinant = jacobian.determinant();
                    result.minimumJacobianDeterminant = std::min(
                        result.minimumJacobianDeterminant,
                        determinant
                    );
                    if (!std::isfinite(determinant) ||
                        determinant <= jacobianFloor)
                    {
                        ++result.nonpositiveJacobianSampleCount;
                        valid = false;
                    }
                }
            }
        }

        if (!valid)
        {
            if (result.firstInvalidCell < 0)
            {
                result.firstInvalidCell = cell;
            }
            ++result.invalidCellCount;
        }
    }

    result.allActiveCellsValid = result.invalidCellCount == 0;
    return result;
}


void compactMesh(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces,
    Eigen::MatrixXd& compactVertices,
    Eigen::MatrixXi& compactFaces
)
{
    if (faces.rows() == 0)
    {
        compactVertices.resize(0, 3);
        compactFaces.resize(0, 3);
        return;
    }

    std::vector<int> oldToNew(
        static_cast<std::size_t>(vertices.rows()),
        -1
    );
    int count = 0;
    for (Eigen::Index index = 0; index < faces.size(); ++index)
    {
        const int oldIndex = faces(index);
        if (oldToNew[static_cast<std::size_t>(oldIndex)] < 0)
        {
            oldToNew[static_cast<std::size_t>(oldIndex)] = count++;
        }
    }

    compactVertices.resize(count, 3);
    for (Eigen::Index oldIndex = 0; oldIndex < vertices.rows(); ++oldIndex)
    {
        const int newIndex = oldToNew[static_cast<std::size_t>(oldIndex)];
        if (newIndex >= 0)
        {
            compactVertices.row(newIndex) = vertices.row(oldIndex);
        }
    }

    compactFaces.resize(faces.rows(), 3);
    for (Eigen::Index index = 0; index < faces.size(); ++index)
    {
        compactFaces(index) = oldToNew[
            static_cast<std::size_t>(faces(index))
        ];
    }
}


bool sameFaceMatrix(
    const Eigen::MatrixXi& first,
    const Eigen::MatrixXi& second
)
{
    return
        first.rows() == second.rows() &&
        first.cols() == second.cols() &&
        (first.size() == 0 || (first.array() == second.array()).all());
}


bool synchronizedInnerGeometry(
    const deformation::DeformedGeometry& geometry,
    const VoxelGrid& voxelGrid,
    double tolerance
)
{
    Eigen::MatrixXi expectedFaces;
    VoidMeshExtractor::extractFaces(voxelGrid, expectedFaces);

    return
        geometry.gridCorners.rows() == voxelGrid.corners().rows() &&
        geometry.gridCorners.cols() == 3 &&
        geometry.gridCorners.allFinite() &&
        geometry.innerVertices.rows() == geometry.gridCorners.rows() &&
        geometry.innerVertices.cols() == 3 &&
        geometry.innerVertices.allFinite() &&
        geometry.innerVertices.isApprox(geometry.gridCorners, tolerance) &&
        sameFaceMatrix(geometry.innerFaces, expectedFaces);
}


CoherentBoundaryMesh makeCoherentBoundary(
    const Eigen::MatrixXd& outerVertices,
    const Eigen::MatrixXi& outerFaces,
    const Eigen::MatrixXd& innerVertices,
    const Eigen::MatrixXi& innerFaces
)
{
    CoherentBoundaryMesh result;
    result.outerVertexCount = static_cast<int>(outerVertices.rows());
    result.outerFaceCount = static_cast<int>(outerFaces.rows());
    result.innerVertexCount = static_cast<int>(innerVertices.rows());
    result.innerFaceCount = static_cast<int>(innerFaces.rows());

    result.vertices.resize(
        outerVertices.rows() + innerVertices.rows(),
        3
    );
    if (outerVertices.rows() > 0)
    {
        result.vertices.topRows(outerVertices.rows()) = outerVertices;
    }
    if (innerVertices.rows() > 0)
    {
        result.vertices.bottomRows(innerVertices.rows()) = innerVertices;
    }

    result.faces.resize(outerFaces.rows() + innerFaces.rows(), 3);
    if (outerFaces.rows() > 0)
    {
        result.faces.topRows(outerFaces.rows()) = outerFaces;
    }
    if (innerFaces.rows() > 0)
    {
        result.faces.bottomRows(innerFaces.rows()) =
            (innerFaces.array() + result.outerVertexCount).matrix();
    }
    return result;
}


std::string escapedJson(const std::string& value)
{
    std::ostringstream output;
    for (const unsigned char character : value)
    {
        switch (character)
        {
            case '\"': output << "\\\""; break;
            case '\\': output << "\\\\"; break;
            case '\b': output << "\\b"; break;
            case '\f': output << "\\f"; break;
            case '\n': output << "\\n"; break;
            case '\r': output << "\\r"; break;
            case '\t': output << "\\t"; break;
            default:
                if (character < 0x20U)
                {
                    output
                        << "\\u"
                        << std::hex
                        << std::setw(4)
                        << std::setfill('0')
                        << static_cast<int>(character)
                        << std::dec;
                }
                else
                {
                    output << character;
                }
        }
    }
    return output.str();
}


const char* jsonBool(bool value)
{
    return value ? "true" : "false";
}


void appendJsonNumber(std::ostringstream& output, double value)
{
    if (std::isfinite(value))
    {
        output << std::setprecision(17) << value;
    }
    else
    {
        output << "null";
    }
}


void appendBoundaryJson(
    std::ostringstream& output,
    const char* name,
    const BoundaryAudit& audit
)
{
    output
        << "  \"" << name << "\": {\n"
        << "    \"empty\": " << jsonBool(audit.empty) << ",\n"
        << "    \"vertices\": " << audit.vertexCount << ",\n"
        << "    \"faces\": " << audit.faceCount << ",\n"
        << "    \"components\": " << audit.faceComponentCount << ",\n"
        << "    \"degenerate_faces\": " << audit.degenerateFaceCount << ",\n"
        << "    \"boundary_edges\": " << audit.boundaryEdgeCount << ",\n"
        << "    \"nonmanifold_edges\": " << audit.nonManifoldEdgeCount << ",\n"
        << "    \"nonmanifold_vertices\": "
        << audit.nonManifoldVertexCount << ",\n"
        << "    \"closed\": " << jsonBool(audit.closed) << ",\n"
        << "    \"consistently_oriented\": "
        << jsonBool(audit.consistentlyOriented) << ",\n"
        << "    \"component_orientations_valid\": "
        << jsonBool(audit.componentOrientationsValid) << ",\n"
        << "    \"self_intersection_pairs\": "
        << audit.selfIntersections.intersectingPairCount << ",\n"
        << "    \"signed_volume\": " << std::setprecision(17)
        << audit.signedVolume << "\n"
        << "  },\n";
}


std::string systemErrorMessage(const std::string& operation)
{
    return operation + ": " + std::strerror(errno);
}


std::vector<unsigned char> readFileBytes(const std::filesystem::path& path)
{
    std::ifstream input(path, std::ios::binary);
    if (!input)
    {
        throw std::runtime_error("failed to reopen transactional output");
    }
    return std::vector<unsigned char>(
        std::istreambuf_iterator<char>(input),
        std::istreambuf_iterator<char>()
    );
}


void writeAtomicBytes(
    const std::string& filename,
    const std::vector<unsigned char>& bytes
)
{
    if (filename.empty())
    {
        throw std::invalid_argument(
            "transactional export requires a nonempty filename"
        );
    }

    static std::atomic<unsigned long long> counter{0};
    const std::filesystem::path destination(filename);
    const std::filesystem::path directory = destination.has_parent_path()
        ? destination.parent_path()
        : std::filesystem::path(".");
    const std::string temporaryName =
        "." + destination.filename().string() +
        ".r8tmp." + std::to_string(static_cast<long long>(::getpid())) +
        "." + std::to_string(counter.fetch_add(1));
    const std::filesystem::path temporary = directory / temporaryName;

    int descriptor = ::open(
        temporary.c_str(),
        O_WRONLY | O_CREAT | O_EXCL,
        0666
    );
    if (descriptor < 0)
    {
        throw std::runtime_error(systemErrorMessage(
            "failed to create transactional temporary file"
        ));
    }

    bool closed = false;
    try
    {
        std::size_t offset = 0;
        while (offset < bytes.size())
        {
            const ssize_t written = ::write(
                descriptor,
                bytes.data() + offset,
                bytes.size() - offset
            );
            if (written < 0)
            {
                if (errno == EINTR)
                {
                    continue;
                }
                throw std::runtime_error(systemErrorMessage(
                    "failed to write transactional temporary file"
                ));
            }
            if (written == 0)
            {
                throw std::runtime_error(
                    "transactional temporary write made no progress"
                );
            }
            offset += static_cast<std::size_t>(written);
        }

        if (::fsync(descriptor) != 0)
        {
            throw std::runtime_error(systemErrorMessage(
                "failed to fsync transactional temporary file"
            ));
        }
        if (::close(descriptor) != 0)
        {
            descriptor = -1;
            throw std::runtime_error(systemErrorMessage(
                "failed to close transactional temporary file"
            ));
        }
        descriptor = -1;
        closed = true;

        if (readFileBytes(temporary) != bytes)
        {
            throw std::runtime_error(
                "transactional output failed byte-for-byte round-trip validation"
            );
        }

        if (::rename(temporary.c_str(), destination.c_str()) != 0)
        {
            throw std::runtime_error(systemErrorMessage(
                "failed to atomically replace export destination"
            ));
        }

        const int directoryDescriptor = ::open(
            directory.c_str(),
            O_RDONLY | O_DIRECTORY
        );
        if (directoryDescriptor >= 0)
        {
            const int syncResult = ::fsync(directoryDescriptor);
            const int savedError = errno;
            ::close(directoryDescriptor);
            if (syncResult != 0)
            {
                errno = savedError;
                throw std::runtime_error(systemErrorMessage(
                    "failed to fsync transactional output directory"
                ));
            }
        }
    }
    catch (...)
    {
        if (!closed && descriptor >= 0)
        {
            ::close(descriptor);
        }
        std::error_code ignored;
        std::filesystem::remove(temporary, ignored);
        throw;
    }
}


std::vector<unsigned char> stringBytes(const std::string& value)
{
    return std::vector<unsigned char>(value.begin(), value.end());
}


std::string serializeObj(const CoherentBoundaryMesh& mesh)
{
    std::ostringstream output;
    output << "# make-it-stand R8 manufacturing boundary\n";
    output << std::setprecision(17);
    for (Eigen::Index row = 0; row < mesh.vertices.rows(); ++row)
    {
        output
            << "v " << mesh.vertices(row, 0)
            << ' ' << mesh.vertices(row, 1)
            << ' ' << mesh.vertices(row, 2) << '\n';
    }
    for (Eigen::Index row = 0; row < mesh.faces.rows(); ++row)
    {
        output
            << "f " << mesh.faces(row, 0) + 1
            << ' ' << mesh.faces(row, 1) + 1
            << ' ' << mesh.faces(row, 2) + 1 << '\n';
    }
    return output.str();
}


Eigen::Vector3d triangleUnitNormal(
    const CoherentBoundaryMesh& mesh,
    Eigen::Index face
)
{
    const Eigen::Vector3d first =
        mesh.vertices.row(mesh.faces(face, 1)).transpose() -
        mesh.vertices.row(mesh.faces(face, 0)).transpose();
    const Eigen::Vector3d second =
        mesh.vertices.row(mesh.faces(face, 2)).transpose() -
        mesh.vertices.row(mesh.faces(face, 0)).transpose();
    Eigen::Vector3d normal = first.cross(second);
    const double norm = normal.norm();
    if (norm > 0.0)
    {
        normal /= norm;
    }
    else
    {
        normal.setZero();
    }
    return normal;
}


std::string serializeAsciiStl(const CoherentBoundaryMesh& mesh)
{
    std::ostringstream output;
    output << "solid make_it_stand_r8\n" << std::setprecision(17);
    for (Eigen::Index face = 0; face < mesh.faces.rows(); ++face)
    {
        const Eigen::Vector3d normal = triangleUnitNormal(mesh, face);
        output
            << "  facet normal " << normal(0) << ' '
            << normal(1) << ' ' << normal(2) << "\n"
            << "    outer loop\n";
        for (int corner = 0; corner < 3; ++corner)
        {
            const Eigen::Index vertex = mesh.faces(face, corner);
            output
                << "      vertex " << mesh.vertices(vertex, 0) << ' '
                << mesh.vertices(vertex, 1) << ' '
                << mesh.vertices(vertex, 2) << "\n";
        }
        output << "    endloop\n  endfacet\n";
    }
    output << "endsolid make_it_stand_r8\n";
    return output.str();
}


void appendLittleEndian16(
    std::vector<unsigned char>& output,
    std::uint16_t value
)
{
    output.push_back(static_cast<unsigned char>(value & 0xffU));
    output.push_back(static_cast<unsigned char>((value >> 8U) & 0xffU));
}


void appendLittleEndian32(
    std::vector<unsigned char>& output,
    std::uint32_t value
)
{
    for (int byte = 0; byte < 4; ++byte)
    {
        output.push_back(static_cast<unsigned char>(
            (value >> (8U * static_cast<unsigned int>(byte))) & 0xffU
        ));
    }
}


void appendLittleEndianFloat(
    std::vector<unsigned char>& output,
    float value
)
{
    std::uint32_t bits = 0;
    static_assert(sizeof(bits) == sizeof(value), "binary STL needs 32-bit float");
    std::memcpy(&bits, &value, sizeof(bits));
    appendLittleEndian32(output, bits);
}


std::vector<unsigned char> serializeBinaryStl(
    const CoherentBoundaryMesh& mesh
)
{
    if (mesh.faces.rows() > std::numeric_limits<std::uint32_t>::max())
    {
        throw std::overflow_error("binary STL triangle count exceeds uint32");
    }

    std::vector<unsigned char> output(80U, 0U);
    const std::string signature = "make-it-stand-r8";
    std::copy(signature.begin(), signature.end(), output.begin());
    appendLittleEndian32(
        output,
        static_cast<std::uint32_t>(mesh.faces.rows())
    );

    for (Eigen::Index face = 0; face < mesh.faces.rows(); ++face)
    {
        const Eigen::Vector3d normal = triangleUnitNormal(mesh, face);
        for (int axis = 0; axis < 3; ++axis)
        {
            appendLittleEndianFloat(output, static_cast<float>(normal(axis)));
        }
        for (int corner = 0; corner < 3; ++corner)
        {
            const Eigen::Index vertex = mesh.faces(face, corner);
            for (int axis = 0; axis < 3; ++axis)
            {
                appendLittleEndianFloat(
                    output,
                    static_cast<float>(mesh.vertices(vertex, axis))
                );
            }
        }
        appendLittleEndian16(output, 0U);
    }
    return output;
}


void validateBinaryFloatGeometry(const ManufacturingGeometry& geometry)
{
    CoherentBoundaryMesh rounded = geometry.boundary;
    double maximumError = 0.0;
    for (Eigen::Index row = 0; row < rounded.vertices.rows(); ++row)
    {
        for (int axis = 0; axis < 3; ++axis)
        {
            const double value = static_cast<double>(
                static_cast<float>(rounded.vertices(row, axis))
            );
            maximumError = std::max(
                maximumError,
                std::abs(value - rounded.vertices(row, axis))
            );
            rounded.vertices(row, axis) = value;
        }
    }

    const BoundaryAudit roundedAudit = auditTriangleBoundary(
        rounded.vertices,
        rounded.faces,
        BoundaryOrientation::AnyNonzero
    );
    if (!roundedAudit.strictBoundaryValid())
    {
        throw std::runtime_error(
            "binary STL float conversion invalidates the manufacturing boundary"
        );
    }

    if (geometry.report.wallClearance.applicable)
    {
        const double slack =
            geometry.report.wallClearance.measuredMinimumClearance -
            geometry.report.wallClearance.requestedThickness;
        if (!(slack > 2.0 * maximumError))
        {
            throw std::runtime_error(
                "binary STL float conversion exhausts wall-clearance slack"
            );
        }
    }
}


} // namespace


bool BoundaryAudit::strictBoundaryValid() const noexcept
{
    return
        !empty &&
        matrixShapeValid &&
        finiteVertices &&
        faceIndicesValid &&
        nondegenerateFaces &&
        unreferencedVertexCount == 0 &&
        edgeManifold &&
        vertexManifold &&
        closed &&
        consistentlyOriented &&
        componentOrientationsValid &&
        selfIntersections.completed &&
        selfIntersections.intersectionFree;
}


BoundaryAudit auditTriangleBoundary(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces,
    BoundaryOrientation requiredOrientation,
    bool allowEmpty,
    double geometricTolerance
)
{
    BoundaryAudit result;
    result.vertexCount = static_cast<int>(vertices.rows());
    result.faceCount = static_cast<int>(faces.rows());
    result.empty = faces.rows() == 0;
    result.matrixShapeValid = vertices.cols() == 3 && faces.cols() == 3;
    result.finiteVertices = result.matrixShapeValid && vertices.allFinite();

    if (result.empty)
    {
        const bool validEmpty = allowEmpty && result.matrixShapeValid &&
            result.finiteVertices && vertices.rows() == 0;
        result.faceIndicesValid = validEmpty;
        result.nondegenerateFaces = validEmpty;
        result.edgeManifold = validEmpty;
        result.vertexManifold = validEmpty;
        result.closed = validEmpty;
        result.consistentlyOriented = validEmpty;
        result.componentOrientationsValid = validEmpty;
        result.selfIntersections.completed = validEmpty;
        result.selfIntersections.intersectionFree = validEmpty;
        return result;
    }

    if (!result.matrixShapeValid || !result.finiteVertices)
    {
        return result;
    }

    const double tolerance = resolvedTolerance(vertices, geometricTolerance);
    const double squaredAreaTolerance =
        tolerance * tolerance * coordinateScale(vertices) * coordinateScale(vertices);

    std::vector<unsigned char> referenced(
        static_cast<std::size_t>(vertices.rows()),
        0U
    );
    std::map<std::pair<int, int>, EdgeUse> edges;
    result.faceIndicesValid = true;

    for (Eigen::Index row = 0; row < faces.rows(); ++row)
    {
        const int indices[3] =
        {
            faces(row, 0), faces(row, 1), faces(row, 2)
        };
        bool faceValid = true;
        for (const int index : indices)
        {
            if (index < 0 || index >= vertices.rows())
            {
                faceValid = false;
                result.faceIndicesValid = false;
            }
            else
            {
                referenced[static_cast<std::size_t>(index)] = 1U;
            }
        }

        if (
            !faceValid ||
            indices[0] == indices[1] ||
            indices[1] == indices[2] ||
            indices[2] == indices[0]
        )
        {
            ++result.degenerateFaceCount;
            continue;
        }

        const Eigen::Vector3d first = vertices.row(indices[0]).transpose();
        const Eigen::Vector3d second = vertices.row(indices[1]).transpose();
        const Eigen::Vector3d third = vertices.row(indices[2]).transpose();
        if ((second - first).cross(third - first).squaredNorm() <=
            squaredAreaTolerance)
        {
            ++result.degenerateFaceCount;
        }

        for (int edge = 0; edge < 3; ++edge)
        {
            const int from = indices[edge];
            const int to = indices[(edge + 1) % 3];
            EdgeUse& use = edges[std::minmax(from, to)];
            if (from < to)
            {
                ++use.forward;
            }
            else
            {
                ++use.reverse;
            }
        }
    }

    result.referencedVertexCount = static_cast<int>(std::count(
        referenced.begin(), referenced.end(), static_cast<unsigned char>(1U)
    ));
    result.unreferencedVertexCount =
        result.vertexCount - result.referencedVertexCount;
    result.nondegenerateFaces = result.degenerateFaceCount == 0;

    if (!result.faceIndicesValid)
    {
        return result;
    }

    for (const auto& entry : edges)
    {
        const int total = entry.second.forward + entry.second.reverse;
        if (total == 1)
        {
            ++result.boundaryEdgeCount;
        }
        else if (total > 2)
        {
            ++result.nonManifoldEdgeCount;
        }

        if (total == 2)
        {
            if (entry.second.forward != 1 || entry.second.reverse != 1)
            {
                ++result.inconsistentEdgeCount;
            }
        }
        else if (entry.second.forward != entry.second.reverse)
        {
            ++result.inconsistentEdgeCount;
        }
    }

    result.edgeManifold = result.nonManifoldEdgeCount == 0;
    result.closed =
        result.boundaryEdgeCount == 0 && result.nonManifoldEdgeCount == 0;
    result.consistentlyOriented = result.inconsistentEdgeCount == 0;

    Eigen::VectorXi manifoldVertices;
    result.vertexManifold = igl::is_vertex_manifold(faces, manifoldVertices);
    result.nonManifoldVertexCount = static_cast<int>(
        (manifoldVertices.array() == 0).count()
    );

    Eigen::VectorXi components;
    result.faceComponentCount = igl::facet_components(faces, components);
    result.componentSignedVolumes.assign(
        static_cast<std::size_t>(result.faceComponentCount),
        0.0
    );

    for (int component = 0;
         component < result.faceComponentCount;
         ++component)
    {
        const Eigen::Index faceCount =
            (components.array() == component).count();
        Eigen::MatrixXi componentFaces(faceCount, 3);
        Eigen::Index destination = 0;
        for (Eigen::Index row = 0; row < faces.rows(); ++row)
        {
            if (components(row) == component)
            {
                componentFaces.row(destination++) = faces.row(row);
            }
        }
        result.componentSignedVolumes[static_cast<std::size_t>(component)] =
            computeSurfaceMassProperties(vertices, componentFaces).mass;
    }

    result.signedVolume = std::accumulate(
        result.componentSignedVolumes.begin(),
        result.componentSignedVolumes.end(),
        0.0
    );

    const double scale = coordinateScale(vertices);
    const double volumeTolerance = std::max(
        tolerance * tolerance * tolerance,
        4096.0 * std::numeric_limits<double>::epsilon() *
            scale * scale * scale
    );
    result.componentOrientationsValid = result.closed &&
        result.consistentlyOriented &&
        std::all_of(
            result.componentSignedVolumes.begin(),
            result.componentSignedVolumes.end(),
            [requiredOrientation, volumeTolerance](double volume)
            {
                if (requiredOrientation == BoundaryOrientation::OutwardPositive)
                {
                    return volume > volumeTolerance;
                }
                if (requiredOrientation == BoundaryOrientation::InwardNegative)
                {
                    return volume < -volumeTolerance;
                }
                return std::abs(volume) > volumeTolerance;
            }
        );

    if (result.nondegenerateFaces)
    {
        result.selfIntersections = auditSelfIntersections(
            vertices,
            faces,
            tolerance
        );
    }
    return result;
}


ManufacturingGeometry prepareManufacturingGeometry(
    const deformation::DeformedGeometry& geometry,
    const VoxelGrid& voxelGrid,
    const ManufacturingOptions& options
)
{
    if (!std::isfinite(options.geometricTolerance) ||
        options.geometricTolerance < 0.0)
    {
        throw std::invalid_argument(
            "manufacturing geometric tolerance must be finite and nonnegative"
        );
    }

    ManufacturingGeometry result;
    ManufacturingReadinessReport& report = result.report;

    Eigen::MatrixXd toleranceVertices(0, 3);
    if (
        geometry.outerVertices.cols() == 3 &&
        geometry.outerVertices.allFinite()
    )
    {
        toleranceVertices = geometry.outerVertices;
    }
    if (
        geometry.gridCorners.rows() > 0 &&
        geometry.gridCorners.cols() == 3 &&
        geometry.gridCorners.allFinite()
    )
    {
        const Eigen::Index oldRows = toleranceVertices.rows();
        toleranceVertices.conservativeResize(
            oldRows + geometry.gridCorners.rows(),
            3
        );
        toleranceVertices.bottomRows(geometry.gridCorners.rows()) =
            geometry.gridCorners;
    }
    const double tolerance = resolvedTolerance(
        toleranceVertices,
        options.geometricTolerance
    );

    report.synchronizedInnerGeometry = synchronizedInnerGeometry(
        geometry,
        voxelGrid,
        tolerance
    );

    bool outerIndexShapeValid =
        geometry.outerVertices.cols() == 3 &&
        geometry.outerFaces.cols() == 3 &&
        geometry.outerVertices.allFinite() &&
        (
            geometry.outerFaces.size() == 0 ||
            (
                geometry.outerFaces.minCoeff() >= 0 &&
                geometry.outerFaces.maxCoeff() < geometry.outerVertices.rows()
            )
        );
    bool innerIndexShapeValid =
        geometry.innerVertices.cols() == 3 &&
        geometry.innerFaces.cols() == 3 &&
        geometry.innerVertices.allFinite() &&
        (
            geometry.innerFaces.size() == 0 ||
            (
                geometry.innerFaces.minCoeff() >= 0 &&
                geometry.innerFaces.maxCoeff() < geometry.innerVertices.rows()
            )
        );

    Eigen::MatrixXd compactOuterVertices(0, 3);
    Eigen::MatrixXi compactOuterFaces(0, 3);
    Eigen::MatrixXd compactInnerVertices(0, 3);
    Eigen::MatrixXi compactInnerFaces(0, 3);

    if (outerIndexShapeValid)
    {
        compactMesh(
            geometry.outerVertices,
            geometry.outerFaces,
            compactOuterVertices,
            compactOuterFaces
        );
    }
    if (innerIndexShapeValid)
    {
        compactMesh(
            geometry.innerVertices,
            geometry.innerFaces,
            compactInnerVertices,
            compactInnerFaces
        );
    }

    result.boundary = makeCoherentBoundary(
        compactOuterVertices,
        compactOuterFaces,
        compactInnerVertices,
        compactInnerFaces
    );

    report.outerBoundary = auditTriangleBoundary(
        compactOuterVertices,
        compactOuterFaces,
        BoundaryOrientation::OutwardPositive,
        false,
        tolerance
    );
    report.innerBoundary = auditTriangleBoundary(
        compactInnerVertices,
        compactInnerFaces,
        BoundaryOrientation::InwardNegative,
        true,
        tolerance
    );
    report.coherentBoundary = auditTriangleBoundary(
        result.boundary.vertices,
        result.boundary.faces,
        BoundaryOrientation::AnyNonzero,
        false,
        tolerance
    );

    report.voxelTopology = analyzeVoxelTopology(voxelGrid);
    report.singleMaterialComponent =
        report.voxelTopology.materialVoxelCount > 0 &&
        report.voxelTopology.materialComponentCount == 1;

    try
    {
        const deformation::GridDeformationDomain domain =
            deformation::GridDeformationDomain::fromVoxelGrid(voxelGrid);
        report.deformedCells = auditDeformedCells(
            domain,
            geometry.gridCorners,
            tolerance
        );
        const deformation::DeformedVoxelMassField field(
            domain,
            geometry.gridCorners
        );
        report.deformedCellsValid =
            field.activeCellCount() == domain.activeCellCount() &&
            report.deformedCells.allActiveCellsValid;

        report.surfaceMassProperties = computeSurfaceMassProperties(
            result.boundary.vertices,
            result.boundary.faces
        );
        report.cellOracleMassProperties = computeSurfaceMassProperties(
            compactOuterVertices,
            compactOuterFaces
        );

        const Eigen::Vector3i dimensions = voxelGrid.dimensions();
        for (int x = 0; x < dimensions(0); ++x)
        {
            for (int y = 0; y < dimensions(1); ++y)
            {
                for (int z = 0; z < dimensions(2); ++z)
                {
                    const Voxel& voxel = voxelGrid(x, y, z);
                    if (voxel.isInner && !voxel.isFilled)
                    {
                        const int cell = domain.cellLinearIndex(x, y, z);
                        report.cellOracleMassProperties -=
                            field.massProperties(cell);
                    }
                }
            }
        }

        report.massRelativeError = std::abs(
            report.surfaceMassProperties.mass -
            report.cellOracleMassProperties.mass
        ) / std::max(
            {
                1.0,
                std::abs(report.surfaceMassProperties.mass),
                std::abs(report.cellOracleMassProperties.mass)
            }
        );
        report.firstMomentRelativeError = (
            report.surfaceMassProperties.firstMoment -
            report.cellOracleMassProperties.firstMoment
        ).norm() / std::max(
            {
                1.0,
                report.surfaceMassProperties.firstMoment.norm(),
                report.cellOracleMassProperties.firstMoment.norm()
            }
        );
        report.massPropertiesConsistent =
            report.massRelativeError <= 1.0e-10 &&
            report.firstMomentRelativeError <= 1.0e-10;
        const bool surfaceCenterValid =
            report.surfaceMassProperties.tryCenterOfMass(
            report.centerOfMass,
            0.0
        );
        const bool oracleCenterValid =
            report.cellOracleMassProperties.tryCenterOfMass(
                report.cellOracleCenterOfMass,
                0.0
            );
        report.centerOfMassRelativeError =
            surfaceCenterValid && oracleCenterValid
                ? (report.centerOfMass - report.cellOracleCenterOfMass).norm() /
                    std::max(
                        {
                            1.0,
                            report.centerOfMass.norm(),
                            report.cellOracleCenterOfMass.norm()
                        }
                    )
                : std::numeric_limits<double>::infinity();
        report.massPropertiesConsistent =
            report.massPropertiesConsistent &&
            report.centerOfMassRelativeError <= 1.0e-10;
    }
    catch (const std::exception& error)
    {
        report.deformedCellsValid = false;
        report.deformedCellFailure = error.what();
    }

    if (
        report.outerBoundary.strictBoundaryValid() &&
        (report.innerBoundary.empty ||
         report.innerBoundary.strictBoundaryValid())
    )
    {
        report.wallClearance = auditWallClearance(
            compactOuterVertices,
            compactOuterFaces,
            compactInnerVertices,
            compactInnerFaces,
            voxelGrid.minimumWallThickness(),
            tolerance
        );
        report.containment = auditContainment(
            compactOuterVertices,
            compactOuterFaces,
            compactInnerVertices,
            compactInnerFaces,
            report.innerBoundary,
            report.wallClearance.intersectionFree
        );
    }
    else
    {
        report.wallClearance.requestedThickness =
            voxelGrid.minimumWallThickness();
    }

    report.signedMaterialVolume =
        report.outerBoundary.signedVolume +
        report.innerBoundary.signedVolume;
    const double scale = coordinateScale(result.boundary.vertices);
    const double materialVolumeTolerance =
        4096.0 * std::numeric_limits<double>::epsilon() *
        scale * scale * scale;
    report.positiveMaterialVolume =
        std::isfinite(report.signedMaterialVolume) &&
        report.signedMaterialVolume > materialVolumeTolerance;

    if (!report.synchronizedInnerGeometry)
    {
        report.failures.push_back("inner_geometry_not_synchronized");
    }
    if (!report.outerBoundary.strictBoundaryValid())
    {
        report.failures.push_back("outer_boundary_invalid");
    }
    else if (report.outerBoundary.faceComponentCount != 1)
    {
        report.failures.push_back("outer_boundary_not_single_component");
    }
    if (!report.innerBoundary.empty &&
        !report.innerBoundary.strictBoundaryValid())
    {
        report.failures.push_back("inner_boundary_invalid");
    }
    if (!report.coherentBoundary.strictBoundaryValid())
    {
        report.failures.push_back("coherent_boundary_invalid");
    }
    if (!report.deformedCellsValid)
    {
        report.failures.push_back("deformed_cells_invalid");
    }
    if (!report.massPropertiesConsistent)
    {
        report.failures.push_back("surface_cell_mass_properties_disagree");
    }
    if (options.requireSingleMaterialComponent &&
        !report.singleMaterialComponent)
    {
        report.failures.push_back("material_not_single_component");
    }
    if (!report.positiveMaterialVolume)
    {
        report.failures.push_back("material_volume_not_positive");
    }
    if (!report.wallClearance.completed)
    {
        report.failures.push_back("wall_clearance_not_completed");
    }
    else if (!report.wallClearance.intersectionFree)
    {
        report.failures.push_back("outer_inner_intersection");
    }
    else if (!report.wallClearance.meetsRequestedThickness)
    {
        report.failures.push_back("deformed_wall_below_requested_thickness");
    }
    if (!report.containment.completed)
    {
        report.failures.push_back("cavity_containment_not_completed");
    }
    else if (!report.containment.allCavitiesInsideOuter)
    {
        report.failures.push_back("cavity_not_contained_by_outer_boundary");
    }

    report.ready = report.failures.empty();
    return result;
}


std::string ManufacturingReadinessReport::toJson() const
{
    std::ostringstream output;
    output << std::boolalpha << "{\n";
    output << "  \"schema\": \"make-it-stand-r8-readiness-v1\",\n";
    output << "  \"ready\": " << jsonBool(ready) << ",\n";
    output << "  \"synchronized_inner_geometry\": "
        << jsonBool(synchronizedInnerGeometry) << ",\n";
    output << "  \"deformed_cells_valid\": "
        << jsonBool(deformedCellsValid) << ",\n";
    output << "  \"single_material_component\": "
        << jsonBool(singleMaterialComponent) << ",\n";
    output << "  \"positive_material_volume\": "
        << jsonBool(positiveMaterialVolume) << ",\n";
    output << "  \"mass_properties_consistent\": "
        << jsonBool(massPropertiesConsistent) << ",\n";
    output << "  \"signed_material_volume\": "
        << std::setprecision(17) << signedMaterialVolume << ",\n";
    output << "  \"mass_relative_error\": ";
    appendJsonNumber(output, massRelativeError);
    output << ",\n";
    output << "  \"first_moment_relative_error\": ";
    appendJsonNumber(output, firstMomentRelativeError);
    output << ",\n";
    output << "  \"center_of_mass_relative_error\": ";
    appendJsonNumber(output, centerOfMassRelativeError);
    output << ",\n";
    output << "  \"center_of_mass\": ["
        << std::setprecision(17) << centerOfMass(0) << ", "
        << centerOfMass(1) << ", " << centerOfMass(2) << "],\n";
    appendBoundaryJson(output, "outer_boundary", outerBoundary);
    appendBoundaryJson(output, "inner_boundary", innerBoundary);
    appendBoundaryJson(output, "coherent_boundary", coherentBoundary);
    output << "  \"wall_clearance\": {\n";
    output << "    \"applicable\": " << jsonBool(wallClearance.applicable)
        << ",\n";
    output << "    \"completed\": " << jsonBool(wallClearance.completed)
        << ",\n";
    output << "    \"intersection_free\": "
        << jsonBool(wallClearance.intersectionFree) << ",\n";
    output << "    \"meets_requested_thickness\": "
        << jsonBool(wallClearance.meetsRequestedThickness) << ",\n";
    output << "    \"requested_thickness\": " << std::setprecision(17)
        << wallClearance.requestedThickness << ",\n";
    output << "    \"measured_minimum_clearance\": ";
    if (std::isfinite(wallClearance.measuredMinimumClearance))
    {
        output << std::setprecision(17)
            << wallClearance.measuredMinimumClearance;
    }
    else
    {
        output << "null";
    }
    output << "\n  },\n";
    output << "  \"containment\": {\n";
    output << "    \"applicable\": " << jsonBool(containment.applicable)
        << ",\n";
    output << "    \"completed\": " << jsonBool(containment.completed)
        << ",\n";
    output << "    \"all_cavities_inside_outer\": "
        << jsonBool(containment.allCavitiesInsideOuter) << ",\n";
    output << "    \"tested_components\": "
        << containment.testedComponentCount << ",\n";
    output << "    \"failed_component\": "
        << containment.failedComponent << ",\n";
    output << "    \"nested_cavity_pairs\": "
        << containment.nestedCavityPairCount << ",\n";
    output << "    \"first_containing_component\": "
        << containment.firstContainingComponent << ",\n";
    output << "    \"first_nested_component\": "
        << containment.firstNestedComponent << "\n";
    output << "  },\n";
    output << "  \"deformed_cells\": {\n";
    output << "    \"completed\": " << jsonBool(deformedCells.completed)
        << ",\n";
    output << "    \"all_active_cells_valid\": "
        << jsonBool(deformedCells.allActiveCellsValid) << ",\n";
    output << "    \"active_cells\": "
        << deformedCells.activeCellCount << ",\n";
    output << "    \"invalid_cells\": "
        << deformedCells.invalidCellCount << ",\n";
    output << "    \"first_invalid_cell\": "
        << deformedCells.firstInvalidCell << ",\n";
    output << "    \"degenerate_boundary_triangles\": "
        << deformedCells.degenerateBoundaryTriangleCount << ",\n";
    output << "    \"nonpositive_jacobian_samples\": "
        << deformedCells.nonpositiveJacobianSampleCount << ",\n";
    output << "    \"minimum_jacobian_determinant\": ";
    if (std::isfinite(deformedCells.minimumJacobianDeterminant))
    {
        output << std::setprecision(17)
            << deformedCells.minimumJacobianDeterminant;
    }
    else
    {
        output << "null";
    }
    output << "\n  },\n";
    output << "  \"voxel_topology\": {\n";
    output << "    \"material_components\": "
        << voxelTopology.materialComponentCount << ",\n";
    output << "    \"material_voxels\": "
        << voxelTopology.materialVoxelCount << ",\n";
    output << "    \"void_components\": "
        << voxelTopology.carvedVoidComponentCount << ",\n";
    output << "    \"void_voxels\": "
        << voxelTopology.carvedVoidVoxelCount << "\n";
    output << "  },\n";
    output << "  \"deformed_cell_failure\": \""
        << escapedJson(deformedCellFailure) << "\",\n";
    output << "  \"failures\": [";
    for (std::size_t index = 0; index < failures.size(); ++index)
    {
        if (index > 0)
        {
            output << ", ";
        }
        output << "\"" << escapedJson(failures[index]) << "\"";
    }
    output << "]\n}\n";
    return output.str();
}


void exportReadyManufacturingMesh(
    const ManufacturingGeometry& geometry,
    const std::string& filename,
    ManufacturingMeshFormat format
)
{
    if (!geometry.report.ready)
    {
        throw std::logic_error(
            "manufacturing mesh export requires a passing readiness report"
        );
    }
    std::vector<unsigned char> bytes;
    if (format == ManufacturingMeshFormat::Obj)
    {
        bytes = stringBytes(serializeObj(geometry.boundary));
    }
    else if (format == ManufacturingMeshFormat::AsciiStl)
    {
        bytes = stringBytes(serializeAsciiStl(geometry.boundary));
    }
    else if (format == ManufacturingMeshFormat::BinaryStl)
    {
        validateBinaryFloatGeometry(geometry);
        bytes = serializeBinaryStl(geometry.boundary);
    }
    else
    {
        throw std::invalid_argument("unsupported manufacturing mesh format");
    }
    writeAtomicBytes(filename, bytes);
}


void writeManufacturingReadinessReport(
    const ManufacturingReadinessReport& report,
    const std::string& filename
)
{
    writeAtomicBytes(filename, stringBytes(report.toJson()));
}


} // namespace manufacturing
