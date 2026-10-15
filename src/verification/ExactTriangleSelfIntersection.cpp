#include "verification/ExactTriangleSelfIntersection.h"

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/properties.h>

#include <algorithm>
#include <map>
#include <vector>

namespace verification {
namespace {

namespace PMP = CGAL::Polygon_mesh_processing;
using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Point = Kernel::Point_3;
using Mesh = CGAL::Surface_mesh<Point>;
using Polygon = std::vector<std::size_t>;

ExactTriangleSelfIntersectionAudit failure(std::string status, std::string diagnostic) {
    ExactTriangleSelfIntersectionAudit result;
    result.status = std::move(status);
    result.diagnostic = std::move(diagnostic);
    return result;
}

} // namespace

ExactTriangleSelfIntersectionAudit auditExactTriangleMeshSelfIntersections(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces) {
    if (vertices.cols() != 3 || faces.cols() != 3 || vertices.rows() < 4 ||
        faces.rows() < 4 || !vertices.allFinite()) {
        return failure("INVALID_OR_OPEN_INPUT", "expected finite 3D vertices and triangle faces");
    }

    std::map<std::pair<int, int>, int> edgeIncidence;
    std::vector<Point> points;
    points.reserve(static_cast<std::size_t>(vertices.rows()));
    for (Eigen::Index row = 0; row < vertices.rows(); ++row) {
        points.emplace_back(vertices(row, 0), vertices(row, 1), vertices(row, 2));
    }
    std::vector<Polygon> polygons;
    polygons.reserve(static_cast<std::size_t>(faces.rows()));
    for (Eigen::Index row = 0; row < faces.rows(); ++row) {
        Polygon polygon;
        polygon.reserve(3);
        for (int column = 0; column < 3; ++column) {
            const int vertex = faces(row, column);
            if (vertex < 0 || vertex >= vertices.rows()) {
                return failure("INVALID_OR_OPEN_INPUT", "face index is outside the vertex array");
            }
            polygon.push_back(static_cast<std::size_t>(vertex));
        }
        if (polygon[0] == polygon[1] || polygon[1] == polygon[2] || polygon[2] == polygon[0]) {
            return failure("INVALID_OR_OPEN_INPUT", "triangle repeats a vertex index");
        }
        for (int edge = 0; edge < 3; ++edge) {
            const int first = static_cast<int>(polygon[static_cast<std::size_t>(edge)]);
            const int second = static_cast<int>(polygon[static_cast<std::size_t>((edge + 1) % 3)]);
            ++edgeIncidence[{std::min(first, second), std::max(first, second)}];
        }
        polygons.push_back(std::move(polygon));
    }
    if (std::any_of(edgeIncidence.begin(), edgeIncidence.end(),
                    [](const auto& item) { return item.second != 2; })) {
        return failure("INVALID_OR_OPEN_INPUT", "input is not a closed edge-manifold triangle mesh");
    }

    const std::size_t originalVertexCount = points.size();
    if (!PMP::orient_polygon_soup(points, polygons) || points.size() != originalVertexCount) {
        return failure("UNKNOWN_NONMANIFOLD_OR_NONORIENTABLE",
                       "CGAL could not orient the soup without splitting geometric vertices");
    }
    if (!PMP::is_polygon_soup_a_polygon_mesh(polygons)) {
        return failure("UNKNOWN_NONMANIFOLD_OR_NONORIENTABLE",
                       "oriented input is not a polygon mesh");
    }

    Mesh mesh;
    PMP::polygon_soup_to_polygon_mesh(points, polygons, mesh);
    if (!CGAL::is_valid_polygon_mesh(mesh) || !CGAL::is_triangle_mesh(mesh) ||
        !CGAL::is_closed(mesh)) {
        return failure("UNKNOWN_NONMANIFOLD_OR_NONORIENTABLE",
                       "CGAL conversion did not preserve a closed triangle mesh");
    }

    ExactTriangleSelfIntersectionAudit result;
    result.checked_face_count = static_cast<int>(faces.rows());
    std::vector<std::pair<Mesh::Face_index, Mesh::Face_index>> intersecting;
    PMP::self_intersections(mesh, std::back_inserter(intersecting));
    const auto faceIndex = get(CGAL::face_index, mesh);
    result.intersecting_face_pairs.reserve(intersecting.size());
    for (const auto& pair : intersecting) {
        result.intersecting_face_pairs.emplace_back(
            static_cast<int>(get(faceIndex, pair.first)),
            static_cast<int>(get(faceIndex, pair.second)));
    }
    result.status = result.intersecting_face_pairs.empty()
        ? "INTERSECTION_FREE" : "INTERSECTIONS_FOUND";
    result.diagnostic = "CGAL exact-predicate triangle self-intersection query; coordinates are exact binary64 values";
    return result;
}

} // namespace verification
