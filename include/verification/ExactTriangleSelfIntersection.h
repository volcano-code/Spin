#pragma once

#include <Eigen/Core>

#include <string>
#include <utility>
#include <vector>

namespace verification {

/** Exact-predicate self-intersection audit of a closed indexed triangle mesh.
 * Coordinates are interpreted as their exact binary64 values. The audit may
 * reverse face winding to construct an oriented polygon mesh, but never moves,
 * welds, drops, or adds geometric vertices. Unsupported topology is UNKNOWN,
 * not intersection-free.
 */
struct ExactTriangleSelfIntersectionAudit {
    std::string status;
    std::string diagnostic;
    int checked_face_count = 0;
    std::vector<std::pair<int, int>> intersecting_face_pairs;
};

ExactTriangleSelfIntersectionAudit auditExactTriangleMeshSelfIntersections(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces);

} // namespace verification
