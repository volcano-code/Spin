#pragma once
#include <Eigen/Core>
#include <array>
#include <string>

namespace verification::detail {
// Internal helper consumes only authenticated, nonempty, original triangles.
// No Nef exposure or floating-mesh reconstruction. Throws on failure.
struct BoundaryDistance {
    std::string squaredExact;
    std::array<std::size_t,2> witness;
    bool atLeastMinimum, contact;
};
BoundaryDistance exactBoundaryDistance(const Eigen::MatrixXd& aVertices,
    const Eigen::MatrixXi& aFaces, const Eigen::MatrixXd& bVertices,
    const Eigen::MatrixXi& bFaces, double minimumClearance);
} // namespace verification::detail
