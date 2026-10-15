#pragma once

#include <Eigen/Core>
#include <array>
#include <cmath>
#include <map>
#include <utility>
#include <vector>
#include <algorithm>

namespace spin::detail {
struct FixedNormalSphere {
    std::vector<Eigen::Vector3d> vertices;
    std::vector<std::array<int,3>> faces;
};

// Shared original generator. Keep arithmetic, traversal and midpoint insertion
// order unchanged: represented rays (not re-normalized copies) are authority.
// The caller must validate the requested subdivision/resource budget first.
inline FixedNormalSphere fixedNormalSphere(int subdivisions) {
    using Face=std::array<int,3>;
    using Edge=std::pair<int,int>;
    const double phi = (1.0 + std::sqrt(5.0)) / 2.0;
    FixedNormalSphere mesh;
    mesh.vertices = {
        {-1.0, phi, 0.0}, {1.0, phi, 0.0},
        {-1.0, -phi, 0.0}, {1.0, -phi, 0.0},
        {0.0, -1.0, phi}, {0.0, 1.0, phi},
        {0.0, -1.0, -phi}, {0.0, 1.0, -phi},
        {phi, 0.0, -1.0}, {phi, 0.0, 1.0},
        {-phi, 0.0, -1.0}, {-phi, 0.0, 1.0}
    };
    for (Eigen::Vector3d& vertex : mesh.vertices) vertex.normalize();
    mesh.faces = {
        {0, 11, 5}, {0, 5, 1}, {0, 1, 7}, {0, 7, 10}, {0, 10, 11},
        {1, 5, 9}, {5, 11, 4}, {11, 10, 2}, {10, 7, 6}, {7, 1, 8},
        {3, 9, 4}, {3, 4, 2}, {3, 2, 6}, {3, 6, 8}, {3, 8, 9},
        {4, 9, 5}, {2, 4, 11}, {6, 2, 10}, {8, 6, 7}, {9, 8, 1}
    };
    for (int level = 0; level < subdivisions; ++level) {
        std::map<Edge, int> midpointIndices;
        std::vector<Face> refined;
        refined.reserve(mesh.faces.size() * 4);
        const auto midpoint = [&](int first, int second) {
            const Edge edge = std::minmax(first, second);
            const auto found = midpointIndices.find(edge);
            if (found != midpointIndices.end()) return found->second;
            Eigen::Vector3d point =
                (mesh.vertices[first] + mesh.vertices[second]).normalized();
            const int index = static_cast<int>(mesh.vertices.size());
            mesh.vertices.push_back(point);
            midpointIndices.emplace(edge, index);
            return index;
        };
        for (const Face& face : mesh.faces) {
            const int ab = midpoint(face[0], face[1]);
            const int bc = midpoint(face[1], face[2]);
            const int ca = midpoint(face[2], face[0]);
            refined.push_back({face[0], ab, ca});
            refined.push_back({face[1], bc, ab});
            refined.push_back({face[2], ca, bc});
            refined.push_back({ab, bc, ca});
        }
        mesh.faces = std::move(refined);
    }
    return mesh;
}
} // namespace spin::detail
