#pragma once
#include "verification/ExactBoundaryAudit.h"
#include <array>
#include <cstdint>
#include <memory>
#include <stdexcept>

#if MAKE_IT_STAND_HAS_SPIN_IT_2014 && MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
namespace spin2014 {
struct SourceSurfaceCageConfig {
    Eigen::Vector3d rootOrigin = Eigen::Vector3d::Zero();
    double rootSide = 1.0;
    int cageLevel = 3; // Supported: 0..20. No automatic level reduction.
    std::size_t maximumSourceNodes = 5000000;
    std::size_t maximumSelectedCells = 1000000;
    std::size_t maximumBoundaryVertices = 1000000;
    std::size_t maximumBoundaryTriangles = 2000000;
};

class SourceSurfaceCageGeometryError : public std::runtime_error {
public:
    explicit SourceSurfaceCageGeometryError(verification::ExactBoundaryAudit audit);
    const verification::ExactBoundaryAudit& boundaryAudit() const { return audit_; }
private:
    verification::ExactBoundaryAudit audit_;
};

/** Complete source-surface cage GEOMETRY prerequisite, independently of Faster.
 * IMPLEMENTATION_DEFAULT: CLOSED cells touching the copied original material
 * boundary count, even at zero-volume contacts. This is a conservative envelope,
 * not an exact material integral or an author-specified tie rule.
 * All topological boundary components (including cavities) and all their grid
 * corners are retained, without repair, smoothing, filling, culling or sampling.
 * Keys are lexicographic (x,y,z); local corners use x bit first. Handles map
 * one-to-one to indexed boundary vertices in that same canonical key order.
 * Triangle diagonals and outward union winding are fixed, inward at cavities.
 * Singular boundaries throw SourceSurfaceCageGeometryError with exact audit;
 * resource exhaustion throws length_error with RESOURCE_LIMIT, never infeasible.
 * No usable partial result, fallback level, BBW, optimization, physics, seal,
 * feasibility or original-reproduction authority is provided.
 * Source identity reuses the unchanged outer-as-design source octree identity.
 */
class SourceSurfaceCage {
public:
    using Key = std::array<std::uint32_t,3>;
    SourceSurfaceCage(const Eigen::MatrixXd& sourceVertices,
                      const Eigen::MatrixXi& sourceFaces,
                      const SourceSurfaceCageConfig& config);
    const Eigen::MatrixXd& sourceVertices() const;
    const Eigen::MatrixXi& sourceFaces() const;
    const SourceSurfaceCageConfig& config() const;
    const std::vector<Key>& selectedCells() const;
    const std::vector<Key>& boundaryVertexKeys() const;
    const Eigen::MatrixXd& boundaryVertices() const;
    const Eigen::MatrixXi& boundaryFaces() const;
    const std::vector<int>& handleVertexIndices() const;
    const verification::ExactBoundaryAudit& boundaryAudit() const;
    std::size_t sourceNodeCount() const;
    const std::string& sourceSha256() const;
    const std::string& configSha256() const;
    const std::string& outputSha256() const;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
};
} // namespace spin2014
#endif
