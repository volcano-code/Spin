#pragma once
#include "physics/RawMassMoments10.h"
#include <array>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#if MAKE_IT_STAND_HAS_SPIN_IT_2014 && MAKE_IT_STAND_HAS_EXACT_GEOMETRY_CHECKS
namespace spin2014 {

struct SourceSurfaceOctreeConfig {
    Eigen::Vector3d rootOrigin = Eigen::Vector3d::Zero();
    double rootSide = 1.0;
    int maximumDepth = 9; // Supported: 0..20.
    std::size_t maximumNodes = 5000000; // Input guard, included in config identity.
};

enum class SourceSurfaceCellKind : std::uint8_t { Internal, Boundary, Full, Empty };
struct SourceSurfaceCellBounds {
    Eigen::Vector3d lower,upper;
};
struct SourceSurfaceCell {
    static constexpr std::size_t noIndex = std::numeric_limits<std::size_t>::max();
    std::size_t index = noIndex, parent = noIndex, firstChild = noIndex;
    std::array<std::uint32_t,3> coordinates{};
    int level = 0;
    SourceSurfaceCellKind kind = SourceSurfaceCellKind::Empty;
    physics::RawMassMoments10 moments;
    double unresolvedVolume = 0;
};
struct SourceSurfaceOctreeStats {
    std::size_t internalNodes=0, boundaryNodes=0, fullNodes=0, emptyNodes=0;
    int maximumVisitedLevel=0;
    std::vector<std::size_t> levels;
};

/** Immutable geometry oracle built exclusively from copied original surfaces.
 * Requires both 2014 and exact-geometry features, independently of Faster.
 * Exact rational audits establish signed material legality and design subset
 * outer. Surface moments integrate the unapproximated original polyhedra with
 * binary64 Algorithm 1 arithmetic, not rational moment arithmetic.
 *
 * IMPLEMENTATION_DEFAULT: finest CLOSED cells touching any design triangle
 * contribute zero carveable moments and their entire volume as uncertainty.
 * The resulting cuboid sum is a conservative geometry approximation. Uncertainty
 * is a volume bound, not a ten-moment relative-error or feasibility certificate.
 * The finest step must strictly exceed the maximum root-endpoint ULP so all
 * represented dyadic planes remain distinct. Equal-ULP inputs are conservatively
 * rejected as numerically unsupported even when a particular alignment works;
 * the implementation does not recenter, modify geometry, or coarsen the request.
 * No occupancy input, optimizer, repair, or mesh/cage export is provided.
 */
class SourceSurfaceOctree {
public:
    SourceSurfaceOctree(const Eigen::MatrixXd& outerVertices,
                        const Eigen::MatrixXi& outerFaces,
                        const Eigen::MatrixXd& designVertices,
                        const Eigen::MatrixXi& designFaces,
                        const SourceSurfaceOctreeConfig& config = {});
    const Eigen::MatrixXd& outerVertices() const;
    const Eigen::MatrixXi& outerFaces() const;
    const Eigen::MatrixXd& designVertices() const;
    const Eigen::MatrixXi& designFaces() const;
    const SourceSurfaceOctreeConfig& config() const;
    const physics::RawMassMoments10& outerSurfaceMoments() const;
    const physics::RawMassMoments10& designSurfaceMoments() const;
    const physics::RawMassMoments10& carveableMoments() const;
    double unresolvedBoundaryVolume() const;
    const std::string& sourceSha256() const;
    const std::string& configSha256() const;
    const std::string& partitionSha256() const;
    const SourceSurfaceOctreeStats& stats() const;
    const std::vector<SourceSurfaceCell>& cells() const;
    /** Stored internal children occupy [firstChild, firstChild+8), x bit first.
     * Descendants of full/empty terminals are analytic virtual cells: index is
     * noIndex, parent identifies the covering stored terminal, no allocation.
     * All other returned cells exactly match stored scientific contributions.
     */
    SourceSurfaceCell cell(int level, std::array<std::uint32_t,3> coordinates) const;
    /** Canonical represented dyadic planes used by geometry classification and
     * moments, including virtual descendants. Each plane is rounded once by
     * FMA so coincident boundaries agree across levels and all consumers.
     */
    SourceSurfaceCellBounds cellBounds(int level,
        std::array<std::uint32_t,3> coordinates) const;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
};
} // namespace spin2014
#endif
