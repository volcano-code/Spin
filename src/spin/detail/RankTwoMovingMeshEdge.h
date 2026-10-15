#pragma once
#include "RankTwoInteriorMesh.h"

namespace spin {
struct RankTwoMovingMeshEdgeRequest {
    std::array<std::size_t,2> pointIndices{};
    CGAL::Gmpq fraction=CGAL::Gmpq(1,2);
};
struct RankTwoMovingMeshEdgeSettings {
    int maximumIsolationBits=512;
    std::uint64_t maximumNativeChecks=1000000;
};
// A dependency-bearing affine point of an ACTUAL unconstrained internal
// mesh edge. Not a native arc sample, a frozen rational witness, or a midpoint
// of a root isolator. The parent mesh is part of the immutable authority.
class RankTwoMovingMeshEdgePoint {
public:
    const std::string& arrangementIdentity()const;
    const std::string& parentMeshIdentity()const;
    const RankTwoMovingMeshEdgeRequest& definition()const;
    const std::array<detail::CertifiedAlgebraicReal,2>& coordinates()const;
    std::size_t registryRegionIndex()const;
    const std::string& canonicalRecord()const;
    const std::string& scientificHash()const;
private:
    struct Impl;
    explicit RankTwoMovingMeshEdgePoint(std::shared_ptr<const Impl>);
    std::shared_ptr<const Impl> impl_;
    friend RankTwoMovingMeshEdgePoint buildRankTwoMovingMeshEdgePoint(
        const RankTwoArrangement&,const RankTwoInteriorMesh&,
        const RankTwoMovingMeshEdgeRequest&,const RankTwoMovingMeshEdgeSettings&);
};
RankTwoMovingMeshEdgePoint buildRankTwoMovingMeshEdgePoint(const RankTwoArrangement&,
    const RankTwoInteriorMesh&,const RankTwoMovingMeshEdgeRequest&,
    const RankTwoMovingMeshEdgeSettings& settings={});
RankTwoInteriorMesh insertRankTwoMovingEdgeNodes(const RankTwoArrangement&,
    const RankTwoInteriorMesh&,const std::vector<RankTwoMovingMeshEdgePoint>&,
    const RankTwoInteriorMeshSettings& settings={});
} // namespace spin
