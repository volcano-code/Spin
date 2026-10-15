#pragma once
#include "RankTwoCoveredFiniteMesh.h"
#include "RankTwoMovingInterior.h"

namespace spin {
class RankTwoMovingMeshEdgePoint;
struct RankTwoInteriorMeshSettings {
    std::uint64_t maximumPredicates=1000000;
};
struct RankTwoInteriorMeshPoint {
    std::string identity;
    std::array<detail::CertifiedAlgebraicReal,2> coordinates;
};
// The original exact constrained mesh is retained. Only incident faces are
// subdivided using privately certified moving points. Original shared edges,
// native ownership and holes are invariant. This is NOT a thickness grant.
class RankTwoInteriorMesh {
public:
    const RankTwoCoveredFiniteMesh& baseMesh()const;
    const std::vector<RankTwoInteriorMeshPoint>& points()const;
    const std::vector<RankTwoFiniteTriangle>& triangles()const;
    const std::vector<RankTwoRegionTriangulation>& regionSummaries()const;
    const std::string& canonicalRecord()const;
    const std::string& scientificHash()const;
private:
    struct Impl;
    explicit RankTwoInteriorMesh(std::shared_ptr<const Impl>);
    std::shared_ptr<const Impl> impl_;
    friend RankTwoInteriorMesh insertRankTwoMovingInteriorNodes(const RankTwoArrangement&,
        const RankTwoCoveredFiniteMesh&,const std::vector<RankTwoMovingInteriorPoint>&,
        const RankTwoInteriorMeshSettings&);
    friend RankTwoInteriorMesh insertRankTwoMovingEdgeNodes(const RankTwoArrangement&,
        const RankTwoInteriorMesh&,const std::vector<RankTwoMovingMeshEdgePoint>&,
        const RankTwoInteriorMeshSettings&);
};
RankTwoInteriorMesh insertRankTwoMovingInteriorNodes(const RankTwoArrangement&,
    const RankTwoCoveredFiniteMesh&,const std::vector<RankTwoMovingInteriorPoint>&,
    const RankTwoInteriorMeshSettings& settings={});
} // namespace spin
